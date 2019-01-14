import sys, io, os, random, logging
import scipy.stats as ST
import numpy as np
from jacks.preprocess import subsample_and_preprocess
from jacks.jacks_io import readControlGeneset, collateTestControlSamples, createSampleSpec, createGeneSpec
from jacks.infer import LOG, inferJACKSGene

LOG.setLevel(logging.WARNING)

if len(sys.argv) != 8 and len(sys.argv) != 9:
    print('Usage: sample_jacks_screen.py condensed_input test_line num_replicates(-1 for all)  num_celllines(-1 for all) outfile num_samples num_guides(-1 for all) job_idx\n')
    print('where, condensed_input = countfile#replicatefile:rep_hdr:sample_hdr:ctrl_sample_or_hdr#guidemappingfile:sgrna_hdr:gene_hdr#ctrl_genes(can be blank)')
else:

    #Minimial checks on this, as this is for a script that is intended for use internally only
    condensed_input = sys.argv[1]
    countfile, replicatestuff, grnastuff, ctrl_genes = condensed_input.split('#')
    replicatefile, rep_hdr, sample_hdr, ctrl_sample_or_hdr = replicatestuff.split(':')
    guidemappingfile, sgrna_hdr, gene_hdr = grnastuff.split(':')
    ctrl_sample_hdr = ctrl_sample_or_hdr if ctrl_sample_or_hdr == 'Control' else None
    sample_spec, ctrl_spec, sample_num_reps = createSampleSpec(countfile, replicatefile, rep_hdr, sample_hdr, ctrl_sample_or_hdr, ctrl_sample_hdr)
    gene_spec = createGeneSpec(guidemappingfile, sgrna_hdr, gene_hdr)
    test_celllines = [sample_id for sample_id in ctrl_spec if ctrl_spec[sample_id] != sample_id]

    ctrl_geneset = readControlGeneset(ctrl_genes) if ctrl_genes is not '' else set()
    normtype = 'median'

    test_line = sys.argv[2]
    num_replicates = eval(sys.argv[3])
    num_celllines = eval(sys.argv[4])
    outfile = sys.argv[5]
    num_samples = eval(sys.argv[6])
    num_guides = eval(sys.argv[7])
    if len(sys.argv) == 9: job_idx = sys.argv[8]
    else: job_idx = os.environ['LSB_JOBINDEX']
    
    for bs in range(num_samples):

        if num_celllines < 0: selected_celllines = [test_line] + test_celllines #All cell lines
        else: selected_celllines = [test_line] + [random.choice(test_celllines) for i in range(num_celllines-1)]

        print(selected_celllines)
        # Read all data from input file
        selected_screens = [(ctrl_spec[x],num_ctrl_reps) for x in selected_celllines] + [(x, num_replicates) for x in selected_celllines]
        data, meta, cell_lines, genes, gene_index = subsample_and_preprocess(selected_screens, sample_spec, gene_spec, ctrl_spec=ctrl_spec, ctrl_geneset=ctrl_geneset, normtype=normtype)
        testdata, ctrldata, _ = collateTestControlSamples(data, cell_lines, ctrl_spec)
        if testdata.shape != ctrldata.shape: raise Exception('Mismatch between test and control shape')
        
        #Run JACKS
        if bs == 0: 
            gene_ws = {gene:[] for gene in gene_index}
            bs_errors = np.zeros((data.shape[0],num_samples))
        bs_jacks_results = {}
        idx = selected_celllines.index(test_line)
        for gene in gene_index:
        
            Ig = gene_index[gene]
            # If restricting the number of guides per gene...
            if num_guides > 0 and len(Ig) >= num_guides:
                Ig = random.sample(Ig, num_guides)

            # Perform inference
            y, tau, x1, x2, w1, w2 = inferJACKSGene(testdata[Ig,:,0], testdata[Ig,:,1], ctrldata[Ig,:,0], ctrldata[Ig,:,1], 50)
            bs_jacks_results[gene] = (None, None, None, None, [w1[idx]], [w2[idx]])

        norm_w1 = 0.0
        for gene in gene_index:
            _, _, _, _, w1, w2 = bs_jacks_results[gene]
            w1, w2 = w1[0], w2[0]
            gene_ws[gene].append((w1-norm_w1, np.sqrt(w2 - w1**2.0), (1.0-ST.norm.sf((w1-norm_w1)/np.sqrt(w2 - w1**2.0)))))

    if not os.path.exists(os.path.dirname(outfile)): os.makedirs(os.path.dirname(outfile))  
    fout = io.open(outfile[:-4] + '_%s_' % job_idx + test_line + outfile[-4:], 'w')
    fout.write(u'Gene\t%s\n' % ('\t'.join(['WMean\tWStd\tPval' for i in range(num_samples)])))
    for gene in gene_ws:
        fout.write(u'%s\t%s\n' % (gene, '\t'.join(['%5e\t%5e\t%5e' % (wmean, wstd, pval) for (wmean, wstd, pval) in gene_ws[gene]])))
    fout.close()

