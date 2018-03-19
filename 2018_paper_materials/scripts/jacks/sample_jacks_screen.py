import sys, io, os, random, logging
from jacks.io_preprocess import subsample_and_preprocess
from jacks.jacks import infer_JACKS_gene, LOG
import numpy as np

if len(sys.argv) != 8 and len(sys.argv) != 9:
    print 'Usage: sample_jacks_screen.py rawcountfile test_line num_replicates(-1 for all) num_celllines(-1 for all) outfile num_samples num_guides(-1 for all) job_idx'
else:

    LOG.setLevel(logging.WARNING)

    inputfile = sys.argv[1]
    test_line = sys.argv[2]
    num_replicates = eval(sys.argv[3])
    num_celllines = eval(sys.argv[4])
    outfile = sys.argv[5]
    num_samples = eval(sys.argv[6])
    num_guides = eval(sys.argv[7])
    if len(sys.argv) == 9: job_idx = sys.argv[8]
    else: job_idx = os.environ['LSB_JOBINDEX']
   
    #Fetch the cell lines
    f = io.open(inputfile)
    full_hdrs = [x for x in f.readline()[:-1].split('\t')[2:]]
    hdrs = [x.split('_')[0] for x in full_hdrs]
    test_celllines = [x for x in set(hdrs) if x != test_line]
    f.close()
 
    for bs in range(num_samples):

        if num_celllines < 0: selected_celllines = [test_line] + test_celllines #All cell lines
        else: selected_celllines = [test_line] + [random.choice(test_celllines) for i in range(num_celllines-1)] #Sample cell lines with replacement
        print selected_celllines

        # Read in and preprocess count data
        data, meta, cell_lines, genes, gene_index, tl_reps = subsample_and_preprocess(inputfile, [('CTRL',-1)] + [(x, num_replicates) for x in selected_celllines])
        ctrldata = data[:,cell_lines.index('CTRL'),:]
        testdata = data[:,[cell_lines.index(x) for x in selected_celllines],:]

        #Run JACKS
        if bs == 0: 
            gene_ws = {x: {gene:[] for gene in gene_index} for x in [test_line] + test_celllines}
        bs_gene_ws = {x: {gene:[] for gene in gene_index} for x in [test_line] + test_celllines}
        for gene in gene_index:
        
            Ig = gene_index[gene]
            # If restricting the number of guides per gene...
            if num_guides > 0 and len(Ig) >= num_guides:
                Ig = random.sample(Ig, num_guides)

            # Perform inference
            y, tau, x1, x2, w1, w2 = infer_JACKS_gene(testdata[Ig,:,0], testdata[Ig,:,1], ctrldata[Ig,0], ctrldata[Ig,1], 50)
                
            for i, cellline in enumerate(selected_celllines):
                bs_gene_ws[cellline][gene].append((w1[i],w2[i]))

        for gene in gene_index:
            for i, cellline in enumerate(selected_celllines):
                for w1, w2 in bs_gene_ws[cellline][gene]:
                    gene_ws[cellline][gene].append((w1,w2))

    if not os.path.exists(os.path.dirname(outfile)): os.makedirs(os.path.dirname(outfile))     
    fout = io.open(outfile[:-4] + '_%s_' % job_idx + test_line + outfile[-4:], 'w')
    num_bs_for_cl = len(gene_ws[test_line][gene_index.keys()[0]])
    fout.write(u'Gene\t%s\n' % ('\t'.join(['WMean\tWStd' for i in range(num_bs_for_cl)])))
    for gene in gene_ws[test_line]:
        fout.write(u'%s\t%s\n' % (gene, '\t'.join(['%5e\t%5e' % (w1, w2) for (w1,w2) in gene_ws[test_line][gene]])))
    fout.close()
