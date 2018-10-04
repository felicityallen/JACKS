import os, io, random, csv, argparse
import numpy as np
from jacks.infer import LOG, inferJACKS
from jacks.preprocess import loadDataAndPreprocess, collateTestControlSamples

REP_HDR_DEFAULT = "Replicate"
SAMPLE_HDR_DEFAULT = "Sample"
COMMON_CTRL_SAMPLE_DEFAULT = "CONTROL"
SGRNA_HDR_DEFAULT = "sgRNA"
GENE_HDR_DEFAULT = "Gene"
OUTPREFIX_DEFAULT = ""
APPLY_W_HP_DEFAULT = False
PICKLE_FILENAME = '_JACKS_results_full.pickle'

def getGeneWs(jacks_results, gene):
    return jacks_results[gene][4]

""" Sort genes in order of mean w1 across cell lines
@param jacks_results: output of infer_JACKS
@return list of genes sorted with largest average effects first
"""
def getSortedGenes(jacks_results):
    #Sort genes by w1
    ordered_genes = [(np.nanmean(getGeneWs(jacks_results, gene)), gene) for gene in jacks_results]
    ordered_genes.sort()
    return ordered_genes

def writeJacksWResults( filename, jacks_results, cell_lines, write_w2=False ):
    #Sort genes by w1
    ordered_genes = getSortedGenes(jacks_results)
    fout = io.open(filename, 'w')
    fout.write(u'Gene\t%s\n' % ('\t'.join(cell_lines)))
    for w1_mean, gene in ordered_genes:
        if write_w2:
            w1s = ['%5e' % np.sqrt(w2 - w1**2.0) for (w1,w2) in zip(jacks_results[gene][4],jacks_results[gene][5])]
        else:
            w1s = ['%5e' % w1 for w1 in jacks_results[gene][4]]
        w1_str = '\t'.join(w1s)
        fout.write(u'%s\t%s\n' % (gene, w1_str))
    fout.close()

def writeJacksXResults( filename, jacks_results, gene_grnas ):
    #Sort genes by mean w1
    ordered_genes = [(np.nanmean(jacks_results[gene][4]),gene) for gene in jacks_results]
    ordered_genes.sort()
    
    fout = io.open(filename, 'w')
    fout.write(u'sgrna\tX1\tX2\n')
    for w1_mean, gene in ordered_genes:
        y, tau, x1, x2, w1, w2 = jacks_results[gene]
        for i, grna in enumerate(gene_grnas[gene]):
            fout.write(u'%s\t%5e\t%5e\n' % (grna, x1[i], x2[i]))
    fout.close()   

def writeFoldChanges(filename, data, meta, sample_ids, write_std=False):
    fout = io.open(filename, 'w')
    fout.write(u'gRNA\tgene\t%s\n' % '\t'.join(sample_ids))
    for i in range(len(meta[:,0])):
        fout.write(u'%s\t%s\t%s\n' % (meta[i,0],meta[i,1],'\t'.join(['%6e' % x for x in data[i,:,int(write_std)]])))
    fout.close()
    
def pickleJacksFullResults( filename, jacks_results, cell_lines, gene_grnas ):    
    import pickle
    full_results = [jacks_results, cell_lines, gene_grnas]
    f = io.open(filename,'wb')
    pickle.dump(full_results, f)
    f.close()
    
def loadJacksFullResultsFromPickle( filename ):
    import pickle
    f = io.open(filename,'rb')
    jacks_results, cell_lines, gene_grnas = pickle.load(f)
    f.close()
    return jacks_results, cell_lines, gene_grnas
    
def prepareFile(filename, hdr):
    # Count any lines before the headers (should be skipped)
    f = io.open(filename)
    skip_lines, line = 0, f.readline()
    while hdr not in line and skip_lines < 100: skip_lines += 1; line = f.readline()
    f.close()

    if skip_lines >= 100:
        raise Exception('Could not find line with header ' + hdr + ' in ' + filename)

    # Check for comma vs tab delimited
    delim = ',' if (filename.split('.')[-1] == 'csv') else '\t'
    
    #Reopen the file and skip to the start of the data
    f = io.open(filename); [f.readline() for i in range(skip_lines)]  
    return f, delim

#output:  {input_filename:[(sample_id, colname)]}
def createSampleSpec(infile, repfile, rep_hdr, sample_hdr, common_ctrl_sample, ctrl_sample_hdr=None):
    f, delim = prepareFile(repfile, rep_hdr)
    rdr = csv.DictReader(f, delimiter=delim)
    sample_spec, sample_num_reps = {infile:[]}, {}
    ctrl_per_sample = not (ctrl_sample_hdr is None)
    if ctrl_per_sample and ctrl_sample_hdr not in rdr.fieldnames:
        raise Exception('Could not find column %s specifying controls in %s' % (ctrl_sample_hdr, repfile))
    ctrl_spec = {}
    for row in rdr:
        sample_id, rep_id = row[sample_hdr], row[rep_hdr]
        if sample_id not in sample_num_reps:
            sample_num_reps[sample_id] = 0
        sample_num_reps[sample_id] += 1
        sample_spec[infile].append((sample_id, rep_id)) 
        if ctrl_per_sample: 
            if row[sample_hdr] in ctrl_spec:
                if ctrl_spec[row[sample_hdr]] != row[ctrl_sample_hdr]:
                    err_msg = '%s vs %s for %s\n' % (ctrl_spec[row[sample_hdr]], row[ctrl_sample_or_hdr], row[sample_hdr])
                    raise Exception(err_msg + 'Different controls for replicates of the sample not supported.')
            else: ctrl_spec[row[sample_hdr]] = row[ctrl_sample_hdr]
        else:
            ctrl_spec[row[sample_hdr]] = common_ctrl_sample
    f.close()
    if not ctrl_per_sample and common_ctrl_sample not in ctrl_spec:
        raise Exception('Could not find control sample %s in %s' % (common_ctrl_sample, repfile))
    return sample_spec, ctrl_spec, sample_num_reps
    
#output:  {grna: gene}
def createGeneSpec(guidemappingfile, sgrna_hdr, gene_hdr):
    f, delim = prepareFile(guidemappingfile, sgrna_hdr)
    rdr = csv.DictReader(f, delimiter=delim)
    gene_spec = {row[sgrna_hdr]: row[gene_hdr] for row in rdr}
    f.close()
    return gene_spec

#input file with set of guide identifiers (one per line)
def readGuideset(filename):
    f = io.open(filename)
    guideset = set([line[:-1] for line in f])
    f.close()
    return guideset

def loadSgrnaReference(filename):
    f = io.open(filename)
    x_ref = {row['sgrna']: row for row in csv.DictReader(f, delimiter='\t')}
    f.close()
    return x_ref

def getJacksParser():
    """
    Create command-line parser for run_jacks program
    """
    run_jacks_doc = ""
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                 description=run_jacks_doc)
    ap.add_argument("countfile",
                    help="Countfile: assumes sgrna label is always in the first column!")
    replicatefile_group = ap.add_argument_group("Replicate file arguments")
    replicatefile_group.add_argument("replicatefile",
                                     help="A CSV or tab delimited file mapping replicates to samples")
    replicatefile_group.add_argument("--rep_hdr",
                                     type=str,
                                     default=REP_HDR_DEFAULT,
                                     help="Column header for column containing the replicate labels")
    replicatefile_group.add_argument("--sample_hdr",
                                     type=str,
                                     default=SAMPLE_HDR_DEFAULT,
                                     help="Column header for column containing the sample labels")
    replicatefile_group.add_argument("--common_ctrl_sample",
                                     type=str,
                                     default=COMMON_CTRL_SAMPLE_DEFAULT,
                                     help="Name of the sample to be used as a control for all other samples")
    replicatefile_group.add_argument("--ctrl_sample_hdr",
                                     type=str,
                                     default=None,
                                     help="Column header for column containing the sample label for the control to be used for each sample")
    guidemapping_group = ap.add_argument_group("Guidemapping file arguments")
    guidemapping_group.add_argument("guidemappingfile",
                                    help="A CSV or tab delimited file mapping guides to genes")
    guidemapping_group.add_argument("--sgrna_hdr",
                                    type=str,
                                    default=SGRNA_HDR_DEFAULT,
                                    help="Column header for the column containing the guide labels")
    guidemapping_group.add_argument("--gene_hdr",
                                    type=str,
                                    default=GENE_HDR_DEFAULT,
                                    help="Column header for the column containing the gene labels")

    ap.add_argument("--outprefix",
                    type=str,
                    default="",
                    help="Output prefix")
    ap.add_argument("--reffile",
                    type=str,
                    default=None,
                    help="Reference file containing pre-trained sgRNA efficacy scores")
    ap.add_argument("--apply_w_hp",
                    action='store_true',
                    default=False,
                    help="Apply hierarchical prior on gene effects (use with caution)")
    return ap

def runJACKS(countfile, replicatefile, guidemappingfile,
              rep_hdr=REP_HDR_DEFAULT, sample_hdr=SAMPLE_HDR_DEFAULT, common_ctrl_sample=COMMON_CTRL_SAMPLE_DEFAULT,
              ctrl_sample_hdr=None, sgrna_hdr=SGRNA_HDR_DEFAULT, gene_hdr=GENE_HDR_DEFAULT, 
              outprefix=OUTPREFIX_DEFAULT, reffile=None, apply_w_hp=APPLY_W_HP_DEFAULT):
    outprefix = outprefix
    if '/' in outprefix and not os.path.exists(os.path.dirname(outprefix)): os.makedirs(os.path.dirname(outprefix))
    outfile_w = outprefix + '_gene_JACKS_results.txt'
    outfile_w2 = outprefix + '_genestd_JACKS_results.txt'
    outfile_x = outprefix + '_grna_JACKS_results.txt'
    outfile_lfc = outprefix + '_logfoldchange_means.txt'
    outfile_lfc_std = outprefix + '_logfoldchange_std.txt'
    outfile_pickle = outprefix + PICKLE_FILENAME

    # Load the specification of samples to include
    LOG.info('Loading sample specification')
    sample_spec, ctrl_spec, sample_num_reps = createSampleSpec(countfile, replicatefile, rep_hdr,
                                                               sample_hdr, common_ctrl_sample, ctrl_sample_hdr)
    # Load the mappings from guides to genes
    LOG.info('Loading gene mappings')
    gene_spec = createGeneSpec(guidemappingfile, sgrna_hdr, gene_hdr)

    sgrna_reference_file = reffile
    if sgrna_reference_file:
        # Load the sgrna reference (precomputed X's)
        LOG.info('Loading sgrna reference values')
        x_ref = loadSgrnaReference(reffile)
        # Check that the data to be loaded have sgrna reference values
        LOG.info('Checking sgrna reference identifiers against gene mappings')
        for guide in gene_spec:
            if guide not in x_ref:
                raise Exception('%s has no sgrna reference in %s' % (guide, sgrna_reference_file))

    # Load the data and preprocess
    LOG.info('Loading data and pre-processing')
    data, meta, sample_ids, genes, gene_index = loadDataAndPreprocess(sample_spec, gene_spec)
    gene_grnas = {gene: [x for x in meta[gene_index[gene], 0]] for gene in gene_index}
    x_reference = None
    if sgrna_reference_file:
        # Create the X reference (in the correct order)
        x_reference = {'X1': np.array([eval(x_ref[x]['X1']) for x in meta[:, 0]]),
                       'X2': np.array([eval(x_ref[x]['X2']) for x in meta[:, 0]])}
    else:
        writeFoldChanges(outfile_lfc, data, meta, sample_ids)
        writeFoldChanges(outfile_lfc_std, data, meta, sample_ids, write_std=True)
        
    #Run all samples against their controls
    LOG.info('Running JACKS inference')
    testdata, ctrldata, test_sample_idxs = collateTestControlSamples(data, sample_ids, ctrl_spec)
    jacks_results = inferJACKS(gene_index, testdata, ctrldata, apply_w_hp=apply_w_hp)

    # Write out the results
    LOG.info('Writing JACKS results')
    sample_ids_without_ctrl = [sample_ids[idx] for idx in test_sample_idxs]
    writeJacksWResults(outfile_w, jacks_results, sample_ids_without_ctrl)
    writeJacksWResults(outfile_w2, jacks_results, sample_ids_without_ctrl, write_w2=True)
    writeJacksXResults(outfile_x, jacks_results, gene_grnas)
    pickleJacksFullResults(outfile_pickle, jacks_results, sample_ids_without_ctrl, gene_grnas)

def runJACKSFromArgs():
    parser = getJacksParser()
    args = parser.parse_args()
    runJACKS(args.countfile, args.replicatefile, args.guidemappingfile,
              args.rep_hdr, args.sample_hdr, args.common_ctrl_sample, 
              args.ctrl_sample_hdr, args.sgrna_hdr, args.gene_hdr, 
              args.outprefix, args.reffile, args.apply_w_hp)
