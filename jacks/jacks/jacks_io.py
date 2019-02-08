import os, io, random, csv, argparse
import numpy as np
from jacks.infer import LOG, inferJACKS
from jacks.preprocess import loadDataAndPreprocess, collateTestControlSamples
import scipy.stats as ST

REP_HDR_DEFAULT = "Replicate"
SAMPLE_HDR_DEFAULT = "Sample"
COMMON_CTRL_SAMPLE_DEFAULT = "CONTROL"
SGRNA_HDR_DEFAULT = "sgRNA"
GENE_HDR_DEFAULT = "Gene"
OUTPREFIX_DEFAULT = ""
APPLY_W_HP_DEFAULT = False
PICKLE_FILENAME = '_JACKS_results_full.pickle'
NORM_TYPE_DEFAULT = 'median'

def getGeneWs(jacks_results, gene):
    return jacks_results[gene][4]

""" Sort genes in order of mean w1 across cell lines
@param jacks_results: output of infer_JACKS
@return list of genes sorted with largest average effects first
"""
def getSortedGenes(jacks_results):
    #Sort genes by w1
    ordered_genes = [(np.nanmean(jacks_results[gene][4]),gene) for gene in jacks_results]
    ordered_genes.sort()
    return ordered_genes

def getMode(vals, allow_bimode=True):
    if len(vals.shape)==2 and vals.shape[1] ==1: vals = vals[:,0]
    kernel, mean_vals, std_vals = ST.gaussian_kde(vals), np.mean(vals), np.nanstd(vals)
    xdata = np.arange(mean_vals-std_vals,mean_vals+std_vals,0.01)
    kdata = kernel(xdata)
    kdata_grad = kdata[1:] - kdata[:-1]
    pos_grads = (kdata_grad>=0)
    mode_idxs = np.where(pos_grads[1:].astype(int) - pos_grads[:-1].astype(int) == -1)[0]
    max_mode_height, max_mode_idx = max([(kdata[idx+1],idx) for idx in mode_idxs])
    mode_idxs = [idx for idx in mode_idxs if kdata[idx+1]/max_mode_height > 0.2]
    modes = [xdata[idx+1] for idx in mode_idxs]
    min_idxs = np.where(pos_grads[1:].astype(int) - pos_grads[:-1].astype(int) == 1)[0]
    mins = [xdata[idx+1] for idx in min_idxs if idx > min(mode_idxs) and idx < max(mode_idxs)]
    if len(modes) == 1: #If one clear mode, use that
        return modes[0]
    elif allow_bimode and len(modes) == 2 and len(mins)==1:  #Sometimes it becomes bimodal around 0, so take the trough between the two modes
        return mins[0]
    else: return xdata[max_mode_idx+1] #Otherwise just take the heighest one

def getFDRGeneSets(jacks_w1_pvals, fdr):
    fdr_gene_sets = []
    a_gene = [x for x in jacks_w1_pvals][0]
    for sample_idx in range(len(jacks_w1_pvals[a_gene])):
        pval_genes = [(jacks_w1_pvals[gene][sample_idx], gene) for gene in jacks_w1_pvals]
        pval_genes.sort()
        below_idxs = [i for i,(pval, gene) in enumerate(pval_genes) if pval < fdr*(i+1)/len(pval_genes)]
        k = max(below_idxs) if len(below_idxs) > 0 else -1
        fdr_gene_sets.append(set([gene for (pval, gene) in pval_genes[:k+1]]))
    return fdr_gene_sets

def getLocalFDRGeneSets(jacks_w1_fdrs,test_sample_ids, fdr):
    fdr_gene_sets = []
    for sample_idx in enumerate(test_sample_ids):
        fdr_gene_sets.append(set([gene for gene in jacks_w1_pvals if jacks_w1_pvals[gene][sample_idx] < fdr ]))
    return fdr_gene_sets

def computeW1PvalsAndFDRs(jacks_results, test_sample_ids, noness_genes = set(), pseudo=False, compute_fdr=False):

    is_noness = lambda x: x in noness_genes
    w1_pvals = {gene: np.zeros(len(test_sample_ids)) for gene in jacks_results}
    w1_fdrs = {gene: np.zeros(len(test_sample_ids)) for gene in jacks_results}
    for idx, sample_id in enumerate(test_sample_ids):
        if compute_fdr:
            w1_all = np.array([jacks_results[gene][4][idx].item() for gene in jacks_results if not np.isnan(jacks_results[gene][4][idx])])
            k_all = ST.gaussian_kde(w1_all)
        noness_w1 = np.array([jacks_results[gene][4][idx].item() for gene in jacks_results if is_noness(gene) and not np.isnan(jacks_results[gene][4][idx])])
        k_noness = ST.gaussian_kde(noness_w1)
        for gene in jacks_results:
            if is_noness(gene) and not pseudo: 	#For non-essential genes, build another kernel that doesn't include the gene itself
                kernel = ST.gaussian_kde(np.array([jacks_results[x][4][idx].item() for x in jacks_results if (is_noness(x) and (x != gene))]))
            else: kernel = k_noness	#Otherwise use all the non-essential (or pseudo) genes
            w1_pvals[gene][idx] = kernel.integrate_box_1d(-100.0, jacks_results[gene][4][idx])
            if compute_fdr: w1_fdrs[gene][idx] = kernel.evaluate(jacks_results[gene][4][idx])/k_all.evaluate(jacks_results[gene][4][idx])
    return w1_pvals, w1_fdrs

def writeJacksWResults( outprefix, jacks_results, cell_lines, write_types=[''], ctrl_geneset=set(), fdr=None, fdr_thresh_type='REGULAR', pseudo=False):
    #Sort genes by w1
    ordered_genes = getSortedGenes(jacks_results)
    fouts = [io.open(outprefix + '_gene%s_JACKS_results.txt' % write_type,'w') for write_type in write_types]
    for fout in fouts: fout.write(u'Gene\t%s\n' % ('\t'.join(cell_lines)))
    if '_fdr' in write_types or '_pval' in write_types:
        LOG.info('Computing P-values')
        jacks_w1_pvals,jacks_w1_fdrs =  computeW1PvalsAndFDRs(jacks_results, cell_lines, noness_genes = ctrl_geneset, pseudo=pseudo, compute_fdr=('_fdr' in write_types))
    
    #Determine threshold sets for fdr cut-offs (blank out non-significant genes)
    if fdr is not None:
        if fdr_thresh_type == 'REGULAR':
            fdr_sets = getFDRGeneSets(jacks_w1_pvals, fdr)
        elif fdr_thresh_type == 'LOCAL_FDR':
            fdr_sets = getLocalFDRGeneSets(jacks_w1_fdrs, fdr)
        else: raise Exception('Unrecognised FDR threshold type (expecting REGULAR or LOCAL_FDR): ', fdr_thresh_type)

    #Write out one line per gene (all cell lines)
    for w1_mean, gene in ordered_genes:
        for write_type,fout in zip(write_types, fouts):
            
            #Determine whether to include the gene for each cell line (if fdr thresholded)
            if fdr is not None:
                sig_gene_flags = [(gene in x) for x in fdr_sets]
            else: 
                sig_gene_flags = [True for x in jacks_results[gene][4]]
            if sum(sig_gene_flags) == 0: continue

            #Write out the values
            if write_type=='_pval':
                w1s = ['%5e' % x for x in jacks_w1_pvals[gene]]
            elif write_type == '_fdr':
                w1s = ['%5e' % x for x in jacks_w1_fdrs[gene]]
            elif write_type == '_std':
                w1s = ['%5e' % np.sqrt(w2 - w1**2.0) for (w1,w2) in zip(jacks_results[gene][4],jacks_results[gene][5])]
            elif write_type == '':
                w1s = [('%5e' % w1) if flag else '' for (w1,flag) in zip(jacks_results[gene][4],sig_gene_flags)]
            else:  raise Exception('Unrecognised write type: %s' % write_type)
            w1_str = '\t'.join(w1s)
            if 'JACKS_PSEUDO_GENE' not in gene:
                fout.write(u'%s\t%s\n' % (gene, w1_str))
    for fout in fouts: fout.close()

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

def writeFoldChanges(filename, testdata, ctrldata, meta, sample_ids, write_std=False, write_raw = False):
    fout = io.open(filename, 'w')
    fout.write(u'gRNA\tgene\t%s\n' % '\t'.join(sample_ids))
    for i in range(len(meta[:,0])):
       if write_std:  vals = [np.sqrt(x**2.0 + y**2.0) for (x,y) in zip(testdata[i,:,1],ctrldata[i,:,1])]
       elif write_raw:  vals = [np.sqrt(x**2.0 + y**2.0) for (x,y) in zip(testdata[i,:,2],ctrldata[i,:,2])]
       else: vals = testdata[i,:,0] - ctrldata[i,:,0]
       fout.write(u'%s\t%s\t%s\n' % (meta[i,0],meta[i,1],'\t'.join(['%6e' % x for x in vals])))
    fout.close()

def writeLogData(filename, data, meta, sample_ids, write_std=False, write_raw = False):
    fout = io.open(filename, 'w')
    fout.write(u'gRNA\tgene\t%s\n' % '\t'.join(sample_ids))
    for i in range(len(meta[:,0])):
       if write_std:  vals = data[i,:,1]
       elif write_raw:  vals = data[i,:,2]
       else: vals = data[i,:,0]
       fout.write(u'%s\t%s\t%s\n' % (meta[i,0],meta[i,1],'\t'.join(['%6e' % x for x in vals])))
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


# output:  {input_filename:[(sample_id, colname)]}
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
def createGeneSpec(guidemappingfile, sgrna_hdr, gene_hdr, ignore_blank_genes=True):
    f, delim = prepareFile(guidemappingfile, sgrna_hdr)
    rdr = csv.DictReader(f, delimiter=delim)
    gene_spec = {row[sgrna_hdr]: row[gene_hdr] for row in rdr if (not ignore_blank_genes or row[gene_hdr] != '')}
    f.close()
    return gene_spec

#Either gene name or input file with set of gene identifiers (one per line)
def readControlGeneset(ctrl_genes, gene_spec):
    known_genes = set([gene_spec[x] for x in gene_spec])
    if os.path.isfile(ctrl_genes):
        f = io.open(ctrl_genes)
        geneset = set([line.split()[0] for line in f if line.split()[0] in known_genes])
        f.close()
        LOG.info('Read %d recognised control genes from %s' % (len(geneset), ctrl_genes))
    else: 
        if ctrl_genes not in known_genes: raise Exception('Not a file or unrecognised control gene: %s' % ctrl_genes) 
        geneset = set([ctrl_genes])
        LOG.info('Using %s as control gene' % (ctrl_genes))
    return geneset

def loadSgrnaReference(filename):
    f = io.open(filename)
    x_ref = {row['sgrna']: row for row in csv.DictReader(f, delimiter='\t')}
    f.close()
    return x_ref

def createPseudoNonessGenes(gene_index, ctrl_geneset, num_psuedo_noness_genes):
    pseudo_gene_index = {}
    ctrl_geneset_lst = [x for x in ctrl_geneset if x in gene_index]
    all_geneset_lst = [x for x in gene_index]
    for i in range(num_psuedo_noness_genes):
        num_guides = len(gene_index[random.choice(all_geneset_lst)])
        pseudo_idxs = [random.choice(gene_index[random.choice(ctrl_geneset_lst)]) for j in range(num_guides)]
        pseudo_gene_index['JACKS_PSEUDO_GENE_%d' % i] = pseudo_idxs
    return pseudo_gene_index

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
    guidemapping_group.add_argument("--ignore_blank_genes",
                                    action='store_true',
                                    default=False,
                                    help="Ignore guides for which the assigned gene is blank")

    ap.add_argument("--outprefix",
                    type=str,
                    default="",
                    help="Output prefix")
    ap.add_argument("--reffile",
                    type=str,
                    default=None,
                    help="Reference file containing pre-trained sgRNA efficacy scores")
    ap.add_argument("--fdr",
                    type=float,
                    default=None,
                    help="False discovery rate threshold at which to accept significant genes (if None output all in order)")
    ap.add_argument("--fdr_thresh_type",
                    type=str,
                    default='REGULAR',
                    help="Method to use to threshold significant genes (REGULAR, LOCAL_FDR)")
    ap.add_argument("--apply_w_hp",
                    action='store_true',
                    default=False,
                    help="Apply hierarchical prior on gene effects (not recommended, use with caution)")
    ap.add_argument("--norm_type",
                    type=str,
                    default=NORM_TYPE_DEFAULT,
                    help="Type of normalisation to use on log count data (median, mode, ctrl_guides - if using ctrl_guides, must provide --ctrl_gene_file)")
    ap.add_argument("--ctrl_genes",
                    type=str,
                    default=None,
                    help="(Required if p-value output is wanted) Either, the name of a gene (as used in sgrnamappingfile) specifying a set of negative control guides (e.g. these could be intergenic, non-targeting etc) OR a text file containing a list (one per line) of genes to use as negative controls. Also used to infer variances in case of single replicate data.")
    ap.add_argument("--n_pseudo",
                    type=int,
                    default=2000,
                    help="Number of pseudo genes to create and infer JACKS results for by randomly sampling guides from the control genes (specified in --ctrl_genes)")
    ap.add_argument("--count_prior",
                    type=int,
                    default=32,
                    help="Prior count to be added to all read counts prior to log operation")
    return ap

def preprocess(countfile, replicatefile, guidemappingfile,
              rep_hdr=REP_HDR_DEFAULT, sample_hdr=SAMPLE_HDR_DEFAULT, common_ctrl_sample=COMMON_CTRL_SAMPLE_DEFAULT,
              ctrl_sample_hdr=None, sgrna_hdr=SGRNA_HDR_DEFAULT, gene_hdr=GENE_HDR_DEFAULT, ignore_blank_genes=False,
              outprefix=OUTPREFIX_DEFAULT, reffile=None):

    # Load the specification of samples to include
    LOG.info('Loading sample specification')
    sample_spec, ctrl_spec, sample_num_reps = createSampleSpec(countfile, replicatefile, rep_hdr,
                                                               sample_hdr, common_ctrl_sample, ctrl_sample_hdr)
    # Load the mappings from guides to genes
    LOG.info('Loading gene mappings')
    gene_spec = createGeneSpec(guidemappingfile, sgrna_hdr, gene_hdr, ignore_blank_genes=ignore_blank_genes)

    sgrna_reference_file = reffile
    x_ref = None
    if sgrna_reference_file:
        # Load the sgrna reference (precomputed X's)
        LOG.info('Loading sgrna reference values')
        x_ref = loadSgrnaReference(reffile)
        # Check that the data to be loaded have sgrna reference values
        LOG.info('Checking sgrna reference identifiers against gene mappings')
        for guide in gene_spec:
            if guide not in x_ref:
                raise Exception('%s has no sgrna reference in %s' % (guide, sgrna_reference_file))
    return sample_spec, ctrl_spec, gene_spec, x_ref


def load_data_and_run(sample_spec, gene_spec, ctrl_spec, sgrna_reference_file, x_ref,
                      outprefix, apply_w_hp=APPLY_W_HP_DEFAULT, norm_type=NORM_TYPE_DEFAULT, 
                      ctrl_genes=None, fdr=None, fdr_thresh_type = 'REGULAR', n_pseudo=0, count_prior=32 ):

    # Load negative control genes (if any)
    ctrl_geneset = readControlGeneset(ctrl_genes, gene_spec) if ctrl_genes is not None else set()

    if '/' in outprefix and not os.path.exists(os.path.dirname(outprefix)): os.makedirs(os.path.dirname(outprefix))
    outfile_x = outprefix + '_grna_JACKS_results.txt'
    outfile_lfc = outprefix + '_logfoldchange_means.txt'
    outfile_lfc_std = outprefix + '_logfoldchange_std.txt'
    outfile_pickle = outprefix + PICKLE_FILENAME

    # Load the data and preprocess
    LOG.info('Loading data and pre-processing')
    data, meta, sample_ids, genes, gene_index = loadDataAndPreprocess(sample_spec, gene_spec,ctrl_spec=ctrl_spec,normtype=norm_type, ctrl_geneset=ctrl_geneset, prior=count_prior)
    gene_grnas = {gene: [x for x in meta[gene_index[gene], 0]] for gene in gene_index}
    testdata, ctrldata, test_sample_idxs = collateTestControlSamples(data, sample_ids, ctrl_spec)
    sample_ids_without_ctrl = [sample_ids[idx] for idx in test_sample_idxs]

    x_reference = None
    if sgrna_reference_file:
        # Create the X reference (in the correct order)
        x_reference = {'X1': np.array([eval(x_ref[x]['X1']) for x in meta[:, 0]]),
                       'X2': np.array([eval(x_ref[x]['X2']) for x in meta[:, 0]])}
    else:
        writeFoldChanges(outfile_lfc, testdata, ctrldata, meta, sample_ids_without_ctrl)
        writeFoldChanges(outfile_lfc_std, testdata, ctrldata, meta, sample_ids_without_ctrl, write_std=True)
        
    #Run all samples against their controls
    LOG.info('Running JACKS inference')
    jacks_results = inferJACKS(gene_index, testdata, ctrldata, apply_w_hp=apply_w_hp, fixed_x=x_reference)

    #Add a set of pseudo genes, created by randomly sampling from guides targeting genes in the control set
    if n_pseudo > 0 and len(ctrl_geneset) > 0:
        LOG.info('Running JACKS inference on %d pseudogenes' % n_pseudo)
        pseudo_gene_index = createPseudoNonessGenes(gene_index, ctrl_geneset, n_pseudo)
        jacks_pseudo_results = inferJACKS(pseudo_gene_index, testdata, ctrldata, apply_w_hp=apply_w_hp)
        writeJacksWResults(outprefix + '_pseudo_noness', jacks_pseudo_results, sample_ids_without_ctrl, write_types=['', '_std'] )
        for gene in jacks_results:
            jacks_pseudo_results[gene] = jacks_results[gene]

    # Write out the results
    LOG.info('Writing JACKS results')
    if len(ctrl_geneset) > 0 and n_pseudo > 0:
        writeJacksWResults(outprefix, jacks_pseudo_results, sample_ids_without_ctrl, ctrl_geneset=set([x for x in jacks_pseudo_results if 'JACKS_PSEUDO_GENE' in x]), write_types=['', '_std', '_pval'], fdr=fdr, pseudo=True, fdr_thresh_type=fdr_thresh_type)
    else:
        writeJacksWResults(outprefix, jacks_results, sample_ids_without_ctrl, ctrl_geneset=ctrl_geneset, write_types=['', '_std'])
    writeJacksXResults(outfile_x, jacks_results, gene_grnas)
    pickleJacksFullResults(outfile_pickle, jacks_results, sample_ids_without_ctrl, gene_grnas)

def runJACKS(countfile, replicatefile, guidemappingfile,
              rep_hdr=REP_HDR_DEFAULT, sample_hdr=SAMPLE_HDR_DEFAULT, common_ctrl_sample=COMMON_CTRL_SAMPLE_DEFAULT,
              ctrl_sample_hdr=None, sgrna_hdr=SGRNA_HDR_DEFAULT, gene_hdr=GENE_HDR_DEFAULT, 
              apply_w_hp=APPLY_W_HP_DEFAULT, norm_type=NORM_TYPE_DEFAULT, 
              ignore_blank_genes=False, ctrl_genes=None, fdr=None, fdr_thresh_type = 'REGULAR',
              outprefix=OUTPREFIX_DEFAULT, reffile=None, n_pseudo=0, count_prior=32):

    sample_spec, ctrl_spec, gene_spec, x_ref = preprocess(countfile, replicatefile, guidemappingfile,
                                                                rep_hdr, sample_hdr, common_ctrl_sample,
                                                                ctrl_sample_hdr, sgrna_hdr, gene_hdr, 
                                                                ignore_blank_genes, outprefix, reffile)
    load_data_and_run(sample_spec, gene_spec, ctrl_spec, reffile, x_ref,
                      outprefix, apply_w_hp=apply_w_hp, norm_type=norm_type, ctrl_genes=ctrl_genes, fdr=fdr, fdr_thresh_type=fdr_thresh_type, n_pseudo=n_pseudo, count_prior=count_prior )
    
def runJACKSFromArgs():
    parser = getJacksParser()
    args = parser.parse_args()
    runJACKS(args.countfile, args.replicatefile, args.guidemappingfile,
              args.rep_hdr, args.sample_hdr, args.common_ctrl_sample, 
              args.ctrl_sample_hdr, args.sgrna_hdr, args.gene_hdr, 
              args.apply_w_hp, args.norm_type, 
              args.ignore_blank_genes, args.ctrl_genes, args.fdr, args.fdr_thresh_type,
              args.outprefix, args.reffile, args.n_pseudo, args.count_prior)

