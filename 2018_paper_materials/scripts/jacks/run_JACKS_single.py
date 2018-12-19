import io, os, sys, csv, random, logging
from jacks.infer import LOG, inferJACKS
from jacks.jacks_io import createPseudoNonessGenes, readControlGeneset, createGeneSpec, createSampleSpec, getJacksParser, collateTestControlSamples, writeJacksWResults
from jacks.preprocess import loadDataAndPreprocess

py_cmd = 'python'

def combineSingleResults(single_jacks_results):
    jacks_results = {}
    for gene in single_jacks_results[0].keys():
        y, tau, x1, x2, w1, w2 = single_jacks_results[0][gene]
        w1 = [x[gene][4] for x in single_jacks_results]
        w2 = [x[gene][5] for x in single_jacks_results]
        jacks_results[gene] = y, tau, x1, x2, w1, w2
    return jacks_results

def filterSampleSpec(sample_spec, cell_line, ctrl_spec):
    new_sample_spec = {}
    for filename in sample_spec:
        for sample_id, colname in sample_spec[filename]:
            if sample_id == cell_line or sample_id == ctrl_spec[cell_line]:
                if filename not in new_sample_spec:
                    new_sample_spec[filename] = []
                new_sample_spec[filename].append((sample_id, colname))
    return new_sample_spec
    
def filterCtrlSpec(ctrl_spec, cell_line):
    new_ctrl_spec = {}
    new_ctrl_spec[cell_line] = ctrl_spec[cell_line]             #Sample
    new_ctrl_spec[ctrl_spec[cell_line]] = ctrl_spec[cell_line]  #Control
    return new_ctrl_spec 

if __name__ == '__main__':

    LOG.setLevel(logging.WARNING)
    parser = getJacksParser()
    parser.add_argument("--cell_line",
                    type=str,
                    default=None,
                    help="cell line to run")
    parser.add_argument("--separate",
                    action='store_true',
                    default=False,
                    help="Run cell lines separately")
    args = parser.parse_args()

    outprefix = args.outprefix
    if '/' in outprefix and not os.path.exists(os.path.dirname(outprefix)): os.makedirs(os.path.dirname(outprefix))

    # Load the specification of samples to include
    LOG.info('Loading sample specification')
    sample_spec, ctrl_spec, sample_num_reps = createSampleSpec(args.countfile, args.replicatefile, args.rep_hdr,
                                                                args.sample_hdr, args.common_ctrl_sample, args.ctrl_sample_hdr)
    if args.cell_line != None:
        sample_spec = filterSampleSpec(sample_spec, args.cell_line, ctrl_spec)
        ctrl_spec = filterCtrlSpec(ctrl_spec, args.cell_line)
        outprefix += ('_' + args.cell_line)

    elif args.separate:
        for cell_line in ctrl_spec:
            if ctrl_spec[cell_line] == cell_line: continue
            cmd = '%s --cell_line=%s' % (' '.join(sys.argv), cell_line)
            os.system('%s %s' % (py_cmd cmd))
        exit()

    # Load the mappings from guides to genes
    LOG.info('Loading gene mappings')
    gene_spec = createGeneSpec(args.guidemappingfile, args.sgrna_hdr, args.gene_hdr, ignore_blank_genes=args.ignore_blank_genes)

    # Load negative control guides (if any)
    ctrl_geneset = readControlGeneset(args.ctrl_genes) if args.ctrl_genes is not None else set()

    # Load the data and preprocess
    LOG.info('Loading data and pre-processing')
    data, meta, sample_ids, genes, gene_index = loadDataAndPreprocess(sample_spec, gene_spec,ctrl_spec=ctrl_spec, normtype=args.norm_type, ctrl_geneset=ctrl_geneset)
    gene_grnas = {gene: [x for x in meta[gene_index[gene], 0]] for gene in gene_index}
    testdata, ctrldata, test_sample_idxs = collateTestControlSamples(data, sample_ids, ctrl_spec)
    sample_ids_without_ctrl = [sample_ids[idx] for idx in test_sample_idxs]

    #Run all samples against their controls
    LOG.info('Running Single JACKS inference')
    single_jacks_results = []
    for ts in range(testdata.shape[1]):
        single_jacks_results.append(inferJACKS(gene_index, testdata[:,[ts],:], ctrldata[:,[ts],:], w_only=True))
    jacks_results = combineSingleResults(single_jacks_results)


    #Add a set of pseudo genes, created by randomly sampling from guides targeting genes in the control set
    if args.n_pseudo > 0 and len(ctrl_geneset) > 0:
        LOG.info('Running Single JACKS inference on %d pseudogenes' % args.n_pseudo)
        pseudo_gene_index = createPseudoNonessGenes(gene_index, ctrl_geneset, args.n_pseudo)
        pseudo_single_results = []
        for ts in range(testdata.shape[1]):
            pseudo_single_results.append(inferJACKS(pseudo_gene_index, testdata[:,[ts],:], ctrldata[:,[ts],:], w_only=True))
        jacks_pseudo_results = combineSingleResults(pseudo_single_results)

        writeJacksWResults(outprefix + '_pseudo_noness', jacks_pseudo_results, sample_ids_without_ctrl, write_types=['', '_std'] )

    # Write out the results
    LOG.info('Writing Single JACKS results')
    if len(ctrl_geneset) > 1:
        writeJacksWResults(outprefix, jacks_results, sample_ids_without_ctrl, ctrl_geneset=ctrl_geneset, write_types=['', '_std', '_pval', '_fdr'], fdr=args.fdr, fdr_thresh_type=args.fdr_thresh_type)
    else:
        writeJacksWResults(outprefix, jacks_results, sample_ids_without_ctrl, ctrl_geneset=ctrl_geneset, write_types=['', '_std'])


    #Write pseudo-normalized pvalue results
    if args.n_pseudo > 0 and len(ctrl_geneset) > 0:
        LOG.info('Writing pseudo-normalized Single JACKS results')
        for gene in jacks_pseudo_results:  jacks_results[gene] = jacks_pseudo_results[gene]
        pseudo_genes = set([gene for gene in jacks_pseudo_results])
        writeJacksWResults(outprefix + '_pseudo_combined', jacks_results, sample_ids_without_ctrl, ctrl_geneset=pseudo_genes, write_types=['', '_std', '_pval', '_fdr'], fdr=args.fdr, fdr_thresh_type=args.fdr_thresh_type)




