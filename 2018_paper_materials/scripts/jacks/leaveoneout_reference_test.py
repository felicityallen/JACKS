import io, sys, os, csv, random, logging
from jacks.infer import LOG, inferJACKS
from jacks.jacks_io import loadSgrnaReference, readControlGeneset, createGeneSpec, createSampleSpec, getJacksParser, collateTestControlSamples, writeJacksWResults, writeJacksXResults
from jacks.preprocess import loadDataAndPreprocess
from run_JACKS_single import filterSampleSpec, filterCtrlSpec
import numpy as np

py_cmd = 'python'

def filterOutSampleSpec(sample_spec, cell_line, ctrl_spec):
    new_sample_spec = {}
    for filename in sample_spec:
        for sample_id, colname in sample_spec[filename]:
            if sample_id == cell_line: continue
            if filename not in new_sample_spec:
                new_sample_spec[filename] = []
            new_sample_spec[filename].append((sample_id, colname))
    return new_sample_spec
    
LOG.setLevel(logging.WARNING)
parser = getJacksParser()
parser.add_argument("--cell_line",
                type=str,
                default=None,
                help="cell line to run")
parser.add_argument("--separate",
                action='store_true',
                default=False,
                help="Run all cell lines separately")
args = parser.parse_args()

outprefix = args.outprefix
if '/' in outprefix and not os.path.exists(os.path.dirname(outprefix)): os.makedirs(os.path.dirname(outprefix))

# Load the specification of samples to include
LOG.info('Loading sample specification')
sample_spec, ctrl_spec, sample_num_reps = createSampleSpec(args.countfile, args.replicatefile, args.rep_hdr,
                                                            args.sample_hdr, args.common_ctrl_sample, args.ctrl_sample_hdr)
if args.cell_line != None:
    single_sample_spec = filterSampleSpec(sample_spec, args.cell_line, ctrl_spec)  
    single_ctrl_spec = filterCtrlSpec(ctrl_spec, args.cell_line)
    ref_sample_spec = filterOutSampleSpec(sample_spec, args.cell_line, ctrl_spec)
    ref_outfile =  outprefix + '_xreference_' + args.cell_line + '.txt'

elif args.separate:
    for cell_line in ctrl_spec:
        if ctrl_spec[cell_line] == cell_line: continue
        cmd = '%s --cell_line=%s' % (' '.join(sys.argv), cell_line)
        os.system('%s %s' % (py_cmd,cmd))
    exit()

# Load the mappings from guides to genes
LOG.info('Loading gene mappings')
gene_spec = createGeneSpec(args.guidemappingfile, args.sgrna_hdr, args.gene_hdr, ignore_blank_genes=args.ignore_blank_genes)

# Load negative control guides (if any)
ctrl_geneset = readControlGeneset(args.ctrl_genes) if args.ctrl_genes is not None else set()

##REFERENCE (to collect X's)

# Load the data and preprocess
LOG.info('Reference: Loading data and pre-processing')
data, meta, sample_ids, genes, gene_index = loadDataAndPreprocess(ref_sample_spec, gene_spec,ctrl_spec=ctrl_spec, normtype=args.norm_type, ctrl_geneset=ctrl_geneset)
gene_grnas = {gene: [x for x in meta[gene_index[gene], 0]] for gene in gene_index}
testdata, ctrldata, test_sample_idxs = collateTestControlSamples(data, sample_ids, ctrl_spec)
sample_ids_without_ctrl = [sample_ids[idx] for idx in test_sample_idxs]

#Run all samples against their controls
LOG.info('Reference: Running JACKS inference')
jacks_results = inferJACKS(gene_index, testdata, ctrldata)
writeJacksXResults(ref_outfile, jacks_results, gene_grnas)

##TEST (using reference)
LOG.info('Test: Loading data and pre-processing')
data, meta, sample_ids, genes, gene_index = loadDataAndPreprocess(single_sample_spec, gene_spec,ctrl_spec=single_ctrl_spec, normtype=args.norm_type, ctrl_geneset=ctrl_geneset)
gene_grnas = {gene: [x for x in meta[gene_index[gene], 0]] for gene in gene_index}
testdata, ctrldata, test_sample_idxs = collateTestControlSamples(data, sample_ids, ctrl_spec)
sample_ids_without_ctrl = [sample_ids[idx] for idx in test_sample_idxs]

x_ref = loadSgrnaReference(ref_outfile)
# Create the X reference (in the correct order)
x_reference = {'X1': np.array([eval(x_ref[x]['X1']) for x in meta[:, 0]]),
                'X2': np.array([eval(x_ref[x]['X2']) for x in meta[:, 0]])}

#Run all samples against their controls
LOG.info('Test: Running JACKS inference using reference')
test_jacks_results = inferJACKS(gene_index, testdata, ctrldata, w_only=True, fixed_x=x_reference)

# Write out the results
LOG.info('Writing JACKS results')
writeJacksWResults(outprefix + '_' + args.cell_line, test_jacks_results, sample_ids_without_ctrl, ctrl_geneset=ctrl_geneset, write_types=[''])
