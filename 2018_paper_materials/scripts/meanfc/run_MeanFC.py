import io, os, csv, random, logging
from jacks.infer import LOG
from jacks.jacks_io import createGeneSpec, createSampleSpec, getJacksParser, collateTestControlSamples, writeJacksWResults
from jacks.preprocess import loadDataAndPreprocess
import scipy as SP

def infer_JACKS_meanfc(gene_index, testdata, ctrldata):
    results = {}
    for gene in gene_index:
        Ig = gene_index[gene]
        y = (testdata[Ig,:,0] - ctrldata[Ig,:,0])
        w1 = SP.nanmean(y,axis=0)
        results[gene] = (y,-1.0,-1.0,-1.0,w1,-1.0)
    return results

LOG.setLevel(logging.WARNING)
parser = getJacksParser()
args = parser.parse_args()

outprefix = args.outprefix
if '/' in outprefix and not os.path.exists(os.path.dirname(outprefix)): os.makedirs(os.path.dirname(outprefix))

outfile_w = outprefix + '_gene_results.txt'
outfile_w2 = outprefix + '_genestd_results.txt'

# Load the specification of samples to include
LOG.info('Loading sample specification')
sample_spec, ctrl_spec, sample_num_reps = createSampleSpec(args.countfile, args.replicatefile, args.rep_hdr,
                                                            args.sample_hdr, args.common_ctrl_sample, args.ctrl_sample_hdr)
# Load the mappings from guides to genes
LOG.info('Loading gene mappings')
gene_spec = createGeneSpec(args.guidemappingfile, args.sgrna_hdr, args.gene_hdr)

# Load the data and preprocess
LOG.info('Loading data and pre-processing')
data, meta, sample_ids, genes, gene_index = loadDataAndPreprocess(sample_spec, gene_spec)
gene_grnas = {gene: [x for x in meta[gene_index[gene], 0]] for gene in gene_index}
        
#Compute MeanFC for all samples against their controls
LOG.info('Running Single JACKS inference')
testdata, ctrldata, test_sample_idxs = collateTestControlSamples(data, sample_ids, ctrl_spec)
jacks_results = infer_JACKS_meanfc(gene_index, testdata, ctrldata)

# Write out the results
LOG.info('Writing Single JACKS results')
sample_ids_without_ctrl = [sample_ids[idx] for idx in test_sample_idxs]
writeJacksWResults(outprefix, jacks_results, sample_ids_without_ctrl, write_types=[''])

