
import io, os, csv, random, sys
from jacks.jacks_io import createSampleSpec, createGeneSpec, getJacksParser, prepareFile
from jacks.infer import LOG

from BAGEL_calc_foldchange import calcBagelFoldchange
from BAGEL import runBAGEL

def addGeneColumn(filename, outfile, gene_spec, gene_hdr):
    delim = ',' if (filename.split('.')[-1] == 'csv') else '\t'
    f = io.open(filename); rdr = csv.DictReader(f, delimiter=delim)
    hdrs = rdr.fieldnames
    if hdrs[1] != gene_hdr:
        fout = io.open(outfile, 'w')
        new_hdrs = [hdrs[0]] + [gene_hdr] + hdrs[1:]
        fout.write(delim.join(new_hdrs) + u'\n')
        for row in rdr:
            grna_id = row[hdrs[0]]
            gene = gene_spec[grna_id] if gene_spec[grna_id] != '' else 'NO_GENE'
            fout.write(delim.join([grna_id] + [gene] + [row[x] for x in hdrs[1:]]) + u'\n')
        fout.close()
        return outfile
    return filename

py_cmd = 'run_python_2.7.sh'
parser = getJacksParser()
parser.add_argument("--sample_id",
                type=str,
                default=None,
                help="Sample id to run BAGEL on")
parser.add_argument("--v10",
                type=str,
                default='',
                help="Data set label")
args = parser.parse_args()

fc_dir = 'fold_change_files'
if not os.path.isdir(fc_dir): os.makedirs(fc_dir)

inputs_dir = 'input_files'
if not os.path.isdir(inputs_dir): os.makedirs(inputs_dir)

# Load the specification of samples to include
LOG.info('Loading sample specification')
sample_spec, ctrl_spec, sample_num_reps = createSampleSpec(args.countfile, args.replicatefile, args.rep_hdr,
                                                            args.sample_hdr, args.common_ctrl_sample, args.ctrl_sample_hdr)
# Load the mappings from guides to genes
LOG.info('Loading gene mappings')
gene_spec = createGeneSpec(args.guidemappingfile, args.sgrna_hdr, args.gene_hdr)

# Sample not specified: re-call self for all samples
if args.sample_id is None:
    for sample_id in ctrl_spec:
        if ctrl_spec[sample_id] == sample_id: continue
        cmd = py_cmd + ' ' + ' '.join(sys.argv) + ' --sample_id="%s"' % sample_id
        os.system(cmd)

#Sample specified - run BAGEL
else:
    sample_id = args.sample_id

    out_dir = 'bagel_single_screens_%s_%s/%s_%s_1' % (args.v10, sample_id, args.v10, sample_id)
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    #Collect sample information
    sample_filenames = set()
    for filename in sample_spec:
        ctrl_col_idxs, sample_colnames = [], []
        outfile = inputs_dir + '/' + filename.split('/')[-1]
        used_filename = addGeneColumn(filename, outfile, gene_spec, args.gene_hdr)
        delim = ',' if (used_filename.split('.')[-1] == 'csv') else '\t'
        f = io.open(used_filename)
        hdrs = csv.DictReader(f, delimiter=delim).fieldnames
        f.close()
        for (spec_sample_id, colname) in sample_spec[filename]:
            if spec_sample_id == sample_id:
                sample_colnames.append(colname)
            elif spec_sample_id == ctrl_spec[sample_id]:
                ctrl_col_idxs.append(hdrs.index(colname)-1)
            else: continue
            sample_filenames.add(used_filename)

    if len(sample_filenames) == 0:
        raise Exception('Could not find sample %s in sample spec' % sample_id)
    elif len(sample_filenames) > 1:
        raise Exception('Multiple input files containing %s - not supported!' % sample_id)

    #Compute the fold change file (will recompute for all samples in case controls are different per sample)
    infile = sample_filenames.pop()
    fcfile = fc_dir + '/' + infile.split('/')[-1][:-4] + '_' + sample_id
    sys.argv = ['BAGEL-calc_foldchange', '-i', infile, '-o', fcfile, '-c', "%s" % str(ctrl_col_idxs) ]
    calcBagelFoldchange()
    fcfile += '.foldchange'

    #Run BAGEL
    f = io.open(fcfile); rdr = csv.DictReader(f, delimiter=delim)
    hdrs = rdr.fieldnames
    f.close()
    sample_col_idxs = [str(hdrs.index(x)-1) for x in sample_colnames]
    sys.argv = ['BAGEL','-i', fcfile, '-o', out_dir + '/' + sample_id + '_1.txt', '-e','../../data/training_essentials.txt','-n','../../data/training_nonessential.txt','-c', ','.join(sample_col_idxs)]
    runBAGEL()
