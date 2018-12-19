
import io, os, csv, random, sys, re
from jacks.jacks_io import createSampleSpec, createGeneSpec, getJacksParser, prepareFile
from jacks.infer import LOG
import pandas as pd

def convertColName(colname):
    return colname
    colname = re.sub('[^0-9a-zA-Z_]', '.', colname)
    if colname[0] in ['0','1','2','3','4','5','6','7','8','9']: colname = 'X' + colname
    return colname


def addGeneColumn(filename, outfile, gene_spec, gene_hdr):
    delim = ',' if (filename.split('.')[-1] == 'csv') else '\t'
    f = io.open(filename); rdr = csv.DictReader(f, delimiter=delim)
    hdrs = rdr.fieldnames
    if hdrs[1] != gene_hdr:
        fout = io.open(outfile, 'w')
        new_hdrs = [hdrs[0]] + [gene_hdr] + hdrs[1:]
        fout.write('\t'.join(new_hdrs) + u'\n')
        for row in rdr:
            grna_id = row[hdrs[0]]
            gene = gene_spec[grna_id] if grna_id in gene_spec and gene_spec[grna_id] != '' else 'NO_GENE'
            fout.write('\t'.join([grna_id] + [gene] + [row[x] for x in hdrs[1:]]) + u'\n')
        fout.close()
        return outfile
    return filename

py_cmd = 'python'
python_local = 'python'
parser = getJacksParser()
parser.add_argument("--sample_id",
                type=str,
                default=None,
                help="Sample id to run PBNPA on")
parser.add_argument("--v10",
                type=str,
                default='',
                help="Data set label")
args = parser.parse_args()

inputs_dir = 'input_files'
if not os.path.isdir(inputs_dir): os.makedirs(inputs_dir)

tmp_dir = 'tmp_files'
if not os.path.isdir(tmp_dir):
    os.mkdir(tmp_dir)

# Load the specification of samples to include
LOG.info('Loading sample specification')
sample_spec, ctrl_spec, sample_num_reps = createSampleSpec(args.countfile, args.replicatefile, args.rep_hdr,
                                                            args.sample_hdr, args.common_ctrl_sample, ctrl_sample_hdr=args.ctrl_sample_hdr)
# Load the mappings from guides to genes
LOG.info('Loading gene mappings')
gene_spec = createGeneSpec(args.guidemappingfile, args.sgrna_hdr, args.gene_hdr)

# Sample not specified: re-call self for all samples
if args.sample_id is None:
    for sample_id in ctrl_spec:
        if ctrl_spec[sample_id] == sample_id: continue
        cmd = py_cmd + ' ' + ' '.join(sys.argv) + ' --sample_id=%s' % sample_id
        os.system(cmd)

#Sample specified - run PBNPA
else:
    sample_id = args.sample_id

    out_dir = 'PBNPA_single_screens_%s_%s/%s_%s_1' % (args.v10, sample_id, args.v10, sample_id)
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    #Collect sample information
    sample_filenames = set()
    for filename in sample_spec:
        ctrl_colnames, sample_colnames = [], []
        outfile = inputs_dir + '/' + filename.split('/')[-1][:-4] + '.txt'
        used_filename = addGeneColumn(filename, outfile, gene_spec, args.gene_hdr)
        delim = ',' if (used_filename.split('.')[-1] == 'csv') else '\t'
        f = io.open(used_filename)
        hdrs = csv.DictReader(f, delimiter=delim).fieldnames
        f.close()
        for (spec_sample_id, colname) in sample_spec[filename]:
            if spec_sample_id == sample_id:
                sample_colnames.append(convertColName(colname))
            elif spec_sample_id == ctrl_spec[sample_id]:
                ctrl_colnames.append(convertColName(colname))
            else: continue
            sample_filenames.add(used_filename)

    if len(sample_filenames) == 0:
        raise Exception('Could not find sample %s in sample spec' % sample_id)
    elif len(sample_filenames) > 1:
        raise Exception('Multiple input files containing %s - not supported!' % sample_id)

    if len(sample_colnames) == 0:
        print(sample_id)
        print(sample_spec[filename])
        raise Exception('Could not find sample columns for %s' % sample_id)
    if len(ctrl_colnames) == 0:
        print(sample_id)
        print(ctrl_spec)
        raise Exception('Could not find ctrl columns for %s' % sample_id)

    infile = sample_filenames.pop()
    infiles= []
    data = pd.read_csv(infile, sep='\t')
    for sample_col in sample_colnames:   
        data['sgRNA'] = data[data.columns[0]]
        data['Gene'] = data[data.columns[1]]
        rep_PBNPA_input_file = tmp_dir + '/%s_%s' % (infile.split('/')[-1], sample_col.replace(' ','_').replace('/','_').replace(',','_'))
        data[['sgRNA','Gene', ctrl_colnames[0], sample_col]].to_csv(rep_PBNPA_input_file, index=False, header=True)
        infiles.append(rep_PBNPA_input_file)
    cmd = 'Rscript PBNPA/R/PBNPA.R "%s" "%s.csv"' % (','.join(infiles), out_dir + '/' + sample_id + '_1')
    print(cmd); os.system(cmd)

