import argparse
import csv
import io
import logging
import os

import numpy as np

from jacks.io_preprocess import load_data_and_preprocess, writeJacksWResults, writeJacksXResults, \
    pickleJacksFullResults, writeFoldChanges
from jacks.jacks import infer_JACKS, LOG


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

    # Reopen the file and skip to the start of the data
    f = io.open(filename)
    [f.readline() for i in range(skip_lines)]
    return f, delim


# output:  {input_filename:[(sample_id, colname)]}
def createSampleSpec(infile, repfile, rep_hdr, sample_hdr, ctrl_sample_or_hdr):
    f, delim = prepareFile(repfile, rep_hdr)
    rdr = csv.DictReader(f, delimiter=delim)
    sample_spec = {infile: []}
    ctrl_per_sample = (ctrl_sample_or_hdr in rdr.fieldnames)
    ctrl_spec = {}
    for row in rdr:
        sample_spec[infile].append((row[sample_hdr], row[rep_hdr]))
        if ctrl_per_sample:
            if row[sample_hdr] in ctrl_spec:
                if ctrl_spec[row[sample_hdr]] != row[ctrl_sample_or_hdr]:
                    err_msg = '%s vs %s for %s\n' % (
                        ctrl_spec[row[sample_hdr]], row[ctrl_sample_or_hdr], row[sample_hdr])
                    raise Exception(err_msg + 'Different controls for replicates of the sample not supported.')
            else:
                ctrl_spec[row[sample_hdr]] = row[ctrl_sample_or_hdr]
    f.close()
    return sample_spec, ctrl_per_sample, ctrl_spec


# output:  {grna: gene}
def createGeneSpec(guidemappingfile, sgrna_hdr, gene_hdr):
    f, delim = prepareFile(guidemappingfile, sgrna_hdr)
    rdr = csv.DictReader(f, delimiter=delim)
    gene_spec = {row[sgrna_hdr]: row[gene_hdr] for row in rdr}
    f.close()
    return gene_spec


def loadSgrnaReference(filename):
    f = io.open(filename)
    x_ref = {row['sgrna']: row for row in csv.DictReader(f, delimiter='\t')}
    f.close()
    return x_ref


def get_run_jacks_parser():
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
    replicatefile_group.add_argument("--rep-hdr",
                                     type=str,
                                     default="Replicate",
                                     help="Column header for columns containing the replicate labels")
    replicatefile_group.add_argument("--sample-hdr",
                                     type=str,
                                     default="Sample",
                                     help="Column header for columns containing the sample labels")
    replicatefile_group.add_argument("--ctrl-sample",
                                     type=str,
                                     default="CONTROL",
                                     help="Name of the control sample")
    guidemapping_group = ap.add_argument_group("Guidemapping file arguments")
    guidemapping_group.add_argument("guidemappingfile",
                                    help="A CSV or tab delimited file mapping guides to genes")
    guidemapping_group.add_argument("--sgrna-hdr",
                                    type=str,
                                    default="Guide",
                                    help="Column headers for the columns containing the guide labels")
    guidemapping_group.add_argument("--gene-hdr",
                                    type=str,
                                    default="Gene",
                                    help="Column headers for the columns containing the gene labels")
    ap.add_argument("--outprefix",
                    type=str,
                    default="",
                    help="Output prefix")
    settings_group = ap.add_mutually_exclusive_group(required=True)
    settings_group.add_argument("--reffile",
                                help="sgRNA_reference_file")
    settings_group.add_argument("--apply_w_hp",
                                action='store_true',
                                default=False,
                                help="Apply with hp")
    return ap


if __name__ == '__main__':

    LOG.setLevel(logging.WARNING)
    parser = get_run_jacks_parser()
    args = parser.parse_args()

    ctrl_sample_or_hdr = args.ctrl_sample
    outprefix = args.outprefix
    if '/' in outprefix and not os.path.exists(os.path.dirname(outprefix)): os.makedirs(os.path.dirname(outprefix))
    outfile_w = outprefix + '_gene_JACKS_results.txt'
    outfile_w2 = outprefix + '_genestd_JACKS_results.txt'
    outfile_x = outprefix + '_grna_JACKS_results.txt'
    outfile_lfc = outprefix + '_logfoldchange_means.txt'
    outfile_lfc_std = outprefix + '_logfoldchange_std.txt'
    outfile_pickle = outprefix + '_JACKS_results_full.pickle'

    # Load the specification of samples to include
    print('Loading sample specification')
    sample_spec, ctrl_per_sample, ctrl_spec = createSampleSpec(args.countfile, args.replicatefile, args.rep_hdr,
                                                               args.sample_hdr, ctrl_sample_or_hdr)

    # Load the mappings from guides to genes
    print('Loading gene mappings')
    gene_spec = createGeneSpec(args.guidemappingfile, args.sgrna_hdr, args.gene_hdr)

    sgrna_reference_file = args.reffile
    if sgrna_reference_file:
        # Load the sgrna reference (precomputed X's)
        print('Loading sgrna reference values')
        x_ref = loadSgrnaReference(args.reffile)
        # Check that the data to be loaded have sgrna reference values
        print('Checking sgrna reference identifiers against gene mappings')
        for guide in gene_spec:
            if guide not in x_ref:
                raise Exception('%s has no sgrna reference in %s' % (guide, sgrna_reference_file))

    # Load the data and preprocess
    print('Loading data and pre-processing')
    data, meta, sample_ids, genes, gene_index = load_data_and_preprocess(sample_spec, gene_spec)
    gene_grnas = {gene: [x for x in meta[gene_index[gene], 0]] for gene in gene_index}
    x_reference = None
    if sgrna_reference_file:
        # Create the X reference (in the correct order)
        x_reference = {'X1': np.array([eval(x_ref[x]['X1']) for x in meta[:, 0]]),
                       'X2': np.array([eval(x_ref[x]['X2']) for x in meta[:, 0]])}
    else:
        writeFoldChanges(outfile_lfc, data, meta, sample_ids)
        writeFoldChanges(outfile_lfc_std, data, meta, sample_ids, write_std=True)

    # Run all samples against their controls
    print('Running JACKS inference')
    if ctrl_per_sample:  # Different control samples specified per test sample
        test_sample_idxs = [i for i, x in enumerate(sample_ids) if ctrl_spec[x] != x]
        testdata = data[:, test_sample_idxs, :]
        ctrldata = data[:, [sample_ids.index(ctrl_spec[sample_ids[idx]]) for idx in test_sample_idxs], :]
    else:  # Same control sample for all tests
        ctrldata = data[:, sample_ids.index(ctrl_sample_or_hdr), :]
        test_sample_idxs = [i for i, x in enumerate(sample_ids) if x != ctrl_sample_or_hdr]
        testdata = data[:, test_sample_idxs, :]
    jacks_results = infer_JACKS(gene_index, testdata, ctrldata, apply_w_hp=args.apply_w_hp, fixed_x=x_reference)

    # Write out the results
    print('Writing JACKS results')
    sample_ids_without_ctrl = [sample_ids[idx] for idx in test_sample_idxs]
    writeJacksWResults(outfile_w, jacks_results, sample_ids_without_ctrl)
    writeJacksWResults(outfile_w2, jacks_results, sample_ids_without_ctrl, write_w2=True)
    writeJacksXResults(outfile_x, jacks_results, gene_grnas)
    pickleJacksFullResults(outfile_pickle, jacks_results, sample_ids_without_ctrl, gene_grnas)
