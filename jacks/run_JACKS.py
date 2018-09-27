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


run_jacks_doc = """

Example: run_JACKS.py countfile replicatefile:rep_hdr:sample_hdr:ctrl_sample_or_hdr guidemappingfile:sgrna_hdr:gene_hdr
 --apply_w_hp --outputprefix prefix
"""
replicatefile_doc = """ 
            "replicatefile:rep_hdr:sample_hdr:ctrl_sample_or_hdr"
            
            where replicatefile is a csv or tab delimited file mapping replicates to samples, 
            rep_hdr and sample_hdr specify the column headers for the columns containing 
            the replicate labels and sample labels respectively, and ctrl_sample specifies 
            the name of the control sample
            """
guidemapping_doc = """
                "guidemappingfile:sgrna_hdr:gene_hdr"
                
                where guidemappingfile is a csv or tab delimited file mapping guides to genes, 
                sgrna_hdr and gene_hdr specify the column headers for the columns containing 
                the guide labels and gene labels respectively.
                """
ap = argparse.ArgumentParser(description=run_jacks_doc)
ap.add_argument("countfile",
                help="Countfile: assumes sgrna label is always in the first column!")
ap.add_argument("replicatefile",
                help=replicatefile_doc)
ap.add_argument("guidemappingfile",
                help=guidemapping_doc)
ap.add_argument("--outprefix",
                type=str,
                default="",
                help="Output prefix")
group = ap.add_mutually_exclusive_group(required=True)
group.add_argument("--reffile",
                   help="sgRNA_reference_file")
group.add_argument("--apply_w_hp",
                   action='store_true',
                   default=False,
                   help="Apply with hp")

if __name__ == '__main__':

    LOG.setLevel(logging.WARNING)
    args = ap.parse_args()

    countfile = args.countfile
    apply_w_hp = args.apply_w_hp
    rep_toks = args.replicatefile.split(':')
    if len(rep_toks) != 4:
        raise Exception(
            'Incorrect replicate file input: expecting "replicatefile:rep_hdr:sample_hdr:ctrl_sample_or_hdr" '
            + replicatefile_doc)
    replicatefile, rep_hdr, sample_hdr, ctrl_sample_or_hdr = rep_toks

    guide_toks = args.guidemappingfile.split(':')
    if len(guide_toks) != 3:
        raise Exception(
            'Incorrect guidemappingfile input: expecting "guidemappingfile:sgrna_hdr:gene_hdr" ' + guidemapping_doc)
    guidemappingfile, sgrna_hdr, gene_hdr = guide_toks

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
    sample_spec, ctrl_per_sample, ctrl_spec = createSampleSpec(countfile, replicatefile, rep_hdr, sample_hdr,
                                                               ctrl_sample_or_hdr)

    # Load the mappings from guides to genes
    print('Loading gene mappings')
    gene_spec = createGeneSpec(guidemappingfile, sgrna_hdr, gene_hdr)

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
    jacks_results = infer_JACKS(gene_index, testdata, ctrldata, apply_w_hp=apply_w_hp, fixed_x=x_reference)

    # Write out the results
    print('Writing JACKS results')
    sample_ids_without_ctrl = [sample_ids[idx] for idx in test_sample_idxs]
    writeJacksWResults(outfile_w, jacks_results, sample_ids_without_ctrl)
    writeJacksWResults(outfile_w2, jacks_results, sample_ids_without_ctrl, write_w2=True)
    writeJacksXResults(outfile_x, jacks_results, gene_grnas)
    pickleJacksFullResults(outfile_pickle, jacks_results, sample_ids_without_ctrl, gene_grnas)
