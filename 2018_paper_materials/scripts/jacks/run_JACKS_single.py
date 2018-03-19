import io, os, csv, random, sys
from jacks.jacks import infer_JACKS
from jacks.io_preprocess import load_data_and_preprocess, writeJacksWResults, writeJacksXResults, pickleJacksFullResults

def prepareFile(filename, hdr):
    #Count any lines before the headers (should be skipped)
    f = io.open(filename)
    skip_lines, line = 0, f.readline()  
    while hdr not in line and skip_lines < 100: skip_lines += 1; line = f.readline()
    f.close()
    
    if skip_lines >= 100:
        raise Exception('Could not find line with header ' + hdr + ' in ' + filename)
    
    #Check for comma vs tab delimited
    delim = ',' if (filename.split('.')[-1] == 'csv') else '\t'
    
    #Reopen the file and skip to the start of the data
    f = io.open(filename); [f.readline() for i in range(skip_lines)]  
    return f, delim

#output:  {input_filename:[(sample_id, colname)]}
def createSampleSpec(infile, repfile, rep_hdr, sample_hdr):
    f, delim = prepareFile(repfile, rep_hdr)
    rdr = csv.DictReader(f, delimiter=delim)
    sample_spec = {infile:[(row[sample_hdr],row[rep_hdr]) for row in rdr]}
    f.close()
    return sample_spec
    
#output:  {grna: gene}
def createGeneSpec(guidemappingfile, sgrna_hdr, gene_hdr):
    f, delim = prepareFile(guidemappingfile, sgrna_hdr)
    rdr = csv.DictReader(f, delimiter=delim)
    gene_spec = {row[sgrna_hdr]:row[gene_hdr] for row in rdr}
    f.close()
    return gene_spec
    
LOG.setLevel(logging.WARNING)

if len(sys.argv) != 5:
    print 'Usage: run_JACKS_single.py countfile replicatefile:rep_hdr:sample_hdr:ctrl_sample guidemappingfile:sgrna_hdr:gene_hdr outputprefix'
else:

    #Parse arguments
    countfile = sys.argv[1] #sgrna label is always in the first column

    rep_toks = sys.argv[2].split(':')
    if len(rep_toks) != 4:
        raise Exception('Incorrect replicate file input: expecting "replicatefile:rep_hdr:sample_hdr:ctrl_sample" where replicatefile is a csv or tab delimited file mapping replicates to samples, rep_hdr and sample_hdr specify the column headers for the columns containing the replicate labels and sample labels respectively, and ctrl_sample specifies the name of the control sample')
    replicatefile, rep_hdr, sample_hdr, ctrl_sample = rep_toks
    
    guide_toks = sys.argv[3].split(':')
    if len(guide_toks) != 3:
        raise Exception('Incorrect guidemappingfile input: expecting "guidemappingfile:sgrna_hdr:gene_hdr" where guidemappingfile is a csv or tab delimited file mapping guides to genes, sgrna_hdr and gene_hdr specify the column headers for the columns containing the guide labels and gene labels respectively.')
    guidemappingfile, sgrna_hdr, gene_hdr = guide_toks   
    
    outprefix = sys.argv[4]
    
    #Load the specification of samples to include
    print 'Loading sample specification'
    sample_spec = createSampleSpec(countfile, replicatefile, rep_hdr, sample_hdr)
    
    #Load the mappings from guides to genes
    print 'Loading gene mappings'
    gene_spec = createGeneSpec(guidemappingfile, sgrna_hdr, gene_hdr)
        
    #Load the data and preprocess (or just load from pickle if we did this already)
    print 'Loading data and pre-processing'
    data, meta, sample_ids, genes, gene_index = load_data_and_preprocess(sample_spec, gene_spec)
    gene_grnas = {gene: [x for x in meta[gene_index[gene],0]] for gene in gene_index}
    
    #Run each sample against the control
    ctrldata = data[:,sample_ids.index(ctrl_sample),:]
    for sample_id in sample_ids:
        if sample_id == ctrl_sample: continue

        outfile_w = outprefix + '_' + sample_id + '_gene_JACKS_results.txt'
        outfile_x = outprefix + '_' + sample_id + '_grna_JACKS_results.txt'
        outfile_pickle = outprefix + '_' + sample_id +'_JACKS_results_full.pickle'

        print 'Running JACKS inference for', sample_id
        testdata = data[:,[sample_ids.index(sample_id)],:]
        jacks_results = infer_JACKS(gene_index, testdata, ctrldata)

        #Write out the results
        print 'Writing JACKS results for', sample_id
        writeJacksWResults( outfile_w, jacks_results, [sample_id])
        writeJacksXResults( outfile_x, jacks_results, gene_grnas )
        pickleJacksFullResults( outfile_pickle, jacks_results, [sample_id], gene_grnas )       
