import io, os, csv, random
import sys; sys.path.append('../..')    #path to jacks
from jacks.code.io_preprocess import load_data_and_preprocess, writeJacksWResults, pickleJacksFullResults

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

def infer_MeanFC(gene_index, testdata, ctrldata):
    results = {}
    for gene in gene_index:
        Ig = gene_index[gene]
        y = (testdata[Ig,:,0].T - ctrldata[Ig,0]).T
        w1 = SP.nanmean(y,axis=0)
        results[gene] = (y,-1.0,-1.0,-1.0,w1,-1.0)
    return results
    
LOG.setLevel(logging.WARNING)

if len(sys.argv) != 5:
    print 'Usage: run_MeanFC.py countfile replicatefile:rep_hdr:sample_hdr:ctrl_sample guidemappingfile:sgrna_hdr:gene_hdr outputprefix'
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
    outfile_w = outprefix + '_gene_meanfc_results.txt'
    outfile_pickle = outprefix + '_meanfc_results_full.pickle'
    
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
    
    #Run all samples against the control
    print 'Running MeanFC inference'
    ctrldata = data[:,sample_ids.index(ctrl_sample),:]
    testdata = data[:,[i for i,x in enumerate(sample_ids) if ctrl_sample not in x],:]
    meanfc_results = infer_MeanFC(gene_index, testdata, ctrldata)

    #Write out the results
    print 'Writing MeanFC results'
    sample_ids_without_ctrl = [x for x in sample_ids if x != ctrl_sample]
    writeJacksWResults( outfile_w, meanfc_results, sample_ids_without_ctrl)
    pickleJacksFullResults( outfile_pickle, meanfc_results, sample_ids_without_ctrl, gene_grnas )       
