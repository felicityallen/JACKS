import io, os, csv, random, sys, logging
from jacks.jacks import infer_JACKS, LOG
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
    

if __name__ == '__main__':

    LOG.setLevel(logging.WARNING)

    if len(sys.argv) != 5 and len(sys.argv) != 6:
        print 'Usage: run_JACKS.py countfile replicatefile:rep_hdr:sample_hdr:ctrl_sample guidemappingfile:sgrna_hdr:gene_hdr outputprefix apply_w_hp'
    else:

        #Parse arguments
        countfile = sys.argv[1] #Note: Assumes sgrna label is always in the first column!!!
        apply_w_hp = False
        if len(sys.argv) == 6: apply_w_hp = eval(sys.argv[5])

        rep_toks = sys.argv[2].split(':')
        if len(rep_toks) != 4:
            raise Exception('Incorrect replicate file input: expecting "replicatefile:rep_hdr:sample_hdr:ctrl_sample" where replicatefile is a csv or tab delimited file mapping replicates to samples, rep_hdr and sample_hdr specify the column headers for the columns containing the replicate labels and sample labels respectively, and ctrl_sample specifies the name of the control sample')
        replicatefile, rep_hdr, sample_hdr, ctrl_sample = rep_toks
        
        guide_toks = sys.argv[3].split(':')
        if len(guide_toks) != 3:
            raise Exception('Incorrect guidemappingfile input: expecting "guidemappingfile:sgrna_hdr:gene_hdr" where guidemappingfile is a csv or tab delimited file mapping guides to genes, sgrna_hdr and gene_hdr specify the column headers for the columns containing the guide labels and gene labels respectively.')
        guidemappingfile, sgrna_hdr, gene_hdr = guide_toks   
        
        outprefix = sys.argv[4]
        if '/' in outprefix and not os.path.isdir(outprefix.split('/')[0]): os.mkdir(outprefix.split('/')[0])
        outfile_w = outprefix + '_gene_JACKS_results.txt'
        outfile_x = outprefix + '_grna_JACKS_results.txt'
        outfile_pickle = outprefix + '_JACKS_results_full.pickle'
        
        #Load the specification of samples to include
        print 'Loading sample specification'
        sample_spec = createSampleSpec(countfile, replicatefile, rep_hdr, sample_hdr)
        
        #Load the mappings from guides to genes
        print 'Loading gene mappings'
        gene_spec = createGeneSpec(guidemappingfile, sgrna_hdr, gene_hdr)
            
        #Load the data and preprocess
        print 'Loading data and pre-processing'
        data, meta, sample_ids, genes, gene_index = load_data_and_preprocess(sample_spec, gene_spec)
        gene_grnas = {gene: [x for x in meta[gene_index[gene],0]] for gene in gene_index}
        
        #Run all samples against the control
        print 'Running JACKS inference'
        ctrldata = data[:,sample_ids.index(ctrl_sample),:]
        testdata = data[:,[i for i,x in enumerate(sample_ids) if ctrl_sample not in x],:]
        jacks_results = infer_JACKS(gene_index, testdata, ctrldata, apply_w_hp=apply_w_hp)

        #Write out the results
        print 'Writing JACKS results'
        sample_ids_without_ctrl = [x for x in sample_ids if x != ctrl_sample]
        writeJacksWResults( outfile_w, jacks_results, sample_ids_without_ctrl)
        writeJacksXResults( outfile_x, jacks_results, gene_grnas )
        pickleJacksFullResults( outfile_pickle, jacks_results, sample_ids_without_ctrl, gene_grnas )       
