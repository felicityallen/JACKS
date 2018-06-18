import io, os, csv, sys, random, logging
from jacks.jacks import infer_JACKS, LOG
from jacks.io_preprocess import load_data_and_preprocess, writeJacksWResults, writeJacksXResults, pickleJacksFullResults
import numpy as np
from run_JACKS import prepareFile, createSampleSpec, createGeneSpec
    
def loadSgrnaReference(filename):
    f = io.open(filename)
    x_ref = {row['sgrna']:row for row in csv.DictReader(f,delimiter='\t')}
    f.close()
    return x_ref

if __name__ == '__main__':
    
    LOG.setLevel(logging.WARNING)
    
    if len(sys.argv) != 6:
        print('Usage: run_JACKS_reference.py countfile replicatefile:rep_hdr:sample_hdr:ctrl_sample_or_hdr guidemappingfile:sgrna_hdr:gene_hdr sgrna_reference_file outputprefix')
    else:

        #Parse arguments
        countfile = sys.argv[1] #sgrna label is always in the first column

        rep_toks = sys.argv[2].split(':')
        if len(rep_toks) != 4:
            raise Exception('Incorrect replicate file input: expecting "replicatefile:rep_hdr:sample_hdr:ctrl_sample" where replicatefile is a csv or tab delimited file mapping replicates to samples, rep_hdr and sample_hdr specify the column headers for the columns containing the replicate labels and sample labels respectively, and ctrl_sample specifies the name of the control sample')
        replicatefile, rep_hdr, sample_hdr, ctrl_sample_or_hdr = rep_toks
        
        guide_toks = sys.argv[3].split(':')
        if len(guide_toks) != 3:
            raise Exception('Incorrect guidemappingfile input: expecting "guidemappingfile:sgrna_hdr:gene_hdr" where guidemappingfile is a csv or tab delimited file mapping guides to genes, sgrna_hdr and gene_hdr specify the column headers for the columns containing the guide labels and gene labels respectively.')
        guidemappingfile, sgrna_hdr, gene_hdr = guide_toks   
        
        sgrna_reference_file = sys.argv[4]
        
        outprefix = sys.argv[5]
        if '/' in outprefix and not os.path.exists(os.path.dirname(outprefix)): os.makedirs(os.path.dirname(outprefix))   
        outfile_w = outprefix + '_gene_JACKS_results.txt'
        outfile_w2 = outprefix + '_genestd_JACKS_results.txt'
        outfile_x = outprefix + '_grna_JACKS_results.txt'
        outfile_pickle = outprefix + '_JACKS_results_full.pickle'
        
        #Load the sgrna reference (precomputed X's)
        print('Loading sgrna reference values')
        x_ref = loadSgrnaReference(sgrna_reference_file)
        
        #Load the specification of samples to include
        print('Loading sample specification')
        sample_spec, ctrl_per_sample, ctrl_spec = createSampleSpec(countfile, replicatefile, rep_hdr, sample_hdr, ctrl_sample_or_hdr)
        
        #Load the mappings from guides to genes
        print('Loading gene mappings')
        gene_spec = createGeneSpec(guidemappingfile, sgrna_hdr, gene_hdr)
        
        #Check that the data to be loaded have sgrna reference values
        print('Checking sgrna reference identifiers against gene mappings')
        for guide in gene_spec:
            if guide not in x_ref: 
                raise Exception('%s has no sgrna reference in %s' % (guide, sgrna_reference_file))
        
        #Load the data and preprocess (or just load from pickle if we did this already)
        print('Loading data and pre-processing')
        data, meta, sample_ids, genes, gene_index = load_data_and_preprocess(sample_spec, gene_spec)
        gene_grnas = {gene: [x for x in meta[gene_index[gene],0]] for gene in gene_index}
            
        #Create the X reference (in the correct order)
        x_reference = {'X1': np.array([eval(x_ref[x]['X1']) for x in meta[:,0]]),
                       'X2': np.array([eval(x_ref[x]['X2']) for x in meta[:,0]])}
            
        #Run all samples against the control
        print('Running JACKS inference')
        if ctrl_per_sample:     #Different control samples specified per test sample
            test_sample_idxs = [i for i,x in enumerate(sample_ids) if ctrl_spec[x] != x]
            testdata = data[:,test_sample_idxs,:]
            ctrldata = data[:,[sample_ids.index(ctrl_spec[sample_ids[idx]]) for idx in test_sample_idxs],:]
        else:                   #Same control sample for all tests
            ctrldata = data[:,sample_ids.index(ctrl_sample_or_hdr),:]
            test_sample_idxs = [i for i,x in enumerate(sample_ids) if x != ctrl_sample_or_hdr]
            testdata = data[:,test_sample_idxs,:]
        jacks_results = infer_JACKS(gene_index, testdata, ctrldata, fixed_x=x_reference)

        #Write out the results
        print('Writing JACKS results')
        sample_ids_without_ctrl = [sample_ids[idx] for idx in test_sample_idxs]
        writeJacksWResults( outfile_w, jacks_results, sample_ids_without_ctrl)
        writeJacksWResults( outfile_w2, jacks_results, sample_ids_without_ctrl, write_w2=True)
        writeJacksXResults( outfile_x, jacks_results, gene_grnas )
        pickleJacksFullResults( outfile_pickle, jacks_results, sample_ids_without_ctrl, gene_grnas )       
