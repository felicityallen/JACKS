import io, os, csv, random
from jacks.jacks import infer_JACKS
from jacks.io_preprocess import load_data_and_preprocess, writeJacksWResults, writeJacksXResults, pickleJacksFullResults

from run_JACKS import prepareFile, createGeneSpec, createSampleSpec
    
LOG.setLevel(logging.WARNING)

def loadSgrnaReference(filename):
    f = io.open(filename)
    x_ref = {row['sgrna']:row for row in csv.DictReader(f,delimiter='\t')}
    f.close()
    return x_ref

def createModifiedRepMapFile(repfile, outrepfile, rep_hdr, sample_hdr, ctrl_sample, sample, sample_only=False):
    f, delim = prepareFile(repfile, rep_hdr)
    fout = io.open(outrepfile, 'w')
    fout.write(u'%s\t%s\n' % (rep_hdr, sample_hdr))
    for row in csv.DictReader(f, delimiter=delim):
        if (row[sample_hdr] == ctrl_sample) or \
           (sample_only and row[sample_hdr] == sample) or \
           (not sample_only and row[sample_hdr] != sample):
                fout.write(u'%s\t%s\n' % (row[rep_hdr], row[sample_hdr]))
    f.close()
    fout.close()
    
if len(sys.argv) != 6:
    print 'Usage: leaveoneout_reference_test.py sample countfile replicatefile:rep_hdr:sample_hdr:ctrl_sample guidemappingfile:sgrna_hdr:gene_hdr outputprefix'
else:

    #Parse arguments
    sample = sys.argv[1]
    countfile = sys.argv[2] #sgrna label is always in the first column

    rep_toks = sys.argv[3].split(':')
    if len(rep_toks) != 4:
        raise Exception('Incorrect replicate file input: expecting "replicatefile:rep_hdr:sample_hdr:ctrl_sample" where replicatefile is a csv or tab delimited file mapping replicates to samples, rep_hdr and sample_hdr specify the column headers for the columns containing the replicate labels and sample labels respectively, and ctrl_sample specifies the name of the control sample')
    replicatefile, rep_hdr, sample_hdr, ctrl_sample = rep_toks
    
    guide_toks = sys.argv[4].split(':')
    if len(guide_toks) != 3:
        raise Exception('Incorrect guidemappingfile input: expecting "guidemappingfile:sgrna_hdr:gene_hdr" where guidemappingfile is a csv or tab delimited file mapping guides to genes, sgrna_hdr and gene_hdr specify the column headers for the columns containing the guide labels and gene labels respectively.')
    guidemappingfile, sgrna_hdr, gene_hdr = guide_toks   
    
    outprefix = sys.argv[5]
    
    #Create a replicate map file containing all but the specified sample
    mod_repfile_loo = outprefix + '_loo_repmap.txt'
    createModifiedRepMapFile(replicatefile, mod_repfile_loo, rep_hdr, sample_hdr, ctrl_sample, sample, sample_only=False)
    
    #Run JACKS on this
    if not os.path.isdir(outprefix + '_loo'): os.mkdir(outprefix + '_loo')
    cmd = 'python run_JACKS.py %s %s:%s:%s:%s %s %s' % (countfile, mod_repfile_loo, rep_hdr, sample_hdr, ctrl_sample, sys.argv[4], outprefix + '_loo' )
    print cmd
    os.system(cmd)
    xref_file = outprefix + '_loo_grna_JACKS_results.txt' 
    
    #Create a replicate map file containing only the specified sample and ctrl sample
    mod_repfile_only = outprefix + '_only_repmap.txt'
    createModifiedRepMapFile(replicatefile, mod_repfile_only, rep_hdr, sample_hdr, ctrl_sample, sample, sample_only=True)
    
    #Now run JACKS reference on the just computed sgrna data
    cmd = 'python run_JACKS_reference.py %s %s:%s:%s:%s %s %s %s' % (countfile, mod_repfile_only, rep_hdr, sample_hdr, ctrl_sample, sys.argv[4], xref_file, outprefix + '_only' )
    print cmd
    os.system(cmd)    
