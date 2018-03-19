import sys, io, os, random, logging
from jacks.io_preprocess import subsample_and_preprocess
from jacks.jacks import infer_JACKS_gene, LOG
import numpy as np

def read_essentiality():
    ess = [l.strip().split("\t")[0] for l in file("../../data/Hart_training_essentials.txt",'r').readlines()[1:]]
    noness = [l.strip().split("\t")[0] for l in file("../../data/Hart_training_nonessential.txt",'r').readlines()[1:]]
    return {False:noness, True:ess}

if len(sys.argv) != 6 and len(sys.argv) != 7:
    print 'Usage: sample_jacks_screen_xs.py rawcountfile num_replicates num_celllines outfile num_bootraps <job_idx - opt, else LSB_JOBINDEX>'
else:

    LOG.setLevel(logging.WARNING)
    
    inputfile = sys.argv[1]
    num_replicates = eval(sys.argv[2])
    num_celllines = eval(sys.argv[3])
    outfile = sys.argv[4]
    num_bootstraps = eval(sys.argv[5])
    if len(sys.argv) == 6: job_idx = os.environ['LSB_JOBINDEX']
    else: job_idx = sys.argv[6]
    
    #Get list of test cell lines
    f = io.open(inputfile)
    hdrs = [x.split('_')[0] for x in f.readline().split('\t')[2:] if 'CTRL' not in x]
    f.close()
    test_celllines = [x for x in set(hdrs)]
    
    ess_genes = set(read_essentiality()[True])
    
    x_values = {}
    for bs in range(num_bootstraps):

        selected_celllines = random.sample(test_celllines, num_celllines)
        print selected_celllines
        
        # Read all data from input file
        data, meta, cell_lines, genes, gene_index, _ = subsample_and_preprocess(inputfile, [('CTRL',-1)] + [(x, num_replicates) for x in selected_celllines])
        ctrldata = data[:,cell_lines.index('CTRL'),:]
        testdata = data[:,[cell_lines.index(x) for x in selected_celllines],:]

        #Run JACKS on essential genes only, record the output x values
        for gene in gene_index:
            if gene not in ess_genes:
                continue
            Ig = gene_index[gene]

            # Perform inference
            y, tau, x1, x2, w1, w2 = infer_JACKS_gene(testdata[Ig,:,0], testdata[Ig,:,1], ctrldata[Ig,0], ctrldata[Ig,1], 50)

            for i,grna in enumerate(meta[Ig,0]):
                if bs == 0:
                    x_values[grna] = []       
                x_values[grna].append((x1[i],x2[i]))
       
    if not os.path.exists(os.path.dirname(outfile)): os.makedirs(os.path.dirname(outfile))     
    fout = io.open(outfile[:-4] + '_%s.txt' % job_idx, 'w')               
    fout.write(u'gRNA\t%s\n' % ('\t'.join(['X1\tX2' for i in range(num_bootstraps)])))
    for grna in x_values:
        fout.write(u'%s\t%s\n' % (grna, '\t'.join(['%.5e\t%5e' % (x1, x2) for (x1,x2) in x_values[grna]])))
    fout.close()
        
