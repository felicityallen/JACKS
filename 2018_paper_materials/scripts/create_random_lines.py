import io
import os
import sys
import csv
import numpy as np
import random

if len(sys.argv) != 6:
    print 'Usage: create_random_lines.py input_file repmap_file:rep_hdr:sample_hdr num_rand rep_per_rand ctrl_sample'
else:

    infile = sys.argv[1]
    outfile = infile[:-4] + '_rand' + infile[-4:]
    repmap_file, rep_hdr, sample_hdr = sys.argv[2].split(':')
    out_repmap_file = repmap_file[:-4] + '_rand' + repmap_file[-4:]
    num_rand = eval(sys.argv[3])
    rep_per_rand = eval(sys.argv[4])
    ctrl_sample = sys.argv[5]
    
    #Load replicate map
    f = io.open(repmap_file)
    rep_delim = ',' if (repmap_file.split('.')[-1] == 'csv') else '\t'
    rep_map = {row[rep_hdr]:row[sample_hdr] for row in csv.DictReader(f, delimiter=rep_delim)}
    f.close()
    
    #Randomly select existing cell lines to randomise
    f = io.open(infile)
    delim = ',' if (infile.split('.')[-1] == 'csv') else '\t'
    rdr = csv.DictReader(f, delimiter=delim)
    ctrl_hdrs = [hdr for hdr in rep_map if rep_map[hdr] == ctrl_sample]
    non_ctrl_hdrs = [hdr for hdr in rep_map if rep_map[hdr] != ctrl_sample]
    rand_spec = {'RAND%d' % idx: random.sample(non_ctrl_hdrs, rep_per_rand) for idx in range(num_rand)}
    rand_hdrs = sum([rand_spec[x] for x in rand_spec], [])
    ctrl_data, rand_data = [], []
    for row in rdr:
        ctrl_data.append([eval(row[hdr]) for hdr in ctrl_hdrs])
        rand_data.append([eval(row[hdr]) for hdr in rand_hdrs])
    ctrl_data = np.log2(np.array(ctrl_data) + 32)
    rand_data = np.log2(np.array(rand_data) + 32)
    f.close()
    
    #Median normalise each column and compute the fold changes
    ctrl_meds, rand_meds = np.nanmedian(ctrl_data,axis=0), np.nanmedian(rand_data,axis=0)
    ctrl_means = np.nanmean(ctrl_data - ctrl_meds,axis=1)
    rand_fcs = ((rand_data- rand_meds).T - ctrl_means).T
  
    #Randomly shuffle the guides within each column of the random samples
    for i in range(rand_fcs.shape[1]):
        np.random.shuffle(rand_fcs[:,i])
        
    #Convert the fold changes back to read counts 
    rand_counts = np.round(2.0**((rand_fcs.T + ctrl_means).T+rand_meds)-32)
    rand_counts *= (rand_counts>=0.0)
    
    #Write the output count file (appending random samples at the end)
    f,fout = io.open(infile), io.open(outfile, 'w')
    rdr = csv.reader(f, delimiter=delim); hdrs = rdr.next()
    rand_hdrs = sum([[rand_id + ' ' + hdr for hdr in rand_spec[rand_id]] for rand_id in rand_spec], []) 
    fout.write(u'%s\n' % delim.join(hdrs + rand_hdrs))
    for i, toks in enumerate(rdr):
        fout.write(u'%s\n' % delim.join(toks + ['%d' % x for x in rand_counts[i,:]]))
    f.close(), fout.close()

    #Write the new rep_map file (appending random samples)
    f, fout = io.open(repmap_file), io.open(out_repmap_file, 'w')
    fout.write(u'%s%s%s\n' % (rep_hdr,rep_delim,sample_hdr))
    for row in csv.DictReader(f, delimiter=rep_delim): 
        fout.write(u'%s%s%s\n' % (row[rep_hdr],rep_delim,row[sample_hdr]))
    for rand_id in rand_spec:
        fout.write(u'%s\n' % '\n'.join(['%s%s%s' % (rand_id + ' ' + hdr, rep_delim, rand_id) for hdr in rand_spec[rand_id]]))
    f.close(), fout.close()
    
