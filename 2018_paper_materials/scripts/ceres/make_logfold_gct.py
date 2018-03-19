import io
import csv
import sys
import numpy as np

if len(sys.argv) != 3:
	raise Exception('Usage: make_logfold_gct.py infile repmapfile:rep_hdr:sample_hdr:ctrl_sample')
else:

	infile = sys.argv[1]
	repmapfile, rep_hdr, sample_hdr,ctrl_sample = sys.argv[2].split(':')
	
	#Load replicate map
	f = io.open(repmapfile)
	rep_delim = ',' if repmapfile[-4:] == '.csv' else '\t'
	rep_map = {row[rep_hdr]:row for row in csv.DictReader(f, delimiter=rep_delim)}
	f.close() 

	num_ctrls = len([x for x in rep_map if rep_map[x][sample_hdr] == ctrl_sample])
	if num_ctrls == 0: raise Exception("No control replicates found for sample " + ctrl_sample)
	if num_ctrls > 2: usebatch = True 
	else: usebatch = False

	#Read in counts
	f = io.open(infile)
	delim = ',' if infile[-4:] == '.csv' else '\t'
	rdr = csv.reader(f, delimiter=delim)
	hdrs = rdr.next()
	has_gene_col = (hdrs[1].lower() == 'gene')
	shdrs = hdrs[1+has_gene_col:]
	guides, genes, counts = [],[],[]
	for toks in rdr:
		guides.append(toks[0])
		if has_gene_col: genes.append(toks[1])
		counts.append([eval(x) for x in toks[1+has_gene_col:]])
	f.close()
	counts = np.array(counts)

	total_counts = np.sum(counts, axis=0)

	#Scale counts to 1 million reads per replicate
	G,L = counts.shape
	counts = counts*1e6/np.tile(total_counts,(G,1)) 

	#Compute log-fold changes from controls
	ctrl_hdrs = [x for x in rep_map if rep_map[x][sample_hdr] == ctrl_sample]; ctrl_hdrs.sort()
	if usebatch:
		ctrl_counts = np.zeros(counts.shape)
		for i, hdr in enumerate(shdrs):
			batch = eval(rep_map[hdr]['Batch'])
			ctrl_counts[:,i] = counts[:,shdrs.index(ctrl_hdrs[batch])]
	else:
		ctrl_counts = np.tile(np.nanmean(counts[:,[shdrs.index(x) for x in ctrl_hdrs]],axis=1),(1,L))
	lognorm_counts = np.log2(counts+1.0)-np.log2(ctrl_counts+1.0)
	
	#Remove ctrl columns amd columns not in the replicate map file
	lognorm_counts = lognorm_counts[:,[i for i,hdr in enumerate(shdrs) if hdr not in ctrl_hdrs and hdr in rep_map]]
	non_ctrl_shdrs = [x for x in shdrs if x not in ctrl_hdrs]
	G,L = lognorm_counts.shape

	#Apply ZMAD
	lognorm_counts -= np.tile(np.nanmedian(lognorm_counts,axis=0),(G,1))
	lognorm_counts = lognorm_counts/np.tile(1.4826*np.nanmedian(abs(lognorm_counts), axis=0),(G,1))

	#Write out the gct
	fout = io.open(infile[:-4] + '_lfnorm.gct','w')
	fout.write(u'#1.2\n%d\t%d\n' % (G,L))
	fout.write(u'Name\tDescription\t%s\n' % '\t'.join(non_ctrl_shdrs))
	for i in range(G):
		fout.write(u'%s\t%s\t%s\n' % (guides[i],guides[i],'\t'.join(['%5e' % x for x in lognorm_counts[i,:]])))
	fout.close()





