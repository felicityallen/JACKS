import os, glob
import io
import random, csv
import numpy as np

""" Ensures monotonicity of x"""
def monotonize(x):
    N = len(x)
    for i in range(1, N):
        if not np.isnan(x[N-i-1]):
            x[N-i-1] = np.nanmax(x[N-i-1:N-i+1])
    return x

""" Apply a mean averaging window across a sample with a shift of 1"""
def window_smooth(x, window=25):
    N = len(x)
    m = np.zeros(N) # mean variance in window
    for k in range(window+1): m[k] = np.mean(x[0:2*window+1])
    for k in range(window+1, N-window-1): m[k] = m[k-1] + (x[k+window+1]-x[k-window])/(2*window+1)
    for k in range(N-window-1, N): m[k] = np.mean(x[N-2*window-1:])
    return m
    
""" Simple posterior estimation - use the average variance from 2*window (default 50) gRNAs with nearest means as the estimate if larger than observed. """
def calc_posterior_sd(data, window=800, do_monotonize=True, window_estimate=True):
    N = len(data)
    
    #If there's only one replicate, set variances undefined (nans)
    if len(data.shape) == 1 or data.shape[1] == 1:
        return np.zeros(N)*np.nan

    #Sort by means, then (in case of equality) variances
    dmean_vars = [(x,y,i) for (x,y,i) in zip(data.mean(axis=1), np.nanstd(data, axis=1)**2, range(data.shape[0]))]
    dmean_vars.sort()
    vars  = np.array([y for (x,y,i) in dmean_vars])
    I = np.array([i for (x,y,i) in dmean_vars])
    
    # window smoothing of variances
    m = window_smooth(vars, window=window)
    # enforce monotonicity (decreasing variance with log count) if desired
    if do_monotonize: m = monotonize(m)   
    
    # construct the posterior estimate
    sd_post = np.zeros(N)
    for i,idx in enumerate(I):
        sd_post[idx] = max(m[i], vars[i])**0.5 # take the larger of window estimate and actual observation; right hand side all [i], since both m and s2s are sorted
        if window_estimate: sd_post[idx] = m[i]**0.5 # ignore if much larger
    return sd_post  

""" Select a specified number of replicates for each of a selection of screens
@param expts List of possible experiment labels to choose from (from the header of the input file)
@param selected_screens A list of tuples, each containing a screen identifier string (cell line) and the number of replicates 
       to select with replacement or -1 to use all replicates available without replacement
"""
def selectSamples(all_screens, selected_screens, norepl):
    col_idxs, Iscreens = [], []
    for screen, exp_per_screen in selected_screens:
        screen_col_idxs = [idx for idx,scr in enumerate(all_screens) if scr == screen]
        if len(screen_col_idxs) == 0: raise KeyError('No columns found matching screen ' + screen)
        if exp_per_screen < 0: 
            selected_col_idxs = screen_col_idxs
        else:
            if norepl:  #Sample from replicates without replacement
                if len(screen_col_idxs) < exp_per_screen:
                    selected_col_idxs = screen_col_idxs
                else:
                    selected_col_idxs = random.sample(screen_col_idxs, exp_per_screen)
            else:       #Sampler from replicates with replacement
                selected_col_idxs = [random.choice(screen_col_idxs) for i in range(exp_per_screen)]
        Iscreens.append(range(len(col_idxs), len(col_idxs)+len(selected_col_idxs)))
        col_idxs.extend(selected_col_idxs)
    return col_idxs, Iscreens
    
""" Function to compute estimated mean and variance values for each sample across replicates"""
def condense_normalised_counts(counts, Iscreens):
    data = np.zeros([counts.shape[0], len(Iscreens), 2]) # mean and variance per gRNA in each sample
    for s, Iscreen in enumerate(Iscreens): # condense to per sample values
        N = len(Iscreen)
        data[:,s,0] = counts[:,Iscreen].mean(axis=1) # mu_hat
        data[:,s,1] = calc_posterior_sd(counts[:,Iscreen]) #sigma_hat
    return data
    
""" Function to calculate additional helper projections"""
def compute_extras(meta, sample_ids):
    genes = np.unique(sorted(meta[:,1]))
    gene_guides = {g:[] for g in genes}
    for i in range(len(meta)): gene_guides[meta[i,1]].append(i)
    return genes, gene_guides
   
""" Function to normalize log counts """ 
def normalizeLogCounts(logcounts, normtype='median'):
    G,L = logcounts.shape
    if normtype == 'median':        
        logcounts -= np.tile(np.nanmedian(logcounts,axis=0),(G,1)) # median-normalize
    elif normtype == 'zmad':        
        logcounts -= np.tile(np.nanmedian(logcounts,axis=0),(G,1)) # median-normalize
        logcounts = logcounts/np.tile(1.4826*np.nanmedian(abs(logcounts), axis=0),(G,1))	#adjust to median absolute deviation = 1
    return logcounts

""" Select a specified number of replicates for each of a selection of screens from a larger file and generate the condensed data,
    the output is the same as for read_condensed in io.py
@param infile The input raw counts file
@param selected_screens A list of tuples, each containing a screen identifier string (cell line) and the number of replicates 
       to select or -1 to use all replicates available without replacement
@param prior The prior for each read count
@param norepl Whether to randomly select samples for each screen with (False) or without (True) replacement""" 
def subsample_and_preprocess(infile, selected_screens, prior=32, normtype='median', norepl=True):
    
    # Select data from the input file
    f = file(infile, 'r')
    expts = [x for x in f.readline().split('\t')]
    all_screens = [x.split('_')[0] for x in expts]  #NOTE: Assumes cell lines is first tok before '_' !!!
    col_idxs, Iscreens = selectSamples(all_screens, selected_screens, norepl)    
    test_line_reps = [expts[col_idxs[i]] for i in Iscreens[1]] #Assumes test line is the second listed (first is control)         

    #Read data from file for selected screens
    counts, meta = [], []
    screens = np.array([all_screens[idx] for idx in col_idxs])
    for toks in csv.reader(f, delimiter='\t'):
        counts.append([np.log2(eval(toks[idx])+prior) for idx in col_idxs])
        meta.append(toks[:2])
    counts = np.array(counts)
    counts = normalizeLogCounts(counts, normtype=normtype)

    #Compile into condensed format 
    cell_lines = [x for (x,y) in selected_screens]
    meta = np.array(meta)
    data = condense_normalised_counts(counts, Iscreens)
    genes, gene_guides = compute_extras(meta, cell_lines)
    return data, meta, cell_lines, genes, gene_guides, test_line_reps

""" Load counts data from a list of one or more files with the same guides (in the same order)

@param sample_spec is {filename:[(sample_id, colname)]} i.e. for each file a list of column names you
   want to load and their sample_ids (to group replicates)"""
def load_data_and_preprocess(sample_spec, gene_spec, normtype='median', prior=32):

    #Load the count data from all files as per sample_spec
    samples, sample_counts, metas = [], [], []
    for filename in sample_spec:
        counts, meta = [], []
        f = io.open(filename)
        delim = ',' if (filename.split('.')[-1] == 'csv') else '\t'
        rdr = csv.DictReader(f, delimiter=delim)
        sgrna_col = rdr.fieldnames[0]   #sgrna label is the first column
        samples.extend([sample_id for sample_id, colname in sample_spec[filename]])
        for row in rdr:
            if row[sgrna_col] not in gene_spec: 
                LOG.warning('No gene mapping found for %s (from file %s)' % (row[sgrna_col],filename))
                continue    #Ignore entries with no gene mapping
            counts.append([np.log2(eval(row[colname])+prior) for sample_id, colname in sample_spec[filename]])
            meta.append([row[sgrna_col], gene_spec[row[sgrna_col]]])
        f.close()   
        counts = np.array(counts)
        counts = normalizeLogCounts(counts, normtype=normtype)
        meta = np.array(meta)
        I = np.argsort(meta[:,0])
        meta, counts = meta[I,:], counts[I,:]
        metas.append(meta); sample_counts.append(counts + 0.0)
        if len(metas) > 1:
            if not (metas[0]==meta).all():
                import pdb; pdb.set_trace()
                raise Exception('Incompatible gRNAs in file %s e.g gRNA naming or order does not match other files' % filename)
    counts = np.concatenate(sample_counts, axis=1)
    
    #Condense counts into means, variances etc per sample
    sample_ids = [x for x in set(samples)]
    Iscreens = [[col_idx for col_idx,x in enumerate(samples) if x == sample_id] for sample_id in sample_ids]
    data = condense_normalised_counts(counts, Iscreens)
    genes, gene_guides = compute_extras(meta, sample_ids)

    return data, meta, sample_ids, genes, gene_guides

def writeJacksWResults( filename, jacks_results, cell_lines, write_w2=False ):
    #Sort genes by w1
    ordered_genes = [(np.nanmean(jacks_results[gene][4]),gene) for gene in jacks_results]
    ordered_genes.sort()
    fout = io.open(filename, 'w')
    fout.write(u'Gene\t%s\n' % ('\t'.join(cell_lines)))
    for w1_mean, gene in ordered_genes:
        if write_w2:
            w1_str = '\t'.join(['%5e' % np.sqrt(w2 - w1**2.0) for (w1,w2) in zip(jacks_results[gene][4],jacks_results[gene][5])])
        else:
            w1_str = '\t'.join(['%5e' % w1 for w1 in jacks_results[gene][4]])
        fout.write(u'%s\t%s\n' % (gene, w1_str))
    fout.close()

def writeJacksXResults( filename, jacks_results, gene_grnas ):
    #Sort genes by mean w1
    ordered_genes = [(np.nanmean(jacks_results[gene][4]),gene) for gene in jacks_results]
    ordered_genes.sort()
    
    fout = io.open(filename, 'w')
    fout.write(u'sgrna\tX1\tX2\n')
    for w1_mean, gene in ordered_genes:
        y, tau, x1, x2, w1, w2 = jacks_results[gene]
        for i, grna in enumerate(gene_grnas[gene]):
            fout.write(u'%s\t%5e\t%5e\n' % (grna, x1[i], x2[i]))
    fout.close()   

def writeFoldChanges(filename, data, meta, sample_ids, write_std=False):
    fout = io.open(filename, 'w')
    fout.write(u'gRNA\tgene\t%s\n' % '\t'.join(sample_ids))
    for i in range(len(meta[:,0])):
        fout.write(u'%s\t%s\t%s\n' % (meta[i,0],meta[i,1],'\t'.join(['%6e' % x for x in data[i,:,int(write_std)]])))
    fout.close()
    
def pickleJacksFullResults( filename, jacks_results, cell_lines, gene_grnas ):    
    import pickle
    full_results = [jacks_results, cell_lines, gene_grnas]
    f = io.open(filename,'wb')
    pickle.dump(full_results, f)
    f.close()
    
def loadJacksFullResultsFromPickle( filename ):
    import pickle
    f = io.open(filename,'rb')
    jacks_results, cell_lines, gene_grnas = pickle.load(f)
    f.close()
    return jacks_results, cell_lines, gene_grnas
    
