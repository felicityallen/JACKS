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
def calc_posterior_sd(data, windowfracN=100, do_monotonize=True, window_estimate=True, guideset_indexs=set()):

    MIN_WINDOW, MAX_WINDOW = 30, 800
    N = data.shape[0]
    
    #If there's only one replicate, set variances undefined (nans)
    if len(data.shape) == 1 or data.shape[1] == 1:
        return np.zeros(N)*np.nan

    #If there's only one 
    if len(guideset_indexs) == 0: guideset_indexs = set([x for x in range(N)])
    window = min(max( int(len(guideset_indexs)/windowfracN), MIN_WINDOW ), MAX_WINDOW)

    #Sort by means, then (in case of equality) variances
    dmean_vars = [(x,y,i,(i in guideset_indexs)) for (x,y,i) in zip(data.mean(axis=1), np.nanstd(data, axis=1)**2, range(data.shape[0]))]
    dmean_vars.sort()
    dmean_vars = np.array(dmean_vars)

    # window smoothing of variances
    guideset_mask = (dmean_vars[:,3] == True)
    m = window_smooth(dmean_vars[guideset_mask,1], window)
    # enforce monotonicity (decreasing variance with log count) if desired
    if do_monotonize: m = monotonize(m)
    if len(guideset_indexs) != N:
        #Interpolate back up to full guide set
        m = np.interp(dmean_vars[:,0], dmean_vars[guideset_mask,0], m)

    # construct the posterior estimate and restore original ordering
    sd_post = np.zeros(N)
    if not window_estimate: m = np.maximum(m, dmean_vars[:,1])
    sd_post[dmean_vars[:,2].astype(int)] = (m**0.5)

    return sd_post  
    
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

def inferMissingVariances(data, meta, sample_ids, ctrl_spec, ctrl_guideset):
    
    #Check for nan variances
    for s, sample_id in enumerate(sample_ids):
        if sum(np.isnan(data[:,s,1]) == 0): continue    
        
        #If this is a control replicate, ignore it, as JACKS will use the data variance for this at model time
        if sample_id in ctrl_spec and ctrl_spec[sample_id] == sample_id: continue

        # If this is a sample, and a ctrl_guideset is specified, use twice the variance between
        # sample and control within this set to infer the mean-variance relationship and apply for others
        if sample_id in ctrl_spec and len(ctrl_guideset) == 0:
            guideset_indexs = [meta[:,0].index(x) for x in ctrl_guideset]
            ctrl_data = data[:,sample_ids.index(ctrl_spec[sample_id]), 1]
            concat_data = np.concatenate((ctrl_data, data[:,s,1]), axis=1)
            data[:,s,1] = 2*calc_posterior_sd(concat_data, guideset_indexs) #sigma_hat

    return data

def collateTestControlSamples(data, sample_ids, ctrl_spec):
    test_sample_idxs = [i for i,x in enumerate(sample_ids) if ctrl_spec[x] != x]
    testdata = data[:,test_sample_idxs,:]
    ctrldata = data[:,[sample_ids.index(ctrl_spec[sample_ids[idx]]) for idx in test_sample_idxs],:]    
    return testdata, ctrldata, test_sample_idxs

""" Load counts data from a list of one or more files with the same guides (in the same order)

@param sample_spec is {filename:[(sample_id, colname)]} i.e. for each file a list of column names you
   want to load and their sample_ids (to group replicates)"""
def loadDataAndPreprocess(sample_spec, gene_spec, ctrl_spec={}, ctrl_guideset=set(), normtype='median', prior=32):

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

    #Infer missing variances (for any cases of single replicates)
    data = inferMissingVariances(data, meta, sample_ids, ctrl_spec, ctrl_guideset)
            
    return data, meta, sample_ids, genes, gene_guides

""" Select a specified number of replicates for each of a selection of screens from a larger file and generate the condensed data,
    the output is the same as for read_condensed in io.py
@param selected_screens A list of tuples, each containing a screen identifier string (cell line) and the number of replicates 
       to select or -1 to use all replicates available without replacement
@param norepl Whether to randomly select samples for each screen with (False) or without (True) replacement
@param All other params as for loadDataAndPreprocess""" 

def subsample_and_preprocess(selected_screens, sample_spec, gene_spec, ctrl_spec={}, ctrl_guideset=set(), norepl=True):
    
    #Collect replicates for relevant samples
    found_replicates = {}
    selected_lookup = {sample_id: num_rep for (sample_id, num_rep) in selected_screens}
    for filename in sample_spec:
        for sample_id, colname in sample_spec[filename]:
            if sample_id not in selected_lookup: continue
            if sample_id not in found_replicates:
                found_replicates[sample_id] = []
            found_replicates[sample_id].append((filename, colname))
    
    #Sub-sample replicates
    selected_sample_spec = {}
    sample_count = {}
    for sample_id, num_rep in selected_screens:
        if sample_id not in found_replicates:
            raise Exception('Could not find data for selected sample %s' % sample_id)
        if num_rep < 0:
            selected_replicates = found_replicates[sample_id]
        elif norepl:
            if num_rep <= len(found_replicates[sample_id]):
                selected_replicates = random.sample(found_replicates[sample_id], num_rep)
            else: selected_replicates = found_replicates[sample_id]
        else:
            selected_replicates = [random.choice(found_replicates[sample_id]) for i in range(num_rep)]

        ctrl_sample_id = ctrl_spec[sample_id]
        if sample_id in sample_count:
            #If using the same sample more than once, give the extra ones a different sample id
            sample_count[sample_id] += 1 
            sample_id += ('_%d' % sample_count[sample_id])
            ctrl_spec[sample_id] = ctrl_sample_id
        elif ctrl_sample_id != sample_id:
            sample_count[sample_id] = 0

        for filename, colname in selected_replicates:
            if filename not in selected_sample_spec:
                selected_sample_spec[filename] = []
            selected_sample_spec[filename].append((sample_id, colname))

    #Load and preprocess data for those replicates
    data, meta, cell_lines, genes, gene_guides = loadDataAndPreprocess(selected_sample_spec, gene_spec, ctrl_spec=ctrl_spec, ctrl_guideset=ctrl_guideset)
    return data, meta, cell_lines, genes, gene_guides


