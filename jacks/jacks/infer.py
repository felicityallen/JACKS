import numpy as np
import random
import scipy as SP
import logging

logging.basicConfig(level=logging.DEBUG, format='[%(asctime)s] %(name)-5s: %(levelname)-8s %(message)s')
LOG = logging.getLogger("jacks")
LOG.setLevel(logging.DEBUG)

""" Main function from which to compute JACKS across a set of specified genes
@param gene_index {gene: [list of row indexes for the gRNAs for that gene in the testdata]}
@param testdata (num_total_grnas x num_conditions x 2) numpy array containing the mean (idx 0) and estimated variance (idx 1 ) 
         of each median normalized log-fold change
@param ctrldata Either (num_total_grnas x num_conditions x 2) numpy array containing the respective control measurements
         for each element of testdata OR a (num_total_grnas x 2) numpy array containing a common control measurement to
         be used for all conditions.
@return {gene: y, tau, x1, x2, w1, w2}, the JACKS inference results for each gene in the input gene_index""" 
def inferJACKS(gene_index, testdata, ctrldata, fixed_x=None, n_iter=50, apply_w_hp=False):
    results = {}
    for gene in gene_index:
        Ig = gene_index[gene]
        if fixed_x is not None:
            gene_fixed_x = {'X1':fixed_x['X1'][Ig], 'X2':fixed_x['X2'][Ig]}
        else: gene_fixed_x = fixed_x
        if testdata.shape == ctrldata.shape: # each line has a matching control:
            results[gene] = inferJACKSGene(testdata[Ig,:,0], testdata[Ig,:,1], ctrldata[Ig,:,0], ctrldata[Ig,:,1], n_iter, fixed_x=gene_fixed_x, apply_w_hp=apply_w_hp)
        else: # shared control
            results[gene] = inferJACKSGene(testdata[Ig,:,0], testdata[Ig,:,1], ctrldata[Ig,0], ctrldata[Ig,1], n_iter, fixed_x=gene_fixed_x, apply_w_hp=apply_w_hp)

    return results

""" Convenience function for matrix-vector and vector-vector dot products, ignoring Nans
@param x1 matrix or vector 1
@param x2 matrix or vector 2
@return x1 dot x2, where the summands that are NaN are ignored """ 
def nandot(x1, x2):
    if len(x1.shape) == 1 and len(x2.shape) == 2:
        x1T = SP.tile(x1, [x2.shape[1],1]).transpose()
        return SP.nansum(SP.multiply(x1T,x2), axis=0)
    elif len(x2.shape) == 1 and len(x1.shape) == 2:
        x2T = SP.tile(x2, [x1.shape[0],1])
        return SP.nansum(SP.multiply(x1,x2T), axis=1)
    elif len(x1.shape) == 1 and len(x2.shape) == 1:
        return SP.nansum(SP.multiply(x1,x2))      
    return None
 
    
""" Run JACKS inference on a single gene
@param data: G x L matrix of guide effects (log2-scale difference of cell line frequency to control, assume correctly normalized); one value for each guide and line
@param yerr: G x L matrix of guide effect standard deviations
@param n_iter: number of iterations of the EM [5]
@param verbose: whether to output debug information [False]
@return length-G vector X of guide efficacies, length-L vector W of cell line effects, GxL posterior variance estimate of Y, GxL reconstruction error (Y - X x W)
"""
def inferJACKSGene(data, data_err, ctrl, ctrl_err, n_iter, tol=0.1, mu0_x=1, var0_x=1.0, mu0_w=0.0, var0_w=1e4, tau_prior_strength=0.5, fixed_x=None, apply_w_hp = False):

    #Adjust estimated variances if needed
    data_err[SP.isnan(data_err)] = 2.0 # very uncertain if a single replicate

    #The control can be specified once for each sample, or common across all cell lines
    if len(ctrl.shape) == 1 or ctrl.shape[1] == 1: 
        #If only 1 control replicate, use mean variance from data across cell lines for that guide 
        ctrl_err[SP.isnan(ctrl_err)] = SP.nanmean(data_err, axis=1)[SP.isnan(ctrl_err)]
        y = (data.T - ctrl).T
        tau_pr_den = tau_prior_strength*1.0*((data_err**2).T + ctrl_err**2 + 1e-2).T
    else:
        #If only 1 control replicate, use data variances for ctrls as well
        ctrl_err[SP.isnan(ctrl_err)] = data_err[SP.isnan(ctrl_err)]
        y = data - ctrl
        tau_pr_den = tau_prior_strength*1.0*(data_err**2 + ctrl_err**2 + 1e-2)

    #Run the inference
    G,L = y.shape
    if fixed_x is None:
        x1 = mu0_x*SP.ones(G)
        x2 = x1**2
    else:
        x1 = fixed_x['X1']
        x2 = fixed_x['X2']

    w1 = np.nanmedian(y, axis=0)
    
    tau = tau_prior_strength*1.0/tau_pr_den
   
    w2 = w1**2
    bound = lowerBound(x1,x2,w1,w2,y,tau)
    LOG.debug("Initially, mean absolute error=%.3f"%(SP.nanmean(abs(y)).mean()))
    LOG.debug("After init, mean absolute error=%.3f, <x>=%.1f <w>=%.1f lower bound=%.1f"%(SP.nanmean(abs(y.T-SP.outer(w1,x1))).mean(), x1.mean(), w1.mean(), bound))
    for i in range(n_iter):
        last_bound = bound
        if fixed_x is None: x1,x2 = updateX(w1,w2,tau,y,mu0_x,var0_x)
        if apply_w_hp and len(w1) > 1: mu0_w, var0_w = w1.mean(), w1.var()*3+1e-4  # hierarchical update on w (to encourage w's together - use with caution!)
        w1,w2 = updateW(x1,x2,tau,y,mu0_w,var0_w)
        tau = updateTau(x1, x2, w1, w2, y, tau_prior_strength, tau_pr_den)
        bound = lowerBound(x1,x2,w1,w2,y,tau)
        LOG.debug("Iter %d/%d. lb: %.1f err: %.3f x:%.2f+-%.2f w:%.2f+-%.2f xw:%.2f"%(i+1, n_iter, bound, SP.nanmean(abs(y.T-SP.outer(w1,x1))).mean(), x1.mean(), SP.median((x2-x1**2)**0.5), w1.mean(), SP.median((w2-w1**2)**0.5), x1.mean()*w1.mean()))
        if abs(last_bound - bound) < tol:
            break
    return y, tau, x1, x2, w1, w2
    
""" Convenience functions following the latex code for variational updates.
"""
def updateX(w1, w2, tau, y, mu0_x, var0_x):
    x1 = (mu0_x/var0_x + nandot((y.T).T*tau,w1))/(nandot(tau,w2) + 1.0/var0_x)
    x2 = x1**2 + 1.0/(nandot(tau,w2)+1.0/var0_x)
    wadj = 0.5/len(x1)
    #Normalize by the median-emphasized mean of x
    x1m = x1.mean() + 2*wadj*np.nanmedian(x1) - wadj*x1.max() - wadj*x1.min()      
    LOG.debug("After X update, <x>=%.1f, mean absolute error=%.3f"%(x1.mean(), SP.nanmean(abs(y.T-SP.outer(w1,x1))).mean()))
    return x1/x1m, x2/x1m/x1m

def updateW(x1, x2, tau, y, mu0_ws, var0_w):
    w1 = (mu0_ws/var0_w + nandot(x1,(y.T).T*tau))/(nandot(x2,tau)+1.0/var0_w)
    w2 = w1**2 + 1.0/(nandot(x2,tau)+1.0/var0_w)
    LOG.debug("After W update, <w>=%.1f, mean absolute error=%.3f"%(w1.mean(), SP.nanmean(abs(y.T-SP.outer(w1,x1))).mean()))
    return w1, w2

def updateTau(x1, x2, w1, w2, y, tau_prior_strength, tau_pr_den):
    b_star = y**2 - 2*y*(SP.outer(x1,w1)) +SP.outer(x2,w2) 
    tau = (tau_prior_strength + 0.5)/(tau_pr_den + 0.5*b_star)   
    return tau

def lowerBound(x1,x2,w1,w2,y,tau):
    xw = SP.outer(x1,w1)
    return SP.nansum(tau*(y**2 + SP.outer(x2,w2) -2*xw*y))
    

    
