# JACKS - joint analysis of CRISPR/Cas9 knock-out screens.
#'@import SummarizedExperiment
#'@importMethodsFrom S4Vectors metadata

library(SummarizedExperiment)

#' Infer JACKS decomposition of gRNA effect into a sample-specific gene effect, and sample-independent gRNA efficacy.
#'
#'\code{infer_jacks()} takes a \link{SummarizedExperiment}, which must pass the \link{validate_logfchange_table} check,
#'and infers the decomposition of JACKS model for all desired genes. The experiment object is endowed with
#'JACKS output: the estimated gene essentialities, their uncertainties, and directional p-values
#'(nominal and FDR corrected) are stored in experiment metadata, the gRNA efficacies with uncertainties in
#'rowData, and prediction of the decomposition and its uncertainty as assays.
#'
#'@param expt \code{\link{SummarizedExperiment}} of replicate-averaged gRNA count log2-fold changes in treatment compared to control.
#'@param target_genes A vector of gene names to run JACKS on if not all genes are required. Default NULL (all genes will be used).
#'@param n_iter An integer number of iterations of JACKS inference performed. Default 5.
#'@param reference_library Name of the gRNA library to be used in inference ("avana","gecko2","yusa_v1o", or the path to a local grna results file). If this is specified, gRNA efficacies are not estimated, which greatly speeds up the computation. Default NULL (recalculate gRNA efficacies).
#'@return \code{\link{SummarizedExperiment}} object endowed with JACKS output.
#'@export
infer_jacks <- function(expt, target_genes=NULL, n_iter=5, reference_library=NA){
    data = assays(expt)
    if (!("logfc_mean" %in% names(data))){ # calculate treatment-control comparison if not made yet
        flog.info("update_jacks: log-fold change not yet calculated, attempting to do so first")
        expt = calc_logfc(expt)
    }
    validate_logfchange_table(expt)

    if(is.null(target_genes)) {
        target_genes = unique(rowData(expt)$gene)
    }
    flog.debug(paste(length(target_genes), "target genes"))

    expt = .initialize_jacks_output(expt, reference_library)

    for (gene in target_genes){ # run JACKS for each gene
        flog.debug(paste("Processing", gene))
        I = (rowData(expt)$gene == gene)
        jacks_result = infer_jacks_gene(as.matrix(data$logfc_mean[I,]), as.matrix(data$logfc_sd[I,]), n_iter, ref_x1=rowData(expt)$jacks_x1[I], ref_x2=rowData(expt)$jacks_x2[I])
        expt = .update_with_jacks_result(expt, I, jacks_result, gene)
    }
    expt = .calc_significance(expt)
    flog.info(paste("Done infer JACKS for", length(target_genes), "genes"))
    return(expt)
}

.initialize_jacks_output <- function(expt, reference_library){
    rowData(expt)$jacks_x1 = rep(NA, dim(rowData(expt))[1]) # add JACKS output to rowdata
    rowData(expt)$jacks_x2 = rep(NA, dim(rowData(expt))[1])
    if (!is.na(reference_library)){ # if precomputed x values used from reference data:
        jacks_x = .read_precomputed_x(reference_library)
        rowData(expt)$jacks_x1 = jacks_x[as.character(rowData(expt)$gRNA),"X1"]
        rowData(expt)$jacks_x2 = jacks_x[as.character(rowData(expt)$gRNA),"X2"]
    }
    assays(expt)$prediction_sd = assays(expt)$logfc_sd*NA # add an assay that stores prediction certainty
    metadata(expt)$jacks_w = data.frame(row.names=rownames(colData(expt)))
    metadata(expt)$jacks_sdw = data.frame(row.names=rownames(colData(expt)))
    return(expt)
}

.update_with_jacks_result <- function(expt, I, result, gene){
    rowData(expt)$jacks_x[I] = result$x1
    rowData(expt)$jacks_sdx[I] = (result$x2 - result$x1**2)**0.5
    assays(expt)$prediction_sd[I,] = (result$tau)**(-0.5)
    w1 = result$w1
    sdw = (result$w2 - w1**2)**0.5
    metadata(expt)$jacks_w[[gene]] = w1
    metadata(expt)$jacks_sdw[[gene]] = sdw
    return(expt)
}

.calc_significance <- function(expt){
    m = metadata(expt)
    m$jacks_neg_pval = data.frame(row.names=rownames(colData(expt)))
    m$jacks_pos_pval = data.frame(row.names=rownames(colData(expt)))
    m$jacks_neg_fdr = data.frame(row.names=rownames(colData(expt)))
    m$jacks_pos_fdr = data.frame(row.names=rownames(colData(expt)))

    for(i in 1:dim(metadata(expt)$jacks_w)[2]){ # for each experiment
        # calculate nominal p-values from z-scores
        m$jacks_neg_pval[,i] = stats::pnorm(m$jacks_w[,i] / m$jacks_sdw[,i])
        m$jacks_pos_pval[,i] = 1. - m$jacks_neg_pval[,i]
        # calculate fdr corrected p-values
        m$jacks_neg_fdr[,i] = stats::p.adjust(m$jacks_neg_pval[,i], 'fdr')
        m$jacks_pos_fdr[,i] = stats::p.adjust(m$jacks_pos_pval[,i], 'fdr')
    }
    for(field in c("jacks_neg_pval", "jacks_pos_pval", "jacks_neg_fdr", "jacks_pos_fdr")){
        colnames(m[[field]]) = colnames(m$jacks_w)
    }
    metadata(expt) <- m
    return(expt)
}

#' Infer JACKS decomposition for a single gene
#'
#'\code{infer_jacks_gene()} takes the mean and standard deviation of gRNA log-fold change for one gene across
#'many experiments, as well as inference parameters (e.g. prior distribution specification),
#'and performs inference.
#'A list of posterior distributions is returned, one each for \emph{w}, \emph{x}, \emph{tau}. the estimated gene essentialities, their uncertainties, and directional p-values
#' are stored in experiment metadata, the gRNA efficacies with uncertainties in
#'rowData, and prediction of the decomposition and its uncertainty as assays.
#'
#'@param data \code{\link{matrix}} of replicate-averaged gRNA count log2-fold changes in treatment compared to control for all the gRNAs and experiments for one gene
#'@param data_err \code{\link{matrix}} standard deviation of the values in \code{data}.
#'@param n_iter An integer number of iterations of JACKS inference performed. Default 5.
#'@param tol Tolerance of model fit. Inference stops when lower bound from consecutive rounds does not improve by at least this much. Default 0.1
#'@param mu0_x Prior mean of gRNA efficacy \emph{x}. 1 means we assume gRNAs to work well on average. Default 1.
#'@param var0_x Prior variance of gRNA efficacy \emph{x}. 1 means we assume gRNAs to have a broad range of efficacies. Default 1.
#'@param mu0_ws Prior mean of gene essentialities \emph{w}. 0 means we assume genes to have no effect on average when knocked out. Default 0.
#'@param var0_w Prior variance of gene essentiality \emph{w}. Large value means we do not constrain the inferred values. Default 10000.
#'@param tau_prior_strength Strength of prior on \emph{tau}. 0.5 means we weight information about uncertainty equally from measurement noise and accuracy of the JACKS decomposition. Default 0.5.
#'@param ref_x1 Reference gRNA efficacy E(\emph{x}). If provided, will be used directly, instead of inferring from data. Default NA.
#'@param ref_x2 Reference gRNA efficacy E(\emph{x}^2). If provided, will be used directly, instead of inferring from data. Default NA.
#'@return \code{List} object endowed with inferred posterior distributions for one gene.
#'@export
infer_jacks_gene <- function(data, data_err, n_iter=5, tol=0.1, mu0_x=1, var0_x=1.0, mu0_ws=0, var0_w=1e4, tau_prior_strength=0.5, ref_x1=NA, ref_x2=NA){
    model_state = .prep_infer_jacks_gene(data, data_err, mu0_x, tau_prior_strength, ref_x1, ref_x2)

    #Run the inference
    model_state = .update_lower_bound(model_state)
    flog.debug(paste("infer_jacks_gene: After init, mean absolute error=",model_state$err,"lower bound=", model_state$lower_bound))

    for (i in 1:n_iter){
        last_bound = model_state$lower_bound
        model_state = .update_x(model_state, mu0_x, var0_x, ref_x1, ref_x2)
        model_state = .update_w(model_state, mu0_ws, var0_w)
        model_state = .update_tau(model_state, tau_prior_strength)
        model_state = .update_lower_bound(model_state)
        flog.debug(paste0("infer_jacks_gene: Iteration ",i, " bound=", model_state$lower_bound))
        if (abs(last_bound - model_state$lower_bound) < tol){
            break
        }
    }
    return(model_state)
}

.prep_infer_jacks_gene <- function(data, data_err, mu0_x, tau_prior_strength, ref_x1=NA, ref_x2=NA){
    tau_prior_b = tau_prior_strength*1.0*(data_err**2+1e-2)
    state = list(y=as.matrix(data), w1=as.vector(apply(data, 2, stats::median)), x1=mu0_x*rep(1, dim(data)[1]), taua=tau_prior_strength*1.0, taub=tau_prior_b)
    state$x2 = state$x1**2
    state$w2 = state$w1**2
    state$tau = state$taua/state$taub
    state$tau_prior_b = tau_prior_b
    if(!is.na(ref_x1) && !is.na(ref_x2)){
        state$x1 = ref_x1
        state$x2 = ref_x2
    }
    return(state)
}

.update_x <- function(state, mu0_x, var0_x, ref_x1, ref_x2){
    if(!is.na(ref_x1) && !is.na(ref_x2)){ # reference library efficacies are given, and informative
        flog.debug("    Using reference X; skipping update")
        return(state)
    }
    state$x1 = as.vector((mu0_x/var0_x + ((state$y*state$tau) %*% state$w1))/((state$tau %*% state$w2) + 1.0/var0_x))
    state$x2 = as.vector(state$x1**2 + 1.0/((state$tau %*% state$w2) + 1.0/var0_x))
    wadj = 0.5/length(state$x1)
    x1m = mean(state$x1) + 2*wadj*stats::median(state$x1) - wadj*max(state$x1) - wadj*min(state$x1)
    state$x1 = state$x1/x1m
    state$x2 = state$x2/x1m/x1m
    flog.debug(paste0("    X update: <x>=",mean(state$x1)))
    return(state)
}

.update_w <- function(state, mu0_ws, var0_w){
    state$w1 = as.vector((mu0_ws/var0_w + (state$x1 %*% (state$y*state$tau))/((state$x2 %*% state$tau) + 1.0/var0_w)))
    state$w2 = as.vector(state$w1**2 + 1.0/((state$x2 %*% state$tau) +1.0/var0_w))
    flog.debug(paste0("    W update: <w>=",mean(state$w1)))
    return(state)
}

.update_tau <- function(state, tau_prior_strength){
    b_star = state$y**2 - 2*state$y*outer(state$x1,state$w1) + outer(state$x2,state$w2)
    state$taua = tau_prior_strength + 0.5
    state$taub = state$tau_prior_b + 0.5*b_star
    state$tau = state$taua/state$taub
    flog.debug(paste0("    Tau update: <t>=",mean(state$tau)))
    return(state)
}

.update_lower_bound <- function(state){
    xw = outer(state$x1, state$w1)
    state$lower_bound = sum(state$tau*(state$y**2 + outer(state$x2, state$w2) -2*xw*state$y))
    #ent = .entropy_normal(state$x1, state$x2) + .entropy_normal(state$w1, state$w2) + .entropy_gamma(state$taua, state$taub)
    #loglik = .logexp_gamma(state$taua, state$taub)
    #loglik = loglik - sum(state$tau*(state$y**2 + outer(state$x2, state$w2) -2*xw*state$y))
    #state$lower_bound = ent + 0.5*loglik
    state$err = mean(abs(state$y-xw))
    return(state)
}

.entropy_gamma <- function(a,b){
    sum(a - log(b) + lgamma(a) + (1-a)*digamma(a))
}

.logexp_gamma <- function(a,b){
    sum(digamma(a) - log(b))
}

.entropy_normal <- function(m1, m2){
    sum(0.5*log(m2-m1^2))
}

.vecprod <- function(m, v){apply(t(m), 1, function(x) x*v)}

#' Calculate replicate-averaged log2-fold changes in gRNA frequency.
#'
#'\code{calc_logfc()} takes a \link{SummarizedExperiment}, which must pass either the \link{validate_logf_table}
#'or \link{validate_count_table} check, and calculates the replicate-averaged log2-fold change of gRNA frequency
#'in treatment condition compared to control. The experiment object is endowed with fields \code{logfc_mean} and
#'\code{logfc_sd} accordingly. If counts are provided, the frequency only averages \code{logf_mean} and \code{logf_sd}
#'are calculated and stored first as well.
#'
#'@param expt \code{\link{SummarizedExperiment}} of raw read counts, or replicate-averaged log2-scale gRNA counts
#'@param reference String corresponding to factor level in colData$Condition that specifies the change baseline. Default "CTRL".
#'@param condition String corresponding to factor level in colData$Condition that specifies the condition for which change is calculated. Default "D21".
#'@return \code{List} object endowed with inferred posterior distributions for one gene.
#'@export
calc_logfc <- function(expt, reference="CTRL", condition="D21"){
    if (!("logf_mean" %in% names(assays(expt)))){ # attempt to calculate log-scale frequencies if not done yet
        flog.info("calc_logfc: log-scale counts not yet calculated, attempting to do so first from assumed counts")
        expt = .collapse_replicates(expt)
    }
    validate_logf_table(expt)
    sample_data = colData(expt)

    if (sum(sample_data$Condition == reference) == 1){ # common reference
        flog.info("Calculating log2-fold changes. Common reference")
    } else{
        flog.info("Calculating log2-fold changes. Individual references")
    }

    meanmat = data.frame(row.names=rownames(assays(expt)$logf_mean))
    sdmat = data.frame(row.names=rownames(assays(expt)$logf_sd))
    good_samples = c()
    for (s in unique(sample_data$Name)){
        flog.debug(paste("Processing sample", s))
        I_cond = (sample_data$Name == s) & (sample_data$Condition == condition)
        if (sum(sample_data$Condition == reference) == 1){ # common reference - use a single sample as reference
            I_ref = sample_data$Condition == reference
        } else{ # per-sample reference
            I_ref = (sample_data$Name == s) & (sample_data$Condition == reference)
        }

        if((sum(I_ref) == 1) && (sum(I_cond) == 1)){ # after collapsing replicates, exactly one reference and sample
            j_ref = which(I_ref)
            j_cond = which(I_cond)
            meanmat[[s]] = assays(expt)$logf_mean[,j_cond] - assays(expt)$logf_mean[,j_ref]
            sdmat[[s]] = (assays(expt)$logf_sd[,j_cond]**2 + assays(expt)$logf_sd[,j_ref]**2)**0.5
            good_samples = c(good_samples, s)
        } else{
            flog.warn(paste("Sample",s,":", sum(I_ref), "reference samples and", sum(I_cond), "condition samples - not precisely one of each as expected. Ignore if this is the control sample, check sample assignment otherwise"))
        }
    }
    result = SummarizedExperiment(assays=list("logfc_mean"=meanmat, "logfc_sd"=sdmat), rowData=rowData(expt), colData=data.frame("Name"=good_samples, row.names=good_samples), metadata=metadata(expt))
    validate_logfchange_table(result)
    flog.info(paste("Done calculating mean log2 fold changes;", length(names), "comparisons"))
    return(result)
}

.collapse_replicates <- function(expt, count_prior=16., window=800) {
    validate_count_table(expt)
    samples = colData(expt)
    counts = assays(expt)$counts

    # create new matrix of log-scale centered and counts
    meanmat = matrix(ncol=0, nrow=dim(counts)[1], dimnames=list(rownames(counts), NULL))
    sdmat = matrix(ncol=0, nrow=dim(counts)[1], dimnames=list(rownames(counts), NULL))
    condensed_samples = c()

    for (i in 1:dim(samples)[1]){ # for each sample
        I = (samples$Name == samples$Name[i]) & (samples$Condition == samples$Condition[i])
        condensed_name = paste0(samples$Name[i], "_", samples$Condition[i])
        if (!(condensed_name %in% condensed_samples)){
            flog.debug(paste("Processing", condensed_name))
            condensed_samples = c(condensed_samples, condensed_name)
            d = .calc_mean_sd(counts[,I], count_prior, window)
            meanmat = cbind(meanmat, as.vector(d$mean))
            sdmat = cbind(sdmat, as.vector(d$sd))
        }
    }
    colnames(meanmat) = condensed_samples
    colnames(sdmat) = condensed_samples
    condensed_meta = t(as.data.frame(strsplit(condensed_samples,"_"), row.names=c("Name","Condition"), col.names=condensed_samples, check.names=FALSE))
    result = SummarizedExperiment(assays=list("logf_mean"=meanmat, "logf_sd"=sdmat), rowData=rowData(expt), colData=condensed_meta, metadata=metadata(expt))
    validate_logf_table(result)
    flog.info(paste("Done condensing replicates; went from", dim(samples)[1], "experiments down to", length(condensed_samples), "samples"))
    return(result)
}

.calc_mean_sd <- function(counts, count_prior=16., window=800) {
    log2c = log2(counts + count_prior)
    log2cm = apply(log2c, 2, function(x) x - stats::median(x))
    result = list(mean=apply(log2cm, 1, mean))
    result$sd = .calc_posterior_sd(log2cm, window=window, do_monotonize=TRUE, window_estimate=TRUE)
    result$sd[is.na(result$sd)] = 2.0 # very uncertain if a single replicate

    return(result)
}

#""" Simple posterior estimation - use the average variance from 2*window (default 50) gRNAs with nearest means as the estimate if larger than observed. """
.calc_posterior_sd <- function(data, window=800, do_monotonize=TRUE, window_estimate=TRUE) {
    flog.debug("    Calculating posterior SD")
    N = dim(data)[1]
    R = dim(data)[2]
    I = order(apply(data, 1, mean),(apply(data, 1, stats::sd)**2)*(R-1)/R)
    vars = (apply(data[I,], 1, stats::sd)**2)*(R-1)/R # observed (biased) variances sorted by corresponding means (NOte: biased for compatibility with Python, since numpy.nanstd produces a biased estimator )
    # window smoothing of variances
    m = .window_smooth(vars, window=window)

    # enforce monotonicity (decreasing variance with log count) if desired
    if (do_monotonize) { m = .monotonize(m) }

    # construct the posterior estimate
    sd_post = rep(0, N)
    for (i in 1:length(I)) {
        idx = I[i]
        sd_post[idx] = max(m[i], vars[i])**0.5 # take the larger of window estimate and actual observation; right hand side all [i], since both m and s2s are sorted
        if (window_estimate) { sd_post[idx] = m[i]**0.5 } # ignore if much larger
    }
    return(sd_post)
}

#""" Apply a mean averaging window across a sample with a shift of 1"""
.window_smooth <- function(x, window=25){
    flog.debug("    Smoothing")
    N = length(x)
    m = rep(0, N) # mean variance in window
    m[1:(window+1)] = mean(x[1:(2*window+1)]) # left shoulder; assume constant start to avoid high variance
    m[(N-window):N] = mean(x[(N-2*window):N]) # right shoulder; assume constant start to avoid high variance
    for (k in (window+2):(N-window-1)) { m[k] = m[k-1] + (x[k+window+1] - x[k-window])/(2*window + 1) } # linear time average
    return(m)
}

# Monotonize a vector so it is strictly decreasing
.monotonize <- function(x){
    flog.debug("    Monotonizing")
    N = length(x)
    for (i in 1:N){
        if (!is.nan(x[N-i+1])){
            x[N-i] = max(x[N-i+1],x[N-i])
        }
    }
    return(x)
}

