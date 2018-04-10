#'@importFrom methods is
#'@importFrom utils read.table

.validate_sample_names <- function(names){
    for (n in names){
        if(nchar(n) - nchar(gsub("_","",n)) != 2) {
            stop(paste0("Sample '", n, " 'does not conform to expected form <sample-name>_<condition>_<replicate>. Please make sure exactly two underscores '_' are present that separate these fields."))
        }
    }
    flog.debug("All sample names successfully validated")
}

#' Validate count table
#'
#'\code{validate_count_table()} takes in an object, and confirms whether it meets the
#'requirements to be used in JACKS. A valid count table is of type \code{\link{SummarizedExperiment}},
#'which includes an assay named \code{'counts'} that contains the raw read counts.
#'To enable gene-level analysis, the \code{rowData} has to contain a column named \code{'gene'}
#'that will be used to collage gRNA information. \code{ColData} of the experiment should
#'be a valid sample table
#'
#'@param counts \code{\link{SummarizedExperiment}} object
#'@return TRUE if JACKS can be run on this object; exception otherwise
#'@export
validate_count_table <- function(counts){
    if(!methods::is(counts, "SummarizedExperiment")) { # check object types
        stop("Counts object is not of type SummarizedExperiment")
    }
    if(!("counts" %in% names(assays(counts)))) {
        stop("Assays of counts object does not include actual counts ('counts' not in assays(obj)). Make sure capitalization is correct.")
    }
    if(!("gene" %in% colnames(rowData(counts)))) {
        stop("gRNA metadata does not include gene names ('gene' not in columns of rowData(obj)). Make sure capitalization is correct")
    }
    validate_sample_table(colData(counts))
    flog.debug("Count table successfully validated")
    return(TRUE)
}

#' Validate sample table
#'
#'\code{validate_sample_table()} takes in an object, and confirms whether it meets the
#'requirements to be used in JACKS. A valid sample table has one row per sample, and columns
#'according to the annotation. Depending on strictness of definition, a valid table must have
#'metadata columns for \code{'Name'} (with capital N; sample name), \code{'Condition'} (condition the sample was
#'measured in, e.g. \code{"CTRL"}, or \code{"Day21-Puromycin"}), and \code{'Replicate'} (replicate index, e.g. \code{'1'},
#' \code{'Rep1'}, or \code{'A'}). Samples with same name and condition will be collated across all unique replicates;
#' different conditions can later be specified to be contrasted in JACKS analysis. If samples are adequately
#' described, they can be automatically analysed from raw counts; if not, the user has to manually
#' calculate log-fold changes before running JACKS.
#'
#'@param samples \link{data.frame} of samples
#'@param require_name A logical scalar whether 'Name' field is required. Default: TRUE
#'@param require_cond A logical scalar whether 'Condition' field is required. Default: TRUE
#'@param require_rep A logical scalar whether 'Replicate' field is required. Default: TRUE
#'@return TRUE if JACKS can be run using this object as sample metadata; exception otherwise
#'@export
validate_sample_table <- function(samples, require_name=TRUE, require_cond=TRUE, require_rep=TRUE){
    c = colnames(samples)
    if(require_name && !("Name" %in% c)){
        stop("Sample metadata does not include sample names ('Name' column). Make sure capitalization is correct")
    }
    if(require_cond && !("Condition" %in% c)){
        stop("Sample metadata does not include treatment condition ('Condition' column). Make sure capitalization is correct")
    }
    if(require_rep && !("Replicate" %in% c)){
        stop("Sample metadata does not include replicate info ('Replicate' column). Make sure capitalization is correct")
    }
    if(require_cond){
        # conditions exist - make sure each sample has exactly one matching control
    }
    return(TRUE)
}

#' Validate log-scale count table
#'
#'\code{validate_logf_table()} takes in an object, and confirms whether it meets the
#'requirements to be used in JACKS. A valid log scale count table is of type SummarizedExperiment,
#'which includes an assay named 'logf_mean' that contains the log-scale read counts, normalised
#'and averaged across replicates, as well as an assay 'logf_sd', which contains the standard deviation
#'of the log-scale counts across replicates. To enable gene-level analysis, the rowData has to contain
#'a column named 'gene' that will be used to collage gRNA information. ColData of the experiment should
#'be a valid sample table. Finally, if desired, it is checked that log-fold changes have zero median.
#'
#'@param logf \code{\link{SummarizedExperiment}} object
#'@param require_0median Whether to require per-sample frequencies to have median of 0. Default TRUE.
#'@return TRUE if JACKS can be run on this object; exception otherwise
#'@export
validate_logf_table <- function(logf, require_0median=FALSE){
    if(!methods::is(logf, "SummarizedExperiment")) { # check object types
        stop("log2 gRNA frequency object is not of type SummarizedExperiment")
    }
    if(!("logf_mean" %in% names(assays(logf))) || !("logf_sd" %in% names(assays(logf)))) {
        stop("Assays of log2 gRNA frequency object do not include means or standard deviations ('logf_mean' or 'logf_sd not in assays(obj)). Make sure capitalization is correct.")
    }
    if(!("gene" %in% colnames(rowData(logf)))) {
        stop("gRNA metadata does not include gene names ('gene' not in columns of rowData(obj)). Make sure capitalization is correct")
    }
    validate_sample_table(colData(logf), require_rep=FALSE)
    if(require_0median){
        logf_mean = assays(logf)$logf_mean
        logf_sd = assays(logf)$logf_sd
        for (i in 1:dim(logf_mean)[2]){
            if(abs(stats::median(logf_mean[,i])) > 1e-8){
                stop(paste0("Sample '", colnames(logf_mean)[i], "' has non-zero median"))
            }
        }
    }
    flog.debug("log2 gRNA frequency table successfully validated")
    return(TRUE)
}

#' Validate log-fold change table
#'
#'\code{validate_logfchange_table()} takes in an object, and confirms whether it meets the
#'requirements to be used in JACKS. A valid log-fold change table is of type SummarizedExperiment,
#'which includes an assay named 'logfc_mean' that contains the log-fold change, combined across replicates,
#'as well as an assay 'logfc_sd', which contains the standard deviation of this estimate.
#'To enable gene-level analysis, the rowData has to contain a column named 'gene'
#'that will be used to collage gRNA information. ColData of the experiment should be a valid sample table.
#'
#'@param logfc \code{\link{SummarizedExperiment}} object
#'@return TRUE if JACKS can be run on this object; exception otherwise
#'@export
validate_logfchange_table <- function(logfc){
    if(!methods::is(logfc, "SummarizedExperiment")) { # check object types
        stop("log2 fold change object is not of type SummarizedExperiment")
    }
    if(!("logfc_mean" %in% names(assays(logfc))) || !("logfc_sd" %in% names(assays(logfc)))) {
        stop("Assays of log2 fold change object do not include means or standard deviations ('logfc_mean' or 'logfc_sd not in assays(obj)). Make sure capitalization is correct.")
    }
    if(!("gene" %in% colnames(rowData(logfc)))) {
        stop("gRNA metadata does not include gene names ('gene' not in columns of rowData(obj)). Make sure capitalization is correct")
    }
    validate_sample_table(colData(logfc), require_cond=FALSE, require_rep=FALSE)
    flog.debug("log2 fold change table successfully validated")
    return(TRUE)
}

.test_inference <- function(){
    N = 5
    M = 30
    data = array(stats::rnorm(N*M,mean=0,sd=1), c(N, M))
    data_var = matrix(stats::rnorm(N*M,mean=0,sd=1), N, M)**2
    ctrl = matrix(stats::rnorm(N*M,mean=0,sd=1), N, M)
    ctrl_var = matrix(stats::rnorm(N*M,mean=0,sd=1), N, M)**2
    infer_jacks_gene(data-ctrl, (data_var+ctrl_var)^0.5, n_iter=10)
}

.test_realdata_inference <- function(test_genes=c("KRAS", "RRM2", "ZNF253")){
    count_file = system.file("extdata", "example_count_data.tab", package="jacks", mustWork=TRUE)
    sample_spec_file = system.file("extdata", "example_repmap.tab", package="jacks", mustWork=TRUE)
    lfc = read_counts_from_spec_files(count_file, sample_spec_file, replicate_col="Replicate", sample_col="Sample", gene_spec_file=count_file, grna_col="sgRNA", gene_col="gene", count_prior=32., normalization='median', window=800, reference_sample="CTRL")
    for(gene in test_genes){
        rjacks =  infer_jacks(lfc, gene)
        plot_jacks(rjacks, gene, do_save=FALSE)
    }
}

.test_samplespec_io <- function(){
    f = system.file("extdata", "avana_head.tab", package="jacks", mustWork=TRUE)
    grna_meta = utils::read.table(f, header=TRUE, stringsAsFactors=FALSE)[,1:2]
    sample_spec = rbind(c(f, "Plasmid_CTRL_1", "ctrl"), c(f, "Plasmid_CTRL_2", "ctrl"), c(f, "Plasmid_CTRL_3", "ctrl"), c(f, "KP-3_D21_1", "KP3"), c(f, "KP-3_D21_2","KP3"))
    read_count_table_from_spec(sample_spec, grna_meta)
}

.test_pyjacks_concordance <- function(mean_tol=1e-2, std_tol=1e-2){
    data = NULL
    data_err = NULL
    pyvals = NULL

    count_file = system.file("extdata", "example_count_data.tab", package="jacks", mustWork=TRUE)
    sample_spec_file = system.file("extdata", "example_repmap.tab", package="jacks", mustWork=TRUE)
    lfc = read_counts_from_spec_files(count_file, sample_spec_file, replicate_col="Replicate", sample_col="Sample", gene_spec_file=count_file, grna_col="sgRNA", gene_col="gene", count_prior=32., normalization='median', window=800, reference_sample="CTRL")
    rjacks = infer_jacks(lfc, c("KRAS"))
    
    #To produce these run 
    pyjacks_genefile = system.file("extdata", "pyjacks_for_R_comparison_gene_JACKS_results.txt", package="jacks", mustWork=TRUE)
    pyjacks_genestdfile = system.file("extdata", "pyjacks_for_R_comparison_genestd_JACKS_results.txt", package="jacks", mustWork=TRUE)
    pyjacks_grnafile = system.file("extdata", "pyjacks_for_R_comparison_grna_JACKS_results.txt", package="jacks", mustWork=TRUE)

    py_wmean <- read.table(file=pyjacks_genefile, header=TRUE)
    py_wstd <- read.table(file=pyjacks_genestdfile, header=TRUE)
    py_x <- read.table(file=pyjacks_grnafile, header=TRUE)    
    
    failed = FALSE
    for(i in 1:length(colnames(rjacks))){
      cline <- colnames(rjacks)[i]
      r_w <- metadata(rjacks)$jacks_w$KRAS[i]
      if( abs(r_w - py_wmean[cline]) > mean_tol ){
        flog.info(paste(cline, "\tpyjacks:", toString(py_wmean[cline]), "\tRjacks:", toString(r_w)))
        failed = TRUE
      }
    }
    for(i in 1:length(colnames(rjacks))){
      cline <- colnames(rjacks)[i]
      r_wstd <- metadata(rjacks)$jacks_sdw$KRAS[i]    
      if( abs(r_wstd - py_wstd[cline]) > std_tol ){
        flog.info(paste(cline, "\tpyjacks:", toString(py_wstd[cline]), "\tRjacks:", toString(r_wstd)))
        failed = TRUE
      }
    }
    if(failed){ stop("Gene KRAS inference not concordant between Rjacks and pyjacks") }
    flog.info(paste("pyjacks concordant with Rjacks on gene KRAS"))
}

.test_pyjacks_reference_concordance <- function(mean_tol=1e-2, std_tol=1e-2){
  data = NULL
  data_err = NULL
  pyvals = NULL
  
  count_file = system.file("extdata", "example_count_data.tab", package="jacks", mustWork=TRUE)
  sample_spec_file = system.file("extdata", "example_repmap.tab", package="jacks", mustWork=TRUE)
  ref_spec_file = system.file("extdata", "example_grna_R_JACKS_results.txt", package="jacks", mustWork=TRUE)
  lfc = read_counts_from_spec_files(count_file, sample_spec_file, replicate_col="Replicate", sample_col="Sample", gene_spec_file=count_file, grna_col="sgRNA", gene_col="gene", count_prior=32., normalization='median', window=800, reference_sample="CTRL")
  rjacks = infer_jacks(lfc, c("KRAS"),reference_library=ref_spec_file)
  
  #To produce these run 
  pyjacks_genefile = system.file("extdata", "pyjacks_for_R_comparison_ref_gene_JACKS_results.txt", package="jacks", mustWork=TRUE)
  pyjacks_genestdfile = system.file("extdata", "pyjacks_for_R_comparison_ref_genestd_JACKS_results.txt", package="jacks", mustWork=TRUE)
  pyjacks_grnafile = system.file("extdata", "pyjacks_for_R_comparison_ref_grna_JACKS_results.txt", package="jacks", mustWork=TRUE)
  
  py_wmean <- read.table(file=pyjacks_genefile, header=TRUE)
  py_wstd <- read.table(file=pyjacks_genestdfile, header=TRUE)
  py_x <- read.table(file=pyjacks_grnafile, header=TRUE)    
  
  failed = FALSE
  for(i in 1:length(colnames(rjacks))){
    cline <- colnames(rjacks)[i]
    r_w <- metadata(rjacks)$jacks_w$KRAS[i]
    if( abs(r_w - py_wmean[cline]) > mean_tol ){
      flog.info(paste(cline, "\tpyjacks:", toString(py_wmean[cline]), "\tRjacks:", toString(r_w)))
      failed = TRUE
    }
  }
  for(i in 1:length(colnames(rjacks))){
    cline <- colnames(rjacks)[i]
    r_wstd <- metadata(rjacks)$jacks_sdw$KRAS[i]    
    if( abs(r_wstd - py_wstd[cline]) > std_tol ){
      flog.info(paste(cline, "\tpyjacks:", toString(py_wstd[cline]), "\tRjacks:", toString(r_wstd)))
      failed = TRUE
    }
  }
  if(failed){ stop("Gene KRAS inference not concordant between Rjacks and pyjacks") }
  flog.info(paste("pyjacks reference concordant with Rjacks on gene KRAS"))
}

#' Run tests
#'
#'\code{run_tests()} runs all defined unit tests. So far, these include inference from simulations,
#'real data, and using a reference gRNA efficacy definition; as well as concordance of
#'python and R versions of JACKS.
#'@return nothing. Execution stops if errors are found.
#'@export
run_tests <- function(){
    .test_inference()
    .test_pyjacks_reference_concordance()
    .test_pyjacks_concordance()
    .test_realdata_inference()
}
