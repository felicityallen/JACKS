#'@import futile.logger
#'@import hashmap
library(utils)
library(hashmap)

#'Read gRNA counts given sample and gene specification files
#'
#'\code{read_counts_from_spec_files()} takes a file containing raw read counts, a replicate to sample specification file, and a gRNA ID to gene specification file, and returns a \code{\link{SummarizedExperiment}}
#'object of normalized log2-scale read counts and standard deviations.
#'
#'@param count_file String path to read count file. See \code{read_count_table_from_spec} for requirements.
#'@param sample_spec_file String path to sample specification file. Has to contain \code{replicate_col} and \code{sample_col} columns.
#'@param replicate_col String Name of the column in \code{sample_spec_file} that contains identifiers for each experimental replicate, as used as the column header for that replicate in the corresponding count file (globally unique identifiers per line, treatment, replicate, etc., not just '1','2', 'A', or 'R1').
#'@param sample_col String Name of the column in \code{sample_spec_file} that contains the sample identifiers. For example, if three columns in the count file are replicate measurements of the same sample, they should all have the same sample id here.
#'@param reference_sample String Sample name of the single reference sample all other samples will be compared against
#'@param gene_spec_file String path to gRNA-gene link specification file. Has to contain \code{grna_col} and \code{gene_col} columns.
#'@param grna_col String Name of the column in \code{gene_spec_file} that contains the gRNA identifiers.
#'@param gene_col String Name of the column in \code{gene_spec_file} that contains the gene identifiers (to identify which gRNAs map to the same gene).
#'@param count_prior Scalar double pseudocounts for each gRNA added to observations. Default: 32.
#'@param normalization String per-sample normalization method. Default: 'median' [sets median log2-scale count to 0]
#'@param window Scalar integer Length of smoothing window used in standard deviation estimation. Default: 800.
#'@return \code{\link{SummarizedExperiment}} object of log2-scale read count changes compared to reference, and estimated standard deviations across replicates
#'@export
read_counts_from_spec_files <- function(count_file, sample_spec_file, replicate_col, sample_col, reference_sample, gene_spec_file, grna_col, gene_col, count_prior=32., normalization='median', window=800){
    sample_spec = .read_sample_spec(count_file, sample_spec_file, replicate_col, sample_col)
    gene_spec = .read_gene_spec(gene_spec_file, grna_col, gene_col)
    logf = read_count_table_from_spec(sample_spec, gene_spec, count_prior=count_prior, normalization=normalization, window=window)
    colData(logf)$Condition = "TEST"
    colData(logf)$Condition = as.character(colData(logf)$Condition)
    colData(logf)$Condition[as.character(colData(logf)$Name) == reference_sample] = "CTRL"
    colData(logf)$Condition = as.factor(colData(logf)$Condition)
    return(calc_logfc(logf, reference="CTRL", condition="TEST"))
}

.read_sample_spec <- function(count_file, sample_spec_file, replicate_col, sample_col){
    sample_spec = utils::read.table(sample_spec_file, header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
    if(!(replicate_col %in% colnames(sample_spec))){
        flog.error(paste("Designated replicate column '", replicate_col, "' not found in columns of sample specification file ", sample_spec_file))
    } else if (!(sample_col %in% colnames(sample_spec))){
        flog.error(paste("Designated sample column '", sample_col, "' not found in columns of sample specification file ", sample_spec_file))
    }
    sample_spec = sample_spec[,c(which(colnames(sample_spec) == replicate_col), which(colnames(sample_spec) == sample_col))]
    sample_spec = cbind(rep(count_file, dim(sample_spec)[1]), sample_spec) # add count file info
    colnames(sample_spec) = c("Filename", "Replicate", "Sample")
    return(sample_spec)
}

.read_gene_spec <- function(gene_spec_file, grna_col, gene_col){
    gene_spec = utils::read.table(gene_spec_file, header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
    if (!(grna_col %in% colnames(gene_spec))){
        flog.error(paste("Designated gRNA ID column '", grna_col, "' not found in columns of gRNA gene specification file ", gene_spec_file))
    } else if (!(gene_col %in% colnames(gene_spec))){
        flog.error(paste("Designated gene column '", gene_col, "' not found in columns of gRNA gene specification file ", gene_spec_file))
    }
    gene_spec = gene_spec[,c(which(colnames(gene_spec) == grna_col), which(colnames(gene_spec) == gene_col))]
    return(gene_spec)
}

#'Read gRNA counts given sample and gene specification
#'
#'\code{read_count_table_from_spec()} takes a sample specification, and returns a \code{\link{SummarizedExperiment}}
#'object of read counts.
#'
#'The sample specification table provides names of count files to read, columns in those files to retain, and the samples these correspond to.
#'For example,
#'\code{}
#'JACKS assumes that each count file has the same number of rows ordered in the same way, corresponding to the same gRNAs.
#'Every count file has one header line, and the rest of the lines give raw gRNA read counts for one gRNA each, with each column being a different measurement.
#'Only columns with headers as given in the "Replicate" (2nd column of the specification) are stored from each count file.
#'Count file columns that map to the same "Sample" (3rd column of the specification) are treated as replicate measurements for that sample.
#'
#'@param sample_spec Data frame of at least three columns: "Filename": input count file, "Replicate": column names in count file, "Sample": corresponding sample ID to allow combining of replicates within samples.
#'@param gene_spec Data frame of at least two columns: gRNA ID (1st col) and corresponding gene (2nd col)
#'@param count_prior Scalar double pseudocounts for each gRNA added to observations. Default: 32.
#'@param normalization String per-sample normalization method. Default: 'median' [sets median log2-scale count to 0]
#'@param window Scalar integer smoothing window used in standard deviation estimation. Default: 800.
#'@return \code{\link{SummarizedExperiment}} object of log-scale sample read counts and estimated standard deviations across replicates
#'@export
read_count_table_from_spec <- function(sample_spec, gene_spec, count_prior=32., normalization='median', window=800){
    all_samples = c()
    all_counts = c()
    genehash = hashmap(as.character(gene_spec[,1]), as.character(gene_spec[,2]))
    for(filename in unique(sample_spec[,1])){ # column 1 of sample spec. is filename
        flog.debug(paste("Reading samples from", filename))
        counts = c()
        meta = c()
        d = .read_sample_counts(filename, sample_spec)
        grnas = rownames(d)
        I = rep(T, dim(d)[1]) # gRNAs to use
        for(i in 1:dim(d)[1]){ I[i] = genehash$has_key(grnas[i]) } # see which gRNAs can be mapped
        meta = cbind(grnas[I], rep(NA, sum(I))) # retain gRNA metadata info - no gene to begin with
        for(i in 1:dim(meta)[1]){ meta[i,2] = genehash[[meta[i,1]]] } # and fill in

        if(is.null(all_counts)) { # if first file, set result
            all_counts = d[I,]
            all_samples = .get_sample_ids(colnames(d), sample_spec, filename)
        } else { # some results already exist - append
            all_counts = cbind(all_counts, d[I,]) # assume all files are aligned
            all_samples = c(all_samples, .get_sample_ids(colnames(d), sample_spec, filename))
        }
    }
    return(.generate_logfc_summarized_experiment(meta, all_counts, all_samples, count_prior, normalization, window))
}

.read_sample_counts <- function(filename, sample_spec){
    sep = "\t"
    if(".csv" %in% filename){
        sep = ","
    }
    d = utils::read.table(filename, sep=sep,header=TRUE, stringsAsFactors=FALSE, row.names=1, check.names=FALSE)
    file_samples = unique(sample_spec[filename == sample_spec[,1], 2]) # Specification column 1 = file name, column 2 = samples to retain from this file
    I_sample = rep(F, dim(d)[2])
    for(j in 1:length(I_sample)){
        I_sample[j] = colnames(d)[j] %in% file_samples
    }
    return(d[,I_sample])
}

.get_sample_ids <- function(colnames, sample_spec, filename){
    sample_ids = c()
    spec = sample_spec[sample_spec[,1] == filename,] # column 1 is filename
    for(i in 1:length(colnames)){
        for(j in 1:dim(spec)[1]){
            if(colnames[i] == spec[j,2]){ # Spec. column 2 is sample name in input file
                sample_ids = c(sample_ids, spec[j,3]) # Spec. column 3 is sample ID
            }
        }
    }
    return(sample_ids)
}

.generate_logfc_summarized_experiment <- function(grna_meta, counts, count_samples, count_prior=32., normalization='median', window=800){
    # typecast vectors of vectors to matrices, add names
    grna_meta = matrix(unlist(grna_meta), ncol=2, dimnames=list(NULL, c("gRNA", "gene")))
    counts = matrix(unlist(counts), nrow=dim(grna_meta)[1], dimnames=list(grna_meta[,1], count_samples))
    samples = unique(count_samples)

    # create new matrix of log-scale centered and counts
    meanmat = matrix(ncol=length(samples), nrow=dim(counts)[1], dimnames=list(grna_meta[,1], samples))
    sdmat = matrix(ncol=length(samples), nrow=dim(counts)[1], dimnames=list(grna_meta[,1], samples))

    for (i in 1:length(samples)){ # for each sample
        flog.debug(paste("Calculating variance estimates for", samples[i]))
        I = (count_samples == samples[i]) # pick corresponding counts from replicates
        d = .calc_mean_sd(counts[,I], count_prior, window) # calculate
        meanmat[,i] = as.vector(d$mean) # and store
        sdmat[,i] = as.vector(d$sd)
    }
    # assemble object - metadata of matrices, row, and column entries
    colnames(grna_meta) = c("gRNA", "gene")
    sample_meta = data.frame(list(Name=samples, Condition=rep("NA", length(samples))), check.names=FALSE)
    rownames(sample_meta) = samples
    result = SummarizedExperiment(assays=list("logf_mean"=meanmat, "logf_sd"=sdmat), rowData=grna_meta, colData=sample_meta)
    validate_logf_table(result)
    flog.info(paste("Done condensing replicates; went from", length(count_samples), "experiments down to", length(samples), "samples"))
    return(result)
}

.read_precomputed_x <- function(library){
    supported = c("avana","gecko2","yusa_v10")
    if( library %in% supported ){ ref = url(paste0("https://raw.githubusercontent.com/felicityallen/JACKS/master/reference_grna_efficacies/", library, "_grna_JACKS_results.txt"))}
    else{ ref = library }
    x = utils::read.table(ref, header=TRUE, sep="\t", check.names=FALSE, stringsAsFactors=FALSE, row.names = 1)
    return(x)
}

#'Extract JACKS output for a gene
#'
#'\code{jacks_w_gene()} takes JACKS output and gene name, and returns a \code{data.frame}
#'that has columns for statistics of JACKS inferred posterior for gene essentiality w - mean
#'and standard deviation. Each row is one cell line, name in row.names.
#'
#'@param expt Output of JACKS inference - SummarizedExperiment endowed with posteriors.
#'@param gene Gene name to extract information for.
#'@return \code{data.table} object estimated means and standard deviations of gene essentiality.
#'@export
jacks_w_gene <- function(expt, gene){
    data.frame(
        row.names = rownames(colData(expt)),
        w = metadata(expt)$jacks_w[[gene]],
        sd_w = metadata(expt)$jacks_sdw[[gene]],
        neg_pval = metadata(expt)$jacks_neg_pval[[gene]],
        pos_pval = metadata(expt)$jacks_pos_pval[[gene]],
        neg_fdr = metadata(expt)$jacks_neg_fdr[[gene]],
        pos_fdr = metadata(expt)$jacks_pos_fdr[[gene]]
    )
}

#'Extract JACKS output for a sample
#'
#'\code{jacks_w_sample()} takes JACKS output and sample name, and returns a \code{data.frame}
#'that has columns for statistics of JACKS inferred posterior for gene essentiality w - mean
#'and standard deviation. Each row is one gene (name in row.names).
#'
#'@param expt Output of JACKS inference - SummarizedExperiment endowed with posteriors.
#'@param sample Sample name to extract information for.
#'@return \code{data.table} object estimated means and standard deviations of gene essentiality.
#'@export
jacks_w_sample <- function(expt, sample){
    i = which(row.names(colData(expt)) == sample)
    m = metadata(expt)
    data.frame(
        row.names = colnames(m$jacks_w),
        w = m$jacks_w[i,],
        sd_w = m$jacks_sdw[i,],
        neg_pval = m$jacks_neg_pval[i,],
        pos_pval = m$jacks_pos_pval[i,],
        neg_fdr = m$jacks_neg_fdr[i,],
        pos_fdr = m$jacks_pos_fdr[i,]
    )
}

#
#write_jacks_output <- function(){
#
#}

#'Pre-computed gRNA efficacy estimates for the Avana library.
#'Two vectors are provided, both sorted according to the gRNA sequence.
#' \itemize{
#'   \item x: Estimated gRNA efficacy for each gRNA. x = 1 means the guide works roughly as an average gRNA would. x = 0 means guide is ineffective.
#'   \item sdx: Standard deviation of x estimate.
#' }
#'
#' @format Two numeric vectors.
#' @source \url{http://www.diamondse.info/}
#' @name avana
#' @docType data
#' @author Leopold Parts \email{lp2@sanger.ac.uk}
#' @keywords data
NULL

#'Test
#' @name example_repmap
#' @docType data
#' @keywords data
NULL


#'Test
#' @name example_count_data
#' @docType data
#' @keywords data
NULL

#'Test
#' @name avana_head
#' @docType data
#' @keywords data
NULL

#'Test
#' @name data
#' @docType data
#' @keywords data
NULL

#'Test
#' @name data_err
#' @docType data
#' @keywords data
NULL

#'Test
#' @name pyvals
#' @docType data
#' @keywords data
NULL

#'Test
#' @name x
#' @docType data
#' @keywords data
NULL

#'Test
#' @name sdx
#' @docType data
#' @keywords data
NULL
