#'@import ggplot2
#'@importFrom GGally ggpairs
library(cowplot)
library(ggplot2)
library(GGally)

.yhist <- function(y, gene, all_logfc){
    if(is.null(all_logfc)){
        return(ggplot(y, aes_string(x="log2fc", y="..density..")) + geom_histogram(fill="red"))
    }

    ycmp_line = data.frame(log2fc=as.vector(as.matrix(y$log2fc)))
    ycmp_line$gRNAs = gene
    ycmp_all = data.frame(log2fc=as.vector(as.matrix(all_logfc)))
    ycmp_all$gRNAs = "All"
    ycmp = rbind(ycmp_line, ycmp_all)
    ggplot(ycmp, aes_string("log2fc", fill="gRNAs")) +
      geom_density(alpha=0.3) +
      theme(legend.position=c(0.6, 0.4), legend.text=element_text(size=8), legend.title=element_text(size=10))
}

.yheat <-function(y, gene){
    ggplot(y, aes_string("CellLine", "gRNA")) +
      ggtitle(paste(gene, "log2 fold change (measured)")) +
      geom_tile(aes_string(fill="log2fc"), colour="black") +
      scale_fill_gradient2(low="steelblue", mid="white", high="red", midpoint=0, limits=c(-3,3)) +
      theme(axis.text.x=element_text(angle=0, size=8)) + #, axis.text.y=element_text(size=8)) +
      theme(legend.title=element_blank(), axis.line=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())
}

.predheat <- function(y){
    ggplot(y, aes_string("CellLine", "gRNA")) +
      ggtitle("Prediction") +
      geom_tile(aes_string(fill="Prediction"), colour="black") +
      scale_fill_gradient2(low="steelblue", mid="white", high="red", midpoint=0, limits=c(-3,3)) +
      theme(axis.text.x=element_text(angle=0, size=8)) + #, axis.text.y=element_text(size=8))
      theme(legend.title=element_blank(), axis.line=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())
}

.xviolin <- function(x1, sds_x, n=10000){
    k = length(x1)
    xs = data.frame(x=matrix(matrix(stats::rnorm(k*n))*sds_x+x1, nrow=k*n), gRNA=as.factor(rep(1:k,n)))

    ggplot(xs, aes_string(x="gRNA", y="x")) +
      geom_violin(aes_string(fill="factor(gRNA)"), draw_quantiles=c(0.5)) +
      geom_hline(yintercept=1, alpha=0.5, linetype="dashed") +
      coord_flip() +
      theme(legend.position="none")
}

.wviolin <- function(w1, sds_w, lines, n=10000){
    l = length(w1)
    ws = data.frame(w=matrix(matrix(stats::rnorm(l*n))*sds_w+w1, nrow=l*n), CellLine=as.factor(rep(lines,n)))

    ggplot(ws, aes_string(x="CellLine", y="w")) +
      geom_violin(aes_string(fill="factor(CellLine)"), draw_quantiles=c(0.5)) +
      geom_hline(yintercept=0, alpha=0.5, linetype="dashed") +
      theme(legend.position="none", plot.margin=margin(l=0.6, r=2.4, unit="cm"))
}

.predscatter <- function(y){
    ggplot(y, aes_string(x="log2fc", y="Prediction", color="CellLine", shape="gRNA")) +
      geom_point() +
      geom_errorbarh(aes_string(xmin="log2fc-MeasurementError", xmax="log2fc+MeasurementError"), alpha=0.2) +
      geom_errorbar(aes_string(ymin="Prediction-PredictionError", ymax="Prediction+PredictionError"), alpha=0.2) +
      geom_abline(slope=1, intercept=0, alpha=0.5, linetype="dashed") + guides(shape=FALSE) +
      xlim(-3,0.5) + ylim(-3,0.5) +
      theme(legend.position=c(0.6, 0.4), legend.text=element_text(size=8), legend.title=element_text(size=12))
}

.prep_y <- function(r, lines, gene){
    line_factor = as.factor(rep(lines, each=length(r$x1)))
    grna_factor = as.factor(rep(paste("gRNA", 1:length(r$x1)), length(r$w1)))
    data.frame(log2fc=as.vector(r$y),
        Prediction=as.vector(outer(r$x1, r$w1)),
        PredictionError=as.vector(1./r$tau**0.5),
        MeasurementError=as.vector(r$err),
        CellLine=line_factor,
        gRNA=grna_factor)
}

# .prep_y_object <- function(d){
#     lines = colnames(d)
#     w1 = metadata(d)$jacks_w
#     sdw = metadata(d)$jacks_sdw
#     x1 = rowData(d)$jacks_x
#     sdx = rowData(d)$jacks_sdx
#     pred = outer(x1, w1)
#
#     line_factor = as.factor(rep(colnames(d), each=length(x1)))
#     grna_factor = as.factor(rep(paste("gRNA", 1:length(x1)), length(w1)))
#     data.frame(log2fc=as.vector(r$y),
#                Prediction=as.vector(outer(r$x1, r$w1)),
#                PredictionError=as.vector(1./r$tau**0.5),
#                MeasurementError=as.vector(r$err),
#                CellLine=line_factor,
#                gRNA=grna_factor)
# }

.prep_y_object_summarizedexperiment <- function(d, gene){
    lines = colnames(d)
    I = rowData(d)$gene == gene
    w1 = metadata(d)$jacks_w[[gene]]
    sdw = metadata(d)$jacks_sdw[[gene]]
    x1 = rowData(d)$jacks_x[I]
    sdx = rowData(d)$jacks_sdx[I]
    pred = outer(x1, w1)
    line_factor = as.factor(rep(colnames(d), each=length(x1)))
    grna_factor = as.factor(rep(paste("gRNA", 1:length(x1)), length(lines)))
    data.frame(log2fc=as.vector(as.matrix(assays(d)$logfc_mean[I,])),
                   Prediction=as.vector(outer(x1, w1)),
                   PredictionError=as.vector(as.matrix(assays(d)$prediction_sd[I,])),
                   MeasurementError=as.vector(as.matrix(assays(d)$logfc_sd[I,])),
                   CellLine=line_factor,
                   gRNA=grna_factor)
}

#' Plot JACKS result for one gene.
#'
#'\code{plot_jacks()} plots the JACKS decomposition of gRNA log2-fold changes for one gene across
#'several experiments into an experiment-specific gene effect, and experiment-independent gRNA efficacy.
#'Input can either be a \link{SummarizedExperiment} object returned by \code{\link{infer_jacks}} or
#'a list of posteriors returned by \code{\link{infer_jacks_gene}}.
#'@param input A \link{SummarizedExperiment} object returned by \code{\link{infer_jacks}} or a list of posteriors returned by \code{\link{infer_jacks_gene}}.
#'@param gene String identifying the gene for which to plot results
#'@param lines Vector of cell line (experiment) names used. Only used when list of posteriors without annotation is given as input. Default NULL.
#'@param all_logfc Vector of all log-fold changes in the experiment. Only used when list of posteriors without annotation is given as input. Default NULL.
#'@param do_save Logical scalar of whether to save a .pdf of the output. Default T.
#'@return nothing. Execution stops if errors are found.
#'@export
plot_jacks <- function(input, gene, lines=NULL, all_logfc=NULL, do_save=TRUE){
    if(class(input)[1] == "SummarizedExperiment"){
        .plot_jacks_summarizedexperiment(input, gene, do_save)
    } else{
        .plot_jacks_generesult(input, lines, gene, do_save=do_save)
    }
}

.plot_jacks_generesult <- function(r, lines, gene, all_logfc=NULL, do_save=TRUE){
    y = .prep_y(r, lines, gene)
    yhist = .yhist(y, gene, all_logfc)
    yheat = .yheat(y, gene)
    predheat = .predheat(y)
    xviol = .xviolin(r$x1, (r$x2-r$x1^2)^0.5)
    wviol = .wviolin(r$w1, (r$w2-r$w1^2)^0.5, lines)
    predscatter = .predscatter(y)
    .plot_jacks_generesult_grid(yhist, wviol, xviol, yheat, predscatter, predheat, gene, do_save=do_save)
}

.plot_jacks_summarizedexperiment <- function(d, gene, do_save=TRUE){
    y <- .prep_y_object_summarizedexperiment(d, gene)
    yhist = .yhist(y, gene, assays(d)$logfc_mean)
    yheat = .yheat(y, gene)
    predheat = .predheat(y)

    I = rowData(d)$gene == gene
    xviol = .xviolin(rowData(d)$jacks_x[I], rowData(d)$jacks_sdx[I])
    wviol = .wviolin(metadata(d)$jacks_w[[gene]], metadata(d)$jacks_sdw[[gene]], colnames(d))
    predscatter = .predscatter(y)

    .plot_jacks_generesult_grid(yhist, wviol, xviol, yheat, predscatter, predheat, gene, do_save=do_save)
}


.plot_jacks_generesult_grid <- function(yhist, wviol, xviol, yheat, predscatter, predheat, gene, do_save=TRUE){
    row1 = cowplot::plot_grid(yhist, wviol, labels=c('',''), align='h', rel_widths = c(1,2))
    row2 = cowplot::plot_grid(xviol, yheat, labels=c('',''), align='h', rel_widths = c(1,2))
    row3 = cowplot::plot_grid(predscatter, predheat, align='h', rel_widths=c(1, 2))
    p = cowplot::plot_grid(row1, row2, row3, labels=rep('',3), ncol=1, rel_heights = c(0.6,1,1))

    if(do_save){
        ggsave(paste0("jacks_",gene,"_plots.pdf"), width=8.24, height=6.74)
    } else{
        print(p)
    }
    return(p)
}

.plot_count_replicates <- function(counts, sample_name, sample_condition){
    my_scatter <- function(data, mapping, ...) {ggplot(data=data, mapping=mapping) + geom_point(..., alpha=0.2) }
    samples = colData(counts)
    I = (samples$Name == sample_name) & (samples$Condition == sample_condition)
    GGally::ggpairs(as.data.frame(log2(assays(counts)$counts + 1)[,I]), lower=list(continuous=my_scatter))
}

.plot_meanvar <- function(counts, condensed, sample, count_prior=16.){
    count_samples = colData(counts)
    I = (count_samples$Name == strsplit(sample,"_")[[1]][1]) & (count_samples$Condition == strsplit(sample,"_")[[1]][2])
    log2c = log2(assays(counts)$counts[,I] + count_prior)
    log2c = apply(log2c, 2, function(x) x - stats::median(x))
    d = data.frame("log2f"=apply(log2c,1,mean), "stddev"=apply(log2c,1,stats::sd), "type"="Data")
    cond_samples = colData(condensed)
    I = (cond_samples$Name == strsplit(sample,"_")[[1]][1]) & (cond_samples$Condition == strsplit(sample,"_")[[1]][2])
    prior_std = assays(condensed)$logf_sd[,I]
    d2 = data.frame("log2f"=assays(condensed)$logf_mean[,I], "stddev"=prior_std, "type"="Prior")
    d = rbind(d, d2)
    ggplot(d, aes_string(x="log2f", y="stddev", color="type")) + geom_point(alpha=0.05) + ylim(0, max(prior_std)+0.2) +
        scale_color_manual(values=c("black", "blue")) + theme(legend.position="none")
}
