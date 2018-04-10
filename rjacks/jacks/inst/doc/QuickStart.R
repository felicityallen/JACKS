## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(jacks)

## ----counts, warning=FALSE, message=FALSE--------------------------------
    count_file = system.file("extdata", "example_count_data.tab", package="jacks", mustWork=TRUE)
    sample_spec_file = system.file("extdata", "example_repmap.tab", package="jacks", mustWork=TRUE)
    lfc = read_counts_from_spec_files(count_file, sample_spec_file, replicate_col="Replicate", sample_col="Sample", gene_spec_file=count_file, grna_col="sgRNA", gene_col="gene", count_prior=32., normalization='median', window=800, reference_sample="CTRL")

## ----JACKS, warning=FALSE, message=FALSE---------------------------------
test_genes = c("KRAS", "RRM2", "ZNF253")
result = infer_jacks(lfc, test_genes)
print(jacks_w_gene(result, "KRAS"), digits=3)

## ----plot, echo=TRUE, warning=FALSE--------------------------------------
p = plot_jacks(result, "KRAS", do_save=FALSE)

## ----jacks_ref, warning=FALSE, message=FALSE-----------------------------
result = infer_jacks(lfc, test_genes, reference_library="yusa_v10")
print(jacks_w_gene(result, "KRAS"), digits=3)

