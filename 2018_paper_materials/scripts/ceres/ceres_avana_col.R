library(ceres)

### Setup

# Edit this line to point to data directory
data_dir <- "../../data"
ceres_data_dir <- "ceres_data"

cn_seg_file <- file.path(ceres_data_dir, "CCLE_copynumber_2013-12-03.seg.txt")
gene_annot_file <- file.path(ceres_data_dir, "CCDS.current.txt")

# Set bowtie index directory. Not needed if $BOWTIE_INDEXES environmental variable is set and includes hg19 index.
Sys.setenv(BOWTIE_INDEXES = file.path(ceres_data_dir, "bowtie_indexes"))

avana_dep_file <- file.path("../../data","Avana_sgrnareplicateslogfcnorm.gct")
avana_rep_map <- file.path(".", "Avana_replicatemap_colon.tsv")


### Run CERES on Gecko data

inputs_dir <- file.path("./avana_ceres_inputs_col", Sys.Date())

prepare_ceres_inputs(inputs_dir=inputs_dir,
                     dep_file=avana_dep_file,
                     cn_seg_file=cn_seg_file,
                     gene_annot_file=gene_annot_file,
                     rep_map_file=avana_rep_map,
                     chromosomes=paste0("chr", c(as.character(1:22),"X","Y")),
                     dep_normalize="zmad")

avana_ceres <-
    wrap_ceres(sg_path=file.path(inputs_dir, "guide_sample_dep.Rds"),
               cn_path=file.path(inputs_dir, "locus_sample_cn.Rds"),
               guide_locus_path=file.path(inputs_dir, "guide_locus.Rds"),
               locus_gene_path=file.path(inputs_dir, "locus_gene.Rds"),
               replicate_map_path=file.path(inputs_dir, "replicate_map.Rds"),
               run_id="Avana",
               params=list(lambda_g=0.561))

avana_ceres_scaled <-
    scale_to_essentials(avana_ceres$gene_essentiality_results$ge_fit)

write.csv(avana_ceres_scaled, file="avana_ceres_col.csv")

