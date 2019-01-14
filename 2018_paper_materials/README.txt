---------------------------------------------
Data and Results: Figshare location
--------------------------------------------

Data and Results files can be found at:

https://figshare.com/articles/Results/6002438

Code can be found at:

https://github.com/felicityallen/JACKS

To generate results (python and R files in scripts or jacks directory):

----------------------------------------
**JACKS results**
----------------------------------------

python run_JACKS.py data/Avana_sgrna_raw_readcounts_matched.csv data/Avana_replicatemap.csv data/Avana_sgrnamapping.csv --outprefix=results/avana_jacks_results/avana_jacks_results --rep_hdr=Replicate --sample_hdr=CellLine --common_ctrl_sample=pDNA --sgrna_hdr=Guide --gene_hdr=Gene --ctrl_genes=data/NEGv1.txt
python run_JACKS.py data/Avana_sgrna_raw_readcounts_rand_matched.csv data/Avana_replicatemap_rand.csv data/Avana_sgrnamapping.csv --outprefix=results/avana_jacks_results_rand/avana_jacks_results_rand --rep_hdr=Replicate --sample_hdr=CellLine --common_ctrl_sample=pDNA --sgrna_hdr=Guide --gene_hdr=Gene --ctrl_genes=data/NEGv1.txt

python run_JACKS.py data/Achilles_raw_GeckoV2.tab data/Achilles_raw_GeckoV2_repmap.tab data/Achilles_raw_GeckoV2.tab --outprefix=results/jacks_results_gecko2/jacks_results_gecko2 --rep_hdr=Replicate --sample_hdr=Sample --common_ctrl_sample=CTRL --sgrna_hdr=sgRNA --gene_hdr=gene --ctrl_genes=data/NEGv1.txt
python run_JACKS.py data/Achilles_raw_GeckoV2_rand.tab data/Achilles_raw_GeckoV2_repmap_rand.tab data/Achilles_raw_GeckoV2_rand.tab --outprefix=results/jacks_results_gecko2_rand/jacks_results_gecko2_rand --rep_hdr=Replicate --sample_hdr=Sample --common_ctrl_sample=CTRL --sgrna_hdr=sgRNA --gene_hdr=gene --ctrl_genes=data/NEGv1.txt

python run_JACKS.py data/fiona_AML_raw_v10.tab data/fiona_AML_raw_v10_repmap.tab data/fiona_AML_raw_v10.tab --outprefix=results/jacks_results_v10/jacks_results_v10 --rep_hdr=Replicate --sample_hdr=Sample --common_ctrl_sample=CTRL --sgrna_hdr=sgRNA --gene_hdr=gene --ctrl_genes=data/NEGv1.txt
python run_JACKS.py data/tko/tko_counts.txt data/tko/tko_repmap.txt data/tko/tko_counts.txt --outprefix=results/jacks_results_tko/jacks_results_tko --rep_hdr=Replicate --sample_hdr=Sample --ctrl_sample_hdr=Control --sgrna_hdr=GENE_CLONE --gene_hdr=GENE --ctrl_genes=data/NEGv1.txt
python run_JACKS.py data/whitehead/Wang2015_2017_merged_counts.txt data/whitehead/whitehead_repmap.txt data/whitehead/Wang2017_1-s2.0-S0092867417300612-mmc2.txt --outprefix=results/jacks_results_whitehead/jacks_results_whitehead --rep_hdr=Replicate --sample_hdr=Sample --ctrl_sample_hdr=Control --sgrna_hdr=sgRNA_ID --gene_hdr=Symbol --ctrl_genes=data/NEGv1.txt

----------------------------------------
**JACKS +HP results**
----------------------------------------

python run_JACKS.py data/Avana_sgrna_raw_readcounts_matched.csv data/Avana_replicatemap.csv data/Avana_sgrnamapping.csv --outprefix=results/avana_jacks_results_hp/avana_jacks_results_hp --apply_w_hp --rep_hdr=Replicate --sample_hdr=CellLine --common_ctrl_sample=pDNA --sgrna_hdr=Guide --gene_hdr=Gene --ctrl_genes=data/NEGv1.txt
python run_JACKS.py data/Avana_sgrna_raw_readcounts_rand_matched.csv data/Avana_replicatemap_rand.csv data/Avana_sgrnamapping.csv --outprefix=results/avana_jacks_results_rand_hp/avana_jacks_results_rand_hp --apply_w_hp --rep_hdr=Replicate --sample_hdr=CellLine --common_ctrl_sample=pDNA --sgrna_hdr=Guide --gene_hdr=Gene --ctrl_genes=data/NEGv1.txt

python run_JACKS.py data/Achilles_raw_GeckoV2.tab data/Achilles_raw_GeckoV2_repmap.tab data/Achilles_raw_GeckoV2.tab --outprefix=results/jacks_results_gecko2_hp/jacks_results_gecko2_hp --apply_w_hp --rep_hdr=Replicate --sample_hdr=Sample --common_ctrl_sample=CTRL --sgrna_hdr=sgRNA --gene_hdr=gene --ctrl_genes=data/NEGv1.txt
python run_JACKS.py data/Achilles_raw_GeckoV2_rand.tab data/Achilles_raw_GeckoV2_repmap_rand.tab data/Achilles_raw_GeckoV2_rand.tab --outprefix=results/jacks_results_gecko2_rand_hp/jacks_results_gecko2_rand_hp --apply_w_hp --rep_hdr=Replicate --sample_hdr=Sample --common_ctrl_sample=CTRL --sgrna_hdr=sgRNA --gene_hdr=gene --ctrl_genes=data/NEGv1.txt

----------------------------------------
**JACKS x results**
----------------------------------------

for <num_celllines> in test cases:
    python scripts/jacks/sample_jacks_screen_xs.py Avana_raw_reformatted_Batch0.txt 2 <num_celllines> avana0_xs/xs_2_<num_celllines>.txt 100
    python scripts/jacks/sample_jacks_screen_xs.py Avana_raw_reformatted_Batch1.txt 2 <num_celllines> avana1_xs/xs_2_<num_celllines>.txt 100

----------------------------------------
**JACKS Colon/Melanoma BRAF Results**
----------------------------------------

python run_JACKS.py data/Avana_sgrna_raw_readcounts_matched.csv data/Avana_replicatemap_colon_nad.csv data/Avana_sgrnamapping.csv --outprefix=results/avana_jacks_results_col_nad/avana_jacks_results_col_nad --rep_hdr=Replicate --sample_hdr=CellLine --common_ctrl_sample=pDNA --sgrna_hdr=Guide --gene_hdr=Gene --ctrl_genes=data/NEGv1.txt
python run_JACKS.py data/Avana_sgrna_raw_readcounts_matched.csv data/Avana_replicatemap_melanoma_nad.csv data/Avana_sgrnamapping.csv --outprefix=results/avana_jacks_results_mel_nad/avana_jacks_results_mel_nad --rep_hdr=Replicate --sample_hdr=CellLine --common_ctrl_sample=pDNA --sgrna_hdr=Guide --gene_hdr=Gene --ctrl_genes=data/NEGv1.txt
python run_JACKS.py data/Avana_sgrna_raw_readcounts_matched.csv data/Avana_replicatemap_colonmelanoma_noampdel.csv data/Avana_sgrnamapping.csv --outprefix=results/avana_jacks_results_melcol_nad/avana_jacks_results_melcol_nad --rep_hdr=Replicate --sample_hdr=CellLine --common_ctrl_sample=pDNA --sgrna_hdr=Guide --gene_hdr=Gene --ctrl_genes=data/NEGv1.txt

----------------------------------------
**JACKS Single Line Results**
----------------------------------------

python run_JACKS_single.py data/fiona_AML_raw_v10.tab data/fiona_AML_raw_v10_repmap.tab data/fiona_AML_raw_v10.tab --outprefix=results/single_jacks_results_v10/single_jacks_results_v10 --rep_hdr=Replicate --sample_hdr=Sample --common_ctrl_sample=CTRL --sgrna_hdr=sgRNA --gene_hdr=gene --ctrl_genes=data/NEGv1.txt --separate
python run_JACKS_single.py data/Avana_sgrna_raw_readcounts_matched.csv data/Avana_replicatemap.csv data/Avana_sgrnamapping.csv --outprefix=results/single_avana_jacks_results/single_avana_jacks_results --rep_hdr=Replicate --sample_hdr=CellLine --common_ctrl_sample=pDNA --sgrna_hdr=Guide --gene_hdr=Gene --ctrl_genes=data/NEGv1.txt --separate
python run_JACKS_single.py data/Achilles_raw_GeckoV2.tab data/Achilles_raw_GeckoV2_repmap.tab data/Achilles_raw_GeckoV2.tab --outprefix=results/single_jacks_results_gecko2/single_jacks_results_gecko2 --rep_hdr=Replicate --sample_hdr=Sample --common_ctrl_sample=CTRL --sgrna_hdr=sgRNA --gene_hdr=gene --ctrl_genes=data/NEGv1.txt --separate
python run_JACKS_single.py data/whitehead/Wang2015_2017_merged_counts.txt data/whitehead/whitehead_repmap.txt data/whitehead/Wang2017_1-s2.0-S0092867417300612-mmc2.txt --outprefix=results/single_jacks_results_whitehead/single_jacks_results_whitehead --rep_hdr=Replicate --sample_hdr=Sample --ctrl_sample_hdr=Control --sgrna_hdr=sgRNA_ID --gene_hdr=Symbol --ctrl_genes=data/NEGv1.txt --separate
python run_JACKS_single.py data/tko/tko_counts.txt data/tko/tko_repmap.txt data/tko/tko_counts.txt --outprefix=results/single_jacks_results_tko/single_jacks_results_tko --rep_hdr=Replicate --sample_hdr=Sample --ctrl_sample_hdr=Control --sgrna_hdr=GENE_CLONE --gene_hdr=GENE --ctrl_genes=data/NEGv1.txt --separate

----------------------------------------
**JACKS Reference Leave-One-Out Results**
----------------------------------------

python leaveoneout_reference_test.py data/fiona_AML_raw_v10.tab data/fiona_AML_raw_v10_repmap.tab data/fiona_AML_raw_v10.tab --outprefix=results/loor_jacks_results_v10/loor_jacks_results_v10 --rep_hdr=Replicate --sample_hdr=Sample --common_ctrl_sample=CTRL --sgrna_hdr=sgRNA --gene_hdr=gene --ctrl_genes=data/NEGv1.txt --separate
python leaveoneout_reference_test.py data/Avana_sgrna_raw_readcounts_matched.csv data/Avana_replicatemap.csv data/Avana_sgrnamapping.csv --outprefix=results/loor_jacks_avana_results/loor_jacks_avana_results --rep_hdr=Replicate --sample_hdr=CellLine --common_ctrl_sample=pDNA --sgrna_hdr=Guide --gene_hdr=Gene --ctrl_genes=data/NEGv1.txt --separate
python leaveoneout_reference_test.py data/Achilles_raw_GeckoV2.tab data/Achilles_raw_GeckoV2_repmap.tab data/Achilles_raw_GeckoV2.tab --outprefix=results/loor_jacks_results_gecko2/loor_jacks_results_gecko2 --rep_hdr=Replicate --sample_hdr=Sample --common_ctrl_sample=CTRL --sgrna_hdr=sgRNA --gene_hdr=gene --ctrl_genes=data/NEGv1.txt --separate
python leaveoneout_reference_test.py data/whitehead/Wang2015_2017_merged_counts.txt data/whitehead/whitehead_repmap.txt data/whitehead/Wang2017_1-s2.0-S0092867417300612-mmc2.txt --outprefix=results/loor_jacks_results_whitehead/loor_jacks_results_whitehead --rep_hdr=Replicate --sample_hdr=Sample --ctrl_sample_hdr=Control --sgrna_hdr=sgRNA_ID --gene_hdr=Symbol --ctrl_genes=data/NEGv1.txt --separate
python leaveoneout_reference_test.py data/tko/tko_counts.txt data/tko/tko_repmap.txt data/tko/tko_counts.txt --outprefix=results/loor_jacks_results_tko/loor_jacks_results_tko --rep_hdr=Replicate --sample_hdr=Sample --ctrl_sample_hdr=Control --sgrna_hdr=GENE_CLONE --gene_hdr=GENE --ctrl_genes=data/NEGv1.txt --separate

----------------------------------------
**Bagel Results**
----------------------------------------

python run_bagel.py data/Avana_sgrna_raw_readcounts_matched.csv data/Avana_replicatemap.csv data/Avana_sgrnamapping.csv --rep_hdr=Replicate --sample_hdr=CellLine --common_ctrl_sample=pDNA --sgrna_hdr=Guide --gene_hdr=Gene --v10=avana
python run_bagel.py data/Achilles_raw_GeckoV2.tab data/Achilles_raw_GeckoV2_repmap.tab data/Achilles_raw_GeckoV2.tab --rep_hdr=Replicate --sample_hdr=Sample --common_ctrl_sample=CTRL --sgrna_hdr=sgRNA --gene_hdr=gene --v10=gecko2
python run_bagel.py data/fiona_AML_raw_v10.tab data/fiona_AML_raw_v10_repmap.tab data/fiona_AML_raw_v10.tab --rep_hdr=Replicate --sample_hdr=Sample --common_ctrl_sample=CTRL --sgrna_hdr=sgRNA --gene_hdr=gene --v10=v10
python run_bagel.py data/whitehead/Wang2015_2017_merged_counts.txt data/whitehead/whitehead_repmap.txt data/whitehead/Wang2017_1-s2.0-S0092867417300612-mmc2.txt --sample_hdr=Sample --rep_hdr=Replicate --ctrl_sample_hdr=Control --sgrna_hdr=sgRNA_ID --gene_hdr=Symbol --v10=whitehead
python run_bagel.py data/tko/tko_counts.txt data/tko/tko_repmap.txt data/tko/tko_counts.txt --rep_hdr=Replicate --sample_hdr=Sample --ctrl_sample_hdr=Control --sgrna_hdr=GENE_CLONE --gene_hdr=GENE --ignore_blank_genes --v10=tko
    
----------------------------------------
**MAGeCK Results**
----------------------------------------

python run_mageck.py data/Avana_sgrna_raw_readcounts_matched.csv data/Avana_replicatemap.csv data/Avana_sgrnamapping.csv --rep_hdr=Replicate --sample_hdr=CellLine --common_ctrl_sample=pDNA --sgrna_hdr=Guide --gene_hdr=Gene --v10=avana
python run_mageck.py data/Achilles_raw_GeckoV2.tab data/Achilles_raw_GeckoV2_repmap.tab data/Achilles_raw_GeckoV2.tab --rep_hdr=Replicate --sample_hdr=Sample --common_ctrl_sample=CTRL --sgrna_hdr=sgRNA --gene_hdr=gene --v10=gecko2
python run_mageck.py data/fiona_AML_raw_v10.tab data/fiona_AML_raw_v10_repmap.tab data/fiona_AML_raw_v10.tab --rep_hdr=Replicate --sample_hdr=Sample --common_ctrl_sample=CTRL --sgrna_hdr=sgRNA --gene_hdr=gene --v10=v10
python run_mageck.py data/whitehead/Wang2015_2017_merged_counts.txt data/whitehead/whitehead_repmap.txt data/whitehead/Wang2017_1-s2.0-S0092867417300612-mmc2.txt --sample_hdr=Sample --rep_hdr=Replicate --ctrl_sample_hdr=Control --sgrna_hdr=sgRNA_ID --gene_hdr=Symbol --v10=whitehead
python run_mageck.py data/tko/tko_counts.txt data/tko/tko_repmap.txt data/tko/tko_counts.txt --rep_hdr=Replicate --sample_hdr=Sample --ctrl_sample_hdr=Control --sgrna_hdr=GENE_CLONE --gene_hdr=GENE --ignore_blank_genes --v10=tko
    
----------------------------------------
**MAGeCK MLE Results**
----------------------------------------

python mageck_mle/run_mageck_mle.py data/yusa_v10/yusa_raw_v10.tab data/yusa_v10/yusa_raw_v10_repmap.tab:Replicate:Sample:CTRL data/yusa_v10/yusa_raw_v10.tab:sgRNA:gene v10_meanfc/v10_meanfc

----------------------------------------
**MeanFC Results**
----------------------------------------

python run_JACKS_meanfc.py data/fiona_AML_raw_v10.tab data/fiona_AML_raw_v10_repmap.tab data/fiona_AML_raw_v10.tab --outprefix=results/meanfc_jacks_results_v10/meanfc_jacks_results_v10 --rep_hdr=Replicate --sample_hdr=Sample --common_ctrl_sample=CTRL --sgrna_hdr=sgRNA --gene_hdr=gene --ctrl_genes=data/NEGv1.txt
python run_JACKS_meanfc.py data/Avana_sgrna_raw_readcounts_matched.csv data/Avana_replicatemap.csv data/Avana_sgrnamapping.csv --outprefix=results/meanfc_avana_results/meanfc_avana_jacks_results --rep_hdr=Replicate --sample_hdr=CellLine --common_ctrl_sample=pDNA --sgrna_hdr=Guide --gene_hdr=Gene --ctrl_genes=data/NEGv1.txt
python run_JACKS_meanfc.py data/Achilles_raw_GeckoV2.tab data/Achilles_raw_GeckoV2_repmap.tab data/Achilles_raw_GeckoV2.tab --outprefix=results/meanfc_results_gecko2/meanfc_results_gecko2 --rep_hdr=Replicate --sample_hdr=Sample --common_ctrl_sample=CTRL --sgrna_hdr=sgRNA --gene_hdr=gene --ctrl_genes=data/NEGv1.txt
python run_JACKS_meanfc.py data/whitehead/Wang2015_2017_merged_counts.txt data/whitehead/whitehead_repmap.txt data/whitehead/Wang2017_1-s2.0-S0092867417300612-mmc2.txt --outprefix=results/meanfc_results_whitehead/meanfc_results_whitehead --rep_hdr=Replicate --sample_hdr=Sample --ctrl_sample_hdr=Control --sgrna_hdr=sgRNA_ID --gene_hdr=Symbol --ctrl_genes=data/NEGv1.txt
python run_JACKS_meanfc.py data/tko/tko_counts.txt data/tko/tko_repmap.txt data/tko/tko_counts.txt --outprefix=results/meanfc_results_tko/meanfc_results_tko --rep_hdr=Replicate --sample_hdr=Sample --ctrl_sample_hdr=Control --sgrna_hdr=GENE_CLONE --gene_hdr=GENE --ctrl_genes=data/NEGv1.txt

----------------------------------------
**ScreenBEAM Results**
----------------------------------------

python run_screenbeam_all.py data/Avana_sgrna_raw_readcounts_matched.csv data/Avana_replicatemap.csv data/Avana_sgrnamapping.csv --rep_hdr=Replicate --sample_hdr=CellLine --common_ctrl_sample=pDNA --sgrna_hdr=Guide --gene_hdr=Gene --v10=avana
python run_screenbeam_all.py data/Achilles_raw_GeckoV2.tab data/Achilles_raw_GeckoV2_repmap.tab data/Achilles_raw_GeckoV2.tab --rep_hdr=Replicate --sample_hdr=Sample --common_ctrl_sample=CTRL --sgrna_hdr=sgRNA --gene_hdr=gene --v10=gecko2
python run_screenbeam_all.py data/fiona_AML_raw_v10.tab data/fiona_AML_raw_v10_repmap.tab data/fiona_AML_raw_v10.tab --rep_hdr=Replicate --sample_hdr=Sample --common_ctrl_sample=CTRL --sgrna_hdr=sgRNA --gene_hdr=gene --v10=v10
python run_screenbeam_all.py data/whitehead/Wang2015_2017_merged_counts.txt data/whitehead/whitehead_repmap.txt data/whitehead/Wang2017_1-s2.0-S0092867417300612-mmc2.txt --sample_hdr=Sample --rep_hdr=Replicate --ctrl_sample_hdr=Control --sgrna_hdr=sgRNA_ID --gene_hdr=Symbol --v10=whitehead
python run_screenbeam_all.py data/tko/tko_counts.txt data/tko/tko_repmap.txt data/tko/tko_counts.txt --rep_hdr=Replicate --sample_hdr=Sample --ctrl_sample_hdr=Control --sgrna_hdr=GENE_CLONE --gene_hdr=GENE --ignore_blank_genes --v10=tko

----------------------------------------
**PBNPA Results**
----------------------------------------

python run_PBNPA_all.py data/Avana_sgrna_raw_readcounts_matched.csv data/Avana_replicatemap.csv data/Avana_sgrnamapping.csv --rep_hdr=Replicate --sample_hdr=CellLine --common_ctrl_sample=pDNA --sgrna_hdr=Guide --gene_hdr=Gene --v10=avana
python run_PBNPA_all.py data/Achilles_raw_GeckoV2.tab data/Achilles_raw_GeckoV2_repmap.tab data/Achilles_raw_GeckoV2.tab --rep_hdr=Replicate --sample_hdr=Sample --common_ctrl_sample=CTRL --sgrna_hdr=sgRNA --gene_hdr=gene --v10=gecko2
python run_PBNPA_all.py data/fiona_AML_raw_v10.tab data/fiona_AML_raw_v10_repmap.tab data/fiona_AML_raw_v10.tab --rep_hdr=Replicate --sample_hdr=Sample --common_ctrl_sample=CTRL --sgrna_hdr=sgRNA --gene_hdr=gene --v10=v10
python run_PBNPA_all.py data/whitehead/Wang2015_2017_merged_counts.txt data/whitehead/whitehead_repmap.txt data/whitehead/Wang2017_1-s2.0-S0092867417300612-mmc2.txt --sample_hdr=Sample --rep_hdr=Replicate --ctrl_sample_hdr=Control --sgrna_hdr=sgRNA_ID --gene_hdr=Symbol --v10=whitehead
python run_PBNPA_all.py data/tko/tko_counts.txt data/tko/tko_repmap.txt data/tko/tko_counts.txt --rep_hdr=Replicate --sample_hdr=Sample --ctrl_sample_hdr=Control --sgrna_hdr=GENE_CLONE --gene_hdr=GENE --ignore_blank_genes --v10=tko

----------------------------------------
**Ceres Results**
----------------------------------------

Rscript ceres_avana.R
Rscript ceres_avana_mel.R
Rscript ceres_avana_col.R
Rscript ceres_avana_melcol.R
Rscript ceres_avana_rand.R
Rscript ceres_gecko2.R
Rscript ceres_gecko2.R

----------------------------------------
**JACKS Sampling Results**
----------------------------------------
for <cell_line, num_celllines, num_guides, num_replicates, avana_sample_outputfilename> in avana tests cases:
    avana_INPUT = ../../data/avana/Avana_sgrna_raw_readcounts_matched.csv#../../data/avana/Avana_replicatemap.csv:Replicate:CellLine:pDNA#../../data/avana/Avana_sgrnamapping.csv:Guide:Gene#
    python sample_jacks_screen.py avana_INPUT cell_line num_replicates(-1 for all) num_celllines(-1 for all) <avana_sample_outputfilename> samples num_guides(-1 for all)

for <cell_line, num_celllines, num_guides, num_replicates, gecko2_sample_outputfilename> in gecko2 tests cases:
    gecko2_INPUT = ../../data/gecko2/Achilles_raw_GeckoV2.tab#../../data/gecko2/Achilles_raw_GeckoV2_repmap.tab:Replicate:Sample:CTRL_pDNA_pXPR003_120K_201406#../../data/gecko2/Achilles_raw_GeckoV2.tab:sgRNA:gene#
    python sample_jacks_screen.py gecko2_INPUT cell_line num_replicates(-1 for all) num_celllines(-1 for all) <avana_sample_outputfilename> samples num_guides(-1 for all)

for <cell_line, num_celllines, num_guides, num_replicates, v10_sample_outputfilename> in yusa_v10 tests cases:
    v10_INPUT = data/yusa_v10/fiona_AML_raw_v10.tab#data/yusa_v10/fiona_AML_raw_v10_repmap.tab:Replicate:Sample:CTRL#data/yusa_v10/fiona_AML_raw_v10.tab:sgRNA:gene#
    python sample_jacks_screen.py v10_INPUT cell_line num_replicates(-1 for all) num_celllines(-1 for all) <avana_sample_outputfilename> samples num_guides(-1 for all)

for <cell_line, num_celllines, num_guides, num_replicates, tko_sample_outputfilename> in tko tests cases:
    tko_INPUT = ../../data/tko/tko_counts.txt#../../data/tko/tko_repmap.txt:Replicate:Sample:Control#../../data/tko/tko_counts.txt:GENE_CLONE:GENE#
    python sample_jacks_screen.py tko_INPUT cell_line num_replicates(-1 for all) num_celllines(-1 for all) <avana_sample_outputfilename> samples num_guides(-1 for all)

for <cell_line, num_celllines, num_guides, num_replicates, whitehead_sample_outputfilename> in whitehead tests cases:
    whitehead_INPUT = ../../data/whitehead/Wang2015_2017_merged_counts.txt#../../data/whitehead/whitehead_repmap.txt:Replicate:Sample:Control#../../data/whitehead/Wang2017_1-s2.0-S0092867417300612-mmc2.txt:sgRNA_ID:Symbol#CTRLGENE
    python sample_jacks_screen.py whitehead_INPUT cell_line num_replicates(-1 for all) num_celllines(-1 for all) <avana_sample_outputfilename> samples num_guides(-1 for all)



