DEPRECATION WARNING: a command line interface for run_jacks.py changed, one should update commands below according to https://github.com/felicityallen/JACKS/blob/master/jacks/README.md to run them

---------------------------------------------
Data and Results: Figshare location
--------------------------------------------

Results files can be found at:

https://figshare.com/articles/Results/6002438

Data files can be found at:

yusa_v10:

https://figshare.com/articles/Yusa_V1_0_Data/6002402

gecko2:

https://figshare.com/articles/GeCKOv2_Data/6002408

avana:

https://figshare.com/articles/Avana_Data/6002417

Instructions for re-generating results can be found at:

https://figshare.com/articles/ReadMe_-_How_to_reproduce_results/6002444

To generate results (python and R files in scripts or jacks directory):

----------------------------------------
**JACKS full results**
----------------------------------------

python run_JACKS.py data/avana/Avana_sgrna_raw_readcounts_matched.csv data/avana/Avana_replicatemap.csv:Replicate:CellLine:pDNA data/Avana_sgrnamapping.csv:Guide:Gene avana_jacks/avana_jacks 0
python run_JACKS.py data/avana/Avana_sgrna_raw_readcounts_rand.csv data/avana/Avana_replicatemap_rand.csv:Replicate:CellLine:pDNA data/Avana_sgrnamapping.csv:Guide:Gene avana_jacks_rand/avana_jacks_rand 0
python run_JACKS.py data/avana/Avana_sgrna_raw_readcounts_rand.csv data/avana/Avana_replicatemap_rand.csv:Replicate:CellLine:pDNA data/Avana_sgrnamapping.csv:Guide:Gene avana_jacks_rand_hp/avana_jacks_rand_hp 1

python run_JACKS.py data/gecko2/Achilles_raw_GeckoV2.tab data/gecko2/Achilles_raw_GeckoV2_repmap.tab:Replicate:Sample:CTRL data/gecko2/Achilles_raw_GeckoV2.tab:sgRNA:gene gecko2_jacks/gecko2_jacks 0
python run_JACKS.py data/gecko2/Achilles_raw_GeckoV2_rand.tab data/gecko2/Achilles_raw_GeckoV2_repmap_rand.tab:Replicate:Sample:CTRL data/gecko2/Achilles_raw_GeckoV2.tab:sgRNA:gene gecko2_jacks_rand/gecko2_jacks_rand 0
python run_JACKS.py data/gecko2/Achilles_raw_GeckoV2_rand.tab data/gecko2/Achilles_raw_GeckoV2_repmap_rand.tab:Replicate:Sample:CTRL data/gecko2/Achilles_raw_GeckoV2.tab:sgRNA:gene gecko2_jacks_rand_hp/gecko2_jacks_rand_hp 1

python run_JACKS.py data/yusa_v10/yusa_raw_v10.tab data/yusa_v10/yusa_raw_v10_repmap.tab:Replicate:Sample:CTRL data/yusa_v10/yusa_raw_v10.tab:sgRNA:gene v10_jacks/v10_jacks 0

----------------------------------------
**JACKS x results**
----------------------------------------

for <num_celllines> in test cases:
    python scripts/jacks/sample_jacks_screen_xs.py Avana_raw_reformatted_Batch0.txt 2 <num_celllines> avana0_xs/xs_2_<num_celllines>.txt 100
    python scripts/jacks/sample_jacks_screen_xs.py Avana_raw_reformatted_Batch1.txt 2 <num_celllines> avana1_xs/xs_2_<num_celllines>.txt 100

----------------------------------------
**JACKS Colon/Melanoma BRAF Results**
----------------------------------------

python run_JACKS.py data/avana/Avana_sgrna_raw_readcounts_matched.csv data/avana/Avana_replicatemap_colon.csv:Replicate:CellLine:pDNA data/Avana_sgrnamapping.csv:Guide:Gene avana_jacks_col/avana_jacks_col 0
python run_JACKS.py data/avana/Avana_sgrna_raw_readcounts_matched.csv data/avana/Avana_replicatemap_melanoma.csv:Replicate:CellLine:pDNA data/Avana_sgrnamapping.csv:Guide:Gene avana_jacks_mel/avana_jacks_mel 0
python run_JACKS.py data/avana/Avana_sgrna_raw_readcounts_matched.csv data/avana/Avana_replicatemap_colonmelanoma.csv:Replicate:CellLine:pDNA data/Avana_sgrnamapping.csv:Guide:Gene avana_jacks_melcol/avana_jacks_melcol 0

----------------------------------------
**JACKS Sampling Results**
----------------------------------------

for <cell_line, num_celllines, num_guides, num_replicates, avana_sample_outputfilename> in avana tests cases:
    python scripts/jacks/sample_jacks_screen.py data/avana/Avana_raw_reformatted.tab cell_line num_replicates(-1 for all) num_celllines(-1 for all) <avana_sample_outputfilename> samples num_guides(-1 for all)
for <cell_line, num_celllines, num_guides, num_replicates, gecko2_sample_outputfilename> in gecko2 tests cases:
    python scripts/jacks/sample_jacks_screen.py data/gecko2/Achilles_raw_GeckoV2.tab cell_line num_replicates(-1 for all) num_celllines(-1 for all) <avana_sample_outputfilename> samples num_guides(-1 for all)
for <cell_line, num_celllines, num_guides, num_replicates, avana_sample_outputfilename> in yusa_v10 tests cases:
    python scripts/jacks/sample_jacks_screen.py data/yusa_v10/yusa_v10_raw.tab cell_line num_replicates(-1 for all) num_celllines(-1 for all) <avana_sample_outputfilename> samples num_guides(-1 for all)

----------------------------------------
**JACKS Single Line Results**
----------------------------------------

python scripts/jacks/run_JACKS_single.py data/avana/Avana_sgrna_raw_readcounts.csv data/avana/Avana_replicatemap.csv:Replicate:CellLine:pDNA data/Avana_sgrnamapping.csv:Guide:Gene avana_single_jacks/avana_single_jacks
python scripts/jacks/combine_single_jacks.py avana_single_jacks

python scripts/jacks/run_JACKS_single.py data/yusa_v10/yusa_raw_v10.tab data/yusa_v10/yusa_raw_v10_repmap.tab:Replicate:Sample:CTRL data/yusa_v10/yusa_raw_v10.tab:sgRNA:gene v10_single_jacks/v10_single_jacks/v
python scripts/jacks/combine_single_jacks.py v10_single_jacks

python scripts/jacks/run_JACKS_single.py data/gecko2/Achilles_raw_GeckoV2.tab data/gecko2/Achilles_raw_GeckoV2_repmap.tab:Replicate:Sample:CTRL data/gecko2/Achilles_raw_GeckoV2.tab:sgRNA:gene gecko2_single_jacks/gecko2_single_jacks
python scripts/jacks/combine_single_jacks.py gecko2_single_jacks

----------------------------------------
**JACKS Reference Leave-One-Out Results**
----------------------------------------

for <idx, cell_line> in avana:
    python scripts/jacks/leaveoneout_reference_test.py <cell_line> data/avana/Avana_sgrna_raw_readcounts_matched.csv data/avana/Avana_replicatemap.csv:Replicate:CellLine:pDNA data/avana/Avana_sgrnamapping.csv:Guide:Gene avana_loor/avana_loor_<idx>

for <idx, cell_line> in gecko2:
    python scripts/jacks/leaveoneout_reference_test.py <cell_line> data/gecko2/Achilles_raw_GeckoV2.tab data/gecko2/Achilles_raw_GeckoV2_repmap.tab:Replicate:Sample:CTRL data/gecko2/Achilles_raw_GeckoV2.tab:sgRNA:gene gecko2_loor/gecko2_loor_<idx>
    
for <idx, cell_line> in yusa_v10:
    python scripts/jacks/leaveoneout_reference_test.py <cell_line> data/yusa_v10/yusa_raw_v10.tab data/yusa_v10/yusa_raw_v10_repmap.tab:Replicate:Sample:CTRL data/yusa_v10/yusa_raw_v10.tab:sgRNA:gene v10_loo/v10_loo_<idx>

----------------------------------------
**Bagel Results**
----------------------------------------

python BAGEL-calc_foldchange.py -i ../data/avana/Avana_raw_reformatted.tab -o ../../data/avana/Avana_raw_reformatted.foldchange -c"[1,2,3,4]"
for <cell_line> in avana:
    python run_bagel.py ../../data/avana/Avana_raw_reformatted.foldchange <cell_line> -1 avana_bagel/avana_bagel_<cell_line>
python collate_bagel_results.py avana_bagel
    
python BAGEL-calc_foldchange.py -i ../data/gecko2/Achilles_raw_GeckoV2.tab -o ../../data/gecko2/Achilles_raw_GeckoV2.foldchange -c"[1]"
for <cell_line> in gecko2:
    python run_bagel.py ../../data/gecko2/Achilles_raw_GeckoV2.foldchange <cell_line> -1 gecko2_bagel/gecko2_bagel_<cell_line>
python collate_bagel_results.py gecko2_bagel

python BAGEL-calc_foldchange.py -i ../data/yusa_v10/yusa_raw_v10.tab -o ../../data/yusa_v10/yusa_raw_v10.foldchange -c"[1,2]"
for <cell_line> in v10:
    python run_bagel.py ../../data/yusa_v10/yusa_raw_v10.foldchange <cell_line> -1 avana_bagel/avana_bagel_<cell_line>
python collate_bagel_results.py v10_bagel    
    
----------------------------------------
**MAGeCK Results**
----------------------------------------
for <cell_line> in avana:
    python run_mageck.py ../../data/avana/Avana_raw_reformatted.tab <cell_line> -1 avana_mageck/avana_mageck_<cell_line>
python collate_mageck_results.py avana_mageck
      
for <cell_line> in gecko2:
    python run_mageck.py ../../data/gecko2/Achilles_raw_GeckoV2.tab <cell_line> -1 gecko2_mageck/gecko2_mageck_<cell_line>
python collate_mageck_results.py gecko2_mageck

for <cell_line> in v10:
    python run_mageck.py ../../data/yusa_v10/yusa_raw_v10.tab <cell_line> -1 avana_mageck/avana_mageck_<cell_line>
python collate_mageck_results.py v10_mageck    

----------------------------------------
**MAGeCK MLE Results**
----------------------------------------

python mageck_mle/run_mageck_mle.py data/avana/Avana_sgrna_raw_readcounts_matched.csv data/avana/Avana_replicatemap.csv:Replicate:CellLine:pDNA data/Avana_sgrnamapping.csv:Guide:Gene avana_meanfc/avana_meanfc
python mageck_mle/run_MeanFC.py data/yusa_v10/yusa_raw_v10.tab data/yusa_v10/yusa_raw_v10_repmap.tab:Replicate:Sample:CTRL data/yusa_v10/yusa_raw_v10.tab:sgRNA:gene v10_meanfc/v10_meanfc
python mageck_mle/run_MeanFC.py data/gecko2/Achilles_raw_GeckoV2.tab data/gecko2/Achilles_raw_GeckoV2_repmap.tab:Replicate:Sample:CTRL data/gecko2/Achilles_raw_GeckoV2.tab:sgRNA:gene gecko2_meanfc/gecko2_meanfc

----------------------------------------
**MeanFC Results**
----------------------------------------

python meanfc/run_MeanFC.py data/avana/Avana_sgrna_raw_readcounts_matched.csv data/avana/Avana_replicatemap.csv:Replicate:CellLine:pDNA data/Avana_sgrnamapping.csv:Guide:Gene avana_meanfc/avana_meanfc
python meanfc/run_MeanFC.py data/yusa_v10/yusa_raw_v10.tab data/yusa_v10/yusa_raw_v10_repmap.tab:Replicate:Sample:CTRL data/yusa_v10/yusa_raw_v10.tab:sgRNA:gene v10_meanfc/v10_meanfc
python meanfc/run_MeanFC.py data/gecko2/Achilles_raw_GeckoV2.tab data/gecko2/Achilles_raw_GeckoV2_repmap.tab:Replicate:Sample:CTRL data/gecko2/Achilles_raw_GeckoV2.tab:sgRNA:gene gecko2_meanfc/gecko2_meanfc

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


