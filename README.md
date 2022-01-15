# larva_scSeq2022

This repo contains most of the information needed to generate the Figures of the paper "A complete series of temporal transcriptions factors in the fly visual system". 

It contains:

- R files for the preparation of the Drosophila and human datasets, as well as a re-analysis of just the neuroepithelial cells in the Drosophila dataset. 
1_preparation_of_larval_pupal_dataset.R
1a_subsetting_NE_and_progenitors.R
2_preparation_of_human_data.R

- Importantly,  all the resulting files from the previous R scripts that are used later for Figures are also provided.
human_snSeq.Rdata
larvaOL.integrated150.rds
medulla_ne_reanalysis.rds
OL.combined_new.annotations.rds (IMPORTANT NOTE: This file contains pupal cluster annotations that were added semi-manually, so it is enriched compared to the file that is generated with the 1_preparation_of_larval_pupal_dataset.R. These annotations can be found here: OL.combined_new.annotations@meta.data$Nikos_Neset_ID3 and they are the most updated up to 02/10/2021.)

- R files that contain the commands used to generate the different Figures.
3_Figure1.R
4_Figure2_S3.R
5_Figure4_5_S7.R
6_Figure6_S10.R
7_FigureS1.R
8_FigureS2.R
9_FigureS11.R

- Supporting functions that are used in the above scripts:
classifier_utils.R
functions.R
LayeredFeaturePlot_white_background_pt.brightness.R

- Tables that are used in the different scripts:
ConversionTable.csv	: Table from Flybase that converts Drosophila gene names to FGgn.
FlyBase_TFs.csv		: All Drosophila TFs, as indicated in Flybase
mart_export_human.txt	: A table that contains the GO terms for all human genes
mart_export.txt		: A table that contains the GO terms for all Drosophila genes
colors.txt		: A vector that color codes the pupal clusters according to their predicted temporal window
colors2.txt		: A vector that color codes non-medulla clusters used in EDF 9.

- Datasets that are used in this project that are available at GEO:
GSE118953_raw_count.tsv	: mouse cortical  scSeq (Telley et al, Science 2019)
P15_10Xv3Adjusted.rds	: scSeq data from Drosophila optic lobes at P15 stage
P30_10Xv3Adjusted.rds	: scSeq data from Drosophila optic lobes at P30 stage
P40_10Xv3Adjusted.rds	: scSeq data from Drosophila optic lobes at P40 stage
P50_10Xv3Adjusted.rds	: scSeq data from Drosophila optic lobes at P50 stage
P70_10Xv3Adjusted.rds	: scSeq data from Drosophila optic lobes at P70 stage
