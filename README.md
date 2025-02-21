# BalticSeaProxyComparison

## Dependencies

```R
library(tidyverse)
library(vegan)
library(knitr)
library(phylotools)

# modeling
library(caret)
library(mgcv)

# taxonomy
library(taxonomizr) # assigning NCBI ID to taxonomy
library(worrms)

#plotting
library(ggpubr)
library(ggExtra)

#color packages
library(nord)
library(Manu)
```
## Necessary metadata in FigShare

| Filname in R Script  | FigShare DOI |
| ------------- | ------------- |
| EMB262_PhytoArk_MUC__detailed_metadata.tsv  | 10.6084/m9.figshare.28457471  |
| Supplementary_Biomarker_data_GoF_EGB.tsv  | 10.6084/m9.figshare.28457492  |
| 8_PhytoArk_Euka02_2024_final__community__MUC.RData | 10.6084/m9.figshare.28457474 | 
| 8_PhytoArk_Euka02_2024_final_community_MUC__mothur_sequences.fasta | 10.6084/m9.figshare.28457495 |
| 8_PhytoArk_Euka02_2024_final_community_MUC__mothur_sequences.0_SSU_mothur__concatinated.wang__tax.tx | 10.6084/m9.figshare.28457483 |


## Script Overview 

```bash
00_R_functions.R
05_Create_SpeciesTable.R
06_new_cleaning_Euka02.Rmd
07_preproccessing_monitoring.R
07_preprocessing_community.R
07_script_unassigned.R
08_MAIN_analysis.Rmd
09_fitting_GAM_and_statistics.Rmd
```


### 00_R_functions.R

R script containing most relevant R functions of the various scripts

### 05_Create_SpeciesTable.R

Input:
- mothur output file from classify.seqs() function

Output:
- 8_PhytoArk_Euka02_2024_final__mothur_sequences.0_SSU_mothur.wang__tax.tx 

###  06_new_cleaning_Euka02.Rmd

Input: 
- EMB262_PhytoArk_MUC__detailed_metadata.tsv
- 8_PhytoArk_Euka02_2024_final__community__MUC.RData

Output:
- 8_PhytoArk_euka_OT4_final__community_cleaned_community_data.RData # community matrix in RData format
- 8_PhytoArk_euka_OT4_final__community_cleaned_matrix.csv           # community matrix in csv format
- 8_PhytoArk_euka_OT4_final__community_cleaned_Cleaning_Stats.csv   # overview ASV and read number statistics for every cleaning step


###  07_preproccessing_monitoring.R

Input: 
- DomePhytoplankton_Data_0420454538.csv

Output:
- ICES_phytodata_with_worms_EMB262MUCs__SprSum.RData
- ICES_phytodata_with_worms_EMB262MUCs.RData
- Suppl_ICES__Month_of_Monitoring_events__march-july__5perc.png
- Suppl_ICES__Month_of_Monitoring_events__april-june__5perc.png

###  07_preprocessing_community.R

Input: 
- 8_PhytoArk_euka_OT4_final__community_cleaned__community_data.RData           #(06_new_cleaning_Euka02.Rmd)
- 8_PhytoArk_Euka02_2024_final__mothur_sequences.0_SSU_mothur.wang__tax.tx     #(05_Create_SpeciesTable.R)
- 8_PhytoArk_Euka02_2024_final__mothur_sequences.fasta                         # mothur fasta input file
- EMB262_PhytoArk_MUC__detailed_metadata.tsv

Output:
- EMB262_PhytoArk_MUC__detailed_metadata_samples.tsv
- sedaDNA_taxonomy_eukaryotes_cleaned.csv
- sedaDNA_taxonomy_eukaryotes_cleaned.RData
- sedaDNA_community_eukaryotes_normed500.RData
- sedaDNA_community_eukaryotes_normed500_transformed.RData
- 02__unassigned_MUC_ASVs__subidivision_unknown.fasta

### 07_script_unassigned.R 

Input: 
- sedaDNA_community_eukaryotes_normed500.RData                               #(07_preprocessing_community.R)
- Euka_Unassigned_Blast_nt240122-1-ex.txt                                    # Megan output

Output:
- BLASTn_NCBI_reassignment_results.png

### 08_MAIN_analysis.Rmd

Input: 
- ICES_phytodata_with_worms_EMB262MUCs__SprSum.RData                         #(07_preproccessing_monitoring.R)
- sedaDNA_community_eukaryotes_normed500_transformed.RData                   #(07_preprocessing_community.R)
- EMB262_MUC_detailed_renaming_metadata.csv
- Supplementary_Biomarker_data_GoF_EGB.tsv

Output:
- MUC__TOC_values_over_time.png
- Monitoring__dinoflagellate_community_composition__orderLevel_yearly.png
- Metabarcoding_community_comp_over_time__interval.png
- community_comp_over_time__interva.tsv
- Metabarocding__community_comp_over_time__interval_abuClass.png
- community_comp_phytoplankton_ALL.png
- Data_allApproach__for_modelling.RData
- Data_allApproach__for_modelling__real.RData
- Supplementary__metabarcoding_results.png

### 09_fitting_GAM_and_statistics.Rmd

Input: 
- Data_allApproach__for_modelling.RData                                      #(08_MAIN_analysis.Rmd)
  
