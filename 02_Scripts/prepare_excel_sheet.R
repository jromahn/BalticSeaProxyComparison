rm(list = ls())
library("xlsx") # excel has to be closest properly
library(gt)
library(tidyverse)

# Set the Java heap size to 2 GB


setwd("/Users/Juliane/Documents/00_Work_SGN/00_PhytoArk/XX_PAPERS/2024_Comparison_Paper__relAbudance/")

# supplementary 1 - metabarcoding metadata
path="03_metadata_final"
filename="Mock_positive_control.csv"
table1 = read.table(file.path(path,filename), header=T, sep=",")
#colnames(table1)[3] <- "DNA concentration (ng/µL)"


# supplementary 2  compoistion positive control
path="03_metadata_final"
filename="EMB262_PhytoArk_MUC__detailed_metadata.tsv"
table2 = read.table(file.path(path,filename), header=T, sep="")
colnames(table2)[20] <- "TOC(%)"

# supplementary 3 range coordinates
path="03_metadata_final"
filename="Monitoring_range_of_coordinates.tsv"
table3 = read.table(file.path(path,filename), header=T, sep="")
colnames(table3) <- c("location", "core.location", "Latitude(DD)", "Longitude(DD)","maxLongitude(DD)", "minLongitude(DD)", "maxLatitude(DD)", "minLatitude(DD)" )

# supplementary 4 biomarker data
path="03_metadata_final"
filename="Supplementary_Biomarker_data_GoF_EGB.tsv"
table4 = read.table(file.path(path,filename), header=T, sep="")
colnames(table4)[3:6] <- c("Age(CE)","Dinosterol & Dinostanol (µg/gTOC)","C25:1 & C27:1 alkene (µg/gTOC)", "depth(cm)" )


# supplementary 5 mothur assignment
path="02_mothur_MUC_Euka02"
filename="8_PhytoArk_Euka02_2024_final_community_MUC__mothur_sequences.0_SSU_mothur__concatinated.wang__tax.tx"
table5 = read.table(file.path(path,filename), header=T, sep=",")
table5[] <- lapply(table5, function(x) if (is.factor(x)) as.character(x) else x)


# supplementary 6 cleaning statistics
path="02_cleaning"
filename="8_PhytoArk_euka_OT4_final__community_cleaned__Cleaning_Stats.csv"
table6 = read.table(file.path(path,filename), header=T, sep=",")

# supplementary 7 GAM stats
path="04_exploration_july24/GAM"
filename="Parameters_GAM_modeling.tsv"
table7 = read.table(file.path(path,filename), header=T, sep="")


# supplementary 8 correlation stats
path="04_exploration_july24"
filename="Correlation_results_as_table_of_gam_n_interval.tsv"
table8 = read.table(file.path(path,filename), header=T, sep="")

Filename="Supplementary_Tables_S1-8.xlsx"

# Write the first data set in a new workbook

write.xlsx(table1, file =Filename,
           sheetName = "S1 - PCR Positive control composition", append = FALSE, row.names = F)
# Add a second data set in a new worksheet
write.xlsx(table2, file = Filename, 
           sheetName="S2 - Metabarcoding & horizon metadata", append=TRUE, row.names = F)

write.xlsx(table3, file = Filename, 
           sheetName="S3 - Biomonitoring coordiante range", append=TRUE, row.names = F)
 
write.xlsx(table4, file = Filename, 
           sheetName="S4 - Biomarker data", append=TRUE, row.names = F)


write.xlsx(table6, file = Filename, 
           sheetName="S6 - Statistics for cleaning steps", append=TRUE, row.names = F)

write.xlsx(table7, file = Filename, 
           sheetName="S7 - Statistics of the GAM modeling", append=TRUE, row.names = F)

write.xlsx(table8, file = Filename, 
           sheetName="S8 - Statistics of the correlations", append=TRUE, row.names = F)


##########
path="03_metadata_final"
output_path=file.path(path, "tex")
if (!dir.exists(output_path)){ dir.create(output_path) }
table1 %>%  gt() %>%
  tab_header(title = md("Composition of the PCR Positive control ")) |> gtsave(filename =file.path(output_path,"Supplemetary_Table_S1.docx"))

colnames(table2)
table2 %>%  select(-c(latitude, longditude, forward_primerseq, reverse_primerseq, water.depth.m.bsl., 
                      latitude.DD.,longditude.DD.,extraction.date,volume.µL., agemodel.version, R1, R2,
                      numeric_age.CE., new_sample_id, sediment.subsample, sample_type))%>% 
  gt() %>%
  tab_header(title = md(" Metabarcoding & core sediment sample metadata")) |> gtsave(filename =file.path(output_path,"Supplemetary_Table_S2.docx"))

table3 %>%  gt() %>%
  tab_header(title = md("Range of Biomonitoring coordiantes")) |> gtsave(filename =file.path(output_path,"Supplemetary_Table_S3.docx"))

table4 %>%  gt() %>%
  tab_header(title = md("Biomarker data")) |> gtsave(filename =file.path(output_path,"Supplemetary_Table_S4.docx"))

table6 %>%  gt() %>%
  tab_header(title = md("Statistics for cleaning steps")) |> gtsave(filename =file.path(output_path,"Supplemetary_Table_S6.docx"))

table7 %>% select(- weights)%>% gt() %>%
  tab_header(title = md("Statistics of the GAM modeling")) |> gtsave(filename =file.path(output_path,"Supplemetary_Table_S7.docx"))

table8 %>% select(- stats_method)%>% gt() %>%
  tab_header(title = md("Statistics of the Spearman correlations")) |> gtsave(filename =file.path(output_path,"Supplemetary_Table_S8.docx"))

