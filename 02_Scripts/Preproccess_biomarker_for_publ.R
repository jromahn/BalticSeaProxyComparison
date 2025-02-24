rm(list=ls())

library(readxl)
library(tidyverse)


path = "/Users/Juliane/Documents/00_Work_SGN/00_PhytoArk/"
output_path= "/Users/Juliane/Documents/00_Work_SGN/00_PhytoArk/XX_PAPERS/2024_Comparison_Paper__relAbudance/03_metadata_upload"
setwd(path)
## Biomarker
biomarker_path <- paste(path,"04_biomarker_jerome", sep="/")

#GoF
biomarker_file.1 <- "EMB262_12-2MUC1_hydrocarbons_for Juliane.xlsx"
biomarker_file.2 <- "EMB262_12-2MUC1_sterols_alkanols_FAs_for Juliane.xlsx"

#EGB
biomarker.file.5 <- "lipid data MUC EGB & Färö.xlsx"


### ## ## ## ## ## ## ## ## ## ## ## Biomarker
## Join Biomarker data for  GOF
carbon_file=file.path(biomarker_path,biomarker_file.1)
carbon_data <- read_excel(carbon_file)
carbon_data$age <- round(carbon_data$age)

carbon_file=file.path(biomarker_path,biomarker_file.2)
sterol_data <- read_excel(carbon_file)

colnames(sterol_data) <- gsub("Results","(ng/gTOC)", colnames(sterol_data))
sterol_data$age <- round(sterol_data$age)
sterol_data[,c("Name", "Data File")] <- NULL

biomarker <- merge(sterol_data, carbon_data, by =c("age","depth"), all = T)

biomarker <- biomarker[,c("age","depth","C25:1 (ng/gTOC)",  "C27:1 (ng/gTOC)","Dinosterl+dinostanols (ng/gTOC)")]


rm(sterol_data,carbon_data)
save(biomarker, file=file.path(output_path,"biomarker_data_GoF.RData"))


############## egb
egb_lipid_file=file.path(biomarker_path,biomarker.file.5)
lipid_data <- read_excel(egb_lipid_file, sheet =3)

# transform data from both location to same unit
biomarker[,grep("ng/gTOC", colnames(biomarker))] <-  biomarker[,grep("ng/gTOC", colnames(biomarker))] / 1000
colnames(biomarker) <- gsub("ng/gTOC", "µg/gTOC",colnames(biomarker))
biomarker <- biomarker %>% 
  mutate( location = rep("Gulf of Finland", length(age)), 
          `C25:1+C27:1 (µg/gTOC)` = `C25:1 (µg/gTOC)`+ `C27:1 (µg/gTOC)`,
          `C25:1 (µg/gTOC)`= NULL ,`C27:1 (µg/gTOC)` = NULL)%>%
  mutate(core = "EMB262/6-28")%>%
  rename("depth (cm)"=depth)


lipid_data <-lipid_data %>% 
  dplyr::rename("age" = `age (yr CE)`, "Dinosterl+dinostanols (µg/gTOC)"= `Dino_sterol+stanol (µg/gTOC)`) %>% 
  mutate( location = "Eastern Gotland Basin")%>%
  mutate(core = "EMB201/7-4")

egb_gof_biomarker <- biomarker %>%
                        bind_rows(lipid_data)%>%
                        select(location,core, age, `Dinosterl+dinostanols (µg/gTOC)`,`C25:1+C27:1 (µg/gTOC)`,`depth (cm)` )

colnames(lipid_data)

write.table(egb_gof_biomarker, file=file.path(output_path,"biomarker_data_GoF_EGB.tsv"), row.names = F, sep = "\t")


