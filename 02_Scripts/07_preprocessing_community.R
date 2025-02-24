rm(list = ls())
library(tidyverse)
library(vegan)
library(phylotools)



setwd("/Users/Juliane/Documents/00_Work_SGN/00_PhytoArk/XX_PAPERS/2025_Comparison_Paper__relAbudance/")
print("Preproccessing starts")

# output data
output_path <- "02_intermediate"
if (!dir.exists(output_path)){ dir.create(output_path) }

## community
path <- "02_intermediate/"
community_file <- "8_PhytoArk_euka_OT4_final__community_cleaned__community_data.RData"

## taxonomy
path_mothur <- "02_mothur_MUC_Euka02"
file_taxo    <- "8_PhytoArk_Euka02_2024_final_community_MUC__mothur_sequences.0_SSU_mothur__concatinated.wang__tax.tx" # mothur taxonomy file
fasta_file <- "8_PhytoArk_Euka02_2024_final_community_MUC__mothur_sequences.fasta" #mothur fasta file

#metadata
meta_path <-"03_metadata_final"
meta_file <- "EMB262_PhytoArk_MUC__detailed_metadata.tsv"
meta_file_sample <- "EMB262_PhytoArk_MUC__detailed_metadata_samples.tsv"

meta_data <- read.table(file.path(meta_path, meta_file), sep="\t", header=T)%>%
  rename("station"="core.location")%>%
  rename("depth"="core.depth.cm.")%>% 
  filter(!grepl("Extraction", core.samplename)) %>%
  mutate(sample_id=gsub("_(R|P)\\d+",  "",new_sample_id, perl=T))%>%
  rename(weight= 'weight.µg.')%>%
  rename(age= 'age.CE.')
write.table(meta_data, file.path(output_path, meta_file_sample), sep="\t", row.names = F)

##read input
seda_taxo <- read.table( file.path(path_mothur, file_taxo), header= T, sep =",")
load(file.path(path, community_file))
seda_community <- final_data
rm(final_data)

# remove asvs which were removed because of the cleaning
asv_list <- names(seda_community)[grep("ASV", names(seda_community))]
seda_taxo <- seda_taxo %>% filter(ID %in% asv_list)%>%
  filter(domain == "Eukaryota")
save(seda_taxo, file=file.path(output_path,"sedaDNA_taxonomy_eukaryotes_cleaned.RData"))
write.csv(seda_taxo, file=file.path(output_path,"sedaDNA_taxonomy_eukaryotes_cleaned.csv"), row.names = F)

#remove ASVs which were removed

seda_community <- seda_community %>% select( sample_id, tag, station, depth, replicate, one_of(seda_taxo$ID))
#seda_community <- cbind(seda_community[,c("sample_id", "tag", "station", "depth","replicate")],   seda_community[,names(seda_community) %in% seda_taxo$ID ])


## combine with meta data
seda_community <- merge(meta_data[,c("tag", "sample_id", "weight", "age")], seda_community, by =c("sample_id", "tag"))

## add location
print(unique(seda_community$station))
seda_community <- seda_community %>% mutate(location = ifelse(station == "EMB262_6_28_MUC", "Eastern Gotland Basin",
                                                       ifelse(station == "EMB262_13_8_MUC","Landsort Deep", "Gulf of Finland")))

## normalize to 500 mg
seda_community[,grep("ASV", names(seda_community))] <-seda_community[,grep("ASV", names(seda_community))]/  as.numeric(seda_community$weight) *500

save(seda_community, file=file.path(output_path,"sedaDNA_community_eukaryotes_normed500.RData"))


############# combine to reduce data
print("Preproccessing starts -- merging taxo & community")
seda_community_sum <-seda_community %>% 
          pivot_longer(!c("sample_id","tag", "weight","station", "location","depth", "age", "replicate" ), names_to = "ID", values_to = "reads")%>%
          left_join(seda_taxo, by ="ID")%>%
          group_by(location, station,  sample_id, replicate, tag ,depth, age, division, subdivision,class, order, family, genus, species )%>%
          dplyr::summarise(ASV_no= n(), reads= sum(reads))

rm(seda_community)
save(seda_community_sum, file=file.path(output_path,"sedaDNA_community_eukaryotes_normed500_transformed.RData"))
rm(seda_community_sum)

################################### unassigned ASVs
# extract asv which are unassigned
print("Preproccessing starts -- unassigned")
unassigned_subdivision <-  seda_taxo %>% filter(is.na(subdivision)) %>%pull(ID)

# fasta and keep unassigned fasta
sequences <- read.fasta(file.path(path_mothur,fasta_file))
sequences_unassigned <- sequences %>% filter(seq.name %in% unassigned_subdivision)

# write fasta
dat2fasta(sequences_unassigned, outfile = file.path(output_path,"02__unassigned_MUC_ASVs__subidivision_unknown.fasta"))
