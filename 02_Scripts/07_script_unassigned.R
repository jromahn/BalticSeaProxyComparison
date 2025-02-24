rm(list = ls())
library(tidyverse)
library(vegan)
library(taxonomizr)


setwd("/Users/Juliane/Documents/00_Work_SGN/00_PhytoArk/XX_PAPERS/2024_Comparison_Paper__relAbudance/")
source("02_Scripts/00_R_functions.R") # read file with functions

output_path <- "04_exploration_july24"
output_path_extra <- file.path(output_path, "suppl")

### load community data
input_path="03_preprocessed_MUC"
load(file.path(input_path,"sedaDNA_community_eukaryotes_normed500.RData"))
     
##load megan data
input_path="04_unassigned"
megan_results<- read.table(file = file.path(input_path, "Euka_Unassigned_Blast_nt240122-1-ex.txt"), sep=",")
colnames(megan_results) <- c("ID", "taxID")

###filter for community data for euka unassigned ASVs
seda_community <- seda_community %>%
  mutate(totalReads = rowSums(select(., starts_with("ASV")), na.rm = TRUE))%>%
  select(-starts_with("ASV"), any_of(megan_results$ID))

##load NCBI taxonomy
unique_NCBIids <- unique(megan_results$taxID)


### get taxonomy
#execute if not downloaded yet
loc_ncbi <- getNamesAndNodes(outDir = ".",url = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz",fileNames = c("names.dmp", "nodes.dmp"))
sqlFile <- 'accessionTaxa.sql'
read.accession2taxid(list.files(loc_ncbi),sqlFile,indexTaxa=TRUE,overwrite=TRUE)
taxaNames<-read.names.sql(loc_ncbi[1],sqlFile)
taxaNodes<-read.nodes.sql(loc_ncbi[2],sqlFile)


ncbi_taxonomy <- function(unique_assigned){
  unique_assigned <- unique(unique_assigned)
  rank_information <- getTaxonomy(unique_assigned, sqlFile, 
                                  desiredTaxa = c("species","genus","family","order","class", "subphylum" ,"phylum", "kingdom","superkingdom"))
  rank_information <- data.frame(rank_information)
  rank_information$taxID <- unique_assigned
  return(rank_information)
}

unassigned_NCBI_taxonomy <- ncbi_taxonomy(unique_NCBIids)%>%select(superkingdom, taxID)%>%
                                right_join(megan_results, by = c("taxID"="taxID"))

### combine taxonomy and community data
seda_community <-seda_community %>%
  select(any_of(c("location", "tag", "age", "totalReads")),starts_with("ASV"))%>%
  pivot_longer(!c(location, tag, age, totalReads), names_to = "ID", values_to = "reads")%>%
  mutate(relAbu= reads/totalReads)%>%
  right_join(unassigned_NCBI_taxonomy, by="ID")%>%
  group_by(location, age, ID, superkingdom, taxID)%>%
  summarise(relAbu=mean(relAbu)*100)%>%
  ungroup()
  

### visualize

# stacked area chart
ggplot(seda_community, aes(x=age, y=relAbu, fill=superkingdom)) + 
  geom_area()+
  facet_wrap(~location, scale="free_x")+ 
  labs(title="Blastn reassignment results of unassigned eukaryotic amplicons",x= "Years CE", y= "Relative Read Abundance(%)", fill= "NCBI Superkingdom")+ 
  figure_theme
ggsave(file = file.path(output_path_extra, "BLASTn_NCBI_reassignment_results.png"), width=14, height=8)
