#setwd("/home/jromahn/PHYTOARK_final_metabarcoding")

library(tidyverse)
library(stringr)

input_dir="02_mothur_MUC_Euka02"
output_dir="02_mothur_MUC_Euka02"

input_files <- list.files(input_dir, pattern="*wang.taxonomy", full.names=TRUE)

file <- input_files[2]


for (file in input_files){
  print(file)
  new_file <- gsub(".taxonomy", "__tax.tx", file)
  new_file2 <- gsub(".taxonomy", "__ORIGINAL_mothur.tx", file)
  
  data <- read.table(file, header=F, sep="\t")
  colnames(data) <- c("ID", "taxonomy")
  
  
  taxo_no <- str_count( data$taxonomy[1],";")
  data_tax <- data.frame()
  i <-3
  for (i in c(1: length(data$ID))){
    tax <- str_split( data$taxonomy[i],";")[[1]]
    tax <- gsub("_unclassified", "", tax) # remove because just copied taxon name which passes the threshold (in case of mothur)
    tax <- unique(tax)
    
    if (taxo_no == 8){
      #kingdom	supergroup	division	class	order	family	genus	species

      data_tax <- rbind(data_tax, data.frame(ID= data$ID[i], kingdom= tax[1],
                                           supergroup = tax[2],division= tax[3],class= tax[4],order= tax[5], family= tax[6],genus= tax[7], species= tax[8]))
    }else{ 
      if(taxo_no == 9){ # after pr2 version 5
        data_tax <- rbind(data_tax, data.frame(ID= data$ID[i], domain= tax[1],
                                             supergroup = tax[2],division =  tax[3],subdivision= tax[4],class= tax[5],order= tax[6], family= tax[7],genus= tax[8], species= tax[9]))
      }
      else{data_tax <- rbind(data_tax, data.frame(ID= data$ID[i], kingdom= tax[1],
                                             phylum= tax[2],class = tax[3], order= tax[4], family= tax[5],genus= tax[6], species= tax[7]))
      }
    }
  }
  write.csv(data_tax, file= new_file2, row.names = F)
  print("separate now")
  if (taxo_no == 8){
    data_tax <- data_tax %>% 
      separate(kingdom, c("kingdom", "sk_success"), sep ="\\(|\\)") %>%
      separate(supergroup, c("taxon", "cl_success"), sep ="\\(|\\)") %>%
      separate(division, c("division", "k_success"), sep ="\\(|\\)") %>%
      separate(class, c("class", "p_success"), sep ="\\(|\\)") %>%
      separate(order, c("order", "o_success"), sep ="\\(|\\)") %>%
      separate(family, c("family", "f_success"), sep ="\\(|\\)") %>%
      separate(genus, c("genus", "g_success"), sep ="\\(|\\)") %>%
      separate(species, c("species", "sp_success"), sep ="\\(|\\)") 
  }else{
    if(taxo_no == 9){
      data_tax <- data_tax %>% 
        separate(domain, c("domain", "dom_success"), sep ="\\(|\\)") %>%
        separate(supergroup, c("supergroup", "sk_success"), sep ="\\(|\\)") %>%
        separate(division, c("division", "subk_success"), sep ="\\(|\\)") %>%
        separate(subdivision, c("subdivision", "ph_success"), sep ="\\(|\\)") %>%
        separate(class, c("class", "cl_success"), sep ="\\(|\\)") %>%
        separate(order, c("order", "o_success"), sep ="\\(|\\)") %>%
        separate(family, c("family", "f_success"), sep ="\\(|\\)") %>%
        separate(genus, c("genus", "g_success"), sep ="\\(|\\)") %>%
        separate(species, c("species", "sp_success"), sep ="\\(|\\)") 
    }
    else{
      data_tax <- data_tax %>%
        separate(kingdom, c("kingdom", "sk_success"), sep ="\\(|\\)") %>%
        separate(phylum, c("taxon", "p_success"), sep ="\\(|\\)") %>%
        separate(class, c("class", "cl_success"), sep ="\\(|\\)") %>%
        separate(order, c("order", "o_success"), sep ="\\(|\\)") %>%
        separate(family, c("family", "f_success"), sep ="\\(|\\)") %>%
        separate(genus, c("genus", "g_success"), sep ="\\(|\\)") %>%
        separate(species, c("species", "sp_success"), sep ="\\(|\\)")
    }
    
   
    
  }
  data_tax[data_tax==""]<-NA
  write.csv(data_tax, file= new_file, row.names = F)
}



