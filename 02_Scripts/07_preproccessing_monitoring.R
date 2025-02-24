rm(list=ls())

library(tidyverse)
library("worrms")
library(tidyverse)

### ## ## ## ## ## ## ## ## ## ## ## # Biomonitoring
setwd("/Users/juliane/Documents/00_Work_SGN/00_PhytoArk/XX_PAPERS/2025_Comparison_Paper__relAbudance/")
source("02_Scripts/00_R_functions.R") # read file with functions

#input original ices data
path = "/Users/juliane/Documents/00_Work_SGN/00_PhytoArk/00_ICES/ICES_Phytoplankton"
ices_file <- "DomePhytoplankton_Data_0420454538.csv"

## output path
output_path <-"02_intermediate"
metadata_path <- "02_intermediate"
output_path_extra <- "02_intermediate"
if (!dir.exists(output_path)){ dir.create(output_path) }
if (!dir.exists(metadata_path)){ dir.create(metadata_path) }


#### core data
MUC_data <- data.frame(location=c("Gulf of Finland"	,"Landsort Deep",	"Eastern Gotland Basin"),
                       short=c("EMB262_12_2_MUC",	"EMB262_13_8_MUC",	"EMB262_6_28_MUC"	),
                       latt =c("59°34.443’","58°38.391’" , "57°17.004'"),
                       long = c("023°36.461’", "018°15.997’","020°07.244'"))

## transform MUC data to digital degree
MUC_data <-MUC_data %>% 
  separate(latt, c('degree', 'minutes', 'sec')) %>%
  mutate(degree= as.numeric(degree), minutes=as.numeric(paste(minutes, sec,sep=".")))%>%
  mutate(latt =round(degree +minutes/60 ,8)) %>% 
  separate(long, c('degree', 'minutes', 'sec')) %>%
  mutate(degree= as.numeric(degree), minutes=as.numeric(paste(minutes, sec,sep=".")))%>%
  mutate(long =round(degree +minutes/60 ,8)) %>%
  select( -degree, -minutes,  -sec)




## Load all and just keep in area from MUC and counting data
ices_data <- read.table(file.path(path,ices_file), header= T, sep ="'")
col_interest <- c("tblSampleID", "STATN", "HELCOM_subbasin", "HELCOM_L4", "Latitude", "Longitude", "MNDEP", "MXDEP", "MYEAR", "DATE", "Year", "Month", "Day",
                  "RLIST", "SPECI", "SPECI_name", "AphiaID", "WoRMS_name","AphiaID_accepted","WoRMS_accepted_name", "SIZCL", "SIZRF", "TRPHY", "STAGE",
                  "Value", "PARAM", "PARAM_desc","MUNIT", "final_value", "final_value_unit" )

ices_data <- ices_data[,col_interest]

#filter for count data
ices_data <- ices_data%>%   filter(PARAM=="ABUNDNR")

ices_data.MUC <- data.frame()

##define radius to extract only monitoring events which where close to sample location
MUC_data <- MUC_data %>%
  mutate(max_long=long *1.05, min_long = long *0.95, 
         max_latt=latt*1.05, min_latt = latt *0.95)
write.table(MUC_data, file=file.path(metadata_path, "Monitoring_range_of_coordinates.tsv"), sep="\t", row.names = F)

##extract close by samples
for ( i in c(1:length(MUC_data$location))){
  small <- ices_data %>% 
    
    filter(Longitude < MUC_data$max_long[i]  , Longitude > MUC_data$min_long[i]) %>%
    filter(Latitude < MUC_data$max_latt[i] , Latitude > MUC_data$min_latt[i] )%>%
    mutate(MUC = rep(MUC_data$location[i], n()))
  
  ices_data.MUC <- rbind(ices_data.MUC, small)

}

#save(ices_data, file=file.path(output_path,"ICES_data__ABUNDNR.RData"))
#save(ices_data.MUC, file=file.path(output_path,"ICES_data__ABUNDNR_EMB262MUCs.RData"))


rm(small,ices_data)

## update species with worms

ices_data <- ices_data.MUC
#extract prokaryotes
ices_data <- ices_data[!grepl("Unicell|MONADER",ices_data$SPECI ),]

rm(ices_data.MUC)


################## check ices taxa without worms id and add to ices_dataset
ices_unknown <- ices_data %>% filter(AphiaID =="")
ices_unknown_species <-ices_data %>% filter(AphiaID =="") %>% pull(SPECI_name) %>% unique()
new_aphias_ICES <- wm_records_taxamatch(ices_unknown_species, marine_only = F)


for ( i in c(1:length(ices_unknown_species))){
  if(length(new_aphias_ICES[[i]]) >0){
    
    ices_rows <- grep(paste(ices_unknown_species[i],"$", sep=""), ices_data$SPECI_name)
    ices_data[ices_rows,"SPECI_name"] <- new_aphias_ICES[[i]]$scientificname
    ices_data[ices_rows,"AphiaID"] <- new_aphias_ICES[[i]]$AphiaID
  }
}

## check worms id
worms_names <- sort(as.numeric(unique(ices_data$AphiaID))) #SPECI_name
worms.tax <- data.frame()

i <- 1
while( i < length(worms_names) ){
  short_ids <- list()
  
  if ( i < length(worms_names) -50){
    short_ids <- wm_record(id = worms_names[c(i:(i+49))], marine_only =F)
  }
  else{
    short_ids <- wm_record(id = worms_names[c(i:length(worms_names))], marine_only =F)
  }
  worms.tax <- rbind(worms.tax,short_ids )
  
  i <- i +49
}
worms.tax <- unique(worms.tax)
rm(short_ids,worms_names)


################## reduce to phytoplankton
worms.tax.all <- worms.tax
#save(worms.tax.all, file=file.path(output_path,"Worms_taxonomy_data_ICES_EMB262MUCs_all.RData"))


#get full worms taxonomy
worms_taxonomy_total <- data.frame()
for (singleID in worms.tax$AphiaID){
  single_entry <- wm_classification(singleID)
  
  single_entry <- single_entry%>% mutate(rank = gsub("Phylum \\(Division\\)", "Phylum", rank),
                                         rank = gsub("Subphylum \\(Subdivision\\)", "Subphylum", rank)) %>%
    select(-AphiaID) %>%
    pivot_wider(names_from = rank, values_from = scientificname) %>% 
    mutate(AphiaID = singleID)
  
  worms_taxonomy_total <- plyr::rbind.fill(worms_taxonomy_total, single_entry)
}
#save(worms_taxonomy_total, file=file.path(output_path,"Worms_FULL_taxonomy_data_ICES_EMB262MUCs_all.RData")) 

colnames(worms_taxonomy_total)
worms_taxonomy_total <- worms_taxonomy_total %>% select(Kingdom,Phylum, Class, Order, Family, Genus, Species, AphiaID )
worms_taxonomy_total <- unique(worms_taxonomy_total)


## reduce dataset to phytoplankton
worms.tax.na <- worms.tax%>%filter(kingdom=="Chromista" & is.na(phylum))

colnames(worms.tax)
worms.tax <- worms.tax %>% filter(phylum=="Myzozoa"| # dinoflagellate
                                    phylum == "Bacillariophyta" |
                                    phylum == "Chlorophyta"|
                                    phylum == "Cryptophyta"|
                                    phylum == "Haptophyta"|
                                    phylum == "Ochrophyta"|
                                    phylum == "Cyanobacteria"|
                                    phylum == "Charophyta") 

#save(worms.tax, file=file.path(output_path,"Worms_taxonomy_data_ICES_EMB262MUCs.RData"))


################## combine ices data and worms taxonomy
#load(file=file.path(output_path,"Worms_taxonomy_data_ICES_EMB262MUCs.RData"))

### units problem
unique(ices_data$MUNIT)
unique(ices_data$final_value_unit)


##transform units
ices_data.m3 <-ices_data %>%ungroup() %>% filter(grepl('/m3', final_value_unit))%>%
      mutate(final_value= as.numeric(final_value) /1000, 
             final_value_unit= gsub("m3", "dm3", final_value_unit))


ices_data.dm3 <- ices_data %>% filter(! grepl("/m3", final_value_unit)) %>%
      mutate(final_value = as.numeric(final_value))
ices_data <- rbind(ices_data.m3, ices_data.dm3)
rm(ices_data.m3, ices_data.dm3)

# sum size classes
ices_data <- merge(ices_data, worms.tax, by="AphiaID") %>%
  mutate(Month = as.numeric(Month), final_value= as.numeric(final_value))%>%
  group_by( PARAM, MUC, tblSampleID,HELCOM_subbasin,HELCOM_L4,Latitude, Longitude ,
            DATE, Year,Month ,Day, 
            phylum, class,order,family,genus,
            SPECI_name,WoRMS_name, AphiaID_accepted,WoRMS_accepted_name,AphiaID)%>%
  summarise(final_value = sum(final_value))

save(ices_data, file=file.path(output_path,"ICES_phytodata_with_worms_EMB262MUCs.RData"))


################# Visualize monitoring frequncy

ices_data_loc <- unique(ices_data[,c(1:10)])
ices_data_loc <- ices_data_loc %>% group_by(MUC,Year,Month)%>% 
  summarise( event_no = n()) %>% mutate( even_pa = ifelse(event_no> 0, 1,0))


plot<- ices_data_loc %>% filter( Month< 7 & Month >3)%>%
  ggplot( aes(x=Year, y= event_no, fill= as.factor(Month))) + 
  geom_bar(position="stack", stat="identity")+facet_wrap( ~MUC, ncol = 1)+
  scale_fill_viridis_d()+
  labs(title= "Abundance of monitoring events in April to June", x= "Years CE", y= "No. of events", fill= "Month in numbers")+
  figure_theme  
print(plot)
ggsave(plot, file=file.path(output_path_extra, "Suppl_ICES__Month_of_Monitoring_events__april-june__5perc.png"),height=8, width=14)

plot<- ices_data_loc %>% filter( Month< 7 & Month >3)%>%
  ggplot( aes(x=Year, y= even_pa, fill= as.factor(Month))) + 
  geom_bar(position="stack", stat="identity")+facet_wrap( ~MUC, ncol = 1)+
  scale_fill_viridis_d()+
  labs(title= "Presence of monitoring events in April to June", x= "Years CE", y= "No. of months", fill= "Month in numbers")+
  figure_theme

print(plot)
ggsave(plot, file=file.path(output_path_extra, "Suppl_ICES__Month_of_Monitoring_events__april-june__5perc.png"),height=8, width=14)

plot<- ices_data_loc %>% filter( Month< 8 & Month >2)%>%
  ggplot( aes(x=Year, y= even_pa, fill= as.factor(Month))) + 
  geom_bar(position="stack", stat="identity")+facet_wrap( ~MUC, ncol = 1)+
  scale_fill_viridis_d()+
  labs(title= "Presence of monitoring events in March to July", x= "Years CE", y= "No. of months", fill= "Month in numbers")+
  figure_theme

print(plot)
ggsave(plot, file=file.path(output_path_extra, "Suppl_ICES__Month_of_Monitoring_events__march-july__5perc.png"),height=8, width=14)


################## extract spring summer month

muc_ICES_data <- ices_data %>%
  group_by(MUC, Year) %>%
  filter(any(Month == 4) & any(Month == 5) & any(Month == 6))%>%
  filter(Month> 3 & Month< 7) %>%
  mutate( phylum = gsub("Myzozoa", "Dinoflagellata", phylum))%>%
  ungroup()%>%
  mutate(across(where(is.numeric), ~ replace_na(., 0)))


save(muc_ICES_data, file=file.path(output_path,"ICES_phytodata_with_worms_EMB262MUCs__SprSum.RData"))
