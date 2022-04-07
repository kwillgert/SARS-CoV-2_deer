#read dated phylogenetic tree
cov_tree<-read.nexus(file="./data/dated_phylo.tre")

#If tip labels contain sampling date, remove sampling date from tip label for consistency with meta data
cov_tree$tip.label<- str_remove(cov_tree$tip.label, pattern="\\|.*$") 
  
#create meta file for samples present in tree
phylo_meta<-data.frame(Sample_ID=cov_tree$tip.label) %>% 
  left_join(combined_meta) 
  
##reorder the meta data to match the order of phylogeny
phylo_meta <- phylo_meta[order(match(phylo_meta$Sample_ID,cov_tree$tip.label)),]

# test the labels match:
identical(cov_tree$tip.label, as.vector(phylo_meta$Sample_ID)) 

#add column for sample names to use in figure
phylo_meta <- phylo_meta %>% 
  mutate(County=replace_na(County, "unknown")) %>% 
  mutate(day_collected=as.numeric(Date_Collected-min(Date_Collected, na.rm=TRUE))) %>%  #using difference in time between collection date of first sample and collection date of the sample (days)
  mutate(spp_ID=str_c(species, County, day_collected, sep="-")) %>%   #add collection date to sequence name
  mutate(spp_ID=paste(spp_ID, data.table::rowid(spp_ID, prefix=""),  sep="-")) #add number if multiple samples from same day for species

