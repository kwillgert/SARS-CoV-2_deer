#read maximum likelihood tree
cov_tree_ml<-ape::read.tree(file="./data/ml_phylo.treefile") #ML phylogenetic tree of SARS-CoV-2 data used in analysis, made in iqtree - not dated
cov_tree_ml$tip.label<- str_remove(cov_tree_ml$tip.label, pattern="\\|.*$") #If tip labels contain sampling date, remove sampling date from tip label for consistency with meta data

#meta data of maximum likelihood phylogenetic tree
phylo_meta_ml<-data.frame(Sample_ID=cov_tree_ml$tip.label) %>% 
  left_join(combined_meta) 

#root-to-tip distance to assess for temporal signal
#need to source script named "functions.R" before running
temp_signal_spp<-plotRootToTip_bySpp(datedTree=cov_tree_ml, meta=phylo_meta_ml, drop_human="no")

temp_signal_spp

ggsave("./figures/root_to_tip_distance_spp.png", width=6.5, height=5.35)

