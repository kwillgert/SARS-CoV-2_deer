#Code for transmission analysis of SARS-CoV-2 in humans and deer

#folders:
#data: folder for all data files
#script: folder for all R script files
#figures: folder where figures are saved


#load packages
source("./script/packages.R")

#read meta data
source("./script/data.R") 

#load functions
source("./script/functions.R")

#set text font to New Times Roman for session
par(family="serif")

#root-to-tip distance to assess for temporal signal
source("./script/root_tip_dist.R")

#read dated phylogenetic tree
source("./script/read_tree.R")

#date randomisation test to assess for temporal signal
source("./script/tip_dating_beast.R")

#skyline population figures
source("./script/skyline.R")

#epidemic curve of samples used in analysis
source("./script/epicurve_samples.R")

#prepare phylogenetic tree for assessment in TransPhylo and visualisation
source("./script/tidy_tree.R")

#visualise phylogenetic tree
source("./script/phylo_figures.R")

#perform transmission analysis in TransPhylo
source("./script/transphylo.R")

#extract medoid tree
source("./script/medoid_tree.R")

#illustrate transmission tree
source("./script/transmission_network.R")

#intermediates between sampled cases
source("./script/intermediates.R")

#estimate generation time between infector and sampled individuals by species
source("./script/generation_times.R")

#estimate time between sampled cases becoming infected and being sampled by species
source("./script/sampling_times.R")

#estimate unsampled cases and case finding
source("./script/unsampled_cases.R")

#calculate geographical distance and transmission events between sampled deer
source("./script/spatial_distance.R")
source("./script/Maps.Rmd")



