#intermediates between sampled cases

# Matrix of how many intermediates were inferred to be in the transmission chain between each pair of individuals
mat=computeMatTDist(res_red, burnin=0)

#rename sample names in matrix for figure
rownames(mat)<-phylo_meta$spp_ID
colnames(mat)<-phylo_meta$spp_ID

#order names by county and date of sampling
mat_order<-phylo_meta %>% 
  arrange(species, County, day_collected) %>% 
  select(spp_ID) 

mat<-mat[mat_order$spp_ID,mat_order$spp_ID, drop=FALSE]

#use traveling salesman problem approach to order sequences according to least intermediate transmission events
library(TSP)

#convert matrix to class tsp
mat_TSP<-as.TSP(mat)
print(mat2)
labels(mat2)

#solve TSP using nearest neighbor method
TSP_sol<-solve_TSP(mat_TSP, method="nn")

#extract order, starting with deer-Allamakee-136-1
TSP_order<-data.frame(orig_loc=cut_tour(TSP_sol, cut="deer-Allamakee-136-1", exclude_cut = FALSE)) %>% 
  rownames_to_column(var = "spp_ID") %>% 
  left_join(phylo_meta %>% select(spp_ID, species)) %>% 
  arrange(species)   #sort by species

#order intermediates matrix accordingly
mat<-mat[TSP_order$spp_ID,TSP_order$spp_ID, drop=FALSE]

#set colours
myCols <- colorRampPalette(c("#01665e", "#35978f","#80cdc1","#c7eae5","#f5f5f5","#f6e8c3","#dfc27d","#bf812d"))

png("./figures/intermediates.png", width=482*4, height=381*4)
lattice::levelplot(mat,xlab='',ylab='', 
                   scales=list(x=list(rot=90, cex=0.7), y=list(cex=0.7)), #cex = edit text size on x and y axis, rot=rotates text 90 degrees
                   col.regions=myCols,
                   colorkey=list(labels=list(cex=4))) #change size of legend text

dev.off()

#convert matrix of number of intermediates to dataframe
mat_df<-data.frame(mat) %>% 
  mutate_if(is.numeric, round, digits=2) %>% 
  na_if(0) #set middle diagonal to NA to not include in summary stats

intermed_summary<-data.frame(
median_intermed=median(as.matrix(mat_df), na.rm=TRUE),
mean_intermed=mean(as.matrix(mat_df), na.rm=TRUE),
min_intermed=min(as.matrix(mat_df), na.rm=TRUE),
max_intermed=max(as.matrix(mat_df), na.rm=TRUE))
