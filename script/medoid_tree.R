#find the medoid transmission tree
med=medTTree(res_red, burnin=0) #no burnin removed as already removed from res_red

med

plot(med)

#rename names of sequences
med[["nam"]]<-phylo_meta$spp_ID

#set colours
myPal<-c("#8C510A", "#BF812D", "#DFC27D", "#80CDC1", "#35978F", "#01665E")

#function to plot the transmission tree, allowing to save it in cowplot
transm_tree<- function(){
  par(family="serif", mar=c(4.5, 0.5, 0.2, 6.1), 
      cex.lab=3, cex.axis=3,
      mgp=c(3, 2.5, 0), 
      yaxs="i") #yaxs="i" --> avoid axes extending beyond origin and max value
  plotCTree_edit(med, showLabels=TRUE, showStars=TRUE, col=myPal, cex=0.95) #size of branch text specified in cex
  }

p2<-cowplot::ggdraw(transm_tree) +
  cowplot::draw_plot_label( #add figure label
    label="B",
    hjust=-2.5, vjust = 1.1,
    size=32,
    fontface = "bold",
    family = "serif")

png("./figures/transmission_tree.png", width=482*2.5, height=381*4)
p2
dev.off()

#summary data frame of medoid tree
info_tree=data.frame(cbind(tree_ID=1:nrow(ttree$ttree), ttree$ttree))  %>%   #extract column 1: time when infected, column 2: time when sampled, column 3: infector of medoid tree
  rename("t_infected"="V2", "t_sampled"="V3", "infector"="V4") %>% 
  mutate(t_infected=round(t_infected), t_sampled=round(t_sampled))

numCases=nrow(info_tree) #total number cases
numSamp=length(ttree$nam)  # number sampled
numUnsamp=numCases-numSamp; # number cases unsampled 

info_tree<-info_tree %>% 
  add_column(Sample_ID=c(ttree$nam, rep(NA, numCases-numSamp)), .after="tree_ID") %>%  #add sample ID, NA for unsampled cases
  left_join(phylo_meta, by=c("Sample_ID"="spp_ID")) %>%  #add meta data
  mutate(sampled=ifelse(!is.na(Sample_ID), 1, 0)) #indicate if node sampled (1=sampled, 0=unsampled)

#number sampled deer and human
info_tree %>% 
  group_by(species) %>% 
  summarise(sampled=n())

#save network data
write.xlsx(info_tree, "./data/deer_CoV_network.xlsx", sheetName="network_meta")

# Plot the matrix of the probability of direct transmission between all pairs of individuals in the outbreak
mat1=computeMatWIW(res_red, burnin=0)

png("./figures/transmission_probability.png", width=482*3, height=381*3)
lattice::levelplot(mat1,xlab='',ylab='', 
                   scales=list(x=list(rot=90)),#scales to rotate x-axis 90 degrees
                   colorkey=list(labels=list(cex=3))) #change size of legend text 

dev.off()

mat1[mat1>0] #check if any direct transmission pairs
mat1[mat1>0.5] #transmission pairs with probability of >0.5

