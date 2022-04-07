#tidy tree for TransPhylo and visualisaion

#check if tree binary
is.binary(cov_tree)

#remove any multifurcations if present:
cov_tree<-multi2di(cov_tree)

#check for negative branch lengths - this should have been avoided by using common ancestor node heights in BEAST
sum(cov_tree$edge.length<0) #any negative branches?
sum(cov_tree$edge.length>0) 
sum(cov_tree$edge.length==0) #any branches of length zero?

#extract time period over which samples collected so that first sample is set as day 0 (if working in days)
length_time<- as.numeric(max(phylo_meta$Date_Collected, na.rm=TRUE)-min(phylo_meta$Date_Collected, na.rm=TRUE))

#convert phylogeny to a ptree object and plot it, aligning the time frame
ptree<- ptreeFromPhylo(cov_tree, dateLastSample=length_time)
