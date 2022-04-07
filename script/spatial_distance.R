#estimate spatial distance and distance in transmission events between each pair of deer

#remove human samples from matrix with number of intermediates between each sampled pair
deer_intermed<-mat_df %>% 
  rownames_to_column(var="spp_ID") %>%  #create column for row names
  select(-contains("human")) %>%  #remove human cases from columns
  filter(str_detect(spp_ID, "deer")) #%>% #only keep deer samples in rows

#set coordinate reference system
r<-raster()
crs(r)<-"+proj=longlat + datum=WGS84"

#create a crs object
longlat<-st_crs(r)
class(longlat)

#create dataset for where and when deer samples collected
deer_location<-meta_deer %>% 
  dplyr::select(Sample_ID, County, Date_Collected, PCR,  Longitude, Latitude) %>% 
  mutate(Sample_ID=as.character(Sample_ID)) %>% 
  filter(PCR==1) %>%  #only include positive cases
  inner_join(y=phylo_meta %>% select(Sample_ID, spp_ID), by="Sample_ID") #only include subset of samples used in transmission analysis

##reorder data to match order of transmission events data
deer_location <- deer_location[order(match(deer_location$spp_ID, deer_intermed$spp_ID)),]
# test the labels match:
identical(deer_intermed$spp_ID, as.vector(deer_location$spp_ID))

#convert to sf object
deer_location<-st_as_sf(deer_location, coords=c("Longitude", "Latitude"), crs=longlat)
distance<-st_distance(deer_location)
distance<- distance/1000 #convert to km

deer_intermed<-deer_intermed %>% 
  column_to_rownames(var="spp_ID") #move column with rownames back into rownames

#plot y axis -> the value of (i,j) from the matrix of number of transmission events
#x axis -> the value for the same pair (i,j) from distance matrix

png("./figures/distance.png", width=414*3, height=480*3) 
par(family="serif", cex=4) #size and font of text
plot(as.vector(distance), as.vector(as.matrix(deer_intermed)), 
     pch = 21, #shape used
     col="lightgrey", bg="#FC8D59", 
     cex=1, #size of points
     ylab="Intermediate transmission events", xlab="Spatial distance (km)")
dev.off()

