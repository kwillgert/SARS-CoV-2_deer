#illustration of transmission network in igraph

#convert network data to igraph opject

#extract edges
transmission_links<-data.frame(from=info_tree$infector, to=info_tree$tree_ID, sampled=info_tree$sampled) %>% 
  mutate(from=replace(from, from==0, max(to)+1)) %>%  #igraph does not allow ID of 0, change to next available ID number
  mutate(arrow_size=ifelse(sampled==1, 0.1, 0)) %>% #set arrow attributes based on if sampled or not
  mutate(arrow_width=ifelse(sampled==1, 0.2, 0)) %>% 
  mutate(weights=ifelse(sampled==1, 0.4, 0.2)) #used in Kamada Kawai layout, placing vertices connected with low weighted edge closer to eachother
  
#create node list
nodes<-info_tree %>%
  select(ID=tree_ID, species, sampled) %>% 
  add_row(ID=max(info_tree$tree_ID)+1, species=NA, sampled=0) %>% #add ID of node previously called 0
  replace_na(list(species="unknown")) 

nodes<-nodes %>% 
  mutate(color=ifelse(nodes$species=="deer", "#fc8d59", ifelse(nodes$species=="human", "#91bfdb", "#fee090"))) %>%  #set node colour to be based on species where known, otherwise unknown
  mutate(node_size=ifelse(nodes$sampled==1, 3, 0.8)) #%>%  #set node size to be based on if sampled or unsampled 

#creat Igraph object  
net<-graph_from_data_frame(d=transmission_links,
                           vertices=nodes,
                           directed=TRUE)

net <- igraph::simplify(net, remove.multiple = F, remove.loops = T) 

#specify transmission network layout
l<-igraph::layout_with_kk(net, 
                          maxiter=400*vcount(net), #maximum number of iterations to perform
                          weights=E(net)$weights) #edge weights, larger values --> longer edges
l <- norm_coords(l, ymin=-1, ymax=1, xmin=-1, xmax=1) #linearly transform coordinates to fit into given limits

png("./figures/covid19_network_igraph.png", width=482*3, height=381*3)

par(family="serif", cex=1.1, mar=c(3, 0, 2, 3)) #edit font, text size and margins

plot(net, 
     edge.color="darkgrey",
     edge.width=0.7, 
     edge.arrow.mode=0, 
     vertex.label=NA,
     vertex.label.color=V(net)$color, #colour based on species and if sampled or not
     vertex.size=V(net)$node_size, #node size based on if sampled or not
     rescale=FALSE, #prevents coordinates to rescale
     layout=l*1.06 #multiply by value to create more compact or spread out layout
     ) 

legend(x=1, y=-0.8, #location in figure
       c("Deer","Human", "Unknown"), #legend labels
       pch=21, #symbol used
       col="#777777", 
       pt.bg=c(nodes[nodes$species=="deer", "color"][1], nodes[nodes$species=="human", "color"][1], nodes[nodes$species=="unknown", "color"][1]),
      pt.cex=c(nodes[nodes$species=="deer", "node_size"][1]+0.2, nodes[nodes$species=="human", "node_size"][1]+0.2, nodes[nodes$species=="unknown", "node_size"][1]+0.5), 
      cex=2.5, 
      bty="n", 
      ncol=1,
      y.intersp=1.2, #space between legend text vertically
      x.intersp=1.8 #space between legend symbol and text
      )

dev.off()
