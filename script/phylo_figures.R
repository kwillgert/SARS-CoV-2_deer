#visualise dated phylogenetic tree

#colour palette
myPal<-colorRampPalette(c("#fc8d59", "#91bfdb"))

#visualise ptree
phylo_tree<- function(){
  par(family="serif", cex=0.95, mar=c(4.5, 0.5, 0.2, 2.1), #set font for plot and label text size, adjust margins
      cex.lab=3, cex.axis=3,#size of axis labels
      mgp=c(3, 2.5, 0)) #move axis tick labels away from axis by increasing mgp
  plot(ptree, show.tip=FALSE) 
  tiplabels(#text="   ", 
    pch=c(21),
    bg=adegenet::fac2col(as.factor(phylo_meta$species), col.pal=myPal), adj=0, cex=2)

  legend(x=0, y=21, #location in figure - needs to be increased if using figure within cowplot
         c("Deer","Human"), #legend labels
         pch=21, #symbol used
         col="black", #outer colour
         pt.bg=c("#fc8d59", "#91bfdb"),
         pt.cex=2.5, #size of symbol
         cex=3, #text size
         ncol=1,
         y.intersp=1.2, #space between symbols on y-axis
         x.intersp=1.5 #space between symbol and text on x-axis
  )
}

p1<-cowplot::ggdraw(phylo_tree) +
  cowplot::draw_plot_label( #add figure label to left top corner 
    label="A",
    hjust=-2, vjust = 1.1,
    size=32,
    fontface = "bold",
    family = "serif")

png("./figures/ptree.png", width=482*2.5, height=381*4)
p1
dev.off()
