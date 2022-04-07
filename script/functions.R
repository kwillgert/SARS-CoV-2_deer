
#######################################################################
# Function to plot root to tip distances against collection date to 
# assess temporal signal in dataset                                 
# Adapted from code written by Andries van Tonder                    
# https://github.com/avantonder/RBCT/blob/main/mainRBCTanalysis.R    
#######################################################################

#Root-to-tip distance by species
plotRootToTip_bySpp <- function(datedTree, meta, drop_human){
  #meta data should include species column with deer and human samples
  #and Date_Collected column with collection date of sample
  
  if(drop_human=="yes"){
  #drop human samples from tree
  datedTree<-drop.tip(phy=datedTree, tip=meta[meta$species=="human","Sample_ID"], trim.internal=TRUE,
                                subtree=FALSE, rooted=is.rooted(cov_tree), collapse.singles=FALSE)
  
  #filter meta data to not include human samples
  meta<-meta %>% 
    filter(species!="human")
  
  }  
  
  sample_dates=meta$Date_Collected #date of samples
  species=meta$species #species samples collected from
  
  # Extract root to tip distances from tree
  RootTotipDistances <- diag(vcv.phylo(datedTree))
  
  # Extract tip labels from tree
  treeLabels <- datedTree$tip.label
 
  #convert date of samples to days since first sample collected
  treeDates <- as.numeric(sample_dates-min(sample_dates))
  
  # Perform linear regression on root to tip distances and tree dates
  treeModel <- lm(RootTotipDistances ~ treeDates)
  
  # Calculate correlation between root to tip distances and tree dates
  treeCorrelation <- cor.test(treeDates, 
                              RootTotipDistances, 
                              method = "pearson", 
                              conf.level = 0.95)
  
  # Summarise linear regression
  modelSummary <- summary(treeModel)
  
  # Add root to tip distances and tree dates to dataframe for plotting
  RootTotipDF <- data.frame(RootTotipDistances,
                            treeDates, species)
  
  # Create root to tip plot and add variables calculated above
  RootTotipPlot <- ggplot(data = RootTotipDF, 
                          aes(treeDates,  
                              RootTotipDistances)) + 
    geom_point(aes(colour=species),
               alpha = 0.70, 
               size = 2) + 
    scale_color_manual(values = brewer.pal(3,"RdYlBu")[c(1,3)])+
    theme_classic() + 
    labs(x = "\nDays since first sample",
         y = "Root to tip distance \n",
         colour="Species" 
         ) +
    geom_smooth(method = 'lm',
                fullrange=T,
                se=T) + 
    ggtitle(paste("First sample: ",
                  min(sample_dates, na.rm=TRUE),
                  "; ",
                  #"TMRCA: ",
                  #round(xIntercept, 0),
                  "Slope: ",
                  formatC(modelSummary$coefficients[2], 
                          format = "e", 
                          digits = 3),
                  ";", 
                  "p-value: ",
                  format(modelSummary$coefficients[8],
                         digits = 2),
                  "; ",
                  "\n", 
                  "Correlation Coefficient: ",
                  round(treeCorrelation$estimate, 3), 
                  ";", 
                  "R^2: ", 
                  format(modelSummary$r.squared,
                         digits = 3)))+
    theme(text=element_text(size=10,  family="serif"), axis.text=element_text(size=16),
          axis.title=element_text(size=16), legend.title=element_text(size=16), legend.text=element_text(size=14))
  return(RootTotipPlot)
}


#specify proportion of MCMC output to be removed as burning
#record=Output from inferTTree
#burning=proportion of MCMC output to be discarded

burnin<-function(record, burnin){
      if (burnin > 0) 
        record = record[(round(length(record) * burnin)+1):length(record)]
}


#sampled and unsampled cases over time by week
#edit of function from TransPhylo package getIncidentCases, 
#min_date=date of first sample
getIncidentCases_weekly<-function (record, burnin, min_date) 
{
  record = record[max(1, round(length(record) * burnin)):length(record)]
  lr = length(record)
  minTime = Inf
  maxTime = -Inf
  for (i in 1:lr) { 
    thisTT <- extractTTree(record[[i]]$ctree)$ttree
    minTime = min(minTime, min(thisTT[, 1]))
    maxTime = max(maxTime, max(thisTT[, 1]))
  }
  
  minDate_week<-floor_date(min_date+minTime, unit="weeks", week_start=1) #date of first day of week of first case
  maxDate_week<-floor_date(min_date+maxTime, unit="weeks", week_start=1) #date of first day of week of last case
  outbreak_weeks<-seq.Date(from=minDate_week, to=maxDate_week, by="week") #create bins for each week up until last sampled case
  
  numBins=length(outbreak_weeks)
  
  breaks = as.numeric(outbreak_weeks-min_date) #specify breaks in days from first sample collected
  sampcounts = matrix(NA, lr, numBins)
  unsampcounts = matrix(NA, lr, numBins)
  for (i in 1:lr) {
    thisTT <- extractTTree(record[[i]]$ctree)$ttree #extract transmission tree for each iteration
    IND = !is.na(thisTT[, 2]) #which cases have been sampled 
    sampcounts[i, ] = tabulate(findInterval(thisTT[IND, 1], #extracts number of sampled cases infected during the time interval specified in bins
                                            breaks), numBins)
    unsampcounts[i, ] = tabulate(findInterval(thisTT[!IND, #extracts number of unsampled cases infected during the time interval specified in bins
                                                     1], breaks), numBins)
  }
  samp = colSums(sampcounts)/lr #average number of cases sampled per each time interval over the iterations performed
  unsamp = colSums(unsampcounts)/lr #average number of cases unsampled per each time interval over the iterations performed
  
  res = list(Time = breaks, sampledCases = samp, unsampCases = unsamp)

    mydata = data.frame(time_day=res$Time, week_start=outbreak_weeks, sampled=res$sampledCases, unsampled=res$unsampCases) %>% 
      pivot_longer(cols=c(sampled, unsampled), names_to = "sampling", values_to = "cases")
    
    unsamp_fig<-ggplot(data=mydata %>% mutate(sampling=factor(sampling, levels=c("unsampled", "sampled"))), 
           aes(x=week_start))+
      geom_col(mapping=aes(y=cases, fill=sampling))+ 
      theme_classic() + 
      theme(text=element_text(size=28,  family="serif"), 
            axis.text.x=element_text(angle=60, hjust=1, size=26),
            axis.text.y=element_text(size=26),
            legend.position = c(1,1), 
            legend.justification = c(1, 1), 
            legend.direction="vertical",
            legend.spacing.y = unit(-0.05, 'cm'),
            legend.background = element_blank(),
            legend.box.background = element_rect(colour = "black")
      )+
      labs(x ="\nDate",    
           y="Weekly cases\n",
           fill=NULL) +
      scale_fill_manual(values = c("#f6e8c3", "#01665e"))+
      #set y axis to start at 0 and remove excessive space below y-axis=0 in expand scale
      scale_y_continuous(limits=c(0, NA), expand = expansion(mult = c(0.01, 0.1)))+
      #remove white space on left and right side of x-axis
      scale_x_date(breaks=seq.Date(from=as.Date(min(outbreak_weeks)), 
                                   to=as.Date(max(outbreak_weeks)), 
                                   by="4 week"), #specify breaks to display on x-axis
                                   expand = c(0.01,0.01), #expand space to left and right of x-axis by 1% 
                                    # set x-axist to weekly
                                    date_labels = "%d-%b")
  
  return(unsamp_fig)
}



#edit of function from TransPhylo package plotCTree
#edits plot of coloured transmisison tree
plotCTree_edit<- function (tree=med, showLabels = TRUE, showStars = TRUE, cols = myPal, 
                           maxTime = NA, cex = 1){
  nam = tree$nam
  tree = tree$ctree
  nsam <- sum(tree[, 2] + tree[, 3] == 0) #sampled
  nh <- nrow(tree) - 3 * nsam + 1 #unsampled
  ntot <- nsam + nh #total number of samples
  oldpar <- par("yaxt", "bty", "xpd")
  on.exit(par(oldpar))
  par(yaxt = "n", bty = "n", xpd = T)
  
  par(family="serif") #set text font and scale size
  plot(0, 0, #produce empty plot
       type = "l", 
       xlim = c(min(tree[, 1]), ifelse(is.na(maxTime), max(tree[, 1]), maxTime)), #min and max x-axis
       ylim = c(-1.2, nsam + 5), #min and max y-axis, add some space above and below for margins
       #xlab = "Days since first sample", #x-axis label
       xlab="",
       ylab = "", #y-axis label
       cex.lab=3, #font size of axis title
       cex.axis=3, #Increase size of axis values
         ) 
  host <- tree[, 4]
  if (ntot > 1) {
    if (is.na(cols[1])) 
      grDevices::palette(grDevices::rainbow(min(1024, ntot)))
    else grDevices::palette(cols)
  }
  root <- which(host == 0)
  ys <- matrix(0, nsam, 1)
  todo <- cbind(root, 0, 0.5, 1)
  while (nrow(todo) > 0) {
    w <- todo[1, 1]
    x <- todo[1, 2]
    y <- todo[1, 3]
    scale <- todo[1, 4]
    if (tree[w, 2] == 0 && tree[w, 3] == 0) {
      ys[w] <- y
    }
    else if (tree[w, 3] == 0) {
      todo <- rbind(todo, cbind(tree[w, 2], tree[w, 1], 
                                y, scale, deparse.level = 0))
    }
    else {
      todo <- rbind(todo, cbind(tree[w, 2], tree[w, 1], 
                                y + scale/2, scale/2, deparse.level = 0), cbind(tree[w, 
                                                                                     3], tree[w, 1], y - scale/2, scale/2, deparse.level = 0))
    }
    todo <- rbind(todo[-1, ])
  }
  ys <- rank(ys)
  for (i in ((nsam + 1):nrow(tree))) {
    children <- c()
    todo <- i
    while (length(todo) > 0) {
      children = c(children, todo[1])
      todo = c(todo[-1], setdiff(tree[todo[1], 2:3], 0))
    }
    ys[i] <- mean(ys[children[which(children <= nsam)]])
  }
  todo <- cbind(root, tree[root, 1])
  while (nrow(todo) > 0) {
    w <- todo[1, 1]
    x <- todo[1, 2]
    y <- ys[w]
    col = host[w]
    if (tree[w, 2] == 0 && tree[w, 3] == 0) {
      lines(c(x, tree[w, 1]), c(y, y), col = col, lwd = 2)
      if (showLabels) 
        text(tree[w, 1], y, nam[w], cex = cex, pos = 4)
    }
    else if (tree[w, 3] == 0) {
      lines(c(x, tree[w, 1]), c(y, y), col = col, lwd = 2)
      todo <- rbind(todo, cbind(tree[w, 2], tree[w, 1]))
    }
    else {
      lines(c(x, tree[w, 1]), c(y, y), col = col, lwd = 2)
      lines(c(tree[w, 1], tree[w, 1]), cbind(ys[tree[w, 
                                                     2]], ys[tree[w, 3]]), col = col, lwd = 2)
      todo <- rbind(todo, cbind(tree[w, 2], tree[w, 1]), 
                    cbind(tree[w, 3], tree[w, 1]))
    }
    todo <- rbind(todo[-1, ])
  }
  todo <- cbind(root, tree[root, 1])
  while (nrow(todo) > 0 && showStars) {
    w <- todo[1, 1]
    x <- todo[1, 2]
    y <- ys[w]
    col = host[w]
    if (tree[w, 2] == 0 && tree[w, 3] == 0) {
    }
    else if (tree[w, 3] == 0) {
      points(tree[w, 1], y, col = "red", pch = 8, cex=0.4) #set shape (pch) and size of shape (cex) highlighting transmissison events
      todo <- rbind(todo, cbind(tree[w, 2], tree[w, 1]))
    }
    else {
      todo <- rbind(todo, cbind(tree[w, 2], tree[w, 1]), 
                    cbind(tree[w, 3], tree[w, 1]))
    }
    todo <- rbind(todo[-1, ])
  }
  return(invisible(tree))
}

