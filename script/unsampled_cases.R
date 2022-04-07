#unsampled cases

#create empty list where can store iterations
sim_case_find<-list()
itr<-0 #record iteration in itr

for(i.itr in seq_along(1:length(res_red))) {
  
  itr<-itr+1 #update iteration
  
  ctree=res_red[[i.itr]]$ctree #select simulation
  
  #for each predicted tree: extract transmission data
  
  # extract the selected transmission tree 
  ttree=extractTTree(ctree)
  
  #summary data frame of selected tree
  info_tree=data.frame(cbind(tree_ID=1:nrow(ttree$ttree), ttree$ttree))  %>%   #extract column 1: time when infected, column 2: time when sampled, column 3: infector
    rename("t_infected"="V2", "t_sampled"="V3", "infector"="V4") %>% 
    mutate(t_infected=round(t_infected), t_sampled=round(t_sampled))

  numCases=nrow(info_tree) #total number cases
  numSamp=length(ttree$nam)  # number sampled
  numUnsamp=numCases-numSamp  # number cases unsampled 
  
cases<-data.frame(numSamp, numUnsamp, numCases) %>%  
 mutate(case_finding=numSamp/numCases) %>% 
  mutate(itr=itr, .before=1)
  
  #add case finding estimates from iteration to list
  sim_case_find[[i.itr]]<- as.data.frame(cases)
  
}


#bind iterations into a dataframe
case_find_sum<-rbindlist(sim_case_find) %>% 
  summarise(mean_unsamp=round(mean(numUnsamp), digits=0), median_unsamp=median(numUnsamp), min_unsamp=min(numUnsamp), max_unsamp=max(numUnsamp), 
            mean_case_find=round(mean(case_finding), digits=4), median_case_find=round(median(case_finding), digits=4), min_case_find=round(min(case_finding), digits=4), max_case_find=round(max(case_finding), digits=4))

#sampled and unsampled cases over time by week
#requires reading R script functions.R before running
unsamp_cases=getIncidentCases_weekly(record=res_red, burnin=0, min_date=min(phylo_meta$Date_Collected, na.rm=TRUE))

unsamp_cases

ggsave("./figures/unsampled_cases.png", width=13, height=10)  

#combine plot of samples collected over time by species and unsampled cases into one figure

cowplot::plot_grid(
  epicurve_samples,
  NULL, #include empty plot to increase space between plots
  unsamp_cases+
    theme(plot.margin = margin(5.5, 35.5, 5.5, 11.5))+ #increase plot margin on left and right side to compensate for removed space from y-axis title
    ylab(NULL), #remove y-axis label
  rel_widths = c(1, 0.06, 1), #specify width of each plot
  labels = c("A", "", "B"),
  label_x=-0.03, #move label to left to prevent overlapping with y-axis labels
  label_y=1.015, #move label up
  label_size=22,
  label_fontfamily="serif",
  nrow=1)+
  theme(plot.margin = unit(c(0.3,0.3,0.3,0.5), "cm")) #add additional margin

ggsave("./figures/sampled_cases.png", width=20, height=10)
