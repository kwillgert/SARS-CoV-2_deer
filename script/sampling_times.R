#estimate sampling times for humans and deer

#create empty list where can store iterations
sim_sampl_time<-list()
itr<-0 #record iteration in itr

for(i.itr in seq_along(1:length(res_red))) {
  
  itr<-itr+1 #update iteration
  
  ctree=res_red[[i.itr]]$ctree #select simulation
  
  #for each predicted transmission tree: extract transmission data
  
  # extract the selected transmission tree  
  ttree=extractTTree(ctree)
  
  #summary data frame of selected tree
  info_tree=data.frame(cbind(tree_ID=1:nrow(ttree$ttree), ttree$ttree))  %>%   #extract column 1: time when infected, column 2: time when sampled, column 3: infector
    rename("t_infected"="V2", "t_sampled"="V3", "infector"="V4") %>% 
    mutate(t_infected=round(t_infected), t_sampled=round(t_sampled))
  
  numCases=nrow(info_tree) #total number cases
  numSamp=length(ttree$nam)  # number sampled
  numUnsamp=numCases-numSamp; # number cases unsampled 
  
  info_tree<-info_tree %>% 
    add_column(Sample_ID=c(ttree$nam, rep(NA, numCases-numSamp)), .after="tree_ID") %>%  #add sample ID, NA for unsampled cases
    left_join(phylo_meta) %>%  #add meta data
    mutate(sampled=ifelse(!is.na(Sample_ID), 1, 0)) #indicate if node sampled (1=sampled, 0=unsampled)
  
  #sampled individuals
  sampl_time<-info_tree %>% 
    filter (sampled==1) %>% 
    mutate(sampl_time=t_sampled-t_infected) %>% 
    mutate(itr=itr, .before=1)
  
  #add sampling time estimates from iteration to list
  sim_sampl_time[[i.itr]]<- as.data.frame(list(sampl_time))
  
}

#bind iterations into a dataframe
sampl_time_df<-rbindlist(sim_sampl_time) 

#summary of median and range of sampling time by species
sampl_time_summary<-sampl_time_df %>% 
  group_by(itr, species) %>% 
  summarise(median_sampl_time=median(sampl_time), min_sampl_time=min(sampl_time), max_sampl_time=max(sampl_time)) %>%  #median, min and max sampling time for an iteration
  ungroup() %>% 
  group_by(species) %>% 
  summarise(median_median_sampl_time=median(median_sampl_time), min_median_sampl_time=min(median_sampl_time), max_median_sampl_time=max(median_sampl_time), #median, min and max of the median sampling time across all iterations
            median_min_sampl_time=median(min_sampl_time), min_min_sampl_time=min(min_sampl_time), max_min_sampl_time=max(min_sampl_time), #median, min and max of the lowest sampling time for all iterations 
            median_max_sampl_time=median(max_sampl_time), min_max_sampl_time=min(max_sampl_time), max_max_sampl_time=max(max_sampl_time)) ##median, min and max of the max sampling time across all iterations

view(sampl_time_summary)

sampl_time_distr <- sampl_time_df %>% 
  group_by(itr, species, sampl_time) %>% 
  summarise(n_sampl_time=n()) %>%  #number of individuals with indicated sampling time for a given iteration and species
  ungroup() %>% 
  group_by(itr, species) %>% 
  mutate(prop_sampl_time=round(n_sampl_time/sum(n_sampl_time), digits=2)) %>%  #proportion of individuals with a given sampling time
  ungroup()


sampl_time_sum <- sampl_time_distr %>% 
  group_by(species, sampl_time) %>% 
  summarise(mean_n=round(mean(n_sampl_time), digits=2), median_n=median(n_sampl_time), min_n=min(n_sampl_time), max_n=max(n_sampl_time), 
            mean_prop=round(mean(prop_sampl_time), digits=2), median_prop=median(prop_sampl_time), min_prop=min(prop_sampl_time), max_prop=max(prop_sampl_time))

#plot figure of sampling time by species
sampl_time<-ggplot(data = sampl_time_sum, 
       aes(x = sampl_time)) +
  geom_ribbon(aes(ymin = min_prop,  # lower edge of ribbon
                  ymax = max_prop, # upper edge of ribbon
                  fill=species), 
              alpha = 0.3,   # make semi-transparent
              color = NA) +     # no border color
  geom_line(aes(y = median_prop, color=species), size=1.2, alpha=1) +  # line for median
  scale_colour_manual(values = c("#fc8d59", "#91bfdb"))+
  scale_fill_manual(values = c("#fc8d59", "#91bfdb"))+
  theme_bw() +  
  theme(text=element_text(size=24,  family="serif"))+
  labs(x="\nSampling time (days)",      
       y="Proportion hosts\n", 
       colour="Species",
       fill="Species")+
  scale_y_continuous(limits=c(0, NA), expand = expansion(mult = c(0.01, 0.1)))+
  scale_x_continuous(limits=c(0, NA), expand = expansion(mult = c(0.00, 0.0)))

sampl_time

ggsave("./figures/sampl_times_prop.png", width=14, height=8)

#combine generation time and sampling time into one figure

# extract the legend from one of the plots to assign single species legend
legend_samplT <- cowplot::get_legend(
  sampl_time + 
    theme(legend.box.margin = margin(0, 0, 0, 0.5)) # create some space to the left of the legend
)

comb_genT_samplT<-cowplot::plot_grid(
  gen_time+
    theme(legend.position = "none"),
  NULL, #include empty plot to increase space between plots
  sampl_time+
    theme(legend.position = "none",
          plot.margin = margin(5.5, 10.5, 5.5, 36.5))+ #remove legend from plot, increase plot margin on left and right side to compensate for removed space from y-axis title
    ylab(NULL), #remove y-axis label
  rel_widths = c(1, 0.04, 1), #specify width of each plot
  labels = c("A", "", "B"),
  label_size=22,
  label_fontfamily="serif",
  nrow=1)

# add the legend to the plot. Give it 0.4 of the width of one plot via rel_widths.
cowplot::plot_grid(comb_genT_samplT, legend_samplT, rel_widths = c(3, 0.4))

ggsave("./figures/genT_sampT_comb.png", width=20, height=8)
