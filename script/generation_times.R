#estimate generation time between infector and sampled individuals by species

#create empty list where can store iterations
sim_gen_time<-list()
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
  numUnsamp=numCases-numSamp; # number cases unsampled 
  
  info_tree<-info_tree %>% 
    add_column(Sample_ID=c(ttree$nam, rep(NA, numCases-numSamp)), .after="tree_ID") %>%  #add sample ID, NA for unsampled cases
    left_join(phylo_meta) %>%  #add meta data
    mutate(sampled=ifelse(!is.na(Sample_ID), 1, 0)) #indicate if node sampled (1=sampled, 0=unsampled)
  
  #sampled individuals
  sampled_hosts<-info_tree %>% 
    filter (sampled==1)
  
  #extract infector of sampled individuals and calculate generation time between infector becoming infected and infecting sampled host
  gen_time<-info_tree %>% 
    right_join(sampled_hosts %>% select(infector), by=c("tree_ID"="infector")) %>% 
    select(ID_source=tree_ID, t_infected) %>% #ID_source = ID of infector
    left_join(sampled_hosts %>% 
                select(ID_inf=tree_ID, t_infect_inf=t_infected, infector, species_inf=species), #ID_inf = ID of sampled hosts, t_infect_inf= time infector infected sampled host, species_inf=species of sampled host
              by=c("ID_source"="infector")) %>% 
    mutate(gen_time=t_infect_inf-t_infected) %>%  #generation time= time from infector becoming infected to infecting sampled host
    mutate(itr=itr, .before=1)
  
  #add generation time estimates from iteration to list
  sim_gen_time[[i.itr]]<- as.data.frame(list(gen_time))
  
}

#bind iterations into a dataframe
gen_time_df<-rbindlist(sim_gen_time) 

#summary of generation time by species across iterations
gen_time_summary<-
  gen_time_df %>% 
  group_by(itr, species_inf) %>% 
  summarise(median_gen_time=median(gen_time), min_gen_time=min(gen_time), max_gen_time=max(gen_time)) %>%  #median, min and max generation time for an iteration
  ungroup() %>% 
  group_by(species_inf) %>% 
  summarise(median_median_gen_time=median(median_gen_time), min_median_gen_time=min(median_gen_time), max_median_gen_time=max(median_gen_time), #median, min and max of the median generation time across all iterations
            median_min_gen_time=median(min_gen_time), min_min_gen_time=min(min_gen_time), max_min_gen_time=max(min_gen_time), #median, min and max of the lowest generation time for all iterations
            median_max_gen_time=median(max_gen_time), min_max_gen_time=min(max_gen_time), max_max_gen_time=max(max_gen_time)) ##median, min and max of the max generation time across all iterations

view(gen_time_summary)  
  
gen_time_distr<-gen_time_df %>% 
  group_by(itr, species_inf, gen_time) %>% 
  summarise(n_gen_time=n()) %>%  #number of individuals with indicated generation time for a given iteration and species
  ungroup() %>% 
  group_by(itr, species_inf) %>% 
  mutate(prop_gen_time=round(n_gen_time/sum(n_gen_time), digits=2)) %>%  #proportion of individuals with a given generation time
  ungroup()

gen_time_sum <- gen_time_distr %>% 
  group_by(species_inf, gen_time) %>% 
  summarise(mean_n=round(mean(n_gen_time), digits=2), median_n=median(n_gen_time), min_n=min(n_gen_time), max_n=max(n_gen_time), 
            mean_prop=round(mean(prop_gen_time), digits=2), median_prop=median(prop_gen_time), min_prop=min(prop_gen_time), max_prop=max(prop_gen_time))

#plot generation time by species
gen_time<- ggplot(data = gen_time_sum, 
       aes(x = gen_time)) +
  geom_ribbon(aes(ymin = min_prop,  # lower edge of ribbon
                  ymax = max_prop, # upper edge of ribbon
                  fill=species_inf), 
              alpha = 0.3,   # make semi-transparent
              color = NA) +     # no border color
  geom_line(aes(y = median_prop, color=species_inf), size=1.2, alpha=1) +       # line for median
  scale_colour_manual(values = c("#fc8d59", "#91bfdb"))+
  scale_fill_manual(values = c("#fc8d59", "#91bfdb"))+
    theme_bw() +  
      theme(text=element_text(size=24,  family="serif"))+
      labs(x="\nGeneration time (days)",      
           y="Proportion hosts\n", 
           colour="Species",
           fill="Species")+
      scale_y_continuous(limits=c(0, NA), expand = expansion(mult = c(0.01, 0.1)))+
      scale_x_continuous(limits=c(0, NA), expand = expansion(mult = c(0.00, 0.0)))
    
gen_time    

ggsave("./figures/gen_times_prop.png", width=14, height=8)