#epidemic curve of samples used in analysis

#by week
meta_sum <- phylo_meta %>% 
  mutate(week=isoweek(Date_Collected)) %>% 
  mutate(week_start=floor_date(Date_Collected, unit="weeks", week_start=1)) %>%   #fist day of week
  group_by(week_start, week, species) %>% 
  summarise(cases=n()) %>% 
  ungroup()

epicurve_samples<-ggplot(data=meta_sum %>% mutate(species=factor(species, levels=c("human", "deer"))), 
       aes(x=week_start))+
  geom_col(mapping=aes(y=cases, fill=species))+ 
  theme_classic() + 
  theme(text=element_text(size=28,  family="serif"), 
        axis.text.x=element_text(angle=60, hjust=1, size=26),
        axis.text.y=element_text(size=26),
        legend.position = c(1,1), 
        legend.justification = c(1, 1), 
        legend.direction="vertical",
        legend.spacing.y = unit(-0.05, 'cm'),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))+
  labs(x ="\nDate of collection",    
       y="Weekly cases\n",
       fill=NULL) +
  scale_fill_manual(values = c("#91bfdb","#fc8d59"))+
  scale_y_continuous(limits=c(0, NA), expand = expansion(mult = c(0.01, 0.1)))+
  scale_x_date(limits=c(as.Date("2019-12-23"), as.Date("2021-02-22")+7), 
               expand = c(0.0,0.0), 
               date_breaks = "4 week", date_labels = "%d-%b")  # set x-axis to 4 weekly 

epicurve_samples

ggsave("./figures/epicurve_samples.png", width=13, height=10)
