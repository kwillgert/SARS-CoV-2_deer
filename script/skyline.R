
#Inferred skyline population plots for human and deer SARS-CoV-2 and human SARS-CoV-2 alone

#read human and deer skyline data generated in BEAST
skyline_HD<-read.csv("./data/skyline_human_deer", 
                     header = TRUE, sep = "\t", skip=1) %>% 
  select(-date, -datetime, -milliseconds) %>% 
  mutate(date_pop=as.Date(date_min)+time, .before=mean) %>% #convert to date
  drop_na(median)

#upper limit of 95% interval
upper_pop<-max(skyline_HD$upper, na.rm=TRUE)

#Human only
phylo_meta_human<-phylo_meta %>% 
  filter(species=="human")

#date of last collected sample
date_max<-max(phylo_meta_human$Date_Collected, na.rm=TRUE)
#date of first collected sample
date_min<-min(phylo_meta_human$Date_Collected, na.rm=TRUE)

#extract time period over which samples collected so that first sample can be set as day 0
length_time_human<- as.numeric(date_max)-as.numeric(date_min)

#read human only skyline data generated in BEAST
skyline_H<-read.csv("./data/skyline_human", 
                    header = TRUE, sep = "\t", skip=1) %>% 
  select(-date, -datetime, -milliseconds) %>% 
  mutate(date_pop=as.Date(date_min)+time, .before=mean) #convert to date

date_min_skyH<-min(skyline_H$date_pop)
date_max_skyH<-max(skyline_H$date_pop)

skyline_plot_H<-ggplot(data = skyline_H, aes(x = date_pop)) +
  geom_ribbon(aes(ymin = lower,  # lower edge of ribbon
                  ymax = upper), # upper edge of ribbon
              alpha = 0.5,   # make semi-transparent
              fill = "#dfc27d", # fill colour
              color = NA) +     # no border color
  geom_line(aes(y = median), alpha=0.75) +       # line for median
  theme_bw() +  
  theme(text=element_text(size=28,  family="serif"), 
        axis.text.x=element_text(angle=60, hjust=1, size=26),
        axis.text.y=element_text(size=26))+
  xlab("Time") +      
  ylab("Population size \n")+ 
  scale_y_continuous(limits=c(0, round(upper_pop)), expand = expansion(mult = c(0.01, 0.1)))+ #set upper limit to same as for other plot
  scale_x_date(limits=c(floor_date(date_min_skyH, unit="weeks", week_start=1), #limits of x-axis
                        floor_date(date_max_skyH, unit="weeks", week_start=1)+7),
    breaks=seq.Date(from=floor_date(min(skyline_H$date_pop), unit="weeks", week_start=1), #breaks of x-axis
                               to=floor_date(max(skyline_H$date_pop), unit="weeks", week_start=1), 
                               by="4 week"), #specify breaks to display on x-axis
               expand = c(0.01,0.01), #expand space to left and right of x-axis by 1% 
               date_labels = "%d-%b")  # set x-axist to weekly 

#human and deer

date_min<-min(phylo_meta$Date_Collected, na.rm=TRUE) #date of first collected sample
date_max<-max(phylo_meta$Date_Collected, na.rm=TRUE) #date of last collected sample

#first date in skyline
date_min_skyHD<-min(skyline_HD$date_pop)
date_max_skyHD<-max(skyline_HD$date_pop)

skyline_plot_HD<-ggplot(data = skyline_HD, aes(x = date_pop)) +
  geom_ribbon(aes(ymin = lower,  # lower edge of ribbon
                  ymax = upper), # upper edge of ribbon
              alpha = 0.5,   # make semi-transparent
              fill = "#dfc27d", # fill
              color = NA) +     # no border color
  geom_line(aes(y = median), alpha=0.75) +       # line for median
  theme_bw() +  
  theme(text=element_text(size=28,  family="serif"), 
        axis.text.x=element_text(angle=60, hjust=1, size=26),
        axis.text.y=element_text(size=26))+
  xlab("Time") +      
  ylab("Population size \n")+ 
  scale_y_continuous(limits=c(0, round(upper_pop)), expand = expansion(mult = c(0.01, 0.1)))+
  scale_x_date(limits=c(floor_date(date_min_skyH, unit="weeks", week_start=1), #set limits of x-axis to same as human only data
                        floor_date(date_max_skyH, unit="weeks", week_start=1)+7),
               breaks=seq.Date(from=floor_date(min(skyline_H$date_pop), unit="weeks", week_start=1), #breaks of x-axis set to match those of human data only for comparison
                               to=floor_date(max(skyline_H$date_pop), unit="weeks", week_start=1),
                               by="4 week"), #specify breaks to display on x-axis
               expand = c(0.01,0.01), #expand space to left and right of x-axis by 1% 
               date_labels = "%d-%b") # set x-axist to weekly 


cowplot::plot_grid(
  skyline_plot_HD,
  NULL, #include empty plot to increase space between plots, size specified in rel_widths
  skyline_plot_H+
    theme(plot.margin = margin(5.5, 35.5, 5.5, 11.5))+ #increase plot margin on left and right side to compensate for removed space from y-axis title
    ylab(NULL), #remove y-axis label
  rel_widths = c(1, 0.06, 1), #specify width of each plot
  labels = c("A", "", "B"),
  label_x=-0.03, #move label to left to prevent overlapping with y-axis labels
  label_y=1.015, #move label up
  label_size=22,
  #labels="AUTO",
  label_fontfamily="serif",
  nrow=1,
  align="v")+ #same width (aligned vertically)
  #Add some space around the edges  
  theme(plot.margin = unit(c(0.3,0.3,0.3,0.5), "cm")) #add additional margin

ggsave("./figures/skyline_pop.png", width=20, height=10)

