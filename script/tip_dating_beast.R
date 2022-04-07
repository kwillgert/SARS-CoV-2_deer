#date randomisation test
#based on tutorial in https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2F1755-0998.12603&file=men12603-sup-0001-SupInfo.pdf

#load package
library(TipDatingBeast)

?RandomDates

#produce 20 repetitions of all sequences with randomly assigned dates
RandomDates(name="./data/reps", reps = 20, writeTrees = T)
#generated xml files are then run in BEAST

#plot comparison of parameter estimates of original datasets with randomised dates, no burnin required as already removed when combining logs in LogCombiner
PlotDRT(name="./data/reps/output/DRTtree.comb", reps=20, burnin=0.0)
#when prompted to specify which parameter of the log files you want DRT to be performed on, type VAR and hit enter
#then type clock.rate

#read summary file of date randomisation test
drt<-read.csv("./data/reps/output/clock.rate.stats.csv") %>% 
  mutate(Dataset=ifelse(drt$calibr==0, "real", "randomized"))

DRT<-ggplot(data=drt %>% 
         mutate(Dataset=factor(Dataset, levels=c("real", "randomized"))),
        aes(x=calibr, y=mean)) +
  geom_point(aes(colour=Dataset))+ #change colour of real dataset to red, leaving date-randomized dataset black
  geom_errorbar(aes(ymin=lowerHPD, ymax=HigherHPD, colour=Dataset))+
  scale_colour_manual(values = c("red", "black"))+
  geom_hline(yintercept=c(drt[drt$calibr==0,"lowerHPD"], drt[drt$calibr==0,"HigherHPD"]), 
             linetype="dashed", col='red')+ #add dotted line
  theme_classic() + 
  theme(text=element_text(size=16,  family="serif"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+ 
  labs(x =NULL,    
       y="Clock rate\n") +
  scale_y_continuous(limits=c(0, NA), expand = expansion(mult = c(0.01, 0.1)))+
  scale_x_continuous(expand = c(0.001,0.01))

DRT

ggsave("./figures/DRT.png", width=6.5, height=5.35)  


#combine figure with root-to-tip distance

comb_temp_signal<-cowplot::plot_grid(
  temp_signal_spp+
    scale_y_continuous(limits=c(0, max_rot_tip_dist), expand = expansion(mult = c(0.01, 0.1))),
    #theme(legend.position = "none"), #remove legend from plot
  NULL, #include empty plot to increase space between plots, size specified in rel_widths
  DRT,
  rel_widths = c(1, 0.04, 1),
  labels = c("A", "", "B"),
  label_fontfamily="serif",
  nrow=1,
  align= "h") #align plots horizontally

comb_temp_signal

ggsave("./figures/temp_signal_comb.png", width=12, height=5.35)

