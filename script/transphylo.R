#TransPhylo

#Specify generation time and sampling time to use
#Use shape and scale of a gamma distribution
#shape = mean^2/sd^2, scale = sd^2/mean

#4 scenarios considered

#Scenario 1: mean = 5.2 days, standard deviation = 1.72
w.shape=(5.2^2)/(1.72^2)
w.scale=(1.72^2)/5.2
ws.shape=(5.2^2)/(1.72^2)
ws.scale=(1.72^2)/5.2

#view distribution
hist(rgamma(100000, shape=w.shape, scale=w.scale), main=NULL)
hist(rgamma(100000, shape=ws.shape, scale=ws.scale), main=NULL)

#Scenario 2:
mult=2 #explore multiples of default mean and standard deviation generation time, specified in mult

w.shape=((mult*5.2)^2)/((1.72)^2)
w.scale=((1.72)^2)/(mult*5.2)
ws.shape=((mult*5.2)^2)/((1.72)^2)
ws.scale=((1.72)^2)/(mult*5.2)

#Scneario 3:
mult=3

w.shape=((mult*5.2)^2)/((1.72)^2)
w.scale=((1.72)^2)/(mult*5.2)
ws.shape=((mult*5.2)^2)/((1.72)^2)
ws.scale=((1.72)^2)/(mult*5.2)

#Scenario 4: possible extended sampling time in deer 
w.shape=(5.2^2)/(1.72^2)
w.scale=(1.72^2)/5.2
ws.shape=(20^2)/(20^2) #standard deviation for sampling time extended to allow for variance in deer
ws.scale=(20^2)/20

#Specify mcmc iterations to run
mcmc_itr<-500000

# set the time observation stopped dateT. 
# Here set to a a little after the last sampling date to avoid lots of unsampled cases towards the tips
# by adding the length of a mean generation time
dateT=length_time + 5.2

#Run TransPhylo MCMC

res<-inferTTree(ptree,mcmcIterations=mcmc_itr,w.shape=w.shape,w.scale=w.scale,
                ws.shape=ws.shape,ws.scale=ws.scale, 
                dateT=dateT, 
                thinning=5) #sample every 5 iterations for manageable output

#Plot MCMC traces
plotTraces(res, burnin = 0, extend = F)

#calculate EES for each parameter, should be at least 100
#at default setting, TransPhylo does not estimate the parameter off.p --> ESS 0.
mcmc<-convertToCoda(res, burnin=0.25)
effectiveSize(mcmc)

#remove burn-in from MCMC output
res_red<-burnin(record=res, burnin=0.25)

#summary of res output after removing burn-in
res_red_sum<-data.frame(as.mcmc.resTransPhylo(res_red, burnin = 0)) %>% 
  summarise_at(vars(pi:off.p), list(mean,min, max)) %>% 
  pivot_longer(cols=everything(), names_to = c("variable", ".value"), names_pattern = "([^_]+)_(.*)") %>% 
  rename(mean="fn1", min="fn2", max="fn3")

#Plot MCMC traces after removing burn-in
plotTraces(res_red, extend = F)

#calculate EES for each parameter after removing burn-in
mcmc<-convertToCoda(res_red, burnin=0)
effectiveSize(mcmc)

