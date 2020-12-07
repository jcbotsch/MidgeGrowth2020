#====Load Packages====
library(tidyverse)
library(lubridate)
library(nlme)
library(broom)

#====Prepare Data for GPP====
#this code is adapted from McCormick et al. (in review) "Intra- and interannual shifts in the partitioning of benthic and pelagic primary production in Lake Mývatn, Iceland"

#====Routine Light Profiles====
#profiles
profiles <- read_csv("raw data/profiles.csv") %>% 
  filter(sta == 33, month(sampledate) %in% c(5:8))

#examine data
yr_list2 = c(2013:2020)
for (i in yr_list2){
  plot =
    ggplot(subset(profiles, year(profiles$sampledate)==i & !is.na(light)), aes(x=sampledepth, y=light)) +
    facet_wrap(~sampledate)+
    geom_point()+scale_y_continuous(trans='log')+
    geom_smooth(method='lm', se=F)+
    ggtitle(paste(i))
  print(plot)
}

#fit a lm to the log light for each date and pull out the slope 
light.atten <- profiles %>% 
  filter(!is.na(light)) %>% 
  nest(-sampledate) %>% 
  mutate(slope = map(data, ~coef(lm(log(light+0.1)~sampledepth, data = .))[2]),
         int = map(data, ~coef(lm(log(light+0.1)~sampledepth, data = .))[1])) %>% 
  unnest(c(slope, int)) %>% 
  mutate(kD = abs(slope)) %>% 
  select(sampledate, kD)

#find any NAs
profiles %>% 
  filter(sampledate %in% {light.atten %>% 
           filter(is.na(kD))}$sampledate) %>% 
  select(sampledate, light, comments)


#look at patterns through time and compare to AM inland waters fig 3
light.atten %>% 
  filter(!is.na(kD)) %>% 
  ggplot(aes(x= sampledate, y  = kD))+
  geom_line()+
  geom_point()+
  facet_wrap(~year(sampledate), scales = "free_x")

# for one instance in 2014, the incubation (and midges and bchl) was done a day after the light profile and pelagic chl; the only reason for alignment is the later analysis of benthic pmax w/ bchl and preceding light. I decided to adjust the date for the profile (and pchl later so that they have the same date)
light.atten = light.atten %>% mutate(sampledate = replace(sampledate, sampledate=="2014-06-18", "2014-06-19")) 


#====Routine Incubations====

#read data
incub_raw <- read_csv("raw data/NEP.csv") 


#prepare data
incub <- incub_raw %>% 
  mutate(sampledate = as.Date(sampledate, format = "%m/%d/%Y")) %>%
  filter(!is.na(sampledate), 
         sta ==33,
         month(sampledate) %in% c(5:8)) %>% 
  mutate(inc_start = as.POSIXct(strptime(paste(sampledate,inc_time_start,sep="-"), 
                                         "%Y-%m-%d-%H:%M", tz="GMT")),
         inc_end = as.POSIXct(strptime(paste(sampledate,inc_time_end,sep="-"), 
                                       "%Y-%m-%d-%H:%M", tz="GMT"))) %>% 
  select(sampledate, sampletime, c(light0_init:dark_light), inc_start, inc_end, do_init, do_final, column_depth, wtemp_final)

#found 2 typos for DO mg/L and 1 for wtemp 
incub <-  incub %>% 
  mutate(do_final=replace(do_final, (sampledate=="2016-07-05" & coreid==3), 7.74), #though the datasheet reads 77.4
         do_init=replace(do_init, (sampledate=="2017-07-17" & coreid==5), 10.89),
         wtemp_final=replace(wtemp_final, (sampledate=="2017-08-20" & coreid==10), 10))


# summary stats on incubation duration
incub %>%  filter(dark_light != '-1') %>% 
  mutate(inc_length = as.numeric(inc_end-inc_start)) %>% 
  summarize(min = min(inc_length, na.rm = T), 
            max = max(inc_length, na.rm = T),
            mean= mean(inc_length, na.rm = T),
            sd = sd(inc_length, na.rm = T)) 

inc2 = incub %>% 
  filter(dark_light != '-1') %>% #remove pitcher readings
  mutate(delta_do_rate = (do_final - do_init)/as.numeric(inc_end-inc_start), #mg O2 L-1 h-1
         do_gm2hr = delta_do_rate * column_depth/100) %>% #converting units to g O2 m-2 hr-1 (1 mg/L = 1g/m3)
  group_by(sampledate) %>%
  mutate(wtemp = mean(wtemp_final, na.rm=T)) %>%
  ungroup() %>%
  select(sampledate, light0_init, light50_init, light0_final, light50_final, light_logger, coreid, dark_light, delta_do_rate, do_gm2hr, wtemp) %>%
  filter(!is.na(delta_do_rate))

inc2 %>% group_by(sampledate, dark_light) %>% tally() %>% print(n=Inf)

inc2 %>% filter(dark_light==0 & do_gm2hr>0) %>% print(n=Inf) #pos dark readings
inc2 %>% filter(dark_light==1 & do_gm2hr<0) %>% print(n=Inf) #neg light readings
inc2 %>% filter(sampledate == "2014-08-09")


inc2wide = spread(inc2, dark_light, do_gm2hr) %>%
  rename(RESP=`0`, NEP=`1`)

# units of RESP and NEP are g O2 m-2 hr-1
inc_agg = inc2wide %>% 
  group_by(sampledate) %>%
  summarize(RESP.mean=mean(RESP, na.rm=T)*-1,
            RESP.var = var(RESP, na.rm=T),
            RESP.n = sum(!is.na(RESP)),
            NEP.mean=mean(NEP, na.rm=T),
            NEP.var = var(NEP, na.rm=T),
            NEP.n = sum(!is.na(NEP)),
            GPP.mean = RESP.mean + NEP.mean, # Use the routine St33 incubations to calculate the mean GPP and 
            GPP.sd = sqrt(RESP.var + NEP.var), #use variances of light and dark reps to calculate sd of GPP
            GPP.se = sqrt(  NEP.var/NEP.n + RESP.var/RESP.n   ), #GPP standard error based on propagated NEP and ER s.e.
            light50_init=mean(light50_init, na.rm=T),
            light50_final=mean(light50_final, na.rm=T),
            light50=(light50_init+light50_final)/2,
            light_logger=mean(light_logger, na.rm=T),
            wtemp=mean(wtemp, na.rm=T)) %>% #wtemp is the mean final temp reading from tubes
  ungroup()



# Note: I put a lot of thought regarding the best way to characterize the light during the incubation. I was hesitant to base light availability on the hobo data because this doesn't account for attentuation of light through the water column, and I couldn't come up with a way to determine a kD during the routine incubations. So I went with using the direct measuremnts of PAR at 0.5 m. The obvious disadvantage is that there are only two snapshots of light availability, with the initial and final readings. But I was more comfortable that it was an in situ measurement

#3 dates were missing a light50 reading (either initial or final) aside from 2019
inc_agg %>% filter(is.na(light50)) %>% as.data.frame()

#use a single light reading (either initial or final) for these 3 dates
inc_agg = inc_agg %>% 
  mutate(light50=replace(light50, (sampledate=="2016-06-05"), 490.00),
         light50=replace(light50, (sampledate=="2017-07-24"), 950.00))


#2019 the light meter was broken. I'll use the light loggers to estimate.
#how reliably can you estimate light using the kd and logger data?

#equation to estimate light at 0.5 m, given kD and light logger at depth = 0
estlight <- function(light_logger, kD){
  y = log(light_logger/54+0.1)-kD*0.5 #divide by 54 to convert lux to par
  y2 = exp(y)-0.1
  return(y2)
}


lightcomp1 <- inc_agg %>% 
  left_join(light.atten) %>% 
  mutate(est50 = estlight(light_logger = light_logger, kD = kD)) %>%  #estimate par from logger and then calculate value at 0.5 depth
  lm(light50~0+est50, data = .) 

summary(lightcomp1)
#r2 = 0.69. Not great, but since that's what we have we'll go ahead and 



inc_agg %>% 
  left_join(light.atten) %>% 
  mutate(est50 = estlight(light_logger = light_logger, kD = kD)) %>%  #estimate par from logger and then calculate value at 0.5 depth
  ggplot(aes(x = est50, y = light50))+
  geom_point()+
  geom_abline(slope=coef(lightcomp1))+
  coord_equal()


#for values where we have estimates for kd, estimate from light logger
inc_agg <- inc_agg %>% 
  left_join(light.atten %>% 
              mutate()) %>% 
  mutate(light50 = ifelse(is.na(light50)&!is.na(kD), estlight(light_logger=light_logger, kD = kD), light50)) 

inc_agg %>% filter(year(sampledate)==2019) %>% 
  select(sampledate, light_logger, light50)

#find still missing
inc_agg %>% filter(is.na(light50)) %>% select(sampledate)

#interpolate kD for period of broken light meter (6/10/2019)
date1 <- as.Date("2019-05-31")
date2 <- as.Date("2019-06-16")
  
kdinterp1 <- light.atten %>% 
  filter(sampledate %in% c(date1, date2)) %>% 
  lm(kD~sampledate, data = .)

#find kD estimate
estkds1 <- inc_agg %>% 
  filter(sampledate >date1, sampledate<date2) %>% 
  select(sampledate) %>% 
  mutate(estkD = predict(kdinterp1, newdata = .))

#interpolate kD for period of broken light meter (6/16-8/1/2019)
date3 <- as.Date("2019-08-01")

kdinterp2 <- light.atten %>% 
  filter(sampledate %in% c(date2, date3)) %>% 
  lm(kD~sampledate, data = .)

estkds2 <- inc_agg %>% 
  filter(sampledate >date2, sampledate<date3) %>% 
  select(sampledate) %>% 
  mutate(estkD = predict(kdinterp2, newdata = .))

estkds <- bind_rows(estkds1, estkds2) 

#estimate light at 50cm for the last few
inc_agg <- inc_agg %>% 
  left_join(estkds) %>% 
  # filter(!is.na(estkD)) #check
  mutate(light50 = ifelse(!is.na(estkD), estlight(light_logger = light_logger, kD = estkD), light50)) 

inc_agg %>% filter(is.na(light50))


inc_agg %>% 
  filter(year(sampledate)==2019) %>% select(sampledate, light_logger, kD, estkD, light50)


#check sampling dates
iaunmatch<- unique(inc_agg$sampledate)[which(!unique(inc_agg$sampledate) %in% unique(light.atten$sampledate))]
inc_agg %>% filter(sampledate %in% iaunmatch)
light.atten %>% filter(sampledate %in% iaunmatch) #misssing because light meter broken


launmatch <- unique(light.atten$sampledate)[which(!unique(light.atten$sampledate) %in% unique(inc_agg$sampledate))]
light.atten %>% filter(sampledate %in% launmatch)

inc_agg %>% filter(sampledate %in% launmatch) #days where incubations were not conducted


# write_csv(inc_agg, "clean data/incubations.csv")

#====Analysis=====
#adapted from McCormick et al. in Review

# The routine incubations should capture benthic GPP max, but probably the ones done on cloudy/low light days may not reach max GPP.
# 1) Fit MM P-I curve with the obs GPP
# 2) using the K parameter and B parameter, rearrange and solve for the 'corrected' Pmax by plugging in the OBSERVED gpp and light, and B parameter from the curve fit to incubation data 

# import the aggregated dataset
d <- inc_agg %>% 
  rename(GPP.obs = GPP.mean) %>%
  arrange(sampledate)


# fit MM P-I curve based on 0.5m par
mm.ben.inc = nls(GPP.obs ~ B*light50/(K+light50) , 
                 data=d %>% filter(!is.na(GPP.obs), !is.na(light50), GPP.obs>0), #non-sensical to have neg GPP; one incubation in 2018 was done when the bloom was really thick -- not enough light to capture prod
                 start=list(B=0.3, K=250))
coef(mm.ben.inc)

K.ben = as.numeric(coef(mm.ben.inc)[2])

# Benthic incubation figure (FIG 1): gpp.obs and light50 (The mean recorded PAR at 0.5 m during incubation)
d %>%  filter(GPP.obs>0) %>% #non-sensical to have neg GPP; one incubation in 2018 was done when the bloom was really thick -- not enough light to capture prod
  ggplot(aes(x=light50, y=GPP.obs)) +
  # scale_shape_manual(values=c(22,15,21,16,24,17,23))+
  geom_point(aes(col=as.factor(year(sampledate))), size=1)+
  geom_smooth(method = "nls", aes(x=light50, y=GPP.obs), formula = y ~ B*x/(K+x) ,
              method.args = list(start = list(K=250,B=0.3)),
              data = d %>% filter(!is.na(GPP.obs), !is.na(light50), GPP.obs>0), color="black", se = FALSE, size=1)+
  ylab(bquote(GPP['obs']~' ( g'~O['2']~ m^-2~hr^-1~')')) +
  xlab(bquote('PAR ('~mu*'mol photons'~ m^-2~s^-1~')')) +
  theme(legend.title = element_blank(), legend.position = "top",
        plot.margin = margin(10,15,10,10),
        legend.text = element_text(size=8), legend.spacing.x = unit(0.075, 'cm'),
        axis.title = element_text(size=10), axis.text = element_text(size=8))

#ggsave("Figure1.pdf", plot = last_plot(),  width = 3.5, height=3.75, units='in', device = "pdf", dpi = 300) 
#ggsave("Figure1.tiff", plot = last_plot(),  width = 6.8, height=7.5, units='cm', device = "tiff", dpi = 1200) 




# Platt inhibition model fit to incubations
ben_platt_inhib = nls(GPP.obs ~ B*(1 - exp(-a*light50/B))*exp(-beta*light50/B),
                      data=d %>% filter(GPP.obs>0, !is.na(GPP.obs), !is.na(light50), month(sampledate)<9),
                      start=list(B=0.3, a=0.005, beta=-9.75e-06))
coef(ben_platt_inhib) # neg beta suggests no photoinhibition



####################################
# Now, we can calc the corrected benthic GPPmax, using B=GPPobs*(K+L)/L
# and propagate the error; using the error for GPPobs that is based on NEP and RESP replicate core errors 
# after talking with Tony on 27 June, have taken out the variance of K from the 'global' model from the error propagation (instead use it for prop estimated in situ GPP later on)

K.ben.var = vcov(mm.ben.inc)[2,2]


# to get partial derivs 
deriv( ~ GPP*(K+light)/light, 
       c("GPP", "K", "light"), 
       function(GPP, K, light) {}) #can remove light from this row b/c don't need partial

#partial derivs are:
# df/dGPPobs = (K+L)/L
# df/dK = GPPobs/L

#add cols to dataset to include the corrected ben_GPPmax and its error
data = d %>% filter(!is.na(GPP.obs), GPP.obs>0) %>%
  mutate(GPP.var = GPP.sd^2) %>% 
  rename(GPPobs.var = GPP.var, GPPobs.sd = GPP.sd, GPPobs.se = GPP.se) %>% #just to clarify these are based on observed GPP
  mutate(
    #This uses the K parameter from the incubation P-I curve to correct GPP
    ben_GPPmax =  GPP.obs*(K.ben+light50)/light50, #'correct' Pmax based on recorded incubation light
    
    ben_GPPmax.sd = sqrt( (  (K.ben+light50)/light50)^2*GPPobs.var  ), #error from K was removed and I included it instead in propagation of est GPP error
    ben_GPPmax.var = ben_GPPmax.sd^2)

data %>% 
  ggplot(aes(y = ben_GPPmax, x = year(sampledate), col = month(sampledate)))+
  geom_jitter(width = 0.15)

data %>% 
  ggplot(aes(y = ben_GPPmax*1000*0.375, x = sampledate))+
  geom_line()+
  facet_wrap(~year(sampledate), scales = "free_x")
#these match figure4.tiff

rhobo_raw <- read_csv("raw data/hobo_clean_20Nov20.csv", col_types = "ccTddcccddclcDtDtTT")


weather = read_csv("raw data/weather_2010_2019.csv") %>% select(Date_Time, RADGL) %>%
  rename(datetime=Date_Time) %>% filter(year(datetime)>=2012) %>%
  filter(month(datetime)>=5, month(datetime)<=8)%>%
  mutate(sampledate=date(datetime),
         hour = hour(datetime),
         hour = as.numeric(hour),
         surf_par_weath = 1.76225*RADGL) %>% #from surfPAR_regressions.R
  select(-c(datetime, RADGL)) 


#hourly irradiance data
incid_irr <- rhobo_raw %>% 
  filter(site == "st33", layer == "out") %>% 
  mutate(sampledate=date(datetime),
         hour = hour(datetime),
         surf_par_hobo = 0.0084519*lux, #from surfPAR_regressions.R
         hour=as.numeric(hour)) %>%
  filter(month(datetime)>=5, month(datetime)<=8)%>%
  select(sampledate, datetime, hour, lux, surf_par_hobo) %>% 
  group_by(sampledate, hour) %>%
  summarize(surf_par_hobo= mean(surf_par_hobo, na.rm=T)) %>%
  ungroup()

incid_irr = full_join(incid_irr, weather) %>% arrange(sampledate) %>%
  mutate(surf_par = ifelse(is.na(surf_par_hobo), surf_par_weath, surf_par_hobo))

expand = expand.grid(sampledate = unique(data$sampledate), hour=c(0:23)) %>% 
  arrange(sampledate) %>% 
  mutate(sampledate = as.Date(sampledate, format="%Y-%m-%d"))



#### 7 day running ave of surface light levels (3 day before and 3 after sampledate)
orig = expand %>% mutate(orig_date = factor(sampledate),
                         set = "original")
full_back1 = orig %>% mutate(sampledate = sampledate - 1, set = "back1")
full_back2 = orig %>% mutate(sampledate = sampledate - 2, set = "back2")
full_back3 = orig %>% mutate(sampledate = sampledate - 3, set = "back3")
full_up1 = orig %>% mutate(sampledate = sampledate + 1, set = "up1")
full_up2 = orig %>% mutate(sampledate = sampledate + 2, set = "up2")
full_up3 = orig %>% mutate(sampledate = sampledate + 3, set = "up3")

full_week = rbind(orig, full_back1, full_back2, full_back3, full_up1, full_up2, full_up3) %>% arrange(sampledate)

full_week = left_join(full_week, incid_irr)

head(full_week)


full_week %>% filter(is.na(surf_par)) # not missing any surface light data
full_week %>% group_by(orig_date) %>% tally() %>% filter(n!=168) #7 full days = 168 entries

#avergae diel light for the 1 week periods
full_week_agg = full_week %>% 
  group_by(orig_date, hour) %>% 
  summarize(surf_par = mean(surf_par)) %>% ungroup() %>%
  mutate(sampledate = as.Date(orig_date, format="%Y-%m-%d")) %>% select(-orig_date)


full3 = left_join(data, full_week_agg)
head(full3)

full3 %>% filter(is.na(surf_par)) # not missing any surface light data
full3 %>% group_by(sampledate) %>% tally() %>% filter(n!=24) #full 24 hr for all days

# 
full.clean = full3 
head(full.clean)


# to get partial derivs in R
deriv( ~ P*light/(K+light), 
       c("P", "K", "light"), 
       function(P, K, light) {}) #can remove light from this row b/c don't need partial

#partial derivs are:
# df/dGPPobs = (K+L)/L
# df/dK = -P*L/((K+L)^2)

# add the ben_GPPest for each hour and its error; ben par: Iz = I0*e^(-k*z)
full.clean = full.clean %>%
  mutate(kL = ifelse(is.na(kD), estkD, kD),
         ben_par = surf_par*exp(-1*kL*3.3),
         ben_GPPest = ben_GPPmax*ben_par/(K.ben+ben_par), #eq gpp=Pmax*L/(K+L); L=e^(-kD*z)
         
         ben_GPPest.sd = sqrt(  (ben_par/(ben_par+K.ben))^2*ben_GPPmax.sd^2 #from Pmax   (df/dPmax)^2 * var of Pmax) 
                                + (-ben_GPPmax*ben_par/((K.ben+ben_par)^2))^2*K.ben.var )) %>% # from K (df/dK)^2 * var for K.ben
  select(-ben_par)


#### aggregate to day
head(full.clean)


daily = full.clean %>% 
  group_by(sampledate) %>%
  summarize(ben_GPPmax = mean(ben_GPPmax), #these are constant for a day
            ben_GPPmax.se = mean(GPPobs.se), #this is just for the error bars on figure 4
            ben_GPPmax.sd = mean(ben_GPPmax.sd),
            
            
            ben_GPPest_day = sum(ben_GPPest), #sum hourly rates for day
            ben_GPPest_day.sd =  sqrt(sum(ben_GPPest.sd^2))) #error is sqrt of the sum of variances
 

daily %>% filter(is.na(ben_GPPest_day))

daily %>% 
  mutate(ben_GPPest_day = ben_GPPest_day*0.375,
         ben_GPPest_day.sd = ben_GPPest_day.sd*0.375) %>% 
  ggplot(aes(x = sampledate, y = ben_GPPest_day))+
  facet_wrap(~year(sampledate), scales = "free_x")+
  geom_point(alpha = 0.5)+
  geom_errorbar(aes(ymin = ben_GPPest_day-ben_GPPest_day.sd, ymax = ben_GPPest_day+ben_GPPest_day.sd), width = 0)+
  geom_line()

daily %>% 
  count(year(sampledate))

unique(daily$sampledate)[which(!unique(daily$sampledate) %in% unique(rmp$sampledate))]

unique(rmp$sampledate)[which(!unique(rmp$sampledate) %in% unique(daily$sampledate))] %>% sort()



rmp%>% 
  filter(sta == "33", year(sampledate)>2012) %>% 
  full_join(daily %>% 
              mutate(ben_GPPest_day = ben_GPPest_day*0.375)) %>% #convert to C
  select(sampledate, ben_GPPest_day, Pdc) %>% 
  gather(taxon, productivity, -sampledate) %>% 
  mutate(taxon = ifelse(taxon=="Pdc", "tanytarsini", "algae")) %>% 
  ggplot(aes(x = sampledate, y = productivity, col = taxon))+
  facet_wrap(~year(sampledate), scales = "free_x")+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_line()+
  geom_point(size = 2)+
  labs(x = "Date",
       y = "Productivity in g C m-2 d-1")

scaleFUN <- function(x) sprintf("%.1f", x)

rmp%>% 
  filter(sta == "33", year(sampledate)>2012) %>% 
  full_join(daily %>% 
              mutate(ben_GPPest_day = ben_GPPest_day*0.375,
                     ben_GPPest_day.sd = ben_GPPest_day.sd*0.375)) %>% 
  select(sampledate, Pdc.se, ben_GPPest_day, ben_GPPest_day.sd, Pdc) %>% 
  ggplot(aes(x = sampledate))+
  facet_wrap(~year(sampledate), scales = "free_x")+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_line(aes(y = Pdc), col = "firebrick4", data = . %>% filter(!is.na(Pdc)))+
  geom_point(aes(y = Pdc, col = "Tanytarsini"), size = 2, data = . %>% filter(!is.na(Pdc)))+
  # geom_errorbar(aes(ymin = Pdc-Pdc.se*2, ymax = Pdc+Pdc.se*2, col = "Tanytarsini"), width = 0, data = . %>% filter(!is.na(Pdc)))+
  geom_line(aes(y = ben_GPPest_day), col = "dodgerblue", data = . %>% filter(!is.na(ben_GPPest_day)))+
  geom_point(aes(y = ben_GPPest_day, col = "Primary Producers"), size= 2, data = . %>% filter(!is.na(ben_GPPest_day)))+
  # geom_errorbar(aes(ymin = ben_GPPest_day-2*ben_GPPest_day.sd, ymax = ben_GPPest_day+2*ben_GPPest_day.sd, col = "Primary Producers"), size = 0.7, width = 0, data = . %>% filter(!is.na(ben_GPPest_day)))+
  labs(x = "Date",
       y = expression("Productivity g C"~m^{-2}~d^{-1}),
       color = "")+
  scale_color_manual(values = c("firebrick4", "dodgerblue"),
                     breaks = c("Tanytarsini", "Primary Producers"))+
  # scale_y_continuous(breaks=c(0.0, 1.0, 2.0), labels=scaleFUN)+
  scale_x_date(date_breaks = "1 month", date_labels = "%b")
  


rmp%>% 
  filter(sta == "33", year(sampledate)>2012) %>% 
  full_join(daily %>% 
              mutate(sampledate = if_else(sampledate == as.Date("2019-07-03"), as.Date("2019-07-04"), sampledate),
                     sampledate = if_else(sampledate == as.Date("2016-08-02"), as.Date("2016-08-01"), sampledate))) %>% 
  ggplot(aes(x = ben_GPPest_day*0.375, y = Pdc, col = (yday(sampledate))))+
  geom_errorbarh(aes(xmin = ben_GPPest_day*0.375-ben_GPPest_day.sd*0.375, xmax = ben_GPPest_day*0.375+ben_GPPest_day.sd*0.375))+
  geom_errorbar(aes(ymin = Pdc-Pdc.se, ymax = Pdc+Pdc.se), width = 0)+
  geom_point()+
  labs(x = expression("Primary Productivity g C"~m^{-2}~d^{-1}),
       y = expression("Tanytarini Productivity g C"~m^{-2}~d^{-1}),
       color = "")+
  geom_abline(slope = 1)+
  scale_color_viridis_c()
