#====Load Packages====
library(tidyverse)
library(lubridate)
library(viridis)
library(gridExtra)
library(lme4)
library(car)

#Define function to give NA if all NAs else give mean
meanna <- function(x){
  if(all(is.na(x))){
    NA
  }
  else{
    mean(x, na.rm = TRUE)
  }
}

#====set aesthetics====
theme_set(theme_bw()+
            theme(panel.grid = element_blank(),
                  strip.background = element_blank(),
                  legend.position = "bottom"))

#====Set Seed====
set.seed(12345)

#====read in files====

meta <- read_csv("clean data/MG_meta.csv") %>% 
  mutate(algae_conc2 = ifelse(algae_conc == 0, 2e-3, algae_conc)) #check
chl <- read_csv("clean data/MG_chl.csv")
om <- read_csv("clean data/MG_om.csv")
ndvi <- read_csv("clean data/MG_ndvi.csv")
ep_raw <- read_csv("clean data/MG_ep.csv")
cm_raw <- read_csv("clean data/MG_cm.csv")
cc_raw <- read_csv("clean data/MG_cc.csv")


#====Figure 1: Chorophyll and organic content in stock sediment====
omplot <- om %>% 
  left_join(meta) %>% 
  ggplot(aes(algae_conc, perc.org))+
  geom_point(size =3)+
  scale_y_continuous(labels = scales::percent_format())+
  labs(x = "Sediment Quality",
       y = "Organic Content")

#corrected for pheophytin
corrchlplot <- chl %>% 
  left_join(meta) %>% 
  ggplot(aes(x = algae_conc, y = chl))+
  geom_point(size = 2, alpha = 0.8)+
  scale_y_continuous(labels = scales::comma_format())+
  labs(x = "",
       y = "Chlorophyll \u03BCg/L") +  #\u03BC is the unicode for mu
  theme(axis.text.x = element_blank())

#uncorrected
uncorrchlplot <- chl %>% 
  left_join(meta) %>% 
  ggplot(aes(x = algae_conc, y = uncorr_chl))+
  geom_point(size = 2, alpha = 0.8)+
  scale_y_continuous(labels = scales::comma_format())+
  labs(x = "",
       y = "Chlorophyll \u03BCg/L") + #\u03BC is the unicode for mu
  theme(axis.text.x = element_blank())


grid.arrange(corrchlplot, omplot, nrow  = 2)
grid.arrange(uncorrchlplot, omplot, nrow  = 2)


#=====Question 1: Effects of Midges on Primary Producers=====
## Midge effect on production
nep <- ep_raw %>% 
  left_join(meta) %>% 
  mutate(coreid = as.character(coreid),
         box = as.character(box))

nep14 <- nep %>% 
  filter(day == 14)

g14log <- lm(gpp~log(algae_conc2)*midge+box, data = nep14) 
summary(g14log)
# plot(g14log)


nep22 <- nep %>% 
  filter(day == 22)

g22log <- lm(gpp~log(algae_conc2)*midge+box, data = nep22)
summary(g22log)
# plot(g22log)


#====Figure 2: GPP====
nep %>% 
  mutate(day = paste("Day", day)) %>% 
  ggplot(aes(x = algae_conc2, y = gpp, fill = midge))+
  facet_wrap(~day)+
  geom_point(size  = 2, alpha = 0.7, shape = 21)+
  geom_smooth(aes(color = midge), method = "lm", size = 0.75, se = FALSE)+
  scale_fill_manual(values = c("black", "gray60"))+
  scale_color_manual(values = c("black", "gray60"))+
  labs(x = "Initial Sediment Quality",
       y = expression(GPP(g~O[2]~m^{-2}~hr^{-1})),
       color = "",
       fill = "")+
  scale_x_log10()


#====Do the slopes differ between the two sample dates====

#parametric bootstrapping verison

#set number of bootstraps
nboot <- 2000

#define function for parametric bootstrapping
para.boot <- function(data1, data2, response.var, lm1, lm2, results, nboot){

    #define data frames for simulated data
    dat1 = data1
    dat2 = data2
    
    #perform bootstrap
    for(i in 1:nboot){
      #simulate results estimates for day 14
      dat1$response.var <- simulate(lm1)[[1]] #simulate data from model
      lm1boot <- update(lm1, response.var ~., data = dat1) #perform lm on simulated data
      
      #extract coefficients 
      results$intercept[i] = coef(summary(lm1boot))[1,"Estimate"]
      results$algae[i] = coef(summary(lm1boot))[2,"Estimate"]
      results$NoMidge[i] = coef(summary(lm1boot))[3,"Estimate"]
      results$box2[i] = coef(summary(lm1boot))[4,"Estimate"]
      results$interaction.algaeNoMidge[i] = coef(summary(lm1boot))[5,"Estimate"]
      
      #add coefficients from algae and its interaction with midges to estimate effect of algae conc on response in absence of midges
      results$algae_x_NoMidge[i] = results$algae[i] + results$interaction.algaeNoMidge[i]
      
      #perform same process for second lm==
      dat2$response.var <- simulate(lm2)[[1]] #simulate data from model
      lm2boot <- update(lm2, response.var~., data = dat2) # fit regression to simulated data
      
      #extract coefficients
      results$intercept[i + nboot] = coef(summary(lm2boot))[1,"Estimate"] 
      results$algae[i + nboot] = coef(summary(lm2boot))[2,"Estimate"]
      results$NoMidge[i + nboot] = coef(summary(lm2boot))[3,"Estimate"]
      results$box2[i + nboot] = coef(summary(lm2boot))[4,"Estimate"]
      results$interaction.algaeNoMidge[i + nboot] = coef(summary(lm2boot))[5,"Estimate"]
      
      #add coefficients from algae and no midge to estimate effect of algae conc on gpp in absence of midge
      results$algae_x_NoMidge[i + nboot] = results$algae[i + nboot] + results$interaction.algaeNoMidge[i + nboot]
      
    }
    return(results)
}


#generate data frame to put all the results
gboot1 <- data.frame(sim = rep(1:nboot, 2), 
                     intercept = NA, 
                     algae = NA, 
                     NoMidge = NA, 
                     box2 = NA, 
                     interaction.algaeNoMidge = NA,
                     algae_x_NoMidge = NA, 
                     day = rep(c("Day 14", "Day 22"), each = nboot))

#peform bootstrapping
gboot1 <- para.boot(data1 = nep14, data2 = nep22,
          response.var = "gpp",
          lm1 = g14log, lm2 = g22log,
          results = gboot1, nboot = nboot)



#plot all coefficients
gboot1 %>% 
  gather(var, val, -sim, -day) %>% 
  group_by(var, day) %>% 
  mutate(mean = mean(val)) %>% 
  ggplot(aes(val, fill = day))+
  geom_histogram(alpha = 0.8, color = "gray50", data = . %>% filter(day == "Day 14"))+
  geom_histogram(alpha = 0.8, color = "gray50", data =. %>% filter(day=="Day 22"))+
  facet_wrap(~var, scales = "free")+
  geom_vline(xintercept = 0, linetype = "dashed")+
  geom_vline(aes(xintercept = mean, col = day), size = 1)


gboot1 %>%
  select(algae, algae_x_NoMidge, day) %>%
  gather(var, val, -day) %>% 
  group_by(var, day) %>% 
  mutate(mean = mean(val),
         var = ifelse(var == 'algae_x_NoMidge', "No Midges", "Midges"),
         var = paste(var, day)) %>% 
  ggplot(aes(val, fill = var))+
  geom_density(alpha = 0.5, color = "gray50")+
  # geom_density(alpha = 0.5, color = "gray50", data =. %>% filter(day=="Day 22"))+
  geom_vline(xintercept = 0, linetype = "dashed")+
  labs(x = "Effect of Inital Algal Concentration on GPP",
       y = "",
       fill = "")+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  scale_fill_brewer(palette = "RdGy")

gboot1 %>%
  select(algae, algae_x_NoMidge, day) %>%
  gather(var, val, -day) %>% 
  group_by(var, day) %>% 
  mutate(mean = mean(val),
         var = ifelse(var == 'algae_x_NoMidge', "No Midges", "Midges")) %>% 
  ggplot(aes(val, fill = var))+
  facet_wrap(~day, scales = "free")+
  geom_density(alpha = 0.5, color = "gray50", data = . %>% filter(var == "Midges"))+
  geom_density(alpha = 0.5, color = "gray50", data =. %>% filter(var == "No Midges"))+
  geom_vline(xintercept = 0, linetype = "dashed")+
  labs(x = "Effect of Inital Algal Concentration on GPP",
       y = "")+
  geom_vline(aes(xintercept = mean, col = var), size = 1)


##Delta Method Version

##NOT WORKING SOME PROBLEMS WITH DEFINITIONS

dmg <- data.frame(day = c("Day 14", "Day 22"), 
           #estimate for the slope of the effect of initial sediment quality on gpp with midges
           algae.estimate = c(coef(g14log)["log(algae_conc2)"], 
                              coef(g22log)["log(algae_conc2)"]), 
           #variance in that estimate
           algae.var = c(vcov(g14log)["log(algae_conc2)", "log(algae_conc2)"],
                         vcov(g22log)["log(algae_conc2)", "log(algae_conc2)"]), 
           #estimate for the interaction between midges and no midges
           interaction.estimate = c(coef(g14log)["log(algae_conc2):midgeNo Midges"],
                                    coef(g22log)["log(algae_conc2):midgeNo Midges"]),
           #variance in that estimate
           interaction.var = c(vcov(g14log)["log(algae_conc2):midgeNo Midges", "log(algae_conc2):midgeNo Midges"],
                               vcov(g22log)["log(algae_conc2):midgeNo Midges", "log(algae_conc2):midgeNo Midges"]), 
           #covariance between the two
           covar = c(vcov(g14log)["log(algae_conc2)", "log(algae_conc2):midgeNo Midges"],
                     vcov(g22log)["log(algae_conc2)", "log(algae_conc2):midgeNo Midges"]))


#not working because variance is being estimated as negative, which can't happen.
dmg <- dmg %>% 
  rowwise() %>% 
  mutate(algae_x_NoMidge = algae.estimate + interaction.estimate,
         algae_x_NoMidge.var = (algae.estimate^2)*algae.var + (interaction.estimate^2)*interaction.var + covar,
         algae.se = sqrt(algae.var),
         interaction.se = sqrt(interaction.var),
         algae_x_NoMidge.se = sqrt(algae_x_NoMidge.var))


dmg %>% 
  select(day, contains("estimate"), algae_x_NoMidge) %>% 
  gather(term, estimate, -day) %>% 
  mutate(term = str_remove(term, ".estimate")) %>% 
  left_join(dmg %>% 
              select(day, contains("se")) %>% 
              gather(term, se, -day) %>% 
              mutate(term = str_remove(term, ".se"))) %>% 
  ggplot(aes(x = term, y = estimate, col = day))+
  geom_hline(yintercept = 0, linetype = "dotted")+
  geom_point(position = position_jitter(width = 0.5, seed = 1), size = 4)+
  geom_pointrange(aes(ymin = estimate-2*se, ymax = estimate+2*se),  position = position_jitter(width = 0.5, seed = 1))+
  coord_flip()
  
#Statistics???
  

#fig 2 with error bars from bootstrapped models



#Supplemental Figure
nep %>% 
  gather(metabolism, value, resp, nep, gpp) %>% 
  mutate(day = day+rnorm(n(), 0.001), #add jitter
         metabolism = toupper(ifelse(metabolism == "resp", "er", metabolism)),
         metabolism = fct_reorder(metabolism, desc(value))) %>%  
  ggplot(aes(x = day, y = value, col = algae_conc2))+
  facet_grid(metabolism~midge, scales = "free")+
  geom_hline(yintercept = 0)+
  geom_point()+
  geom_line(aes( group = coreid), alpha = 0.4)+
  geom_smooth(method = "lm", col = "black", se = FALSE)+
  theme(legend.position = "bottom")+
  labs(y = expression("Change in" ~O[2] (g~m^{-2}~hr^{-1})),
       x = "Incubation Day",
       color = "Initial Sediment Quality")+
  scale_color_viridis_c(trans = "log10")+
  scale_x_continuous(breaks = c(14, 22))

nep %>% 
  mutate(day = paste("Day", day)) %>% 
  gather(metabolism, value, resp, nep, gpp) %>% 
  mutate(metabolism = toupper(metabolism),
         metabolism = ifelse(metabolism == "RESP", "ER", metabolism),
         metabolism = fct_reorder(metabolism, -value)) %>% 
  ggplot(aes(x = algae_conc2, y = value, fill = midge))+
  geom_hline(yintercept = 0, linetype = "dotted")+
  facet_grid(metabolism~day, scale = "free_y", switch = "y")+
  geom_point(size  = 2, alpha = 0.7, shape = 21)+
  geom_smooth(aes(color = midge), method = "lm", size = 0.75, se = FALSE)+
  scale_fill_manual(values = c("black", "gray60"))+
  scale_color_manual(values = c("black", "gray60"))+
  labs(x = "Initial Sediment Quality",
       y = expression(g~O[2]~m^{-2}~hr^{-1}),
       color = "",
       fill = "")+
  scale_x_log10()+
  theme(strip.placement = "outside")

nep %>% 
  ggplot(aes(x = -resp, y = gpp, col = algae_conc2))+
  geom_point()+
  facet_wrap(~day)+
  geom_abline(slope = 1)+
  coord_equal()

##Midge Effects on chlorophyll

#This is especially sloppy because I might remove all of this. 
lm(NDVI_R~log(algae_conc2)*midge+box, data = ndvi %>% left_join(meta) %>% filter(day ==14)) %>% summary()
lm(NDVI_R~log(algae_conc2)*midge+box, data = ndvi %>% left_join(meta) %>% filter(day ==22)) %>% summary()

ndvi %>% 
  left_join(meta) %>% 
  ggplot(aes(x = algae_conc2, y = NDVI_R, fill = midge))+
  geom_point(size = 2, shape = 21)+
  facet_wrap(~paste("Day", day))+
  # geom_smooth(aes(col = midge), method = "lm", se = FALSE)+
  scale_fill_manual(values = c("black", "gray60"))+
  scale_color_manual(values = c("black", "gray60"))+
  labs(x = "Initial Sediment Quality",
       y = "NDVI")+
  scale_x_log10()

#supplement
ndvi %>%
  mutate(coreid = as.character(coreid)) %>% 
  left_join(nep) %>% 
  ggplot(aes(x = gpp, y = NDVI_R, col = algae_conc2))+
  geom_point()+
  facet_grid(midge~day)+
  geom_smooth(method = "lm", se = FALSE, col = "black", size = 0.7)+
  scale_color_viridis_c(trans = "log10")+
  labs(x = expression(GPP~(g^{2}~O[2]~m^{-2}~hr^{-1})),
       y = "NDVI",
       color = "Initial Sediment Quality")

ndvi %>% 
  left_join(meta) %>% 
  ggplot(aes(x = factor(day), y = NDVI_R, col = algae_conc2, group = factor(algae_conc)))+
  geom_jitter(width = 0.2)+
  facet_wrap(~midge)+
  geom_smooth(method = "lm", se = FALSE, size = 0.7)+
  scale_color_viridis_c(trans = "log10")+
  labs(x = "Day",
       y = "NDVI",
       color = "Initial Sediment Quality")



#====Question 2: Sediment Effects on Midges =====
##====Q2.1: Number of Midges Present====
cc <- cc_raw %>% 
  left_join(meta) %>% 
  mutate(coreid = as.character(coreid),
         box = as.character(box))

cc14 <- cc %>% 
  filter(day == 14)

nm14log <- glm(live_tt~log(algae_conc2)*midge+box, 
               data = cc14,
               family = quasipoisson())

summary(nm14log)
# plot(nm14log)

cc22 <- cc %>% 
  filter(day == 22)

nm22log <- glm(live_tt~log(algae_conc2)*midge+box, 
               data = cc22,
               family = quasipoisson())

summary(nm22log)
# plot(nm22log)

#====Figure 3: Number of Midges====
cc %>% 
  left_join(nep) %>% 
  filter(day %in% c(14, 22)) %>% 
  mutate(day = paste("Day", day)) %>% 
  ggplot(aes(x = algae_conc2, y = live_tt, fill = gpp))+
  facet_grid(midge~day, scales = "free")+
  geom_jitter(shape = 21, size = 2, alpha = 0.6, width = 0.0005, height = 0.3)+
  geom_smooth(method = "lm", se = FALSE, size = 0.6, col = "black")+
  scale_fill_viridis_c()+
  labs(y = "Live Tanytarsini",
       x = "Initial Sediment Quality",
       fill = expression(GPP(g~O[2]~m^{-2}~hr^{-1})~" "))+
  scale_x_log10()

cc %>% 
  left_join(nep) %>% 
  filter(day %in% c(14, 22)) %>% 
  mutate(day = paste("Day", day)) %>% 
  ggplot(aes(x = algae_conc2, y = live_tt, col = gpp, shape = day, linetype = day))+
  facet_wrap(~midge, scales = "free")+
  geom_jitter(size = 2, alpha = 0.6, width = 0.0005, height = 0.3)+
  geom_smooth(method = "lm", se = FALSE, size = 0.6, col = "black")+
  scale_color_viridis_c(trans = "log10")+
  scale_shape_manual(values = c(21,16))+
  scale_linetype_manual(values = c("dashed", "solid"))+
  labs(y = "Live Tanytarsini",
       x = "Initial Sediment Quality",
       fill = expression(GPP(g~O[2]~m^{-2}~hr^{-1})~" "))+
  scale_x_log10()+
  lims(y = c(0, NA))

#====Q2.1.a Bootstrapping Coefficients====
#generate data frame to put all the results
nmboot <- data.frame(sim = rep(1:nboot, 2), 
                     intercept = NA, 
                     algae = NA, 
                     NoMidge = NA, 
                     box2 = NA, 
                     interaction.algaeNoMidge = NA,
                     algae_x_NoMidge = NA, 
                     day = rep(c("Day 14", "Day 22"), each = nboot))

#unable to bootstrap quasi-models. The estimates should be the same if I use family = nonquasi, so we'll do that
nm142 <- update(nm14log, .~., family = poisson())
nm222 <- update(nm22log, .~., family = poisson())


#peform bootstrapping
nmboot <- para.boot(data1 = cc14, data2 = cc22, 
          response.var = "live_tt", 
          lm1 = nm142, lm2 = nm222, 
          results = nmboot, nboot = nboot)


#plot all coefficients
nmboot %>% 
  gather(var, val, -sim, -day) %>% 
  group_by(var, day) %>% 
  mutate(mean = mean(val)) %>% 
  ggplot(aes(val, fill = day))+
  geom_histogram(alpha = 0.8, color = "gray50", data = . %>% filter(day == "Day 14"))+
  geom_histogram(alpha = 0.8, color = "gray50", data =. %>% filter(day=="Day 22"))+
  facet_wrap(~var, scales = "free")+
  geom_vline(xintercept = 0, linetype = "dashed")+
  geom_vline(aes(xintercept = mean, col = day), size = 1)


nmboot %>%
  select(algae, algae_x_NoMidge, day) %>%
  gather(var, val, -day) %>% 
  group_by(var, day) %>% 
  mutate(mean = mean(val),
         var = ifelse(var == 'algae_x_NoMidge', "No Midges", "Midges")) %>% 
  ggplot(aes(val, fill = day))+
  facet_wrap(~var, scales = "free")+
  geom_density(alpha = 0.5, color = "gray50", data = . %>% filter(day == "Day 14"))+
  geom_density(alpha = 0.5, color = "gray50", data =. %>% filter(day=="Day 22"))+
  geom_vline(xintercept = 0, linetype = "dashed")+
  labs(x = "Effect of Inital Algal Concentration on Number of Live Midges",
       y = "")+
  geom_vline(aes(xintercept = mean, col = day), size = 1)


nmboot %>%
  select(algae, algae_x_NoMidge, day) %>%
  gather(var, val, -day) %>% 
  group_by(var, day) %>% 
  mutate(mean = mean(val),
         var = ifelse(var == 'algae_x_NoMidge', "No Midges", "Midges")) %>% 
  ggplot(aes(val, fill = var))+
  facet_wrap(~day, scales = "free")+
  geom_density(alpha = 0.5, color = "gray50", data = . %>% filter(var == "Midges"))+
  geom_density(alpha = 0.5, color = "gray50", data =. %>% filter(var == "No Midges"))+
  geom_vline(xintercept = 0, linetype = "dashed")+
  labs(x = "Effect of Inital Algal Concentration on Number of Live Midges",
       y = "")+
  geom_vline(aes(xintercept = mean, col = var), size = 1)



##====Q2.2: Midge Development====
cm <- cm_raw %>% 
  left_join(meta) %>% 
  mutate(coreid = as.character(coreid),
         box = as.character(box),
         body_size = body_size/1000) #convert from um to mm


#find starting midge stages
startingprop <- {cm %>% 
    #only tanytarsini
  filter(species_name == "tt") %>% 
  #set 1 if instar>2
  mutate(instar = ifelse(instar=="3"|instar=="4", 1, 0)) %>% 
  group_by(coreid, day, algae_conc, instar, midge) %>% 
  #count number of greater than 2nd instar
  summarise(count = n()) %>% 
  group_by(coreid) %>% 
  #get proportion within a mesocosm
  mutate(prop  = count/sum(count),
         live_tt = sum(count)) %>% 
  filter(instar == 1, 
         is.na(algae_conc))}$prop


p3 <- cm %>% 
  filter(species_name == "tt") %>% 
  mutate(instar = ifelse(instar=="3"|instar=="4", 1, 0)) %>% 
  group_by(coreid, day, algae_conc, algae_conc2, instar, midge, box) %>% 
  summarise(count = n()) %>% 
  group_by(coreid) %>% 
  mutate(prop  = count/sum(count),
         live_tt = sum(count)) %>% 
  filter(instar == 1)

#====Figure 4: Instar====
p3 %>% 
  filter(!is.na(algae_conc)) %>% 
  mutate(day = paste("Day", day)) %>% 
  ggplot(aes(y = prop, x = algae_conc2, fill = live_tt))+
  geom_point(shape = 21, size = 2, alpha = 0.7)+
  geom_hline(yintercept = startingprop)+
  facet_grid(midge~day)+
  scale_fill_viridis_c(option = "plasma")+
  lims(y = c(0,1))+
  labs(x = expression(GPP~(g^{2}~O[2]~m^{-2}~hr^{-1})),
       y = "Proportion of Larval Midges\n3rd instar or higher",
       fill = "Number of Larvae")+
  theme(legend.position = "bottom")+
  scale_x_log10()

#subset for day 14
p314 <- p3 %>% 
  filter(day == 14)

prop14log <- glm(prop~log(algae_conc2)*midge+box, 
                 data = p314, 
                 family = quasibinomial())

summary(prop14log)
# plot(prop14log)

#subset for day 22
p322 <- p3 %>% 
  filter(day == 22)

prop22log <- glm(prop~log(algae_conc2)*midge+box, 
                 data = p322, 
                 family = quasibinomial())
summary(prop22log)
# plot(prop22log)

#====Q2.2 Bootstrap====

##PROBLEMS!!!

#generate data frame to put all the results
propboot <- data.frame(sim = rep(1:nboot, 2), 
                     intercept = NA, 
                     algae = NA, 
                     NoMidge = NA, 
                     box2 = NA, 
                     interaction.algaeNoMidge = NA,
                     algae_x_NoMidge = NA, 
                     day = rep(c("Day 14", "Day 22"), each = nboot))

#update quasi-models
prop142 <- update(prop14log, .~., family = binomial())
prop222 <- update(prop22log, .~., family = binomial())

#peform bootstrapping
propboot <- para.boot(data1 = p314, data2 = p322, 
          response.var = "prop", 
          lm1 = prop142, lm2 = prop222, 
          results = propboot, nboot = nboot)


#plot all coefficients
propboot %>% 
  gather(var, val, -sim, -day) %>% 
  group_by(var, day) %>% 
  mutate(mean = mean(val)) %>% 
  ggplot(aes(val, fill = day))+
  geom_histogram(alpha = 0.8, color = "gray50", data = . %>% filter(day == "Day 14"))+
  geom_histogram(alpha = 0.8, color = "gray50", data =. %>% filter(day=="Day 22"))+
  facet_wrap(~var, scales = "free")+
  geom_vline(xintercept = 0, linetype = "dashed")+
  geom_vline(aes(xintercept = mean, col = day), size = 1)


propboot %>%
  select(algae, algae_x_NoMidge, day) %>%
  gather(var, val, -day) %>% 
  group_by(var, day) %>% 
  mutate(mean = mean(val),
         var = ifelse(var == 'algae_x_NoMidge', "No Midges", "Midges")) %>% 
  ggplot(aes(val, fill = day))+
  facet_wrap(~var, scales = "free")+
  geom_density(alpha = 0.5, color = "gray50", data = . %>% filter(day == "Day 14"))+
  geom_density(alpha = 0.5, color = "gray50", data =. %>% filter(day=="Day 22"))+
  geom_vline(xintercept = 0, linetype = "dashed")+
  labs(x = "Effect of Inital Algal Concentration on GPP",
       y = "")+
  geom_vline(aes(xintercept = mean, col = day), size = 1)


propboot %>%
  select(algae, algae_x_NoMidge, day) %>%
  gather(var, val, -day) %>% 
  group_by(var, day) %>% 
  mutate(mean = mean(val),
         var = ifelse(var == 'algae_x_NoMidge', "No Midges", "Midges")) %>% 
  ggplot(aes(val, fill = var))+
  facet_wrap(~day, scales = "free")+
  geom_density(alpha = 0.5, color = "gray50", data = . %>% filter(var == "Midges"))+
  geom_density(alpha = 0.5, color = "gray50", data =. %>% filter(var == "No Midges"))+
  geom_vline(xintercept = 0, linetype = "dashed")+
  labs(x = "Effect of Inital Algal Concentration on GPP",
       y = "")+
  geom_vline(aes(xintercept = mean, col = var), size = 1)



##====Q2.3: Midge Length====
#subset to 3rd instar tanytarsini and to not include the stocked midges
l3 <- cm %>% 
  filter(instar == 3,
         day!=0)

#subset day 14
l314 <- l3 %>% 
  filter(day == 14)

bl14log <- lmer(body_size~log(algae_conc2)*midge+box+
                  (1|coreid), 
                data = l314) 

summary(bl14log)
# plot(bl14log)

Anova(bl14log, type = "3", test.statistic = "F")
Anova(bl14log, type = "2", test.statistic = "F")

#subset day 22
l322 <- l3 %>% 
  filter(day == 22)

bl22log <- lmer(body_size~log(algae_conc2)*midge+box +
                  (1|coreid), 
                data = l322)
summary(bl22log)
# plot(bl22log)

Anova(bl22log, type = "3", test.statistic = "F")
Anova(bl22log, type = "2", test.statistic = "F")


#====Figure 5: Body Length====
l3 %>% 
  left_join(nep) %>%
  mutate(day = paste("Day", day)) %>% 
  ggplot(aes(x = algae_conc2, y = body_size, fill = gpp))+
  facet_grid(midge~day, scales = "free_x")+
  geom_point(alpha = 0.2, shape = 21, size = 2)+
  geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.7)+
  scale_fill_viridis_c()+
  labs(x = "Initial Sediment Conditions",
       y = "Body Length of \nThird Instar Tanytarsini (mm)",
       fill = expression(GPP~(g^{2}~O[2]~m^{-2}~hr^{-1})~" "))+
  theme(legend.position = "bottom")+
  scale_x_log10()


#====Q2.3a: Bootstrap====
#generate data frame to put all the results
blboot <- data.frame(sim = rep(1:nboot, 2), 
                     intercept = NA, 
                     algae = NA, 
                     NoMidge = NA, 
                     box2 = NA, 
                     interaction.algaeNoMidge = NA,
                     algae_x_NoMidge = NA, 
                     day = rep(c("Day 14", "Day 22"), each = nboot))

#peform bootstrapping
blboot <- para.boot(data1 = l314 %>% filter(!is.na(body_size)), data2 = l322 %>% filter(!is.na(body_size)), 
          response.var = "live_tt", 
          lm1 = bl14log, lm2 = bl22log, 
          results = blboot, nboot = nboot)


#plot all coefficients
blboot %>% 
  gather(var, val, -sim, -day) %>% 
  group_by(var, day) %>% 
  mutate(mean = mean(val)) %>% 
  ggplot(aes(val, fill = day))+
  geom_histogram(alpha = 0.8, color = "gray50", data = . %>% filter(day == "Day 14"))+
  geom_histogram(alpha = 0.8, color = "gray50", data =. %>% filter(day=="Day 22"))+
  facet_wrap(~var, scales = "free")+
  geom_vline(xintercept = 0, linetype = "dashed")+
  geom_vline(aes(xintercept = mean, col = day), size = 1)


blboot %>%
  select(algae, algae_x_NoMidge, day) %>%
  gather(var, val, -day) %>% 
  group_by(var, day) %>% 
  mutate(mean = mean(val),
         var = ifelse(var == 'algae_x_NoMidge', "No Midges", "Midges")) %>% 
  ggplot(aes(val, fill = day))+
  facet_wrap(~var, scales = "free")+
  geom_density(alpha = 0.5, color = "gray50", data = . %>% filter(day == "Day 14"))+
  geom_density(alpha = 0.5, color = "gray50", data =. %>% filter(day=="Day 22"))+
  geom_vline(xintercept = 0, linetype = "dashed")+
  labs(x = "Effect of Inital Algal Concentration on Body Length",
       y = "")+
  geom_vline(aes(xintercept = mean, col = day), size = 1)


blboot %>%
  select(algae, algae_x_NoMidge, day) %>%
  gather(var, val, -day) %>% 
  group_by(var, day) %>% 
  mutate(mean = mean(val),
         var = ifelse(var == 'algae_x_NoMidge', "No Midges", "Midges")) %>% 
  ggplot(aes(val, fill = var))+
  facet_wrap(~day, scales = "free")+
  geom_density(alpha = 0.5, color = "gray50", data = . %>% filter(var == "Midges"))+
  geom_density(alpha = 0.5, color = "gray50", data =. %>% filter(var == "No Midges"))+
  geom_vline(xintercept = 0, linetype = "dashed")+
  labs(x = "Effect of Inital Algal Concentration on Body Length",
       y = "")+
  geom_vline(aes(xintercept = mean, col = var), size = 1)




##====Q2.3b Midge Length Alternate====
#3rd instars averaged
bl14log2 <- lm(body_size~log(algae_conc2)*midge+box,
               weights = live_tt,
                data = l314 %>% 
                  group_by(algae_conc2, midge, box) %>% 
                summarise(body_size = meanna(body_size)) %>% 
                 left_join(cc)) 

summary(bl14log2)
plot(bl14log2)


bl22log2 <- lm(body_size~log(algae_conc2)*midge+box, 
               weights = live_tt,
                data = l322 %>% 
                 group_by(algae_conc2, midge, box) %>% 
                 summarise(body_size = meanna(body_size)) %>% 
                 left_join(cc))
summary(bl22log2)
plot(bl22log2)


l3 %>% 
  left_join(nep) %>%
  group_by(coreid, gpp, algae_conc2, day, midge) %>% 
  summarise(bs_se = sd(body_size, na.rm = TRUE), 
            body_size = meanna(body_size)) %>% 
  mutate(day = paste("Day", day)) %>% 
  ggplot(aes(x = algae_conc2, y = body_size, fill = gpp))+
  facet_grid(midge~day, scales = "free_x")+
  geom_point(alpha = 0.7, shape = 21, size = 2)+
  geom_errorbar(aes(ymin = body_size-bs_se, ymax = body_size + bs_se), 
                col = "black", width = 0)+
  geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.7)+
  scale_fill_viridis_c()+
  labs(x = "Initial Sediment Conditions",
       y = "Body Length of \nThird Instar Tanytarsini (mm)",
       fill = expression(GPP~(g^{2}~O[2]~m^{-2}~hr^{-1})~" "))+
  theme(legend.position = "bottom")+
  scale_x_log10()

#body size for all instars averaged
bl14log3 <- lm(body_size~log(algae_conc2)*midge+box, 
               weights = live_tt,
               data = cm %>% 
                 filter(day == 14) %>% 
                 group_by(algae_conc2, midge, box) %>% 
                 summarise(body_size = meanna(body_size)) %>% 
                 left_join(cc)) 

summary(bl14log3)
plot(bl14log3) 


bl22log3 <- lm(body_size~log(algae_conc2)*midge+box, 
               weights = live_tt,
               data = cm %>% 
                 filter(day == 22) %>% 
                 group_by(algae_conc2, midge, box) %>% 
                 summarise(body_size = meanna(body_size)) %>% 
                 left_join(cc))
summary(bl22log3)
plot(bl22log3)

cm %>% 
  filter(day!=0) %>% 
  left_join(nep) %>%
  group_by(coreid, gpp, algae_conc2, day, midge) %>% 
  summarise(bs_se = sd(body_size, na.rm = TRUE), 
            body_size = meanna(body_size)) %>% 
  mutate(day = paste("Day", day)) %>% 
  ggplot(aes(x = algae_conc2, y = body_size, fill = gpp))+
  facet_grid(midge~day, scales = "free_x")+
  geom_point(alpha = 0.7, shape = 21, size = 2)+
  geom_errorbar(aes(ymin = body_size-bs_se, ymax = body_size + bs_se), 
                col = "black", width = 0)+
  geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.7)+
  scale_fill_viridis_c()+
  labs(x = "Initial Sediment Conditions",
       y = "Body Length of \nThird Instar Tanytarsini (mm)",
       fill = expression(GPP~(g^{2}~O[2]~m^{-2}~hr^{-1})~" "))+
  theme(legend.position = "bottom")+
  scale_x_log10()

#Separating midge no midge
lm(body_size~log(algae_conc2)*day+box, 
   weights = live_tt,
   data = cm %>% 
     filter(midge == "Midges") %>% 
     group_by(algae_conc2, midge, box) %>% 
     summarise(body_size = meanna(body_size)) %>% 
     left_join(cc)) %>% summary()

lm(body_size~log(algae_conc2)*day+box, 
   weights = live_tt,
   data = cm %>% 
     filter(midge == "No Midges") %>% 
     group_by(algae_conc2, midge, box) %>% 
     summarise(body_size = meanna(body_size)) %>% 
     left_join(cc)) %>% summary()




##====Q2 SUPPLEMENT1 : Did midges increase in body length?====
growth <- lm(body_size~day, data = cm %>% filter(species_name == "tt"))
summary(growth)
plot(growth)

growth3 <- lm(body_size~day*factor(instar), data = cm %>% filter(species_name == "tt"))
summary(growth3)
plot(growth3)

cm %>% group_by(day) %>% 
  summarise(mean = mean(body_size, na.rm = TRUE),
            se = sd(body_size, na.rm = TRUE)) %>% 
  mutate(xmin = day-2,
         xmax = day+2,
         ymin = mean-se,
         ymax = mean+se) %>% 
  ggplot(aes(x = day, y = mean))+
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), alpha = 0.5)+
  geom_jitter(aes(y = body_size, col = factor(instar)), width = 2, height = 0, alpha = 0.5, data = cm %>% filter(species_name == "tt", !is.na(instar)))+
  geom_smooth(aes(y = body_size), size = 0.7, color = "red", method = "lm", se = FALSE, data = cm)+
  geom_segment(aes(yend = mean, x = xmin, xend = xmax), size = 1, col = "blue")+
  geom_text(y = 5.5, x = 0, label = paste0("y = ", round(coef(growth)[2], 2), "x + ", round(coef(growth)[1], 2), "\n R2 = ", round(summary(growth)$r.squared, 2)))+
  labs(x = "Day",
       y = "Tanytarsini Length (mm)",
       color = "Instar")

growth4 <- lm(weight(body_size)~day, data = cm %>% filter(species_name == "tt"))
summary(growth4)
plot(growth4)


#====Q2.4: Midge Movement====
mmv1 <- glm(live_tt~log(algae_conc2)*factor(day) + box, 
    data = cc %>% filter(midge=="No Midges"), 
    family = quasipoisson())
summary(mmv1)

# mmv2 <- glm(live_tt~log(gpp)*factor(day) + box, 
#             data = cc %>% filter(midge=="No Midges") %>% left_join(nep), 
#             family = quasipoisson())
# summary(mmv2)
# 
# 
# anova(mmv1, mmv2)

#====Question 3: How Much do midges rely on algae?====
#convert midges to biomass following Lindegaard and Jonasson 1979
#AFDW in mg^(1/3) = 0.0442 + 0.0879 * length in mm 
weight <- function(l){
  (0.442 + 0.0879 * l )^3 
}


#calculate carbon using Amanda's Isotope Data
dwc = 0.5 #plae holder


#increment Summation method to estimate production 
incsumprod <- function(n1, n2, wt1, wt2, deltaT){
  g = wt2-wt1 #individual growth
  Nbar = (n1 + n2)/2 #Average abundance
  intp = Nbar*g #production over the interval
  Pd = intp/deltaT #daily production
  list(g = g, Nbar = Nbar, intp = intp, Pd = Pd)
}

#Calculate standard errors on midge productivity

#step 1: calculate partial derivatives for delta method
deriv(~((n1+n2)/2)*(wt2-wt1)/deltaT,
      c("n1", "wt2", "wt1")) #we assume perfect knowledge of nt and t
#step 2:  
#Delta method to propogate error:
#Sy = sqrt((dy/dx)^2*s^2x + (dz/dx)^2*s^2z + .... + cov)
prodse <- function(n1, n2, wt1, wt2, deltaT, n1se, w1se, w2se){
  expr2 = (n1+n2)/2
  expr3 = wt2-wt1
  expr9 = expr2/deltaT
  dn1 = 1/2*expr3/deltaT
  dw2 = expr9
  dw1 = expr9
  P.se = sqrt(dw2^2*w2se^2 + dw1^2*w2se^2)
  return(P.se)
}



#calculate midge productivity
midgeprod <- cm %>% 
  left_join(cc) %>% 
  filter(species_name == "tt",
         day!=0) %>% 
  group_by(coreid) %>% 
  mutate(w = weight(body_size)) %>% #estimate dry mass
  group_by(coreid, day, box, algae_conc) %>% 
  summarise(wt = meanna(w), #average biomass of a midge (mg)
            nt = unique(live_tt), #number of individuals in the mesocosm
            wt.se = sd(w, na.rm = TRUE)/sqrt(nt)) %>% 
  # combine with initial values
  bind_cols(cm %>% 
              filter(species_name == "tt",
                     day == 0) %>% 
              group_by(coreid) %>% 
              add_count() %>% 
              mutate(w = weight(body_size),
                     c = w*dwc,
                     live_tt = n) %>%
              ungroup %>% 
              group_by(sampledate) %>%
              summarise(n= n(),
                        wt1 = mean(w), #weight at t = 0
                        n1 = mean(live_tt), #number at t = 0
                        wt1.se = sd(w)/sqrt(n),
                        nt1.se = sd(live_tt)/sqrt(3)) %>% 
              select(-sampledate)) %>% 
  ungroup %>% 
  mutate(Pd = incsumprod(n1 = n1, n2 = nt, wt1 = wt1, wt2 = wt, deltaT = day)$Pd, #production in mg,
         g = incsumprod(n1 = n1, n2 = nt, wt1 = wt1, wt2 = wt, deltaT = day)$g,
         Pd.se = prodse(n1 = n1, n2 = nt, wt1 = wt1, wt2 = wt, deltaT = day, n1se = nt1.se, w1se = wt1.se, w2se = wt.se),
         Pdc.se = Pd.se*dwc*1000,
         Pdc = Pd*dwc*1000) #convert from weight to c (mg) then convert to micrograms 
 



#compare measurements of growth (lm) to these estimates
midgeprod %>% 
  group_by(day) %>% 
  mutate(meanw = mean((wt-wt1)/day)) %>% 
  ggplot(aes((wt-wt1)/day, fill = factor(day)))+
  geom_histogram(bins = 15, alpha = 0.7)+
  geom_vline(xintercept = 0, linetype = "dotted")+
  geom_vline(xintercept = coef(growth4)[2], linetype = "dashed")+
  geom_vline(aes(xintercept = meanw, col = factor(day)))

##Estimate daily production of primary producers in C
#converting g O2--> C
#1:1 molar mass
12/32

pq <- 0.375

#calculate area of mesocosm
mesor <- 15 #mm
area <- mesor^2*pi #mm^2

#calculate algal productivity
algaeprod <- nep %>%
  mutate(gpp_c = gpp*pq, #g m^-2 h^-1
         #multiply by 1e6 to convert g to ug and divide by 1e6 to get mm^-2 multiply by 24 to get hourly 
         gpp_cm = gpp_c*area, #hourly primary production of C (micrograms) in a mesocosm
         resp_cm = resp*pq*area, #hourly respiration
         nep_cm = nep*pq*area, #hourly nep 
         net_daily = 18*nep_cm + 6*resp_cm, #nep over course of day
         gpp_daily = 18*gpp_cm, #daily GPP
         resp_daily = 24*resp_cm) #daily respiration


apm <- algaeprod %>% 
  filter(midge == "Midges")

#Set values for transfer efficiencies
AEa = 1 #arbitrary value set for the assimilation of fixed carbon by algae
NPEm = 0.54 #Net Production efficiency for T. gracilentus at Myvatn Lindegaard 1994
FA = 1 #Fraction of assimilated material in midge attributed to algae
AEm = c(0.1, 0.5, 1) #Assimilation efficiency


#Im = (Pm*FA)/(NPE*AE), Where Im = ingestion by midges, Pm = midge productivity, FA = Fraction of algae assimilated, NPE = Net production efficiency, AE = assimilation efficiency
#when Im = Pa, Where Pa = algal productivity 
#Pa/(Pm*FA) = NPE*AE

#confirming the slope = NPE*AE
data.frame(GPP = 0:max(apm$gpp_daily)) %>% 
  crossing(AE = c(0.1, 0.5, 1)) %>%
  mutate( NPP = GPP*AEa,
          Pdc = (NPP*NPEm*AE)/FA) %>% 
  lm(Pdc~0+NPP:factor(AE), data = .)



#====Figure 6: Midge productivity and primary production=====
#combine estimates of midge production with algae production
production <- midgeprod %>%  
  full_join(algaeprod %>% 
              group_by(coreid, algae_conc2, midge) %>% 
              # select(-day) %>% 
              mutate(NPP = gpp_daily*AEa,
                     NPP.mean = meanna(NPP)) %>% #multiply GPP by AEa (the efficiency at which algae convert fixed carbon into new biomass)
              summarise(NPP.sd = sd(NPP),
                        NPP = meanna(NPP))) #average by day


#slopes and lines for I/P figure
IPs = data.frame(AEm = AEm, 
                 NPEm = NPEm, 
                 x = max(production$NPP)+5,
                 slope = AEm*NPEm,
                 label = paste("AE[m]==", format(AEm, nsmall = 1)),
                 label2 = "I[m]/P[a] == 1") %>% 
  mutate(y = x*NPEm*AEm)


production %>% 
  mutate(Pdc = ifelse(is.na(Pdc), 0, Pdc),
         Pdc.se = ifelse(is.na(Pdc), 0, Pdc.se)) %>% 
  ggplot(aes(x = NPP, y = Pdc))+
  facet_wrap(~midge)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_errorbar(aes(ymin = Pdc-Pdc.se, ymax = Pdc+Pdc.se, color = algae_conc))+
  geom_point(aes(color = algae_conc, shape = factor(day)), fill = "white")+
  geom_segment(aes(x = 0, y = 0, xend = max(production$NPP)+5, yend =  slope * (max(production$NPP)+5)), data = IPs)+
  geom_text(aes(x = x, y = y, label = label), parse = TRUE, hjust = 0, vjust = 1, color = "black", data = IPs)+
  geom_text(aes(x = x, y = y, label = label2), parse = TRUE, hjust = 0,vjust = 0, color = "black", data = IPs)+
  labs(x = expression("Algal Productivity \u03BCg C "~d^{-1}),
       y = expression("Midge Productivity \u03BCg C "~d^{-1}),
       color = "Initial Sediment Quality",
       shape = "Day")+
  scale_color_gradient(low = "goldenrod2", high = "darkgreen")+
  scale_shape_manual(values = c(21, 16))+
  lims(x = c(0, 200),
       y = c(NA, 100))+
  coord_equal()+
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         shape = guide_legend(title.position="top", title.hjust = 0.5))


production %>% 
  mutate(growth = g*dwc*1000,
         gd = growth/day) %>% 
  # filter(midge == "Midges") %>% 
  # lm(gd~NPP, data = .) %>% summary()
  ggplot(aes(x = NPP, y = gd))+
  geom_point()+
  facet_wrap(~midge)+
  geom_smooth(method = "lm")+
  geom_hline(yintercept = 0)


midgeprod %>% 
  mutate(growth = g*dwc*1000,
         gd = growth/day) %>% 
  left_join(algaeprod %>% 
              select(coreid, algae_conc2, midge, day, gpp_daily) %>% 
              group_by(coreid, algae_conc2) %>% 
              mutate(NPP = gpp_daily*AEa,
                     day = paste0("NPP_", day)) %>%
              select(-gpp_daily) %>% 
              spread(day, NPP) %>% 
              rowwise() %>% 
              mutate(NPP = meanna(c(NPP_14, NPP_22)), 
                     deltaNPP = NPP_22-NPP_14)) %>% 
  filter(!is.na(NPP_22)) %>% 
  ggplot(aes(x = gd, y = deltaNPP))+
  facet_wrap(~midge, scales = "free")+
  geom_point()

#in progress: looking at changes from first incubation to second 
cm %>% 
  left_join(cc) %>% 
  filter(species_name == "tt",
         day!=22) %>% 
  group_by(coreid) %>% 
  mutate(w = weight(body_size)) %>% #estimate dry mass
  group_by(coreid, day, box, algae_conc) %>% 
  summarise(wt = meanna(w), #average biomass of a midge (mg)
            nt = unique(live_tt), #number of individuals in the mesocosm
            wt.se = sd(w, na.rm = TRUE)/sqrt(nt)) %>% 
  # combine with initial values
  bind_cols(cm %>% 
              left_join(cc) %>% 
              filter(species_name == "tt",
                     day == 14) %>% 
              mutate(w = weight(body_size),
                     c = w*dwc) %>%
              group_by(algae_conc, midge) %>%
              summarise(samps = n_distinct(coreid),
                        n= unique(live_tt),
                        wt1 = meanna(w), #weight at t = 0
                        n1 = meanna(n), #number at t = 0
                        wt1.se = sd(w)/sqrt(n),
                        nt1.se = sd(n)/sqrt(samps)) %>% 
              select(-sampledate)) %>% 
  ungroup %>% 
  mutate(Pd = incsumprod(n1 = n1, n2 = nt, wt1 = wt1, wt2 = wt, deltaT = day)$Pd, #production in mg,
         g = incsumprod(n1 = n1, n2 = nt, wt1 = wt1, wt2 = wt, deltaT = day)$g,
         Pd.se = prodse(n1 = n1, n2 = nt, wt1 = wt1, wt2 = wt, deltaT = day, n1se = nt1.se, w1se = wt1.se, w2se = wt.se),
         Pdc.se = Pd.se*dwc*1000,
         Pdc = Pd*dwc*1000) #convert from weight to c (mg) then convert to micrograms 



midgeprod %>% 
  left_join(algaeprod %>% 
              select(coreid, algae_conc2, day, gpp_daily) %>% 
              group_by(coreid, algae_conc2) %>% 
              mutate(NPP = gpp_daily*AEa,
                     day = paste0("NPP_", day)) %>%
              select(-gpp_daily) %>% 
              spread(day, NPP) %>% 
              rowwise() %>% 
              mutate(NPP = meanna(c(NPP_14, NPP_22)))) %>% 
  ggplot(aes(y = Pdc, x = NPP))+
  geom_errorbar(aes(xmin = NPP_14, xmax = NPP_22, col = algae_conc2))+
  geom_errorbar(aes(ymin = Pdc-Pdc.se, ymax = Pdc.se+Pdc, col = algae_conc2))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_point(aes(col = algae_conc2))+
  geom_segment(aes(x = 0, y = 0, xend = max(production$NPP)+5, yend =  slope * (max(production$NPP)+5)), data = IPs)+
  labs(x = expression("Algal Productivity \u03BCg C "~d^{-1}),
       y = expression("Midge Productivity \u03BCg C "~d^{-1}),
       color = "Initial Sediment Quality",
       shape = "Day")+
  scale_color_gradient(low = "goldenrod2", high = "darkgreen", trans = "log10")+
  scale_fill_gradient(low = "goldenrod2", high = "darkgreen", trans = "log10", guide = FALSE)+
  facet_wrap(~day)+
  coord_equal()+
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         shape = guide_legend(title.position="top", title.hjust = 0.5))



#What is the AE if they only eat algae??
prodmod1 <- lm(Pdc~NPP, 
              data = production)

coef(prodmod1)[2]/NPEm

prodmod2<- lm(Pdc~NPP, 
               weights = 1/Pdc.se^2,
               data = production)
coef(prodmod2)[2]/NPEm

prodmod3 <- lm(Pdc~NPP:factor(day), 
   data = production)

coef(prodmod3)[2]/NPEm
coef(prodmod3)[3]/NPEm

prodmod4 <- lm(Pdc~NPP:factor(day), 
               weights = 1/Pdc.se^2,
               data = production)

coef(prodmod4)[2]/NPEm
coef(prodmod4)[3]/NPEm


production %>% 
  mutate(lower = Pdc-2*Pdc.se) %>% 
  filter(lower/(NPEm*1)>NPP)

#====Q3: Supplement====
# using day of incubation
sim%>% 
  ggplot(aes(x = NPP, y = Pdc))+
  facet_wrap(~paste("AEm = ", AE), ncol = 5)+
  stat_contour(aes(z = ip2, fill = ip2, linetype = factor(..level..)), col = "black", binwidth = 0.25, show.legend = FALSE)+
  geom_hline(yintercept = 0, col = "black")+
  #add points of observed values
  geom_point(aes(shape = as.character(day), col = algae_conc), alpha = 0.7, fill = "black", 
             data = midgeprod %>%  
               left_join(algaeprod %>% 
                           mutate(NPP = gpp_daily*AEa)))+
  geom_errorbar(aes(ymin = Pdc-Pdc.se, ymax = Pdc+Pdc.se, col = algae_conc), width = 0,
                data = midgeprod %>%
                  left_join(algaeprod %>%
                              mutate(NPP = gpp_daily*AEa,
                                     ip2 = NA)))+
  labs(x = expression("Algal Productivity \u03BCg C "~d^{-1}),
       y = expression("Midge Productivity \u03BCg C "~d^{-1}),
       fill = expression(I[m]/P[a]),
       color = "Initial Sediment Quality",
       shape = "Day")+
  scale_color_viridis_c()+
  # theme(legend.position = "right")+
  scale_linetype_manual(values = c("dotted", "dashed", "solid"))


#same plot but using only gpp values from day 14
sim%>% 
  ggplot(aes(x = NPP, y = Pdc, fill = ip2))+
  facet_wrap(~paste("AEm = ", AE), ncol = 5)+
  stat_contour(aes(z = ip2, linetype = factor(..level..)), col = "black", binwidth = 0.25, show.legend = FALSE)+
  geom_hline(yintercept = 0, col = "black")+
  #add points of observed values
  geom_point(aes(shape = as.character(day), col = algae_conc), alpha = 0.7, fill = "black", 
             data = midgeprod %>%  
               left_join(algaeprod %>% filter(day != 22) %>% select(-day) %>% 
                           mutate(NPP = gpp_daily*AEa)))+
  labs(x = expression("Algal Productivity \u03BCg C "~d^{-1}),
       y = expression("Midge Productivity \u03BCg C "~d^{-1}),
       fill = expression(I[m]/P[a]),
       color = "Initial Sediment Quality",
       shape = "Day")+
  scale_color_viridis_c()+
  # theme(legend.position = "right")+
  scale_linetype_manual(values = c("dotted", "dashed", "solid"))




