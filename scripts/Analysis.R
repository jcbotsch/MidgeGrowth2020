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

g14 <- lm(gpp~algae_conc*midge+box, data = nep14) 
summary(g14)
plot(g14)

nep22 <- nep %>% 
  filter(day == 22)

g22 <- lm(gpp~algae_conc*midge+box, data = nep22)
summary(g22)
plot(g22)

#Figure 2
nep %>% 
  mutate(day = paste("Day", day)) %>% 
  ggplot(aes(x = algae_conc, y = gpp, fill = midge))+
  facet_wrap(~day)+
  geom_point(size  = 2, alpha = 0.7, shape = 21)+
  geom_smooth(aes(color = midge), method = "lm", size = 0.75, se = FALSE)+
  scale_fill_manual(values = c("black", "gray60"))+
  scale_color_manual(values = c("black", "gray60"))+
  labs(x = "Initial Sediment Quality",
       y = expression(GPP(g~O[2]~m^{-2}~hr^{-1})),
       color = "",
       fill = "")



#Supplemental Figure
nep %>% 
  gather(metabolism, value, resp, nep, gpp) %>% 
  mutate(day = day+rnorm(n(), 0.001), #add jitter
         metabolism = toupper(ifelse(metabolism == "resp", "er", metabolism)),
         metabolism = fct_reorder(metabolism, desc(value))) %>%  
  ggplot(aes(x = day, y = value, col = algae_conc))+
  facet_grid(metabolism~midge, scales = "free")+
  geom_hline(yintercept = 0)+
  geom_point()+
  geom_line(aes( group = coreid), alpha = 0.4)+
  geom_smooth(method = "lm", col = "black", se = FALSE)+
  theme(legend.position = "bottom")+
  labs(y = expression("Change in" ~O[2] (g~m^{-2}~hr^{-1})),
       x = "Incubation Day",
       color = "Initial Sediment Quality")+
  scale_color_viridis_c()+
  scale_x_continuous(breaks = c(14, 22))


##Midge Effects on chlorophyll

#This is especially sloppy because I might remove all of this. 
lm(NDVI_R~algae_conc*midge+box, data = ndvi %>% left_join(meta) %>% filter(day ==14)) %>% summary()
lm(NDVI_R~algae_conc*midge+box, data = ndvi %>% left_join(meta) %>% filter(day ==22)) %>% summary()

ndvi %>% 
  left_join(meta) %>% 
  ggplot(aes(x = algae_conc2, y = NDVI_R, fill = midge))+
  geom_point(size = 2, shape = 21)+
  facet_wrap(~day)+
  geom_smooth(aes(col = midge), method = "lm", se = FALSE)+
  scale_fill_manual(values = c("black", "gray60"))+
  scale_color_manual(values = c("black", "gray60"))+
  labs(x = "Initial Sediment Quality",
       y = "NDVI")

#supplement
ndvi %>%
  mutate(coreid = as.character(coreid)) %>% 
  left_join(nep) %>% 
  ggplot(aes(x = gpp, y = NDVI_R, col = algae_conc))+
  geom_point()+
  facet_grid(midge~day)+
  geom_smooth(method = "lm", se = FALSE, col = "black", size = 0.7)+
  scale_color_viridis_c()+
  labs(x = expression(GPP~(g^{2}~O[2]~m^{-2}~hr^{-1})),
       y = "NDVI",
       color = "Sediment\nTreatment")

ndvi %>% 
  left_join(meta) %>% 
  ggplot(aes(x = factor(day), y = NDVI_R, col = algae_conc, group = factor(algae_conc)))+
  geom_jitter(width = 0.2)+
  facet_wrap(~midge)+
  geom_smooth(method = "lm", se = FALSE, size = 0.7)+
  scale_color_viridis_c(trans = "log1p")+
  labs(x = "Day",
       y = "NDVI",
       color = "Sediment\nTreatment")



#====Question 2: Sediment Effects on Midges =====
## Number of Midges Present
cc <- cc_raw %>% 
  left_join(meta) %>% 
  mutate(coreid = as.character(coreid),
         box = as.character(box))

cc14 <- cc %>% 
  filter(day == 14)

nm14 <- glm(live_tt~algae_conc*midge+box, 
    data = cc14,
    family = quasipoisson())

summary(nm14)
plot(nm14)

cc22 <- cc %>% 
  filter(day == 22)

nm22 <- glm(live_tt~algae_conc*midge+box, 
    data = cc22,
    family = quasipoisson())

summary(nm22)
plot(nm22)

#Figure 3
cc %>% 
  left_join(nep) %>% 
  filter(day %in% c(14, 22)) %>% 
  mutate(day = paste("Day", day)) %>% 
  ggplot(aes(x = algae_conc, y = live_tt, fill = gpp))+
  facet_grid(midge~day, scales = "free")+
  geom_jitter(shape = 21, size = 2, alpha = 0.6, width = 0.0005, height = 0.3)+
  scale_fill_viridis_c()+
  labs(y = "Live Tanytarsini",
       x = "Initial Sediment Quality",
       fill = expression(GPP(g~O[2]~m^{-2}~hr^{-1})~" "))

##Midge Development
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

#figure 3
p3 %>% 
  filter(!is.na(algae_conc)) %>% 
  mutate(day = paste("Day", day)) %>% 
  ggplot(aes(y = prop, x = algae_conc, fill = live_tt))+
  geom_point(shape = 21, size = 2, alpha = 0.7)+
  geom_hline(yintercept = startingprop)+
  facet_grid(midge~day)+
  scale_fill_viridis_c(option = "plasma")+
  lims(y = c(0,1))+
  labs(x = expression(GPP~(g^{2}~O[2]~m^{-2}~hr^{-1})),
       y = "Proportion of Larval Midges\n3rd instar or higher",
       fill = "Number of Larvae")+
  theme(legend.position = "bottom")


p314 <- p3 %>% 
  filter(day == 14)

prop14 <- glm(prop~algae_conc*midge+box, data = p314, family = quasibinomial())

summary(prop14)
plot(prop14)


p322 <- p3 %>% 
  filter(day == 22)

prop22 <- glm(prop~algae_conc*midge+box, data = p322, family = quasibinomial())

summary(prop22)
plot(prop22)


##Length of Midges
l3 <- cm %>% 
  filter(instar == 3,
         day!=0)

l314 <- l3 %>% 
  filter(day == 14)

bl14 <- lmer(body_size~algae_conc*midge+box+(1|coreid), data = l314) 

summary(bl14)
plot(bl14)

Anova(bl14, type = "3", test.statistic = "F")
Anova(bl14, type = "2", test.statistic = "F")



l322 <- l3 %>% 
  filter(day == 22)

bl22 <- lmer(body_size~algae_conc*midge+box +(1|coreid), data = l322)
summary(bl22)
plot(bl22)

Anova(bl22, type = "3", test.statistic = "F")
Anova(bl22, type = "2", test.statistic = "F")


#Figure 4
l3 %>% 
  left_join(nep) %>%
  mutate(day = paste("Day", day)) %>% 
  ggplot(aes(x = algae_conc, y = body_size, fill = gpp))+
  facet_grid(midge~day, scales = "free_x")+
  geom_point(alpha = 0.2, shape = 21, size = 2)+
  geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.7)+
  scale_fill_viridis_c()+
  labs(x = "Initial Sediment Conditions",
       y = "Body Length  (mm)",
       fill = expression(GPP~(g^{2}~O[2]~m^{-2}~hr^{-1})~" "))+
  theme(legend.position = "bottom")


##SUPPLEMENT: Did midges increase in body length?
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




#====Question 3: How Much do midges rely on algae?====
#====Q3 Version 1====
#convert midges to biomass following Lindegaard and Jonasson 1979
#AFDW in mg^(1/3) = 0.0442 + 0.0879 * length in mm 
weight <- function(l){
  (0.442 + 0.0879 * l )^3 
}

#calculate carbon using Amanda's Isotope Data
dwc = 0.5


#increment Summation method to estimate growth 
incsumprod <- function(n1, n2, wt1, wt2, deltaT){
  g = wt2-wt1 #individual growth
  Nbar = (n1 + n2)/2 #Average abundance
  intp = Nbar*g #production over the interval
  Pd = intp/deltaT #daily production
  list(g = g, Nbar = Nbar, intp = intp, Pd = Pd)
}

#calculate midge productivity
midgeprod <- cm %>% 
  left_join(cc) %>% 
  filter(species_name == "tt",
         day!=0,
         midge == "Midges") %>% 
  group_by(coreid) %>% 
  mutate(w = weight(body_size)) %>% #estimate dry mass
  group_by(coreid, day, box, algae_conc) %>% 
  summarise(wt = meanna(w), #average biomass of a midge (mg)
            nt = unique(live_tt)) %>% #number of individuals in the mesocosm
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
              summarise(wt1 = mean(w), #weight at t = 0
                        n1 = mean(live_tt)) %>% #number at t = 0
              select(-sampledate)) %>% 
  ungroup %>% 
  mutate(Pd = incsumprod(n1 = n1, n2 = nt, wt1 = wt1, wt2 = wt, deltaT = day)$Pd, #production in mg,
         g = incsumprod(n1 = n1, n2 = nt, wt1 = wt1, wt2 = wt, deltaT = day)$g,
         Pdc = Pd*dwc*1000) #convert from weight to c (mg) then convert to micrograms 

#compare measurements of growth (lm) to these estimates
midgeprod %>% 
  group_by(day) %>% 
  mutate(meanw = mean(wt-wt1)) %>% 
  ggplot(aes(wt-wt1, fill = factor(day)))+
  geom_histogram(bins = 15, alpha = 0.7)+
  geom_vline(xintercept = 0, linetype = "dotted")+
  geom_vline(xintercept = coef(growth)[2], linetype = "dashed")+
  geom_vline(xintercept = weight(coef(growth)[2]))+
  geom_vline(aes(xintercept = meanw, col = factor(day)))

##Estimate daily production of primary producers in C
#converting g O2--> C
#1:1 molar mass
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


#Set values for transfer efficiencies
AEa = 1 #arbitrary value set for the assimilation of fixed carbon by algae
NPEm = 0.54 #Net Production efficiency for T. gracilentus at Myvatn Lindegaard 1994

#simulate parameter space following migge and algal productivity values
sim <- data.frame(Pdc = 0:max(midgeprod$Pdc)) %>% 
  crossing(GPP = min(algaeprod$gpp_daily):max(algaeprod$gpp_daily)) %>% 
  crossing(AE = seq(0.1, 1, by = 0.1)) %>% 
  #I = Assimilation attributed to algae/GPE
  mutate(ingestion = Pdc/(NPEm*AE), #this assumes 100% of production attributed to algae
         NPP = GPP*AEa, #converting C fixation to biomass accumulation
         ip = ingestion/NPP) %>%  #ingestion by midges/production of resource
  mutate(ip2 = ifelse(ip>1, NA, ip)) 

#exploring midge productivity and algal productivity across the parameter space
sim%>% 
  ggplot(aes(x = NPP, y = Pdc, fill = ip2))+
  facet_wrap(~paste("AEm = ", AE))+
  geom_tile()+
  #add points of observed values
  geom_point(aes(shape = as.character(day), col = algae_conc), alpha = 0.7, fill = "black", 
             data = midgeprod %>%  
               left_join(algaeprod %>% 
                           mutate(NPP = gpp_daily*AEa)))+
  labs(x = expression("Algal Productivity \u03BCg C "~d^{-1}),
       y = expression("Midge Productivity \u03BCg C "~d^{-1}),
       fill = expression(I[m]/P[a]),
       color = "Initial Sediment Quality",
       shape = "Day")+
  scale_fill_gradient2(low = "darkgreen", high = "goldenrod", midpoint = 0.5)+
  scale_color_viridis_c()+
  theme(legend.position = "right")

#how many mesocosms could survive on algae
midgeprod %>%  
  filter(Pdc>0) %>% 
  left_join(algaeprod %>% 
              mutate(NPP = gpp_daily*AEa)) %>% 
  select(coreid, day, box, algae_conc, Pdc, gpp_daily, NPP) %>% 
  crossing(AE = seq(0.1,1,by = 0.1)) %>% 
  mutate(ingestion = Pdc/(NPEm*AE),
         ip = ingestion/NPP) %>% 
  group_by(AE) %>% 
  count(algae_only = ip<1) %>% 
  ggplot(aes(fill = algae_only, y = n, x= factor(AE)))+
  geom_col()+
  labs(fill = "Possible for Midges to Survive only on Algae",
       y = "Number of Mesocosms",
       x = "AEm")+
  scale_fill_manual(values = c("firebrick3", "dodgerblue"))


#What AE could midges survive on algae?
midgeprod %>%  
  filter(Pdc>0) %>% 
  left_join(algaeprod %>% 
              mutate(NPP = gpp_daily*AEa)) %>% 
  select(coreid, day, box, algae_conc, Pdc, gpp_daily, NPP) %>% 
  crossing(AE = seq(0.1,1,by = 0.1)) %>% 
  mutate(ingestion = Pdc/(NPEm*AE),
         ip = ingestion/NPP) %>%
  ggplot(aes(x = AE, y = ip, col = NPP))+
  facet_wrap(~round(algae_conc, 3), nrow = 2)+
  lims(y = c(0,1),
       x = c(0,1))+
  geom_point()+
  geom_line(aes(group = coreid))+
  geom_text(aes(label = n), color = "black",
             data = midgeprod %>%  
               filter(Pdc>0) %>% 
               left_join(algaeprod %>% 
                           mutate(NPP = gpp_daily*AEa)) %>% 
               select(coreid, day, box, algae_conc, Pdc, gpp_daily, NPP) %>% 
               crossing(AE = seq(0.1,1,by = 0.1)) %>% 
               mutate(ingestion = Pdc/(NPEm*AE),
                      ip = ingestion/NPP) %>%
               filter(AE==1, ip>1) %>% 
               group_by(algae_conc) %>% 
               count() %>% 
               mutate(AE = 0.1,
                      ip = 0.9))+
  scale_color_viridis_c()+
  coord_equal()+
  labs(x = expression(AE[m]),
       y = expression(I[m]/P[a]))+
  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5))

#how many mesocosms are not in the plot?
midgeprod %>%  
  filter(Pdc>0) %>% 
  left_join(algaeprod %>% 
              mutate(NPP = gpp_daily*AEa)) %>% 
  select(coreid, day, box, algae_conc, Pdc, gpp_daily, NPP) %>% 
  crossing(AE = seq(0.1,1,by = 0.1)) %>% 
  mutate(ingestion = Pdc/(NPEm*AE),
         ip = ingestion/NPP) %>%
  filter(AE==1, ip>1) %>% 
  group_by(algae_conc) %>% 
  count()
  
#====Q3 Version 2====

algaeprod %>% 
  select(coreid, day, gpp_daily) %>% 
  mutate(day = paste0("day",day)) %>% 
  spread(day, gpp_daily) %>% 
  mutate(avg.gpp = (day22+day14)/2,
         change = (day22-day14)/avg.gpp) %>% 
  left_join(midgeprod %>% select(coreid, Pdc)) %>% 
  filter(!is.na(Pdc)) %>% 
  ggplot(aes(x = change, y = Pdc, col = avg.gpp))+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  geom_point()+
  scale_color_viridis_c()


algaeprod %>% 
  select(coreid, day, gpp_daily) %>% 
  mutate(day = paste0("day",day)) %>% 
  spread(day, gpp_daily) %>% 
  mutate(avg.gpp = (day22+day14)/2,
         change = (day22-day14)/avg.gpp) %>% 
  left_join(midgeprod %>% select(coreid, Pdc)) %>% 
  filter(!is.na(Pdc), !is.na(change)) %>% 
  lm(Pdc~change, data = .) %>% 
  summary()

#=====Q3 Version 3====

gppchange <- lm(gpp*18~midge*day, data = nep) #daily GPP

growth2 <- lm(body_size~day*midge, data = cm %>% filter(species_name == "tt"))


nd <- crossing(midge = c("Midges", "No Midges"), day = 1:22)
nd$gpp <- predict(gppchange, newdata = nd)
nd$gppse <- predict(gppchange, newdata = nd, se.fit = TRUE)$se.fit

nd$length <- predict(growth2, newdata = nd)-coef(growth2)[1]
nd$lengthse <- predict(growth2, newdata = nd, se.fit = TRUE)$se.fit


nd %>% 
  gather(var, val, gpp, length) %>%
  gather(term, se, contains("se")) %>% 
  filter(term == paste0(var, "se")) %>% 
  ggplot(aes(x = day, y = val, col = var, fill = var))+
  facet_wrap(~midge)+
  geom_ribbon(aes(ymin = val-se, ymax = val+se), alpha = 0.5)+
  geom_line()+
  labs(x = "Day",
       y = "Change",
       color = "",
       fill = "")

#same exercise but with C

nppchange <- lm(npp~midge*day, data = algaeprod %>% mutate(npp = gpp_daily*AEa)) 

growthw2 <- lm(Pdc~day, data = midgeprod)


nd <- crossing(midge = c("Midges", "No Midges"), day = 1:22)
nd$npp <- predict(nppchange, newdata = nd)
nd$nppse <- predict(nppchange, newdata = nd, se.fit = TRUE)$se.fit

nd$midgegrowth <- predict(growthw2, newdata = nd)
nd$midgegrowthse <- predict(growthw2, newdata = nd, se.fit = TRUE)$se.fit

coef(nppchange)
coef(growthw2)

nd %>% 
  gather(var, val, npp, midgegrowth) %>%
  gather(term, se, contains("se")) %>% 
  filter(term == paste0(var, "se")) %>% 
  mutate(var = ifelse(var == "npp", "Algae", "Midge"),
         val = ifelse(var == "Midge"& midge == "No Midges", NA, val)) %>% 
  ggplot(aes(x = day, y = val, col = var, fill = var))+
  facet_wrap(~midge)+
  geom_hline(yintercept = 0)+
  geom_ribbon(aes(ymin = val-se, ymax = val+se), alpha = 0.5)+
  geom_line()+
  labs(x = "Day",
       y = expression("Production \u03BCg C"~d^{-1}),
       col = "",
       fill = "")

#same as the last one but by algal concentration
nppchange2 <- lm(npp~midge*day*algae_conc, data = algaeprod %>% mutate(npp = gpp_daily*AEa)) 

growthw3 <- lm(Pdc~day*algae_conc, data = midgeprod)


nd <- crossing(midge = c("Midges", "No Midges"), day = 1:22, algae_conc = 1/10^(0:5))
nd$npp <- predict(nppchange2, newdata = nd)
nd$nppse <- predict(nppchange2, newdata = nd, se.fit = TRUE)$se.fit

nd$midgegrowth <- predict(growthw3, newdata = nd)
nd$midgegrowthse <- predict(growthw3, newdata = nd, se.fit = TRUE)$se.fit

coef(nppchange)
coef(growthw2)

nd %>% 
  gather(var, val, npp, midgegrowth) %>%
  gather(term, se, contains("se")) %>% 
  filter(term == paste0(var, "se")) %>% 
  mutate(var = ifelse(var == "npp", "Algae", "Midge"),
         val = ifelse(var == "Midge"& midge == "No Midges", NA, val)) %>% 
  filter(midge == "Midges") %>% 
  ggplot(aes(x = day, y = val, col = var, fill = var))+
  facet_wrap(~algae_conc)+
  geom_hline(yintercept = 0)+
  geom_ribbon(aes(ymin = val-se, ymax = val+se), alpha = 0.5)+
  geom_line()+
  labs(x = "Day",
       y = expression("Production \u03BCg C"~d^{-1}),
       col = "",
       fill = "")



#########################
#====All Repeated using logged values=====

#Repeat using Logged algae concentration

#Q1 LOG

g14log <- lm(gpp~log(algae_conc2)*midge+box, data = nep14) 
summary(g14log)
plot(g14log)

g22log <- lm(gpp~log(algae_conc2)*midge+box, data = nep22)
summary(g22log)
plot(g22log)

#Figure 2 (LOG)
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


##Q2 LOG

## Number of Midges Present
nm14log <- glm(live_tt~log(algae_conc2)*midge+box, 
            data = cc14,
            family = quasipoisson())

summary(nm14log)
plot(nm14log)

nm22log <- glm(live_tt~log(algae_conc2)*midge+box, 
            data = cc22,
            family = quasipoisson())

summary(nm22log)
plot(nm22log)

#Figure 3
cc %>% 
  left_join(nep) %>% 
  filter(day %in% c(14, 22)) %>% 
  mutate(day = paste("Day", day)) %>% 
  ggplot(aes(x = algae_conc2, y = live_tt, fill = gpp))+
  facet_grid(midge~day, scales = "free")+
  geom_jitter(shape = 21, size = 2, alpha = 0.6, width = 0.0005, height = 0.3)+
  scale_fill_viridis_c()+
  labs(y = "Live Tanytarsini",
       x = "Initial Sediment Quality",
       fill = expression(GPP(g~O[2]~m^{-2}~hr^{-1})~" "))+
  scale_x_log10()

#figure 4
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



prop14log <- glm(prop~log(algae_conc2)*midge+box, data = p314, family = quasibinomial())

summary(prop14log)
plot(prop14log)


prop22log <- glm(prop~log(algae_conc2)*midge+box, data = p322, family = quasibinomial())

summary(prop22log)
plot(prop22log)


##Length of Midges
bl14log <- lmer(body_size~log(algae_conc2)*midge+box+(1|coreid), data = l314) 

summary(bl14log)
plot(bl14log)

Anova(bl14log, type = "3", test.statistic = "F")
Anova(bl14log, type = "2", test.statistic = "F")


bl22log <- lmer(body_size~log(algae_conc2)*midge+box +(1|coreid), data = l322)
summary(bl22log)
plot(bl22log)

Anova(bl22log, type = "3", test.statistic = "F")
Anova(bl22log, type = "2", test.statistic = "F")


#Figure 4
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





