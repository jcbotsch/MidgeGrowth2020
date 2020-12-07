#====Load Packages====
library(tidyverse)
library(lubridate)

sumna <- function(x){
  if(all(is.na(x))){
    NA
  }
  else{
    sum(x, na.rm = TRUE)
  }
}

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

#====Secondary Production====
#====Routine Data====
# read data
rcc_raw<- read_csv("raw data/routine_cc_9Nov20.csv", col_types = "ccDDccddddd?????_")
rcm_raw <- read_csv("raw data/routine_cm_9Nov20.csv", col_types = "ccDDcccddddc__")
emg_raw <- read_csv("raw data/emergence_clean.csv")


st33emg <- emg_raw %>% filter(site == 3) %>% 
  group_by(site, setdate, coldate) %>% 
  summarise(adults = meanna(smch_tot)/0.17, #0.17 m2 dia trap
            weight = 0.1*adults, #mg (Dreyer 2015)
            weight.se = 0.5*adults) %>% 
  ungroup %>% 
  filter(year(setdate) >2012) %>% 
  select(-site) %>% 
  rename(sampledate = coldate)

#kajak core area
corearea <- 2.5^2*pi/10000 #convert core to m2

#instar data
mg_instars = c(0,4.3,7.3,12, 19.5)/55*1000 #divide by ocular micrometer units to get mm then 1000 to get micrometers

#====Prepare Midge Count Data====
#prep count data
rcc <- rcc_raw %>% 
  filter(month(sampledate)%in%5:8) %>% 
  group_by(sta, coreid, fract_count) %>% 
  gather(species, count, tanyt, ortho, chiro, tanyp, contains("pupae"), tubifex) %>% 
  mutate(count = count/fract_count) %>% 
  group_by(sta, sampledate, coreid, species) %>% 
  summarise(sum = sumna(count)) %>% 
  spread(species, sum) %>% 
  select(sta, sampledate, coreid, tanyt, chiro, ortho, tanyp, everything()) %>% 
  ungroup

#====Preview Midge Measurement Data====
#preview body length
rcm_raw %>% 
  filter(species_name == "tt") %>% 
  group_by(sta, coreid, sampledate) %>% 
  # filter(sta == "33") %>%
  mutate(bl = ifelse(is.na(body_length), 0,1)) %>% 
  filter(!is.na(body_length)) %>% 
  ggplot(aes(body_length/body_length_units_to_mm, fill  = factor(sampledate)))+
  facet_wrap(~year(sampledate))+
  geom_histogram()

#preview head size data
rcm_raw %>% 
  filter(species_name == "tt") %>% 
  ggplot(aes(head_size/head_size_units_to_mm*1000))+
  geom_histogram(bins = 50)

rcm_raw %>% 
  filter(species_name == "tt") %>% 
  ggplot(aes(head_size/head_size_units_to_mm*1000))+
  geom_histogram(bins = 50)+
  geom_vline(xintercept = mg_instars)

#preview head size and body length data
rcm_raw %>% 
  filter(species_name == "tt") %>% 
  ggplot(aes(x = head_size/head_size_units_to_mm*1000, y = body_length/body_length_units_to_mm*50, col = sta))+
  geom_vline(xintercept = mg_instars)+
  geom_point()


#====Prepare Midge Measurement Data 1=====
rcm_2 <- rcm_raw %>% 
  filter(species_name == "tt") %>% 
  mutate(conversion = head_size_units_to_mm, 
         #convert ocular micrometer units to mm then micrometers
         head_size = head_size/conversion*1000, 
         #assign instar based on head capsule widths
         instar = cut(head_size, 
                      breaks = mg_instars, 
                      include.lowest =  FALSE, 
                      labels = c(1,2,3,4)),
         #set instar to NA for non-Tanytarsini
         instar = ifelse(species_name!="tt", NA, instar),
         #convert body lengths to mm
         body_size = body_length/body_length_units_to_mm *1000) %>% 
  select(coreid, sta, sampledate, species_name, head_size, body_size, instar, comments)



#===body length by instars====
rcm_bl <- rcm_2%>% 
  filter(!is.na(body_size)) 

#preview body length by instar
rcm_bl %>% 
  ggplot(aes(x = instar, y = body_size/1000, group = instar))+
  geom_boxplot(outlier.shape = NA)+
  # facet_wrap(~sampledate)+
  geom_jitter(width = 0.4, alpha = 0.2, aes(col = sta))

rcm_bl %>% 
  count(instar)

#is there a difference in size through time?
rcm_bl %>% 
  filter(instar!=1) %>% 
  ggplot(aes(body_size/1000, fill = factor(sampledate)))+
  geom_density(alpha  = 0.5)+
  facet_wrap(~instar)+
  geom_point(aes(col = factor(sampledate)), y = 0, alpha = 0.2)


lengths <- rcm_bl %>% 
  group_by(instar) %>% 
  summarise(mean = mean(body_size),
            sd = sd(body_size),
            q25 = quantile(body_size, 0.25),
            q75 = quantile(body_size, 0.75),
            median = median(body_size),
            samples = n()) %>% 
  filter(!is.na(instar))



weight <- function(l){
  (0.442 + 0.0879 * l )^3 
}

AE = 1

deriv(~log(wt/wtm1)/deltaT * ((Bt+Btm1)/2),
      c("wt", "wtm1", "deltaT", "Bt", "Btm1"))

# expression({
#   .expr1 <- wt/wtm1
#   .expr2 <- log(.expr1)
#   .expr3 <- .expr2/deltaT
#   .expr5 <- (Bt + Btm1)/2
#   .expr22 <- .expr3 * (1/2)
#   .value <- .expr3 * .expr5
#   .grad <- array(0, c(length(.value), 5L), list(NULL, c("wt", 
#                                                         "wtm1", "deltaT", "Bt", "Btm1")))
#   .grad[, "wt"] <- 1/wtm1/.expr1/deltaT * .expr5
#   .grad[, "wtm1"] <- -(wt/wtm1^2/.expr1/deltaT * .expr5)
#   .grad[, "deltaT"] <- -(.expr2/deltaT^2 * .expr5)
#   .grad[, "Bt"] <- .expr22
#   .grad[, "Btm1"] <- .expr22
#   attr(.value, "gradient") <- .grad
#   .value
# })

# OR

deriv(~log(wt/wtm1)/deltaT * ((nt*wt+ntm1*wtm1)/2),
      c("wt", "wtm1", "deltaT", "nt", "ntm1"))

# expression({
#   .expr1 <- wt/wtm1
#   .expr2 <- log(.expr1)
#   .expr3 <- .expr2/deltaT
#   .expr7 <- (nt * wt + ntm1 * wtm1)/2
#   .value <- .expr3 * .expr7
#   .grad <- array(0, c(length(.value), 5L), list(NULL, c("wt", 
#                                                         "wtm1", "deltaT", "nt", "ntm1")))
#   .grad[, "wt"] <- 1/wtm1/.expr1/deltaT * .expr7 + .expr3 * 
#     (nt/2)
#   .grad[, "wtm1"] <- .expr3 * (ntm1/2) - wt/wtm1^2/.expr1/deltaT * 
#     .expr7
#   .grad[, "deltaT"] <- -(.expr2/deltaT^2 * .expr7)
#   .grad[, "nt"] <- .expr3 * (wt/2)
#   .grad[, "ntm1"] <- .expr3 * (wtm1/2)
#   attr(.value, "gradient") <- .grad
#   .value
# })



#Delta method to propogate error:
#Sy = sqrt((dy/dx)^2*s^2x + (dz/dx)^2*s^2z + .... + cov)



rmp <- rcm_2 %>% 
  filter(month(sampledate) %in% c(5:8)) %>% 
  left_join(lengths) %>% 
  group_by(sta, sampledate, species_name, coreid) %>% 
  summarise(l = mean(median)/1000) %>% #average estimated length of a midge in this sample
  mutate(w = weight(l)) %>% 
  group_by(sta, sampledate) %>% 
  summarise(n = n(),
            wse = sd(w, na.rm = TRUE)/sqrt(n), 
            w = meanna(w)) %>% 
  full_join(rcc %>% 
              select(sta, sampledate, coreid, tanyt) %>% 
              mutate(nt = tanyt/corearea) %>% 
              group_by(sta, sampledate) %>% 
              summarise(n = n(),
                        ntse = sd(nt)/sqrt(n),
                        nt = meanna(nt))) %>% 
  group_by(sta) %>% 
  mutate(nt = nt+1,
         w = ifelse(nt == 1, 0, w), #mg
         ntm1 = lag(nt),
         ntm1se = lag(ntse),
         wt = w+0.03, #add 10% of 1 1st instar weight to prevent 0s
         wtm1 = lag(wt),
         deltaT = as.numeric(sampledate - lag(sampledate)),
         Bt = nt*wt,
         Bt.se = sqrt(wt^2*wse^2 + nt^2*ntse^2),
         Btm1 = lag(Bt),
         Btm1.se = lag(Bt.se),
         barB = (Bt + Btm1)/2) %>% 
  mutate(g= log(wt/wtm1)/deltaT, #daily instantaneous growth
         Pd = g*barB, #Daily midge productivity mg AFDM m-2 d-1
         PdE1 = wt/wtm1,
         PdE2 = log(PdE1),
         PdE3 = PdE2/deltaT,
         PdE7 = (nt*wt+ntm1*wtm1)/2,
         Pdwt = 1/wtm1/PdE1/deltaT*PdE7+PdE3*(nt/2), #partial derivative of Pd with respect to wt
         Pdwtm1 = PdE3 * (ntm1/2) - wt/wtm1^2/PdE1/deltaT * PdE7, #partial derivative of Pd with respect to wt-1
         Pdnt = PdE3 * (wt/2), #partial derivative of Pd with respect to nt
         Pdntm1 = PdE3 * (wtm1/2), #partial derivative of Pd with respect to nt-1,
         Pd.se = sqrt(Pdwt^2*wse^2 + Pdwtm1^2*wse^2 + Pdnt^2*ntse^2 + Pdntm1^2*ntse^2), #estimated std error on Pd
         Pdc.se = (Pd.se*0.5)/1000,
         Pdc = (Pd*0.5)/1000) %>%  #Daily midge productivity g C m-2 d-1) %>% 
  ungroup


rmp %>% 
  ggplot(aes(x = sampledate, y = Bt))+
  geom_point()+
  geom_errorbar(aes(ymin = Bt-Bt.se, ymax = Bt+Bt.se))+
  facet_wrap(~year(sampledate), scales = "free")+
  geom_line()

#spot checks

skimr::skim(rmp)
skimr::skim(rcm_2)


rmp %>% filter(is.infinite(g))
rmp %>% filter(is.na(wse)) %>% left_join(rcm_2) %>% 
  count(sta, sampledate, coreid) %>% as.data.frame()

rmp %>% filter(is.na(wse)) %>% left_join(rcm_2) %>% 
  count(sta, sampledate, coreid) %>% as.data.frame()

sample_n(rmp, 10) 


rcc %>% filter(sta == "33") %>% left_join()





rmp%>% 
  filter(sta == "33", year(sampledate)>2012) %>% 
  filter(!is.na(Pd.se)) %>% 
  ggplot(aes(x = sampledate, y = Pdc, fill = log10(nt+1), group = sta))+
  facet_wrap(~year(sampledate), scales = "free_x")+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_line()+
  geom_point(shape = 21, size = 2)+
  geom_errorbar(aes(ymin = Pdc-Pdc.se, ymax = Pdc.se+Pdc), width = 0)+
  labs(x = "Date",
       y = "Tanytarsini Productivity in g C m-2 d-1")+
  scale_fill_viridis_c(option = "plasma")

rmp %>% 
  crossing(AE = c(0.1, 0.5, 1)) %>% 
  mutate( NPE = 0.54, #net production efficiency (Lindegaard)
          AE = AE,  
          ingestion  = Pdc/(NPE*AE)) %>% 
  filter(AE ==0.1) %>% 
  filter(sta == "33", year(sampledate)!=2012) %>% 
  ggplot(aes(x = yday(sampledate), y = ingestion, fill = log10(nt+1), group = sta))+
  facet_wrap(~year(sampledate), scales = "free_y")+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_line()+
  geom_point(shape = 21, size = 2)+
  labs(x = "Date",
       y = "Daily Ingestion in g C m-2 d-1",
       title = "AE = 0.1")+
  scale_fill_viridis_c(option = "plasma")

rmp %>% 
  crossing(AE = c(0.1, 0.5, 1)) %>% 
  mutate( NPE = 0.54, #net production efficiency (Lindegaard)
          AE = AE,  
          ingestion  = Pdc/(NPE*AE)) %>% 
  filter(AE ==0.5) %>% 
  filter(sta == "33", year(sampledate)!=2012) %>% 
  ggplot(aes(x = yday(sampledate), y = ingestion, fill = log10(nt+1), group = sta))+
  facet_wrap(~year(sampledate), scales = "free_y")+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_line()+
  geom_point(shape = 21, size = 2)+
  labs(x = "Date",
       y = "Daily Ingestion in g C m-2 d-1",
       title = "AE = 0.5")+
  scale_fill_viridis_c(option = "plasma")

rmp %>% 
  crossing(AE = c(0.1, 0.5, 1)) %>% 
  mutate( NPE = 0.54, #net production efficiency (Lindegaard)
          AE = AE,  
          ingestion  = Pdc/(NPE*AE)) %>% 
  filter(AE ==1) %>% 
  filter(sta == "33", year(sampledate)!=2012) %>% 
  ggplot(aes(x = yday(sampledate), y = ingestion, fill = log10(nt+1), group = sta))+
  facet_wrap(~year(sampledate), scales = "free_y")+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_line()+
  geom_point(shape = 21, size = 2)+
  labs(x = "Date",
       y = "Daily Ingestion in g C m-2 d-1",
       title = "AE = 1")+
  scale_fill_viridis_c(option = "plasma")






