#====Load Packages====
library(tidyverse)
library(lubridate)

#====Metadata====
algae <- read_csv("raw data/Algal_Treatments29Oct20.csv")
tubes <- read_csv("raw data/Tube_Treatments29Oct20.csv")

meta <- left_join(tubes, algae)

#====Preliminary Measurements on Sediment Treatments====
## chlorophyll
chl_raw <- read_csv("raw data/chl29Oct20.csv")

chl <- chl_raw %>% 
  #first I need to calculate the uncorrected chlorophyll
  mutate(uncorr_chl = (vol_sample+vol_met)/vol_sample * (vol_met_dilu+extracted_vol)/extracted_vol*reading_pre,
         #then I need to calculate the phaeophytin
         pheo = (vol_sample+vol_met)/vol_sample * (vol_met_dilu+extracted_vol)/ extracted_vol * reading_post *(1 + 0.03), #phaeophytin (micrograms/L)
         #finally calculate corrected chlorophyll
         chl = uncorr_chl-pheo, #micrograms/L
         day = 0) %>%  
  select(algae, rep, sampledate, day, uncorr_chl, pheo, chl)


## Organic Content
loi_raw <- read_csv("raw data/LOI29OCt20.csv")

#check comments
loi_raw %>% 
  filter(!is.na(comments))

loi <- loi_raw %>% 
  select(-comments) %>% 
  mutate(ww = ww_tray-trayweight,
         dw = dw_tray-trayweight,
         comb = comb_tray-trayweight,
         perc.org = (dw-comb)/dw,
         day = 0) %>% 
  select(sampledate, day, algae, ww, dw, comb, perc.org)

##NDVI 
#these may want to be altered in the future
NDVI1 <- read_csv("raw data/MidgeGrowth_NDVI_R_20Aug20.csv")
NDVI2 <- read_csv("raw data/MidgeGrowth_NDVI_R_29Aug20.csv")

ndvi <- NDVI1 %>% 
  mutate(day = 14) %>% 
  select(coreid, day, NDVI_R) %>% 
  bind_rows(NDVI2 %>% 
              mutate(day = 22) %>% 
              select(coreid, day, NDVI_R))

#====Metabolism Measurements====
nep_raw <- read_csv("raw data/DO29Oct20.csv")

#check flags
nep_raw %>% 
  filter(flag!=0)

#check comments
nep_raw %>% 
  filter(!is.na(comments)) %>% 
  select(coreid, dark_light, date_start, comments)


#account for dilutions by adding water
wateradditions <- nep_raw %>% 
  filter(coreid == 0) %>%
  group_by(dark_light, date_start) %>% 
  filter(time_start == max(time_start)) %>% 
  select(dark_light, date_start, contains("init"), -contains("probe"), -contains("water"), -contains("barom"))

colnames(wateradditions)[3:5] <- paste("added", str_remove(colnames(wateradditions)[3:5], "_init"), sep = "_")

#water column height is 6.8 cm
height = 6.8 #cm
height_m = height/100 #m
nep1 <- nep_raw %>% 
  select(-flag, -comments) %>% 
  filter(coreid!=0) %>% #remove pitcher readings
  left_join(wateradditions) %>% #add dilutions added
  mutate(dttm_start = as_datetime(paste(date_start, #convert to date time. 
                                        paste0(time_start, ":00"))), #need to add seconds to times
         dttm_end = as_datetime(paste(date_end, 
                                      paste0(time_end, ":00"))),
         delta_time = as.numeric(dttm_end-dttm_start), #change in time (hours)
         proportion = (column_vol - water_added_init)/column_vol, #proportion of initial DO from tube 
         wtemp_init_cor = wtemp_init*proportion+added_wtemp*(1-proportion), #correct water temp for adding water
         do_init_cor = do_init*proportion+added_do*(1-proportion), #correct do for adding water
         delta_do = do_final - do_init_cor, #change in do (mg/L)
         delta_do_rate = delta_do/delta_time, #change in do (mg/L)/hour
         #multiply by 1000 to get m3 divide by 1000 to get g
         delta_do_rate = delta_do_rate*height_m, #multiply by height in m to get g (O2/m2)/h
         date = if_else(date_start == "2020-08-19", as.Date("2020-08-18"),
                        if_else(date_start == "2020-08-28", as.Date("2020-08-27"), as.Date(date_start))))  

#final cleaning step
nep <- nep1 %>%
  select(coreid, date, dark_light, delta_do_rate) %>% 
  spread(dark_light, delta_do_rate) %>% 
  rename(nep = light, resp = dark, startdate = date) %>% 
  mutate(gpp = nep  - resp,
         day = ifelse(startdate == "2020-08-18", 14, 22)) %>% 
  arrange(startdate, coreid) %>% 
  select(coreid, startdate, day, resp, nep, gpp)


#====Midge Measurements====
#read data
cm_raw <- read_csv("raw data/midge_measures29Oct20.csv")

#check flags and comments
cm_raw %>% 
  filter(!is.na(comments)) %>% 
  select(coreid, sampledate, comments) %>% 
  unique()

cm_raw %>% 
  filter(flag!=0) %>% 
  select(coreid, sampledate, flag, comments)

#find instars
mg_instars = c(0,4.3,7.3,12, 15)/55*1000 #divide by ocular micrometer units to get mm then 1000 to get micrometers
#Instars from Lindegaard: 75 125-140 225-250 325-350 (I'm skeptical of the veracity of the 4th instar measures)

#check instars
cm_raw %>% 
  mutate(head_size = ifelse(scope == "Wild", head_size/55*1000, head_size/51*1000)) %>% 
  filter(species_name == "tt") %>% 
  ggplot(aes(head_size, fill = sampledate))+
  geom_histogram(bins = 40)+
  geom_vline(xintercept = mg_instars)+
  labs(x = "Head Size \u03BCm",
       fill = "Date",
       y = "Count")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_blank())


cm <- cm_raw %>% 
        #account for differences in microscopes used
  mutate(conversion = ifelse(scope == "Wild", 55, 
                             ifelse(scope == "Zeiss",51,NA)), 
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
         body_size_mag = body_size_mag/50, #account for differing magnifications during measurments
         body_size = (body_size/body_size_mag)/conversion *1000,
         day = if_else(sampledate == "2020-08-21", 14,
                        if_else(sampledate == "2020-08-07", 0, 22))) %>% 
  select(coreid, sampledate, day, species_name, head_size, body_size, instar, comments)


#====Midge Counts====
midge_raw <- read_csv("raw data/midge_counts29Oct20.csv")

#check comments
midge_raw %>% 
  filter(!is.na(comments))

#check flags
midge_raw %>% 
  filter(flag!=0)

cc <- midge_raw %>% 
  select(-comments) %>% 
  mutate(day = ifelse(sampledate == "2020-08-21", 14, 22)) %>% 
  select(coreid, sampledate, day, live_tt, live, dead)

#add counts for stocked midges
cc <- cc %>% 
  bind_rows(cm %>% 
              filter(day == 0) %>% 
              group_by(coreid) %>% 
              add_count() %>% 
              rename(live_tt = n) %>% 
              mutate(live = live_tt,
                     dead = 0) %>% 
              select(coreid, sampledate, day, live_tt, live, dead) %>% 
              unique())

#=====Write to CSV====


# write_csv(meta, "clean data/MG_meta.csv")
# write_csv(chl, "clean data/MG_chl.csv")
# write_csv(loi, "clean data/MG_om.csv")
# write_csv(ndvi, "clean data/MG_ndvi.csv")
# write_csv(nep, "clean data/MG_ep.csv")
# write_csv(cm, "clean data/MG_cm.csv")
# write_csv(cc, "clean data/MG_cc.csv")
