#====Load Packages====
library(tidyverse)
library(lubridate)
source("MG_Functions.R")

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
  mutate(ww = ww_tray-trayweight, #account for weight of aluminum tray sample weighed in
         dw = dw_tray-trayweight,
         comb = comb_tray-trayweight,
         perc.org = (dw-comb)/dw,
         day = 0) %>% 
  select(sampledate, day, algae, ww, dw, comb, perc.org)

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

#extract timing for Fig S3
inc_timing <- nep1 %>% 
  group_by(dark_light, date) %>% 
  summarise(dttm_start = min(dttm_start),
            dttm_end = max(dttm_end))


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

#final cm edits
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
  filter(head_size>70,#first instar head size is 75. there is one individual with a lower head size that was damaged.
         species_name == "tt") %>%  # there are two non-Tanytarsus
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


#====Hobo Data====
hobo_raw <- read_csv("raw data/MGhobo.csv", skip = 1) %>% 
  janitor::clean_names() 

#fix column names and trim to experiment dates
hobo <- hobo_raw %>% 
  rename(date_time = 2,
         temp = 3,
         lux = 4) %>% 
  select(date_time, temp, lux) %>% 
  mutate(date_time = mdy_hms(date_time)) %>% 
  filter(date_time>as_datetime("2020-08-07 15:00:00"),
         date_time<as_datetime("2020-08-29 10:50:00"))

#mean and SD temperature
hobo %>% 
  filter(date_time>as.Date("2020-08-10")) %>% 
  summarise(tempmean = mean(temp),
            tempsd = sd(temp))


hobo %>% 
  filter(date_time<as_datetime("2020-08-16 23:00:00"),
         lux!=0) %>% 
  summarise(tempmean = mean(temp),
            tempsd = sd(temp),
            luxmean = mean(lux),
            luxsd = sd(lux))

hobo %>% 
  filter(date_time>as_datetime("2020-08-16 23:00:00"),
         lux!=0) %>% 
  summarise(tempmean = mean(temp),
            tempsd = sd(temp),
            luxmean = mean(lux),
            luxsd = sd(lux))
  
#====Number of midges at sample site 3 August 2020====
corearea <- 2.5^2*pi/10000 #kajak core area in m2

data.frame(fract_count = rep(1/8,3), tanyt = c(20, 12, 21)) %>% 
  mutate(tanyt_tot = tanyt/fract_count) %>% 
  summarise(mean = mean(tanyt_tot)/corearea,
            sd = sd(tanyt_tot)/corearea)
  
#=====Write to CSV====
# write_csv(hobo, "clean data/MG_hobo.csv")
# write_csv(inc_timing, "clean data/MG_inc_timing.csv")
# write_csv(meta, "clean data/MG_meta.csv")
# write_csv(chl, "clean data/MG_chl.csv")
# write_csv(loi, "clean data/MG_om.csv")
# write_csv(nep, "clean data/MG_ep.csv")
# write_csv(cm, "clean data/MG_cm.csv")
# write_csv(cc, "clean data/MG_cc.csv")
