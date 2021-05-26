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

ggpreview <- function(...) {
  fname <- tempfile(fileext = ".png")
  ggsave(filename = fname, ...)
  system2("open", fname)
  invisible(NULL)
}

#====set aesthetics====
theme_set(theme_bw()+
            theme(panel.grid = element_blank(),
                  strip.background = element_blank(),
                  legend.position = "bottom"))

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

# add temp (unsaved. )
nep %>% 
  left_join(nep1 %>% 
              select(coreid, date, dark_light, contains("wtemp")) %>% 
              mutate(wtemp = (wtemp_init_cor+wtemp_final)/2) %>% 
              select(coreid, date, dark_light, wtemp) %>% 
              spread(dark_light, wtemp) %>% 
              rename(resp_temp = dark, nep_temp = light)) 


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
mg_instars = c(0,4.3,7.3,12, 18)/55*1000 #divide by ocular micrometer units to get mm then 1000 to get micrometers
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
  filter(head_size>70) %>%  #first instar head size is 75. there is one individual with a lower head size that was damaged.
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

hobo <- hobo_raw %>% 
  rename(date_time = 2,
         temp = 3,
         lux = 4) %>% 
  select(date_time, temp, lux) %>% 
  mutate(date_time = mdy_hms(date_time)) %>% 
  filter(date_time>as_datetime("2020-08-07 15:00:00"),
         date_time<as_datetime("2020-08-29 10:50:00"))


hobo %>% 
  ggplot(aes(x = date_time, y = lux))+
  geom_line()




hobo %>% 
  gather(var, val, -date_time) %>% 
  bind_cols(nep1 %>% 
              group_by(date) %>% 
              summarise(dttm_start = min(dttm_start),
                        dttm_end = max(dttm_end)) %>% 
              gather(type, time, -date) %>% 
              mutate(incub = ifelse(date == "2020-08-18", 1, 2),
                     type = paste0(type, incub)) %>% 
              select(type, time) %>% 
              spread(type, time)) %>% 
  filter(date_time>dttm_start1&date_time<dttm_end1|
           date_time>dttm_start2 & date_time<dttm_end2) %>%
  mutate(group = ifelse(date_time<dttm_end1, "Day 14", "Day 22")) %>% 
  ggplot()+
  facet_wrap(var~group, scales = "free", strip.position = "left")+
  geom_rect(aes(xmin = dttm_start, xmax =dttm_end, fill = dark_light, ymin = -Inf, ymax = Inf), 
            data = nep1 %>% 
              group_by(dark_light, date) %>% 
              summarise(dttm_start = min(dttm_start),
                        dttm_end = max(dttm_end)) %>% 
              mutate(group = ifelse(date == "2020-08-18", "Day 14", "Day 22")))+
  geom_line(aes(x = date_time, y = val))+
  theme(strip.placement = "outside")+
  labs(x = "Date",
       y = NULL,
       fill = "")





hobo %>% 
  gather(var, val, -date_time) %>%
  mutate(var = ifelse(var == "temp", "Air Temperature (°C)",
                      "Lux")) %>% 
  ggplot()+
  facet_wrap(~var, scales = "free_y", strip.position = "left", ncol = 1)+
  geom_rect(aes(xmin = dttm_start, xmax =dttm_end, fill = dark_light, ymin = -Inf, ymax = Inf), 
            data = nep1 %>% 
              group_by(dark_light, date) %>% 
              summarise(dttm_start = min(dttm_start),
                        dttm_end = max(dttm_end)),
            alpha = 0.8)+
  geom_line(aes(x = date_time, y = val))+
  theme(strip.placement = "outside",
        legend.position = "none")+
  labs(x = "Date",
       y = NULL,
       fill = "")+
  scale_fill_viridis_d()

# ggpreview(plot = last_plot(), dpi = 650, width = 120, height = 80, units = "mm")



hobo %>% 
  filter(date_time>as.Date("2020-08-10")) %>% 
  summarise(tempmean = mean(temp),
            tempsd = sd(temp))


hobo %>% 
  group_by(hour = hour(date_time)) %>% 
  summarise(temp = mean(temp),
            lux = mean(lux)) %>%
  gather(var, val, -hour) %>% 
  ggplot(aes(x = hour, y = val))+
  facet_wrap(~var, scales = "free")+
  geom_line()


hobo %>% 
  mutate(hour = hour(date_time),
         date = as.Date(date_time)) %>% 
  gather(var, val, -hour, -date, -date_time) %>% 
  filter(! date %in% unique(c(as.Date(nep1$dttm_start), as.Date(nep1$dttm_end)))) %>% 
  group_by(date, hour, var) %>% 
  summarise(val = mean(val)) %>% 
  mutate(var = ifelse(var == "lux", "Irradiance (Lux)", "Temperature (°C)"),
         box = ifelse(date>="2020-08-20", "Top", "Bottom")) %>% 
  ggplot(aes(x = hour, y = val, group = date, col = date))+
  facet_grid(var~box, scales = "free_y", switch = "y")+
  geom_line(alpha = 0.5, size = 0.7)+
  scale_color_gradient2(low = "dodgerblue", high = "firebrick", trans = "date", mid = "darkorchid", midpoint = as.numeric(mean(as.Date(hobo$date_time))))+
  labs(y = "",
       x = "Hour",
       color = NULL)+
  theme(
        strip.placement = "outside")
ggpreview(plot = last_plot(), dpi = 650, width = 120, height = 80, units = "mm")


hobo %>% 
  mutate(hour = hour(date_time),
         date = as.Date(date_time)) %>% 
  group_by(date) %>% 
  gather(var, val, -hour, -date, -date_time) %>% 
  filter(! date %in% unique(c(as.Date(nep1$dttm_start), as.Date(nep1$dttm_end)))) %>% 
  group_by(date, var) %>% 
  summarise(max = max(val),
            min = min(val),
            starthr = min(hour),
            endhr = max(hour)) %>% 
  filter(starthr == 0, endhr == 23) %>% 
  group_by(var) %>% 
  summarise(day = mean(max),
            night = mean(min),
            sdday = sd(max),
            sdnight = sd(min))
  


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
data.frame(fract_count = rep(1/8,3), tanyt = c(20, 12, 21)) %>% 
  mutate(tanyt_tot = tanyt/fract_count) %>% 
  summarise(mean = mean(tanyt_tot)/corearea,
            sd = sd(tanyt_tot)/corearea)
  



#=====Write to CSV====


# write_csv(meta, "clean data/MG_meta.csv")
# write_csv(chl, "clean data/MG_chl.csv")
# write_csv(loi, "clean data/MG_om.csv")
# write_csv(ndvi, "clean data/MG_ndvi.csv")
# write_csv(nep, "clean data/MG_ep.csv")
# write_csv(cm, "clean data/MG_cm.csv")
# write_csv(cc, "clean data/MG_cc.csv")
