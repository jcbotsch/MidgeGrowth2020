#====Load Packages====
library(tidyverse)
library(lubridate)
library(rsample)

#====Load Data=====
meta <- read_csv("clean data/MG_meta.csv") %>% 
  mutate(algae_conc2 = ifelse(algae_conc == 0, 2e-3, algae_conc)) #half of lowest value
ep_raw <- read_csv("clean data/MG_ep.csv")
cm_raw <- read_csv("clean data/MG_cm.csv")
cc_raw <- read_csv("clean data/MG_cc.csv")

#====Prep Data====
nep <- ep_raw %>% 
  left_join(meta) %>% 
  mutate(coreid = as.character(coreid),
         box = as.character(box))


cc <- cc_raw %>% 
  left_join(meta) %>% 
  mutate(coreid = as.character(coreid),
         box = as.character(box))

cm <- cm_raw %>% 
  left_join(meta) %>% 
  mutate(coreid = as.character(coreid),
         box = as.character(box),
         body_size = body_size/1000) #convert from um to mm

start_cm <- cm %>% 
  filter(day == 0) %>% 
  select(coreid, sampledate, day, species_name, head_size, body_size, instar) 

exp_cm <- cm %>% 
  filter(day!=0,
         !is.na(body_size),
         species_name == "tt") 

growth_ests <- bind_rows(start_cm, exp_cm)
  

#Lindegaard et al (1979) formula for converting length to biomass (tab.21)
#AFDW in mg^(1/3) = 0.0442 + 0.0879 * length in mm 
weight <- function(l){
  (0.0442 + 0.0879 * l )^3 
}

#increment Summation method to estimate production (Benke and Huryn 2017)
incsumprod <- function(n1, n2, wt1, wt2, deltaT){
  g = wt2-wt1 #average individual growth
  Nbar = (n1 + n2)/2 #Average abundance
  intp = Nbar*g #production over the interval (mg)
  Pd = intp/deltaT #daily production
  list(g = g, Nbar = Nbar, intp = intp, Pd = Pd)
}

start_cm %>% 
  group_by(coreid) %>% 
  add_count() %>% 
  ungroup %>% 
  sample_n(size = n, replace = TRUE) %>% 
  summarise(body_size = mean(body_size))


exp_cm %>%
  group_by(day, algae_conc2) %>% 
  add_count() %>% 
  sample_n(size = n, replace = TRUE) 

startboot <- list()

set.seed(1234)
nboot = 1000
for(i in 1:nboot){
  startboot[[i]] <- start_cm %>% 
    group_by(day) %>% 
    sample_n(size = 20, replace = TRUE) %>% 
    summarise(body_size = mean(body_size)) %>% 
    mutate(bootstrap = i)
}

sb <- bind_rows(startboot)

expboot <- list()

for(i in 1:nboot){
  expboot[[i]] <- exp_cm %>% 
    group_by(day, algae_conc2) %>% 
    add_count() %>% 
    sample_n(size = n, replace = TRUE) %>% 
    summarise(body_size = mean(body_size)) %>%  
    mutate(bootstrap = i) 
}

eb <- bind_rows(expboot)


estimated_growth <- sb %>% 
  full_join(meta %>% select(algae_conc2) %>% unique(), by = character()) %>% 
  full_join(eb) 



estimated_growth%>% 
  ggplot(aes(x = factor(day), y = body_size))+
  geom_violin()+
  facet_wrap(~algae_conc2)


nts <- cc %>% 
  bind_rows (cm %>% 
               filter(day ==0) %>% 
               group_by(coreid, day) %>% 
               count(name = "live_tt") %>% 
               full_join(meta %>% select(algae_conc2) %>% unique(), by = character())) %>% 
  group_by(algae_conc2, day) %>% 
  summarise(nt = mean(live_tt))




production_boot <- estimated_growth %>% 
  mutate(wt = weight(body_size)) %>% 
  left_join(nts) %>% 
  group_by(bootstrap, algae_conc2) %>% 
  mutate(ntlag = lag(nt),
         wtlag = lag(wt),
         deltaT = day-lag(day),
         Pd = incsumprod(n1 = ntlag, n2 = nt, wt1 = wtlag, wt2 = wt, deltaT)$Pd,
         g = incsumprod(n1 = ntlag, n2 = nt, wt1 = wtlag, wt2 = wt, deltaT)$g)

 
ggplot()+
  geom_violin(aes(x = algae_conc2, group = interaction(algae_conc2, day), y = Pd, fill = day),
              data = production_boot %>% 
                filter(!is.na(Pd)) %>% 
                mutate(day = ifelse(day == 14, "Day 0-14", "Day 14-22")))+
  scale_x_continuous(trans = "log", breaks = c(0.01, 0.1, 1))+
  labs(x = "Sediment Treatment",
       y = "Secondary Production",
       fill = element_blank())+
  geom_rect(aes(xmin = 0.0015, xmax = 1.2, ymin = -Inf, ymax = 0), fill = "white", alpha = 0.7)+
  geom_hline(yintercept = 0, linetype = "dashed")

ggplot()+
  geom_violin(aes(x = algae_conc2, group = interaction(algae_conc2, day), y = g, fill = day),
              data = production_boot %>% 
                filter(!is.na(Pd)) %>% 
                mutate(day = ifelse(day == 14, "Day 0-14", "Day 14-22")))+
  scale_x_continuous(trans = "log", breaks = c(0.01, 0.1, 1))+
  labs(x = "Sediment Treatment",
       y = "Average Individual Growth (mg)",
       fill = element_blank())+
  geom_rect(aes(xmin = 0.0015, xmax = 1.2, ymin = -Inf, ymax = 0), fill = "white", alpha = 0.7)+
  geom_hline(yintercept = 0, linetype = "dashed")
  
production_boot %>% 
  group_by(day, algae_conc2) %>% 
  ggplot(aes(x = algae_conc2, group = interaction(algae_conc2, day), y = wt, fill = factor(day)))+
    geom_violin()+
    geom_segment(aes(y = lag_avg_wt, yend = avg_wt, x = algae_conc21, xend = algae_conc22), data = . %>% 
                   group_by(day, algae_conc2) %>% 
                   summarise(avg_wt = mean(wt)) %>%
                   group_by(algae_conc2) %>% 
                   mutate(lag_avg_wt = lag(avg_wt),
                          algae_conc21 = ifelse(day == 14, algae_conc2-0.2*algae_conc2, algae_conc2),
                          algae_conc22 = ifelse(day == 22, algae_conc2+0.25*algae_conc2, algae_conc2)))+
    scale_x_continuous(trans = "log", breaks = c(0.01, 0.1, 1))+
  geom_point(aes(y = wt, x = algae_conc2, fill = factor(day)), 
              size = 1.5, alpha = 0.75, 
             position = position_jitterdodge(dodge.width = 0.6),
             data = cm %>% 
               filter(midge == "Midges") %>% 
               full_join(start_cm %>% 
                           full_join(meta %>% select(algae_conc2) %>% unique, by = character())) %>% 
               group_by(algae_conc2, day, coreid) %>% 
               summarise(wt = weight(meanna(body_size)))) +
    labs(x = "Sediment Treatment",
         y = "Average Individual AFDM (mg)",
         fill = "Day")


production_boot %>% 
  mutate(day = paste("Day", day, sep = " ")) %>% 
  group_by(day, algae_conc2) %>% 
  ggplot(aes(x = algae_conc2, group = interaction(algae_conc2, day), y = wt, fill = factor(day)))+
  geom_violin()+
  geom_segment(aes(y = lag_avg_wt, yend = avg_wt, x = algae_conc21, xend = algae_conc22), data = . %>% 
                 group_by(day, algae_conc2) %>% 
                 summarise(avg_wt = mean(wt)) %>%
                 group_by(algae_conc2) %>% 
                 mutate(lag_avg_wt = lag(avg_wt),
                        algae_conc21 = ifelse(day == "Day 14", algae_conc2-0.2*algae_conc2, algae_conc2),
                        algae_conc22 = ifelse(day == "Day 22", algae_conc2+0.25*algae_conc2, algae_conc2)) %>% 
                 filter(!is.na(lag_avg_wt)))+
  scale_x_continuous(trans = "log", breaks = c(0.01, 0.1, 1))+
  labs(x = "Sediment Treatment",
       y = "Average Individual AFDM (mg)",
       fill = NULL)+
  theme(legend.position = c(0.3,0.9),
        legend.direction = "horizontal")+
  scale_fill_brewer(type = "seq", palette = 3)
# ggpreview(plot = last_plot(), dpi = 650, width = 120, height = 80, units = "mm")




production_boot %>% 
  mutate(day = paste("Day", day, sep = " ")) %>% 
  group_by(day, algae_conc2) %>% 
  ggplot(aes(group = interaction(algae_conc2, day), fill = factor(day)))+
  ggridges::geom_density_ridges2(aes(y = algae_conc2, x = wt), alpha = 0.75)+
  scale_y_continuous(trans = "log", breaks = c(0.01, 0.1, 1), limits = c(NA, 4))+
  labs(y = "Sediment Treatment",
       x = "Average Individual AFDM (mg)",
       fill = NULL)+
  theme(legend.position = c(0.5,0.95),
        legend.direction = "horizontal")+
  scale_fill_brewer(type = "seq", palette = 3)

ggpreview(plot = last_plot(), dpi = 650, width = 80, height =120, units = "mm")




cm %>% 
  group_by(coreid) %>% 
  mutate(wt = weight(body_size)) %>% 
  summarise(avg_wt = meanna(wt)) %>% 
  right_join(cc) %>% 
  mutate(est_resp = live_tt*avg_wt*7.7/1e6) %>% #convert to grams to match respiration
  full_join(nep) %>% 
  mutate(midge_contr = (est_resp/-resp)*100) %>% #percent of respiration likely caused by midges
  select(midge_contr) %>% 
  summary()
  
nep %>% 
  group_by(day, midge, algae_conc2) %>% 
  mutate(prod = (gpp*18)*(12/32*1.2)+resp*24*(12/32*0.8)) %>% #photosynthetic quotient of 
  add_count() %>% 
  mutate(n = unique(n),
            mean = mean(prod),
            sd = sd(prod)) %>% 
  ggplot()+
  scale_x_continuous(trans = "log", breaks = c(0.01, 0.1, 1))+
  geom_violin(aes(x = algae_conc2, group = interaction(algae_conc2, day), y = Pd*.5, fill = factor(day)),
              data = production_boot %>% mutate(midge = "Midges") %>% 
                filter(!is.na(Pd)))+
  geom_pointrange(aes(x = algae_conc2, y = mean, ymin = mean-sd, ymax = mean+sd, shape = midge, fill = factor(day)), position = position_dodge(width = 0.2))+
  scale_shape_manual(values = c(21, 22))+
  # geom_point(aes(x = algae_conc2, y = prod, fill = factor(day)), shape = 21, position = position_dodge(width = 0.2))+
  labs(y = "Production (mg C d^-1)",
       x = "Sediment Treatment")

ftube_r = 0.03/2 #30 mm/ 2 to get radius /1000 to convert to m
ftube_area = ftube_r^2*pi
pq = 1

prod1 <- nep %>% 
  group_by(day, midge, algae_conc2) %>% 
  mutate(prod = (gpp*18)*(12/32*pq)) %>% #daily production of algae in g C m^2d^-1
  add_count() %>% 
  summarise(n = unique(n),
         mean_pp = mean(prod),
         sd_pp = sd(prod)) %>% 
  full_join(production_boot %>% 
              # filter(Pd>0) %>% 
              mutate(Pdc = ((Pd*0.5)/1000)/ftube_area,
                     midge = "Midges") %>%  #0.5 g C/ g AFDM
              group_by(day, algae_conc2, midge) %>% 
              summarise(mean_sp = mean(Pdc),
                        sd_sp = sd(Pdc))) %>% 
  filter(day!=0,
         midge == "Midges") 


prod1 %>% lm(mean_sp~mean_pp:factor(day), data = .) %>% summary()


prod1 %>% 
  ggplot(aes(x = mean_pp, y = mean_sp))+
  facet_wrap(~day)+
  geom_errorbarh(aes(xmin = mean_pp-sd_pp, xmax = mean_pp+sd_pp))+
  geom_errorbar(aes(ymin = mean_sp-sd_sp, ymax = mean_sp+sd_sp), width  = NA)+
  geom_point(aes(fill = algae_conc2), shape = 21, size = 2)+
  viridis::scale_fill_viridis(trans = "log", breaks = c(0.01, 0.1, 1))+
  coord_cartesian(xlim = c(0, NA),
                  ylim = c(0, NA))
  

prod1 %>% 
  ungroup %>% 
  select(day, algae_conc2, contains("mean")) %>% 
  gather(production, val, contains("mean")) %>% 
  ggplot(aes(x = algae_conc2, y = val, fill = production))+
  # geom_errorbar(ymin = val, ymax = val+sd)+
  facet_wrap(~day)+
  geom_col(position = "dodge")+
  scale_x_continuous(trans = "log")
