#====Load Packages====
library(tidyverse)
library(lubridate)
library(lme4)
library(car)
library(phytools)
source("scripts/MG_Functions.R")

options(dplyr.summarise.inform = FALSE)
#====read in files====
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

cm <- cm_raw %>% 
  left_join(meta) %>% 
  mutate(coreid = as.character(coreid),
         box = as.character(box),
         body_size = body_size/1000) #convert from um to mm


cc <- cc_raw %>% 
  left_join(meta) %>% 
  mutate(coreid = as.character(coreid),
         box = as.character(box))

#number of individuals on a given day in a given treatment
nts <- cc %>% 
  filter(midge == "Midges") %>% 
  bind_rows (cm %>% 
               filter(day ==0) %>% 
               group_by(coreid, day) %>% 
               count(name = "live_tt") %>% 
               full_join(meta %>% select(algae_conc2) %>% unique(), by = character()) %>% 
               mutate(midge = "Midges")) %>% 
  group_by(algae_conc2, day) %>% 
  summarise(nt = mean(live_tt))



#====Plotting boimass====

cm %>% 
  filter(day>0) %>% 
  mutate(wt = weight(body_size)) %>% #mg
  group_by(day, coreid) %>% 
  summarise(Bt = ((sumna(wt) *0.5)/1000)/ftube_area) %>% 
  left_join(nep %>% 
              mutate(gpp = gpp*((12/32)*pq))) %>% 
  group_by(day, midge, algae_conc2) %>% 
  arrange(algae_conc2) %>% 
  # summarise(Bt= mean(Bt, na.rm = T),
  #           gpp = mean(gpp, na.rm = T)) %>% 
  ggplot(aes(x = Bt, y = gpp))+
  # geom_path()+
  facet_grid(midge~day)+
  geom_rug(length = unit(0.3, "lines"), y = NA, data = cm %>% 
             filter(day == 0) %>% 
             mutate(wt = weight(body_size)) %>% 
             group_by(coreid) %>% 
             summarise(Bt = sumna(wt)))+
  geom_point(aes(fill = algae_conc2), size = 2, shape = 21)+
  scale_fill_viridis_c(trans = "log", breaks = c(0.01, 0.1, 1))

#=====Calculate Secondary Production=====
set.seed(12345)
#subset starting midges and experimental midges
start_cm <- cm %>% 
  filter(day == 0) %>% 
  select(coreid, sampledate, day, species_name, head_size, body_size, instar) 

exp_cm <- cm %>% 
  filter(day!=0,
         midge == "Midges",
         !is.na(body_size),
         species_name == "tt") 



# perform non-parametric bootstrapping
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

#number of individuals on a given day in a given treatment
nts <- cc %>% 
  filter(midge == "Midges") %>% 
  bind_rows (cm %>% 
               filter(day ==0) %>% 
               group_by(coreid, day) %>% 
               count(name = "live_tt") %>% 
               full_join(meta %>% select(algae_conc2) %>% unique(), by = character()) %>% 
               mutate(midge = "Midges")) %>% 
  group_by(algae_conc2, day) %>% 
  summarise(nt = mean(live_tt))

#calculate production and growth
day0 <- estimated_growth %>% 
  mutate(wt = weight(body_size)) %>% 
  left_join(nts) %>% 
  filter(day == 0) %>% 
  rename(body_size_day = body_size, nt_day0 = nt, wt_day0 = wt) %>% 
  select(bootstrap, algae_conc2, contains("day0"))

production_boot2 <- estimated_growth %>% 
  filter(day!=0) %>% 
  mutate(wt = weight(body_size)) %>% 
  left_join(nts) %>% 
  left_join(day0) %>% 
  group_by(bootstrap, algae_conc2) %>% 
  mutate(Pd = incsumprod(n1 = nt_day0, n2 = nt, wt1 = wt_day0, wt2 = wt, day)$Pd,
         g = incsumprod(n1 = nt_day0, n2 = nt, wt1 = wt_day0, wt2 = wt, day)$g, #average growth in mg/ind
         gd = g/day) 


total_sp <- eb %>% 
  filter(day == 22) %>% 
  left_join(meta) %>% 
  full_join(sb %>% select(-day) %>% rename(body_size1 = body_size)) %>% 
  mutate(wt2 = weight(body_size),
         wt1 = weight(body_size1),
         deltaT = 22) %>% 
  left_join(nts %>% filter(day == 22)) %>% 
  mutate(nt1 = unique({nts %>% filter(day == 0)}$nt)) %>% 
  mutate(Pd = incsumprod(n1 = nt1, n2 = nt, wt1 = wt1, wt2 = wt2, deltaT)$Pd,
         g = incsumprod(n1 = nt1, n2 = nt, wt1 = wt1, wt2 = wt2, deltaT)$g,
         gd = g/deltaT) %>% 
  group_by(algae_conc2, midge, day) %>% 
  mutate(Pdc = ((Pd*0.5)/1000)/ftube_area) %>% 
  summarise(mean_sp = mean(Pdc),
            sd_sp = sd(Pdc))


prod2 <- nep %>% 
  group_by(day, midge, algae_conc2) %>% 
  mutate(prod = (gpp*18)*(12/32*pq) ) %>%#daily production of algae in g C m^2d^-1
  add_count() %>% 
  summarise(n = unique(n),
            mean_pp = mean(prod),
            sd_pp = sd(prod)/sqrt(n)) %>% 
  full_join(production_boot2 %>% 
              # filter(Pd>0) %>% 
              mutate(Pdc = ((Pd*0.5)/1000)/ftube_area,
                     gdc = gd*0.5*1000, #convert to ug
                     midge = "Midges") %>%  #0.5 g C/ g AFDM
              group_by(day, algae_conc2, midge) %>% 
              summarise(mean_sp = mean(Pdc),
                        sd_sp = sd(Pdc),
                        Bt = mean(nt*wt),
                        mean_gdc = mean(gdc),
                        sd_gdc = sd(gdc))) %>% 
  filter(day!=0,
         midge == "Midges") %>% 
  ungroup %>% 
  mutate(unique_id =1:n())



#====Fit Measurement Error Model====
cov.PP <- matrix(nrow  = 20, ncol = 20, 0)

diag(cov.PP) <- (prod2$sd_pp)^2

cov.SP <- matrix(nrow  = 20, ncol = 20, 0)

diag(cov.SP) <- (prod2$sd_sp)^2

cov.base <- matrix(nrow  = 20, ncol = 20, 0)

#create a phylogeny with no structure
startree <- starTree(prod2$unique_id, branch.lengths = rep(0.000000000001, 20))

starTree(prod2$unique_id) %>% multi2di() %>% plot.phylo()

#convert phylogeny to being rooted
startree <- multi2di(startree)

#check
plot.phylo(startree)

prodmod <- pgls.Ives(startree, X = prod2$mean_pp, y = prod2$mean_sp, Vx = cov.PP, Vy = cov.SP, Cxy = cov.base)

prodmod$beta

-prodmod$beta[1]/prodmod$beta[2] #=x at y = 0

#====Figure: Primary and Secondary Production=====
estfit <- data.frame(mean_pp = seq(min(prod2$mean_pp-prod2$sd_pp), max(prod2$mean_pp + prod2$sd_pp), by = 0.001)) %>% 
  mutate(mean_sp = prodmod$beta[1] + prodmod$beta[2]*mean_pp)

prodfig <- prod2 %>% 
  mutate(day = paste("Day", day)) %>% 
  ggplot(aes(x = mean_pp, y = mean_sp))+
  geom_line(data = estfit)+
  # geom_abline(slope = prodmod$beta[2], intercept = prodmod$beta[1])+
  # geom_ribbon(aes(ymin = psp-psp_se, y = psp, ymax = psp+psp_se, group = day), fill = "gray20", alpha = 0.2, data = nd %>% mutate(day = paste("Day", day)))+
  geom_errorbarh(aes(xmin = mean_pp-sd_pp, xmax = mean_pp+sd_pp), alpha = 0.75)+
  geom_errorbar(aes(ymin = mean_sp-sd_sp, ymax = mean_sp+sd_sp), alpha = 0.75, width  = NA)+
  geom_point(aes(shape = day, fill = algae_conc2), size = 2)+
  viridis::scale_fill_viridis(trans = "log", breaks = c(0.01, 0.1, 1))+
  # geom_line(aes(y = psp, linetype = day), data = nd %>% mutate(day = paste("Day", day)))+
  coord_cartesian(xlim = c(0, NA))+
  scale_shape_manual(values = c(21, 22))+
  labs(x = expression("Primary Production g C"~m^{-2}~d^{-1}),
       y = expression("Secondary Production g C"~m^{-2}~d^{-1}),
       fill = "Sediment Treatment",
       shape = element_blank(),
       linetype = element_blank())+
  guides(fill = guide_colorbar(title.position = "top"))+
  theme(legend.box = "vertical",
        legend.spacing = unit(0,units = 'points'),
        legend.box.spacing = unit(0, units = "points"))
plot(prodfig)

# ggpreview(plot = prodfig, dpi = 650, width = 80, height = 100, units = "mm")

###################################################
#=====Reproduce with midge growth====

prod_cm <- nep %>% 
  group_by(day, midge, algae_conc2) %>% 
  mutate(prod = (gpp*18)*(12/32*pq),#daily production of algae in g C m^2d^-1
         prod = prod*100) %>% #convert g C m^-2d^-1 to ug C cm^-2 d^-1
  add_count() %>% 
  summarise(n = unique(n),
            mean_pp = mean(prod),
            sd_pp = sd(prod)/sqrt(n)) %>% 
  full_join(production_boot2 %>% 
              # filter(Pd>0) %>% 
              mutate(Pdc = ((Pd*0.5)/1000)/ftube_area,
                     gdc = gd*0.5*1000, #convert to ug
                     midge = "Midges") %>%  #0.5 g C/ g AFDM
              group_by(day, algae_conc2, midge) %>% 
              summarise(mean_sp = mean(Pdc),
                        sd_sp = sd(Pdc),
                        Bt = mean(nt*wt),
                        mean_gdc = mean(gdc),
                        sd_gdc = sd(gdc))) %>% 
  filter(day!=0,
         midge == "Midges") %>% 
  ungroup %>% 
  mutate(unique_id =1:n())

#====Fit Measurement Error Model====
cov.PPc <- matrix(nrow  = 20, ncol = 20, 0)

diag(cov.PPc) <- (prod_cm$sd_pp)^2

cov.MG <- matrix(nrow  = 20, ncol = 20, 0)

diag(cov.MG) <- (prod_cm$sd_gdc)^2

cov.base <- matrix(nrow  = 20, ncol = 20, 0)

#create a phylogeny with no structure
startree <- starTree(prod_cm$unique_id, branch.lengths = rep(0.000000000001, 20))

starTree(prod_cm$unique_id) %>% multi2di() %>% plot.phylo()

#convert phylogeny to being rooted
startree <- multi2di(startree)

#check
plot.phylo(startree)

gmod <- pgls.Ives(startree, X = prod_cm$mean_pp, y = prod_cm$mean_gdc, Vx = cov.PPc, Vy = cov.MG, Cxy = cov.base)

gmod$beta

-gmod$beta[1]/gmod$beta[2] #=x at y = 0

#====Figure: Primary and Secondary Production=====
estfit <- data.frame(mean_pp = seq(min(prod_cm$mean_pp-prod_cm$sd_pp), max(prod_cm$mean_pp + prod_cm$sd_pp), by = 0.001)) %>% 
  mutate(mean_gdc = gmod$beta[1] + gmod$beta[2]*mean_pp)



growfig <- prod_cm%>% 
  mutate(day = paste("Day", day)) %>% 
  ggplot(aes(x = mean_pp, y = mean_gdc))+
  geom_line(data = estfit)+
  geom_errorbarh(aes(xmin = mean_pp-sd_pp, xmax = mean_pp+sd_pp), alpha = 0.5)+
  geom_errorbar(aes(ymin = mean_gdc-sd_gdc, ymax = mean_gdc+sd_gdc), alpha = 0.5, width  = NA)+
  geom_point(aes(shape = day, fill = algae_conc2), alpha = 0.75, size = 2)+
  viridis::scale_fill_viridis(trans = "log", breaks = c(0.01, 0.1, 1))+
  coord_cartesian(xlim = c(0, NA))+
  scale_y_continuous(labels = scales::comma_format())+
  scale_shape_manual(values = c(21, 22))+
  labs(x = expression("Primary Production \u03BCg C"~cm^{-2}~d^{-1}),
       y = expression("Midge Growth \u03BCg C "~ind^{-1}~d^{-1}),
       fill = "Sediment Treatment",
       shape = element_blank(),
       linetype = element_blank())+
  guides(fill = guide_colorbar(title.position = "top"))+
  theme(legend.box = "vertical",
        legend.spacing = unit(0,units = 'points'),
        legend.box.spacing = unit(0, units = "points"))

ggpreview(plot = growfig, dpi = 650, width = 80, height = 100, units = "mm")
