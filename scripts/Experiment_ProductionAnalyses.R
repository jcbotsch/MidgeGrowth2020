#====Load Packages====
library(tidyverse)
library(lubridate)
library(lme4)
library(car)
library(phytools)
source("scripts/MG_Functions.R")
library(shades)

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
         box = as.character(box),
         gpp = gpp*1000) #convert GPP from gm-2h-1 to mgm-2h-1

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
  #add initial densities
  bind_rows (cm %>% 
               filter(day ==0) %>% 
               group_by(coreid, day) %>% 
               count(name = "live_tt") %>% 
               full_join(meta %>% select(algae_conc2) %>% unique(), by = character()) %>% 
               mutate(midge = "Midges")) %>% 
  group_by(algae_conc2, day) %>% 
  summarise(nt = mean(live_tt))

#=====Calculate Midge Production=====
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

# initials for non-parametric bootstrapping
startboot <- list()
expboot <- list() 
nboot = 1000

#perform non-parametric bootstrapping
#bootstrap initial samples
for(i in 1:nboot){
  startboot[[i]] <- start_cm %>% 
    group_by(day) %>% 
    sample_n(size = 20, replace = TRUE) %>% 
    summarise(body_size = mean(body_size)) %>% 
    mutate(bootstrap = i)
}
sb <- bind_rows(startboot)

#bootstrap samples during experiment
for(i in 1:nboot){
  expboot[[i]] <- exp_cm %>% 
    group_by(day, algae_conc2) %>% 
    add_count() %>% 
    sample_n(size = n, replace = TRUE) %>% 
    summarise(body_size = mean(body_size)) %>%  
    mutate(bootstrap = i) 
}
eb <- bind_rows(expboot)

#combine bootstraps
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

#extract initial bootstrapped sizes
day0 <- estimated_growth %>% 
  mutate(wt = weight(body_size)) %>% 
  left_join(nts) %>% 
  filter(day == 0) %>% 
  rename(body_size_day = body_size, nt_day0 = nt, wt_day0 = wt) %>% 
  select(bootstrap, algae_conc2, contains("day0"))

# calculate growth and production on the bootstraps
production_boot <- estimated_growth %>% 
  filter(day!=0) %>% 
  mutate(wt = weight(body_size)) %>% 
  left_join(nts) %>% 
  left_join(day0) %>% 
  group_by(bootstrap, algae_conc2) %>% 
  mutate(Pd = incsumprod(n1 = nt_day0, n2 = nt, wt1 = wt_day0, wt2 = wt, day)$Pd,
         g = incsumprod(n1 = nt_day0, n2 = nt, wt1 = wt_day0, wt2 = wt, day)$g, #average growth in mg/ind
         gd = g/day) 

#daily secondary production of midges on day 22
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

#=====Prepare Data====
prod_cm <- nep %>% 
  group_by(day, midge, algae_conc2) %>% 
  mutate(prod = gpp_omgm2h_to_cugcm2d(gpp)) %>% 
  add_count() %>% 
  summarise(n = unique(n),
            mean_pp = mean(prod), #average GPP for a given treatment midge, day for each initial algal abundance
            sd_pp = sd(prod)/sqrt(n)) %>% #calculate standard error of GPP
  full_join(production_boot %>% #join with bootstrapped production data
              mutate(Pdc = ((Pd*0.5)/1000)/ftube_area, #midge production in C m^-2
                     gdc = gd*0.5*1000, #convert to ug
                     midge = "Midges") %>%  #0.5 g C/ g AFDM
              group_by(day, algae_conc2, midge) %>% 
              summarise(mean_sp = mean(Pdc), #secondary production
                        sd_sp = sd(Pdc), #standard error of secondary production
                        Bt = mean(nt*wt), #calculate biomass 
                        mean_gdc = mean(gdc), #calculate growth in C
                        sd_gdc = sd(gdc))) %>% #calculate standard error of growth
  filter(day!=0,
         midge == "Midges") %>% 
  ungroup %>% 
  mutate(unique_id =1:n())

#====Fit Measurement Error Model====
#setup covariance matrices
cov.PPc <- matrix(nrow  = 20, ncol = 20, 0)
diag(cov.PPc) <- (prod_cm$sd_pp)^2
cov.MG <- matrix(nrow  = 20, ncol = 20, 0)
diag(cov.MG) <- (prod_cm$sd_gdc)^2
cov.base <- matrix(nrow  = 20, ncol = 20, 0)

#create a phylogeny with no structure
startree <- starTree(prod_cm$unique_id, branch.lengths = rep(0.000000000001, 20))

#convert phylogeny to being rooted
startree <- multi2di(startree)

#check
plot.phylo(startree)

#fit measurement model
gmod2 <- pgls.Ives(startree, y = prod_cm$mean_pp, X = prod_cm$mean_gdc, Vy = cov.PPc, Vx = cov.MG, Cxy = cov.base)
gmod2$beta #slope and intercept

#====Figure: Primary and Secondary Production=====
estfit <- data.frame(mean_gdc = seq(min(prod_cm$mean_gdc-prod_cm$sd_gdc), max(prod_cm$mean_gdc + prod_cm$sd_gdc), by = 0.001)) %>% 
  mutate(mean_pp = gmod2$beta[1] + gmod2$beta[2]*mean_gdc)

#plot figure
growfig <- prod_cm %>% 
  mutate(day = paste("Day", day)) %>% 
  ggplot(aes(x = mean_gdc, y = mean_pp))+
  geom_line(data = estfit)+
  geom_errorbar(aes(ymin = mean_pp-sd_pp, ymax = mean_pp+sd_pp), alpha = 0.5)+
  geom_errorbar(aes(xmin = mean_gdc-sd_gdc, xmax = mean_gdc+sd_gdc), alpha = 0.5)+
  geom_point(aes(shape = day, fill = algae_conc2), alpha = 1, size = 2)+
  algae_fill+
  scale_shape_manual(values = c(21, 22))+
  labs(y = expression("Primary Production \u03BCg C"~cm^{-2}~d^{-1}),
       x = expression("Midge Growth \u03BCg C "~ind^{-1}~d^{-1}),
       fill = "Initial Algal Abundance",
       shape = element_blank())+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, barheight = 0.5, barwidth = 8))+
  theme(legend.box = "vertical")

ggpreview(plot = growfig, dpi = 800, width = 3, height = 3.5, units = "in")
# ggsave(plot = growfig, filename = "figures/Botsch_MidgeGrowth_fig2.pdf", dpi = 800, width = 3, height = 3.5, units = "in")
