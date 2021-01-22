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


fig1 <- grid.arrange(corrchlplot, omplot, nrow  = 2)

# ggpreview(plot = fig1, dpi = 650, width = 80, height = 100, units = "mm")
# ggsave(plot = fig1, filename = "Botsch_MG_Fig1.pdf", dpi = 650, device = "pdf", width = 80, height = 80, units = "mm")



grid.arrange(uncorrchlplot, omplot, nrow  = 2, )


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


#====Do the slopes differ between the two sample dates====


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
# gboot1 <- para.boot(data1 = nep14, data2 = nep22,
#           response.var = "gpp",
#           lm1 = g14log, lm2 = g22log,
#           results = gboot1, nboot = nboot)
# 


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

#Figure 2: GPP 
nepv2.1 <- nep %>% 
  mutate(day = paste("Day", day)) %>% 
  ggplot(aes(x = algae_conc2, y = gpp, fill = midge))+
  facet_wrap(~day, ncol = 1)+
  geom_point(size  = 2, alpha = 0.7, shape = 21)+
  geom_smooth(aes(color = midge), method = "lm", size = 0.75, se = FALSE)+
  scale_fill_manual(values = c("black", "gray60"))+
  scale_color_manual(values = c("black", "gray60"))+
  labs(x = "Initial Sediment Quality",
       y = expression(GPP(g~O[2]~m^{-2}~hr^{-1})),
       color = "",
       fill = "")+
  scale_x_log10()+
  theme(legend.position = c(0.2,0.95),
        legend.text.align = 0,
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA, color = NA))

gppbootfig <- gboot1 %>%
  select(algae, algae_x_NoMidge, day) %>%
  gather(var, val, -day) %>% 
  group_by(var, day) %>% 
  mutate(mean = mean(val),
         var = ifelse(var == 'algae_x_NoMidge', "No Midges", "Midges"),
         label.x = median(val),
         label.y = ifelse(day == "Day 14", 1000, 600)) %>% 
  ggplot(aes(val, fill = day))+
  geom_density(alpha = 0.5, color = "gray50")+
  # geom_density(alpha = 0.5, color = "gray50", data =. %>% filter(day=="Day 22"))+
  geom_vline(xintercept = 0, linetype = "dashed")+
  facet_wrap(~var)+
  geom_text(aes(label = day, y = label.y, x = label.x, color = day))+
  labs(x = "Effect of Inital Sediment Quality on GPP",
       y = "",
       fill = "")+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")+
  scale_fill_manual(values = c("#F4A582", "#CA0020"))+
  scale_color_manual(values = c("#F4A582", "#CA0020"))+
  scale_x_continuous(breaks = c(0.003,0, -0.003))

cowplot::plot_grid(nepv2.1, gppbootfig, nrow = 2, rel_heights = c(2/3,1/3))

# ggpreview(plot = last_plot(), dpi = 650, width = 80, height = 120, units = "mm")
# ggsave(plot = last_plot(), filename = "Botsch_MG_Fig2.pdf", device = "pdf", dpi = 650, width = 80, height = 120, units = "mm")

gboot1 %>% 
  filter(day == "Day 22") %>%
  count(algae_x_NoMidge<0) %>% 
  mutate(p = n/2000)


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
         metabolism = ifelse(metabolism == "RESP", "R", metabolism),
         metabolism = fct_reorder(metabolism, -value),
         value = ifelse(metabolism == "R", -value, value)) %>% 
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

#how many midges out of the total stocked (1,000) survived?
cc %>% 
  filter(!coreid %in% c(108, 119, 152)) %>% 
  summarise(sum(live_tt))

cc %>% 
  filter(coreid %in% c(108, 119, 152)) %>% 
  summarise(stocked = sum(live_tt)/3) %>% 
  mutate(stocked.tot = stocked*50)

646/883

#====Figure 3: Number of Midges====
cc %>% 
  left_join(nep) %>% 
  filter(day %in% c(14, 22)) %>% 
  mutate(day = paste("Day", day)) %>% 
  ggplot(aes(x = algae_conc2, y = live_tt, fill = gpp))+
  facet_grid(midge~day, scales = "free")+
  geom_jitter(shape = 21, size = 2, alpha = 0.6, width = 0.0005, height = 0.3)+
  geom_smooth(method = "lm", se = FALSE, size = 0.6, col = "black")+
  scale_fill_viridis_c( breaks = c(0, 0.025, 0.05))+
  labs(y = "Live Tanytarsini",
       x = "Initial Sediment Quality",
       fill = expression(GPP(g~O[2]~m^{-2}~hr^{-1})~" "))+
  scale_x_log10()

fig3.v1 <- cc %>% 
  filter(day %in% c(14, 22)) %>% 
  mutate(day = paste("Day", day)) %>% 
  ggplot(aes(x = algae_conc2, y = live_tt, fill = midge))+
  facet_wrap(~day, ncol = 1)+
  geom_point(size  = 2, alpha = 0.7, shape = 21)+
  geom_hline(yintercept = 20, linetype = "dashed", alpha = 0.5)+
  geom_smooth(aes(color = midge), method = "lm", size = 0.75, se = FALSE)+
  scale_fill_manual(values = c("black", "gray60"))+
  scale_color_manual(values = c("black", "gray60"))+
  labs(x = "Initial Sediment Quality",
       y = "Live Tanytarsini",
       color = "",
       fill = "")+
  scale_x_log10()+
  theme(legend.position = c(0.5,0.95),
        legend.text.align = 0,
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA, color = NA),
        legend.direction = "horizontal")



fig3.v2 <- cc %>% 
  filter(day %in% c(14, 22)) %>% 
  mutate(day = paste("Day", day)) %>% 
  ggplot(aes(x = algae_conc2, y = live_tt, fill = day))+
  facet_wrap(~midge, ncol = 1, scales = "free")+
  geom_jitter(size  = 2, alpha = 0.5, shape = 21, width = 0.0005, height = 0)+
  geom_smooth(aes(color = day), method = "lm", size = 0.75, alpha = 0.5, se = FALSE)+
  scale_fill_manual(values = c("#F4A582", "#CA0020"))+
  scale_color_manual(values = c("#F4A582", "#CA0020"))+
  labs(x = "Initial Sediment Quality",
       y = "Live Tanytarsini",
       color = "",
       fill = "")+
  lims(y = c(0,NA))+
  scale_x_log10()+
  theme(legend.position = c(0.3,0.95),
        legend.text.align = 0,
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA, color = NA),
        legend.direction = "horizontal")

# ggpreview(plot = last_plot(), dpi = 650, width = 80, height = 120, units = "mm")
# ggsave(plot = last_plot(), filename = "Botsch_MG_Fig3.pdf", device = "pdf", dpi = 650, width = 80, height = 120, units = "mm")


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
  scale_linetype_manual(values = c("dashed", "solid"), guide = NULL)+
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
# nmboot <- para.boot(data1 = cc14, data2 = cc22, 
#           response.var = "live_tt", 
#           lm1 = nm142, lm2 = nm222, 
#           results = nmboot, nboot = nboot)


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
         var = ifelse(var == 'algae_x_NoMidge', "No Midges", "Midges"),
         var = paste(day, var),
         label.x = median(val),
         label.y = ifelse(day == "Day 14", ifelse(var == "Midges", 12, 4), ifelse(var == "Midges", 15, 6))) %>% 
  ggplot(aes(val, fill = var))+
  geom_density(alpha = 0.5, color = "gray50")+
  geom_vline(xintercept = 0, linetype = "dashed")+
  labs(x = "Effect of Inital Sediment Quality \non Number of Midges",
       y = "",
       fill = "")+
  lims(y = c(0, 16))+
  scale_x_continuous(breaks = c(0, 0.5, 1.0))

nmboot %>%
  select(algae, algae_x_NoMidge, day) %>%
  gather(var, val, -day) %>% 
  group_by(var, day) %>% 
  summarise(mean = mean(val),
            se = sd(val))


nmbootfig <- nmboot %>%
  select(algae, algae_x_NoMidge, day) %>%
  gather(var, val, -day) %>% 
  group_by(var, day) %>% 
  mutate(mean = mean(val),
         var = ifelse(var == 'algae_x_NoMidge', "No Midges", "Midges"),
         label.x = median(val),
         label.y = ifelse(day == "Day 14", ifelse(var == "Midges", 12, 4), ifelse(var == "Midges", 15, 6))) %>% 
  ggplot(aes(val, fill = day))+
  geom_density(alpha = 0.5, color = "gray50")+
  geom_vline(xintercept = 0, linetype = "dashed")+
  facet_wrap(~var)+
  geom_text(aes(label = day, y = label.y, x = label.x, color = day))+
  labs(x = "Effect of Inital Sediment Quality \non Number of Midges",
       y = "",
       fill = "")+
  lims(y = c(0, 16))+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")+
  scale_fill_manual(values = c("#F4A582", "#CA0020"))+
  scale_color_manual(values = c("#F4A582", "#CA0020"))+
  scale_x_continuous(breaks = c(0, 0.5, 1.0))

cowplot::plot_grid(fig3.v1, nmbootfig, nrow = 2, rel_heights = c(2/3,1/3))



# ggpreview(plot = last_plot(), dpi = 650, width = 80, height = 120, units = "mm")
# ggsave(plot = last_plot(), filename = "Botsch_MG_Fig3.pdf", device = "pdf", dpi = 650, width = 80, height = 120, units = "mm")


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
  labs(x = "Initial Sediment Quality",
       y = "Proportion of Larval Midges\n3rd instar or higher",
       fill = "Number of Larvae")+
  theme(legend.position = "bottom")+
  scale_x_log10()


p3 %>% 
  filter(!is.na(algae_conc)) %>% 
  mutate(day = paste("Day", day)) %>% 
  ggplot(aes(y = prop, x = algae_conc2, fill = midge, size = live_tt))+
  geom_point(shape = 21, alpha = 0.5)+
  geom_hline(yintercept = startingprop)+
  facet_grid(midge~day)+
  lims(y = c(0,1))+
  labs(x = "Initial Sediment Quality",
       y = "Proportion of Larval Midges\n3rd instar or higher",
       size = "",
       fill = "")+
  theme(legend.position = "bottom")+
  scale_fill_manual(values = c("black", "gray60"), guide = NULL)+
  scale_x_log10()

# ggpreview(plot = last_plot(), dpi = 650, width = 80, height = 120, units = "mm")
# ggsave(plot = last_plot(), filename = "Botsch_MG_Fig4.pdf", device = "pdf", dpi = 650, width = 80, height = 120, units = "mm")


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

l3 %>% 
  filter(day %in% c(14, 22)) %>% 
  mutate(day = paste("Day", day)) %>% 
  ggplot(aes(x = algae_conc2, y = body_size, fill = midge))+
  facet_wrap(~day, ncol = 1)+
  geom_point(size  = 1.5, alpha = 0.7, shape = 21, stroke = 0.2, position = position_jitterdodge(dodge.width = 0.15, jitter.width = 0.1))+
  geom_smooth(aes(color = midge), method = "lm", size = 0.65, se = FALSE)+
  scale_fill_manual(values = c("black", "gray60"))+
  scale_color_manual(values = c("black", "gray60"))+
  labs(y = "Body Length of \nThird Instar Tanytarsini (mm)",
       x = "Initial Sediment Quality",
       color = "",
       fill = "")+
  scale_x_log10()+
  theme(legend.position = c(0.5,0.95),
        legend.text.align = 0,
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA, color = NA),
        legend.direction = "horizontal")

# ggpreview(plot = last_plot(), dpi = 650, width = 80, height = 80, units = "mm")
# ggsave(plot = last_plot(), filename = "Botsch_MG_Fig5.pdf", device = "pdf", dpi = 650, width = 80, height = 80, units = "mm")



##====Q2.3b Midge Length Alternate====
#3rd instars averaged
bl14log2 <- lm(body_size~log(algae_conc2)*midge+box,
               weights = live_tt,
                data = l314 %>% 
                  group_by(algae_conc2, midge, box) %>% 
                summarise(body_size = meanna(body_size)) %>% 
                 left_join(cc)) 

summary(bl14log2)
# plot(bl14log2)


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
# plot(bl14log3) 


bl22log3 <- lm(body_size~log(algae_conc2)*midge+box, 
               weights = live_tt,
               data = cm %>% 
                 filter(day == 22) %>% 
                 group_by(algae_conc2, midge, box) %>% 
                 summarise(body_size = meanna(body_size)) %>% 
                 left_join(cc))
summary(bl22log3)
# plot(bl22log3)

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

lmer(body_size~log(algae_conc2)*factor(day)+box+
       (1|coreid), 
   data = cm %>% 
     filter(midge == "Midges", instar  == 3)) %>% Anova(type = "3", test.statistic = "F")


lm(body_size~log(algae_conc2)*factor(day)+box, 
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
# plot(growth)

growth3 <- lm(body_size~day*factor(instar), data = cm %>% filter(species_name == "tt"))
summary(growth3)
plot(growth3)


cm %>% 
  filter(day!=0) %>% 
  group_by(day, algae_conc2) %>% 
  summarise(n = n(), 
            length = meanna(body_size), 
            length.se = sd(body_size, na.rm = TRUE)/sqrt(n), 
            avginstar = meanna(instar), 
            instar.se = sd(instar, na.rm = TRUE)/sqrt(n)) %>%
  bind_rows(cm %>% 
              filter(day == 0) %>% 
              summarise(n = n(), 
                        length = meanna(body_size), 
                        length.se = sd(body_size, na.rm = TRUE)/sqrt(n), 
                        avginstar = meanna(instar), 
                        instar.se = sd(instar, na.rm = TRUE)/sqrt(n)) %>% 
              crossing(day = 0, algae_conc2 = unique(meta$algae_conc2))) %>% 
  ggplot(aes(x = algae_conc2, y = length, ymin = length-length.se, ymax = length+length.se, col = factor(day)))+
  geom_pointrange(position = position_dodge(width = 0.1))+
  scale_x_log10()+
  labs(x = "Initial Sediment Quality",
       y = "Average Midge Length (mm)",
       color = "Day")


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
 

midgeprod2 <- cm %>% 
  left_join(cc) %>% 
  filter(species_name == "tt",
         day==22) %>% 
  group_by(coreid) %>% 
  mutate(w = weight(body_size)) %>% #estimate dry mass
  group_by(coreid, day, box, algae_conc) %>% 
  summarise(wt = meanna(w), #average biomass of a midge (mg)
            nt = unique(live_tt), #number of individuals in the mesocosm
            wt.se = sd(w, na.rm = TRUE)/sqrt(nt)) %>% 
  left_join(cm %>% 
              left_join(cc) %>% 
              filter(species_name == "tt",
                     day==14) %>% 
              group_by(coreid) %>% 
              mutate(w = weight(body_size)) %>% #estimate dry mass
              group_by(box, algae_conc) %>% 
              summarise(w0 = meanna(w), #average biomass of a midge (mg)
                        n0 = unique(live_tt), #number of individuals in the mesocosm
                        w0.se = sd(w, na.rm = TRUE)/sqrt(n0))) %>% 
  ungroup %>% 
  mutate(Pd = incsumprod(n1 = n0, n2 = nt, wt1 = w0, wt2 = wt, deltaT = 8)$Pd, #production in mg,
         g = incsumprod(n1 = n0, n2 = nt, wt1 = w0, wt2 = wt, deltaT = 8)$g,
         Pd.se = prodse(n1 = n0, n2 = nt, wt1 = w0, wt2 = wt, deltaT = 8, n1se = nt1.se, w1se = wt1.se, w2se = wt.se),
         Pdc.se = Pd.se*dwc*1000,
         Pdc = Pd*dwc*1000) #convert from weight to c (mg) then convert to micrograms 



#compare measurements of growth (lm) to these estimates
midgeprod %>% 
  group_by(day) %>% 
  mutate(meanw = mean((wt-wt1)/day)) %>% 
  ggplot(aes((wt-wt1)/day, fill = factor(day)))+
  geom_histogram(bins = 15, alpha = 0.7)+
  geom_vline(xintercept = 0, linetype = "dotted")+
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
  lm(Pdc~0+NPP:factor(AE), data = .) %>% 
  summary()



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
                 x = max(production$NPP)-5,
                 slope = AEm*NPEm,
                 label = paste("AE[m]==", format(AEm, nsmall = 1)),
                 label2 = "I[m]/P[a] == 1") %>% 
  mutate(y = x*NPEm*AEm)


proda <- production %>% 
  filter(!is.na(day), midge == "Midges") %>% 
  mutate(Pdc = ifelse(is.na(Pdc), 0, Pdc),
         Pdc.se = ifelse(is.na(Pdc), 0, Pdc.se)) %>% 
  ggplot(aes(x = NPP, y = Pdc))+
  # facet_wrap(~midge)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_errorbar(aes(ymin = Pdc-Pdc.se, ymax = Pdc+Pdc.se, color = algae_conc2))+
  geom_point(aes(color = algae_conc, shape = factor(day)), fill = "white")+
  geom_segment(aes(x = 0, y = 0, xend = max(production$NPP)-5, yend =  slope * (max(production$NPP)+5)), data = IPs)+
  geom_text(aes(x = x, y = y, label = label), parse = TRUE, hjust = 0, vjust = 1, color = "black", data = IPs, size = 3)+
  geom_text(aes(x = x, y = y, label = label2), parse = TRUE, hjust = 0,vjust = 0, color = "black", data = IPs, size = 3)+
  labs(x = expression(~P[a]~" \u03BCg C "~d^{-1}),
       y = expression(~P[m]~" \u03BCg C "~d^{-1}),
       color = "Initial\n Sediment Quality",
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

production %>% 
  filter(midge == "Midges") %>% 
  mutate(ingestion = Pdc/NPEm) %>% 
  lm(Pdc~NPP, data = .) %>% 
  summary()


prodb <- midgeprod2 %>% 
  mutate(gd = g*dwc*1000/8) %>% 
  left_join(algaeprod %>% 
              select(coreid, algae_conc2, midge, day, gpp_daily) %>% 
              group_by(coreid, algae_conc2) %>% 
              mutate(GPP = gpp_daily,
                     day = paste0("NPP_", day)) %>%
              select(-gpp_daily) %>% 
              spread(day, GPP) %>% 
              rowwise() %>% 
              mutate(GPP = meanna(c(NPP_14, NPP_22)), 
                     deltaGPP = NPP_22-NPP_14)) %>% 
  filter(!is.na(NPP_22)) %>% 
  ggplot(aes(x = gd, y = deltaGPP))+
  facet_wrap(~midge, scales = "free")+
  geom_hline(yintercept = 0, col = "gray80")+
  geom_vline(xintercept = 0, col = "gray80")+
  geom_point()+
  labs(y = expression("Change in GPP \u03BCg C "~d^{-1}),
       x = expression("Midge Growth \u03BCg C "~d^{-1}))


cowplot::plot_grid(proda, prodb, ncol = 1, rel_heights = c(2/3,1/3))
# ggpreview(plot = last_plot(), dpi = 650, width = 80, height = 120, units = "mm")
# ggsave(plot = last_plot(), filename = "Botsch_MG_Fig6.pdf", device = "pdf", dpi = 650, width = 80, height = 80, units = "mm")



midgeprod2 %>% 
  left_join(algaeprod %>% 
              select(coreid, algae_conc2, midge, day, gpp_daily) %>% 
              group_by(coreid, algae_conc2) %>% 
              mutate(GPP = gpp_daily,
                     day = paste0("NPP_", day)) %>%
              select(-gpp_daily) %>% 
              spread(day, GPP) %>% 
              rowwise() %>% 
              mutate(GPP = meanna(c(NPP_14, NPP_22)), 
                     deltaGPP = NPP_22-NPP_14)) %>% 
  filter(!is.na(NPP_22), midge == "No Midges") %>%
  lm(deltaGPP~Pdc, data = .) %>% 
  summary()

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
  filter(!is.na(NPP_22),
         midge == "Midges") %>%
  lm(deltaNPP~gd, data = .) %>% 
  summary()
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




