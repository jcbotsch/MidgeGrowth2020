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
  ggplot(aes(algae_conc2, perc.org))+
  geom_point(size =3)+
  scale_y_continuous(labels = scales::percent_format())+
  labs(x = "Resource Availability",
       y = "Organic Content")

#corrected for pheophytin
corrchlplot <- chl %>% 
  left_join(meta) %>% 
  ggplot(aes(x = algae_conc2, y = chl/1000))+
  geom_point(size = 2, alpha = 0.8)+
  scale_y_continuous(labels = scales::comma_format())+
  labs(x = "",
       y = "Chlorophyll mg/L") +  #\u03BC is the unicode for mu
  theme(axis.text.x = element_blank())

#uncorrected
uncorrchlplot <- chl %>% 
  left_join(meta) %>% 
  ggplot(aes(x = algae_conc2, y = uncorr_chl/1000))+
  geom_point(size = 2, alpha = 0.8)+
  scale_y_continuous(labels = scales::comma_format())+
  labs(x = "",
       y = "Chlorophyll mg/L") + #\u03BC is the unicode for mu
  theme(axis.text.x = element_blank())


fig1 <- grid.arrange(corrchlplot, omplot, nrow  = 2)
plot(fig1)

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

# #peform bootstrapping
# gboot1 <- para.boot(data1 = nep14, data2 = nep22,
#           response.var = "gpp",
#           lm1 = g14log, lm2 = g22log,
#           results = gboot1, nboot = nboot)

#estimate for effect of initial food availability on gpp day 14
axn_14 <- as.numeric(g14log$coefficients[5] + g14log$coefficients[2])

#standard error of day 14 estimate
axn_14se <- sqrt(vcov(g14log)[2,2] + vcov(g14log)[5,5] + 2*vcov(g14log)[2,5])

#print
c(estimate = axn_14, se = axn_14se)
# exttract the same info from the bootstrap
gboot1 %>% filter(day == "Day 14") %>% summarise(mean = mean(algae_x_NoMidge),
                                                             sd = sd(algae_x_NoMidge))

axn_22 <- as.numeric(g22log$coefficients[5] + g22log$coefficients[2])

axn_22se <- sqrt(vcov(g22log)[2,2] + vcov(g22log)[5,5] + 2*vcov(g22log)[2,5])

c(estimate = axn_22, se = axn_22se)


gppbootfig <- data.frame(estimate = c(axn_14, axn_22, g14log$coefficients["log(algae_conc2)"], g22log$coefficients["log(algae_conc2)"]), 
           se = c(axn_14se, axn_22se, coef(summary(g14log))["log(algae_conc2)", "Std. Error"], coef(summary(g22log))["log(algae_conc2)", "Std. Error"]), 
           day = c("Day 14", "Day 22", "Day 14", "Day 22"), 
           midge = c("No Midges", "No Midges", "Midges", "Midges")) %>% 
  mutate(midge = fct_relevel(midge, c("No Midges", "Midges"))) %>% 
  ggplot(aes(x = estimate, xmin = estimate-2*se, xmax = estimate+2*se, midge, col = day))+
  geom_vline(xintercept = 0)+
  geom_pointrange(position = position_dodge(width = 0.4))+
  geom_text(aes(label = day, x = estimate-2*se), hjust = 0, data = . %>% filter(midge == "No Midges"))+
  scale_color_manual(values = c("#F4A582", "#CA0020"))+
  labs(y = "", 
       x = "Effect of Initial Resource Availability")+
  theme(legend.position = "none")





gboot1 %>% filter(day == "Day 22") %>% summarise(mean = mean(algae_x_NoMidge), sd = sd(algae_x_NoMidge))

ggplot(data.frame(algae_x_NoMidge = c(-0.005, 0.005)), aes(x = algae_x_NoMidge))+
  geom_density(data = gboot1 %>% filter(day == "Day 14"), fill = "blue")+
  stat_function(fun = dnorm, args = list(axn_14, axn_14se), col = "red", size = 2)

ggplot(data.frame(algae_x_NoMidge = c(-0.005, 0.005)), aes(x = algae_x_NoMidge))+
  geom_density(data = gboot1 %>% filter(day == "Day 22"), fill = "blue")+
  stat_function(fun = dnorm, args = list(axn_22, axn_22se), col = "red", size = 2)

ggplot(data.frame(x = c(-0.0035, 0.004)), aes(x = x))+
  stat_function(geom = "area", aes(fill = "Day 22", color = "No Midges"), fun = dnorm, args = list(axn_22, axn_22se), size = 1, alpha = 0.5)+
  stat_function(geom = "area",aes(fill = "Day 14", color = "No Midges"), fun = dnorm, args = list(axn_14, axn_14se),  size = 1, alpha = 0.5)+
  stat_function(geom = "area",aes(fill = "Day 22", color = "Midges"), fun = dnorm, args = list(g14log$coefficients[2], coef(summary(g14log))["log(algae_conc2)", "Std. Error"]), size = 1, alpha = 0.5)+
  stat_function(geom = "area",aes(fill = "Day 14", color = "Midges"), fun = dnorm, args = list(coef(g22log)["log(algae_conc2)"], coef(summary(g22log))["log(algae_conc2)", "Std. Error"]),  size = 1, alpha = 0.5)+
  scale_color_manual(values = c("gray", "black"))+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank())


ggplot(data.frame(x = c(-0.0035, 0.004)), aes(x = x))+
  stat_function(geom = "area", aes(fill = "Day 22", color = "No Midges"), fun = dnorm, args = list(axn_22, axn_22se), size = 1, alpha = 0.5)+
  stat_function(geom = "area",aes(fill = "Day 14", color = "No Midges"), fun = dnorm, args = list(axn_14, axn_14se),  size = 1, alpha = 0.5)+
  stat_function(geom = "area",aes(fill = "Day 22", color = "Midges"), fun = dnorm, args = list(g14log$coefficients[2], coef(summary(g14log))["log(algae_conc2)", "Std. Error"]), size = 1, alpha = 0.5)+
  stat_function(geom = "area",aes(fill = "Day 14", color = "Midges"), fun = dnorm, args = list(coef(g22log)["log(algae_conc2)"], coef(summary(g22log))["log(algae_conc2)", "Std. Error"]),  size = 1, alpha = 0.5)+
  scale_color_manual(values = c("gray", "black"))+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank())

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
  labs(x = "Initial Resource Availability",
       y = expression("GPP "~(g~O[2]~m^{-2}~hr^{-1})),
       color = "",
       fill = "")+
  scale_x_continuous(trans = "log", breaks = c(0.01, 0.1, 1))+
  theme(legend.position = c(0.2,0.95),
        legend.text.align = 0,
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA, color = NA))

#get densities
gdens <- data.frame(x = seq(-0.005, 0.005, by = 0.0001)) %>% 
  mutate(density = dnorm(x, axn_14, axn_14se),
         label.x = axn_14,
         day = "Day 14",
         midge = "No Midges") %>% 
  bind_rows(data.frame(x = seq(-0.005, 0.005, by = 0.0001)) %>% 
              mutate(density = dnorm(x, axn_22, axn_22se),
                     label.x = axn_22,
                     day = "Day 22",
                     midge = "No Midges")) %>% 
  bind_rows(data.frame(x = seq(-0.005, 0.005, by = 0.0001)) %>% 
              mutate(density = dnorm(x, g14log$coefficients["log(algae_conc2)"], coef(summary(g14log))["log(algae_conc2)", "Std. Error"]),
                     label.x = g14log$coefficients["log(algae_conc2)"],
                     day = "Day 14",
                     midge = "Midges")) %>% 
  bind_rows(data.frame(x = seq(-0.005, 0.005, by = 0.0001)) %>% 
              mutate(density = dnorm(x, g22log$coefficients["log(algae_conc2)"], coef(summary(g22log))["log(algae_conc2)", "Std. Error"]),
                     label.x = g22log$coefficients["log(algae_conc2)"],
                     day = "Day 22",
                     midge = "Midges")) %>% 
  group_by(x) %>% 
  mutate(total_density = sum(density)) %>% 
  filter(total_density>1) %>% 
  group_by(day, midge) %>% 
  mutate(label.y = max(density))
  ungroup



gppbootfig <- gdens %>% 
  ggplot(aes(x, y = density, fill = day))+
  geom_area(position = "identity", alpha = 0.5, color = "gray50")+
  # geom_density(alpha = 0.5, color = "gray50", data =. %>% filter(day=="Day 22"))+
  geom_vline(xintercept = 0, linetype = "dashed")+
  facet_wrap(~midge)+
  geom_text(aes(label = day, y = label.y, x = label.x, color = day))+
  labs(x = "Effect of Initial Resource Availability on GPP",
       y = "",
       fill = "")+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")+
  scale_fill_manual(values = c("#F4A582", "#CA0020"))+
  scale_color_manual(values = c("#F4A582", "#CA0020"))+
  scale_x_continuous(breaks = c(0.003,0, -0.003))

gdens %>% 
  ggplot()+
  ggridges::geom_density_ridges(aes(x = x, y = midge, fill = day, height = density, group = interaction(midge, day)), stat = "identity", alpha = 0.5)+
  geom_vline(xintercept = 0, linetype = "dashed")+
  labs(x = "Effect of Initial Resource Availability on GPP",
       y = "",
       fill = "")+
  theme(legend.position = c(0.1, 0.8))+
  scale_fill_manual(values = c("#F4A582", "#CA0020"))+
  scale_color_manual(values = c("#F4A582", "#CA0020"))+
  scale_x_continuous(breaks = c(0.003,0, -0.003))

cowplot::plot_grid(nepv2.1, gppbootfig, nrow = 2, rel_heights = c(2/3,1/3))

ggpreview(plot = last_plot(), dpi = 650, width = 80, height = 120, units = "mm")
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
       color = "Initial Resource Availability")+
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
  labs(x = "Initial Resource Availability",
       y = expression(g~O[2]~m^{-2}~hr^{-1}),
       color = "",
       fill = "")+
  scale_x_continuous(trans = "log", breaks = c(0.01, 0.1, 1))+
  theme(strip.placement = "outside")

r14log <- lm(resp~ log(algae_conc2)*midge + box, data = nep14)
summary(r14log)

r22log <- lm(resp~ log(algae_conc2)*midge + box, data = nep22)
summary(r22log)

nep %>% 
  ggplot(aes(x = -resp, y = gpp, col = algae_conc2))+
  geom_point()+
  facet_wrap(~day)+
  geom_abline(slope = 1)+
  coord_equal()


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
       x = "Initial Resource Availability",
       fill = expression(GPP(g~O[2]~m^{-2}~hr^{-1})~" "))+
  scale_x_continuous(trans = "log", breaks = c(0.01, 0.1, 1))

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
  labs(x = "Initial Resource Availability",
       y = "Live Tanytarsini",
       color = "",
       fill = "")+
  scale_x_continuous(trans = "log", breaks = c(0.01, 0.1, 1))+
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
  labs(x = "Initial Resource Availability",
       y = "Live Tanytarsini",
       color = "",
       fill = "")+
  lims(y = c(0,NA))+
  scale_x_continuous(trans = "log", breaks = c(0.01, 0.1, 1))+
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
       x = "Initial Resource Availability",
       fill = expression(GPP(g~O[2]~m^{-2}~hr^{-1})~" "))+
  scale_x_continuous(trans = "log", breaks = c(0.01, 0.1, 1))+
  lims(y = c(0, NA))





cc %>%
  ungroup %>% 
  filter(as.numeric(coreid)<=100) %>% 
  arrange(as.numeric(coreid)) %>% 
  mutate(column = rep(1:5, times = 20),
         row = rep(1:20, each = 5),
         diff = ifelse(midge == "Midges", live_tt-20, live_tt)) %>% 
  ggplot(aes(x = row, y = column, shape = midge, col = midge))+
  facet_grid(day~box, scales = "free")+
  geom_tile(aes(fill = algae_conc2), col = "white")+
  geom_label(aes(label = live_tt), size = 5, alpha = 1)+
  scale_shape_manual(values = c(16,18), )+
  scale_fill_viridis(trans = "log10")+
  scale_color_manual(values = c("firebrick4", "dodgerblue"))+
  scale_x_continuous(breaks = c(1:20))+
  theme(axis.title = element_blank())


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
  labs(x = "Effect of Initial Resource Availability \non Number of Midges",
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
  labs(x = "Effect of Initial Resource Availability \non Number of Midges",
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
  labs(x = "Initial Resource Availability",
       y = "Proportion of Larval Midges\n3rd instar or higher",
       fill = "Number of Larvae")+
  theme(legend.position = "bottom")+
  scale_x_continuous(trans = "log", breaks = c(0.01, 0.1, 1))


p3 %>% 
  filter(!is.na(algae_conc)) %>% 
  mutate(day = paste("Day", day)) %>% 
  ggplot(aes(y = prop, x = algae_conc2, fill = midge, size = live_tt))+
  geom_point(shape = 21, alpha = 0.5)+
  geom_hline(yintercept = startingprop)+
  facet_grid(midge~day)+
  lims(y = c(0,1))+
  labs(x = "Initial Resource Availability",
       y = "Proportion of Larval Midges\n3rd instar or higher",
       size = "",
       fill = "")+
  theme(legend.position = "bottom")+
  scale_fill_manual(values = c("black", "gray60"), guide = NULL)+
  scale_x_continuous(trans = "log", breaks = c(0.01, 0.1, 1))

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
#subset to not include the stocked midges
l3 <- cm %>% 
  filter(day!=0)

#subset day 14
l314 <- cm %>% 
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

#average body legnth in initial midges
start_length <- {cm %>% 
  filter(day == 0) %>% 
  summarise(length = meanna(body_size))}$length

l3 %>% 
  filter(day %in% c(14, 22)) %>% 
  mutate(day = paste("Day", day)) %>% 
  ggplot(aes(x = algae_conc2, y = body_size))+
  facet_wrap(~day, ncol = 1)+
  geom_hline(yintercept = start_length, linetype = 2, size = 0.5)+
  geom_point(aes(fill = midge), size  = 1.5, alpha = 0.7, shape = 21, stroke = 0.2, position = position_jitterdodge(dodge.width = 0.15, jitter.width = 0.1))+
  geom_smooth(aes(color = midge), method = "lm", size = 0.65, se = FALSE)+
  scale_fill_manual(values = c("black", "gray60"))+
  scale_color_manual(values = c("black", "gray60"))+
  labs(y = "Midge Body Length (mm)",
       x = "Initial Resource Availability",
       color = "",
       fill = "")+
  scale_x_continuous(trans = "log", breaks = c(0.01, 0.1, 1))+
  theme(legend.position = c(0.5,0.95),
        legend.text.align = 0,
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA, color = NA),
        legend.direction = "horizontal")

ggpreview(plot = last_plot(), dpi = 650, width = 80, height = 80, units = "mm")
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
            body_size = meanna(body_size),
            n = n()) %>% 
  mutate(day = paste("Day", day)) %>% 
  ggplot(aes(x = algae_conc2, y = body_size, fill = midge))+
  facet_wrap(~day, scales = "free_x")+
  geom_point(aes(size = n), alpha = 0.7, shape = 21, position = position_dodge())+
  geom_smooth(aes(col = midge ), method = "lm", se = FALSE, size = 0.7)+
  labs(x = "Initial Sediment Conditions",
       y = "Body Length of \nThird Instar Tanytarsini (mm)",
       fill = "")+
  theme(legend.position = "bottom")+
  scale_x_continuous(trans = "log", breaks = c(0.01, 0.1, 1))

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
       y = "Midge Body Length (mm)",
       fill = expression(GPP~(g^{2}~O[2]~m^{-2}~hr^{-1})~" "))+
  theme(legend.position = "bottom")+
  scale_x_continuous(trans = "log", breaks = c(0.01, 0.1, 1))

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
# plot(growth3)


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
  scale_x_continuous(trans = "log", breaks = c(0.01, 0.1, 1))+
  labs(x = "Initial Resource Availability",
       y = "Average Midge Length (mm)",
       color = "Day")+
  scale_color_brewer(palette = "Set1")

#=====Supplement: Looking at Spatial patterns within boxes====


#====Spatial Patterns of GPP====
nep %>%
  filter(day == 14) %>% 
  ungroup %>% 
  filter(as.numeric(coreid)<=100) %>% 
  arrange(as.numeric(coreid)) %>% 
  mutate(column = rep(1:5, times = 20),
         row = rep(1:20, each = 5)) %>% 
  ggplot(aes(x = row, y = column)) +
  facet_wrap(~box, scales = "free")+
  geom_point(aes(fill = algae_conc2, size = gpp, shape = midge), col = "black")+
  scale_shape_manual(values = c(21,22))+
  scale_color_manual(values = c("firebrick4", "dodgerblue"))+
  scale_x_continuous(breaks = c(1:20))+
  scale_fill_viridis_c(trans = "log10")+
  theme(axis.title = element_blank())

# ggpreview(plot = last_plot(), dpi  =1000, width = 10, height= 5, units = "in")



nep %>%
  filter(day == 22) %>% 
  ungroup %>% 
  filter(as.numeric(coreid)<=100) %>% 
  left_join(nep %>%
              filter(day == 14) %>% 
              ungroup %>% 
              filter(as.numeric(coreid)<=100) %>% 
              arrange(as.numeric(coreid)) %>% 
              mutate(column = rep(1:5, times = 20),
                     row = rep(1:20, each = 5)) %>% 
              select(coreid, column, row)) %>% 
  
  ggplot(aes(x = row, y = column)) +
  facet_wrap(~box, scales = "free")+
  geom_point(aes(fill = algae_conc2, size = gpp, shape = midge), col = "black")+
  scale_shape_manual(values = c(21,22))+
  scale_color_manual(values = c("firebrick4", "dodgerblue"))+
  scale_x_continuous(breaks = c(1:20))+
  scale_fill_viridis_c(trans = "log10")+
  theme(axis.title = element_blank())

#====Movement patterns====
cc %>%
  ungroup %>% 
  filter(as.numeric(coreid)<=100) %>% 
  arrange(as.numeric(coreid)) %>% 
  mutate(column = rep(1:5, times = 20),
         row = rep(1:20, each = 5),
         diff = ifelse(midge == "Midges", live_tt-18, live_tt),
         diff = ifelse(diff>0, paste0("+", diff), diff)) %>% 
  ggplot(aes(x = row, y = column, shape = midge, col = midge))+
  facet_wrap(~box, scales = "free")+
  geom_tile(aes(fill = algae_conc2), col = "white")+
  geom_label(aes(label = diff), size = 5, alpha = 1)+
  scale_shape_manual(values = c(16,18), )+
  scale_fill_viridis(trans = "log10")+
  scale_color_manual(values = c("firebrick4", "dodgerblue"))+
  scale_x_continuous(breaks = c(1:20))+
  theme(axis.title = element_blank())

# ggpreview(plot = last_plot(), dpi  =1000, width = 10, height= 5, units = "in")
