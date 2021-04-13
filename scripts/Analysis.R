#====Load Packages====
library(tidyverse)
library(lubridate)
library(viridis)
library(gridExtra)
library(lme4)
library(car)
library(broom)

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
                  legend.position = "bottom",
                  text = element_text(size = 10),
                  axis.title = element_text(size = 10.5)))

#====read in files====

meta <- read_csv("clean data/MG_meta.csv") %>% 
  mutate(algae_conc2 = ifelse(algae_conc == 0, 2e-3, algae_conc)) #check
chl <- read_csv("clean data/MG_chl.csv")
om <- read_csv("clean data/MG_om.csv")
ndvi <- read_csv("clean data/MG_ndvi.csv")
ep_raw <- read_csv("clean data/MG_ep.csv")
cm_raw <- read_csv("clean data/MG_cm.csv")
cc_raw <- read_csv("clean data/MG_cc.csv")


#====create data frame to predict model response variables ====

topredict <- data.frame(algae_conc2 = rep(unique(meta$algae_conc2), 4),
                        midge = rep(rep(c("Midges", "No Midges"), each = 10), 2),
                        box = rep(c("1", "2"), each = 20))

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
  filter(day == 14) %>% 
  mutate(midge = fct_relevel(midge, c("No Midges", "Midges")))

g14log <- lm(gpp~log(algae_conc2)*midge+box, data = nep14) 
summary(g14log) 
# plot(g14log)


nep22 <- nep %>% 
  filter(day == 22) %>% 
  mutate(midge = fct_relevel(midge, c("No Midges", "Midges")))

g22log <- lm(gpp~log(algae_conc2)*midge+box, data = nep22)
summary(g22log)
# plot(g22log)


predicted14 <- topredict %>% 
  mutate(gpp = predict(g14log, ., type = "response")) %>% 
  group_by(algae_conc2, midge) %>% 
  summarise(gpp = mean(gpp)) %>% 
  mutate(day = "Day 14")

predicted22 <- topredict %>% 
  mutate(gpp = predict(g22log, ., type = "response")) %>% 
  group_by(algae_conc2, midge) %>% 
  summarise(gpp = mean(gpp)) %>% 
  mutate(day = "Day 22")

gpredict <- rbind(predicted14, predicted22)



#====Figure 2: GPP====
nepv2.1 <- nep %>% 
  mutate(day = paste("Day", day)) %>% 
  ggplot(aes(x = algae_conc2, y = gpp, fill = midge))+
  facet_wrap(~day, ncol = 1)+
  geom_point(size  = 2, alpha = 0.75, shape = 21)+
  geom_line(aes(color = midge), size = 0.7, data = gpredict)+
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


#estimate for effect of initial food availability on gpp day 14
axn_14 <- as.numeric(g14log$coefficients[5] + g14log$coefficients[2])

#standard error of day 14 estimate
axn_14se <- sqrt(vcov(g14log)[2,2] + vcov(g14log)[5,5] + 2*vcov(g14log)[2,5])

#print
c(estimate = axn_14, se = axn_14se)

axn_22 <- as.numeric(g22log$coefficients[5] + g22log$coefficients[2])

axn_22se <- sqrt(vcov(g22log)[2,2] + vcov(g22log)[5,5] + 2*vcov(g22log)[2,5])

c(estimate = axn_22, se = axn_22se)


gppbootfig <- data.frame(estimate = c(axn_14, axn_22, g14log$coefficients["log(algae_conc2)"], g22log$coefficients["log(algae_conc2)"]), 
           se = c(axn_14se, axn_22se, coef(summary(g14log))["log(algae_conc2)", "Std. Error"], coef(summary(g22log))["log(algae_conc2)", "Std. Error"]), 
           day = c("Day 14", "Day 22", "Day 14", "Day 22"), 
           midge = c("No Midges", "No Midges", "Midges", "Midges")) %>% 
  mutate(midge = fct_relevel(midge, c("No Midges", "Midges"))) %>% 
  ggplot(aes(x = estimate, xmin = estimate-se, xmax = estimate+se, midge, col = day))+
  geom_vline(xintercept = 0)+
  geom_pointrange(position = position_dodge(width = 0.4))+
  geom_text(aes(label = day, x = estimate-se), hjust = 0, vjust = 2, data = . %>% filter(midge == "No Midges"), size = (5/14)*10*0.8, family = "", fontface = "plain", position = position_dodge(width = 0.4))+
  scale_color_manual(values = c("#F4A582", "#CA0020"))+
  labs(y = "", 
       x = "Effect of Initial Resource Availability")+
  theme(legend.position = "none")





cowplot::plot_grid(nepv2.1, gppbootfig, nrow = 2, rel_heights = c(2/3,1/3))

# ggpreview(plot = last_plot(), dpi = 650, width = 80, height = 120, units = "mm")
# ggsave(plot = last_plot(), filename = "Botsch_MG_Fig2.pdf", device = "pdf", dpi = 650, width = 80, height = 120, units = "mm")


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

nep14log <- lm(nep~ log(algae_conc2)*midge + box, data = nep14)
summary(nep14log)

nep22log <- lm(nep~ log(algae_conc2)*midge + box, data = nep22)
summary(nep22log)

nep %>% 
  mutate(day = paste("Day", day,  sep = " ")) %>% 
  ggplot(aes(x = -resp, y = gpp, fill = algae_conc2))+
  geom_point(shape = 21, alpha = 0.5)+
  facet_grid(midge~day)+
  geom_abline(slope = 1)+
  coord_equal()+
  scale_fill_viridis_c(trans = "log", breaks = c(0.01, 0.1, 1))+
  labs(x = "Respiration",
       y = "GPP",
       fill = "Initial Resource Availability")


#====Question 2: Sediment Effects on Midges =====
##====Q2.1: Number of Midges Present====
cc <- cc_raw %>% 
  left_join(meta) %>% 
  mutate(coreid = as.character(coreid),
         box = as.character(box))

cc14 <- cc %>% 
  filter(day == 14) %>% 
  mutate(midge = fct_relevel(midge, c("Midges", "No Midges")))


nm14log <- glm(live_tt~log(algae_conc2)*midge+box, 
               data = cc14,
               family = quasipoisson())

summary(nm14log)
# plot(nm14log)

cc22 <- cc %>% 
  filter(day == 22) %>% 
  mutate(midge = fct_relevel(midge, c("Midges", "No Midges")))

nm22log <- glm(live_tt~log(algae_conc2)*midge+box, 
               data = cc22,
               family = quasipoisson())

summary(nm22log)
# plot(nm22log)


#extract fits

predicted14 <- topredict %>% 
  mutate(fit = predict(nm14log, ., type = "link")) %>% 
  group_by(algae_conc2, midge) %>% 
  summarise(live_tt = exp(mean(fit))) %>% 
  mutate(day = "Day 14")

predicted22 <- topredict %>% 
  mutate(fit = predict(nm22log, ., type = "link")) %>% 
  group_by(algae_conc2, midge) %>% 
  summarise(live_tt = exp(mean(fit))) %>% 
  mutate(day = "Day 22")

nmpredict <- rbind(predicted14, predicted22)

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

fig3.v1 <- cc %>% 
  filter(day %in% c(14, 22)) %>% 
  mutate(day = paste("Day", day)) %>% 
  ggplot(aes(x = algae_conc2, y = live_tt, fill = midge))+
  facet_wrap(~day, ncol = 1)+
  geom_point(size  = 2, alpha = 0.7, shape = 21)+
  geom_hline(yintercept = 20, linetype = "dashed", alpha = 0.5)+
  # geom_smooth(aes(color = midge), method = "glm", size = 0.75, se = FALSE, method.args = list(family = quasipoisson()) )+
  geom_line(aes(color = midge), data = nmpredict)+
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


#estimate for effect of initial food availability on gpp day 14
nm_axn_14 <- as.numeric(nm14log$coefficients[5] + nm14log$coefficients[2])

#standard error of day 14 estimate
nm_axn_14se <- sqrt(vcov(nm14log)[2,2] + vcov(nm14log)[5,5] + 2*vcov(nm14log)[2,5])

#print
c(estimate = nm_axn_14, se = nm_axn_14se)

nm_axn_22 <- as.numeric(nm22log$coefficients[5] + nm22log$coefficients[2])

nm_axn_22se <- sqrt(vcov(nm22log)[2,2] + vcov(nm22log)[5,5] + 2*vcov(nm22log)[2,5])

c(estimate = nm_axn_22, se = nm_axn_22se)


nmbootfig <- data.frame(estimate = c(nm_axn_14, nm_axn_22, nm14log$coefficients["log(algae_conc2)"], nm22log$coefficients["log(algae_conc2)"]), 
                         se = c(nm_axn_14se, nm_axn_22se, coef(summary(nm14log))["log(algae_conc2)", "Std. Error"], coef(summary(nm22log))["log(algae_conc2)", "Std. Error"]), 
                         day = c("Day 14", "Day 22", "Day 14", "Day 22"), 
                         midge = c("No Midges", "No Midges", "Midges", "Midges")) %>% 
  mutate(midge = fct_relevel(midge, c("No Midges", "Midges"))) %>% 
  ggplot(aes(x = estimate, xmin = estimate-se, xmax = estimate+se, midge, col = day))+
  geom_vline(xintercept = 0)+
  geom_pointrange(position = position_dodge(width = 0.4))+
  geom_text(aes(label = day, x = estimate-se), hjust = 0, vjust = 2, data = . %>% filter(midge == "No Midges"), size = (5/14)*10*0.8, family = "", fontface = "plain", position = position_dodge(width = 0.4))+
  scale_color_manual(values = c("#F4A582", "#CA0020"))+
  labs(y = "", 
       x = "Effect of Initial Resource Availability")+
  theme(legend.position = "none")


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
    filter(species_name == "tt",
           day == 0) %>%  #remove the four non-tanytarsini
    add_count(coreid, name = "measured") %>% 
    group_by(coreid, day, midge, algae_conc2, box, instar, measured) %>% 
    count(instar, name = "s2") %>% 
    filter(instar == 2) %>% 
    mutate(prop = s2/measured)}$prop


p2 <- cm %>% 
  filter(species_name == "tt",
         day!= 0) %>%  #remove the four non-tanytarsini
  add_count(coreid, name = "measured") %>% 
  group_by(coreid, day, midge, algae_conc2, box, instar, measured) %>% 
  count(instar, name = "count") %>% 
  mutate(instar = paste0("s", instar)) %>% 
  spread(instar, count, fill = 0) %>% 
  mutate(prop = s2/measured,
         others = s3+s4)

#subset for day 14
p314 <- p2 %>% 
  filter(day == 14) %>% 
  mutate(midge = factor(midge, levels = c("Midges", "No Midges")))

prop14log <- glm(cbind(s2, others) ~log(algae_conc2)*midge+box, 
                 data = p314, 
                 family = quasibinomial())

summary(prop14log)
# plot(prop14log)



#subset for day 22
p322 <- p2 %>% 
  filter(day == 22) %>% 
  mutate(midge = factor(midge, levels = c("Midges", "No Midges")))

prop22log <- glm(cbind(s2, others) ~log(algae_conc2)*midge+box, 
                 data = p322, 
                 family = quasibinomial())
summary(prop22log)
# plot(prop22log)

predicted14 <- topredict %>% 
  mutate(fit = predict(prop14log, ., type = "link")) %>% 
  group_by(algae_conc2, midge) %>% 
  summarise(prop = exp(mean(fit))) %>% 
  mutate(day = "Day 14")

predicted22 <- topredict %>% 
  mutate(fit = predict(prop22log, ., type = "link")) %>% 
  group_by(algae_conc2, midge) %>% 
  summarise(prop = exp(mean(fit))) %>% 
  mutate(day = "Day 22")

propredict <- rbind(predicted14, predicted22)



#====Figure 4: Instar====

fig4 <- p2 %>% 
  mutate(day = paste("Day", day)) %>% 
  ggplot(aes(y = prop, x = algae_conc2))+
  geom_point(aes(fill = midge), shape = 21, alpha = 0.75, size = 2, position = position_jitter(width = 0.15, height = 0))+
  geom_line(aes(color = midge), data = propredict)+
  geom_hline(yintercept = startingprop, alpha = 0.75)+
  facet_wrap(~day, nrow = 2)+
  lims(y = c(0,1))+
  labs(x = "Initial Resource Availability",
       y = "Proportion of Larval Midges in 2nd Instar",
       color = NULL,
       fill = NULL)+
  scale_color_manual(values = c("black", "gray60"))+
  scale_fill_manual(values = c("black", "gray60"))+
  theme(legend.position = "bottom")+
  scale_x_continuous(trans = "log", breaks = c(0.01, 0.1, 1))+
  theme(legend.position = c(0.20,0.95),
        legend.text.align = 0,
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA, color = NA),
        legend.direction = "vertical")

fig4
ggpreview(plot = last_plot(), dpi = 650, width = 80, height = 120, units = "mm")
# ggsave(plot = last_plot(), filename = "Botsch_MG_Fig4.pdf", device = "pdf", dpi = 650, width = 80, height = 120, units = "mm")


#estimate for effect of initial food availability on gpp day 14
prop_axn_14 <- as.numeric(prop14log$coefficients[5] + prop14log$coefficients[2])

#standard error of day 14 estimate
prop_axn_14se <- sqrt(vcov(prop14log)[2,2] + vcov(prop14log)[5,5] + 2*vcov(prop14log)[2,5])

#print
c(estimate = prop_axn_14, se = prop_axn_14se)

prop_axn_22 <- as.numeric(prop22log$coefficients[5] + prop22log$coefficients[2])

prop_axn_22se <- sqrt(vcov(prop22log)[2,2] + vcov(prop22log)[5,5] + 2*vcov(prop22log)[2,5])

c(estimate = prop_axn_22, se = prop_axn_22se)


propbootfig <- data.frame(estimate = c(prop_axn_14, prop_axn_22, prop14log$coefficients["log(algae_conc2)"], prop22log$coefficients["log(algae_conc2)"]), 
                        se = c(prop_axn_14se, prop_axn_22se, coef(summary(prop14log))["log(algae_conc2)", "Std. Error"], coef(summary(prop22log))["log(algae_conc2)", "Std. Error"]), 
                        day = c("Day 14", "Day 22", "Day 14", "Day 22"), 
                        midge = c("No Midges", "No Midges", "Midges", "Midges")) %>% 
  mutate(midge = fct_relevel(midge, c("No Midges", "Midges"))) %>% 
  ggplot(aes(x = estimate, xmin = estimate-se, xmax = estimate+se, midge, col = day))+
  geom_vline(xintercept = 0)+
  geom_pointrange(position = position_dodge(width = 0.4))+
  geom_text(aes(label = day, x = estimate-se), hjust = 0, vjust = 2, data = . %>% filter(midge == "No Midges"), size = (5/14)*10*0.8, family = "", fontface = "plain", position = position_dodge(width = 0.4))+
  scale_color_manual(values = c("#F4A582", "#CA0020"))+
  labs(y = "", 
       x = "Effect of Initial Resource Availability")+
  theme(legend.position = "none")

cowplot::plot_grid(fig4, propbootfig, nrow = 2, rel_heights = c(2/3,1/3))
# ggpreview(plot = last_plot(), dpi = 650, width = 80, height = 120, units = "mm")


##====Q2.3: Midge Length====
#subset to not include the stocked midges
l3 <- cm %>% 
  filter(day!=0)

#subset day 14
l314 <- cm %>% 
  filter(day == 14) %>% 
  mutate(midge = factor(midge, levels = c("Midges", "No Midges")))


bl14log <- lmer(body_size~log(algae_conc2)*midge+box+
                  (1|coreid), 
                data = l314) 

summary(bl14log)
# plot(bl14log)

Anova(bl14log, type = "3", test.statistic = "F")
Anova(bl14log, type = "2", test.statistic = "F")

#subset day 22
l322 <- l3 %>% 
  filter(day == 22) %>% 
  mutate(midge = factor(midge, levels = c("Midges", "No Midges")))


bl22log <- lmer(body_size~log(algae_conc2)*midge+box +
                  (1|coreid), 
                data = l322)
summary(bl22log)
# plot(bl22log)

Anova(bl22log, type = "3", test.statistic = "F") 
Anova(bl22log, type = "2", test.statistic = "F")


predicted14 <- topredict %>% 
  mutate(fit = predict(bl14log, newdata = ., re.form = NA, type = "response")) %>% 
  group_by(algae_conc2, midge) %>% 
  summarise(body_size = mean(fit)) %>% 
  mutate(day = "Day 14")

predicted22 <- topredict %>% 
  mutate(fit = predict(bl22log, newdata = ., re.form = NA, type = "response")) %>% 
  group_by(algae_conc2, midge) %>% 
  summarise(body_size = mean(fit)) %>% 
  mutate(day = "Day 22")

blpredict <- rbind(predicted14, predicted22)


#====Figure 5: Body Length====

#average body legnth in initial midges
start_length <- {cm %>% 
  filter(day == 0) %>% 
  summarise(length = meanna(body_size))}$length

blfig <- l3 %>% 
  filter(day %in% c(14, 22)) %>% 
  mutate(day = paste("Day", day)) %>% 
  ggplot(aes(x = algae_conc2, y = body_size))+
  facet_wrap(~day, ncol = 1)+
  geom_hline(yintercept = start_length, linetype = 2, size = 0.5)+
  geom_point(aes(fill = midge), size  = 1.5, alpha = 0.7, shape = 21, stroke = 0.2, position = position_jitterdodge(dodge.width = 0.15, jitter.width = 0.1))+
  geom_line(aes(col = midge), size = 0.7, data = blpredict)+
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


#estimate for effect of initial food availability on gpp day 14
blaxn_14 <- as.numeric(fixef(bl14log)[5] + fixef(bl14log)[2])

#standard error of day 14 estimate
blaxn_14se <- sqrt(vcov(bl14log)[2,2] + vcov(bl14log)[5,5] + 2*vcov(bl14log)[2,5])

#print
c(estimate = blaxn_14, se = blaxn_14se)

blaxn_22 <- as.numeric(fixef(bl22log)[5] + fixef(bl22log)[2])

blaxn_22se <- sqrt(vcov(bl22log)[2,2] + vcov(bl22log)[5,5] + 2*vcov(bl22log)[2,5])

c(estimate = blaxn_22, se = blaxn_22se)


bl_slopes <- data.frame(estimate = c(blaxn_14, blaxn_22, fixef(bl14log)["log(algae_conc2)"], fixef(bl22log)["log(algae_conc2)"]), 
                         se = c(blaxn_14se, blaxn_22se, coef(summary(bl14log))["log(algae_conc2)", "Std. Error"], coef(summary(bl22log))["log(algae_conc2)", "Std. Error"]), 
                         day = c("Day 14", "Day 22", "Day 14", "Day 22"), 
                         midge = c("No Midges", "No Midges", "Midges", "Midges")) %>% 
  mutate(midge = fct_relevel(midge, c("No Midges", "Midges"))) %>% 
  ggplot(aes(x = estimate, xmin = estimate-se, xmax = estimate+se, midge, col = day))+
  geom_vline(xintercept = 0)+
  geom_pointrange(position = position_dodge(width = 0.4))+
  geom_text(aes(label = day, x = estimate-se), hjust = 0, vjust = 2, data = . %>% filter(midge == "No Midges"), size = (5/14)*10*0.8, family = "", fontface = "plain", position = position_dodge(width = 0.4))+
  scale_color_manual(values = c("#F4A582", "#CA0020"))+
  labs(y = "", 
       x = "Effect of Initial Resource Availability")+
  theme(legend.position = "none")

cowplot::plot_grid(blfig, bl_slopes, nrow = 2, rel_heights = c(2/3,1/3))


# ggpreview(plot = last_plot(), dpi = 650, width = 80, height = 100, units = "mm")
# ggsave(plot = last_plot(), filename = "Botsch_MG_Fig5.pdf", device = "pdf", dpi = 650, width = 80, height = 80, units = "mm")





#====Print all Model outputs to a csv======
#extract LMs (LMERs more complicated)
# tidy(g14log) %>% mutate(response = "GPP", day = "14") %>%
#   rbind(tidy(g22log) %>% mutate(response = "GPP", day = "22")) %>%
#   rbind(tidy(nm14log) %>% mutate(response = "# Midges", day = "14")) %>%
#   rbind(tidy(nm22log) %>% mutate(response = "# Midges", day = "22")) %>%
#   rbind(tidy(prop14log) %>% mutate(response = "proportion 2nd instar", day = "14")) %>%
#   rbind(tidy(prop22log) %>% mutate(response = "proportion 2nd instar", day = "22")) %>%
#   mutate(estimate = format(estimate, digits = 2),
#          std.error = format(std.error, digits = 2),
#          statistic = round(statistic, digits = 2),
#          p.value = format(p.value, digits = 2)) %>%
#   write_csv(. , "lmoutput.csv")

  

#===============================================

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
