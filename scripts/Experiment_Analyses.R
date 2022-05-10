#====Load Packages====
library(tidyverse)
library(lubridate)
library(lme4)
library(car)
library(cowplot)
library(grid)
library(gridExtra)
source("scripts/MG_Functions.R")


#====read in files====

meta <- read_csv("clean data/MG_meta.csv") %>% 
  mutate(algae_conc2 = ifelse(algae_conc == 0, 2e-3, algae_conc)) #half of lowest value
chl <- read_csv("clean data/MG_chl.csv")
om <- read_csv("clean data/MG_om.csv")
ndvi <- read_csv("clean data/MG_ndvi.csv")
ep_raw <- read_csv("clean data/MG_ep.csv")
cm_raw <- read_csv("clean data/MG_cm.csv")
cc_raw <- read_csv("clean data/MG_cc.csv")


#====create data frame to predict model response variables ====

topredict <- data.frame(algae_conc2 = rep(unique(meta$algae_conc2), 4),
                        midge = rep(rep(c("Midges", "No Midges"), each = 10), 2),
                        box = rep(c("1", "2"), each = 20)) %>% 
  mutate(midge = fct_relevel(midge, "No Midges", "Midges"))

#model matrix used for all regressions
mm <- model.matrix(~log(algae_conc2)*midge + box, topredict) 
mm[,"box2"] <- 0.5
mm <- unique(mm)


#====Figure 1: Chorophyll and organic content in stock sediment====
omplot <- om %>% 
  left_join(meta) %>% 
  ggplot(aes(algae_conc2, perc.org))+
  geom_point(size =3)+
  scale_y_continuous(labels = scales::percent_format(),
                     limits = c(0, NA))+
  labs(x = "Sediment Treatment",
       y = "Organic Content")

#corrected for pheophytin
corrchlplot <- chl %>% 
  left_join(meta) %>% 
  ggplot(aes(x = algae_conc2, y = chl/1000))+
  geom_point(size = 2, alpha = 0.8)+
  scale_y_continuous(labels = scales::comma_format())+
  labs(x = "",
       y = expression("Chlorophyll mg"~L^{-1})) +  
  theme(axis.text.x = element_blank())

#uncorrected
uncorrchlplot <- chl %>% 
  left_join(meta) %>% 
  ggplot(aes(x = algae_conc2, y = uncorr_chl/1000))+
  geom_point(size = 2, alpha = 0.8)+
  scale_y_continuous(labels = scales::comma_format())+
  labs(x = "",
       y = "Chlorophyll mg/L") + 
  theme(axis.text.x = element_blank())


fig1 <- plot_grid(corrchlplot, omplot, nrow  = 2)
plot(fig1)

# ggpreview(plot = fig1, dpi = 650, width = 80, height = 100, units = "mm")
# ggsave(plot = fig1, filename = "Botsch_MG_Fig1.pdf", dpi = 650, device = "pdf", width = 80, height = 80, units = "mm")


#=====Question 1: Effects of Midges on Primary Producers=====
## Midge effect on production
nep <- ep_raw %>% 
  left_join(meta) %>% 
  mutate(coreid = as.character(coreid),
         box = as.character(box))

nep14 <- nep %>% 
  filter(day == 14) %>% 
  mutate(midge = fct_relevel(midge, c("No Midges", "Midges")))

g14log <- lm(log(gpp)~log(algae_conc2)*midge+box, data = nep14) 
summary(g14log) 
summary(update(g14log, .~. -log(algae_conc2):midge))
# plot(g14log)

#get predictions and standard errors
predicted14 <- mod_predict(model.matrix = mm, model = g14log) %>% 
  mutate(day = "Day 14")

nep22 <- nep %>% 
  filter(day == 22)  %>% 
  mutate(midge = fct_relevel(midge, c("No Midges", "Midges")))

g22log <- lm(log(gpp)~log(algae_conc2)*midge+box, data = nep22)
summary(g22log)
summary(update(g22log, .~. -log(algae_conc2):midge))
# plot(g22log)

predicted22 <- mod_predict(model.matrix = mm, model = g22log) %>% 
  mutate(day = "Day 22")

gpredict <- rbind(predicted14, predicted22) %>%
  mutate(gpp = exp(estimate),
         upper = exp(estimate+se),
         lower = exp(estimate-se))

#====Figure 2: GPP====
nepfig <- nep %>% 
  mutate(day = paste("Day", day)) %>% 
  ggplot(aes(x = algae_conc2, y = gpp, fill = midge))+
  facet_wrap(~day, ncol = 1)+
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = midge, linetype = midge), color = "black", size = 0.2, alpha = 0.5, data = gpredict, show.legend = FALSE)+
  geom_point(size  = 2, alpha = 0.75, shape = 21, position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.1, jitter.height = 0))+
  geom_line(aes(linetype = midge), size = 0.7, data = gpredict)+
  midge_color+
  midge_fill+ 
  midge_lines+
  labs(x = "Sediment Treatment",
       y = expression("GPP "~(g~O[2]~m^{-2}~hr^{-1})),
       color = "",
       fill = "",
       linetype = "")+
  scale_y_continuous(trans = "log", breaks = c(0.005, 0.015, 0.05))+
  scale_x_continuous(trans = "log", breaks = c(0.01, 0.1, 1))+
  theme(legend.position = c(0.5,0.95),
        legend.direction = "horizontal",
        legend.text.align = 0,
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA, color = NA))

plot(nepfig)
# ggpreview(plot = nepfig, dpi = 650, width = 80, height = 80, units = "mm")

nepfig2 <- nep %>% 
  mutate(day = paste("Day", day)) %>% 
  ggplot(aes(x = algae_conc2, y = gpp, fill = midge))+
  facet_wrap(~day, ncol = 2)+
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = midge, linetype = midge), color = "black", size = 0.2, alpha = 0.5, data = gpredict, show.legend = FALSE)+
  geom_point(size  = 2, alpha = 0.75, shape = 21, position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.1, jitter.height = 0))+
  geom_line(aes(linetype = midge), size = 0.7, data = gpredict)+
  midge_color+
  midge_fill+ 
  midge_lines+
  labs(x = "Sediment Treatment",
       y = expression("GPP "~(g~O[2]~m^{-2}~hr^{-1})),
       color = "",
       fill = "",
       linetype = "")+
  scale_y_continuous(trans = "log", breaks = c(0.005, 0.015, 0.05))+
  scale_x_continuous(trans = "log", breaks = c(0.01, 0.1, 1))+
  theme(legend.position = "none",
        legend.direction = "horizontal",
        legend.text.align = 0,
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA, color = NA))

#====Question 2: Sediment Effects on Midges =====
##====Q2.1: Number of Midges Present====
cc <- cc_raw %>% 
  left_join(meta) %>% 
  mutate(coreid = as.character(coreid),
         box = as.character(box))

cc14 <- cc %>% 
  filter(day == 14) %>% 
  mutate(midge = fct_relevel(midge, c("No Midges", "Midges")))


nm14log <- glm(live_tt~log(algae_conc2)*midge+box, 
               data = cc14,
               family = quasipoisson())

summary(nm14log)
# plot(nm14log)

summary(update(nm14log, .~.-log(algae_conc2):midge))



cc22 <- cc %>% 
  filter(day == 22) %>% 
  mutate(midge = fct_relevel(midge, c("No Midges", "Midges")))

nm22log <- glm(live_tt~log(algae_conc2)*midge+box, 
               data = cc22,
               family = quasipoisson())

summary(nm22log)
# plot(nm22log)

summary(update(nm22log, .~. -log(algae_conc2):midge))

#extract fits

predicted14 <- mod_predict(mm, nm14log) %>% 
  mutate(day = "Day 14")

predicted22 <- mod_predict(mm, nm22log) %>% 
  mutate(day = "Day 22")

nmpredict <- rbind(predicted14, predicted22) %>% 
  mutate(live_tt = exp(estimate),
         upper = exp(estimate+se),
         lower = exp(estimate-se))

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

nmidge_fig <- cc %>% 
  filter(day %in% c(14, 22)) %>% 
  mutate(day = paste("Day", day)) %>% 
  ggplot(aes(x = algae_conc2, y = live_tt, fill = midge))+
  facet_wrap(~day, ncol = 1)+
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = midge, linetype = midge), color = "black", size = 0.2, alpha = 0.5, data = nmpredict, show.legend = FALSE)+
  geom_point(size  = 2, alpha = 0.7, shape = 21, position = position_jitter(width = 0.1, height = 0))+
  geom_hline(yintercept = 20, linetype = "dashed", alpha = 0.5)+
  # geom_smooth(aes(color = midge), method = "glm", size = 0.75, se = FALSE, method.args = list(family = quasipoisson()) )+
  geom_line(aes(linetype = midge), data = nmpredict)+
  midge_color+
  midge_fill+ 
  labs(x = "Sediment Treatment",
       y = "Number of Larvae",
       color = "",
       fill = "",
       linetype = "")+
  scale_x_continuous(trans = "log", breaks = c(0.01, 0.1, 1))+
  theme(legend.position = c(0.5,0.95),
        legend.text.align = 0,
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA, color = NA),
        legend.direction = "horizontal")

plot(nmidge_fig)
# ggpreview(plot = nmidge_fig, dpi = 650, width = 80, height = 80, units = "mm")


nm_fig2 <- cc %>% 
  filter(day %in% c(14, 22)) %>% 
  mutate(day = paste("Day", day)) %>% 
  ggplot(aes(x = algae_conc2, y = live_tt, fill = midge))+
  facet_wrap(~day, ncol = 2)+
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = midge, linetype = midge), color = "black", size = 0.2, alpha = 0.5, data = nmpredict, show.legend = FALSE)+
  geom_point(size  = 2, alpha = 0.7, shape = 21, position = position_jitter(width = 0.1, height = 0))+
  geom_hline(yintercept = 20, linetype = "dashed", alpha = 0.5)+
  # geom_smooth(aes(color = midge), method = "glm", size = 0.75, se = FALSE, method.args = list(family = quasipoisson()) )+
  geom_line(aes(linetype = midge), size = 0.7, data = nmpredict)+
  midge_color+
  midge_fill+ 
  labs(x = "Sediment Treatment",
       y = "Number of Larvae",
       color = "",
       fill = "",
       linetype = "")+
  scale_x_continuous(trans = "log", breaks = c(0.01, 0.1, 1))+
  theme(legend.position = "none",
        legend.text.align = 0,
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA, color = NA),
        legend.direction = "horizontal")

##====Q2.2: Midge Development (SUPPLEMENT)====
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
  mutate(midge = factor(midge, levels = c("No Midges", "Midges")))

prop14log <- glm(cbind(s2, others) ~ log(algae_conc2)*midge+box, 
                 data = p314, 
                 family = quasibinomial())

summary(prop14log)
# plot(prop14log)
summary(update(prop14log, .~. -log(algae_conc2):midge))


#subset for day 22
p322 <- p2 %>% 
  filter(day == 22) %>% 
  mutate(midge = factor(midge, levels = c("No Midges", "Midges")))

prop22log <- glm(cbind(s2, others) ~log(algae_conc2)*midge+box, 
                 data = p322, 
                 family = quasibinomial())
summary(prop22log)
# plot(prop22log)
summary(update(prop22log, .~. -log(algae_conc2):midge))



predicted14 <- mod_predict(mm, prop14log) %>% 
  mutate(day = "Day 14")

predicted22 <- mod_predict(mm, prop22log) %>% 
  mutate(day = "Day 22")

propredict <- rbind(predicted14, predicted22) %>% 
  mutate(prop = exp(estimate),
         upper = exp(estimate+se),
         lower = exp(estimate-se))



#====Figure 4: Instar====

prop_fig <- p2 %>% 
  mutate(day = paste("Day", day)) %>% 
  ggplot(aes(y = prop, x = algae_conc2))+
  geom_hline(yintercept = startingprop, alpha = 0.75, color = "orange")+
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = midge, linetype = midge), color = "black", size = 0.2, alpha = 0.5, data = propredict, show.legend = FALSE)+
  geom_point(aes(fill = midge), shape = 21, alpha = 0.75, size = 2, position = position_jitterdodge(dodge.width = 0.1, jitter.width = 0.1, jitter.height = 0))+
  geom_line(aes(linetype = midge), data = propredict)+
  facet_wrap(~day, nrow = 2)+
  labs(x = "Sediment Treatment",
       y = "Proportion of Larvae in 2nd Instar",
       color = NULL,
       fill = NULL,
       linetype = NULL)+
  midge_color+
  midge_fill+ 
  midge_lines+
  theme(legend.position = "bottom")+
  scale_x_continuous(trans = "log", breaks = c(0.01, 0.1, 1))+
  theme(legend.position = c(0.8,0.35),
        legend.text.align = 0,
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA, color = NA),
        legend.direction = "vertical")

prop_fig
ggpreview(plot = prop_fig, dpi = 650, width = 80, height = 80, units = "mm")
# ggsave(plot = last_plot(), filename = "Botsch_MG_Fig4.pdf", device = "pdf", dpi = 650, width = 80, height = 120, units = "mm")

#====Q2.3: Midge Length====
#subset day 14
l314 <- cm %>% 
  filter(day == 14) %>% 
  mutate(midge = factor(midge, levels = c("No Midges", "Midges")),
         instar = factor(instar))


bl14log <- lmer(body_size~log(algae_conc2)*midge+box+
                  (1|coreid), 
                data = l314) 

summary(bl14log)
# plot(bl14log)

Anova(bl14log, type = "3", test.statistic = "F")

summary(update(bl14log, .~.-log(algae_conc2):midge))
Anova(update(bl14log, .~.-log(algae_conc2):midge), type = "2", test.statistic = "F")

#subset day 22
l322 <- cm %>% 
  filter(day == 22) %>% 
  mutate(midge = factor(midge, levels = c("No Midges", "Midges")),
         instar = factor(instar))


bl22log <- lmer(body_size~log(algae_conc2)*midge+box +
                  (1|coreid), 
                data = l322)
summary(bl22log)
# plot(bl22log)

#F Tests
Anova(bl22log, type = "3", test.statistic = "F") #includes interaction 

summary(update(bl22log, .~. -log(algae_conc2):midge))
Anova(update(bl22log, .~. -log(algae_conc2):midge), type = "2", test.statistic = "F") #drops interaction


predicted14 <- mod_predict(mm, bl14log, mixed.effects.model = TRUE) %>% 
  mutate(day = "Day 14")

predicted22 <- mod_predict(mm, bl22log, mixed.effects.model = TRUE) %>% 
  mutate(day = "Day 22")

blpredict <- rbind(predicted14, predicted22) %>% 
  mutate(body_size = estimate, 
         lower = estimate-se,
         upper = estimate+se)


#====Figure 5: Body Length====

#average body legnth in initial midges
start_length <- {cm %>% 
    filter(day == 0) %>% 
    summarise(length = meanna(body_size))}$length

blfig <- cm %>% 
  filter(day %in% c(14, 22)) %>% 
  mutate(day = paste("Day", day)) %>% 
  ggplot(aes(x = algae_conc2, y = body_size))+
  facet_wrap(~day, ncol = 1)+
  geom_hline(yintercept = start_length, linetype = 2, size = 0.5)+
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = midge, linetype = midge), color = "black", size = 0.2, alpha = 0.5, data = blpredict, show.legend = FALSE)+
  geom_point(aes(fill = midge), size  = 1.5, alpha = 0.6, shape = 21, stroke = 0.2, position = position_jitterdodge(dodge.width = 0.4, jitter.width = 0.2, jitter.height = 0))+
  geom_line(aes(linetype = midge), size = 0.7, data = blpredict)+
  midge_color+
  midge_fill+ 
  midge_lines+
  labs(y = "Body Length (mm)",
       x = "Sediment Treatment",
       color = "",
       fill = "",
       linetype = "")+
  scale_x_continuous(trans = "log", breaks = c(0.01, 0.1, 1))+
  theme(legend.position = c(0.5,0.95),
        legend.text.align = 0,
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA, color = NA),
        legend.direction = "horizontal")

plot(blfig)
# ggpreview(plot = blfig, dpi = 650, width = 80, height = 80, units = "mm")



blfig2 <- cm %>% 
  filter(day %in% c(14, 22),
         !is.na(body_size)) %>% 
  mutate(day = paste("Day", day)) %>% 
  ggplot(aes(x = algae_conc2, y = body_size))+
  facet_wrap(~day, ncol = 2)+
  geom_hline(yintercept = start_length, linetype = 2, size = 0.5, alpha = 0.5)+
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = midge, linetype = midge), color = "black", size = 0.2, alpha = 0.5, data = blpredict, show.legend = FALSE)+
  geom_point(aes(fill = midge), size  = 1.5, alpha = 0.6, shape = 21, stroke = 0.2, position = position_jitterdodge(dodge.width = 0.4, jitter.width = 0.2, jitter.height = 0))+
  geom_line(aes(linetype = midge), size = 0.7, data = blpredict)+
  midge_color+
  midge_fill+ 
  midge_lines+
  labs(y = "Body Length (mm)",
       x = "Sediment Treatment",
       color = "",
       fill = "",
       linetype = "")+
  scale_x_continuous(trans = "log", breaks = c(0.01, 0.1, 1))+
  theme(legend.position = "none",
        legend.text.align = 0,
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA, color = NA),
        legend.direction = "horizontal")

# ggpreview(plot = blfig, dpi = 650, width = 80, height = 80, units = "mm")
# ggsave(plot = last_plot(), filename = "Botsch_MG_Fig5.pdf", device = "pdf", dpi = 650, width = 80, height = 80, units = "mm")



#====Combine All Plots into one figure====

figcomb <- plot_grid(nepfig2 + theme(axis.title.x = element_blank()), 
                     nm_fig2+ theme(axis.title.x = element_blank(), axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 15))), 
                     blfig2+theme(axis.title.x = element_blank(), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 15))), ncol = 1, labels = c("A", "B", "C"))

# ggpreview(figcomb, width = 6, height = 6, units = "in", dpi = 650)

leg <- get_legend(nepfig+theme(legend.direction = "horizontal", legend.position = "top"))

figcomb2 <- plot_grid(leg, figcomb, ncol = 1, rel_heights = c(0.05, 0.95))


x.grob <- textGrob("Sediment Treatment")

figcomb3 <- grid.arrange(figcomb2, bottom = x.grob)


# ggpreview(figcomb3, width = 5, height = 6, units = "in", dpi = 650)
