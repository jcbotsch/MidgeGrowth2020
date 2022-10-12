#====Load Packages====
library(tidyverse)
library(lubridate)
library(lme4)
library(car)
library(cowplot)
library(grid)
library(gridExtra)
library(shades)
source("scripts/MG_Functions.R")

#====read in files====
meta <- read_csv("clean data/MG_meta.csv") %>% 
  mutate(algae_conc2 = ifelse(algae_conc == 0, 2e-3, algae_conc)) #half of lowest value
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

#=====Question 1: Effects of Midges on Primary Producers=====
## prep data
nep <- ep_raw %>% 
  left_join(meta) %>% 
  mutate(coreid = as.character(coreid),
         box = as.character(box),
         gpp = gpp*1000) #convert GPP from gm-2h-1 to gcm-2h-1

#extract day 14 NEPs
nep14 <- nep %>% 
  filter(day == 14) %>% 
  mutate(midge = fct_relevel(midge, c("No Midges", "Midges")))

#fit model
g14log <- lm(gpp~log(algae_conc2)*midge+box, data = nep14) 
summary(g14log) #type III
summary(update(g14log, .~. -log(algae_conc2):midge)) #type II

#extract day 22 NEPs
nep22 <- nep %>% 
  filter(day == 22)  %>% 
  mutate(midge = fct_relevel(midge, c("No Midges", "Midges")))

#fit model
g22log <- lm(gpp~log(algae_conc2)*midge+box, data = nep22)
summary(g22log) #type III
summary(update(g22log, .~. -log(algae_conc2):midge)) #type II

#get predictions and standard errors
predicted14 <- mod_predict(model.matrix = mm, model = g14log) %>% 
  mutate(day = "Day 14")
predicted22 <- mod_predict(model.matrix = mm, model = g22log) %>% 
  mutate(day = "Day 22")

#join predictions and standard errors
gpredict <- rbind(predicted14, predicted22) %>%
  mutate(gpp = (estimate),
         upper = (estimate+se),
         lower = (estimate-se))

#====Figure 1A====
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
  labs(x = "Initial Algal Abundance",
       y = expression("GPP "~(mg~O[2]~m^{-2}~hr^{-1})),
       color = "",
       fill = "",
       linetype = "")+
  scale_x_continuous(trans = "log", breaks = c(0.01, 0.1, 1))+
  theme(legend.position = "none")

#====Question 2: Sediment Effects on Midges =====
##====Q2.1: Number of Midges Present====
#prep data
cc <- cc_raw %>% 
  left_join(meta) %>% 
  mutate(coreid = as.character(coreid),
         box = as.character(box))

#how many midges out of the total stocked (1,000) survived?
cc %>% 
  filter(!coreid %in% c(108, 119, 152)) %>% 
  summarise(sum(live_tt))

cc %>% 
  filter(coreid %in% c(108, 119, 152)) %>% 
  summarise(stocked = sum(live_tt)/3) %>% 
  mutate(stocked.tot = stocked*50)

646/883

#extract day 14 midge counts
cc14 <- cc %>% 
  filter(day == 14) %>% 
  mutate(midge = fct_relevel(midge, c("No Midges", "Midges")))

#fit model
nm14log <- glm(live_tt~log(algae_conc2)*midge+box, 
               data = cc14,
               family = quasipoisson())

summary(nm14log) #type III 
summary(update(nm14log, .~.-log(algae_conc2):midge)) #type II


#extract day 22 midge counts
cc22 <- cc %>% 
  filter(day == 22) %>% 
  mutate(midge = fct_relevel(midge, c("No Midges", "Midges")))

#fit model
nm22log <- glm(live_tt~log(algae_conc2)*midge+box, 
               data = cc22,
               family = quasipoisson())

summary(nm22log) #type III
summary(update(nm22log, .~. -log(algae_conc2):midge)) #type II

# get preditions
predicted14 <- mod_predict(mm, nm14log) %>% 
  mutate(day = "Day 14")
predicted22 <- mod_predict(mm, nm22log) %>% 
  mutate(day = "Day 22")

nmpredict <- rbind(predicted14, predicted22) %>% 
  mutate(live_tt = exp(estimate),
         upper = exp(estimate+se),
         lower = exp(estimate-se))

#====Figure 1B: Number of Midges====
nm_fig2 <- cc %>% 
  filter(day %in% c(14, 22)) %>% 
  mutate(day = paste("Day", day)) %>% 
  ggplot(aes(x = algae_conc2, y = live_tt, fill = midge))+
  facet_wrap(~day, ncol = 2)+
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = midge, linetype = midge), color = "black", size = 0.2, alpha = 0.5, data = nmpredict, show.legend = FALSE)+
  geom_point(size  = 2, alpha = 0.7, shape = 21, position = position_jitter(width = 0.1, height = 0))+
  geom_hline(yintercept = 20, linetype = "dashed", alpha = 0.5)+
  geom_line(aes(linetype = midge), size = 0.7, data = nmpredict)+
  midge_color+
  midge_fill+ 
  labs(x = "Initial Algal Abundance",
       y = "Number of Larvae",
       color = "",
       fill = "",
       linetype = "")+
  scale_x_continuous(trans = "log", breaks = c(0.01, 0.1, 1))+
  theme(legend.position = "none")

#====Q2.2: Midge Length====
#prep data
cm <- cm_raw %>% 
  left_join(meta) %>% 
  mutate(coreid = as.character(coreid),
         box = as.character(box),
         body_size = body_size/1000) #convert from um to mm

#subset day 14
l314 <- cm %>% 
  filter(day == 14) %>% 
  mutate(midge = factor(midge, levels = c("No Midges", "Midges")),
         instar = factor(instar))

#fit model
bl14log <- lmer(body_size~log(algae_conc2)*midge+box+
                  (1|coreid), 
                data = l314) 
#singular fit because coreid explains no variation.

summary(bl14log) #coefficients 
Anova(bl14log, type = "3", test.statistic = "F") # Type III
Anova(update(bl14log, .~.-log(algae_conc2):midge), type = "2", test.statistic = "F") # Type II

#subset day 22
l322 <- cm %>% 
  filter(day == 22) %>% 
  mutate(midge = factor(midge, levels = c("No Midges", "Midges")),
         instar = factor(instar))

#fit model
bl22log <- lmer(body_size~log(algae_conc2)*midge+box +
                  (1|coreid), 
                data = l322)

summary(bl22log) #coefficients
Anova(bl22log, type = "3", test.statistic = "F") #Type III
Anova(update(bl22log, .~. -log(algae_conc2):midge), type = "2", test.statistic = "F") #Type II

#extract predictions
predicted14 <- mod_predict(mm, bl14log, mixed.effects.model = TRUE) %>% 
  mutate(day = "Day 14")
predicted22 <- mod_predict(mm, bl22log, mixed.effects.model = TRUE) %>% 
  mutate(day = "Day 22")

#combine predictions
blpredict <- rbind(predicted14, predicted22) %>% 
  mutate(body_size = estimate, 
         lower = estimate-se,
         upper = estimate+se)

#====Figure 1C: Body Length====
#average body legnth in initial midges
start_length <- {cm %>% 
    filter(day == 0) %>% 
    summarise(length = meanna(body_size))}$length

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
       x = "Initial Algal Abundance",
       color = "",
       fill = "",
       linetype = "")+
  scale_x_continuous(trans = "log", breaks = c(0.01, 0.1, 1))+
  theme(legend.position = "none")

#====Combine Figure 1====
figcomb <- plot_grid(nepfig2 + theme(axis.title.x = element_blank()), 
                     nm_fig2+ theme(axis.title.x = element_blank()), 
                     blfig2+theme(axis.title.x = element_blank()), ncol = 1, labels = c("A", "B", "C"),
                     align = "v")

leg <- get_legend(nepfig2+theme(legend.direction = "horizontal", legend.position = "top"))

figcomb2 <- plot_grid(leg, figcomb, ncol = 1, rel_heights = c(0.05, 1))


x.grob <- textGrob("Initial Algal Abundance", gp = gpar(fontsize = 11))

figcomb3 <- grid.arrange(figcomb2, bottom = x.grob)

ggpreview(figcomb3, width = 5, height = 6, units = "in", dpi = 800)
# ggsave(figcomb3, filename = "figures/Botsch_MidgeGrowth_fig1.pdf", device = "pdf", width = 5, height = 6, units = "in", dpi = 800)



