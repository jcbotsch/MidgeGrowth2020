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

#Lindegaard et al (1979) formula for converting length to biomass (tab.21)
#AFDW in mg^(1/3) = 0.0442 + 0.0879 * length in mm 
weight <- function(l){
  (0.0442 + 0.0879 * l )^3 
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
  mutate(algae_conc2 = ifelse(algae_conc == 0, 2e-3, algae_conc), #half of lowest value
         box = factor(box),
         midge = factor(midge, levels = c("No Midges", "Midges"))) 
chl <- read_csv("clean data/MG_chl.csv")
om <- read_csv("clean data/MG_om.csv")
ndvi <- read_csv("clean data/MG_ndvi.csv")
ep_raw <- read_csv("clean data/MG_ep.csv")
cm_raw <- read_csv("clean data/MG_cm.csv")
cc_raw <- read_csv("clean data/MG_cc.csv")


#====create data frame to predict model response variables ====
topredict <- crossing(algae_conc2 = unique(meta$algae_conc2),
         midge = c("Midges", "No Midges"),
         box = c("1", "2"),
         day = c("14", "22")) %>% 
  mutate(midge = factor(midge, levels = c("No Midges", "Midges")))

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



#=====Question 1: Effects of Midges on Primary Producers=====
## Midge effect on production
nep <- ep_raw %>% 
  left_join(meta) %>% 
  mutate(coreid = as.character(coreid),
         day = factor(day))

nep2 <- nep %>% left_join(chl %>% group_by(algae) %>% summarise(chl = mean(chl)))


nep2 %>% 
  gather(metric, value, chl, algae_conc2) %>% 
  ggplot(aes(y = gpp, x = log(value), color = midge))+
  facet_wrap(metric~day, scales = "free")+
  geom_smooth(method = "lm")+
  geom_point()

nep14 <- nep %>% 
  filter(day == 14) %>% 
  mutate(midge = fct_relevel(midge, c("No Midges", "Midges")))

g14log <- lm(gpp~log(algae_conc2)*midge+box, data = nep14) 
summary(g14log) 

summary(update(g14log, .~. -log(algae_conc2):midge)) 

# plot(g14log)


nep22 <- nep %>% 
  filter(day == 22) %>% 
  mutate(midge = fct_relevel(midge, c("No Midges", "Midges")))

g22log <- lm(gpp~log(algae_conc2)*midge+box, data = nep22)
summary(g22log)
plot(g22log)

gpredict <- topredict %>% 
  mutate(gpp = if_else(day == "14", 
                       predict(g14log, ., type = "response"), 
                       predict(g22log, ., type = "response")),
         day = paste("Day", day, sep = " ")) %>% 
  group_by(algae_conc2, midge, day) %>% 
  summarise(gpp = mean(gpp)) 



glog <- lmer(log(gpp)~log(algae_conc2)*day*midge+box + (1|coreid), data = nep)

qqnorm(residuals(glog)); qqline(residuals(glog))
plot(glog)

summary(glog)
Anova(glog, test.statistic = "F", type = "3")
Anova(update(glog, .~. -log(algae_conc2):day:midge), test.statistic = "F", type = "3")
Anova(update(glog, .~. -log(algae_conc2):day:midge - log(algae_conc2):day - log(algae_conc2):midge -day:midge), test.statistic = "F", type = "2")


#====Figure 2: GPP====
 nep %>% 
  mutate(day = paste("Day", day)) %>% 
  ggplot(aes(x = algae_conc2, y = gpp, fill = midge))+
  facet_wrap(~day, ncol = 1)+
  geom_point(size  = 2, alpha = 0.75, shape = 21)+
  # geom_line(aes(color = midge), size = 0.7, data = gpredict)+
  geom_smooth(method = "lm")+
  scale_fill_manual(values = c("gray60", "black"))+
  scale_color_manual(values = c("gray60", "black"))+
  labs(x = "Sediment Treatment",
       y = expression("GPP "~(g~O[2]~m^{-2}~hr^{-1})),
       color = "",
       fill = "")+
  scale_y_continuous(trans = "log", breaks = c(0.005, 0.01, 0.02))+
  scale_x_continuous(trans = "log", breaks = c(0.01, 0.1, 1))+
  theme(legend.position = c(0.5,0.95),
        legend.direction = "horizontal",
        legend.text.align = 0,
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA, color = NA))



ggpreview(plot = nepv2, dpi = 650, width = 80, height = 80, units = "mm")
# ggsave(plot = last_plot(), filename = "Botsch_MG_Fig2.pdf", device = "pdf", dpi = 650, width = 80, height = 120, units = "mm")


#Supplemental Figure
nep %>% 
  mutate(day = paste("Day", day)) %>% 
  gather(metabolism, value, resp, nep, gpp) %>% 
  mutate(metabolism = toupper(metabolism),
         metabolism = ifelse(metabolism == "RESP", "ER", metabolism),
         metabolism = fct_reorder(metabolism, -value),
         value = ifelse(metabolism == "ER", -value, value)) %>% 
  ggplot(aes(x = algae_conc2, y = value, fill = midge))+
  geom_hline(yintercept = 0, linetype = "dotted")+
  facet_grid(metabolism~day, scale = "free_y", switch = "y")+
  geom_point(size  = 2, alpha = 0.7, shape = 21)+
  geom_smooth(aes(color = midge), method = "lm", size = 0.75, se = FALSE)+
  scale_fill_manual(values = c("black", "gray60"))+
  scale_color_manual(values = c("black", "gray60"))+
  labs(x = "Sediment Treatment",
       y = expression(g~O[2]~m^{-2}~hr^{-1}),
       color = "",
       fill = "")+
  scale_x_continuous(trans = "log", breaks = c(0.01, 0.1, 1))+
  theme(strip.placement = "outside")

# ggpreview(plot = last_plot(), dpi = 650, width = 120, height = 120, units = "mm")

#====Question 2: Sediment Effects on Midges =====
##====Q2.1: Survival ====
cc <- cc_raw %>% 
  left_join(meta) %>% 
  mutate(coreid = as.character(coreid),
         day = factor(day))

nmlog <- glm(live_tt~log(algae_conc2)*midge*day+box, 
             data = cc,
             family = quasipoisson())

summary(nmlog)
Anova(nmlog, test.statistic = "F", type = "3")
Anova(update(nmlog, .~. - log(algae_conc2):midge:day), test.statistic = "F", type = "3")
Anova(update(nmlog, .~. -log(algae_conc2):midge:day - log(algae_conc2):midge -log(algae_conc2):day - day:midge))


#extract fits

nmpredict <- topredict %>% 
  mutate(live_tt = predict(nmlog, ., type = "response"),
         day = paste("Day", day, sep = " ")) %>% 
  group_by(algae_conc2, midge, day) %>% 
  summarise(live_tt = mean(live_tt)) 


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
  labs(x = "Sediment Treatment",
       y = "Number of Larvae",
       color = "",
       fill = "")+
  scale_x_continuous(trans = "log", breaks = c(0.01, 0.1, 1))+
  theme(legend.position = c(0.5,0.95),
        legend.text.align = 0,
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA, color = NA),
        legend.direction = "horizontal")

#====General Info====
#how many midges out of the total stocked (1,000) survived?

cc %>% 
  filter(day!=0) %>% 
  summarise(count = sum(live_tt))

cc %>% 
  filter(coreid %in% c(108, 119, 152)) %>% 
  summarise(stocked = sum(live_tt)/3) %>% 
  mutate(stocked.tot = stocked*50)

646/883


cc %>% 
  group_by(day) %>% 
  summarise(count = sum(live_tt)) %>% 
  mutate(midges_added = ifelse(day == 14, 20*20,
                               ifelse(day == 22, 30*20, 3*20)),
         proportion = count/midges_added) 


##====Q2: Midge Growth====
cm <- cm_raw %>% 
  left_join(meta) %>% 
  mutate(coreid = as.character(coreid),
         day = factor(day),
         body_size = body_size/1000, #convert from um to mm
         wt = weight(body_size)) 

#===
llog <- lmer(body_size~log(algae_conc2)*day*midge+box + (1|coreid), data = cm %>% filter(day!=0)) 
summary(llog)
Anova(llog, test.statistic = "F", type = "3")
Anova(update(llog, .~. -log(algae_conc2):day:midge), test.statistic = "F", type = "3")
Anova(update(llog, .~. -log(algae_conc2):day:midge - log(algae_conc2):midge - log(algae_conc2):day - midge:day), test.statistic = "F", type = "2")

lmlog <- lmer(body_size~log(algae_conc2)*day+box + (1|coreid), data = cm %>% filter(day!=0, midge == "Midges")) 

summary(lmlog)
Anova(lmlog, test.satistic = "F", type = "3")
Anova(update(lmlog, .~. - log(algae_conc2):day), test.satistic = "F", type = "3")


lnmlog <- lmer(body_size~log(algae_conc2)*day+box + (1|coreid), data = cm %>% filter(day!=0, midge == "No Midges")) 

summary(lnmlog)
Anova(lnmlog, test.satistic = "F", type = "3")
Anova(update(lnmlog, .~. - log(algae_conc2):day), test.satistic = "F", type = "3")
#===

l14log <- lmer(body_size~log(algae_conc2)*midge+box + (1|coreid), data = cm %>% filter(day==14)) 

summary(l14log)
Anova(l14log, test.satistic = "F", type = "3")
Anova(update(l14log, .~. - log(algae_conc2):midge), test.satistic = "F", type = "3")


l22log <- lmer(body_size~log(algae_conc2)*midge+box + (1|coreid), data = cm %>% filter(day==22)) 

summary(l22log)
Anova(l22log, test.satistic = "F", type = "3")
Anova(update(l22log, .~. - log(algae_conc2):midge), test.satistic = "F", type = "3")

#===
wtlog <- lmer(wt~log(algae_conc2)*day*midge+box + (1|coreid), data = cm %>% filter(day!=0)) 
summary(llog)
Anova(wtlog, test.statistic = "F", type = "3")
Anova(update(wtlog, .~. -log(algae_conc2):day:midge), test.statistic = "F", type = "3")
Anova(update(wtlog, .~. -log(algae_conc2):day:midge - log(algae_conc2):midge - log(algae_conc2):day - midge:day), test.statistic = "F", type = "2")


wtlog <- lmer(wt~log(algae_conc2)*day+box + (1|coreid), data = cm %>% filter(day!=0, midge == "Midges")) 
summary(wtmlog)
Anova(wtmlog, type = "3", test.statistic = "F")
Anova(update(wtmlog, .~.-log(algae_conc2):day), type = "2", test.statistic = "F")

wtlog <- lmer(wt~log(algae_conc2)*day+box + (1|coreid), data = cm %>% filter(day!=0, midge == "No Midges")) 
summary(wtnmlog)
Anova(wtnmlog, type = "3", test.statistic = "F")
Anova(update(wtnmlog, .~.-log(algae_conc2):day), type = "2", test.statistic = "F")


plot(wtmlog)

qqnorm(residuals(wtmlog));qqline(residuals(wtmlog))

lm(wt~log(algae_conc2)*day*midge+box, 
   weight = n,
   data = cm %>% 
     group_by(coreid, midge, day, algae_conc2, box) %>% 
     summarise(n = n(), wt = meanna(wt)) %>% mutate(midge = fct_relevel(midge, c("Midges", "No Midges")))) %>% summary() 


summary(wtlog)
Anova(wtlog, test.statistic = "F", type = "3")
Anova(update(wtlog,.~. -log(algae_conc2):day:midge), test.statistic = "F", type = "3")
Anova(update(wtlog,.~. -log(algae_conc2):day:midge - log(algae_conc2):day - log(algae_conc2):midge - day:midge), test.statistic = "F", type = "2")

anova(wtlog, 
      update(wtlog,.~. -log(algae_conc2):day:midge),
      update(wtlog,.~. -log(algae_conc2):day:midge - log(algae_conc2):day - log(algae_conc2):midge - day:midge), 
      test = "F")

wtpredict <- topredict %>% 
  mutate(fit = predict(wtlog, newdata = ., re.form = NA, type = "response")) %>% 
  group_by(day, algae_conc2, midge) %>% 
  summarise(wt = mean(fit)) 

#===Plot Growth====

cm4plot <- cm %>% 
  filter(day !=0) %>% 
  bind_rows(cm %>% filter(day == 0) %>% 
              mutate(midge = "Midges") %>% 
              select(-algae_conc2) %>% 
              full_join(meta %>% 
                          select(algae_conc2) %>% 
                          unique(), 
                        by = character())) %>% 
  mutate(ac = ifelse(day ==0, algae_conc2-0.15*algae_conc2, ifelse(day == 14, algae_conc2, algae_conc2+0.2*algae_conc2)))

wtplot <- wtpredict %>% 
  full_join(cm %>% 
              filter(day == 0) %>% 
              group_by(day) %>% 
              summarise(wt = mean(wt)) %>% 
              full_join(meta %>% select(algae_conc2, midge) %>% unique(), by = character())) %>% 
  arrange(algae_conc2, midge, day) %>% 
  group_by(algae_conc2, midge) %>% 
  mutate(ac1 = ifelse(day == 14, algae_conc2 - 0.15*algae_conc2, algae_conc2),
         ac2 = ifelse(day == 22, algae_conc2 + 0.2 * algae_conc2, algae_conc2),
         ystart = lag(wt)) 


cm4plot %>% 
  ggplot(aes(x = algae_conc2, y = wt, group = interaction(day, algae_conc2)))+
  geom_line(aes(color = day, group = day), data = wtplot)+
  geom_point(aes(x = ac, fill = day), shape = 21, position = position_jitter(width = 0.05, height = 0), alpha = 0.5)+
  scale_x_continuous(trans = "log", breaks = c(0.01, 0.1, 1))+
  facet_wrap(~midge, nrow = 2)+
  geom_segment(aes(group = algae_conc2, x = ac1, xend = ac2, y = ystart, yend = wt), size = 0.9, color = "black", data = wtplot)+
  labs(x = "Sediment Treatment",
       y = "Midge Ash Free Dry Mass (mg)")+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")


#===

lpredict <- topredict %>% 
  mutate(fit = predict(llog, newdata = ., re.form = NA, type = "response")) %>% 
  group_by(day, algae_conc2, midge) %>% 
  summarise(body_size = mean(fit)) 

lplot <-  lpredict %>% 
  full_join(cm %>% 
              filter(day == 0) %>% 
              group_by(day) %>% 
              summarise(body_size = mean(body_size)) %>% 
              full_join(meta %>% select(algae_conc2, midge) %>% unique(), by = character())) %>% 
  arrange(algae_conc2, midge, day) %>% 
  group_by(algae_conc2, midge) %>% 
  mutate(ac1 = ifelse(day == 14, algae_conc2 - 0.15*algae_conc2, algae_conc2),
         ac2 = ifelse(day == 22, algae_conc2 + 0.2 * algae_conc2, algae_conc2),
         ystart = lag(body_size))

cm4plot %>% 
  ggplot(aes(x = algae_conc2, y = body_size, group = interaction(day, algae_conc2)))+
  geom_line(aes(color = day, group = day), data = lplot)+
  geom_point(aes(x = ac, fill = day), shape = 21, position = position_jitter(width = 0.05, height = 0), alpha = 0.5)+
  scale_x_continuous(trans = "log", breaks = c(0.01, 0.1, 1))+
  facet_wrap(~midge, nrow = 2)+
  geom_segment(aes(group = algae_conc2, x = ac1, xend = ac2, y = ystart, yend = body_size), size = 0.9, color = "black", data = lplot)+
  labs(x = "Sediment Treatment",
       y = "Midge Size (mm)")+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")

