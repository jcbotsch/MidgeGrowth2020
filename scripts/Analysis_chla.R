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


midge_color <- scale_color_manual(values = c("Midges" = "black", "No Midges" = "gray60"))
midge_fill <- scale_fill_manual(values = c("Midges" = "black", "No Midges" = "gray60"))



#====read in files====

meta <- read_csv("clean data/MG_meta.csv") %>% 
  mutate(algae_conc2 = ifelse(algae_conc == 0, 2e-3, algae_conc)) #half of lowest value
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
  scale_y_continuous(labels = scales::percent_format(),
                     limits = c(0, NA))+
  labs(x = "Fraction of Surface Sediment",
       y = "Organic Content")


om %>% 
  left_join(meta) %>% 
  ggplot(aes(algae_conc2, perc.org*100))+
  geom_point(size =3)+
  scale_y_continuous(trans = "log")+
  scale_x_continuous(trans = "log", breaks = c(0.01, 0.1, 1))+
  labs(x = "Fraction of Surface Sediment",
       y = "Organic Content")

#corrected for pheophytin
corrchlplot <- chl %>% 
  group_by(algae) %>% 
  mutate(avg.chl = mean(chl)/1000) %>% 
  left_join(meta) %>% 
  ggplot(aes(x = algae_conc, y = chl/1000))+
  geom_point(size = 2, alpha = 0.8)+
  geom_point(aes(y = avg.chl), color = "red", size = 2, shape = 0, stroke = 1)+
  scale_y_continuous(labels = scales::comma_format())+
  labs(x = "",
       y = expression("Chlorophyll mg"~L^{-1})) +  
  theme(axis.text.x = element_blank())

chl %>% 
  gather(var, val, uncorr_chl:chl) %>% 
  mutate(val= val/1000) %>% 
  group_by(algae) %>% 
  left_join(meta) %>% 
  ggplot(aes(x = algae_conc2, y = val))+
  facet_wrap(~var, scales = "free")+
  geom_point(size = 2, alpha = 0.8)+
  scale_y_continuous(trans = "log", labels = scales::comma_format())+
  scale_x_continuous(trans = "log", breaks = c(0.01, 0.1, 1))+
  labs(x = "Sediment Treatment",
       y = expression("mg"~L^{-1})) 



#uncorrected
uncorrchlplot <- chl %>% 
  left_join(meta) %>% 
  ggplot(aes(x = algae_conc, y = uncorr_chl/1000))+
  geom_point(size = 2, alpha = 0.8)+
  scale_y_continuous(labels = scales::comma_format())+
  labs(x = "",
       y = "Chlorophyll mg/L") + 
  theme(axis.text.x = element_blank())


fig1 <- grid.arrange(corrchlplot, omplot, nrow  = 2)
plot(fig1)

# ggpreview(plot = fig1, dpi = 650, width = 80, height = 100, units = "mm")
# ggsave(plot = fig1, filename = "Botsch_MG_Fig1.pdf", dpi = 650, device = "pdf", width = 80, height = 80, units = "mm")



grid.arrange(uncorrchlplot, omplot, nrow  = 2)

#====Meta====
meta <- meta %>% 
  left_join(chl %>% 
              group_by(algae) %>% 
              summarise(init.chl = mean(chl)/1000))


#====create data frame to predict model response variables ====

topredict <- data.frame(init.chl = rep(unique(meta$init.chl), 4),
                        midge = rep(rep(c("Midges", "No Midges"), each = 10), 2),
                        box = rep(c("1", "2"), each = 20))


#=====Question 1: Effects of Midges on Primary Producers=====
## Midge effect on production
nep <- ep_raw %>% 
  left_join(meta) %>% 
  mutate(coreid = as.character(coreid),
         box = as.character(box))

nep14 <- nep %>% 
  filter(day == 14) %>% 
  mutate(midge = fct_relevel(midge, c("No Midges", "Midges")))

g14log <- lm(log(gpp)~log2(init.chl)*midge+box, data = nep14) 
summary(g14log) 
summary(update(g14log, .~. -log2(init.chl):midge))
# plot(g14log)


nep22 <- nep %>% 
  filter(day == 22)  %>% 
  mutate(midge = fct_relevel(midge, c("No Midges", "Midges")))

g22log <- lm(log(gpp)~log2(init.chl)*midge+box, data = nep22)
summary(g22log)
summary(update(g22log, .~. -log2(init.chl):midge))

# plot(g22log)


predicted14 <- topredict %>% 
  mutate(gpp = predict(g14log, ., type = "response")) %>% 
  group_by(init.chl, midge) %>% 
  summarise(gpp = exp(mean(gpp))) %>% 
  mutate(day = "Day 14")

predicted22 <- topredict %>% 
  mutate(gpp = predict(g22log, ., type = "response")) %>% 
  group_by(init.chl, midge) %>% 
  summarise(gpp = exp(mean(gpp))) %>% 
  mutate(day = "Day 22")

gpredict <- rbind(predicted14, predicted22)

gpredict %>% 
  mutate(day = str_replace(day, " ", "_")) %>% 
  spread(day, gpp) %>% 
  ggplot(aes(x = Day_14, y = Day_22, col = midge))+
  geom_line(aes(group = factor(init.chl)), color = "black")+
  geom_point()+
  geom_abline(slope = 1)+
  lims(x = c(0,0.026),
       y = c(0,0.026))+
  coord_equal()



#====Figure 2: GPP====
nepfig <- nep %>% 
  mutate(day = paste("Day", day)) %>% 
  ggplot(aes(x = init.chl, y = gpp, fill = midge))+
  facet_wrap(~day, ncol = 1)+
  geom_point(size  = 2, alpha = 0.75, shape = 21)+
  geom_line(aes(color = midge), size = 0.7, data = gpredict)+
  midge_color+
  midge_fill+ 
  labs(x = "Average Initial Chlorophyll (mg/L)",
       y = expression("GPP "~(g~O[2]~m^{-2}~hr^{-1})),
       color = "",
       fill = "")+
  scale_y_continuous(trans = "log2", breaks = c(0.005, 0.015, 0.05))+
  scale_x_continuous(trans = "log2", breaks = c(1, 3, 9))+
  theme(legend.position = c(0.5,0.95),
        legend.direction = "horizontal",
        legend.text.align = 0,
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA, color = NA))

plot(nepfig)
# ggpreview(plot = nepfig, dpi = 650, width = 80, height = 80, units = "mm")

#Supplemental Figure
nep %>% 
  gather(metabolism, value, resp, nep, gpp) %>% 
  mutate(day = day+rnorm(n(), 0.001), #add jitter
         metabolism = toupper(ifelse(metabolism == "resp", "er", metabolism)),
         metabolism = fct_reorder(metabolism, desc(value))) %>%  
  ggplot(aes(x = day, y = value, col = init.chl))+
  facet_grid(metabolism~midge, scales = "free")+
  geom_hline(yintercept = 0)+
  geom_point()+
  geom_line(aes( group = coreid), alpha = 0.4)+
  geom_smooth(method = "lm", col = "black", se = FALSE)+
  theme(legend.position = "bottom")+
  labs(y = expression("Change in" ~O[2] (g~m^{-2}~hr^{-1})),
       x = "Incubation Day",
       color = "Average Initial Chlorophyll (mg/L)")+
  scale_color_viridis_c(trans = "log10")+
  scale_x_continuous(breaks = c(14, 22))

nep %>% 
  mutate(day = paste("Day", day)) %>% 
  gather(metabolism, value, resp, nep, gpp) %>% 
  mutate(metabolism = toupper(metabolism),
         metabolism = ifelse(metabolism == "RESP", "ER", metabolism),
         metabolism = fct_reorder(metabolism, -value),
         value = ifelse(metabolism == "ER", -value, value)) %>% 
  ggplot(aes(x = init.chl, y = value, fill = midge))+
  geom_hline(yintercept = 0, linetype = "dotted")+
  facet_grid(metabolism~day, scale = "free_y", switch = "y")+
  geom_point(size  = 2, alpha = 0.7, shape = 21)+
  geom_smooth(aes(color = midge), method = "lm", size = 0.75, se = FALSE)+
  midge_color+
  midge_fill+ 
  labs(x = "Average Initial Chlorophyll (mg/L)",
       y = expression(g~O[2]~m^{-2}~hr^{-1}),
       color = "",
       fill = "")+
    scale_x_continuous(trans = "log2", breaks = c(1,3,9))+
  theme(strip.placement = "outside")


# ggpreview(plot = last_plot(), dpi = 650, width = 120, height = 120, units = "mm")

r14log <- lm(log(-resp)~ log2(init.chl)*midge + box, data = nep14)
summary(r14log)

r22log <- lm(log(-resp)~ log2(init.chl)*midge + box, data = nep22)
summary(r22log)


nep %>% 
  mutate(day = paste("Day", day,  sep = " ")) %>% 
  ggplot(aes(x = -resp, y = gpp, fill = init.chl))+
  geom_point(shape = 21, alpha = 0.5, size = 2)+
  facet_grid(midge~day)+
  geom_abline(slope = 1)+
  coord_equal()+
  # lims(x = c(0,0.06),
  #      y = c(0,0.06))+
  scale_fill_viridis_c(trans = "log2", breaks = c(0.01, 0.1, 1))+
  labs(x = expression("ER" ~ (g~O[2]~m^{-2}~hr^{-1})),
       y = expression("GPP"~(g~O[2]~m^{-2}~hr^{-1})),
       fill = "Average Initial Chlorophyll (mg/L)")

#====Question 2: Sediment Effects on Midges =====
##====Q2.1: Number of Midges Present====
cc <- cc_raw %>% 
  left_join(meta) %>% 
  mutate(coreid = as.character(coreid),
         box = as.character(box))

cc14 <- cc %>% 
  filter(day == 14) %>% 
  mutate(midge = fct_relevel(midge, c("No Midges", "Midges")))


nm14log <- glm(live_tt~log2(init.chl)*midge+box, 
               data = cc14,
               family = quasipoisson())

summary(nm14log)
# plot(nm14log)

summary(update(nm14log, .~.-log2(init.chl):midge))



cc22 <- cc %>% 
  filter(day == 22) %>% 
  mutate(midge = fct_relevel(midge, c("No Midges", "Midges")))

nm22log <- glm(live_tt~log2(init.chl)*midge+box, 
               data = cc22,
               family = quasipoisson())

summary(nm22log)
# plot(nm22log)

summary(update(nm22log, .~. -log2(init.chl):midge))

#extract fits

predicted14 <- topredict %>% 
  mutate(fit = predict(nm14log, ., type = "link")) %>% 
  group_by(init.chl, midge) %>% 
  summarise(live_tt = exp(mean(fit))) %>% 
  mutate(day = "Day 14")

predicted22 <- topredict %>% 
  mutate(fit = predict(nm22log, ., type = "link")) %>% 
  group_by(init.chl, midge) %>% 
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

nmidge_fig <- cc %>% 
  filter(day %in% c(14, 22)) %>% 
  mutate(day = paste("Day", day)) %>% 
  ggplot(aes(x = init.chl, y = live_tt, fill = midge))+
  facet_wrap(~day, ncol = 1)+
  geom_point(size  = 2, alpha = 0.7, shape = 21)+
  geom_hline(yintercept = 20, linetype = "dashed", alpha = 0.5)+
  # geom_smooth(aes(color = midge), method = "glm", size = 0.75, se = FALSE, method.args = list(family = quasipoisson()) )+
  geom_line(aes(color = midge), data = nmpredict)+
  midge_color+
  midge_fill+ 
  labs(x = "Average Initial Chlorophyll (mg/L)",
       y = "Number of Larvae",
       color = "",
       fill = "")+
  scale_x_continuous(trans = "log2", breaks = c(1,3,9))+
  theme(legend.position = c(0.5,0.95),
        legend.text.align = 0,
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA, color = NA),
        legend.direction = "horizontal")

plot(nmidge_fig)
ggpreview(plot = nmidge_fig, dpi = 650, width = 80, height = 80, units = "mm")


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
    group_by(coreid, day, midge, init.chl, box, instar, measured) %>% 
    count(instar, name = "s2") %>% 
    filter(instar == 2) %>% 
    mutate(prop = s2/measured)}$prop


p2 <- cm %>% 
  filter(species_name == "tt",
         day!= 0) %>%  #remove the four non-tanytarsini
  add_count(coreid, name = "measured") %>% 
  group_by(coreid, day, midge, init.chl, box, instar, measured) %>% 
  count(instar, name = "count") %>% 
  mutate(instar = paste0("s", instar)) %>% 
  spread(instar, count, fill = 0) %>% 
  mutate(prop = s2/measured,
         others = s3+s4)

#subset for day 14
p314 <- p2 %>% 
  filter(day == 14) %>% 
  mutate(midge = factor(midge, levels = c("No Midges", "Midges")))

prop14log <- glm(cbind(s2, others) ~ log2(init.chl)*midge+box, 
                 data = p314, 
                 family = quasibinomial())

summary(prop14log)
# plot(prop14log)
summary(update(prop14log, .~. -log2(init.chl):midge))


#subset for day 22
p322 <- p2 %>% 
  filter(day == 22) %>% 
  mutate(midge = factor(midge, levels = c("No Midges", "Midges")))

prop22log <- glm(cbind(s2, others) ~log2(init.chl)*midge+box, 
                 data = p322, 
                 family = quasibinomial())
summary(prop22log)
# plot(prop22log)
summary(update(prop22log, .~. -log2(init.chl):midge))



predicted14 <- topredict %>% 
  mutate(fit = predict(prop14log, ., type = "link")) %>% 
  group_by(init.chl, midge) %>% 
  summarise(prop = exp(mean(fit))) %>% 
  mutate(day = "Day 14")

predicted22 <- topredict %>% 
  mutate(fit = predict(prop22log, ., type = "link")) %>% 
  group_by(init.chl, midge) %>% 
  summarise(prop = exp(mean(fit))) %>% 
  mutate(day = "Day 22")

propredict <- rbind(predicted14, predicted22)



#====Figure 4: Instar====

prop_fig <- p2 %>% 
  mutate(day = paste("Day", day)) %>% 
  ggplot(aes(y = prop, x = init.chl))+
  geom_point(aes(fill = midge), shape = 21, alpha = 0.75, size = 2)+
  geom_line(aes(color = midge), data = propredict)+
  geom_hline(yintercept = startingprop, alpha = 0.75, linetype = "dashed")+
  facet_wrap(~day, nrow = 2)+
  lims(y = c(0,1))+
  labs(x = "Average Initial Chlorophyll (mg/L)",
       y = "Proportion of Larvae in 2nd Instar",
       color = NULL,
       fill = NULL)+
  midge_color+
  midge_fill+ 
  theme(legend.position = "bottom")+
    scale_x_continuous(trans = "log2", breaks = c(1,3,9))+
  theme(legend.position = c(0.8,0.35),
        legend.text.align = 0,
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA, color = NA),
        legend.direction = "vertical")

prop_fig
# ggpreview(plot = prop_fig, dpi = 650, width = 80, height = 80, units = "mm")
# ggsave(plot = last_plot(), filename = "Botsch_MG_Fig4.pdf", device = "pdf", dpi = 650, width = 80, height = 120, units = "mm")

#====Q2.3: Midge Length====
#subset day 14
l314 <- cm %>% 
  filter(day == 14) %>% 
  mutate(midge = factor(midge, levels = c("No Midges", "Midges")),
         instar = factor(instar))


bl14log <- lmer(body_size~log2(init.chl)*midge+box+
                  (1|coreid), 
                data = l314) 

summary(bl14log)
# plot(bl14log)

Anova(bl14log, type = "3", test.statistic = "F")
Anova(update(bl14log, .~.-log2(init.chl):midge), type = "2", test.statistic = "F")

#subset day 22
l322 <- cm %>% 
  filter(day == 22) %>% 
  mutate(midge = factor(midge, levels = c("No Midges", "Midges")),
         instar = factor(instar))


bl22log <- lmer(body_size~log2(init.chl)*midge+box +
                  (1|coreid), 
                data = l322)
summary(bl22log)
# plot(bl22log)

#F Tests
Anova(bl22log, type = "3", test.statistic = "F") #includes interaction 
Anova(update(bl22log, .~. -log2(init.chl):midge), type = "2", test.statistic = "F") #drops interaction


predicted14 <- topredict %>% 
  mutate(fit = predict(bl14log, newdata = ., re.form = NA, type = "response")) %>% 
  group_by(init.chl, midge) %>% 
  summarise(body_size = mean(fit)) %>% 
  mutate(day = "Day 14")

predicted22 <- topredict %>% 
  mutate(fit = predict(bl22log, newdata = ., re.form = NA, type = "response")) %>% 
  group_by(init.chl, midge) %>% 
  summarise(body_size = mean(fit)) %>% 
  mutate(day = "Day 22")

blpredict <- rbind(predicted14, predicted22)


#====Figure 5: Body Length====

#average body legnth in initial midges
start_length <- {cm %>% 
    filter(day == 0) %>% 
    summarise(length = meanna(body_size))}$length

blfig <- cm %>% 
  filter(day %in% c(14, 22)) %>% 
  mutate(day = paste("Day", day)) %>% 
  ggplot(aes(x = init.chl, y = body_size))+
  facet_wrap(~day, ncol = 1)+
  geom_hline(yintercept = start_length, linetype = 2, size = 0.5)+
  geom_point(aes(fill = midge), size  = 1.5, alpha = 0.75, shape = 21, stroke = 0.2, position = position_jitterdodge(dodge.width = 0.01, jitter.width = 0.01, jitter.height = 0))+
  geom_line(aes(col = midge), size = 0.7, data = blpredict)+
  midge_color+
  midge_fill+ 
  labs(y = "Body Length (mm)",
       x = "Average Initial Chlorophyll (mg/L)",
       color = "",
       fill = "")+
    scale_x_continuous(trans = "log2", breaks = c(1,3,9))+
  theme(legend.position = c(0.5,0.95),
        legend.text.align = 0,
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA, color = NA),
        legend.direction = "horizontal")

plot(blfig)
# ggpreview(plot = blfig, dpi = 650, width = 80, height = 80, units = "mm")



# ggpreview(plot = blfig, dpi = 650, width = 80, height = 80, units = "mm")
# ggsave(plot = last_plot(), filename = "Botsch_MG_Fig5.pdf", device = "pdf", dpi = 650, width = 80, height = 80, units = "mm")

cm %>% 
  filter(!is.na(instar)) %>% 
  ggplot(aes(x = factor(instar), y = body_size))+
  geom_jitter(height = 0, shape = 21, alpha = 0.5, fill = "gray50")+
  scale_fill_viridis(trans = "log2", breaks = c(0.01, 0.1, 1))+
  labs(y = "Body Length (mm)",
       x = "Instar")
# ggpreview(plot = last_plot(), dpi = 650, width = 80, height = 80, units = "mm")


