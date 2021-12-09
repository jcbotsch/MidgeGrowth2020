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


#====create data frame to predict model response variables ====

topredict <- data.frame(algae_conc2 = rep(unique(meta$algae_conc2), 4),
                        midge = rep(rep(c("Midges", "No Midges"), each = 10), 2),
                        box = rep(c("1", "2"), each = 20))

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


fig1 <- grid.arrange(corrchlplot, omplot, nrow  = 2)
plot(fig1)

# ggpreview(plot = fig1, dpi = 650, width = 80, height = 100, units = "mm")
# ggsave(plot = fig1, filename = "Botsch_MG_Fig1.pdf", dpi = 650, device = "pdf", width = 80, height = 80, units = "mm")



grid.arrange(uncorrchlplot, omplot, nrow  = 2)


#=====Question 1: Effects of Midges on Primary Producers=====
## Midge effect on production
nep <- ep_raw %>% 
  left_join(meta) %>% 
  mutate(coreid = as.character(coreid),
         box = box-1)

nep14 <- nep %>% 
  filter(day == 14) %>% 
  mutate(midge = fct_relevel(midge, c("No Midges", "Midges")))

g14log <- lm(log(gpp)~log(algae_conc2)*midge+box, data = nep14) 
summary(g14log) 
summary(update(g14log, .~. -log(algae_conc2):midge))
# plot(g14log)


nep22 <- nep %>% 
  filter(day == 22)  %>% 
  mutate(midge = fct_relevel(midge, c("No Midges", "Midges")))

g22log <- lm(log(gpp)~log(algae_conc2)*midge+box, data = nep22)
summary(g22log)
summary(update(g22log, .~. -log(algae_conc2):midge))

# plot(g22log)


predicted14 <- topredict %>% 
  mutate(gpp = predict(g14log, ., type = "response")) %>% 
  group_by(algae_conc2, midge) %>% 
  summarise(gpp = exp(mean(gpp))) %>% 
  mutate(day = "Day 14")

predicted22 <- topredict %>% 
  mutate(gpp = predict(g22log, ., type = "response")) %>% 
  group_by(algae_conc2, midge) %>% 
  summarise(gpp = exp(mean(gpp))) %>% 
  mutate(day = "Day 22")

gpredict <- rbind(predicted14, predicted22)

gpredict %>% 
  mutate(day = str_replace(day, " ", "_")) %>% 
  spread(day, gpp) %>% 
  ggplot(aes(x = Day_14, y = Day_22, col = midge))+
  geom_line(aes(group = factor(algae_conc2)), color = "black")+
  geom_point()+
  geom_abline(slope = 1)+
  lims(x = c(0,0.026),
       y = c(0,0.026))+
  coord_equal()



#====Figure 2: GPP====
nepfig <- nep %>% 
  mutate(day = paste("Day", day)) %>% 
  ggplot(aes(x = algae_conc2, y = gpp, fill = midge))+
  facet_wrap(~day, ncol = 1)+
  geom_point(size  = 2, alpha = 0.75, shape = 21, position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.1, jitter.height = 0))+
  geom_line(aes(color = midge), size = 0.7, data = gpredict)+
  midge_color+
  midge_fill+ 
  labs(x = "Sediment Treatment",
       y = expression("GPP "~(g~O[2]~m^{-2}~hr^{-1})),
       color = "",
       fill = "")+
  scale_y_continuous(trans = "log", breaks = c(0.005, 0.015, 0.05))+
  scale_x_continuous(trans = "log", breaks = c(0.01, 0.1, 1))+
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
  ggplot(aes(x = day, y = value, col = algae_conc2))+
  facet_grid(metabolism~midge, scales = "free")+
  geom_hline(yintercept = 0)+
  geom_point()+
  geom_line(aes( group = coreid), alpha = 0.4)+
  geom_smooth(method = "lm", col = "black", se = FALSE)+
  theme(legend.position = "bottom")+
  labs(y = expression("Change in" ~O[2] (g~m^{-2}~hr^{-1})),
       x = "Incubation Day",
       color = "Sediment Treatment")+
  scale_color_viridis_c(trans = "log10")+
  scale_x_continuous(breaks = c(14, 22))

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
  midge_color+
  midge_fill+ 
  labs(x = "Sediment Treatment",
       y = expression(g~O[2]~m^{-2}~hr^{-1}),
       color = "",
       fill = "")+
  scale_x_continuous(trans = "log", breaks = c(0.01, 0.1, 1))+
  theme(strip.placement = "outside")


# ggpreview(plot = last_plot(), dpi = 650, width = 120, height = 120, units = "mm")

r14log <- lm(log(-resp)~ log(algae_conc2)*midge + box, data = nep14)
summary(r14log)

r22log <- lm(log(-resp)~ log(algae_conc2)*midge + box, data = nep22)
summary(r22log)


nep %>% 
  mutate(day = paste("Day", day,  sep = " ")) %>% 
  ggplot(aes(x = -resp, y = gpp, fill = algae_conc2))+
  geom_point(shape = 21, alpha = 0.5, size = 2)+
  facet_grid(midge~day)+
  geom_abline(slope = 1)+
  coord_equal()+
  # lims(x = c(0,0.06),
  #      y = c(0,0.06))+
  scale_fill_viridis_c(trans = "log", breaks = c(0.01, 0.1, 1))+
  labs(x = expression("ER" ~ (g~O[2]~m^{-2}~hr^{-1})),
       y = expression("GPP"~(g~O[2]~m^{-2}~hr^{-1})),
       fill = "Sediment Treatment")

#====Question 2: Sediment Effects on Midges =====
##====Q2.1: Number of Midges Present====
cc <- cc_raw %>% 
  left_join(meta) %>% 
  mutate(coreid = as.character(coreid),
         box = box)

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

nmidge_fig <- cc %>% 
  filter(day %in% c(14, 22)) %>% 
  mutate(day = paste("Day", day)) %>% 
  ggplot(aes(x = algae_conc2, y = live_tt, fill = midge))+
  facet_wrap(~day, ncol = 1)+
  geom_point(size  = 2, alpha = 0.7, shape = 21, position = position_jitter(width = 0.1, height = 0))+
  geom_hline(yintercept = 20, linetype = "dashed", alpha = 0.5)+
  # geom_smooth(aes(color = midge), method = "glm", size = 0.75, se = FALSE, method.args = list(family = quasipoisson()) )+
  geom_line(aes(color = midge), data = nmpredict)+
  midge_color+
  midge_fill+ 
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

plot(nmidge_fig)
# ggpreview(plot = nmidge_fig, dpi = 650, width = 80, height = 80, units = "mm")


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

prop_fig <- p2 %>% 
  mutate(day = paste("Day", day)) %>% 
  ggplot(aes(y = prop, x = algae_conc2))+
  geom_point(aes(fill = midge), shape = 21, alpha = 0.75, size = 2, position = position_jitterdodge(dodge.width = 0.1, jitter.width = 0.1, jitter.height = 0))+
  geom_line(aes(color = midge), data = propredict)+
  geom_hline(yintercept = startingprop, alpha = 0.75, linetype = "dashed")+
  facet_wrap(~day, nrow = 2)+
  lims(y = c(0,1))+
  labs(x = "Sediment Treatment",
       y = "Proportion of Larvae in 2nd Instar",
       color = NULL,
       fill = NULL)+
  midge_color+
  midge_fill+ 
  theme(legend.position = "bottom")+
  scale_x_continuous(trans = "log", breaks = c(0.01, 0.1, 1))+
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


bl14log <- lmer(body_size~log(algae_conc2)*midge+box+
                  (1|coreid), 
                data = l314) 

summary(bl14log)
# plot(bl14log)

Anova(bl14log, type = "3", test.statistic = "F")
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
Anova(update(bl22log, .~. -log(algae_conc2):midge), type = "2", test.statistic = "F") #drops interaction


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

blfig <- cm %>% 
  filter(day %in% c(14, 22)) %>% 
  mutate(day = paste("Day", day)) %>% 
  ggplot(aes(x = algae_conc2, y = body_size))+
  facet_wrap(~day, ncol = 1)+
  geom_hline(yintercept = start_length, linetype = 2, size = 0.5)+
  geom_point(aes(fill = midge), size  = 1.5, alpha = 0.75, shape = 21, stroke = 0.2, position = position_jitterdodge(dodge.width = 0.4, jitter.width = 0.2, jitter.height = 0))+
  geom_line(aes(col = midge), size = 0.7, data = blpredict)+
  midge_color+
  midge_fill+ 
  labs(y = "Body Length (mm)",
       x = "Sediment Treatment",
       color = "",
       fill = "")+
  scale_x_continuous(trans = "log", breaks = c(0.01, 0.1, 1))+
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
  scale_fill_viridis(trans = "log", breaks = c(0.01, 0.1, 1))+
  labs(y = "Body Length (mm)",
       x = "Instar")
# ggpreview(plot = last_plot(), dpi = 650, width = 80, height = 80, units = "mm")



#=====Primary and Secondary Production=====

#subset starting midges and experimental midges
start_cm <- cm %>% 
  filter(day == 0) %>% 
  select(coreid, sampledate, day, species_name, head_size, body_size, instar) 

exp_cm <- cm %>% 
  filter(day!=0,
         midge == "Midges",
         !is.na(body_size),
         species_name == "tt") 

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
production_boot <- estimated_growth %>% 
  mutate(wt = weight(body_size)) %>% 
  left_join(nts) %>% 
  group_by(bootstrap, algae_conc2) %>% 
  mutate(ntlag = lag(nt),
         wtlag = lag(wt),
         deltaT = day-lag(day),
         Pd = incsumprod(n1 = ntlag, n2 = nt, wt1 = wtlag, wt2 = wt, deltaT)$Pd,
         g = incsumprod(n1 = ntlag, n2 = nt, wt1 = wtlag, wt2 = wt, deltaT)$g,
         gd = g/deltaT) 


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
         g = incsumprod(n1 = nt_day0, n2 = nt, wt1 = wt_day0, wt2 = wt, day)$g,
         gd = g/day) 


  

ftube_r = 0.03/2 #30 mm/ 2 to get radius /1000 to convert to m
ftube_area = ftube_r^2*pi
pq = 1 #photosynthetic quotient

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

#combine with estimates of primary production
prod1 <- nep %>% 
  group_by(day, midge, algae_conc2) %>% 
  mutate(prod = (gpp*18)*(12/32*pq)) %>% #daily production of algae in g C m^2d^-1
  add_count() %>% 
  summarise(n = unique(n),
            mean_pp = mean(prod),
            sd_pp = sd(prod)/sqrt(n)) %>% 
  full_join(production_boot %>% 
              # filter(Pd>0) %>% 
              mutate(Pdc = ((Pd*0.5)/1000)/ftube_area,
                     midge = "Midges") %>%  #0.5 g C/ g AFDM
              group_by(day, algae_conc2, midge) %>% 
              summarise(mean_sp = mean(Pdc),
                        sd_sp = sd(Pdc),
                        Bt = mean(nt*wt))) %>% 
  filter(day!=0,
         midge == "Midges") 


prod2 <- nep %>% 
  group_by(day, midge, algae_conc2) %>% 
  mutate(prod = (gpp*18)*(12/32*pq)) %>% #daily production of algae in g C m^2d^-1
  add_count() %>% 
  summarise(n = unique(n),
            mean_pp = mean(prod),
            sd_pp = sd(prod)/sqrt(n)) %>% 
  full_join(production_boot2 %>% 
              # filter(Pd>0) %>% 
              mutate(Pdc = ((Pd*0.5)/1000)/ftube_area,
                     midge = "Midges") %>%  #0.5 g C/ g AFDM
              group_by(day, algae_conc2, midge) %>% 
              summarise(mean_sp = mean(Pdc),
                        sd_sp = sd(Pdc),
                        Bt = mean(nt*wt))) %>% 
  filter(day!=0,
         midge == "Midges") %>% 
  ungroup %>% 
  mutate(unique_id =1:n())

nep %>% 
  filter(midge == "Midges") %>% 
  group_by(midge, algae_conc2) %>% 
  mutate(prod = (gpp*18)*(12/32*pq)) %>% #daily production of algae in g C m^2d^-1
  add_count() %>% 
  summarise(n = unique(n),
            mean_pp = mean(prod),
            sd_pp = sd(prod)/sqrt(n)) %>% 
  left_join(total_sp)%>% 
  ggplot(aes(x = mean_pp, y = mean_sp))+
  geom_errorbarh(aes(xmin = mean_pp-sd_pp, xmax = mean_pp+sd_pp), height = NA)+
  geom_errorbar(aes(ymin = mean_sp-sd_sp, ymax = mean_sp+sd_sp), width  = NA)+
  geom_point(aes(fill = algae_conc2), shape = 21, size = 2)+
  viridis::scale_fill_viridis(trans = "log", breaks = c(0.01, 0.1, 1))+
  coord_cartesian(xlim = c(0, NA),
                  ylim = c(0, NA))+
  labs(x = expression("Primary Production g C"~m^{-2}~d^{-1}),
       y = expression("Secondary Production g C"~m^{-2}~d^{-1}),
       fill = "Sediment\nTreatment")

library(phytools)

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
#====Figure 6: Primary and Secondary Production=====
prodfig <- prod1 %>% 
  mutate(day = paste("Day", day)) %>% 
  ggplot(aes(x = mean_pp, y = mean_sp))+
  facet_wrap(~day, ncol = 1)+
  geom_errorbarh(aes(xmin = mean_pp-sd_pp, xmax = mean_pp+sd_pp))+
  geom_errorbar(aes(ymin = mean_sp-sd_sp, ymax = mean_sp+sd_sp), width  = NA)+
  geom_point(aes(fill = algae_conc2), shape = 21, size = 2)+
  viridis::scale_fill_viridis(trans = "log", breaks = c(0.01, 0.1, 1))+
  coord_cartesian(xlim = c(0, NA),
                  ylim = c(0, NA))+
  labs(x = expression("Primary Production g C"~m^{-2}~d^{-1}),
       y = expression("Secondary Production g C"~m^{-2}~d^{-1}),
       fill = "Sediment\nTreatment")

nd = data.frame(mean_pp = rep(seq(0, 0.22, by = 0.01), 2)) %>% crossing(day = c(14, 22))

nd$psp <- predict(prodmod, nd)

nd$psp_se <- predict(prodmod, nd, se = TRUE)$se.fit

prodfig2 <- prod2 %>% 
  mutate(day = paste("Day", day)) %>% 
  ggplot(aes(x = mean_pp, y = mean_sp))+
  facet_wrap(~day, ncol = 1)+
  geom_errorbarh(aes(xmin = mean_pp-sd_pp, xmax = mean_pp+sd_pp))+
  geom_errorbar(aes(ymin = mean_sp-sd_sp, ymax = mean_sp+sd_sp), width  = NA)+
  geom_point(aes(fill = algae_conc2), shape = 21, size = 2)+
  viridis::scale_fill_viridis(trans = "log", breaks = c(0.01, 0.1, 1))+
  coord_cartesian(xlim = c(0, NA))+
  labs(x = expression("Primary Production g C"~m^{-2}~d^{-1}),
       y = expression("Secondary Production g C"~m^{-2}~d^{-1}),
       fill = "Sediment\nTreatment")

prod2 %>% 
  mutate(day = paste("Day", day)) %>% 
  ggplot(aes(x = mean_pp, y = mean_sp))+
  geom_abline(slope = prodmod$beta[2], intercept = prodmod$beta[1])+
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


# ggpreview(plot = last_plot(), dpi = 650, width = 80, height = 100, units = "mm")


prod2 %>% 
  mutate(day = paste("Day", day),
         mean_pp = mean_pp*algae_conc2,
         sd_pp = sd_pp*algae_conc2) %>% 
  ggplot(aes(x = mean_pp, y = mean_sp))+
  geom_ribbon(aes(ymin = psp-psp_se, y = psp, ymax = psp+psp_se, group = day), fill = "gray20", alpha = 0.2, data = nd %>% mutate(day = paste("Day", day)))+
  geom_errorbarh(aes(xmin = mean_pp-sd_pp, xmax = mean_pp+sd_pp), alpha = 0.75)+
  geom_errorbar(aes(ymin = mean_sp-sd_sp, ymax = mean_sp+sd_sp), alpha = 0.75, width  = NA)+
  geom_point(aes(shape = day, fill = algae_conc2), size = 2)+
  viridis::scale_fill_viridis(trans = "log", breaks = c(0.01, 0.1, 1))+
  geom_line(aes(y = psp, linetype = day), data = nd %>% mutate(day = paste("Day", day)))+
  # coord_cartesian(xlim = c(0, NA))+
  scale_x_continuous(trans = "log1p")+
  scale_shape_manual(values = c(21, 22))+
  labs(x = expression("Primary Production g C"~m^{-2}~d^{-1}),
       y = expression("Secondary Production g C"~m^{-2}~d^{-1}),
       fill = "Sediment\nTreatment",
       shape = element_blank(),
       linetype = element_blank())

prod2 %>% 
  mutate(day = paste0("Day_",day)) %>% 
  left_join(left_join(meta, chl %>% 
                        group_by(algae) %>% 
                        summarise(chl = mean(chl)))) %>% 
  mutate(pp_norm = mean_pp/chl) %>% 
  ggplot(aes(y = mean_sp, x = pp_norm, col = algae_conc2, shape = paste("Day", day)))+
  geom_point(size = 3)+
  scale_x_continuous(labels = scales::comma_format())+
  scale_color_viridis_c(trans = "log", breaks = c(0.001, 0.01, 0.1, 1))

prod2 %>% 
    mutate(day = paste0("Day_",day)) %>% 
  left_join(left_join(meta, chl %>% 
                        group_by(algae) %>% 
                        summarise(init_chl = mean(chl)))) %>% 
    ggplot(aes(y = mean_sp, x = mean_pp/init_chl, col = algae_conc2, shape = paste("Day", day)))+
    geom_point(size = 3)+
    scale_x_continuous(labels = scales::comma_format())+
    scale_color_viridis_c(trans = "log", breaks = c(0.001, 0.01, 0.1, 1))
  
nep %>% 
  ggplot(aes(x = midge, y = gpp/algae_conc2))+
  facet_wrap(~day)+
  geom_jitter(aes(color = algae_conc2), height = 0, width = 0.2)+
  scale_color_viridis_c(trans = "log", breaks = c(0.001, 0.01, 0.1, 1))


nep %>% 
  left_join(left_join(meta %>% select(algae_conc2, algae), chl %>% 
                        group_by(algae) %>% 
                        summarise(init_chl = mean(chl)))) %>% 
  unique() %>% 
  ggplot(aes(x = midge, y = gpp/init_chl))+
  facet_wrap(~paste("Day", day))+
  geom_jitter(aes(color = algae_conc2), height = 0, width = 0.2)+
  scale_color_viridis_c(trans = "log", breaks = c(0.001, 0.01, 0.1, 1))

nep %>% 
  select(coreid, algae_conc2, midge, day, gpp) %>% 
  mutate(day = paste0("gpp_day_", day)) %>% 
  spread(day, gpp) %>% 
  mutate(gpp_change = (gpp_day_22-gpp_day_14)/((gpp_day_22+gpp_day_14)/2)) %>% 
  filter(!is.na(gpp_change)) %>% 
  left_join(cc %>% select(coreid, live_tt)) %>% 
  ggplot(aes(y = gpp_change, x = live_tt, color = algae_conc2))+
  geom_jitter(width = 0.1, height = 0)+
  labs(x = "live Tanytarsini found at day 22",
       y = "Change in GPP/average GPP")+
  scale_color_viridis_c(trans = "log", breaks = c(0.001, 0.01, 0.1, 1))

nep %>% 
  select(coreid, algae_conc2, midge, day, gpp) %>% 
  mutate(day = paste0("gpp_day_", day)) %>% 
  spread(day, gpp) %>% 
  mutate(gpp_change = (gpp_day_22-gpp_day_14)/((gpp_day_22+gpp_day_14)/2)) %>% 
  filter(!is.na(gpp_change)) %>% 
  left_join(cc %>% select(coreid, live_tt)) %>% 
  ggplot(aes(y = gpp_change, x = midge, color = algae_conc2))+
  geom_jitter(width = 0.3, height = 0)+
  labs(x = "Midge Treatment",
       y = "Change in GPP/average GPP")+
  scale_color_viridis_c(trans = "log", breaks = c(0.001, 0.01, 0.1, 1))


nep %>% 
  select(coreid, algae_conc2, midge, day, gpp) %>% 
  mutate(day = paste0("gpp_day_", day)) %>% 
  spread(day, gpp) %>% 
  mutate(gpp_change = (gpp_day_22-gpp_day_14)) %>% 
  filter(!is.na(gpp_change)) %>% 
  left_join(cc %>% select(coreid, live_tt)) %>% 
  ggplot(aes(y = gpp_change, x = midge, color = algae_conc2))+
  geom_jitter(width = 0.3, height = 0)+
  labs(x = "Midge Treatment",
       y = "Change in GPP")+
  scale_color_viridis_c(trans = "log", breaks = c(0.001, 0.01, 0.1, 1))

#=======

startbt <- {day0 %>% summarise(Btlag = mean(ntlag*wtlag))}$Btlag


data <- prod2 %>% 
  filter(day == 22) %>% 
  left_join(chl %>% 
              left_join(meta) %>% 
              group_by(algae_conc2) %>% 
              summarise(chl = mean(chl)/1000)) %>% 
  mutate(Btlag = startbt)


data2 <- prod1 %>% 
  select(algae_conc2, mean_pp, mean_sp, Bt) %>% 
  group_by(algae_conc2) %>% 
  mutate(Btlag = lag(Bt),
         pplag = lag(mean_pp)) %>% 
  filter(day == 22)

obs <- cc %>% 
  select(algae_conc2, coreid, day, live_tt) %>% 
  full_join(cm %>% 
              filter(species_name == "tt") %>% 
              group_by(coreid) %>% 
              summarise(wt = meanna(weight(body_size)))) %>% 
  mutate(Bt = ifelse(live_tt ==0, 0, wt*live_tt)) %>%  #biomass in mg
  left_join(chl %>%
              left_join(meta %>% select(algae, algae_conc2)) %>% 
              group_by(algae_conc2) %>% 
              summarise(chl = mean(chl))) %>% 
  full_join(nep %>% select(midge, coreid, day, algae_conc2, gpp)) %>% 
  mutate(coreid = as.numeric(coreid)) %>% 
  arrange(coreid, day)

# write_csv(data, "data.csv")
# write_csv(data2, "data2.csv")
# write_csv(obs, "microcosm_biomass.csv")

#===
growth2 <- nep %>% 
  group_by(day, midge, algae_conc2) %>% 
  mutate(prod = (gpp*18)*(12/32*pq)) %>% #daily production of algae in g C m^2d^-1
  add_count() %>% 
  summarise(n = unique(n),
            mean_pp = mean(prod),
            sd_pp = sd(prod)/sqrt(n)) %>% 
  full_join(production_boot2 %>% 
              # filter(Pd>0) %>% 
              mutate(g = ((g*0.5)/1000)/ftube_area,
                     midge = "Midges") %>%  #0.5 g C/ g AFDM
              group_by(day, algae_conc2, midge) %>% 
              summarise(mean_g = mean(gd),
                        sd_g = sd(gd))) %>% 
  filter(day!=0,
         midge == "Midges") 



growth1 <- nep %>% 
  group_by(day, midge, algae_conc2) %>% 
  mutate(prod = (gpp*18)*(12/32*pq)) %>% #daily production of algae in g C m^2d^-1
  add_count() %>% 
  summarise(n = unique(n),
            mean_pp = mean(prod),
            sd_pp = sd(prod)/sqrt(n)) %>% 
  full_join(production_boot %>% 
              # filter(Pd>0) %>% 
              mutate(g = ((g*0.5)/1000)/ftube_area,
                     midge = "Midges") %>%  #0.5 g C/ g AFDM
              group_by(day, algae_conc2, midge) %>% 
              summarise(mean_g = mean(gd),
                        sd_g = sd(gd))) %>% 
  filter(day!=0,
         midge == "Midges") 

growth1 %>% 
  mutate(day = ifelse(day == 14, "Day 0 - 14", "Day 14 - 22")) %>% 
  ggplot(aes(x = mean_pp, y = mean_g*1000))+
  facet_wrap(~day, nrow = 1)+
  geom_errorbarh(aes(xmin = mean_pp-sd_pp, xmax = mean_pp+sd_pp), height = NA)+
  geom_errorbar(aes(ymin = mean_g*1000-sd_g*1000, ymax = mean_g*1000+sd_g*1000), width  = NA)+
  geom_point(aes(fill = algae_conc2), shape = 21, size = 2)+
  viridis::scale_fill_viridis(trans = "log", breaks = c(0.01, 0.1, 1))+
  coord_cartesian(xlim = c(0, NA))+
  labs(x = expression("Primary Production g C"~(m^{-2}~d^{-1})),
       y = expression("Growth"~(mu~g~ind^{-1}~d^{-1})),
       fill = "Sediment Treatment")



growth1 %>% 
  mutate(day = ifelse(day == 14, "Day 0 - 14", "Day 14 - 22")) %>% 
  ggplot(aes(x = mean_pp, y = mean_g*1000))+
  facet_wrap(~day, nrow = 1)+
  geom_errorbarh(aes(xmin = mean_pp-sd_pp, xmax = mean_pp+sd_pp), height = NA)+
  geom_errorbar(aes(ymin = mean_g*1000-sd_g*1000, ymax = mean_g*1000+sd_g*1000), width  = NA)+
  geom_point(aes(fill = algae_conc2), shape = 21, size = 2)+
  viridis::scale_fill_viridis(trans = "log", breaks = c(0.01, 0.1, 1))+
  coord_cartesian(xlim = c(0, NA),
                  ylim = c(0, NA))+
  labs(x = expression("Primary Production g C"~(m^{-2}~d^{-1})),
       y = expression("Growth"~(mu~g~ind^{-1}~d^{-1})),
       fill = "Sediment Treatment")


#====Print all Model outputs to a csv======
#extract LMs (LMERs more complicated)
tidy(g14log) %>% mutate(response = "GPP", day = "14") %>%
  rbind(tidy(g22log) %>% mutate(response = "GPP", day = "22")) %>%
  rbind(tidy(nm14log) %>% mutate(response = "# Midges", day = "14")) %>%
  rbind(tidy(nm22log) %>% mutate(response = "# Midges", day = "22")) %>%
  rbind(tidy(prop14log) %>% mutate(response = "proportion 2nd instar", day = "14")) %>%
  rbind(tidy(prop22log) %>% mutate(response = "proportion 2nd instar", day = "22")) %>%
  mutate(estimate = format(estimate, digits = 2),
         std.error = format(std.error, digits = 2),
         statistic = round(statistic, digits = 2),
         p.value = format(p.value, digits = 2)) 
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
# plot(bl22log2)


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
  labs(x = "Sediment Treatment",
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
