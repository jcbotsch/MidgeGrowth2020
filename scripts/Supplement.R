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
source("scripts/MG_Model.R") #this might take some time

#====read in files====
meta <- read_csv("clean data/MG_meta.csv") %>% 
  mutate(algae_conc2 = ifelse(algae_conc == 0, 2e-3, algae_conc)) #half of lowest value
chl <- read_csv("clean data/MG_chl.csv")
om <- read_csv("clean data/MG_om.csv")
hobo <- read_csv("clean data/MG_hobo.csv")
inc_timing <- read_csv("clean data/MG_inc_timing.csv")
cm_raw <- read_csv("clean data/MG_cm.csv")
ltreb <- read_csv("clean data/ltreb_measurements.csv") # https://doi.org/10.6084/m9.figshare.14095697.v3

#====create data frame to predict model response variables ====
topredict <- data.frame(algae_conc2 = rep(unique(meta$algae_conc2), 4),
                        midge = rep(rep(c("Midges", "No Midges"), each = 10), 2),
                        box = rep(c("1", "2"), each = 20)) %>% 
  mutate(midge = fct_relevel(midge, "No Midges", "Midges"))

#model matrix used for all regressions
mm <- model.matrix(~log(algae_conc2)*midge + box, topredict) 
mm[,"box2"] <- 0.5
mm <- unique(mm)


#====Figure S2====
#plot organic content
omplot <- om %>% 
  left_join(meta %>% select(algae, algae_conc2)) %>% 
  unique() %>% 
  ggplot(aes(algae_conc2, perc.org))+
  geom_point(size =2, alpha = 0.6)+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, NA))+
  labs(x = "Initial Algal Abundance",
       y = "Organic Content")

#plot chlorophyll a
corrchlplot <- chl %>% 
  left_join(meta %>% select(algae, algae_conc2)) %>% 
  unique() %>% 
  ggplot(aes(x = algae_conc2, y = chl/1000))+
  geom_point(size = 2, alpha = 0.6)+
  scale_y_continuous(labels = scales::comma_format())+
  labs(x = "",
       y = expression("Chlorophyll mg"~L^{-1})) +  
  theme(axis.text.x = element_blank())

fig1 <- plot_grid(corrchlplot, omplot, nrow  = 2, align = "v")
plot(fig1)

# ggpreview(plot = fig1, dpi = 650, width = 3, height = 4, units = "in")
# ggsave(plot = fig1, filename = "Botsch_MG_Fig1.pdf", dpi = 650, device = "pdf", width = 80, height = 80, units = "mm")

#====Figure S3: Air temperature and Lux====
env <- hobo %>% 
  gather(var, val, -date_time) %>%
  mutate(var = ifelse(var == "temp", "Air Temperature (°C)",
                      "Lux")) %>% 
  ggplot()+
  facet_wrap(~var, scales = "free_y", strip.position = "left", ncol = 1)+
  geom_rect(aes(xmin = dttm_start, xmax =dttm_end, fill = dark_light, ymin = -Inf, ymax = Inf), 
            data = inc_timing,
            alpha = 0.8)+
  geom_line(aes(x = date_time, y = val))+
  theme(strip.placement = "outside",
        legend.position = "none")+
  labs(x = "Date",
       y = NULL,
       fill = "")+
  scale_fill_viridis_d()

plot(env)
# ggpreview(plot = env, dpi = 650, width = 5, height = 3, units = "in")

#====Figure S3: Routine Instars====
labels <- data.frame(breaks = mg_instars) %>% 
  mutate(lag = lag(breaks)) %>% 
  rowwise() %>% 
  mutate(pos = mean(c(breaks, lag))) %>% 
  ungroup %>% 
  mutate(instar = c( NA, "I", "II", "III", "IV"),
         pos = ifelse(instar == "I", 65, pos))


routine_instars <- ltreb %>% 
  filter(species_name == "tt") %>% 
  mutate(head_size_um = head_size/head_size_units_to_mm * 1000) %>%
  filter(head_size_um<max(mg_instars)) %>% 
  ggplot(aes(head_size_um))+
  geom_histogram(bins = 40)+
  geom_vline(xintercept = mg_instars[2:4])+
  geom_text(aes(label = instar, x = pos, y = 3200), data = labels)+
  labs(x = "Head Width (\u03BCm)",
       y = "count")

ggpreview(plot = routine_instars, dpi = 650, width = 3, height = 3, units = "in")

#====Figure S4: Midge Development====
#prep data
cm <- cm_raw %>% 
  left_join(meta) %>% 
  mutate(coreid = as.character(coreid),
         box = as.character(box),
         body_size = body_size/1000) #convert from um to mm

#find starting midge stages
startingprop <- {cm %>% 
    #only tanytarsini
    filter(species_name == "tt",
           day == 0) %>%  
    add_count(coreid, name = "measured") %>% 
    group_by(coreid, day, midge, algae_conc2, box, instar, measured) %>% 
    count(instar, name = "s2") %>% 
    filter(instar == 2) %>% 
    mutate(prop = s2/measured)}$prop

#proportion of larvae in instar
p2 <- cm %>% 
  filter(species_name == "tt",
         day!= 0) %>%  
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

#fit model
prop14log <- glm(cbind(s2, others) ~ log(algae_conc2)*midge+box, 
                 data = p314, 
                 family = quasibinomial())

summary(prop14log) #Type III
summary(update(prop14log, .~. -log(algae_conc2):midge)) #Type II

#subset for day 22
p322 <- p2 %>% 
  filter(day == 22) %>% 
  mutate(midge = factor(midge, levels = c("No Midges", "Midges")))

#fit model
prop22log <- glm(cbind(s2, others) ~log(algae_conc2)*midge+box, 
                 data = p322, 
                 family = quasibinomial())

summary(prop22log) #Type III
summary(update(prop22log, .~. -log(algae_conc2):midge)) #Type II

#get predictions
predicted14 <- mod_predict(mm, prop14log) %>% 
  mutate(day = "Day 14")
predicted22 <- mod_predict(mm, prop22log) %>% 
  mutate(day = "Day 22")

#combine predictions
propredict <- rbind(predicted14, predicted22) %>% 
  mutate(prop = exp(estimate),
         upper = exp(estimate+se),
         lower = exp(estimate-se))

#====Figure S5====
prop_fig <- p2 %>% 
  mutate(day = paste("Day", day)) %>% 
  ggplot(aes(y = prop, x = algae_conc2))+
  geom_hline(yintercept = startingprop, alpha = 0.75, color = "orange")+
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = midge, linetype = midge), color = "black", size = 0.2, alpha = 0.5, data = propredict, show.legend = FALSE)+
  geom_point(aes(fill = midge), shape = 21, alpha = 0.75, size = 2, position = position_jitterdodge(dodge.width = 0.1, jitter.width = 0.1, jitter.height = 0))+
  geom_line(aes(linetype = midge), data = propredict)+
  facet_wrap(~day, nrow = 2)+
  labs(x = "Initial Algal Abundance",
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
# ggpreview(plot = prop_fig, dpi = 650, width = 80, height = 80, units = "mm")
# ggsave(plot = last_plot(), filename = "Botsch_MG_Fig4.pdf", device = "pdf", dpi = 650, width = 80, height = 120, units = "mm")

#====Figure S6====
#simulate
sim <- data.frame(x.init = unique(data$x0),
                  r = exp(par[1]),
                  K = exp(par[2]),
                  a = exp(par[3]),
                  c = exp(par[4])/exp(par[3]), 
                  Tmax = 22) %>% 
  crossing(y.init = unique(data$y0)) %>% 
  mutate(paramset = 1:n())

#join with data
obsmod <- trajectory(sim) %>% 
  left_join(sim) %>% 
  left_join(init.data %>% 
              select(-coreid) %>% 
              unique() %>% 
              rename(x.init = x0, y.init=y0) %>% 
              mutate(midge = ifelse(y.init>0, "Midges", "No Midges"))) %>% 
  convert_to_exp_units()

#plot
mod_data <- obsmod %>% 
  ggplot(aes(x = t, y = val, col = var, linetype = midge, shape = midge))+
  facet_wrap(~round(algae_conc2, 3), ncol = 3)+
  geom_line(aes(y = wt*500, col = "Midge"))+
  geom_line(aes(y = gpp, col = "GPP"))+
  geom_point(aes(y = val, color = taxon, x = day), alpha = 0.75, position = position_dodge(width = 2),  data = data %>% rename(Midge = wt, GPP = gpp) %>% mutate(Midge = 500*Midge) %>%   gather(taxon, val, Midge, GPP)%>% filter(!(taxon == "Midge" & midge == "No Midges")))+
  geom_point(aes(y = val, color = taxon, x = day, alpha = ifelse(taxon == "GPP", "a", "m")),position = position_dodge(width = 2),  data = data  %>% mutate(day = 0, Midge = y0*meany*500, GPP = (x0*gpp.rate)*meanx) %>% select(algae_conc2, midge, day, Midge, GPP) %>% unique() %>%  gather(taxon, val, Midge, GPP)%>% filter(!(taxon == "Midge" & midge == "No Midges")))+
  scale_linetype_manual(values = c("solid", "dotted"))+
  scale_alpha_manual(values = c(0.25, 0.75), guide = "none")+
  scale_shape_manual(values = c(16,21))+
  scale_y_continuous(sec.axis = sec_axis(trans = ~./500, name = "Average Midge Mass (mg)"))+
  scale_color_manual(values = c("black", "red"), guide = "none")+
  labs(linetype = element_blank(),
       color = element_blank(),
       shape = element_blank(),
       x = "Time",
       y = expression("GPP "~(mg~O[2]~m^{-2}~hr^{-1})))+
  theme(axis.title.y.right = element_text(color = "red"),
        axis.text.y.right = element_text(color = "red"),
        axis.ticks.y.right = element_line(color = "red"))

# ggpreview(plot = mod_data, dpi = 650, width = 5, height = 6, units = "in")

#====Figure S7====
arange14 <- obsmoda %>% 
  filter(t == 14,
         y.init>0) %>% 
  ggplot(aes(x = a, y = gdc, fill = algae_conc2, group = x.init))+
  geom_vline(xintercept = aest)+
  geom_line(aes(color = algae_conc2), size = 1)+
  labs(x = "Consumption Rate",
       y = expression("Midge Growth \u03BCg C "~ind^{-1}~d^{-1}),
       color = "Initial Algal Abundance")+
  guides(color = guide_colorbar(title.position = "top", title.hjust = 0.5, barheight = 0.5, barwidth = 8))+
  algae_color
# ggpreview(plot = arange14, width = 3, height = 4, units = "in", dpi = 650)
