#=====load packages====
library(tidyverse)
library(cowplot)
source("scripts/MG_Functions.R")
library(gridExtra)
library(grid)
library(shades)

options(dplyr.summarise.inform = FALSE)



#############################################################
#====Read Data===
meta <- read_csv("clean data/MG_meta.csv") %>% 
  mutate(algae_conc2 = ifelse(algae_conc == 0, 2e-3, algae_conc)) #half of lowest value
chl_raw <- read_csv("clean data/MG_chl.csv")
ep_raw <- read_csv("clean data/MG_ep.csv")
cm_raw <- read_csv("clean data/MG_cm.csv")
cc_raw <- read_csv("clean data/MG_cc.csv")

#====Prepare Data====
#chlorophll data
chl <- chl_raw %>% 
  left_join(meta) %>% 
  group_by(algae_conc2) %>% 
  summarise(chl = mean(chl, na.rm = TRUE)/1000) 

#primary production
gpp <- ep_raw %>% 
  select(coreid, day, gpp) %>% 
  mutate(gpp = gpp*1000) #convert GPP from mgm-2h-1 to gm-2h-1


#midge data
midge <- cc_raw %>% 
  select(coreid, day, live_tt) %>%
  left_join(meta %>% select(coreid, algae_conc2)) %>% 
  full_join(cm_raw %>% 
              filter(species_name == "tt") %>% 
              group_by(coreid) %>% 
              summarise(wt = meanna(weight(body_size/1000)))) %>% 
  arrange(day, coreid) %>% 
  mutate(Bt = ifelse(live_tt ==0, 0, wt*live_tt)) 

#Join 
data <- gpp %>% 
  full_join(meta %>% 
              select(coreid, midge, algae_conc2)) %>% 
  full_join(chl) %>% 
  left_join(midge) %>% 
  select(algae_conc2, midge, day, coreid, everything())

# midge weights at time 0
init.midges <- midge %>% 
  filter(day == 0)
y.init.mean <- {init.midges %>% 
    summarise(y = sum(live_tt*wt)/sum(live_tt))}$y

# prune to experimental falcon tubes
data <- data %>% 
  filter(day>0) %>% 
  arrange(day, coreid)

# standardize variables
data <- data %>% 
  mutate(meanx = mean(gpp, na.rm = TRUE),
         meanx0 = mean(chl, na.rm = TRUE),
         meany = mean(wt, na.rm  = TRUE),
         x = gpp/meanx,
         x0 = chl/meanx0,
         y = wt/meany,
         y0 = ifelse(midge == "Midges", y.init.mean/mean(wt, na.rm = TRUE), 0)) %>% 
  arrange(coreid)

#====estimate parameters====
# start values
r <- .5
K <- .9
a <- .2
ac <- .03
# m <- 0
gpp.rate <- 10^0

#get initial data
init.data <- data %>% 
  select(coreid, algae_conc2, x0, y0, contains("mean")) %>% 
  unique()

x.init <- init.data$x0
y.init <- init.data$y0

SS.weight <- 30

par.full <- log(c(r, K, a, ac, gpp.rate))
par.fixed <- c(NA, NA, NA, NA, NA)
par <- par.full[is.na(par.fixed)]

#====Fit====
set.seed(12345)
SSmin <- 10^10
nrep <- 10

for (i.rep in 1:nrep){
  opt <- optim(par = par, fn=SSfit, data=data, x.init = x.init, y.init = y.init, SS.weight = SS.weight, par.fixed=par.fixed, method = "SANN")
  opt <- optim(par = opt$par, fn=SSfit, data=data, x.init = x.init, y.init = y.init, SS.weight = SS.weight, par.fixed=par.fixed, method = "Nelder-Mead")
  if(opt$value < SSmin){
    SSmin <- opt$value
    opt.opt <- opt
  }
  par <- exp(rnorm(n=length(par), sd=1)) * opt.opt$par
  show(c(opt$value, SSfit(par=opt$par, data=data, x.init = x.init, y.init=y.init, SS.weight = SS.weight, par.fixed=par.fixed, tofit=F)))
}
par.temp <- par.fixed
par.temp[is.na(par.fixed)] <- opt.opt$par
par <- par.temp

#extract fitted parameters
r <- exp(par[1])
K <- exp(par[2])
a <- exp(par[3])
ac <- exp(par[4])
c <- ac/a
# m <- exp(par[5])
# gpp.rate <- exp(par[6])
gpp.rate <- exp(par[5])

round(c(r=r, K=K, a=a, c=c, m=m, gpp.rate=gpp.rate, SS=opt.opt$value), digits=4)
# r        K        a        c        m       gpp.rate       SS 
# 0.1725   0.9290   0.0447   0.1788   0.0015 0.2088       66.8795 

#====Simulate to compare with observations====
#simulate model as is 
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
              mutate(midge = ifelse(y.init>0, "Midges", "No Midges")))
#plot
mod_data <- obsmod %>% 
  mutate(gpp = (x*gpp.rate)*unique(data$meanx),
         wt = y*unique(data$meany)) %>% 
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

obsmod %>% 
  filter(midge == "Midges") %>% 
  mutate(gpp = (x*gpp.rate)*unique(data$meanx),
         wt = y*unique(data$meany)) %>% 
  ggplot(aes(x = t, y = val, col = var, linetype = midge, shape = midge))+
  facet_wrap(~round(algae_conc2, 3), ncol = 3)+
  geom_line(aes(y = wt*500, col = "Midge"))+
  geom_line(aes(y = gpp, col = "GPP"))+
  # geom_point(aes(y = val, color = taxon, x = day), alpha = 0.75, position = position_dodge(width = 2),  data = data %>% rename(Midge = wt, GPP = gpp) %>% mutate(Midge = 500*Midge) %>%   gather(taxon, val, Midge, GPP)%>% filter(!(taxon == "Midge" & midge == "No Midges")))+
  # geom_point(aes(y = val, color = taxon, x = day, alpha = ifelse(taxon == "GPP", "a", "m")),position = position_dodge(width = 2),  data = data  %>% mutate(day = 0, Midge = y0*meany*500, GPP = (x0*gpp.rate)*meanx) %>% select(algae_conc2, midge, day, Midge, GPP) %>% unique() %>%  gather(taxon, val, Midge, GPP)%>% filter(!(taxon == "Midge" & midge == "No Midges")))+
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

#====Varying attack rates====
#parameters to simulate
sima <- data.frame( a = c(seq(0, 0.4, by = 0.001), a),
                    K = K,
                    r = r,
                    c = c, 
                    Bm = 0,
                    Tmax = 22) %>% 
  crossing(x.init = unique(data$x0)) %>% 
  crossing(y.init = unique(data$y0)) %>% 
  mutate(paramset = 1:n())

#simulate
obsmoda <- trajectory(sima) %>% 
  left_join(sima) %>% 
  left_join(init.data %>% 
              select(-coreid) %>% 
              unique() %>% 
              rename(x.init = x0, y.init=y0) %>% 
              mutate(midge = ifelse(y.init>0, "Midges", "No Midges")))

aest <- a
rest <- r

#plot 
arange <- obsmoda %>% 
  filter(t == Tmax,
         y.init>0) %>% 
  convert_to_exp_units() %>% 
  ggplot(aes(x = a, y = gdc, fill = algae_conc2, group = x.init))+
  geom_vline(xintercept = aest)+
  geom_line(size = 1.1)+
  geom_line(aes(color = algae_conc2))+
  labs(x = "Consumption Rate",
       y = expression("Midge Growth \u03BCg C "~ind^{-1}~d^{-1}),
       color = "Sediment Treatment")+
  guides(color = guide_colorbar(title.position = "top", title.hjust = 0.5))+
  scale_color_viridis_c(trans = "log", breaks = c(0.01, 0.1, 1))

# ggpreview(plot = arange, width = 3, height = 4, units = "in", dpi = 650)
# ggpreview(plot = arange2, width = 3, height = 5, units = "in", dpi = 650)

arange14 <- obsmoda %>% 
  filter(t == 14,
         y.init>0) %>% 
  convert_to_exp_units() %>% 
  ggplot(aes(x = a, y = gdc, fill = algae_conc2, group = x.init))+
  geom_vline(xintercept = aest)+
  geom_line(size = 1.1)+
  geom_line(aes(color = algae_conc2))+
  labs(x = "Consumption Rate",
       y = expression("Midge Growth \u03BCg C "~ind^{-1}~d^{-1}),
       color = "Sediment Treatment")+
  guides(color = guide_colorbar(title.position = "top", title.hjust = 0.5))+
  scale_color_viridis_c(trans = "log", breaks = c(0.01, 0.1, 1))

# ggpreview(plot = arange14, width = 3, height = 4, units = "in", dpi = 650)


##Plot Trajectories of Optimal Attack Rates
optima <- obsmoda %>% 
  filter(t == Tmax,
         y.init>0,
         y-y.init == max(y-y.init)) %>%
  pull(a)


optims <- obsmoda %>% 
  filter(t == Tmax,
         y.init>0) %>% 
  group_by(x.init) %>% 
  filter(y-y.init == max(y-y.init)) %>% 
  select(x.init, a)

#plot trajectories 
obsmoda %>% 
  filter(a == optima,
         y.init != 0) %>% 
  convert_to_exp_units() %>% 
  mutate(wt = wt*200) %>% 
  gather(var, val, gpp, wt) %>% 
  ggplot(aes(x = t, y = val, color = var))+
  facet_wrap(~round(init.chl, digits = 2))+
  geom_line()+
  lims(y = c(0, NA))+
  labs(y = "Biomass",
       color = element_blank(),
       linetype = element_blank())

# ggpreview(plot = last_plot(), width  = 6, height = 6.5, units = "in", dpi = 650)

#====Optimal attack rates under different resource growth rates ====

simr1 <- data.frame( a = a,
                    K = K,
                    r =  c(seq(0, 0.5, by = 0.001), rest),
                    c = c,
                    Bm = 0,
                    Tmax = 22) %>% 
  crossing(x.init = unique(data$x0)) %>% 
  crossing(y.init = unique(data$y0)) %>% 
  mutate(paramset = 1:n())

#simulate
obsmodr <- trajectory(simr1) %>% 
  left_join(simr1) %>% 
  left_join(init.data %>% 
              select(-coreid) %>% 
              unique() %>% 
              rename(x.init = x0, y.init=y0) %>% 
              mutate(midge = ifelse(y.init>0, "Midges", "No Midges")))



curve1 <- obsmodr %>% 
  filter(t == Tmax,
         y.init>0) %>% 
  convert_to_exp_units() %>% 
  gather(var, val, gdc, gpp) %>% 
  mutate(var = ifelse(var == "gdc", 
                      "Midge~Growth~Rate~ind^{-1}",
                      "Primary~Production~cm^{-2}")) %>% 
  ggplot(aes(x = r, y = val, fill = algae_conc2, group = x.init))+
  facet_wrap(~var, scales = "free_y", labeller = label_parsed)+
  geom_vline(xintercept = c(0.5*r, 1.5*r, 2.5*r, r))+
  geom_line(size = 1.1)+
  geom_line(aes(color = algae_conc2))+
  labs(x = "Resource Growth Rate",
       y = expression("\u03BCg C "~d^{-1}),
       color = "Sediment Treatment")+
  guides(color = guide_colorbar(title.position = "top", title.hjust = 0.5))+
  scale_color_viridis_c(trans = "log", breaks = c(0.01, 0.1, 1))


curve2 <- obsmoda %>% 
  filter(t == Tmax,
         y.init>0) %>% 
  convert_to_exp_units() %>% 
  gather(var, val, gdc, gpp) %>% 
  mutate(var = ifelse(var == "gdc", 
                      "Midge~Growth~Rate~ind^{-1}",
                      "Primary~Production~cm^{-2}")) %>% 
  ggplot(aes(x = a, y = val, group = x.init))+
  facet_wrap(~var, scales = "free_y", labeller = label_parsed) +
  geom_vline(xintercept = c(a, 0.2))+
  geom_line(size = 1.1)+
  geom_vline(aes(xintercept = a, color = algae_conc2), data = optims %>% left_join(init.data %>% select(algae_conc2, x.init = x0)))+
  geom_line(aes(color = algae_conc2))+
  labs(x = "Consumption Rate",
       y = expression("\u03BCg C "~d^{-1}),
       color = "Sediment Treatment")+
  guides(color = guide_colorbar(title.position = "top", title.hjust = 0.5))+
  scale_color_viridis_c(trans = "log", breaks = c(0.01, 0.1, 1))

grid.arrange(curve1+theme(legend.position = "none"), curve2 + theme(legend.position = "none"))

#====Plot biomasses with different attack rates====
plota <- obsmoda %>% 
  filter(t %in% c(Tmax),
         y.init>0,
         a %in% c(aest, 0.25)) %>% 
  bind_rows(obsmoda %>% filter(t == 22, y.init != 0) %>% right_join(optims)) %>%
  mutate(al = ifelse(a == aest, 
                     "a == a[est]", 
      ifelse(a == 0.25, paste("a==", a), "a == a[opt]")),
         al = fct_reorder(al, a)) %>% 
  convert_to_exp_units() %>% 
  ggplot(aes(y = gpp, x = gdc, fill = algae_conc2))+
  facet_wrap(~al, scales = "free", labeller = label_parsed, ncol = 1)+
  geom_path(aes(group = t))+
  labs(y = expression("Primary Production \u03BCg C"~cm^{-2}~d^{-1}),
       x = expression("Midge Growth \u03BCg C "~ind^{-1}~d^{-1}),
       fill = "Sediment Treatment",
       shape = element_blank())+
  geom_point(shape = 22, size = 2)+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, barheight = 0.5, barwidth = 8))+
  scale_fill_viridis_c(trans = "log", breaks = c(0.01, 0.1, 1))+
  scale_x_continuous(n.breaks= 4)+
  scale_y_continuous(n.breaks = 5)

plotas <- obsmoda %>% 
  filter(t %in% c(14),
         y.init>0,
         a %in% c(aest, 0.25)) %>% 
  bind_rows(obsmoda %>% filter(t == 14, y.init != 0) %>% right_join(optims)) %>%
  mutate(al = ifelse(a == aest, 
                     "a == a[est]", 
                     ifelse(a == 0.25, paste("a==", a), "a == a[opt]")),
         al = fct_reorder(al, a)) %>% 
  ggplot(aes(y = x, x = y, fill = x.init))+
  facet_wrap(~al, scales = "free", labeller = label_parsed, ncol = 1)+
  geom_path(aes(group = t))+
  labs(y = "Resource Biomass",
       x = "Consumer Biomass",
       fill = "Initial Resource Biomass",
       shape = element_blank())+
  geom_point(shape = 21, size = 2)+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, barheight = 0.5, barwidth = 8))+
  scale_fill_viridis_c(trans = "log", breaks = c(0.25, 1,4))+
  scale_x_continuous(n.breaks= 4)+
  scale_y_continuous(n.breaks = 5)
# ggpreview(plot = last_plot(), width = 3, height = 4, units = "in", dpi = 650)

atraj <-obsmoda %>% 
  filter(y.init>0,
         a %in% c(aest, 0.25),
         x.init %in% c(min(x.init), max(x.init))) %>% 
  bind_rows(obsmoda %>% filter(y.init != 0, x.init %in% c(min(x.init), max(x.init))) %>% right_join(optims %>% ungroup %>% filter(x.init %in% c(min(x.init), max(x.init))))) %>%
  mutate(al = ifelse(a == aest, 
                     "a == a[est]", 
                     ifelse(a == 0.25, paste("a==", a), "a == a[opt]")),
         al = fct_reorder(al, a)) %>% 
  convert_to_exp_units() %>% 
  # filter(x.init == min(x.init)) %>% 
  ggplot(aes(x = t, y = val, linetype = var, group = interaction(x.init), shape = midge, col = x.init))+
  facet_wrap(~al, labeller = label_parsed, ncol = 1)+
 
  # geom_line(aes(y = wt*200), linetype = "solid", alpha = 0.5, color = "black", size = 1.1)+
  # geom_line(aes(y = gpp), linetype = "solid", alpha = 0.5, color = "black", size = 1.1)+
  geom_line(aes(y = wt*200, linetype = "Midge"))+
  geom_line(aes(y = gpp, linetype = "GPP"))+
  scale_linetype_manual(values = c("solid", "dashed"))+
  scale_alpha_manual(values = c(0.25, 0.75), guide = "none")+
  scale_shape_manual(values = c(16,21))+
  scale_y_continuous(sec.axis = sec_axis(trans = ~./200, name = "Average Midge Mass (mg)"))+
  # scale_color_manual(values = c("black", "red"), guide = "none")+
  scale_color_viridis_c(guide = "none")+
  labs(linetype = element_blank(),
       color = element_blank(),
       shape = element_blank(),
       x = "Day",
       y = expression("GPP "~(mg~O[2]~m^{-2}~hr^{-1})))+
  theme(strip.placement = "outside")

atraj
  
#trajectories
obsmoda %>% 
  filter(y.init>0,
         a %in% c(aest, optima, 0.2)) %>% 
  rename(consumer = y, resource = x) %>% 
  filter(x.init %in% unique(.$x.init)[c(1,6,10)]) %>%
  gather(var, val, consumer, resource) %>%
  mutate(al = ifelse(a == aest, 
                     "a == a[est]", 
                     ifelse(a == 0.2, paste("a==", a), "a == a[opt]")),
         al = fct_reorder(al, a)) %>% 
  mutate(val = ifelse(var == "consumer", val, val)) %>% 
  ggplot(aes(x = t, y = val))+
  facet_grid(al~paste( "x0 == ", round(x.init,1)), labeller = label_parsed, scales = "free_y")+
  geom_line(aes(col = var, group = interaction(var, x.init)))+
  # geom_point(aes(fill = x.init), shape = 21)+
  labs(x = "t",
       y = "Biomass",
       fill = "Initial Resource Biomass",
       color = "")+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))+
  scale_fill_viridis_c(trans = "log", breaks = c(0.25, 1,4))

# ggpreview(plot = last_plot(), width = 6, height  = 6, units = "in", dpi = 650)


plot_grid(plota, atraj+theme(axis.title = element_blank()))

ggpreview(plot = last_plot(), width = 3, height  = 6, units = "in", dpi = 650)

#====Varying Resource Growth Rates====
#set ranges 
rrange <- c(0.5, 1.5, 3)

#parameters to simulate
simr <- data.frame( a = aest,
                    K = K,
                    rrange = rrange,
                    r = rrange*rest,
                    c = c, 
                    
                    Bm = 0,
                    Tmax = 22) %>% 
  crossing(x.init = unique(data$x0)) %>% 
  crossing(y.init = max(data$y0)) %>% 
  mutate(paramset = 1:n())

#simulate
rmod <- trajectory(simr) %>% 
  left_join(simr) %>% 
  left_join(init.data %>% rename(x.init = x0, y.init = y0)) %>% 
  convert_to_exp_units()

#====Plot Biomasses across resource growth rates====
plot1 <- rmod %>% 
  filter(t%in% c(Tmax)) %>% 
  ggplot(aes(y = gpp, x = gdc, fill = algae_conc2))+
  facet_wrap(~paste("r = ", rrange, "\u00D7 r"), scales = "free", ncol = 1)+
  geom_path(aes(group  = t))+
  labs(y = expression("Primary Production \u03BCg C"~cm^{-2}~d^{-1}),
       x = expression("Midge Growth \u03BCg C "~ind^{-1}~d^{-1}),
       fill = "Sediment Treatment",
       shape = element_blank())+
  geom_point(shape = 22, size = 2)+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, barheight = 0.5, barwidth = 8))+
  scale_fill_viridis_c(trans = "log", breaks = c(0, 0.01,1))+
  scale_x_continuous(n.breaks= 4)+
  scale_y_continuous(n.breaks = 5)

rtraj <- rmod %>% 
  filter(x.init %in% c(min(x.init), max(x.init))) %>% 
  ggplot(aes(x = t, group = interaction(x.init), color = algae_conc2))+
  facet_wrap(~paste("r = ", rrange, "\u00D7 r"), scales = "free_y", ncol = 1)+
  
  # geom_line(aes(y = wt*200), linetype = "solid", alpha = 0.5, color = "black", size = 1.1)+
  # geom_line(aes(y = gpp), linetype = "solid", alpha = 0.5, color = "black", size = 1.1)+
  geom_line(aes(y = wt*200, linetype = "Midge"))+
  geom_line(aes(y = gpp, linetype = "GPP"))+
  scale_linetype_manual(values = c("solid", "dashed"))+
  scale_alpha_manual(values = c(0.25, 0.75), guide = "none")+
  scale_shape_manual(values = c(16,21))+
  scale_y_continuous(sec.axis = sec_axis(trans = ~./200, name = "Average Midge Mass (mg)"))+
  # scale_color_manual(values = c("black", "red"), guide = "none")+
  scale_color_viridis_c(guide = "none")+
  labs(linetype = element_blank(),
       color = element_blank(),
       shape = element_blank(),
       x = "Day",
       y = expression("GPP "~(mg~O[2]~m^{-2}~hr^{-1})))+
  theme(strip.placement = "outside")

#T = 14
plot1s <- rmod %>% 
  filter(t%in% c(14)) %>% 
  ggplot(aes(y = x, x = y, fill = x.init))+
  facet_wrap(~paste("r = ", rrange, "\u00D7 r"), scales = "free", ncol = 1)+
  geom_path(aes(group  = t))+
  labs(y = "Resource Biomass",
       x = "Consumer Biomass",
       fill = "Initial Resource Biomass")+
  geom_point(shape = 22, size = 2)+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, barheight = 0.5, barwidth = 8))+
  scale_fill_viridis_c(trans = "log", breaks = c(0.25, 1,4))+
  scale_x_continuous(n.breaks= 4)+
  scale_y_continuous(n.breaks = 5)


#====Simulate optimal attack rate using different resource growth rates====
ar <- data.frame( a = c(seq(0, 0.4, by = 0.001), a),
                  K = K,
                  y.init = max(data$y0),
                  c = c, 
                  
                  Bm = 0,
                  Tmax = 22) %>% 
  crossing(x.init = unique(data$x0)) %>% 
  crossing(rrange =c( 0.5,1,  1.5, 3)) %>% 
  mutate(r= rest*rrange) %>% 
  mutate(paramset = 1:n())

simar <- trajectory(ar) %>% 
  left_join(ar) %>% 
  left_join(init.data %>% rename(x.init = x0, y.init = y0)) %>% 
  convert_to_exp_units()

arangev3 <- simar %>% 
  filter(t == Tmax,
  ) %>% 
  ggplot(aes(x = a, y = gdc, fill = x.init, group = x.init))+
  facet_wrap(~paste(rrange, "\u00D7 r"), scales = "free_y", ncol = 1)+
  geom_line(size = 1.1)+
  geom_line(aes(color = algae_conc2))+
  geom_vline(xintercept = aest)+
  labs(x = "Consumption Rate",
       y = expression("Midge Growth \u03BCg C "~ind^{-1}~d^{-1}),
       color = "Sediment Treatment")+
  guides(color = guide_colorbar(title.position = "top", title.hjust = 0.5))+
  scale_color_viridis_c(trans = "log", breaks = c(0.01, 0.1, 1))

arangev3.1 <- simar %>% 
  filter(t == Tmax,
  ) %>% 
  ggplot(aes(x = a, y = gdc, fill = x.init, group = x.init))+
  facet_wrap(~paste(rrange, "\u00D7 r"), scales = "free_y")+
  geom_line(size = 1.1)+
  geom_line(aes(color = algae_conc2))+
  geom_vline(xintercept = aest)+
  labs(x = "Consumption Rate",
       y = expression("Midge Growth \u03BCg C "~ind^{-1}~d^{-1}),
       color = "Sediment Treatment")+
  guides(color = guide_colorbar(title.position = "top", title.hjust = 0.5))+
  scale_color_viridis_c(trans = "log", breaks = c(0.01, 0.1, 1))


ggpreview(plot = arangev3.1, width = 3, height = 4, units = "in", dpi = 650)

arangev3.5 <- simar %>% 
  filter(t == Tmax,
  ) %>% 
  ggplot(aes(x = a, y = gppc, fill = x.init, group = x.init))+
  facet_wrap(~paste(rrange, "\u00D7 r"), scales = "free_y", ncol = 1)+
  geom_line(size = 1.1)+
  geom_line(aes(color = algae_conc2))+
  geom_vline(xintercept = aest)+
  labs(x = "Consumption Rate",
       y = expression("Primary Production \u03BCg C "~m^{-2}~d^{-1}),
       color = "Sediment Treatment")+
  guides(color = guide_colorbar(title.position = "top", title.hjust = 0.5))+
  scale_color_viridis_c(trans = "log", breaks = c(0.01, 0.1, 1))

arangexgrob <- textGrob("Consumption Rate")
sedtrtleg <- get_legend(arangev3 + theme(legend.key.height = unit (1, "lines")))


ap1 <- plot_grid(arangev3 +theme(legend.position = "none", axis.title.x = element_blank()), 
          arangev3.5 +theme(legend.position = "none", axis.title.x = element_blank()))

ap2 <- plot_grid(sedtrtleg, ap1, ncol = 1, rel_heights = c(0.15, 0.85))

apf <- grid.arrange(ap2, bottom = arangexgrob)

ggpreview(plot = apf, width = 3, height = 6, units = "in", dpi = 650)


simar %>% 
  # filter(rrange == 2.5,
  #        a<0.05,
  #        a>0.025) %>% 
  filter(t == Tmax)%>% 
  ggplot(aes(x = a, y = gdc, fill = x.init, group = x.init))+
  facet_wrap(~paste(rrange, "\u00D7 Resource Growth Rate"), scales = "free")+
  geom_line(size = 1.1)+
  geom_line(aes(color = algae_conc2))+
  geom_vline(xintercept = aest)+
  labs(x = "Consumption Rate",
       y = expression("Midge Growth \u03BCg C "~ind^{-1}~d^{-1}),
       color = "Sediment Treatment")+
  guides(color = guide_colorbar(title.position = "top", title.hjust = 0.5))+
  scale_color_viridis_c(trans = "log", breaks = c(0.01, 0.1, 1))


#====Varying initial x biomass and resource growth====
#parameters
simequi <- data.frame( a = aest,
                    K = K,
                    r = seq(0, 0.5, by = 0.01),
                    c = c, 
                    
                    Bm = 0,
                    Tmax = 22) %>% 
  crossing(x.init = seq(0,8, by = 0.1)) %>% 
  crossing(y.init = max(data$y0)) %>% 
  mutate(paramset = 1:n())

#simulate
equi2 <- trajectory(simequi) %>% 
  left_join(simequi) %>% 
  bind_cols(init.data %>% select(contains("mean")) %>% unique()) %>% 
  convert_to_exp_units()

#plotting 
nb <- 10 #number of breaks
logbreaks <- round(exp(seq(log(0.1), log(max(equi2$gdc)), length.out = nb)), 2)[2:nb-1]

#plot
equip <- equi2 %>% 
  filter(t == Tmax) %>% 
  ggplot(aes(x = init.chl, y = r))+
  geom_tile(aes(fill = gdc))+
  geom_contour(aes(z = gdc, colour = ..level..), breaks = logbreaks)+
  geom_rug(aes(x = chl), data = data)+
  geom_hline(aes(yintercept = rest), linetype = "dashed")+
  scale_fill_gradient2(high = "black", low = "dodgerblue", mid = "white", midpoint = 0)+
  labs(x = "Initial Chlorophyll-a",
       y = "Resource Growth Rate",
       fill = "Midge\nGrowth")+
  scale_color_viridis_c(option = "plasma", trans = "log")

equipf <- directlabels::direct.label(equip, list("far.from.others.borders", "calc.boxes", "enlarge.box", 
                          box.color = NA, fill = "transparent", "draw.rects", hjust = 1, vjust = 1, cex = 0.75))


nb2 <- 8
logbreaks2 <- round(exp(seq(log(2), log(max(equi2$gpp)), length.out = nb2)), 2)[2:nb2-1]

equip2 <- equi2 %>% 
  filter(t == Tmax) %>% 
  ggplot(aes(x = init.chl, y = r))+
  geom_tile(aes(fill = gpp))+
  geom_contour(aes(z = gpp, col = ..level..), breaks = logbreaks2)+
  geom_rug(aes(x = chl), data = data)+
  geom_hline(aes(yintercept = rest), linetype = "dashed")+
  scale_fill_gradient2(high = "black", low = "dodgerblue", mid = "white", midpoint = 0)+
  labs(x = "Initial Chlorophyll-a",
       y = "Resource Growth Rate",
       fill = "GPP")+
  scale_color_viridis_c(trans = "log", guide = NULL)

equipf2 <- directlabels::direct.label(equip2, list("far.from.others.borders", "calc.boxes", "enlarge.box", 
                                       box.color = NA, fill = "transparent", "draw.rects", hjust = 1, vjust = 1, cex = 0.75))

#====Combine Initial Biomass and resource growth plots====
equipcomb <- plot_grid(equipf + ggtitle("A")+ theme(axis.title = element_blank(), plot.title = element_text(face = "bold"), legend.position = "right"), 
                       equipf2 + ggtitle("B")+ theme(axis.title = element_blank(),plot.title = element_text(face = "bold"), legend.position = "right"), 
                       ncol = 1)

y.grob1 <- textGrob("Resource Growth Rate", rot=90)
x.grob1 <- textGrob("Initial Chlorophyll-a")

equipfinal <- grid.arrange(equipcomb, left = y.grob1, bottom = x.grob1)
ggpreview(plot = equipfinal, width = 3, height = 5, units = "in", dpi = 650)

#====Combine PP and SP plots====
zplot <- plot_grid( plota +theme(axis.title = element_blank(), legend.position = "none"),
                    plot1 + theme(axis.title = element_blank(), legend.position = "none"), 
                    nrow = 1,
                    align = "h")

leg <- get_legend(plota)
plotfin <- plot_grid(leg, zplot, ncol = 1, rel_heights = c(0.15, 0.85))

y.grob <- textGrob(expression("Primary Production \u03BCg C"~cm^{-2}~d^{-1}),
                   rot=90)
x.grob <- textGrob(expression("Midge Growth \u03BCg C "~ind^{-1}~d^{-1}))

comb <- grid.arrange(plotfin, left = y.grob, bottom = x.grob)
ggpreview(plot = comb, width = 4, height = 5, units = "in", dpi = 650)


leg2 <- get_legend(atraj)
tplot <- plot_grid(atraj + theme(axis.title.x = element_blank(), axis.title.y.right = element_blank(), legend.position = "none"),
          rtraj + theme(axis.title.x = element_blank(), axis.title.y.left = element_blank(), legend.position = "none"), align = "h")

x.grob2 <- textGrob("Day")

plott <- plot_grid(leg2, tplot, ncol = 1, rel_heights = c(0.15, 0.85))

combt <- grid.arrange(plott, bottom = x.grob2)


combf <- plot_grid(comb, combt, rel_widths = c(1, 1.2), labels = c("A", "B"), align = "h")

ggpreview(plot = combf, width = 6, height = 4, units = "in", dpi = 650)

zsplot <- plot_grid( plotas +theme(axis.title = element_blank(), legend.position = "none"),
                    # plot2s+ theme(axis.title = element_blank(), legend.position = "none"), 
                    plot1s + theme(axis.title = element_blank(), legend.position = "none"), 
                    labels = c("A", "B", "C"),
                    nrow = 1)

zs2 <- plot_grid(leg, zsplot, ncol = 1, rel_heights = c(0.15, 0.85))


combs <- grid.arrange(zs2, left = y.grob, bottom = x.grob)




# ggpreview(plot = combs, width = 3, height = 5, units = "in", dpi = 650)

