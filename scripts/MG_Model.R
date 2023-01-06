#=====load packages====
library(tidyverse)
library(cowplot)
library(gridExtra)
library(grid)
library(shades)
source("scripts/MG_Functions.R")


options(dplyr.summarise.inform = FALSE)

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

#====fit model====
# starting values for optimization
r <- .5
K <- .9
a <- .2
ac <- .03
gpp.rate <- 10^0

#get initial data
init.data <- data %>% 
  select(coreid, algae_conc2, x0, y0, contains("mean")) %>% 
  unique()

x.init <- init.data$x0
y.init <- init.data$y0

SS.weight <- 30 

par.full <- log(c(r, K, a, ac, gpp.rate))
par.fixed <- c(NA, NA, NA, NA, NA) #estimate all parameters
par <- par.full[is.na(par.fixed)]


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
gpp.rate <- exp(par[5])

#print fitted parameters
round(c(r=r, K=K, a=a, c=c, gpp.rate=gpp.rate, SS=opt.opt$value), digits=4)
# r        K        a        c         gpp.rate       SS 
# 0.1706   0.9326   0.0438   0.1656      0.2061  67.0157


#====Simulate Varying attack rates====
#parameters to simulate
sima <- data.frame( a = c(seq(0, 0.4, by = 0.001), a),
                    K = K,
                    r = r,
                    c = c,
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
              mutate(midge = ifelse(y.init>0, "Midges", "No Midges"))) %>% 
  convert_to_exp_units()

aest <- a
rest <- r


#Pull optima
optims <- obsmoda %>% 
  filter(t == Tmax,
         y.init>0) %>% 
  group_by(x.init) %>% 
  filter(y-y.init == max(y-y.init)) %>% 
  select(x.init, a) %>% 
  ungroup

#====Figure 3==== 
arange <- obsmoda %>% 
  filter(t == Tmax,
         y.init>0) %>% 
  ggplot(aes(x = a, y = gdc, fill = algae_conc2, group = x.init))+
  geom_vline(xintercept = aest)+
  geom_line(aes(color = algae_conc2), size = 1)+
  labs(x = "Consumption Rate",
       y = expression("Midge Growth \u03BCg C "~ind^{-1}~d^{-1}),
       color = "Initial Algal Abundance")+
  guides(color = guide_colorbar(title.position = "top", title.hjust = 0.5, barheight = 0.5, barwidth = 8))+
  algae_color

ggpreview(plot = arange, width = 945, height = 1102, units = "px", dpi = 300)
ggsave(plot = arange, filename = "figures/Botsch_MidgeGrowth_fig3.pdf", device = "pdf", width = 945, height = 1102, units = "px", dpi = 300)

arange_pres <- arange+geom_vline(xintercept = aest, color = "white")+theme_black()
ggpreview(plot = arange_pres, width = 945, height = 1102, units = "px", dpi = 300)


arange_pres <- obsmoda %>% 
  filter(t == Tmax,
         y.init>0) %>% 
  ggplot(aes(x = a, fill = algae_conc2, group = x.init))+
  geom_vline(xintercept = aest)+
  geom_line(aes(y = gdc, color = algae_conc2), size = 1)+
  geom_rug(aes(color = algae_conc2), size = 0.8, data = optims %>% left_join(obsmoda %>% select(x.init, algae_conc2) %>% unique()))+
  # geom_rug(x = aest, size = 1)+
  geom_vline(xintercept = aest, size = 0.8)+
  geom_rug(x = 0.25, size = 0.8)+
  labs(x = "Consumption Rate",
       y = expression("Midge Growth \u03BCg C "~ind^{-1}~d^{-1}),
       color = "Initial Algal Abundance")+
  guides(color = guide_colorbar(title.position = "top", title.hjust = 0.5, barheight = 0.5, barwidth = 8))+
  algae_color


obsmoda %>% 
  filter(t == Tmax,
         x.init == max(x.init),
         y.init>0) %>% 
  ggplot(aes(x = a, y = gdc))+
  # geom_vline(xintercept = aest)+
  geom_line(size = 1)+
  labs(x = "Consumption Rate",
       y = "Consumer Growth Rate")+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())
ggpreview(plot = last_plot(), width = 945, height = 945, units = "px", dpi = 300)


#parameters to simulate
presa <- data.frame( a = c(seq(0, 0.4, by = 0.001), a),
                    K = 0.9,
                    r = 0.5,
                    c = 0.03/0.2,
                    Tmax = 22) %>% 
  crossing(x.init = max(data$x0)) %>% 
  crossing(y.init = max(data$y0)) %>% 
  mutate(paramset = 1:n())

#simulate
pres_a <- trajectory(presa) %>% 
  left_join(presa) %>% 
  left_join(init.data %>% 
              select(-coreid) %>% 
              unique() %>% 
              rename(x.init = x0, y.init=y0) %>% 
              mutate(midge = ifelse(y.init>0, "Midges", "No Midges"))) %>% 
  convert_to_exp_units()

pa <- pres_a %>% filter(t == Tmax, y == max(y)) %>% pull(a)


aselect <- data.frame(a = c(0.015, pa, 0.25), v = 1:3)

plot1 <- pres_a %>% 
  filter(t == Tmax) %>% 
  ggplot(aes(x = a, y = gdc))+
  geom_line(size = 1)+
  geom_point(size = 4, data = . %>% filter(a %in% aselect$a))+
  geom_text(aes(label = v), hjust = -1, vjust = -0.25, size = 4, data = . %>% inner_join(aselect))+
  labs(x = "Consumption Rate",
       y = "Consumer Growth Rate")+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())
  


plot2 <- pres_a %>% 
  inner_join(aselect) %>% 
  mutate(y = y/max(y),
         x = x/max(x)) %>% 
  ggplot(aes(x = t))+
  geom_line(aes(y = y, color = "Consumer"), size = 1)+
  geom_line(aes(y = x, color = "Resource"), linetype = "dashed", size = 1)+
  scale_color_manual(values = c("red", "blue"))+
  geom_text(aes(label = v, x = 2, y = 0.9))+
  facet_wrap(~v, ncol = 1)+
  labs(linetype = element_blank(),
       color = element_blank(),
       x = "Time",
       y = "Biomass")+
  theme(strip.placement = "outside",
        # legend.key.width = unit(0, "lines"),
        # legend.key.height = unit(0.1, "lines"),
        # legend.spacing.y = unit(0.5, "lines"),
        # legend.box = "vertical",
        legend.position = "right",
        strip.text = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank())

plot_grid(plot1, plot2, rel_widths = c(1.2,1))

ggpreview(plot = last_plot(), dpi = 650, width = 7, height = 4)

#====Figure 5a midge growth and pp under different attack rates====
plota <- obsmoda %>% 
  #extract only estimated a, high a (0.25)
  filter(t %in% c(Tmax),
         y.init>0,
         a %in% c(aest, 0.25)) %>% 
  #add optimum attack rates
  bind_rows(obsmoda %>% 
              filter(t == 22, y.init != 0) %>% 
              right_join(optims)) %>%
  #add labels for figure panels
  mutate(al = ifelse(a == aest, 
                     "a == a[est]", 
      ifelse(a == 0.25, paste("a==", a), "a == a[opt]")),
         al = fct_reorder(al, a)) %>% 
  #plot
  ggplot(aes(y = gppc, x = gdc, fill = algae_conc2))+
  facet_wrap(~al, scales = "free", labeller = label_parsed, ncol = 1)+
  geom_path(aes(group = t))+
  labs(y = expression("Primary Production \u03BCg C"~cm^{-2}~d^{-1}),
       x = expression("Midge Growth \u03BCg C "~ind^{-1}~d^{-1}),
       fill = "Initial Algal Abundance",
       shape = element_blank())+
  geom_point(shape = 22, size = 2)+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, barheight = 0.5, barwidth = 8))+
  algae_fill+
  scale_x_continuous(n.breaks= 4)+
  scale_y_continuous(n.breaks = 5)+
  theme(legend.position = "none")


ggpreview(plot = plota, width = 945, height = 1102, units = "px", dpi = 300)


#====Figure 4 biomass trajectories under different attack rates====
atraj <-obsmoda %>% 
  #extract only the estimated and high attack rates, with midges and the highest Xt=0
  filter(y.init>0,
         a %in% c(aest, 0.25),
         x.init %in% c(min(x.init), max(x.init))) %>% 
  #add optimum
  bind_rows(obsmoda %>% 
              filter(y.init != 0, 
                     x.init %in% c(min(x.init), max(x.init))) %>% 
              right_join(optims %>% 
                           ungroup %>% 
                           filter(x.init %in% c(min(x.init), max(x.init))))) %>%
  #add label information
  mutate(al = ifelse(a == aest, 
                     "a == a[est]", 
                     ifelse(a == 0.25, paste("a==", a), "a == a[opt]")),
         al = fct_reorder(al, a),
         algae_conc_d = ifelse(algae_conc2 == 1, "Ambient", "Low"),
         algae_conc_d  = fct_reorder(algae_conc_d, algae_conc2)) %>% 
  #plot
  ggplot(aes(x = t, y = val, linetype = var, group = interaction(x.init), shape = midge, col = algae_conc_d))+
  facet_wrap(~al, labeller = label_parsed, ncol = 1)+
  geom_line(aes(y = wtc/4, linetype = "Midge"), size = 1)+ #midge weight
  geom_line(aes(y = gppc, linetype = "GPP"), size = 1)+ #algal production
  scale_linetype_manual(values = c("solid", "dashed"))+
  scale_y_continuous(sec.axis = sec_axis(trans = ~.*4, name = "Average Midge Mass (\u03BCg C)"))+ #add second axis for midge C
  saturation(scale_color_viridis_d(end = 0.9), scalefac(5))+
  guides(linetype = guide_legend(order = 1), color =  guide_legend(title.position = "top", order = 2))+
  labs(linetype = element_blank(),
       color = "Initial Algal Abundance",
       x = "Day",
       y = expression("Primary Production \u03BCg C "~cm^{-2}~d^{-1}))+
  theme(strip.placement = "outside",
        legend.key.width = unit(2, "lines"),
        legend.key.height = unit(0.1, "lines"),
        legend.spacing.y = unit(0.5, "lines"),
        legend.box = "vertical")

ggpreview(plot = atraj, width = 945, height = 1575, units = "px", dpi = 300)
ggsave(plot = atraj, filename = "figures/Botsch_MidgeGrowth_fig4.pdf", device = "pdf", width = 3, height = 5, units = "in", dpi = 800)


atraj + guides(linetype = guide_legend(override.aes = list(color = "white"))) +theme_black() + theme(legend.box = "vertical")
ggpreview(plot = last_plot(), width = 945, height = 1575, units = "px", dpi = 300)

#examples
obsmoda %>% 
  #extract only the estimated and high attack rates, with midges and the highest Xt=0
  filter(y.init>0,
         a %in% c(aest, 0.25),
         x.init %in% c(min(x.init), max(x.init))) %>% 
  #add optimum
  bind_rows(obsmoda %>% 
              filter(y.init != 0, 
                     x.init %in% c(min(x.init), max(x.init))) %>% 
              right_join(optims %>% 
                           ungroup %>% 
                           filter(x.init %in% c(min(x.init), max(x.init))))) %>%
  #add label information
  mutate(al = ifelse(a == aest, 
                     "Below", 
                     ifelse(a == 0.25, "Above", "Optimum")),
         al = fct_reorder(al, a),
         algae_conc_d = ifelse(algae_conc2 == 1, "High Initial Algal Abundance", "Low Initial Algal Abundance"),
         algae_conc_d  = fct_reorder(algae_conc_d, algae_conc2)) %>% 
  #plot
  ggplot(aes(x = t, y = val, linetype = var, group = interaction(x.init), shape = midge, col = algae_conc_d))+
  facet_grid(al~algae_conc_d, switch = "y")+
  geom_line(aes(y = wtc/4, linetype = "Midge", color = "Midge"), size = 1)+ #midge weight
  geom_line(aes(y = gppc, linetype = "GPP", color = "GPP"), size = 1)+ #algal production
  scale_linetype_manual(values = c("dashed", "solid"))+
  scale_y_continuous(sec.axis = sec_axis(trans = ~.*4, name = "Average Midge Mass (\u03BCg C)"))+ #add second axis for midge C
  scale_color_manual(values = c("blue", "red"))+
  # guides(linetype = guide_legend(order = 1), color =  guide_legend(title.position = "top", order = 2))+
  labs(linetype = element_blank(),
       color = element_blank(),
       x = "Day",
       y = expression("Primary Production \u03BCg C "~cm^{-2}~d^{-1}))+
  theme(strip.placement = "outside",
        legend.key.width = unit(2, "lines"),
        legend.key.height = unit(0.1, "lines"),
        legend.spacing.y = unit(0.5, "lines"),
        legend.box = "vertical",
        legend.position = "none",
        axis.title.y.right = element_text(color  = "red"),
        axis.text.y.right = element_text(color = "red"),
        axis.title.y.left = element_text(color  = "blue"),
        axis.text.y.left = element_text(color = "blue"))
  
obsmoda %>% 
  #extract only the estimated and high attack rates, with midges and the highest Xt=0
  filter(y.init>0,
         a %in% c(aest, 0.25),
         x.init %in% c(max(x.init))) %>% 
  #add optimum
  bind_rows(obsmoda %>% 
              filter(y.init != 0, 
                     x.init %in% c(max(x.init))) %>% 
              right_join(optims %>% 
                           ungroup %>% 
                           filter(x.init %in% c(min(x.init), max(x.init))))) %>%
  #add label information
  mutate(al = ifelse(a == aest, 
                     "Below", 
                     ifelse(a == 0.25, "Above", "Optimum")),
         al = fct_reorder(al, a),
         algae_conc_d = ifelse(algae_conc2 == 1, "High Initial Algal Abundance", "Low Initial Algal Abundance"),
         algae_conc_d  = fct_reorder(algae_conc_d, algae_conc2)) %>% 
  filter(!is.na(algae_conc_d)) %>% 
  #plot
  ggplot(aes(x = t, y = val, linetype = var, group = interaction(x.init), shape = midge, col = algae_conc_d))+
  facet_grid(al~algae_conc_d, switch = "y")+
  geom_line(aes(y = wtc/4, linetype = "Midge", color = "Midge"), size = 1)+ #midge weight
  geom_line(aes(y = gppc, linetype = "GPP", color = "GPP"), size = 1)+ #algal production
  scale_linetype_manual(values = c("dashed", "solid"))+
  scale_y_continuous(sec.axis = sec_axis(trans = ~.*4, name = "Average Midge Mass (\u03BCg C)"))+ #add second axis for midge C
  scale_color_manual(values = c("blue", "red"))+
  # guides(linetype = guide_legend(order = 1), color =  guide_legend(title.position = "top", order = 2))+
  labs(linetype = element_blank(),
       color = element_blank(),
       x = "Day",
       y = expression("Primary Production \u03BCg C "~cm^{-2}~d^{-1}))+
  theme(strip.placement = "outside",
        legend.key.width = unit(2, "lines"),
        legend.key.height = unit(0.1, "lines"),
        legend.spacing.y = unit(0.5, "lines"),
        legend.box = "vertical",
        legend.position = "none",
        axis.title.y.right = element_text(color  = "red"),
        axis.text.y.right = element_text(color = "red"),
        axis.title.y.left = element_text(color  = "blue"),
        axis.text.y.left = element_text(color = "blue"))

ggpreview(plot = last_plot(), width = 2.75, height = 4, units = "in", dpi = 600)
slowgrowth<- data.frame(t = unique(obsmoda$t)) %>% 
  mutate(y = max(obsmoda$y.init) + 0.005*t)

obsmoda %>% 
  filter(x.init == max(x.init),
         a == max(a),
         midge == "Midges") %>% 
  ggplot(aes(x = t))+
  geom_line(aes(y = y, color = "Resource"), size = 1)+ #midge weight
  geom_line(aes(y = y, color = "Consumer"), size = 1, data = slowgrowth)+ #algal production
  scale_linetype_manual(values = c("solid", "dashed"))+
  scale_color_manual(values = c("red", "blue"))+
  guides(linetype = guide_legend(order = 1), color =  guide_legend(title.position = "top", order = 2))+
  labs(linetype = element_blank(),
       color = element_blank(),
       x = "Time",
       y = "Biomass")+
  # theme_black()+
  theme(strip.placement = "outside",
        legend.key.width = unit(2, "lines"),
        legend.key.height = unit(0.1, "lines"),
        legend.spacing.y = unit(0.5, "lines"),
        legend.box = "vertical",
        axis.text.x = element_blank(),
        axis.text.y = element_blank())


presa <- data.frame( a = c(seq(0, 0.4, by = 0.001), a),
            K = K,
            r = r,
            c = c,
            Tmax = 22) %>% 
  crossing(x.init = unique(data$x0)) %>% 
  crossing(y.init = unique(data$y0)) %>% 
  mutate(paramset = 1:n()) %>% 
  filter(x.init == max(x.init),
         y.init == max(y.init))

#simulate
pres_a <- trajectory(presa) %>% 
  left_join(presa) %>% 
  left_join(init.data %>% 
              select(-coreid) %>% 
              unique() %>% 
              rename(x.init = x0, y.init=y0) %>% 
              mutate(midge = ifelse(y.init>0, "Midges", "No Midges"))) %>% 
  convert_to_exp_units()

#====Simulate Varying Resource Growth Rates====
#set ranges 
rrange <- c(0.5, 1.5, 3)

#parameters to simulate
simr <- data.frame( a = aest,
                    K = K,
                    rrange = rrange,
                    r = rrange*rest,
                    c = c,
                    Tmax = 22) %>% 
  crossing(x.init = unique(data$x0)) %>% 
  crossing(y.init = max(data$y0)) %>% 
  mutate(paramset = 1:n())

#simulate
rmod <- trajectory(simr) %>% 
  left_join(simr) %>% 
  left_join(init.data %>% rename(x.init = x0, y.init = y0)) %>% 
  convert_to_exp_units()

#====Figure 5b midge growth and pp under different per capita resource growth rates====
plotr <- rmod %>% 
  filter(t%in% c(Tmax)) %>% #get t = 22
  ggplot(aes(y = gppc, x = gdc, fill = algae_conc2))+
  facet_wrap(~paste("r = ", rrange, "\u00D7 r"), scales = "free", ncol = 1)+
  geom_path(aes(group  = t))+
  labs(y = expression("Primary Production \u03BCg C"~cm^{-2}~d^{-1}),
       x = expression("Midge Growth \u03BCg C "~ind^{-1}~d^{-1}),
       fill = "Initial Algal Abundance",
       shape = element_blank())+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, barheight = 0.5, barwidth = 8))+
  geom_point(shape = 22, size = 2)+
  algae_fill+
  scale_x_continuous(n.breaks= 4)+
  scale_y_continuous(n.breaks = 5)


rmod %>% 
  filter(t%in% c(Tmax)) %>% #get t = 22
  ggplot(aes(y = gppc, x = gdc, fill = algae_conc2))+
  facet_wrap(~paste("r = ", rrange, "\u00D7 r"), scales = "free", ncol = 1)+
  geom_path(aes(group  = t), color = "white")+
  labs(y = expression("Primary Production \u03BCg C"~cm^{-2}~d^{-1}),
       x = expression("Midge Growth \u03BCg C "~ind^{-1}~d^{-1}),
       fill = "Initial Algal Abundance",
       shape = element_blank())+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, barheight = 0.5, barwidth = 8))+
  geom_point(shape = 22, size = 2, color = "white")+
  algae_fill+
  scale_x_continuous(n.breaks= 4)+
  scale_y_continuous(n.breaks = 5)+
  theme_black()

#====Combine panels for Figure 4====
#4A
zplot <- plot_grid( plota +theme(axis.title = element_blank(), legend.position = "none", plot.margin = margin(5,7,5,5)),
                    plotr + theme(axis.title = element_blank(), legend.position = "none", plot.margin = margin(5,7,5,5)), 
                    nrow = 1,
                    align = "h")

leg <- get_legend(plota) # extract legend
plotfin <- plot_grid(leg, zplot, ncol = 1, rel_heights = c(0.15, 0.85)) #plot with legend

#add axis labels
y.grob <- textGrob(expression("Primary Production \u03BCg C"~cm^{-2}~d^{-1}), rot=90, gp = gpar(fontsize = 11))
x.grob <- textGrob(expression("Midge Growth \u03BCg C "~ind^{-1}~d^{-1}), gp = gpar(fontsize = 11))
#finished 4A
comb <- grid.arrange(plotfin, left = y.grob, bottom = x.grob)


#preview and save
ggpreview(plot = comb, width = 945, height = 1260, units = "px", dpi = 300)
ggsave(plot = comb, filename = "figures/Botsch_MidgeGrowth_fig5.pdf",  device = "pdf", width = 3, height = 4, units = "in", dpi = 800)
