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

# ggpreview(plot = arange, width = 3, height = 3.5, units = "in", dpi = 800)
# ggsave(plot = arange, filename = "figures/Botsch_MidgeGrowth_fig3.pdf", device = "pdf", width = 3, height = 3.5, units = "in", dpi = 800)

#====Figure 4Aa midge growth and pp under different attack rates====
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
  scale_y_continuous(n.breaks = 5)

#====Figure 4Ba biomass trajectories under different attack rates====
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
         al = fct_reorder(al, a)) %>% 
  #plot
  ggplot(aes(x = t, y = val, linetype = var, group = interaction(x.init), shape = midge, col = x.init))+
  facet_wrap(~al, labeller = label_parsed, ncol = 1)+
  geom_line(aes(y = wtc/4, linetype = "Midge"), size = 1)+ #midge weight
  geom_line(aes(y = gppc, linetype = "GPP"), size = 1)+ #algal production
  scale_linetype_manual(values = c("solid", "dashed"))+
  scale_y_continuous(sec.axis = sec_axis(trans = ~.*4, name = "Average Midge Mass (\u03BCg C)"))+ #add second axis for midge C
  algae_color+
  guides(color = "none")+
  labs(linetype = element_blank(),
       x = "Day",
       y = expression("Primary Production \u03BCg C "~cm^{-2}~d^{-1}))+
  theme(strip.placement = "outside",
        legend.key.width = unit(2, "lines"))

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

#====Figure 4Ab midge growth and pp under different per capita resource growth rates====
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


#====Figure 4Ba biomass trajectories under different per capita resource growth rates====
rtraj <- rmod %>% 
  #only min and max Xt=0
  filter(x.init %in% c(min(x.init), max(x.init))) %>% 
  ggplot(aes(x = t, group = interaction(x.init), color = algae_conc2))+
  facet_wrap(~paste("r = ", rrange, "\u00D7 r"), scales = "free_y", ncol = 1)+
    geom_line(aes(y = wtc/4, linetype = "Midge"), size = 1)+
  geom_line(aes(y = gppc, linetype = "GPP"), size = 1)+
  scale_linetype_manual(values = c("solid", "dashed"))+
  scale_alpha_manual(values = c(0.25, 0.75), guide = "none")+
  scale_shape_manual(values = c(16,21))+
  algae_color+
  guides(color = "none")+
  scale_y_continuous(sec.axis = sec_axis(trans = ~.*4, name = "Average Midge Mass (\u03BCg C)"))+
  labs(linetype = element_blank(),
       x = "Day",
       y = expression("GPP \u03BCg C"~cm^{-2}~d^{-1}))+
  theme(strip.placement = "outside",
        legend.key.width = unit(2, "lines"))

#====Combine panels for Figure 4====
#4A
zplot <- plot_grid( plota +theme(axis.title = element_blank(), legend.position = "none"),
                    plotr + theme(axis.title = element_blank(), legend.position = "none"), 
                    nrow = 1,
                    align = "h")

leg <- get_legend(plota) # extract legend
plotfin <- plot_grid(leg, zplot, ncol = 1, rel_heights = c(0.15, 0.85)) #plot with legend

#add axis labels
y.grob <- textGrob(expression("Primary Production \u03BCg C"~cm^{-2}~d^{-1}), rot=90, gp = gpar(fontsize = 11))
x.grob <- textGrob(expression("Midge Growth \u03BCg C "~ind^{-1}~d^{-1}), gp = gpar(fontsize = 11))
#finished 4A
comb <- grid.arrange(plotfin, left = y.grob, bottom = x.grob)

#4B
tplot <- plot_grid(atraj + theme(axis.title.x = element_blank(), axis.title.y.right = element_blank(), legend.position = "none"),
          rtraj + theme(axis.title.x = element_blank(), axis.title.y.left = element_blank(), legend.position = "none"), align = "h")
leg2 <- get_legend(atraj) #extract legend for 4B

#add axis labels
x.grob2 <- textGrob("Day", gp = gpar(fontsize = 11))
plott <- plot_grid(leg2, tplot, ncol = 1, rel_heights = c(0.15, 0.85))
#finished 4B
combt <- grid.arrange(plott, bottom = x.grob2)

#combine 4A and 4B
combf <- plot_grid(comb, combt, rel_widths = c(1, 1.2), labels = c("A", "B"), align = "h")

#preview and save
# ggpreview(plot = combf, width = 6, height = 4, units = "in", dpi = 800)
# ggsave(plot = combf, filename = "figures/Botsch_MidgeGrowth_fig4.pdf",  device = "pdf",
#        width = 6.5, height = 4, units = "in", dpi = 800)
