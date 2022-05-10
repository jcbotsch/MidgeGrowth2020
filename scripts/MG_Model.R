#=====load packages====
library(tidyverse)
library(cowplot)
source("scripts/MG_Functions.R")
library(gridExtra)
library(grid)

options(dplyr.summarise.inform = FALSE)



#############################################################
#====Read Data===
meta <- read_csv("clean data/MG_meta.csv") %>% 
  mutate(algae_conc2 = ifelse(algae_conc == 0, 2e-3, algae_conc)) #half of lowest value
chl_raw <- read_csv("clean data/MG_chl.csv")
ep_raw <- read_csv("clean data/MG_ep.csv")
cm_raw <- read_csv("clean data/MG_cm.csv")
cc_raw <- read_csv("clean data/MG_cc.csv")
ndvi_raw <- read_csv("clean data/MG_ndvi.csv")

#====Prepare Data====
chl <- chl_raw %>% 
  left_join(meta) %>% 
  group_by(algae_conc2) %>% 
  summarise(chl = mean(chl, na.rm = TRUE)/1000) 


gpp <- ep_raw %>% 
  select(coreid, day, gpp)

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
  mutate(x = gpp/mean(gpp, na.rm = TRUE),
         x0 = chl/mean(chl, na.rm = TRUE),
         y = wt/mean(wt, na.rm = TRUE),
         y0 = ifelse(midge == "Midges", y.init.mean/mean(wt, na.rm = TRUE), 0)) %>% 
  arrange(coreid)


#====estimate parameters====

# start values
r <- .5
K <- .9
a <- .2
ac <- .03
m <- .001
gpp.rate <- 10^0

#get initial data
init.data <- data %>% 
  select(coreid, algae_conc2, x0, y0) %>% 
  unique()

x.init <- init.data$x0
y.init <- init.data$y0

SS.weight <- 30

par.full <- log(c(r, K, a, ac, m, gpp.rate))
par.fixed <- c(NA, NA, NA, NA, NA, NA)
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

r <- exp(par[1])
K <- exp(par[2])
a <- exp(par[3])
ac <- exp(par[4])
c <- ac/a
m <- exp(par[5])
gpp.rate <- exp(par[6])
round(c(r=r, K=K, a=a, c=c, m=m, gpp.rate=gpp.rate, SS=opt.opt$value), digits=4)
# r        K        a        c        m       gpp.rate       SS 
# 0.1725   0.9290   0.0447   0.1788   0.0015 0.2088       66.8795 


#====Plot Model and Data used to estimate parameters====
#simulate model as is 
sim <- data.frame(x.init = unique(data$x0),
                  r = exp(par[1]),
                  K = exp(par[2]),
                  a = exp(par[3]),
                  c = exp(par[4])/exp(par[3]), 
                  m = exp(par[5]),
                  Ny = 1,
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
  ggplot(aes(x = t, y = val, col = var, linetype = midge, shape = midge))+
  facet_wrap(~round(algae_conc2, 3), ncol = 3, scales = "free_y")+
  geom_line(aes(y = y, col = "Midge"))+
  geom_line(aes(y = x*gpp.rate, col = "Algae"))+
  geom_point(aes(y = val, color = taxon, x = day), alpha = 0.5, position = position_dodge(width = 2),  data = data %>% rename(Midge = y, Algae = x) %>%  gather(taxon, val, Midge, Algae))+
  geom_point(aes(y = val, color = taxon), alpha = 0.5, position = position_dodge(width = 2),  data = init.data %>% rename(Midge = y0, Algae = x0) %>% mutate(t = 1, midge = ifelse(Midge>0, "Midges", "No Midges")) %>%   gather(taxon, val, Midge, Algae))+
  # scale_y_continuous("Algae", sec.axis = sec_axis(~./5, "Midges"))+
  scale_linetype_manual(values = c("solid", "dotted"))+
  scale_shape_manual(values = c(16,21))+
  scale_color_manual(values = c("black", "red"))+
  labs(linetype = element_blank(),
       color = element_blank(),
       shape = element_blank(),
       x = "Time",
       y = "Scaled Biomass")
# ggpreview(plot = mod_data, dpi = 650, width = 5, height = 6, units = "in")

#====Vary attack rates====
#parameters to simulate
sima <- data.frame( a = c(seq(0, 0.4, by = 0.001), a),
                    K = K,
                    r = r,
                    c = c, 
                    m = m,
                    Bm = 0,
                    Ny = 1,
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

#plot 
arange <- obsmoda %>% 
  filter(t == Tmax,
         y.init>0) %>% 
  ggplot(aes(x = a, y = y-y.init, fill = x.init, group = x.init))+
  geom_vline(xintercept = aest)+
  geom_line(size = 1.1)+
  geom_line(aes(color = x.init))+
  labs(x = "\u03b1",
       y = "Standardized Consumer Growth",
       color = "Initial Resource Biomass")+
  guides(color = guide_colorbar(title.position = "top", title.hjust = 0.5))+
  scale_color_viridis_c(trans = "log", breaks = c(0.25, 1,4))

# ggpreview(plot = arange, width = 3, height = 4, units = "in", dpi = 650)

#extract optima
optima <- obsmoda %>% 
  filter(t == Tmax,
         y.init>0) %>% 
  filter(y-y.init==max(y-y.init)) %>% 
  pull(a)


optims <- obsmoda %>% 
  filter(t == Tmax,
         y.init>0) %>% 
  group_by(x.init) %>% 
  filter(y-y.init==max(y-y.init)) %>% 
  select(x.init, a)

#plot variation in optima 
optims %>% 
  bind_rows(crossing(a = optima, x.init = optims$x.init)) %>% 
  left_join(obsmoda %>% filter(y.init>0)) %>% 
  rename(resource = x, consumer = y) %>% 
  gather(var, val, resource, consumer) %>% 
  mutate(al = ifelse(x.init < 4 & a != optima, "optimum", "max optimum")) %>% 
  ggplot(aes(x = t, y = val, color = var))+
  facet_wrap(~round(x.init, digits = 3), scales = "free_y")+
  geom_line(aes( linetype = al))+
  geom_text(aes(y = 0.1, x = 15, label = a), color = "gray40", data = optims)+
  lims(y = c(0, NA))+
  labs(y = "Biomass",
       color = element_blank(),
       linetype = element_blank())

# ggpreview(plot = last_plot(), width  = 6, height = 6.5, units = "in", dpi = 650)


plota <- obsmoda %>% 
  filter(t %in% c(Tmax),
         y.init>0,
         a %in% c(aest, 0.2)) %>% 
  bind_rows(optims %>% 
              left_join(obsmoda %>% 
                          filter(y.init>0,
                                 t %in% c(22)))) %>% 
  mutate(al = ifelse(a == aest, 
                     "\u03b1 = \u03b1 est", 
      ifelse(a == 0.2, paste("\u03b1 = ", a), "\u03b1 = \u03b1 opt")),
         al = fct_reorder(al, a, min)) %>% 
  ggplot(aes(x = exp(r)*x^K, y = y-y.init, fill = x.init))+
  facet_wrap(~al, scales = "free", labeller = label_value, ncol = 1)+
  geom_path(aes(group = t))+
  labs(x = "Standardized Primary Production",
       y = "Standardized Consumer Growth",
       fill = "Initial Resource Biomass",
       shape = element_blank())+
  geom_point(shape = 22, size = 2)+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, barheight = 0.5, barwidth = 8))+
  scale_fill_viridis_c(trans = "log", breaks = c(0.25, 1,4))+
  scale_x_continuous(n.breaks= 3)+
  scale_y_continuous(n.breaks = 3)

plotas <- obsmoda %>% 
  filter(t %in% c(14),
         y.init>0,
         a %in% c(aest, 0.2)) %>% 
  bind_rows(optims %>% 
              left_join(obsmoda %>% 
                          filter(y.init>0,
                                 t %in% c(14)))) %>% 
  mutate(al = ifelse(a == aest, 
                     "\u03b1 = \u03b1 est", 
                     ifelse(a == 0.2, paste("\u03b1 = ", a), "\u03b1 = \u03b1 opt")),
         al = fct_reorder(al, a, min)) %>% 
  ggplot(aes(x = exp(r)*x^K, y = y-y.init, fill = x.init))+
  facet_wrap(~al, scales = "free", labeller = label_value, ncol = 1)+
  geom_path(aes(group = t))+
  labs(x = "Standardized Primary Production",
       y = "Standardized Consumer Growth",
       fill = "Initial Resource Biomass",
       shape = element_blank())+
  geom_point(shape = 21, size = 2)+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, barheight = 0.5, barwidth = 8))+
  scale_fill_viridis_c(trans = "log", breaks = c(0.25, 1,4))+
  scale_x_continuous(n.breaks= 3)+
  scale_y_continuous(n.breaks = 3)


# ggpreview(plot = last_plot(), width = 3, height = 4, units = "in", dpi = 650)


#trajectories
obsmoda %>% 
  filter(y.init>0,
         a %in% c(aest, optima, 0.2)) %>% 
  rename(consumer = y, resource = x) %>% 
  filter(x.init %in% unique(.$x.init)[c(1,6,10)]) %>%
  gather(var, val, consumer, resource) %>%
  mutate(facetlab = paste("\u03b1 = ", a),
         facetlab = fct_reorder(facetlab, a)) %>% 
  mutate(val = ifelse(var == "consumer", val, val)) %>% 
  ggplot(aes(x = t, y = val))+
  facet_grid(facetlab~paste( "x0 = ", round(x.init,1)), scales = "free_y")+
  geom_line(aes(col = var, group = interaction(var, x.init)))+
  # geom_point(aes(fill = x.init), shape = 21)+
  labs(x = "t",
       y = "Biomass",
       fill = "Initial Resource Biomass",
       color = "")+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))+
  scale_fill_viridis_c(trans = "log", breaks = c(0.25, 1,4))

# ggpreview(plot = last_plot(), width = 6, height  = 6, units = "in", dpi = 650)

#====Varying Resource Growth Rates====
rrange <- c(0.5, 1.5, 2.5)

simr <- data.frame( a = aest,
                    K = K,
                    rrange = rrange,
                    r = rrange*r,
                    c = c, 
                    m = m,
                    Bm = 0,
                    Ny = 1,
                    Tmax = 22) %>% 
  crossing(x.init = unique(data$x0)) %>% 
  crossing(y.init = max(data$y0)) %>% 
  mutate(paramset = 1:n())


rmod <- trajectory(simr) %>% 
  left_join(simr)

plot1 <- rmod %>% 
  filter(t%in% c(Tmax)) %>% 
  ggplot(aes(x = exp(r)*x^K, y = y-y.init, fill = x.init))+
  facet_wrap(~paste("r = ", rrange, "\u00D7 r"), scales = "free", ncol = 1)+
  geom_path(aes(group  = t))+
  labs(x = "Standardized Primary Production",
       y = "Standardized Consumer Growth",
       fill = "Initial Resource Biomass")+
  geom_point(shape = 22, size = 2)+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, barheight = 0.5, barwidth = 8))+
  scale_fill_viridis_c(trans = "log", breaks = c(0.25, 1,4))+
  scale_x_continuous(n.breaks= 3)+
  scale_y_continuous(n.breaks = 3)


plot1s <- rmod %>% 
  filter(t%in% c(14)) %>% 
  ggplot(aes(x = exp(r)*x^K, y = y-y.init, fill = x.init))+
  facet_wrap(~paste("r = ", rrange, "\u00D7 r"), scales = "free", ncol = 1)+
  geom_path(aes(group  = t))+
  labs(x = "Standardized Primary Production",
       y = "Standardized Consumer Growth",
       fill = "Initial Resource Biomass")+
  geom_point(shape = 21, size = 2)+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, barheight = 0.5, barwidth = 8))+
  scale_fill_viridis_c(trans = "log", breaks = c(0.25, 1,4))+
  scale_x_continuous(n.breaks= 3)+
  scale_y_continuous(n.breaks = 3)



#plot trajectories
rmod %>% 
  rename(consumer = y, resource = x) %>%
  gather(var, val, consumer, resource) %>% 
  filter(x.init %in% unique(.$x.init)[c(1,6,10)]) %>%
  ggplot(aes(x = t, y = val))+
  facet_grid(paste("change in growth rate = ",rrange)~paste("x0 =",round(x.init, 2)), scales = "free_y")+
  geom_line(aes(color = var, group = interaction(var, r)))+
  labs(x = "t",
       y = "Biomass",
       fill = "Initial Resource Availability",
       color = "")+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))+
  scale_fill_viridis_c(option = "plasma")

# ggpreview(plot = last_plot(), width = 6, height  = 6, units = "in", dpi = 650)

# simulate across a range of attack rates
ar <- data.frame( a = c(seq(0, 0.2, by = 0.001), a),
                  K = K,
                  y.init = max(data$y0),
                  c = c, 
                  m = m,
                  Bm = 0,
                  Ny = 1,
                  Tmax = 22) %>% 
  crossing(x.init = unique(data$x0)) %>% 
  crossing(rrange =c(0.5, 1.5)) %>% 
  mutate(r= r*rrange) %>% 
  mutate(paramset = 1:n())



simar <- trajectory(ar) %>% 
  left_join(ar)

simar %>% 
  filter(t == Tmax,
  ) %>% 
  ggplot(aes(x = a, y = y-y.init, fill = x.init, group = x.init))+
  facet_wrap(~paste(rrange, "\u00D7 Resource Growth Rate"), scales = "free")+
  geom_line(size = 1.1)+
  geom_line(aes(color = x.init))+
  labs(x = "\u03b1",
       y = "Standardized Consumer Growth",
       color = "Initial Resource Biomass")+
  guides(color = guide_colorbar(title.position = "top", title.hjust = 0.5))+
  scale_color_viridis_c(trans = "log", breaks = c(0.25, 1,4))


#====Vary Consumer Density====

simxy <- data.frame( a = aest,
            K = K,
            r = r,
            c = c, 
            m = m,
            Bm = 0,
            Ny = c(0.5, 5, 30),
            Tmax = 22) %>% 
  crossing(x.init = c(unique(data$x0))) %>% 
  crossing(y.init = max(data$y0)) %>% 
  mutate(paramset = 1:n())


xyres <- trajectory(simxy) %>% 
  left_join(simxy)



plot2 <- xyres %>% 
  filter(t%in% c(Tmax)) %>% 
  mutate(facetlab = paste("Ny = ", Ny),
         facetlab = fct_reorder(facetlab, Ny)) %>% 
  ggplot(aes(x = exp(r)*x^K, y = y-y.init, fill = x.init))+
  facet_wrap(~facetlab, scales = "free", ncol = 1)+
  geom_path(aes(group = t))+
  labs(x = "Standardized Primary Production",
       y = "Standardized Consumer Growth",
       fill = "Initial Resource Biomass")+
  geom_point(shape = 22, size = 2)+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, barheight = 0.5, barwidth = 8))+
  scale_fill_viridis_c(trans = "log", breaks = c(0.25, 1,4))+
  scale_x_continuous(n.breaks= 3)+
  scale_y_continuous(n.breaks = 3)

plot2s <- xyres %>% 
  filter(t%in% c(14)) %>% 
  mutate(facetlab = paste("Ny = ", Ny),
         facetlab = fct_reorder(facetlab, Ny)) %>% 
  ggplot(aes(x = exp(r)*x^K, y = y-y.init, fill = x.init))+
  facet_wrap(~facetlab, scales = "free", ncol = 1)+
  geom_path(aes(group = t))+
  labs(x = "Standardized Primary Production",
       y = "Standardized Consumer Growth",
       fill = "Initial Resource Biomass")+
  geom_point(shape = 21, size = 2)+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, barheight = 0.5, barwidth = 8))+
  scale_fill_viridis_c(trans = "log", breaks = c(0.25, 1,4))+
  scale_x_continuous(n.breaks= 3)+
  scale_y_continuous(n.breaks = 3)


#trajectories
xyres %>% 
  rename(consumer = y, resource = x) %>% 
  filter(x.init %in% unique(xyres$x.init)[c(1,6,10)]) %>%
  gather(var, val, consumer, resource) %>%
  mutate(facetlab = paste("Ny = ", Ny),
         facetlab = fct_reorder(facetlab, Ny)) %>% 
  ggplot(aes(x = t, y = val))+
  facet_grid(facetlab~paste( "x0 = ", round(x.init,1)), scales = "free_y")+
  geom_line(aes(col = var, group = interaction(var, x.init)))+
  # geom_point(aes(fill = x.init), shape = 21)+
  labs(x = "t",
       y = "Biomass",
       fill = "Initial Resource Biomass",
       color = "")+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))+
  scale_fill_viridis_c(trans = "log", breaks = c(0.25, 1,4))

# ggpreview(plot = last_plot(), width = 6, height  = 6, units = "in", dpi = 650)

#simulate across a range of attack rates
an <- data.frame( a = c(seq(0, 0.05, by = 0.001), a),
                  K = K,
                  r = r,
                  c = c, 
                  m = m,
                  Bm = 0,
                  y.init = max(y.init),
                  Tmax = 22) %>% 
  crossing(x.init = unique(data$x0)) %>% 
  crossing(Ny = c(5)) %>% 
  mutate(paramset = 1:n())

anx <- data.frame( a = c(seq(0, 0.05, by = 0.001), a),
                  K = K,
                  r = r,
                  c = c, 
                  m = m,
                  Bm = 0,
                  y.init = max(y.init),
                  Tmax = 22) %>% 
  crossing(x.init = unique(data$x0)) %>% 
  crossing(Ny = c(30)) %>% 
  mutate(paramset = 1:n())



siman <- trajectory(an) %>% 
  left_join(an) %>% 
  bind_rows(trajectory(anx) %>% 
              left_join(anx))
#plot
siman %>% 
  filter(t == Tmax,
  ) %>% 
  ggplot(aes(x = a, y = y-y.init,  group = x.init))+
  facet_wrap(~paste("Consumer Density =", Ny), scales = "free")+
  geom_line(size = 1.1)+
  geom_line(aes(color = x.init))+
  labs(x = "\u03b1",
       y = "Standardized Consumer Growth",
       color = "Initial Resource Biomass")+
  guides(color = guide_colorbar(title.position = "top", title.hjust = 0.5))+
  scale_color_viridis_c(trans = "log", breaks = c(0.25, 1,4))


#====Combine PP and SP plots====
zplot <- plot_grid( plota +theme(axis.title = element_blank(), legend.position = "none"),
                    plot2+ theme(axis.title = element_blank(), legend.position = "none"), 
                    plot1 + theme(axis.title = element_blank(), legend.position = "none"), 
                    labels = c("A", "B", "C"),
                    nrow = 1)
leg <- get_legend(plota)

plotfin <- plot_grid(leg, zplot, ncol = 1, rel_heights = c(0.15, 0.85))

y.grob <- textGrob("Standardized Consumer Growth", rot=90)

x.grob <- textGrob("Standardized Primary Production")

comb <- grid.arrange(plotfin, left = y.grob, bottom = x.grob)

# ggpreview(plot = comb, width = 4, height = 5, units = "in", dpi = 650)



zsplot <- plot_grid( plotas +theme(axis.title = element_blank(), legend.position = "none"),
                    plot2s+ theme(axis.title = element_blank(), legend.position = "none"), 
                    plot1s + theme(axis.title = element_blank(), legend.position = "none"), 
                    labels = c("A", "B", "C"),
                    nrow = 1)

zs2 <- plot_grid(leg, zsplot, ncol = 1, rel_heights = c(0.15, 0.85))


combs <- grid.arrange(zs2, left = y.grob, bottom = x.grob)

# ggpreview(plot = combs, width = 4, height = 5, units = "in", dpi = 650)

