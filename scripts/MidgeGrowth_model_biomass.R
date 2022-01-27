#=====load packages====
library(tidyverse)
library(cowplot)
options(dplyr.summarise.inform = FALSE)

theme_set(theme_bw()+
            theme(panel.grid = element_blank(),
                  strip.background = element_blank(),
                  legend.position = "bottom"))



meanna <- function(x){
  if(all(is.na(x))){
    NA
  }
  else{
    mean(x, na.rm = TRUE)
  }
}


#Lindegaard et al (1979) formula for converting length to biomass (tab.21)
#AFDW in mg^(1/3) = 0.0442 + 0.0879 * length in mm 
weight <- function(l){
  (0.0442 + 0.0879 * l )^3 
}

ggpreview <- function(...) {
  fname <- tempfile(fileext = ".png")
  ggsave(filename = fname, ...)
  system2("open", fname)
  invisible(NULL)
}

#====Model=====
theta <- 3/4 #allometric scaling factor (currently using 3/4 on the basis of the WBE model, Riveros and Enquist 2011, )

trajectory <- function(simdat){ # a dataframe containing columns x.init, y.init,a, K, c, m, Bm, Tmax 
  #list
  trajlist <- list()
  
  for(i in 1:nrow(simdat)){
    x <- simdat$x.init[i]
    y <- simdat$y.init[i]
    r <- simdat$r[i]
    a <- simdat$a[i]
    K <- simdat$K[i]
    c <- simdat$c[i]
    m <- simdat$m[i]
    Tmax <- simdat$Tmax[i]

    out <- data.frame(paramset = i, t = 1:Tmax,x = NA, y = NA)
    out$x[1] = x
    out$y[1] = y
    for(t in 2:Tmax){
      
      out$x[t] <- max(0,exp(r)*out$x[t-1]^K - a*out$x[t-1]*out$y[t-1]) #growth of algae ##CHECK 
      out$y[t] <- max(0, out$y[t-1] + c*a*out$x[t-1]*out$y[t-1]- m*out$y[t-1]) #midge growth
      
    }	
    
    trajlist[[i]] <- out
  }
  return(bind_rows(trajlist))
}


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
    summarise(y = mean(Bt))}$y

# prune to experimental falcon tubes
data <- data %>% 
  filter(day>0) %>% 
  arrange(day, coreid)





#############################################################
#=====fit with biomass specific loss rates====



SSfit2 <- function(par, data, x.init, y.init, SS.weight, par.fixed=par.fixed, tofit = T){
  
# standardize variables
data <- data %>% 
  mutate(x = gpp/mean(gpp, na.rm = TRUE),
         x0 = chl/mean(chl, na.rm = TRUE),
         y = Bt/mean(Bt, na.rm = TRUE),
         y0 = ifelse(midge == "Midges", y.init.mean/mean(Bt, na.rm = TRUE), 0)) %>% 
  arrange(coreid)

  par.temp <- par.fixed
  par.temp[is.na(par.fixed)] <- par
  par <- par.temp
  
  r <- exp(par[1])
  K <-  exp(par[2])
  a <-  exp(par[3])
  ac <-  exp(par[4])
  Bm <-  exp(par[5])
  gpp.rate <-  exp(par[6])
  
  x <- x.init
  y <- y.init
  
  Tmax <- 22
  
  SS1 <- 0
  SS2 <- 0
  for(t in 1:Tmax){
    
    x.next <- exp(r)*x^K - a*x*y
    y <- y + ac*x*y - Bm*y
    x <- x.next
    
    if(t == 14){
      dif.x <- log(gpp.rate*x) - log(data$x[data$day == 14])
      SS1 <- SS1 + sum(dif.x^2)
      
      sampled.data <- (data$day == 14) & !is.na(data$y) & (data$midge == "Midges")
      sampled.sim <- data$coreid[sampled.data]
      dif.y <- log(y[sampled.sim]) - log(data$y[sampled.data])
      #dif.y <- y[sampled.sim] - data$wt[sampled.data]
      SS2 <- SS2 + sum(dif.y^2)
    }
    
    if(t == 22){
      sampled.data <- data$day == 22
      sampled.sim <- data$coreid[sampled.data]
      dif.x <- log(gpp.rate*x[sampled.sim]) - log(data$x[sampled.data])
      SS1 <- SS1 + sum(dif.x^2)
      
      sampled.data <- (data$day == 22) & (data$midge == "Midges") & !is.na(data$y)
      sampled.sim <- data$coreid[sampled.data]
      dif.y <- log(y[sampled.sim]) - log(data$y[sampled.data])
      #dif.y <- y[sampled.sim] - data$wt[sampled.data]
      SS2 <- SS2 + sum(dif.y^2)
    }
  }
  if(tofit){
    if(any(is.nan(x))) {
      return(10^8)
    }else{
      return(SS1 + SS.weight*SS2)
    }
  }else{
    return(c(SS1,SS.weight*SS2, exp(par), mean(x), mean(y)))
  }
}

# start values
r <- .5
K <- .9
a <- .5
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

par.full2 <- log(c(r, K, a, ac, m, gpp.rate))
par.fixed2 <- c(NA, NA, NA, NA, NA, NA)
par2 <- par.full2[is.na(par.fixed)]


#====Fit2====
SSmin <- 10^10
nrep <- 10
for (i.rep in 1:nrep){
  opt <- optim(par = par, fn=SSfit2, data=data, x.init = x.init, y.init = y.init, SS.weight = SS.weight, par.fixed=par.fixed, method = "SANN")
  opt <- optim(par = opt$par, fn=SSfit2, data=data, x.init = x.init, y.init = y.init, SS.weight = SS.weight, par.fixed=par.fixed, method = "Nelder-Mead")
  if(opt$value < SSmin){
    SSmin <- opt$value
    opt.opt <- opt
  }
  par2 <- exp(rnorm(n=length(par), sd=1)) * opt.opt$par
  show(c(opt$value, SSfit(par=opt$par, data=data, x.init = x.init, y.init=y.init, SS.weight = SS.weight, par.fixed=par.fixed, tofit=F)))
}
par.temp2 <- par.fixed2
par.temp2[is.na(par.fixed)] <- opt.opt$par
par2 <- par.temp2

r2 <- exp(par2[1])
K2 <- exp(par2[2])
a2 <- exp(par2[3])
ac2 <- exp(par2[4])
c2 <- ac2/a2
Bm <- exp(par2[5])
gpp.rate <- exp(par2[6])
round(c(r=r2, K=K2, a=a2, c=c2, m=Bm, gpp.rate=gpp.rate, SS=opt.opt$value), digits=8)

# r        K        a        c        m gpp.rate       SS 
# 0.216    0.919    0.034    0.671    0.071    0.161  274.263 




obssim <- data.frame(r = r2,
                     K = K2,
                     a = a2,
                     c = c,
                     m = Bm,
                     Tmax = 22,
                     y.init = unique(data$y0)) %>% 
  crossing(x.init = unique(data$x0)) %>% 
  mutate(paramset = 1:n())

chck <- trajectory(obssim) %>% 
  left_join(obssim)

chck %>% 
  mutate(x = x*gpp.rate,
         midge = ifelse(y.init>0, "Midges", "No Midges")) %>% 
  gather(var, val, x, y) %>% 
  ggplot(aes(x = t, y = val, color = var, linetype = midge))+
  facet_grid(midge~round(x.init, 3))+
  geom_line()+
  geom_point( alpha = 0.7, position=position_dodge(width = 2), data = data %>% rename(t = day) %>%  gather(var, val, x, y))+
  geom_point( alpha = 0.7, position= position_dodge(width = 2), data = data %>% select(algae_conc2, midge, x0, y0) %>% mutate(x.init = x0) %>%  rename(x = x0, y = y0) %>%  gather(var, val, x, y) %>% mutate(t = 0))+
  
  scale_shape_manual(values = c(16,21))
  


simbmdat <- data.frame(r = r2,
                       Ny = 1,
                       K = K2,
                       a = seq(0, 0.2, by = 0.01),
                       c = c2,
                       m = Bm,
                       y.init = max(data$y0),
                       Tmax = 22) %>%
  crossing(x.init = unique(data$x0)) %>%
  mutate(paramset = 1:n())

simbm <- trajectory(simbmdat) %>%
  left_join(simbmdat)


simbm %>% 
  filter(t %in% c(22),
         y.init>0) %>%
  ggplot(aes(x = a, y = y-y.init, group = x.init, fill = x.init))+
  geom_vline(xintercept = a2)+
  scale_x_continuous(breaks = c(0, 0.1, 0.2))+
  geom_line()+
  geom_point(shape = 21, alpha = 0.7)+
  labs(x = "Attack Rate",
       y = "Standardized Secondary Production", 
       fill = "Initial Resource Biomass")+
  scale_fill_viridis_c()+
  guides(fill = guide_colorbar(title.position = "top"))

# ggpreview(plot = last_plot(), dpi = 650, width = 3, height = 5, units = "in")

simbm %>% 
  filter(t %in% c(22),
         y.init>0) %>%
  mutate(t2 = ifelse(t == 5000, "Equilibrium", paste("t =", t)),
         t3 = fct_reorder(t2, t, .desc = FALSE) ) %>% 
  ggplot(aes(x = a, y = y-y.init, group = x.init, fill = x.init))+
  geom_vline(xintercept = a2)+
  # facet_wrap(~t3, scales = "free_y", ncol = 1)+
  scale_x_continuous(breaks = c(0, 0.1, 0.2))+
  geom_line()+
  geom_point(shape = 21, alpha = 0.7)+
  labs(x = "Attack Rate",
       y = "Standardized Secondary Production", 
       fill = "Initial Resource Biomass")+
  scale_fill_viridis_c()+
  guides(fill = guide_colorbar(title.position = "top"))

# ggpreview(plot = last_plot(), width = 3, height = 4, units = "in", dpi = 650)



simbm %>% 
  filter(t %in% c(22),
         y.init>0, 
         a %in% c(0.02, 0.03, 0.08, 0.15)) %>%
  mutate(t2 = ifelse(t == 5000, "Equilibrium", paste("t =", t)),
         t3 = fct_reorder(t2, t, .desc = FALSE)) %>% 
  ggplot(aes(x = exp(r)*x^K, y = y-y.init, fill = x.init))+
  # geom_vline(xintercept = a2)+
  facet_wrap(~paste("\u03b1 =", a), scales = "free", ncol = 2)+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
  
  geom_path()+
  geom_point(shape = 21, alpha = 0.7)+
  labs(x = "Primary Production",
       y = "Secondary Production", 
       fill = "Initial Resource Biomass")+
  scale_fill_viridis_c()+
  guides(fill = guide_colorbar(title.position = "top"))

# ggpreview(plot = last_plot(), width = 3, height = 4, units = "in", dpi = 650)



simbm %>% 
  filter(t %in% c(65),
         y.init>0, 
         a %in% c(0.01, 0.05, 0.08, 0.15)) %>%
  mutate(t2 = ifelse(t == 5000, "Equilibrium", paste("t =", t)),
         t3 = fct_reorder(t2, t, .desc = FALSE)) %>% 
  ggplot(aes(x = exp(r)*x^K, y = y-y.init, fill = x.init))+
  # geom_vline(xintercept = a2)+
  facet_wrap(~paste("\u03b1 =", a), scales = "free", ncol = 2)+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3))+
  geom_path()+
  geom_point(shape = 21, alpha = 0.7)+
  labs(x = "Primary Production",
       y = "Secondary Production", 
       fill = "Initial Resource Biomass")+
  scale_fill_viridis_c()+
  guides(fill = guide_colorbar(title.position = "top"))

# ggpreview(plot = last_plot(), width = 3, height = 4, units = "in", dpi = 650)
