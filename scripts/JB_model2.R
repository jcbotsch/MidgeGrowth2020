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
#model the dynamics through time
# trajectory.Ives <- function(z.init, r, K, a, c, m, Tmax){
#   
#   for(i in 1:nrow(z.init)){
#     x <- z.init[i,1]
#     y <- z.init[i,2]
#     
#     X <- matrix(NA, 2, Tmax)
#     for(t in 1:Tmax){
#       #x.next <- exp(r*(1 - x/K))*x - a*x*y
#       x.next <- exp(r)*x^K - a*x*y
#       y <- y + c*a*x*y - m
#       x <- x.next
#       y <- max(y,0)
#       
#       X[,t] <- matrix(c(x,y),2,1)
#     }	
#     if(i == 1) maxX <- max(X)
#     plot(X[1,], typ="l", ylim=c(0,maxX))
#     lines(X[2,], col="red")
#   }
# }
# 

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
    Bm <- simdat$Bm[i] #lgl (0 or 1)
    Tmax <- simdat$Tmax[i]
    Ny <- simdat$Ny[i]
    
    out <- data.frame(paramset = i, t = 1:Tmax,x = NA, y = NA)
    out$x[1] = x
    out$y[1] = y
    for(t in 2:Tmax){
      
      out$x[t] <- max(0,exp(r)*out$x[t-1]^K - Ny*(a*out$x[t-1]*out$y[t-1])) #growth of algae ##CHECK 
      out$y[t] <- max(0, out$y[t-1] + c*a*out$x[t-1]*out$y[t-1] - (1-Bm)*m - Bm*m*out$y[t-1] ) #midge growth
      
    }	
    
    trajlist[[i]] <- out
  }
  return(bind_rows(trajlist))
}





# data.frame(x.init = 1,
#            y.init = 0.5,
#            a = seq(0, 0.2, by = 0.01),
#            r = 0.2,
#            K = 0.9,
#            c = 0.2,
#            m = 0.0015,
#            Bm = 0,
#            Tmax = 5000,
#            Ny = 1) %>% 
#   traject_equil()
#   
# simdatequil <- data.frame(x.init = 1,
#                  y.init = 0.5,
#                  a = seq(0, 0.2, by = 0.01),
#                  # r = 0.2,
#                  # K = 0.9,
#                  c = 0.2,
#                  m = 0.0015,
#                  Bm = 0,
#                  Tmax = 5000,
#                  Ny = 1) %>% 
#   crossing(K = c(0.8, 0.9, 0.99)) %>% 
#   crossing(r = c(0.1, 0.2, 0.4)) %>% 
#   mutate(paramset = 1:n())
# 
# 
# simequil <- trajectory(simdatequil) %>% 
#   left_join(simdatequil) %>% 
#   tibble()
# 
# simequil %>% 
#   group_by(paramset) %>% 
#   filter(round(y, 1)!=round(lag(y),1)) %>% 
#   filter(t == max(t)) %>% 
#   ungroup %>% 
#   count(t) %>% 
#   arrange(desc(t))
# #there are still a fair number of parameter sets not reaching equilibrium
# 
# simequil %>% 
#   ggplot(aes(x = t, y = y, color = factor(r)))+
#   facet_grid(paste("K =", K)~ paste("a =", a))+
#   geom_line()
# 
# 
# simequil %>% 
#   filter(t%in% c(22, Tmax)) %>% 
#   ggplot(aes(x = a, y = y, color = factor(r)))+
#   facet_wrap(paste("t = ", t)~paste("K =", K), scales = "free_y")+
#   geom_line()
# 
# 
# #for these values, midges consumed all algae 
# 
# simequil %>% 
#   select(paramset, K, r, a, t, x, y) %>% 
#   group_by(paramset) %>% 
#   filter(x == 0) %>% 
#   filter(t == min(t)) %>% 
#   print(n = 100)
# #this ocurred only when the density dependence of algae was weak and their growth rates were high.
# 
# 
# simequil %>% 
#   filter(t == Tmax,
#          a>0,
#          x>0) %>% 
#   ggplot(aes(x = exp(r)*x^K, y = y-y.init, col = r))+
#   facet_wrap(~K, scales = "free")+
#     geom_line(aes(group = r))+
#   geom_point()
# 
# 
# simequil %>% 
#   filter(a>0) %>% 
#   filter(t == Tmax) %>%
#   group_by(r, K) %>% 
#   mutate(ydiff = y-max(y)) %>% 
#   # filter(x == 0)
#   ggplot(aes(x = exp(r)*x^K, y = y-y.init, shape = ydiff==0))+
#   facet_wrap(paste("r=", r)~paste("K=",K), scales = "free")+
#   # geom_line(aes(group = K))+
#   geom_point(alpha = 0.5)+
#   scale_shape_manual(values = c(21, 16))
# 
# simequil %>% 
#   filter(a>0) %>% 
#   filter(t == 22) %>%
#   group_by(r, K) %>% 
#   mutate(ydiff = y-max(y)) %>% 
#   # filter(x == 0)
#   ggplot(aes(x = exp(r)*x^K, y = y-y.init, color = a, shape = ydiff==0))+
#   facet_wrap(paste("r=", r)~paste("K=",K), scales = "free")+
#   # geom_line(aes(group = K))+
#   geom_point(alpha = 0.5)+
#   scale_shape_manual(values = c(21, 16))


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


# r        K        a        c        m gpp.rate       SS 
# 0.187    0.927    0.045    0.163    0.002    0.183   66.817 


obssimul <- data.frame(r = 0.187,
                       # K = 0.927,
                       # a = 0.045,
                       c = 0.163,
                       m = 0.0016,
                       Bm =0,
                       Ny = 1,
                       x.init = unique(data$x0),
                       y.init = max(data$y0),
                       Tmax = 2000) %>% 
  crossing(a = seq(0.005, 0.1, by = 0.005)) %>% 
  crossing(K = seq(0.8, 0.95, by = 0.05)) %>% 
  mutate(paramset = 1:n())

obssimv <- trajectory(obssimul) %>% 
  left_join(obssimul) %>% 
  tibble()


obssimv %>% 
  filter(t == 22) %>%
  ggplot(aes(x = a, y = y, color = x.init, group = interaction(x.init, K)))+
  geom_line()+
  geom_point()


obssimv %>% 
  filter(a %in% c(0.005, 0.045, 0.1)) %>% 
  ggplot(aes(x = t, y = y, color = K, group = interaction(x.init, K)))+
  facet_wrap(~a, scales = "free")+
  geom_line()

  


 obssimv %>% 
   filter(t == 22,
          a %in% c(0.005, 0.045, 0.1)) %>%
   group_by(a, K, t) %>% 
   mutate(ygrowth = scale(y-y.init),
          xgrowth = scale(exp(r)*x^K)) %>%
  ggplot(aes(y = ygrowth, x = xgrowth, fill = x.init))+
   facet_grid(paste("\u03b1", "=", a) ~paste("K=", K), scales = 'free')+
   geom_path()+
  geom_point(shape = 21, size =2)+
   labs(x = "Standardized Primary Production",
        y = "Standardized Secondary Production",
        fill = "Starting Resource Density")+
   scale_fill_viridis_c()
 # ggpreview(plot= last_plot(), width = 6, height = 5, units = "in", dpi = 650)
 
 obssimv %>% 
   filter(t %in% c(22, Tmax),
          a %in% c(0.005, 0.045, 0.1),
          x.init == max(x.init)) %>%
   group_by(a, t) %>%
   mutate(t = ifelse(t == Tmax, "Equilibrium", paste("t =", t)),
          t = fct_rev(t),
          ygrowth = scale(y-y.init),
          xgrowth = scale(exp(r)*x^K)) %>%
   ggplot(aes(y = ygrowth, x = xgrowth, fill = K))+
   facet_grid(paste("\u03b1", "=", a)~t, scales = 'free')+
   geom_path()+
   geom_point(shape = 21, size =2)+
   labs(x = "Standardized Primary Production",
        y = "Standardized Secondary Production",
        fill = "K")+
   scale_fill_viridis_c(option = "plasma")
 
 # ggpreview(plot= last_plot(), width = 4, height = 5, units = "in", dpi = 650)
 
 
 
 obssimul2 <- data.frame( K = 0.927,
                        # r = 0.187,
                        # a = 0.045,
                        c = 0.163,
                        m = 0.0016,
                        Bm =0,
                        Ny = 1,
                        x.init = unique(data$x0),
                        y.init = max(data$y0),
                        Tmax = 2000) %>% 
   crossing(a = seq(0.005, 0.1, by = 0.005)) %>% 
   crossing(r = seq(0.05, 0.3, by = 0.05)) %>% 
   mutate(paramset = 1:n())
 
 obssimv2 <- trajectory(obssimul2) %>% 
   left_join(obssimul2) %>% 
   tibble()
 
 
 
 
 obssimv2 %>% 
   filter(t %in% c(22, Tmax)) %>%
   mutate(t = ifelse(t == Tmax, "Equilibrium", paste("t =", t)),
          t = fct_rev(t)) %>%
   ggplot(aes(x = a, y = y-y.init, group = r, color  = r))+
   facet_wrap(~t, scales = "free_y")+
   # facet_grid(paste("\u03b1", "=", a)~t, scales = 'free')+
   geom_path()+
   # geom_point(shape = 21, size =2)+
   labs(x = "",
        y = "Secondary Production",
        fill = "Resource Growth Rate")+
   scale_color_viridis_c(option = "plasma")+
   theme(panel.background = element_rect(fill = "gray80"))
 
 
 
 
 obssimv2 %>% 
   filter(t == 22,
          a %in% c(0.005, 0.045, 0.1)) %>%
   group_by(a, r, t) %>% 
   mutate(ygrowth = scale(y-y.init),
          xgrowth = scale(exp(r)*x^K)) %>%
   ggplot(aes(y = ygrowth, x = xgrowth, fill = x.init))+
   facet_grid(paste("r=", r)~ paste("\u03b1", "=", a), scales = 'free')+
   geom_path()+
   geom_point(shape = 21, size =2)+
   labs(x = "Standardized Primary Production",
        y = "Standardized Secondary Production",
        fill = "Starting Resource Density")+
   scale_fill_viridis_c()
 # ggpreview(plot= last_plot(), width = 6, height = 7, units = "in", dpi = 650)
 
 obssimv2 %>% 
   filter(t %in% c(22, Tmax),
          a %in% c(0.005, 0.045, 0.1),
          x.init == max(x.init)) %>%
   group_by(a, t) %>%
   mutate(t = ifelse(t == Tmax, "Equilibrium", paste("t =", t)),
          t = fct_rev(t),
          ygrowth = scale(y-y.init),
          xgrowth = scale(exp(r)*x^K)) %>%
   ggplot(aes(y = ygrowth, x = xgrowth, fill = r))+
   facet_grid(paste("\u03b1", "=", a)~t, scales = 'free')+
   geom_path()+
   geom_point(shape = 21, size =2)+
   labs(x = "Standardized Primary Production",
        y = "Standardized Secondary Production",
        fill = "Resource Growth Rate")+
   scale_fill_viridis_c(option = "plasma")
 
 # ggpreview(plot= last_plot(), width = 4, height = 5, units = "in", dpi = 650)
 # 
 
 
 
 #########################
 obssimul3 <- data.frame( K = 0.927,
                          # r = 0.187,
                          # a = 0.045,
                          c = 0.163,
                          m = 0.0016,
                          Bm =0,
                          Ny = 1,
                          x.init = mean(data$x0),
                          y.init = max(data$y0),
                          Tmax = 10000) %>% 
   crossing(a = seq(0.0005, 0.01, by = 0.0005)) %>% 
   crossing(r = seq(0.05, 0.3, by = 0.05)) %>% 
   mutate(paramset = 1:n())
 
 obssim_long3 <- trajectory(obssimul3) %>% 
   left_join(obssimul3) %>% 
   tibble()
 
 obssim_long3 %>% 
   ggplot(aes(x = t, y = y, color = a, group = factor(a)))+
   facet_wrap(~r, scales = "free")+
   geom_line()
 
 
 obssimul4 <- data.frame( K = 0.927,
                          # r = 0.187,
                          # a = 0.045,
                          c = 0.163,
                          m = 0.0016,
                          Bm =0,
                          Ny = 1,
                          x.init = unique(data$x0),
                          y.init = max(data$y0),
                          Tmax = 22) %>% 
   crossing(a = seq(0.01, 0.2, by = 0.01)) %>% 
   crossing(r = seq(0.05, 0.3, by = 0.05)) %>% 
   mutate(paramset = 1:n())
 
 obssim_short3 <- trajectory(obssimul4) %>% 
   left_join(obssimul4) %>% 
   tibble()
 
 
obssim_full <-  bind_rows(obssim_long3 %>% mutate(type = "equil"), obssim_short3 %>% mutate(type = "transient"))
 
 
 
 obssim_full %>% 
   filter(t ==Tmax) %>%
   ggplot(aes(x = a, y = y-y.init, group = interaction(r, x.init), color  = r))+
   facet_wrap(~type, scales = "free")+
   # facet_grid(paste("\u03b1", "=", a)~t, scales = 'free')+
   geom_line()+
   # geom_point(shape = 21, size =2)+
   labs(x = "Attack Rate",
        y = "Secondary Production",
        fill = "Resource Growth Rate")+
   scale_color_viridis_c(option = "plasma")+
   theme(panel.background = element_rect(fill = "gray80"))
 
 obssim_short3 %>% 
   filter(t ==Tmax) %>%
   ggplot(aes(x = a, y = y-y.init, group = interaction(r, x.init), fill  = x.init))+
   facet_wrap(~paste("r =", r), scales = "free_y", ncol = 2)+
   # facet_grid(paste("\u03b1", "=", a)~t, scales = 'free')+
   geom_line()+
   geom_point(shape = 21, size =1)+
   labs(x = "Attack Rate",
        y = "Secondary Production",
        fill = "Initial Resource Density")+
   scale_fill_viridis_c()+
   scale_x_continuous(breaks = c(0, 0.1, 0.2))+
   guides(fill = guide_colorbar(title.position = "top"))
 # ggpreview(plot = last_plot(), width = 3, height = 5, units = "in", dpi = 650)
 
 
 
 obssim_short3 %>% 
   group_by(x.init, r) %>%
   filter(t == Tmax) %>% 
   filter(y == max(y)) %>% 
   select(x.init, a, r, y, x, K, y.init) %>% 
   group_by(x.init, r) %>% 
   mutate(ygrowth = (y-y.init),
          xgrowth = (exp(r)*x^K)) %>%
   ggplot(aes(y = ygrowth, x = xgrowth, fill = r))+
   # facet_wrap(~a)+
   labs(x = "Primary Production",
        y = "Secondary Production",
        fill = "Resource\nGrowth Rate")+
   geom_line(aes(group = x.init), alpha = 0.5)+
   geom_point(shape = 21, size = 2)+
   scale_fill_viridis_c(option = "plasma", breaks = c(0.1, 0.2, 0.3))
 
 # ggpreview(plot = last_plot(), width = 3, height = 3.5, units = "in", dpi = 650)
 
 
 obssim_long3 %>% 
   group_by(x.init, r) %>%
   filter(t == Tmax) %>% 
   filter(y == max(y)) %>% 
   select(x.init, a, r, y, x, K, y.init) %>% 
   group_by(x.init, r) %>% 
   mutate(ygrowth = (y-y.init),
          xgrowth = (exp(r)*x^K)) %>%
   ggplot(aes(y = ygrowth, x = xgrowth, fill = r))+
   # facet_wrap(~a)+
   geom_line(aes(group = x.init), alpha = 0.5)+
   geom_point(shape = 21)+
   scale_fill_viridis_c()
 
   
 
 
 obssim_short3 %>% 
   filter(t == 22,
          a %in% c(0.01, 0.05, 0.2)) %>%
   group_by(a, r, t) %>% 
   mutate(ygrowth = scale(y-y.init),
          xgrowth = scale(exp(r)*x^K)) %>%
   ggplot(aes(y = ygrowth, x = xgrowth, fill = x.init))+
   facet_grid(paste("r=", r)~ paste("\u03b1", "=", a), scales = 'free')+
   geom_path()+
   geom_point(shape = 21, size =2)+
   labs(x = "Standardized Primary Production",
        y = "Standardized Secondary Production",
        fill = "Initial Resource Density")+
   scale_fill_viridis_c()
 # ggpreview(plot= last_plot(), width = 6, height = 7, units = "in", dpi = 650)
 
 obssimv2 %>% 
   filter(t %in% c(22, Tmax),
          a %in% c(0.005, 0.045, 0.1),
          x.init == max(x.init)) %>%
   group_by(a, t) %>%
   mutate(t = ifelse(t == Tmax, "Equilibrium", paste("t =", t)),
          t = fct_rev(t),
          ygrowth = scale(y-y.init),
          xgrowth = scale(exp(r)*x^K)) %>%
   ggplot(aes(y = ygrowth, x = xgrowth, fill = r))+
   facet_grid(paste("\u03b1", "=", a)~t, scales = 'free')+
   geom_path()+
   geom_point(shape = 21, size =2)+
   labs(x = "Standardized Primary Production",
        y = "Standardized Secondary Production",
        fill = "Resource Growth Rate")+
   scale_fill_viridis_c(option = "plasma")
 
 # ggpreview(plot= last_plot(), width = 4, height = 5, units = "in", dpi = 650)
 # 