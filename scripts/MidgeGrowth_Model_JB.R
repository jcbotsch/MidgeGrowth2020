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


#############################################################
#====Fig. 1 Concept====
#parameters over which to simulate data
simdat1 <- data.frame(y.init = 1,
                       K = 0.9,
                       c = 0.2,
                       m = 0.0015,
                       r = 0.2,
                       Bm = 1,
                       Tmax = 2000,
                       Ny = 1) %>% 
  crossing(a = seq(0,0.05, by = 0.0005),
           x.init = 5) %>% 
  mutate(paramset = 1:n()) 

#run model
traj1 <- trajectory(simdat1) %>% 
  left_join(simdat1) %>% 
  tibble()


aopt <- traj1 %>% 
  filter(t == Tmax) %>% 
  filter(y == max(y)) %>% 
  pull(a)

ymax <- traj1 %>% 
  filter(t == Tmax) %>% 
  filter(y == max(y)) %>% 
  pull(y)
  
traj1 %>% 
  filter(t == Tmax) %>% 
  ggplot(aes(x = a, y = y))+
  geom_vline(xintercept = aopt, alpha = 0.5)+
  geom_hline(yintercept = ymax, alpha = 0.5)+
  geom_text(x = aopt+0.003, y = 1, label = expression(~alpha["opt"]))+
  geom_text(y = ymax-2, x = 0.045, label = expression(~y["max"]))+
  geom_line()+
  labs(x = "Attack Rate of Consumer on Resource",
       y = "Consumer Biomass")+
  theme(axis.text = element_blank())

# ggpreview(plot = last_plot(), width = 3, height = 3, dpi = 650, units = "in")


#====Fig 2. Optimum across parameters====
datax <- data.frame(y.init = 1,
                    x.init = 2,
           K = 0.9,
           c = 0.2,
           m = 0.0015,
           # r = 0.2,
           Bm = 1,
           Tmax = 2000,
           Ny = 1) %>% 
  crossing(a = seq(0,0.025, by = 0.0001),
           r = c(0.1, 0.2, 0.3, 0.4)) %>% 
  mutate(paramset = 1:n()) 

#run model
traj2 <- trajectory(datax) %>% 
  left_join(datax) %>% 
  tibble()


opts <- traj2 %>% 
  group_by(r) %>% 
  filter(t == Tmax,
         y == max(y))

traj2 %>% 
  filter(t == Tmax) %>% 
  ggplot(aes(x = a, y = y, color = r, group = factor(r)))+
  # geom_vline(aes(xintercept = a, color = r), alpha = 0.5, dat = opts)+
  # geom_hline(aes(yintercept = y, color = r), alpha = 0.5, data = opts)+
  geom_point(data = opts)+
  geom_line()+
  labs(x = "Attack Rate of Consumer on Resource",
       y = "Consumer Biomass",
       color = "r")+
  scale_color_gradient(low = "dodgerblue", high = "firebrick4")+
  theme(axis.text = element_blank())
ggpreview(plot = last_plot(), width = 3, height = 3, dpi = 650, units = "in")



#=====Equilibrium====
#find equilibrium values
trajequil <- data.frame(y.init = 1,
                        x.init = 5,
                        K = 0.9,
                        c = 0.2,
                        m = 0.0015,
                        # r = 0.2,
                        Bm = 1,
                        Tmax = 2000,
                        Ny = 1) %>% 
  crossing(a = seq(0,0.5, by = 0.01),
           r = c(0.01, 0.1, 0.2, 0.45, 0.6, 0.8)) %>% 
  mutate(paramset = 1:n()) 

trajequil <- trajectory(trajequil) %>% 
  left_join(trajequil) %>% 
  tibble()

trajequil %>%
  gather(var, val, x, y) %>% 
  filter(a %in% c(0.01, 0.1,0.2, 0.4)) %>% 
  ggplot(aes(x = t, y = val, color = factor(a)))+
  facet_wrap(var~r)+
  geom_line()

trajequil


equils <- trajequil %>% 
  filter(t == Tmax) %>% 
  select(y.init, K, c, m, r, Bm, Ny, a, xequil = x, yequil = y)


varxy <- equils %>% 
  crossing(frac_equily = c(0.1,0.5,0.75,1)) %>% 
  crossing(frac_equilx = c(0.1, 0.5, 0.75, 1)) %>% 
  mutate(Tmax = 20,
         x.init = frac_equilx*xequil,
         y.init = frac_equily*yequil, 
         paramset = 1:n())

trajxy <- trajectory(varxy) %>% 
  left_join(varxy) %>% 
  tibble()
  
trajxy %>% 
  filter(t == Tmax) %>% 
  ggplot(aes(x = a, y = y, color = factor(frac_equilx)))+
  facet_grid(frac_equily~r, scales = "free_y")+
  geom_line()


trajxy %>% 
  filter(t == Tmax,
         a!=0) %>% 
  ggplot(aes(x = exp(r)*x^K, y = y-y.init, color = r))+
  facet_grid(frac_equilx~frac_equily)+
  geom_point()

#====Illustrate optimum in midge feeding====
# Calculate the midge growth as a function of their consumption rate

# x <- simdat$x.init[i]
# y <- simdat$y.init[i]
# r <- simdat$r[i]
# a <- simdat$a[i]
# K <- simdat$K[i]
# c <- simdat$c[i]
# m <- simdat$m[i]
# Bm <- simdat$Bm[i] #lgl (0 or 1)
# Tmax <- simdat$Tmax[i]

#test
simdat1 <- data.frame(paramset = 1,
                      x.init = 0.2,
                      y.init = 1,
                      r = 0.2,
                      a = 0.1,
                      K = 0.9,
                      c = 0.3,
                      m = 0.5,
                      Bm = 1,
                      Tmax = 100,
                      Ny = 2) #0 states that metabolic costs are independent of biomass

trajectory(simdat1) %>% 
  gather(var, val, x, y) %>% 
  ggplot(aes(x = t, y = val, col = var))+
  facet_wrap(~paramset)+
  geom_line()






#run model
trajax <- trajectory(simdatax) %>% 
  left_join(simdatax) %>% 
  tibble()

trajax %>% 
  filter(t == 20) %>% 
  gather(var, val, x, y) %>% 
  ggplot(aes(x = a, y = val, group = x.init, color = x.init))+
  facet_wrap(~var, scales = "free")+
  geom_line()+
  scale_color_viridis_c(option = "plasma")

# ggpreview(plot = last_plot(), width = 3, height = 3, units = "in")


#plots
trajax %>% 
  filter(t %in% c(20,50,100,200)) %>% 
  group_by(t) %>% 
  mutate(z.x = scale(log(x))) %>% 
  ggplot(aes(x = a, y = x.init, fill = z.x))+
  facet_wrap(~t)+
  geom_tile()+
  scale_fill_viridis_c()


trajax %>% 
  filter(t %in% c(20,50,100,200)) %>% 
  group_by(t) %>% 
  mutate(z.y = scale(log(y))) %>% 
  ggplot(aes(x = a, y = x.init, fill = z.y))+
  facet_wrap(~t)+
  geom_tile()+
  scale_fill_viridis_c()


trajax %>% 
  filter(t %in% c(20, 100, 1000, 2000)) %>% 
  filter(a<0.2) %>%
  gather(var, val, x, y) %>% 
  ggplot(aes(x = a, y = val, group = var, color = var))+
  facet_wrap(~t, scales = "free")+
  geom_point()+
  geom_line()
  # scale_color_viridis_c(option = "plasma")


trajax %>% 
  filter(t %in% c(20, 100, 1000, 2000)) %>% 
  filter(a<0.2) %>% 
  ggplot(aes(x = a, y = y, group = x.init, color = x.init))+
  facet_wrap(~t, scales = "free")+
  geom_line()+
  scale_color_viridis_c(option = "plasma")



trajax %>% 
  filter(x.init == 4.6) %>%
  filter(a<0.2) %>% 
  gather(var, val, x, y) %>% 
  ggplot(aes(x = t, y = val, group = a, color = a))+
  facet_wrap(~var, scales= "free_y")+
  geom_line()+
  scale_color_viridis_c(option = "plasma")

# ggpreview(plot = last_plot(), width = 3, height = 3.5, units = "in")

trajax %>% 
  filter(t %in% c(20, 50, 100, 200)) %>% 
  ggplot(aes(x = a, y = x, group = x.init, color = x.init))+
  facet_wrap(~t, scales = "free")+
  geom_line()+
  scale_color_viridis_c(option = "plasma")
# ggpreview(plot = last_plot(), width = 3, height = 3.5, units = "in")


trajax %>% 
  filter(a %in% c(0.01, 0.05, 0.1, 0.2),
         x.init %in% c(0.1,1.1,4.6)) %>% 
  gather(var, val, x, y) %>% 
  ggplot(aes(x = t))+
  facet_grid(var~x.init, scales = "free_y")+
  geom_line(aes(y = val, group = a, color = factor(a)))+
  scale_color_viridis_d()
# ggpreview(plot = last_plot(), dpi = 650, width = 6, height = 3.5, units = "in")



equils <- trajax %>% 
  filter(t == Tmax) %>% 
  select(y.init, K, c, m, r, Bm, Ny, a, xequil = x, yequil = y)


simdatax <- data.frame(K = 0.9,
                       c = 0.3,
                       m = 0.01,
                       r = 0.2,
                       Bm = 0,
                       Tmax = 20,
                       Ny = 1) %>% 
  crossing(a = seq(0,0.4, by = 0.01),
           frac_equily = c(0.1,0.5,0.75, 1, 10)) %>% 
  left_join(equils) %>%
  crossing(frac_equilx = c(0.1, 0.5, 0.75, 1, 10)) %>% 
  mutate(x.init = frac_equilx*xequil,
         y.init = frac_equily*yequil) %>% 
  mutate(paramset = 1:n()) 

#run model
trajax <- trajectory(simdatax) %>% 
  left_join(simdatax) %>% 
  tibble()



trajax %>% 
  left_join(equils) %>% 
  filter(t == 20) %>% 
  gather(var, val, x, y) %>% 
  ggplot(aes(x = a, y = val, color = factor(frac_equilx), group = factor(frac_equilx)))+
  facet_grid(var~frac_equily, scales = "free")+
  geom_line()

trajax %>% 
  left_join(equils) %>% 
  filter(t == 20) %>% 
  ggplot(aes(x = r*x, y = y-y.init, color = a))+
  facet_grid(paste("x =", frac_equilx)~paste( "y = ", frac_equily), scales = "free_y")+
  geom_point()


trajax %>% 
  left_join(equils) %>% 
  filter(t == 20) %>% 
  ggplot(aes(x = exp(r)*x^K, y = y-y.init))+
  # facet_wrap(~a, scales = "free")+
  # facet_grid(paste("x =", frac_equilx)~paste( "y = ", frac_equily), scales = "free_y")+
  geom_line(aes(color = a, group = interaction(a, frac_equilx)), alpha = 0.5)+
  geom_point(aes(fill = frac_equily), shape = 21, alpha = 0.5)+
  geom_smooth(method = "lm")+
  scale_color_viridis_c(option = "plasma")+
  scale_fill_viridis_c()

trajax %>% 
  left_join(equils) %>% 
  filter(t == 20,
         frac_equilx|frac_equily == 1) %>% 
  ggplot(aes(x = exp(r)*x^K, y = y-y.init))+
  # facet_wrap(~a, scales = "free")+
  # facet_grid(paste("x =", frac_equilx)~paste( "y = ", frac_equily), scales = "free_y")+
  geom_line(aes(color = a, group = interaction(a, frac_equilx)), alpha = 0.5)+
  geom_point(aes(fill = frac_equily), shape = 21, alpha = 0.5)+
  geom_smooth(method = "lm")+
  scale_color_viridis_c(option = "plasma")+
  scale_fill_viridis_c()

unique(trajax$y.init)

trajax %>% 
  filter(round(y.init,1)==0.5) %>% 
  count(a, c, )

trajax %>% 
  left_join(equils) %>% 
  filter(t == 20) %>% 
  ggplot(aes(x = exp(r)*x^K, y = y-y.init, group = interaction(a, frac_equilx)))+
  # facet_wrap(~a, scales = "free")+
  # facet_grid(paste("x =", frac_equilx)~paste( "y = ", frac_equily), scales = "free_y")+
  geom_line(aes(color = a), alpha = 0.5)+
  facet_wrap(~frac_equilx)+
  geom_point(aes(fill = frac_equily), shape = 21, alpha = 0.5)+
  scale_color_viridis_c(option = "plasma")+
  scale_fill_viridis_c()




trajax %>% 
  left_join(equils) %>% 
  filter(t == 20) %>% 
  ggplot(aes(x = exp(r)*x^K, y = y-y.init, group = interaction(a, frac_equily)))+
  # facet_wrap(~a, scales = "free")+
  # facet_grid(paste("x =", frac_equilx)~paste( "y = ", frac_equily), scales = "free_y")+
  geom_line(aes(color = a), alpha = 0.5)+
  geom_point(aes(fill = frac_equilx), shape = 21, alpha = 0.5)+
  scale_color_viridis_c(option = "plasma")+
  scale_fill_viridis_c()

trajax %>% 
  left_join(equils) %>% 
  filter(t == 20) %>% 
  ggplot(aes(x = exp(r)*x^K, y = y-y.init, group = interaction(a, frac_equily)))+
  # facet_wrap(~a, scales = "free")+
  # facet_grid(paste("x =", frac_equilx)~paste( "y = ", frac_equily), scales = "free_y")+
  facet_wrap(~frac_equily)+
  geom_line(aes(color = a), alpha = 0.5)+
  geom_point(aes(fill = frac_equilx), shape = 21, alpha = 0.5)+
  scale_color_viridis_c(option = "plasma")+
  scale_fill_viridis_c()

  

trajax %>% 
  filter(a %in% c(0.01, 0.1, 0.2, 0.4)) %>% 
  left_join(equils) %>% 
  # gather(var, val, x, y) %>% 
  ggplot(aes(x = t, y = y, color = a, group = factor(a)))+
  facet_wrap(frac_equilx~frac_equily, scales = "free_y")+
  geom_line()


######

#find equilibrium values
trajequil <- data.frame(y.init = 1,
           x.init = 5,
           K = 0.9,
           c = 0.2,
           m = 0.002,
           Bm = 1,
           Tmax = 2000,
           Ny = 1) %>% 
  crossing(a = seq(0,0.5, by = 0.01),
           r = c(0.01, 0.1, 0.2, 0.45, 0.6, 0.8)) %>% 
  mutate(paramset = 1:n()) 

trajequil <- trajectory(trajequil) %>% 
  left_join(trajequil) %>% 
  tibble()

trajequil %>%
  gather(var, val, x, y) %>% 
  filter(a %in% c(0.01, 0.1,0.2, 0.4)) %>% 
  ggplot(aes(x = t, y = val, color = factor(a)))+
  facet_wrap(var~r)+
  geom_line()

########3333

trajsim <- data.frame(y.init = 1,
           K = 0.9,
           c = 0.2,
           m = 0.002,
           Bm = 1,
           Tmax = 20,
           Ny = 1) %>% 
  crossing(a = seq(0,0.5, by = 0.01),
           x.init = c(0.01, 0.1, 1, 5),
           r = c(0.01, 0.1, 0.2, 0.45, 0.6, 0.8)) %>% 
  mutate(paramset = 1:n()) 

trajsim <- trajectory(trajsim) %>% 
  left_join(trajsim) %>% 
  tibble()



trajsim %>% 
  filter(t == Tmax,
         y.init == 1) %>% 
  rename(Resource = x, Consumer = y) %>% 
  gather(var, val, Resource, Consumer) %>% 
  ggplot(aes(x = a, y = val, group = x.init, color = x.init))+
  facet_wrap(var~paste("r =", r), scales = "free_y", nrow = 2)+
  geom_line()+
  labs(x = "Attack Rate of Consumer on Resource",
       y = "Biomass",
       color = "Initial Resource Biomass")+
  geom_vline(xintercept = 0.045)+
  scale_color_viridis_c(option = "plasma")


trajsim %>% 
  filter(r %in% c(0.01, 0.1, 0.45, 0.8),
         a %in% seq(0, 0.5, by = 0.1)) %>% 
  mutate(r = paste("r = ", r), 
         a = paste("a = ", a)) %>% 
  rename(Resource = x, Consumer = y) %>% 
  gather(var, val, Resource, Consumer) %>% 
  group_by(a, r, x.init, var) %>%
  # mutate(val = val/mean(val)) %>%
  ggplot(aes(x = t, y = val, group = interaction(var, x.init), color = x.init))+
  facet_wrap(r~a, scales = "free_y", ncol = 5)+
  geom_line(aes(linetype = var))+
  labs(y = "Biomass",
       linetype = element_blank(),
       color = "Initial Resource Biomass")+
  scale_color_viridis_c(option= "plasma")

trajsim %>% 
  mutate(r = paste("r = ", r), 
         a = paste("a = ", a)) %>% 
  ggplot(aes(x = t, y = y, group = x.init, color = x.init))+
  facet_wrap(r~a, scales = "free_y", ncol = 7)+
  geom_line()+
  scale_color_viridis_c(option= "plasma", trans = "log")

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




#############################################################

#======Function to fit model=======
SSfit <- function(par, data, x.init, y.init, SS.weight, par.fixed=par.fixed, tofit = T){
	
	par.temp <- par.fixed
	par.temp[is.na(par.fixed)] <- par
	par <- par.temp

	r <- exp(par[1])
	K <-  exp(par[2])
	a <-  exp(par[3])
	ac <-  exp(par[4])
	m <-  exp(par[5])
	gpp.rate <-  exp(par[6])

	x <- x.init
	y <- y.init
	
	Tmax <- 22
	
	SS1 <- 0
	SS2 <- 0
	for(t in 1:Tmax){
	  
		x.next <- exp(r)*x^K - a*x*y
		y <- y + ac*x*y - m
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
round(c(r=r, K=K, a=a, c=c, m=m, gpp.rate=gpp.rate, SS=opt.opt$value), digits=5)
       # r        K        a        c        m gpp.rate       SS 
   # 0.187    0.927    0.045    0.163    0.002    0.183   66.817 
 
   
# figures

#====plot model=====

sim <- data.frame(x.init = unique(data$x0),
                  r = exp(par[1]),
                  K = exp(par[2]),
                  a = exp(par[3]),
                  c = exp(par[4])/exp(par[3]), 
                  m = exp(par[5]),
                  Bm = 0,
                  Ny = 1,
                  Tmax = 22) %>% 
  crossing(y.init = unique(data$y0)) %>% 
  mutate(paramset = 1:n())

obsmod <- trajectory(sim) %>% 
  left_join(sim) %>% 
  left_join(init.data %>% 
              select(-coreid) %>% 
              unique() %>% 
              rename(x.init = x0, y.init=y0) %>% 
              mutate(midge = ifelse(y.init>0, "Midges", "No Midges")))



obsmod %>% 
  ggplot(aes(x = t, y = val, col = var, linetype = midge, shape = midge))+
  facet_wrap(~round(algae_conc2, 3), ncol = 2)+
  geom_line(aes(y = y, col = "Midge"))+
  geom_line(aes(y = x*gpp.rate, col = "Algae"))+
  geom_point(aes(y = val, color = taxon, x = day), alpha = 0.5, position = position_dodge(width = 2),  data = data %>% rename(Midge = y, Algae = x) %>%  gather(taxon, val, Midge, Algae))+
  # scale_y_continuous("Algae", sec.axis = sec_axis(~./5, "Midges"))+
  scale_linetype_manual(values = c("solid", "dotted"))+
  scale_shape_manual(values = c(16,21))+
  scale_color_manual(values = c("black", "red"))+
  labs(linetype = element_blank(),
       color = element_blank(),
       shape = element_blank(),
       x = "Time",
       y = "Scaled Biomass")
# ggpreview(plot = last_plot(), dpi = 650, width = 5, height = 6, units = "in")




obsmod %>% 
  ggplot(aes(x = t, y = val, col = var, linetype = midge, shape = midge))+
  facet_wrap(~round(algae_conc2, 3), ncol = 2)+
  geom_line(aes(y = y, col = "Midge"))+
  geom_line(aes(y = x*gpp.rate, col = "Algae"))+
  scale_linetype_manual(values = c("solid", "dotted"))+
  scale_shape_manual(values = c(16,21))+
  scale_color_manual(values = c("black", "red"))+
  labs(linetype = element_blank(),
       color = element_blank(),
       shape = element_blank(),
       x = "Time",
       y = "Scaled Biomass")


obsmod %>% 
  mutate(x = x*gpp.rate) %>% 
  gather(var, val, x, y) %>% 
  filter(val!=0) %>% 
  ggplot(aes(x = t, y = val, col = algae_conc2, group = algae_conc2, shape = midge))+
  facet_wrap(var~midge, scales = "free_y", ncol = 2)+
  geom_line()+
  scale_color_viridis_c(trans = "log", breaks = c(0.01, 0.1, 1))+
  labs(linetype = element_blank(),
       color = "Sediment Treatment",
       shape = element_blank(),
       x = "Time",
       y = "Scaled Biomass")
# ggpreview(plot = last_plot(), dpi = 650, width = 3, height = 4, units = "in")




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

obsmoda <- trajectory(sima) %>% 
  left_join(sima) %>% 
  left_join(init.data %>% 
              select(-coreid) %>% 
              unique() %>% 
              rename(x.init = x0, y.init=y0) %>% 
              mutate(midge = ifelse(y.init>0, "Midges", "No Midges")))


aest <- a
obsmoda %>% 
  filter(t == Tmax,
         y.init>0) %>% 
  gather(var, val, x, y) %>% 
  ggplot(aes(x = a, y = val, color = x.init, group = x.init))+
  geom_vline(xintercept = aest)+
  facet_wrap(~var, scales = 'free')+
  geom_line()+
  scale_color_gradient(trans = "log2", low = "goldenrod", high = "darkgreen")



#
sim_tsa <- data.frame(x.init = unique(data$x0), 
                      y.init = max(unique(data$y0)),
                      r = exp(par[1]),
                      K = exp(par[2]),
                      # a = exp(par[3]),
                      c = exp(par[4])/exp(par[3]), 
                      m = exp(par[5]),
                      Bm = 0,
                      Ny = 1,
                      Tmax = 5000) %>% 
  crossing(a = seq(0, 0.3, by = 0.01)) %>% 
  mutate(paramset = 1:n())
  

sim_tsa <- trajectory(sim_tsa) %>% 
  left_join(sim_tsa) %>% 
  tibble()

sim_tsa %>% 
  filter(t %in% c(22, 65, 5000),
       y.init>0) %>% 
  rename(Midges = y, Algae = x) %>%
  mutate(t2 = ifelse(t == 5000, "Equilibrium", paste("t =", t)),
         t3 = fct_reorder(t2, t, .desc = FALSE) ) %>% 
  ggplot(aes(x = a, y = Midges-y.init, group = x.init))+
  geom_vline(xintercept = aest)+
  facet_wrap(~t3, scales = 'free_y', ncol = 1)+
  scale_x_continuous(breaks = c(0, 0.15, 0.3))+
  geom_line(alpha = 0.5)+
  geom_point(aes(fill = x.init), alpha = 0.7, shape = 21)+
  labs(x = "Attack Rate",
       y = "Secondary Production", 
       fill = "Initial Resource Biomass")+
  scale_fill_viridis_c()+
  guides(fill = guide_colorbar(title.position = "top"))
# ggpreview(plot = last_plot(), width = 3, height = 4, units = "in", dpi = 650)


sim_tsa %>% 
  filter(t %in% c(22),
         y.init>0) %>% 
  rename(Midges = y, Algae = x) %>%
  mutate(t2 = ifelse(t == 5000, "Equilibrium", paste("t =", t)),
         t3 = fct_reorder(t2, t, .desc = FALSE) ) %>% 
  ggplot(aes(x = a, y = Midges-y.init, group = x.init))+
  geom_vline(xintercept = aest)+
  scale_x_continuous(breaks = c(0, 0.15, 0.3))+
  geom_line(alpha = 0.5)+
  geom_point(aes(fill = x.init), alpha = 0.7, shape = 21)+
  labs(x = "Attack Rate",
       y = "Secondary Production", 
       fill = "Initial Resource Biomass")+
  scale_color_gradient(low = "goldenrod", high = "darkgreen")+
  scale_fill_viridis_c()+
  guides(fill = guide_colorbar(title.position = "top"))

# ggpreview(plot = last_plot(), width = 3, height = 4, units = "in", dpi = 650)




sim_tsa %>% 
  filter(t %in% c(22),
         y.init>0, 
         a %in% c(0.01, 0.05, 0.08, 0.15)) %>%
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


sim_tsa %>% 
  filter(t %in% c(65),
         y.init>0, 
         a %in% c(0.01, 0.05, 0.08, 0.15)) %>%
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

ggpreview(plot = last_plot(), width = 3, height = 4, units = "in", dpi = 650)


round(aest, 2)

sim_tsa %>% 
  filter(t<100) %>% 
  group_by(t) %>% 
  filter(y == max(y)) %>% 
  filter(a == 0.04) %>% 
  ungroup %>% 
  filter(y == max(y))



#longer timescale
sim_ts <- data.frame(x.init = unique(data$x0),
                     r = r,
                     K = K,
                     a = a,
                     c = c, 
                     m = m,
                     Bm = 0,
                     Ny = 1,
                     Tmax = 10000) %>% 
  crossing(y.init = unique(data$y0)) %>% 
  mutate(paramset = 1:n())


sim_ts <- trajectory(sim_ts) %>% 
  left_join(sim_ts) %>% 
  left_join(init.data %>% 
              select(-coreid) %>% 
              unique() %>% 
              rename(x.init = x0, y.init=y0) %>% 
              mutate(midge = ifelse(y.init>0, "Midges", "No Midges"))) 


sim_ts %>% 
  gather(var, val, x, y) %>% 
  filter(!(var=="y"& midge == "No Midges")) %>% 
  ggplot(aes(x = t, y = val, group = algae_conc2, col = algae_conc2))+
  facet_wrap(var~midge, scales = "free_y", ncol=2)+
  geom_line()+
  scale_linetype_manual(values = c("solid", "dotted"))+
  scale_color_viridis_c(trans = "log", breaks = c(0.01, 0.1, 1))+
  labs(linetype = element_blank(),
       color = element_blank(),
       shape = element_blank(),
       x = "Time",
       y = "Scaled Biomass")
# ggpreview(plot = last_plot(), dpi = 650, width = 3, height = 4, units = "in")

sim_ts %>% 
  filter(t == Tmax) %>% 
  rename(y.equil = y,
         x.equil = x) %>% 
  mutate(y.frac = y.init/y.equil,
         x.frac = x.init/x.equil) %>% 
  ungroup %>% 
  select(algae_conc2, contains("x."), contains("y.")) %>% 
  mutate_all(round, 3)

obsmod

#===Add Midges====

N.range <- c(1,5,10)

sim_N <- data.frame(x.init = unique(data$x0),
                    y.init = 0.734, #scaled initial midge weight
           r = r,
           K = K,
           a = a,
           c = c, 
           m = m,
           Bm = 0,
           Tmax = 22) %>% 
  crossing(Ny = N.range) %>% 
  mutate(paramset = 1:n())


sim_N <- trajectory(sim_N) %>% 
  left_join(sim_N) %>% 
  left_join(init.data %>% 
              select(algae_conc2, x0) %>% 
              unique() %>% 
              rename(x.init = x0)) 

sim_N %>% 
  mutate(x = x*gpp.rate) %>% 
  rename(Algae = x, Midges = y) %>% 
  gather(midge, val, Algae, Midges) %>% 
  ggplot(aes(x = t, y = val, group = algae_conc2, col = algae_conc2))+
  facet_wrap(midge~Ny, scales = "free_y")+
  geom_line()+
  scale_color_viridis_c(trans = "log", breaks = c(0.01, 0.1, 1))+
  labs(color = "Initial Algal Density",
       x = "Time",
       y = "Scaled Biomass")

# ggpreview(plot = last_plot(), dpi = 650, width = 6, height = 4, units = "in")



# Same but longer
sim_Nlong <- data.frame(x.init = unique(data$x0),
                    y.init = 0.734, #scaled initial midge weight
                    r = r,
                    K = K,
                    a = a,
                    c = c, 
                    m = m,
                    Bm = 0,
                    Tmax = 200) %>% 
  crossing(Ny = N.range) %>% 
  mutate(paramset = 1:n())


sim_Nlong <- trajectory(sim_Nlong) %>% 
  left_join(sim_Nlong) %>% 
  left_join(init.data %>% 
              select(algae_conc2, x0) %>% 
              unique() %>% 
              rename(x.init = x0)) 

sim_Nlong %>% 
  mutate(x = x*gpp.rate) %>% 
  rename(Algae = x, Midges = y) %>% 
  gather(midge, val, Algae, Midges) %>% 
  ggplot(aes(x = t, y = val, group = algae_conc2, col = algae_conc2))+
  facet_wrap(midge~Ny, scales = "free_y")+
  geom_line()+
  scale_color_viridis_c(trans = "log", breaks = c(0.01, 0.1, 1))+
  labs(color = "Initial Algal Density",
       x = "Time",
       y = "Scaled Biomass")
# ggpreview(plot = last_plot(), dpi = 650, width = 6, height = 4, units = "in")



#====Midge Density and Engineering=====
##Add Engineering
#engineering effect should be saturating
#herren found for every 100 midges in a mesocosm gpp increased by 1.5% 
#phillips found that at high densities consumption exceeds algal growth so r should saturate


#per capita engineering effect as a proportion of the r

eesat <- function(N, maxee){
  maxee*(N-1)/(5+N-1)
}
eesat(20, 0.5)

sim_Nee <- data.frame(x.init = unique(data$x0),
                    y.init = 0.734, #scaled initial midge weight
                    K = K,
                    a = a,
                    c = c, 
                    m = m,
                    Bm = 0,
                    Tmax = 22) %>% 
  crossing(Ny = N.range) %>% 
  mutate(r = r+eesat(Ny, 0.2)) %>% 
  mutate(paramset = 1:n())


sim_Nee <- trajectory(sim_Nee) %>% 
  left_join(sim_Nee) %>% 
  left_join(init.data %>% 
              select(algae_conc2, x0) %>% 
              unique() %>% 
              rename(x.init = x0)) 

sim_Nee %>% 
  mutate(x = x*gpp.rate) %>% 
  rename(Algae = x, Midges = y) %>% 
  gather(midge, val, Algae, Midges) %>% 
  ggplot(aes(x = t, y = val, group = algae_conc2, col = algae_conc2))+
  facet_wrap(midge~Ny, scales = "free")+
  geom_line()+
  scale_color_viridis_c(trans = "log", breaks = c(0.01, 0.1, 1))+
  labs(color = "Initial Algal Density",
       x = "Time",
       y = "Scaled Biomass")

# ggpreview(plot = last_plot(), dpi = 650, width = 6, height = 4, units = "in")


ee_comp <- sim_N %>% 
  select(Ny, x.init, t, x, y, algae_conc2) %>% 
  rename(x.ne = x,
         y.ne = y) %>% 
  left_join(sim_Nee %>% select(Ny, x.init, t, x, y) %>% 
              rename(x.ee = x, y.ee = y))




ee_comp %>% 
  ggplot(aes(x = x.ne, y = x.ee))+
  facet_wrap(~round(x.init,3))+
  geom_point(aes(col = factor(Ny)))

ee_comp %>% 
  ggplot(aes(x = y.ne, y = y.ee))+
  facet_wrap(~round(x.init,3), scales = "free")+
  geom_point(aes(col = factor(Ny)))
  

ee_comp %>% 
  mutate(x.diff = (x.ee - x.ne)/((x.ee+x.ne)/2),
         y.diff = (y.ee-y.ne)/((y.ee+y.ne)/2)) %>% 
  gather(var, diff, contains("diff")) %>% 
  filter(Ny>1) %>% 
  ggplot(aes(x = t, y = diff, color = algae_conc2))+
  facet_grid(var~Ny, scales = "free")+
  geom_line(aes(group = algae_conc2))+
  scale_color_viridis_c(trans = "log", breaks = c(0.01, 0.1, 1))


ee_comp2 <- ee_comp %>% 
  gather(var, val, x.ee, y.ee, x.ne, y.ne) %>% 
  mutate(ee = substr(var, 3,4),
         var = substr(var,1,1)) 


ee_comp2%>% 
  filter(Ny>1) %>% 
  ggplot(aes(x = t, y = val, linetype = ee, color = algae_conc2, group = interaction(ee, algae_conc2)))+
  facet_grid(var~Ny, scales = "free")+
  geom_line()+
  scale_color_viridis_c(trans = "log", breaks = c(0.01, 0.1, 1))+
  scale_shape_manual(values = c(16, 21))

# ggpreview(plot = last_plot(), dpi = 650, width = 5, height = 4, units = "in")


#=====Dynamics??=====

#egg weight (from Lindegaard)
egg <- 0.003/mean(data$wt, na.rm = TRUE)

#fourth instar weight (from my analyses of the station 33 data)
adult <- 0.23/mean(data$wt, na.rm = TRUE)



gents <- data.frame(x.init = max(data$x0),
           y.init = egg, #scaled initial midge weight
           r = r,
           Ny = 1,
           K = K,
           a = a,
           c = c, 
           m = m,
           Bm = 0,
           Tmax = 200) 

genlength <- trajectory(gents) 

genlength %>% 
  filter(y>adult) %>% 
  filter(t == min(t))

genlength%>% 
  ggplot(aes(x = t, y = y))+
  geom_line()+
  geom_hline(yintercept = adult)


genlength %>% 
  ggplot(aes(x = t, y = x))+
  geom_line()+
  geom_vline(xintercept = 180)


trajectory(data.frame(x.init = 0.375,
           y.init = egg, #scaled initial midge weight
           r = r,
           Ny = 1,
           K = K,
           a = a,
           c = c, 
           m = m,
           Bm = 0,
           Tmax = 200) ) %>% 
  gather(var, val, x, y) %>% 
  ggplot(aes(x = t, y = val))+
  facet_wrap(~var)+
  geom_line()+
  geom_hline(yintercept = adult)+
  geom_vline(xintercept = 180)

maxlength <- 1000

traj_gen <- function(simdat){ # a dataframe containing columns x.init, y.init,a, K, c, m, Bm, Tmax 
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
    
    out <- data.frame(paramset = i, t = 1:Tmax,x = NA, y = NA, gen = NA, gent = NA)
    out$x[1] = x
    out$y[1] = y
    out$gen[1] = 1
    out$gent[1] = 1
    for(t in 2:Tmax){
      out$x[t] <- max(0,exp(r)*out$x[t-1]^K - Ny*(a*out$x[t-1]*out$y[t-1])) #growth of algae ##CHECK 
      
      if((round(out$gent[t-1]/maxlength))==(out$gent[t-1]/maxlength)){
        out$gent[t] <- 1
        out$gen[t] <- out$gen[t-1]+1
        out$y[t] <- ifelse(out$y[t-1]>=adult, egg, 
                           ifelse(out$y[t-1]==0, egg, 0))
        
      }
      else{
        if(out$y[t-1]>=adult){
          out$y[t] <- egg
          out$gent[t] <- 1
          out$gen[t] <- out$gen[t-1]+1
        }
        else{
          out$y[t] <- max(0, out$y[t-1] + c*a*out$x[t-1]*out$y[t-1] - (1-Bm)*m - Bm*m*out$y[t-1] )
          out$gent[t] <- out$gent[t-1]+1
          out$gen[t] <- out$gen[t-1]
        }
    
        
      }
      

    }	
    
    trajlist[[i]] <- out
  }
  return(bind_rows(trajlist))
}


dyn_v1 <- traj_gen(data.frame(x.init = 1,
           y.init = egg, #scaled initial midge weight
           r = r,
           Ny = 1,
           K = K,
           a = a,
           c = c, 
           m = m,
           Bm = 0,
           Tmax = 2000)) 

max(dyn_v1$gen)

dyn_v1 %>% 
  group_by(gen) %>% 
  filter(gent==max(gent))

dyn_v1 %>% 
  mutate(x = x*gpp.rate) %>% 
  gather(var, val, x, y) %>% 
  ggplot(aes(x = t, val, color = var))+
  geom_line()+
  scale_color_manual(values = c("black","red"))
# ggpreview(plot = last_plot(), dpi = 650, width = 4, height = 4, units = "in")



#do initial conditions alter the dynamics?

dyn_ic <- data.frame(r = r,
           Ny = 1,
           K = K,
           a = a,
           c = c, 
           m = m,
           Bm = 0,
           Tmax = 2000) %>% 
  crossing(x.init = c(min(init.data$x0), max(init.data$x0)),
           y.init = c(egg, max(init.data$y0), adult/2)) %>% 
  mutate(paramset = 1:n())

dyn_ic <- traj_gen(dyn_ic) %>% 
  left_join(dyn_ic)


dyn_ic %>% 
  mutate(x = x*gpp.rate) %>% 
  gather(var, val, x, y) %>% 
  ggplot(aes(x = t, y = val, col = var, ))+
  geom_hline(yintercept = adult)+
  facet_grid(round(y.init,3)~round(x.init,3))+
  geom_line()+
  scale_color_manual(values = c("black", "red"))
# ggpreview(plot = last_plot(), dpi = 650, width = 5, height = 4, units = "in")


dyn_ic %>% 
  group_by(gen, x.init, y.init) %>% 
  filter(gent == max(gent)) %>% 
  group_by(x.init, y.init) %>% 
  filter(gen !=max(gen)) %>% 
  count(gent)

#====altered growth rates====

dyn_v1 <- traj_gen(data.frame(x.init = 1,
                              y.init = egg, #scaled initial midge weight
                              r = r,
                              Ny = 1,
                              K = K,
                              a = a,
                              c = c, 
                              m = m,
                              Bm = 0,
                              Tmax = 2000)) 

max(dyn_v1$gen)

dyn_v1 %>% 
  group_by(gen) %>% 
  filter(gent==max(gent))

dyn_v1 %>% 
  mutate(x = x*gpp.rate) %>% 
  gather(var, val, x, y) %>% 
  ggplot(aes(x = t, val, color = var))+
  geom_line()+
  scale_color_manual(values = c("black","red"))
# ggpreview(plot = last_plot(), dpi = 650, width = 4, height = 4, units = "in")



#do initial conditions alter the dynamics?

dyn_icr <- data.frame(r = c(r*0.5, r),
                     Ny = 1,
                     K = K,
                     a = a,
                     c = c, 
                     m = m,
                     Bm = 0,
                     Tmax = 2000) %>% 
  crossing(x.init = c(min(init.data$x0), max(init.data$x0)),
           y.init = c(max(init.data$y0))) %>% 
  mutate(paramset = 1:n())

dyn_icr <- traj_gen(dyn_icr) %>% 
  left_join(dyn_icr)


dyn_icr %>% 
  mutate(x = x*gpp.rate) %>% 
  gather(var, val, x, y) %>% 
  ggplot(aes(x = t, y = val, col = var, ))+
  geom_hline(yintercept = adult)+
  facet_grid(paste0("r = ",round(r,3))~paste0("x.init = ", round(x.init,3)))+
  geom_line()+
  scale_color_manual(values = c("black", "red"))
# ggpreview(plot = last_plot(), dpi = 650, width = 5, height = 4, units = "in")

dyn_icr %>% 
  count(y==0)


dyn_icr %>% 
  group_by(r, x.init) %>% 
  summarise(max.x = max(x),
            min.x = min(x),
            max.y = max(y),
            min.y = min(y))

dyn_icr %>% 
  group_by(r, x.init) %>% 
  summarise(max.x = max(x),
            min.x = min(x))

dyn_icr %>% 
  filter(t>500, t<800) %>% 
  ggplot(aes(x = t, y = y))+
  facet_grid(paste0("r = ",round(r,3))~paste0("x.init = ", round(x.init,3)))+
  geom_line()

dyn_icr %>% 
  filter(t>700, t<800) %>% 
  ggplot(aes(x = t, y = y))+
  facet_grid(paste0("r = ",round(r,3))~paste0("x.init = ", round(x.init,3)))+
  geom_line()+
  geom_vline(xintercept = 750)

######

#====Spatially Explicit IBM=====
#simulate a landscape with resource variation

#landscape size
landsize = 5

res0 <- tibble(loc_x = 1:landsize) %>% 
  crossing(loc_y = 1:landsize) %>% 
  mutate(x = runif(n(), min = 0, max = 6),
         t = 1) #randomly generate variability in algal abundance through space (uniformly)

#simulate individuals

#number of individuals to start
ny.start <- 25
mid0 <- tibble(y = rep.int(egg, times = ny.start), #all start as eggs
               loc_x = sample(1:landsize, ny.start, replace = TRUE),
               loc_y = sample(1:landsize, ny.start, replace = TRUE),
               t.born = 1,
               t = 1) %>% 
  mutate(y.id = paste(t.born, row_number(), sep = "."))
  


#plot initial state of environment
res0 %>% 
  left_join(mid0 %>% 
              group_by(loc_x, loc_y) %>% 
              arrange(loc_x, loc_y) %>% 
              count(name = "Ny")) %>% 
  ggplot(aes(x = loc_x, y = loc_y, fill = x)) + 
  geom_tile() +
  geom_point(aes(size = Ny))+
  coord_equal()+
  scale_fill_viridis_c(option = "plasma")
# ggpreview(plot = last_plot(), dpi = 650, width = 4, height = 4, units = "in")



#====Function Describing Births====
#simulate births as a poisson process
birth.lambda = 4 #choose lambda alue

repro <- function(ydf, lambda, t){
  #get number of individuals large enough to reproduce
  repro_size <- ydf %>% 
    filter(y >= adult) %>% 
    nrow()
  
  if(repro_size>0){
    #get number of new individuals
    newY <- sum(rpois(repro_size, birth.lambda))
    
    #create data frame with characteristics
    newMid <- tibble(y = rep.int(egg, times = newY),
                     loc_x = sample(1:landsize, newY, replace = TRUE),
                     loc_y = sample(1:landsize, newY, replace = TRUE),
                     t = t, 
                     t.born = t) %>% 
      mutate(y.id = paste(t.born, row_number(), sep = "."))
    
    #join to old data frame
    full <- bind_rows(ydf, newMid)
    return(full)
  }
  else{
    return(ydf)
  }
  
}


#====Function Describing Deaths====
lifespan <- 300
deathprob <- 0.9


death <- function(ydf, lifespan){
  survivors <- ydf %>% 
    rowwise() %>% 
    mutate(age = t-t.born,
           death = ifelse(age>lifespan|y>adult, 
                          rbinom(n = 1, size =1, prob= deathprob), 0)) %>% 
    filter(death<1) %>% 
    select(-age, -death)
  
  return(survivors)
  
}


#I will not simulate movement


#====Simulate Dynamics====

#for each pixel on each timestep calculate the loss of algae as the sum of the attack constant * xi * y
# for each individual midge growth is a function of the algae at their spot at the previous timestep
#for each pixel algal growth is calculated as above with the loss 
# for each timestep births and deaths must be calculated before the time progresses
#additionally at each timestep I want to know the total biomass of midges at a pixel and the total number of adults 



xlist <- list()
xlist[[1]] <- res0

ylist <- list()
ylist[[1]] <- mid0

Tmax <- 2000
for(i in 2:Tmax){
  xin <- xlist[[i-1]]
  yin <- ylist[[i-1]]
  
  xout <- xin %>% 
    full_join(yin %>% select(loc_x, loc_y, y), by = c("loc_x", "loc_y")) %>%   
    mutate(xlag = x,
           xgrowth = exp(r)*xlag^K, 
           t = i,
           y = ifelse(is.na(y), 0, y)) %>% 
    group_by(loc_x, loc_y, t) %>% 
    summarise(xloss = sum(a*xlag*y),
              x = xgrowth-xloss) %>% 
    select(loc_x, loc_y, t, x) %>% 
    unique() %>% 
    ungroup

  
  yout <- yin %>% 
    left_join(xin %>% select(loc_x, loc_y, x), by = c("loc_x", "loc_y")) %>%
    mutate(ylag = y, 
           y = max(0, ylag + c*a*x*ylag - m), 
           t = i) %>% 
    select(loc_x, loc_y, y.id, t.born, t, y)
  yout <- death(repro(yout, birth.lambda, i), lifespan = lifespan)
  
  xlist[[i]] <- xout
  ylist[[i]] <- yout
  print(i)
}
# xlist
# ylist

y_ibm <- bind_rows(ylist) 

x_ibm <- bind_rows(xlist)




Ny_ibm <- y_ibm %>% 
  group_by(loc_x, loc_y, t) %>% 
  add_count(name = "Ny") %>% 
  group_by(loc_x, loc_y, t, Ny) %>% 
  summarise(avg.born = mean(t.born),
            avg.size = mean(y))

#====Plot dynamics of system====
y_ibm %>% 
  group_by(t) %>% 
  count(name = "Ny") %>% 
  ggplot(aes(x = t, y = Ny/10))+
  geom_line(aes(col = "Number of midges"))+
  geom_line(aes(y = avg.x, color = "Average Algal Biomass"), data = x_ibm %>% group_by(t) %>% summarise(avg.x = mean(x)))+
  scale_color_manual(values = c("black", "red"))+
  scale_y_continuous(name = "Average Algal Biomass", sec.axis = sec_axis(trans = ~.*10, name ="Number of Midges"))+
  theme(axis.title.y.right = element_text(color = "red"))+
  labs(color = element_blank())

# ggpreview(plot = last_plot(), dpi = 650, width = 5, height = 3.5, units = "in")



y_ibm %>% 
  group_by(t) %>% 
  summarise(Bt = sum(y)) %>% 
  ggplot(aes(x = t, y = Bt))+
  geom_line(aes(col = "Midge"))+
  geom_line(aes(y = tot.x, color = "Algae"), data = x_ibm %>% group_by(t) %>% summarise(tot.x = sum(x)))+
  scale_color_manual(values = c("black", "red"))+
  scale_y_continuous(name = "Total Biomass")+
  theme(axis.title.y.right = element_text(color = "red"))+
  labs(color = element_blank())
# ggpreview(plot = last_plot(), dpi = 650, width = 5, height = 3.5, units = "in")



#===Dynamics in individual systems====
x_ibm %>% 
  full_join(y_ibm) %>% 
  ggplot(aes(x = t, y = x))+
  facet_grid(loc_x~loc_y)+
  geom_line()+
  geom_line(aes(y = y*2, group = y.id, col = "midge"), col = "red")


x_ibm %>% 
  full_join(y_ibm ) %>% 
  group_by(t, loc_x, loc_y, x) %>% 
  summarise(By = sum(y, na.rm = TRUE)) %>% 
  ggplot(aes(x = t, y = x))+
  facet_grid(loc_x~loc_y)+
  geom_line()+
  geom_line(aes(y = By*1.5, col = "midge"), col = "red")+
  labs(y = "Scaled Biomass")
# ggpreview(plot = last_plot(), dpi = 650, width = 10, height = 6, units = "in")


x_ibm %>% 
  full_join(y_ibm) %>% 
  group_by(t, loc_x, loc_y, x) %>% 
  count(name = "Ny") %>% 
  ggplot(aes(x = t, y = x))+
  facet_grid(loc_x~loc_y)+
  geom_line()+
  geom_line(aes(y = Ny*2, col = "midge"), col = "red")+
  scale_y_continuous(name = "Average Algal Biomass", sec.axis = sec_axis(trans = ~./2, name ="Number of Midges"))+
  theme(strip.placement = "outside")
# ggpreview(plot = last_plot(), dpi = 650, width = 10, height = 6, units = "in")


y_ibm %>% 
  group_by(y.id) %>% 
  mutate(adult = ifelse(y>6.5, 1, 0)) %>% 
  group_by(t) %>% 
  summarise(n.adults = sum(adult)) %>% 
  ggplot(aes(x = t, y = n))+
  geom_line()
