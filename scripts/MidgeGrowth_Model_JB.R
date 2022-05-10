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
    Tmax <- simdat$Tmax[i]

    out <- data.frame(paramset = i, t = 1:Tmax,x = NA, y = NA)
    out$x[1] = x
    out$y[1] = y
    for(t in 2:Tmax){
      
      out$x[t] <- max(0,exp(r)*out$x[t-1]^K - a*out$x[t-1]*out$y[t-1]) #growth of algae ##CHECK 
      out$y[t] <- max(0, out$y[t-1] + c*a*out$x[t-1]*out$y[t-1] - m) #midge growth
      
    }	
    
    trajlist[[i]] <- out
  }
  return(bind_rows(trajlist))
}

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
 
   
# figures

#====plot model=====

sim <- data.frame(x.init = unique(data$x0),
                  r = exp(par[1]),
                  K = exp(par[2]),
                  a = exp(par[3]),
                  c = exp(par[4])/exp(par[3]), 
                  m = exp(par[5]),
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
# ggpreview(plot = last_plot(), dpi = 650, width = 5, height = 6, units = "in")


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
  ggplot(aes(x = a, y = val, fill = x.init, group = x.init))+
  geom_vline(xintercept = aest)+
  facet_wrap(~var, scales = 'free')+
  geom_line(size = 1.1)+
  geom_line(aes(color = x.init))+
  scale_color_viridis_c(trans = "log", breaks = c(0.25, 1,4))


obsmoda %>% 
  filter(t == Tmax,
         y.init>0) %>% 
  ggplot(aes(x = a, y = y-y.init, fill = x.init, group = x.init))+
  geom_vline(xintercept = aest)+
  geom_line(size = 1.1)+
  geom_line(aes(color = x.init))+
  labs(x = "\u03b1",
       y = "Standardized Secondary Production",
       color = "Iniital Resource Biomass")+
  guides(color = guide_colorbar(title.position = "top", title.hjust = 0.5))+
  scale_color_viridis_c(trans = "log", breaks = c(0.25, 1,4))

# ggpreview(plot = last_plot(), width = 3, height = 4, units = "in", dpi = 650)

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



optims %>% 
  left_join(obsmoda %>% 
              filter(y.init>0,
                     t ==22)) %>% 
  ggplot(aes(x = exp(r)*x^K, y = y-y.init, fill = x.init))+
  geom_point()+
  geom_path()+
  lims(y = c(0, NA))+
  labs(y = "Biomass",
       color = element_blank(),
       linetype = element_blank())
  
  
obsmoda %>% 
  filter(t == Tmax,
         y.init>0,
         a %in% c(0.02, aest, optima, 0.2)) %>% 
  mutate(al = ifelse(a == aest, 
                        "\u03b1 = \u03b1 est", 
                        ifelse(a == optima, 
                               "\u03b1 = \u03b1 opt", 
                               paste("\u03b1 = ", a))),
         al = fct_reorder(al, a, min)) %>% 
  ggplot(aes(x = exp(r)*x^K, y = y-y.init, fill = x.init))+
  facet_wrap(~al, scales = "free", labeller = label_value)+
  geom_path()+
  labs(x = "Standardized Primary Production",
       y = "Standardized Secondary Production",
       fill = "Initial Resource Biomass")+
  geom_point(shape = 21, size = 2)+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))+
  scale_fill_viridis_c(trans = "log", breaks = c(0.25, 1,4))
# ggpreview(plot = last_plot(), width = 3, height = 4, units = "in", dpi = 650)

