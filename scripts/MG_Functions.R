#====Universal Functions====
#Define function to give NA if all NAs else give mean
meanna <- function(x){
  if(all(is.na(x))){
    NA
  }
  else{
    mean(x, na.rm = TRUE)
  }
}

#Define function to give NA if all NAs else give sum
sumna <- function(x){
  if(all(is.na(x))){
    NA
  }
  else{
    sum(x, na.rm = TRUE)
  }
}

#preview figures prior to saving
ggpreview <- function(...) {
  fname <- tempfile(fileext = ".png")
  ggsave(filename = fname, ...)
  system2("open", fname)
  invisible(NULL)
}

#====Experiment Analysis Functions=====

#find instars
mg_instars = c(0,4.3,7.3,12, 18)/55*1000 #divide by ocular micrometer units to get mm then 1000 to get micrometers
#Instars from Lindegaard: 75 125-140 225-250 325-350 (I'm skeptical of the veracity of the 4th instar measures)


#fix predicted names
fix.names <- function(preddata){
  preddata %>% 
    select(-1) %>% 
    rename(algae_conc2 = 1, midge = 2, box = 3, algae_conc2_x_midge = 4) %>% 
    mutate(midge = ifelse(midge == 1, "Midges", "No Midges"),
           algae_conc2 = exp(algae_conc2))
}

#get estimates and standard errors
mod_predict <- function(model.matrix, model, mixed.effects.model = FALSE){
  
  if(!mixed.effects.model){
    est <- model.matrix %*% coef(model)
    vcov <- as.matrix(vcov(model))
    se <- apply(model.matrix, 1, function(x){sqrt(t(x) %*% vcov %*% x) })
    return(fix.names(data.frame(model.matrix, estimate = est, se = se)))
  }
  if(mixed.effects.model){
    est <- model.matrix %*% fixef(model)
    vcov <- as.matrix(vcov(model))
    se <- apply(model.matrix, 1, function(x){sqrt(t(x) %*% vcov %*% x) })
    return(fix.names(data.frame(model.matrix, estimate = est, se = se)))
  }
  
}

#====Production Functions====

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


#Define variables to scale biomass to same units
ftube_r = 0.03/2 #30 mm/ 2 to get radius /1000 to convert to m
ftube_area = ftube_r^2*pi
pq = 1 #photosynthetic quotient

#convert gpp from O2 mg m^-2 hr-1 to micrograms of C cm^-2 hr^-1
gpp_omgm2h_to_cugcm2d <- function(x){
  y <- (x*18)*(12/32*pq) #daily production of algae in g C m^2d^-1
  y <- y/10 #convert mg C cm^-2d^-1 to ug C cm^-2 d^-1
  return(y)
}

#====Model Functions====
#function to simulate dynamics
trajectory <- function(simdat){ # a dataframe containing columns x.init, y.init,a, K, c, m, Bm, Tmax 
  #list to collect multiple simulations
  trajlist <- list()
  
  for(i in 1:nrow(simdat)){
    x <- simdat$x.init[i]
    y <- simdat$y.init[i]
    r <- simdat$r[i]
    a <- simdat$a[i]
    K <- simdat$K[i]
    c <- simdat$c[i]
    Tmax <- simdat$Tmax[i]
    
    out <- data.frame(paramset = i, t = 1:Tmax,x = NA, y = NA)
    out$x[1] = x
    out$y[1] = y
    for(t in 2:Tmax){
      
      out$x[t] <- max(0,exp(r)*out$x[t-1]^K - a*out$x[t-1]*out$y[t-1]) #growth of algae ##CHECK 
      out$y[t] <- max(0, out$y[t-1] + c*a*out$x[t-1]*out$y[t-1]) #midge growth
      
    }	
    
    trajlist[[i]] <- out
  }
  return(bind_rows(trajlist))
}

#function to fit model to data
SSfit <- function(par, data, x.init, y.init, SS.weight, par.fixed=par.fixed, tofit = T){
  
  par.temp <- par.fixed
  par.temp[is.na(par.fixed)] <- par
  par <- par.temp
  
  r <- exp(par[1])
  K <-  exp(par[2])
  a <-  exp(par[3])
  ac <-  exp(par[4])
  gpp.rate <- exp(par[5])
  
  x <- x.init
  y <- y.init
  
  Tmax <- 22
  
  SS1 <- 0
  SS2 <- 0
  for(t in 1:Tmax){
    
    x.next <- exp(r)*x^K - a*x*y
    y <- y + ac*x*y
    x <- x.next
    
    if(t == 14){
      dif.x <- log(gpp.rate*x) - log(data$x[data$day == 14])
      SS1 <- SS1 + sum(dif.x^2)
      
      sampled.data <- (data$day == 14) & !is.na(data$y) & (data$midge == "Midges")
      sampled.sim <- data$coreid[sampled.data]
      dif.y <- log(y[sampled.sim]) - log(data$y[sampled.data])
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

#converting units from standardized model units to experimental units
convert_to_exp_units <- function(dataframe){
  dataframe %>% 
    mutate(init.chl = (x.init)*unique(meanx0),
           gpp = (x*gpp.rate)*unique(meanx),
           wt = y*unique(meany),
           wtc = wt*0.5*1000,
           gd = (y*meany - y.init*meany)/t,
           gdc = gd*0.5*1000,
           gppc = gpp_omgm2h_to_cugcm2d(gpp)) %>% 
    return()
}

#====set aesthetics====
theme_set(theme_bw()+
            theme(panel.grid = element_blank(),
                  strip.background = element_blank(),
                  legend.position = "top",
                  text = element_text(size = 10),
                  axis.title = element_text(size = 10.5),
                  legend.spacing = unit(0,units = 'points'),
                  legend.margin = margin(c(1,5,5,5)),
                  legend.box.spacing = unit(0, units = "points")))

theme_black = function(base_size = 12, base_family = "") {
  
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    
    theme(
      # Specify axis options
      axis.line = element_blank(),  
      axis.text.x = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),  
      axis.text.y = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),  
      axis.ticks = element_line(color = "white", size  =  0.2),  
      axis.title.x = element_text(size = base_size, color = "white", margin = margin(0, 10, 0, 0)),  
      axis.title.y = element_text(size = base_size, color = "white", angle = 90, margin = margin(0, 10, 0, 0)),  
      axis.ticks.length = unit(0.3, "lines"),   
      # Specify legend options
      legend.background = element_rect(color = NA, fill = "black"),  
      legend.key = element_rect(color = "black",  fill = "black"),  
      legend.key.size = unit(1.2, "lines"),
      legend.key.height = NULL,  
      legend.key.width = NULL,      
      legend.text = element_text(size = base_size*0.8, color = "white"),  
      legend.title = element_text(size = base_size*0.8, face = "bold", hjust = 0, color = "white"),  
      legend.position = "top",  
      legend.text.align = NULL,  
      legend.title.align = NULL,  
      legend.direction = "horizontal",  
      legend.box = NULL, 
      legend.spacing = unit(0, "points"),
      # Specify panel options
      panel.background = element_rect(fill = "black", color  =  NA),  
      panel.border = element_rect(fill = NA, color = "white"),  
      panel.grid.major = element_blank(),  
      panel.grid.minor = element_blank(),  
      panel.margin = unit(0.5, "lines"),   
      # Specify facetting options
      strip.background = element_blank(),  
      strip.text.x = element_text(size = base_size*0.8, color = "white"),  
      strip.text.y = element_text(size = base_size*0.8, color = "white",angle = -90),  
      # Specify plot options
      plot.background = element_rect(color = "black", fill = "black"),  
      plot.title = element_text(size = base_size*1.2, color = "white"),  
      plot.margin = unit(rep(1, 4), "lines"),
      
    )
  
}


midge_color <- scale_color_manual(values = c("Midges" = "black", "No Midges" = "gray60"))
midge_fill <- scale_fill_manual(values = c("Midges" = "black", "No Midges" = "gray60"))
midge_lines <- scale_linetype_manual(values = c("Midges" = "solid", "No Midges" = "dashed"))
algae_color <- saturation(scale_color_viridis_c(trans = "log", breaks = c(0.01, 0.1, 1), end = 0.9), scalefac(5))
algae_fill <-  saturation(scale_fill_viridis_c(trans = "log", breaks = c(0.01, 0.1, 1), end = 0.9), scalefac(5))

midge_color_dark <- scale_color_manual(values = c("Midges" = "yellow", "No Midges" = "red"))
midge_fill_dark <- scale_fill_manual(values = c("Midges" = "yellow", "No Midges" = "red"))
