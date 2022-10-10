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

#preview
ggpreview <- function(...) {
  fname <- tempfile(fileext = ".png")
  ggsave(filename = fname, ...)
  system2("open", fname)
  invisible(NULL)
}

#====Experiment Analysis Functions=====

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
#function to fit dynamics
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
    # m <- simdat$m[i]
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
  # m <-  exp(par[5])
  # gpp.rate <-  exp(par[6])
  gpp.rate <- exp(par[5])
  
  x <- x.init
  y <- y.init
  
  Tmax <- 22
  
  SS1 <- 0
  SS2 <- 0
  for(t in 1:Tmax){
    
    x.next <- exp(r)*x^K - a*x*y
    # y <- y + ac*x*y - m
    y <- y + ac*x*y
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

convert_to_exp_units <- function(dataframe){
  dataframe %>% 
    mutate(init.chl = (x.init)*unique(meanx0),
           gpp = (x*gpp.rate)*unique(meanx),
           wt = y*unique(meany),
           gd = (y*meany - y.init*meany)/t,
           gdc = gd*0.5*1000,
           gppc = x) %>% 
    return()
}

#====set aesthetics====
theme_set(theme_bw()+
            theme(panel.grid = element_blank(),
                  strip.background = element_blank(),
                  legend.position = "bottom",
                  text = element_text(size = 10),
                  axis.title = element_text(size = 10.5)))


midge_color <- scale_color_manual(values = c("Midges" = "black", "No Midges" = "gray60"))
midge_fill <- scale_fill_manual(values = c("Midges" = "black", "No Midges" = "gray60"))
midge_lines <- scale_linetype_manual(values = c("Midges" = "solid", "No Midges" = "dashed"))

