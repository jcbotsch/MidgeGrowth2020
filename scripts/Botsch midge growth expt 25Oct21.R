trajectory <- function(z.init, r, K, a, c, m, Tmax){

	for(i in 1:nrow(z.init)){
		x <- z.init[i,1]
		y <- z.init[i,2]
	
		X <- matrix(NA, 2, Tmax)
		for(t in 1:Tmax){
			#x.next <- exp(r*(1 - x/K))*x - a*x*y
			x.next <- exp(r)*x^K - a*x*y
			y <- y + c*a*x*y - m
			x <- x.next
			y <- max(y,0)

			X[,t] <- matrix(c(x,y),2,1)
		}	
		if(i == 1) maxX <- max(X)
		plot(X[1,], typ="l", ylim=c(0,maxX))
		lines(X[2,], col="red")
	}
}


growth <- function(z.init, r, K, a, c, m, Tmax){
	
	X <- matrix(NA, nrow(z.init), 2)
	for(i in 1:nrow(z.init)){
		x <- z.init[i,1]
		y <- z.init[i,2]
		for(t in 1:Tmax){
			#x.next <- exp(r*(1 - x/K))*x - a*x*y
			x.next <- exp(r)*x^K - a*x*y
			y <- y + c*a*x*y - m
			x <- x.next
		}
		X[i,] <- c(x,y)
	}	
	
	return(cbind(x.growth=exp(r)*X[,1]^K, 
	             delta.y=X[,2] - z.init[,2]))
}

#############################################################
# Calculate the midge growth as a function of their consumption rate

x.init = .5
y.init = 1
r = .2
K = 0.9
q = 0
c = 0.3
m = .0
Tmax = 22

a.range <- .02*(0:20)
x.init.range <- c(.5, 4)
df <- data.frame(a = rep(a.range, times=length(x.init.range)), x.init=rep(x.init.range, each=length(a.range)))
for(x.init in x.init.range) for(a in a.range){
	X <- growth(z.init = cbind(x.init, y.init), r = r, K = K, a = a, c = c, m = m, Tmax = Tmax)
	df$x.growth[df$a == a & df$x.init == x.init] <- X[1]
	df$delta.y[df$a == a & df$x.init == x.init] <- X[2]
}

par(mfrow=c(1,1))
plot(delta.y ~ a, data=df[df$x.init==x.init.range[2],], typ="l")
lines(delta.y ~ a, data=df[df$x.init==x.init.range[1],], col="red")

#############################################################
# Look at data
read.csv("data.csv")
read.csv("data2.csv")

data <- read.csv("microcosm_biomass.csv")
names(data)

# midge weights at time 0
init.midges <- data[is.na(data$algae_conc2),]
y.init.mean <- sum(init.midges$live_tt * init.midges$wt)/sum(init.midges$live_tt)
c(min(init.midges$wt), y.init.mean,max(init.midges$wt))

# midge weights at time 14
midge.wt <- data$wt[data$day == 14 & data$midge == "Midges"]
c(min(midge.wt, na.rm=T), mean(midge.wt, na.rm=T), max(midge.wt, na.rm=T))

# midge weights at time 22
midge.wt <- data$wt[data$day == 22 & data$midge == "Midges"]
c(min(midge.wt, na.rm=T), mean(midge.wt, na.rm=T), max(midge.wt, na.rm=T))

# prune to experimental falcon tubes
data <- data[!is.na(data$algae_conc2),]
data <- data[order(data$day),]
data <- data[order(data$coreid),]
data$midge.trt <- data$midge == "Midges"

################################
# algae
plot(chl ~ algae_conc2, data=data, ylim=c(0,1.2*10^4))
lines(c(0, 1), c(0, 1.2*10^4))

par(mfrow=c(1,2))
plot(gpp ~ chl, data=data[data$day == 14,], col=midge.trt+1)
plot(gpp ~ chl, data=data[data$day == 22,], col=midge.trt+1)

summary(lm(gpp ~ chl + day*midge.trt, data=data))

# algal trajectories
gpp.rate <- 10^-6

par(mfrow=c(2,3))
plot(c(0,22),c(0.0001, .06), typ="n", log="y")
for(i.core in unique(data$coreid)){
	d <- data[data$coreid == i.core,]
	df <- data.frame(day=c(0,d$day), gpp=c(gpp.rate*d$chl[nrow(d)], d$gpp))
	lines(gpp ~ day, data = df, col=d$midge.trt[1]+1)
}
plot(c(0,22),c(0.0001, .06), typ="n", log="y")
for(i.core in unique(data$coreid)){
	d <- data[data$coreid == i.core,]
	df <- data.frame(day=c(0,d$day), gpp=c(gpp.rate*d$chl[nrow(d)], d$gpp))
	if(d$midge.trt[1]) lines(gpp ~ day, data = df, col=d$midge.trt[1]+1)
}
plot(c(0,22),c(0.0001, .06), typ="n", log="y")
for(i.core in unique(data$coreid)){
	d <- data[data$coreid == i.core,]
	df <- data.frame(day=c(0,d$day), gpp=c(gpp.rate*d$chl[nrow(d)], d$gpp))
	if(!d$midge.trt[1]) lines(gpp ~ day, data = df, col=d$midge.trt[1]+1)
}

gpp.rate <- 2*10^-2

plot(c(0,22),c(0.0001, .06), typ="n", log="y")
for(i.core in unique(data$coreid)){
	d <- data[data$coreid == i.core,]
	df <- data.frame(day=c(0,d$day), gpp=c(gpp.rate*d$algae_conc2[nrow(d)], d$gpp))
	lines(gpp ~ day, data = df, col=d$midge.trt[1]+1)
}
plot(c(0,22),c(0.0001, .06), typ="n", log="y")
for(i.core in unique(data$coreid)){
	d <- data[data$coreid == i.core,]
	df <- data.frame(day=c(0,d$day), gpp=c(gpp.rate*d$algae_conc2[nrow(d)], d$gpp))
	if(d$midge.trt[1]) lines(gpp ~ day, data = df, col=d$midge.trt[1]+1)
}
plot(c(0,22),c(0.0001, .06), typ="n", log="y")
for(i.core in unique(data$coreid)){
	d <- data[data$coreid == i.core,]
	df <- data.frame(day=c(0,d$day), gpp=c(gpp.rate*d$algae_conc2[nrow(d)], d$gpp))
	if(!d$midge.trt[1]) lines(gpp ~ day, data = df, col=d$midge.trt[1]+1)
}

################################
# midges
par(mfrow=c(1,2))
d <- data[data$midge.trt,]
plot(wt ~ chl, data=d[d$day == 14,], ylim =c(0,.06))
plot(wt ~ chl, data=d[d$day == 22,], ylim =c(0,.06))

#############################################################
# fit model
#============================================================
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
		#x.next <- exp(r*(1 - x/K))*x - a*x*y
		x.next <- exp(r)*x^K - a*x*y
		y <- y + ac*x*y - m
		x <- x.next
		
		# if(!tofit) show(c(mean(x), mean(y)))
				
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

#============================================================
# midge.trt.row <- data$coreid[data$midge == "Midges" & data$day==14]

# standardize variables
data$x <- data$gpp/mean(data$gpp, na.rm=T)
data$x0 <- data$chl/mean(data$chl, na.rm=T)
data$y <- data$wt/mean(data$wt, na.rm=T)
data$y0 <- 0
data$y0[data$midge.trt] <- y.init.mean/mean(data$wt, na.rm=T)

# fit model
r <- .5
K <- .9
a <- .5
ac <- .03
m <- .001
gpp.rate <- 10^0

x.init <- data$x0[!is.na(data$x0)]
y.init <- data$y0[!is.na(data$x0)]
# cbind(x.init, y.init, data$midge[!is.na(data$x0)])

SS.weight <- 30

par.full <- log(c(r, K, a, ac, m, gpp.rate))
par.fixed <- c(NA, NA, NA, NA, NA, NA)
#par.fixed <- c(NA, log(K), NA, NA, NA, NA)
par <- par.full[is.na(par.fixed)]

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
round(c(r=r, K=K, a=a, c=c, m=m, gpp.rate=gpp.rate, SS=opt.opt$value), digits=3)
       # r        K        a        c        m gpp.rate       SS 
   # 0.187    0.927    0.045    0.163    0.002    0.183   66.817 
 
   
# figures

x.sim <- unique(x.init[order(x.init, decreasing=T)])
y.sim <- unique(y.init[order(y.init)])
z.sim <- cbind(rep(x.sim, times=2), rep(y.sim, each=length(x.sim)))
Tmax <- 22

# trajectories
par(mfrow=c(4,5), mai=c(.1,.1,.1,.1))
trajectory(z.init = z.sim, r = r, K = K, a = a, c = c, m = m, Tmax = Tmax)

# Calculate the midge growth as a function of their consumption rate
a.range <- a*(0:40)/10
x.sim <- max(x.init)
# x.init.range <- c(min(x.init),mean(x.init),max(x.init))
y.init.sim <- y.init.mean/mean(data$wt, na.rm=T)
y.init.range <- c(1,3,5)*y.init.sim

df <- data.frame(a = rep(a.range, times=length(y.init.range)), y.init=rep(y.init.range, each=length(a.range)))
for(y.sim in y.init.range) for(a.sim in a.range){
	X <- growth(z.init = cbind(x.sim, y.sim), r = r, K = K, a = a.sim, c = c, m = m, Tmax = Tmax)
	df$x.growth[df$a == a.sim & df$y.init == y.sim] <- X[1]
	df$delta.y[df$a == a.sim & df$y.init == y.sim] <- X[2]
}

par(mfrow=c(1,1), mai=c(.8,.8,.1,.1))
plot(delta.y ~ a, data=df[df$y.init==y.init.range[1],], typ="l", ylim=c(0, max(delta.y)))
lines(delta.y ~ a, data=df[df$y.init==y.init.range[2],], col="red")
lines(delta.y ~ a, data=df[df$y.init==y.init.range[3],], col="blue")
lines(c(a,a),c(0,10), col="orange")
lines(c(0,10),c(0,0), lty=2)
