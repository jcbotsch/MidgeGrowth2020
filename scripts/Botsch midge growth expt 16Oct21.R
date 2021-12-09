trajectory <- function(x.init, y.init, r, K, q, a, c, m, Tmax){
	
	x <- x.init
	y <- y.init
	
	X <- matrix(NA, 2, Tmax)
	for(t in 1:Tmax){
		x.next <- exp(r*(1 - x/K))*x - (1 - q)*a*max(0, x)*y^(2/3)
		y <- y + c*(1 - q)*a*max(0, x)*y^(2/3) - m
		x <- x.next
		
		X[,t] <- matrix(c(x,y),2,1)
	}	
	plot(X[1,], typ="l", ylim=c(0,max(X)))
	lines(X[2,], col="red")
}
trajectory <- Vectorize(trajectory)

growth <- function(x.init, y.init, r, K, q, a, c, m, Tmax){
	
	x <- x.init
	y <- y.init
	
	for(t in 1:Tmax){
		x.next <- exp(r*(1 - x/K))*x - (1 - q)*a*max(0, x)*y^(2/3)
		# x.next <- exp(r*(1 - x/K))*x - (1 - q)*a*max(0, x)*y
		y <- y + c*(1 - q)*a*max(0, x)*y^(2/3) - m
		x <- x.next
	}	
	
	return(c(x.growth=exp(r*(1 - x/K))*x, delta.y=y - y.init))
}
growth <- Vectorize(growth)

#############################################################
# Calculate the midge growth as a function of their consumption rate

x.init = .5
y.init = 1
r = .2
K = 20
q = 0
c = 0.3
m = .0
Tmax = 22

a.range <- .02*(0:20)
x.init.range <- c(.5, 4)
df <- data.frame(a = rep(a.range, times=length(x.init.range)), x.init=rep(x.init.range, each=length(a.range)))
for(x.init in x.init.range) for(a in a.range){
	X <- growth(x.init = x.init, y.init = y.init, r = r, K = K, q = q, a = a, c = c, m = m, Tmax = Tmax)
	df$x.growth[df$a == a & df$x.init == x.init] <- X[1,1]
	df$delta.y[df$a == a & df$x.init == x.init] <- X[2,1]
}

par(mfrow=c(1,1))
plot(delta.y ~ a, data=df[df$x.init==x.init.range[2],], typ="l")
lines(delta.y ~ a, data=df[df$x.init==x.init.range[1],], col="red")

############################################################
# Calculate the effect of initial algae concentrations

a = .07
x.init = c(.5*(1:10))

# trajectories
par(mfrow=c(2,6))
trajectory(x.init = x.init, y.init = y.init, r = r, K = K, q = q, a = a, c = c, m = m, Tmax = Tmax)

par(mfrow=c(2,2))

X <- growth(x.init = x.init, y.init = y.init, r = r, K = K, q = q, a = a, c = c, m = m, Tmax = Tmax)

plot(X[1,]/x.init, typ="l", ylim=c(0, max(X)), xlab="init.x")
lines(X[2,], col="red")

plot(X[1,]/x.init, X[2,], xlim=c(0,max(X[1,])), ylim=c(0,max(X[2,])))

plot(X[1,], typ="l", xlab="init.x", ylim=c(0,max(X)))
lines(X[2,], col="red")

plot(t(X), xlim=c(0,max(X[1,])), ylim=c(0,max(X[2,])))
