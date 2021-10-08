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
		y <- y + c*(1 - q)*a*max(0, x)*y^(2/3) - m
		x <- x.next
	}	
	
	return(c(x.growth=exp(r*(1 - x/K))*x, delta.y=y - y.init))
}
growth <- Vectorize(growth)
#############################################################

x.init = 1
y.init = .5
r = .1
K = 200
q = 0
a = .04
c = .3
m = .01
Tmax = 20

# trajectories
par(mfrow=c(2,5))
trajectory(x.init = 1:10, y.init = y.init, r = r, K = K, q = q, a = a, c = c, m = m, Tmax = Tmax)


X <- growth(x.init = 1:10, y.init = y.init, r = r, K = K, q = q, a = a, c = c, m = m, Tmax = Tmax)

par(mfrow=c(1,2))
plot(X[1,], typ="l", ylim=c(0, max(X)), xlab="init.x")
lines(X[2,], col="red")

plot(t(X), xlim=c(0,max(X[1,])), ylim=c(0,max(X[2,])))