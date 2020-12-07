
data {
  int nobs; // number of observations
  int initobs; //number of initial observations
  real l[nobs]; //length of midges in a mesocosm
  real linit[initobs]; //starting lengths
  real wla_mean; //mean intercept for length to weight conversion from lindegaard ms
  real wlb_mean; //mean slope for length ot weight conversion from lindegaard ms 
  real wla_se; //standard error for intercept for length to weight conversion from lindegaard ms
  real wlb_se; // standard error for slope for length ot weight conversion from lindegaard ms 
  matrix[nobs, 2] N; //number of midges in a mesocosm at time t and t0 
  real t; // sample event
}

transformed data{
  vector[nobs] Nbar; //average abundance in mesocosm between time t and t-1
  
  Nbar = (N[,1]+N[,2])/2; //calculate Nbar
}


parameters {
  real<lower=0> sdl; //estimated standard error on length
  real wla; //estimated length:weight intercept
  real wlb; // estimated length:weight slope
  vector[nobs] Pm; //productivity of midges
  real<lower=0> sdlinit; //estimated standard error 
  real<lower=0> linit_mean; //average initial length of a midge

}

transformed parameters {
    vector[nobs] w; // weight (mg)
    vector[nobs] lpred; //calculated mean length of a midge
    real winit; //average initial weight of a midge
    
    
    winit = (linit_mean*wlb+wla)^(3); //calculate average initial weight of a midge
    for (i in 1:nobs){
          w[i] = (Pm[i] * t)/Nbar[i] + winit; //calculate weight from production (increment summation)
          lpred[i] = (w[i]^(1/3)-wla)/wlb; //estimate weight
    }
}

model {
  sdl ~ exponential(0.5); 
  sdlinit ~ exponential(0.5);
  l ~ normal(lpred, sdl); // length is normally distributed around the mean and sd of length
  wla ~ normal(wla_mean, wla_se); //distribution for w:l intercept
  wlb ~ normal(wlb_mean, wlb_se); //distribution for w:l slope
  Pm ~ normal(0, 0.1); //implicitly allows Pm to be negative
  linit_mean ~ exponential(0.5); //mm 
  linit ~ normal(linit_mean, sdlinit); //initial body lengths 
}

generated quantities {
  //for calculations later...
}
