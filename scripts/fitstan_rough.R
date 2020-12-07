library(tidyverse)
library(rstan)

meanna <- function(x){
  if(all(is.na(x))){
    NA
  }
  else{
    mean(x, na.rm = TRUE)
  }
}


cm_raw <- read_csv("clean data/MG_cm.csv")

meta <- read_csv("clean data/MG_meta.csv") %>% 
  mutate(algae_conc2 = ifelse(algae_conc == 0, 2e-3, algae_conc)) #check


cm <- cm_raw %>% 
  left_join(meta) %>% 
  mutate(coreid = as.character(coreid),
         box = as.character(box),
         body_size = body_size/1000) #convert from um to mm

avgd <- cm %>% 
  filter(species_name == "tt") %>% 
  group_by(coreid, day) %>% 
  summarise(n = n(),
            l = meanna(body_size))

size_init <- avgd %>% 
  filter(day == 0)

size14 <- avgd %>% 
  filter(day == 14,
         !is.na(l)) %>% 
  bind_cols(avgd %>% 
              ungroup %>% 
              filter(day == 0) %>% 
              summarise(n0 =mean(n)) %>% 
              select(n0))

N <- data.frame(n = size14$n, n0 = size14$n0)
  
data_list <- list(nobs = nrow(N),
     initobs = 3,
     l = size14$l,
     linit = size_init$l,
     wla_mean = 0.0442,
     wla_se = 0.004, #not in lindegaard
     wlb_mean = 0.0879,
     wlb_se = 0.008, #not in lindegaard
     N = N,
     t = 14)

model_path = "scripts/midgeproductivity.stan"
chains = 4
iter = 8000

P14_fit <- stan(file = model_path, data = data_list, seed=1, chains = chains, iter = iter, cores = 1, verbose = TRUE)
