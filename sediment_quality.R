library(tidyverse)

#c molar mass
cmm <- 12.011 * #g/mol
  1e6 #microgram/g
#n molar mass
nmm <- 14.0067 * #g/mol
  1e6 #micrgram/g


sed_iso <- read_csv("raw data/sediment_layeranalysis_2019.csv") %>% 
  mutate(c = mass_c_ug/cmm,
         n = mass_n_ug/nmm,
         cn = c/n,
         layer = factor(depth_cm, unique(.$depth_cm), labels = c("top", "mid", "bot")),
         site_depth = factor(ifelse(sta == "E5", 2.5, 4.3))) 



sed_iso %>% 
  ggplot(aes(x = depth_cm, y = cn, col = sta, group= interaction(sta, core)))+
  geom_point()+
  geom_line()



lm(cn~layer*sta, data = sed_iso) %>% 
car::Anova()


 
  
lm(cn~layer, data = sed_iso) %>% car::Anova()



boxplot(cn~layer, data = sed_iso)


lm(cn~depth_cm*sta, data = sed_iso) %>% car::Anova()


kruskal.test(cn~layer, data =sed_iso) 


sed_iso %>% 
  ggplot(aes(x = depth_cm, y = d15N, col = sta, group= interaction(sta, core)))+
  geom_point()+
  geom_line()


sed_iso %>% 
  gather(var, val, contains("d1")) %>% 
  ggplot(aes(x = depth_cm, y = val, col = site_depth, group= interaction(sta, core)))+
  facet_wrap(~var, scales = "free")+
  geom_point(size = 2)+
  geom_line(size = 1)+
  labs(x = "Depth of Layer in Sediment (cm)",
       color = "Site Depth (m)",
       y = element_blank()) +
  scale_color_manual(values = c("dodgerblue", "black"))+
  # scale_color_viridis_d(direction = -1)+
  theme_bw()


sed_iso %>% 
  mutate(site_depth = factor(ifelse(sta == "E5", 2.5, 4.3))) %>% 
  ggplot(aes(x = depth_cm, y = cn, col = site_depth, group= interaction(sta, core)))+
  geom_point(size = 2)+
  geom_line(size = 1)+
  scale_color_manual(values = c("dodgerblue", "black"))+
  labs(x = "Depth of Layer in Sediment (cm)",
       color = "Site Depth (m)",
       y = "C:N") +
  # scale_color_viridis_d(direction = -1)+
  theme_bw()
