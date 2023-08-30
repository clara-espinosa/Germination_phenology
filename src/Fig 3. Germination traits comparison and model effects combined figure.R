library(FD); library(vegan); library(FactoMineR); library(emmeans);
library(tidyverse); library(ggrepel); library(cowplot);library(ggpubr);
library (binom);library (ggsignif);library (rstatix); library (stringr);
library(plyr);library(patchwork)
Sys.setlocale("LC_ALL","English")
##### effect size fig 3 A #####
library(viridis)
x11()
read.csv("data/test_effectsize.csv", sep=";") %>%
  convert_as_factor(Trait, community, terms) %>%
  mutate (Trait = fct_relevel(Trait, "Environmental heat sum", "t50",
                              "Germination in winter conditions", "Summer germination", 
                              "Spring germination","Autumn germination", "Total germination" )) -> effect_size

ggplot(effect_size, aes(x= Trait, y =effect_size, ymin = L95, ymax = U95, color = Trait))+
  geom_point( size = 3) +
  geom_errorbar (width = 0.2, size =1.2) +
  facet_wrap (~ community, ncol = 2, nrow =7, scales = "free_x") +
  scale_color_manual (values = c("#AC1926", "#891171", "#33407D", "#077395", "#00BC7F", "#AADB41", "#FDE333")) +
  scale_x_discrete(labels = function(Trait) str_wrap(Trait, width = 13)) +
  scale_y_continuous (limits = c(-5,5), breaks = seq (-4, 4, by = 2)) +
  geom_hline(yintercept = 0, linetype = "dashed", size =1, color = "black") +
  coord_flip() +
  labs(title = "A", y = "Effect size") + 
  theme_classic(base_size = 14) +
  theme(plot.title = element_text (size = 22),
        strip.text = element_text(face = "bold", size = 22),
        legend.position = "none",
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.x = element_text (size=18)) +
  annotate("text", x= 0.70, y = -2.5, label = "Higher in Fellfield", size = 3.5, color = "chocolate2" , fontface ="bold" ) +
  annotate ("segment", x = 0.5, xend = 0.5, y= -0.1, yend=-4.9, colour= "chocolate2", size = 2, arrow = arrow()) +
  annotate("text", x= 0.70, y = 2.6, label = "Higher in Snowbed", size = 3.5, colour= "deepskyblue3", fontface ="bold" ) +
  annotate ("segment", x = 0.5, xend = 0.5, y= 0.1, yend= 4.9, colour= "deepskyblue3", size = 2, arrow = arrow()) ->effect_size_graph

ggsave(filename = "results/Fig3Av4.png", effect_size_graph, path = NULL, 
       scale = 1, width = 360, height = 360, units = "mm", dpi = 600)
#### mean value graph fig 3 B  ####
#Total germination  
detach(package:plyr)
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  group_by(species, code, incubator, petridish) %>%
  summarise(total_germ = sum(germinated)) %>% # 
  merge(viables) %>%
  mutate(germPER = (total_germ/viable) *100, # 
         germPER = round (germPER, digit =2))%>%
  filter(viablePER>25)%>%
  filter(germPER>0)%>%
  merge(species) %>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>% 
  group_by (incubator, community) %>%
  summarise (total_germ = sum(total_germ),
             viable = sum (viable))%>%
  mutate (binom.confint(total_germ, viable, methods = "wilson")) -> totalgerm_graph
# AUTUMN (test representation) 
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  group_by (species, code, incubator, petridish, time)%>%
  summarise(seeds_germ = sum(germinated)) %>%
  filter (between(time, 1, 105)) %>%## 105 = 12/11 last check before winter
  group_by (species, code, incubator, petridish) %>%
  summarise(seeds_germ = sum(seeds_germ))%>%
  merge(viables)%>%
  filter(viablePER>25)%>%
  filter(germPER>0)%>% 
  merge(species, by = c("code", "species")) %>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>% 
  select (code, species, community, incubator, petridish, seeds_germ, viable)%>%
  group_by (incubator, community)%>%
  summarise (seeds_germ = sum (seeds_germ),
             viable = sum(viable)) %>%
  mutate (binom.confint(seeds_germ, viable, methods = "wilson")) %>% 
  mutate (Season = "Autumn") -> autumn_graph
# SPRING (test representation) 
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  spread(date, germinated, fill = 0) %>% # wide format for dates, and fill Na with 0
  gather ("date", "germinated", 8: last_col() )%>% # back in long format frrom 8th colum to the last
  arrange (species, accession, code, incubator, petridish, date)%>% # sort row observations this way
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  group_by (species, code, incubator, petridish, time)%>% 
  summarise(seeds_germ = sum(germinated)) %>%
  mutate(cumulative = cumsum(seeds_germ)) %>%
  filter (between(time, 106, 322)) %>% # amount of days calculated from the dates 106 = 13/11 and 322 = 16/6
  ungroup () %>%
  group_by (species, code, incubator, petridish) %>%
  summarise (seeds_germ = sum(seeds_germ),
             cumulative = last(cumulative)) %>%
  merge(viables)%>%
  filter(viablePER>25)%>%
  filter(germPER>0)%>%
  merge(species) %>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>% 
  group_by (incubator, community)%>%
  summarise (seeds_germ = sum (seeds_germ),
             viable = sum(viable)) %>%
  mutate (binom.confint(seeds_germ, viable, methods = "wilson")) %>% 
  mutate (Season = "Spring")-> spring_graph  
# END-SUMMER (test representation) 
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))  %>%
  spread(date, germinated, fill = 0) %>% # wide format for dates, and fill Na with 0
  gather ("date", "germinated", 8: last_col() )%>% # back in long format frrom 8th colum to the last
  arrange (species, accession, code, incubator, petridish, date)%>% # sort row observations this way
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  group_by (species, code, incubator, petridish, time)%>% 
  summarise(seeds_germ = sum(germinated)) %>%
  mutate(cumulative = cumsum(seeds_germ)) %>%
  filter (between(time, 323, 430)) %>% # amount of days calculated from the dates 323 = 17/6 and 430 = 19/09
  ungroup () %>%
  group_by (species, code, incubator, petridish) %>%
  summarise (seeds_germ = sum(seeds_germ),
             cumulative = last(cumulative)) %>%
  merge(viables)%>%
  filter(viablePER>25)%>%
  filter(germPER>0)%>%
  select(species, code, incubator, petridish, seeds_germ, cumulative, viable) %>%
  merge(species, by = c("code", "species")) %>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>% 
  group_by (incubator, community)%>%
  summarise (seeds_germ = sum (seeds_germ),
             viable = sum(viable)) %>%
  mutate (binom.confint(seeds_germ, viable, methods = "wilson")) %>% 
  mutate (Season = "Summer")-> summer_graph  
# GERMINATION WINTER CONDITIONS: NON-GROWING season germination <3ºC
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))  %>%
  spread(date, germinated, fill = 0) %>% # wide format for dates, and fill Na with 0
  gather ("date", "germinated", 8: last_col()) %>% # back in long format from 8th colum to the last
  arrange (species, accession, code, incubator, petridish, date)%>% # sort row observations this way
  group_by (species, code, incubator, petridish, date)%>% 
  summarise(seeds_germ = sum(germinated)) %>%
  mutate(cumulative = cumsum(seeds_germ)) %>%
  filter (date %in% c("2022-04-04", # # first check after winter before Tmin 2ºC
                      "2022-05-11","2022-05-12", "2022-05-18", "2022-05-19", 
                      "2022-05-25", "2022-05-26")) %>% #undersnow + first checks Tmin 2ºC
  ungroup () %>%
  group_by (species, code, incubator, petridish) %>%
  summarise (seeds_germ = sum(seeds_germ)) %>%
  merge (viables) %>%
  filter(viablePER>25)%>%
  filter(germPER>0)%>%
  select(species, code, incubator, petridish, seeds_germ, viable)%>%
  arrange(species) %>%
  merge(species, by = c("code", "species")) %>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>% 
  group_by (community, incubator)%>%
  summarise (seeds_germ = sum (seeds_germ),
             viable = sum(viable)) %>%
  mutate (binom.confint(seeds_germ, viable, methods = "wilson")) %>% 
  mutate (Temperature = "< 2ºC")-> nogrow_graph 

# T50 
t50model %>%
  merge(viables) %>%
  filter(viablePER>25)%>%
  filter(germPER>0)%>%
  merge(species, by = c("code", "species"))-> t50_graph

# Heat sum 
t50_dates%>%
  group_by (species, code, incubator, petridish) %>%
  do (heatsum(.)) %>%
  merge(species, by = c("code", "species"))-> heatsum_graph

### individual 

library("plyr")
x11()
read.csv("data/meanvalues_graph.csv", sep = ";")%>%
  convert_as_factor(trait, community, incubator) %>%
  mutate (trait = fct_relevel(trait, "Total germination", "Autumn germination","Spring germination", "Summer germination", 
                              "Germination in winter conditions",  "t50","Environmental heat sum")) -> mean_values
ggplot(mean_values, aes(incubator, mean, fill=incubator))+
  geom_col(position = position_dodge(0.7), width = 0.75) +
  facet_grid (trait ~ community, scales = "free_y", labeller = label_wrap_gen (width = 13)) +
  geom_errorbar(aes(incubator, mean, ymin = lower, ymax = upper), color = "black",width = 0.2, size =1) +
  scale_fill_manual (name= "Incubator", values = c ("Fellfield"= "chocolate2", "Snowbed" ="deepskyblue3")) +
  labs(title = "B", x = "Incubator", y = "Degres (ºC)         Time (days)                                                                           Germination proportion                                ") + 
  theme_classic(base_size = 14) +
  theme(plot.title = element_text (size = 22),
        strip.text.x = element_text(face = "bold", size = 22),
        strip.text.y = element_text(size = 14, angle = 360),
        legend.position = "none",
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title.y = element_text (size=16), 
        axis.title.x = element_text (size=16)) -> mean_values_graph

ggsave(filename = "results/Fig3Bv3.png", mean_values_graph, path = NULL, 
       scale = 1, width = 240, height = 360, units = "mm", dpi = 600)

#combine both graphs 
Fig3 <- effect_size_graph + mean_values_graph

ggsave(filename = "results/Fig3 no sig.png", Fig3, path = NULL, 
       scale = 1, width = 400, height = 360, units = "mm", dpi = 600)


