library(FD); library(vegan); library(FactoMineR); library(emmeans);
library(tidyverse); library(ggrepel); library(cowplot);
library (binom);library (ggsignif);library (rstatix);
theme_set(theme_cowplot(font_size = 10)) 

## GERMINATION EXPERIMENTS
#### dataframe  AND VISUALIZATION for temperature programs ####
temp <- read.csv("data/date_temp.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y")) %>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  mutate(incubator = as.factor(incubator))
temp <- as.data.frame(temp)
str(temp)
summary(temp)
# visualization 1
# mulitple lines with shaded areas in between
x11()
ggplot(temp, aes ()) +
  geom_line (aes (x=time, y=Tmin, colour =incubator), size =1.25) +
  geom_line (aes (x=time, y=Tmax, colour = incubator), size =1.25) +
  geom_ribbon (data = temp, aes (x=time, ymin =Tmin, ymax=Tmax, colour = incubator, fill = incubator), alpha =0.5) +
  scale_fill_manual (values =c("chocolate2", "deepskyblue3")) +
  scale_color_manual (values =c("chocolate2", "deepskyblue3")) +
  labs (title = "Move-along temperature regimes", y= "Temperature ºC", x = "Time (days)") + 
  theme (axis.title.y = element_text (size=14), 
         axis.title.x = element_text (size=14), 
         plot.title = element_text (size =16), 
         legend.title = element_text(size = 14),
         legend.text = element_text (size =12)) +
  geom_vline(xintercept = 122, linetype = "dashed", size =1.25) +
  annotate (geom ="text", x= 137, y = 22, label ="Winter begins", colour = "black", size = 4,  fontface ="bold") +
  geom_vline(xintercept = 240, linetype = "dashed", size =1.25, color = "chocolate2") +
  annotate (geom ="text", x= 255, y = 22, label ="Fellfield\nwinter ends", colour = "chocolate2", size = 4,  fontface ="bold") +
  geom_vline(xintercept = 291, linetype = "dashed", size =1.25,color = "deepskyblue3") +
  annotate (geom ="text", x= 307, y = 22, label ="Snowbed\nwinter ends", colour = "deepskyblue3", size = 4,  fontface ="bold") +
  geom_hline(yintercept=0, linetype ="dashed", size =1, colour = "red")
  
########### GERMINATION RATE ###################################################
#### germination peaks  #####
data.frame(macroclimate = c("Mediterranean","Mediterranean","Mediterranean","Temperate", "Temperate"),
           incubator  = c("Fellfield", "Fellfield","Snowbed", "Fellfield", "Fellfield"),
           time = c(105, 247, 105, 105, 247),
           germinated = c(0,0,0,0,0)) -> data_vis
x11()
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date)))) %>%
  merge (species) %>%
  group_by(macroclimate, incubator, time) %>%
  summarise(germinated = sum(germinated)) %>%
  bind_rows (data_vis)%>%
  merge(viables_macroclimate) %>%
  mutate(germination = germinated/viable) %>%
  ggplot(aes(time, germination, color = incubator, fill = incubator)) +
  geom_line(size = 2) +
  scale_color_manual (name= "Incubator", values = c ("Fellfield"= "chocolate2", "Snowbed" ="deepskyblue3")) +
  facet_wrap(~ macroclimate, scales = "free_x", ncol = 2) +
  labs(x = "Time (days)", y = "Germination proportion") +
  theme(strip.text = element_text(face = "italic", size = 26), 
        legend.title = element_text(size=20), 
        legend.text = element_text(size=16),
        axis.title.y = element_text (size=20), 
        axis.title.x = element_text (size=20)) + 
  geom_vline(xintercept = 122, linetype = "dashed", size= 2) +
  geom_vline(xintercept = 240, linetype = "dashed", size= 2, color = "chocolate2") +
  geom_vline(xintercept = 291, linetype = "dashed", size= 2, color = "deepskyblue3")


#### cumulative germination  ####
# tidyverse modification to have the accumulated germination along the whole experiment + ggplot
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date)))) %>%
  merge (species) %>%
  group_by(macroclimate, incubator, time) %>%
  bind_rows (data_vis)%>%
  summarise(germinated = sum(germinated)) %>%
  mutate(germinated = cumsum(germinated)) %>%
  merge(viables_macroclimate) %>%
  mutate(germination = germinated/viable) %>%
  ggplot(aes(time, germination, color = incubator, fill = incubator)) +
  geom_line(size = 2) +
  scale_color_manual (name= "Incubator", values = c ("Fellfield"= "chocolate2", "Snowbed" ="deepskyblue3")) +
  facet_wrap(~ macroclimate, scales = "free_x", ncol = 2) +
  labs(title = "", x = "Time (days)", y = "Germination proportion") +
  theme(strip.text = element_text(face = "italic", size = 26), 
        legend.title = element_text(size=20), 
        legend.text = element_text(size=16),
        axis.title.y = element_text (size=20), 
        axis.title.x = element_text (size=20)) +
  geom_vline(xintercept = 122, linetype = "dashed", size= 2) +
  geom_vline(xintercept = 240, linetype = "dashed", size= 2, color = "chocolate2") +
  geom_vline(xintercept = 291, linetype = "dashed", size= 2, color = "deepskyblue3") 

########### TOTAL GERMINATIONv + TIMING  #######################################
# Total germination (test representation) ####
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  group_by(species, code, incubator, petridish) %>%
  summarise(total_germ = sum(germinated)) %>% # 
  merge(viables) %>%
  mutate(germPER = (total_germ/viable) *100, # 
         germPER = round (germPER, digit =2))%>%
  filter(viablePER>25)%>%
  filter(germPER>0)%>%
  mutate(species= str_replace(species, "Cerastium sp", "Cerastium pumilum"))%>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>% 
  merge(species) %>%  
  group_by (macroclimate, incubator) %>%
  summarise (total_germ = sum(total_germ),
             viable = sum (viable))%>%
  mutate (binom.confint(total_germ, viable, methods = "wilson")) -> totalgerm_graph

x11()
ggplot(data=totalgerm_graph, aes(incubator, mean, color=macroclimate))+
  geom_point(size =2) +
  #facet_grid (.~macroclimate) +
  geom_signif(comparisons = list(c("Fellfield", "Snowbed")),annotations = "***", y_position = 0.95, tip_length = 0.05, color = "black", size = 1, textsize = 9) +
  #geom_signif(comparisons = list(c("Mediterranean", "Temperate")),annotations = "***", y_position = 0.85, tip_length = 0.2, color = "black", size = 1, textsize = 9) +
  geom_errorbar(data= totalgerm_graph, aes(incubator, mean, ymin = lower, ymax = upper,
                                      color = macroclimate),width = 0.3, size =1.2) +
  geom_segment (aes(x=1.5, y =0.72, xend=1.5, yend =0.82), color = "black", size = 1)+
  geom_segment (aes(x=1.45, y = 0.72, xend =1.55, yend =0.72), color = "black", size = 1) +
  geom_segment (aes(x=1.45, y = 0.82, xend =1.55, yend =0.82), color = "black", size = 1) +
  annotate (geom ="text", x= 1.60, y = 0.79, label =".", colour = "black", size = 9) +
  #scale_color_manual (name= "Incubator", values = c ("Fellfield"= "chocolate2", "Snowbed" ="deepskyblue3")) +
  scale_color_manual (name = "Macroclimate", values = c ("Mediterranean"= "darkgoldenrod1" , "Temperate" = "forestgreen")) +
  scale_y_continuous (limits = c(0,1), breaks = seq (0, 1, by= 0.2)) +
  labs(title = "Total germination", x = "Incubator", y = "Germination proportion") +
  theme_classic(base_size = 16)+
  theme(plot.title = element_text (hjust = 0.5, size = 30),
      strip.text = element_text(face = "italic", size = 24), 
      legend.position = "right",
      legend.title = element_text(size=20), 
      legend.text = element_text(size=16),
      legend.background = element_rect(fill="transparent",colour=NA),
      axis.title.y = element_text (size=18), 
      axis.title.x = element_text (size=18)) 

# AUTUMN (test representation) ####
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  group_by (species, code, incubator, petridish, time)%>%
  summarise(seeds_germ = sum(germinated)) %>%
  filter (between(time, 1, 105)) %>%## 105 = 12/11 last check before winter
  group_by (species, code, incubator, petridish) %>%
  summarise(seeds_germ = sum(seeds_germ))%>%
  merge(viables_germ)%>%
  filter(viablePER>25)%>%
  filter(germPER>0)%>%
  mutate(species= str_replace(species, "Cerastium sp", "Cerastium pumilum"))%>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>% 
  merge(species, by = c("code", "species")) %>%
  select (code, species, macroclimate, incubator, petridish, seeds_germ, viable)%>%
  group_by (macroclimate, incubator)%>%
  summarise (seeds_germ = sum (seeds_germ),
             viable = sum(viable)) %>%
  mutate (binom.confint(seeds_germ, viable, methods = "wilson")) %>% 
  mutate (Season = "Autumn") -> autumn_graph

x11()
ggplot(data=autumn_graph, aes(incubator, mean, color=macroclimate))+
  geom_point(size =2) +
  geom_signif(comparisons = list(c("Fellfield", "Snowbed")),annotations = "***", y_position = 0.65, tip_length = 0.04, color = "black", size = 1, textsize = 9) +
  #geom_signif(comparisons = list(c("Mediterranean", "Temperate")),annotations = "***", y_position = 0.65, tip_length = 0.04, color = "black", size = 1, textsize = 9) +
  #facet_grid (.~macroclimate) +
  geom_errorbar(data= autumn_graph, aes(incubator, mean, ymin = lower, ymax = upper, 
                                        color = macroclimate), width = 0.3, size =1.2) +
  geom_segment (aes(x=1.5, y =0.24, xend=1.5, yend =0.4), color = "black", size = 1)+
  geom_segment (aes(x=1.45, y = 0.24, xend =1.55, yend =0.24), color = "black", size = 1) +
  geom_segment (aes(x=1.45, y = 0.4, xend =1.55, yend =0.4), color = "black", size = 1) +
  annotate (geom ="text", x= 1.62, y = 0.3, label ="**", colour = "black", size = 9) +
  #scale_color_manual (name= "Incubator", values = c ("Fellfield"= "chocolate2", "Snowbed" ="deepskyblue3")) +
  scale_color_manual (name = "Macroclimate", values = c ("Mediterranean"= "darkgoldenrod1" , "Temperate" = "forestgreen")) +
  scale_y_continuous (limits = c(0,1), breaks = seq (0, 1, by= 0.2)) +
  labs(title = "Autumn germination", x = "Incubator", y = "Germination proportion") +
  theme_classic(base_size = 16)+
  theme(plot.title = element_text (hjust = 0.5, size = 30),
        strip.text = element_text(face = "italic", size = 24), 
        legend.position = "right",
        legend.title = element_text(size=20), 
        legend.text = element_text(size=16),
        legend.background = element_rect(fill="transparent",colour=NA),
        axis.title.y = element_text (size=18), 
        axis.title.x = element_text (size=18)) 
   
# SPRING (test representation) ####
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
  merge(viables_germ)%>%
  filter(viablePER>25)%>%
  filter(germPER>0)%>%
  mutate(species= str_replace(species, "Cerastium sp", "Cerastium pumilum"))%>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>% 
  merge(species) %>%
  group_by (macroclimate, incubator)%>%
  summarise (seeds_germ = sum (seeds_germ),
             viable = sum(viable)) %>%
  mutate (binom.confint(seeds_germ, viable, methods = "wilson")) %>% 
  mutate (Season = "Spring")-> spring_graph  
x11()
ggplot(data= spring_graph, aes(incubator, mean, color=macroclimate))+
  geom_point(size =2) +
  geom_signif(comparisons = list(c("Fellfield", "Snowbed")), annotations = "***", y_position = 0.5, tip_length = 0.04, color = "black", size = 1, textsize = 9) +
  #facet_grid (.~macroclimate) +
  geom_errorbar(data=  spring_graph, aes(incubator, mean, ymin = lower, ymax = upper, 
                                        color = macroclimate), width = 0.3, size =1.2) +
  #scale_color_manual (name= "Incubator", values = c ("Fellfield"= "chocolate2", "Snowbed" ="deepskyblue3")) +
  scale_color_manual (name = "Macroclimate", values = c ("Mediterranean"= "darkgoldenrod1" , "Temperate" = "forestgreen")) +
  scale_y_continuous (limits = c(0,1), breaks = seq (0, 1, by= 0.2)) +
  labs(title = "Spring germination", x = "Incubator", y = "Germination proportion") +
  theme_classic(base_size = 16)+
  theme(plot.title = element_text (hjust = 0.5, size = 30),
        strip.text = element_text(face = "italic", size = 24), 
        legend.position = "right",
        legend.title = element_text(size=20), 
        legend.text = element_text(size=16),
        legend.background = element_rect(fill="transparent",colour=NA),
        axis.title.y = element_text (size=18), 
        axis.title.x = element_text (size=18))
# END-SUMMER (test representation) ####
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
  merge(viables_germ)%>%
  filter(viablePER>25)%>%
  filter(germPER>0)%>%
  select(species, code, incubator, petridish, seeds_germ, cumulative, viable) %>%
  mutate(species= str_replace(species, "Cerastium sp", "Cerastium pumilum"))%>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>% 
  merge(species, by = c("code", "species")) %>%
  group_by (macroclimate, incubator)%>%
  summarise (seeds_germ = sum (seeds_germ),
             viable = sum(viable)) %>%
  mutate (binom.confint(seeds_germ, viable, methods = "wilson")) %>% 
  mutate (Season = "Summer")-> summer_graph  

x11()
ggplot(data= summer_graph , aes(incubator, mean, color=macroclimate))+
  geom_point(size =2) +
  #geom_signif(comparisons = list(c("Mediterranean", "Temperate")),annotations = "***", y_position = 0.45, tip_length = 0.05, color = "black", size = 1, textsize = 9) +
  geom_signif(comparisons = list(c("Fellfield", "Snowbed")), annotations = "***", y_position = 0.45, tip_length = 0.04, color = "black", size = 1, textsize = 9) +
  #facet_grid (.~macroclimate) +
  geom_errorbar(data=  summer_graph , aes(incubator, mean, ymin = lower, ymax = upper, 
                                         color = macroclimate), width = 0.3, size =1.2) +
  #geom_segment (aes(x=1.5, y =0.08, xend=1.5, yend =0.13), color = "black", size = 1)+
  #geom_segment (aes(x=1.45, y = 0.08, xend =1.55, yend =0.08), color = "black", size = 1) +
  #geom_segment (aes(x=1.45, y = 0.13, xend =1.55, yend =0.13), color = "black", size = 1) +
  #annotate (geom ="text", x= 1.62, y = 0.1, label ="**", colour = "black", size = 9) +
  #scale_color_manual (name= "Incubator", values = c ("Fellfield"= "chocolate2", "Snowbed" ="deepskyblue3")) +
  scale_color_manual (name = "Macroclimate", values = c ("Mediterranean"= "darkgoldenrod1" , "Temperate" = "forestgreen")) +
  scale_y_continuous (limits = c(0,1), breaks = seq (0, 1, by= 0.2)) +
  labs(title = "Summer germination", x = "Incubator", y = "Germination proportion") +
  theme_classic(base_size = 16)+
  theme(plot.title = element_text (hjust = 0.5, size = 30),
        strip.text = element_text(face = "italic", size = 24), 
        legend.position = "right",
        legend.title = element_text(size=20), 
        legend.text = element_text(size=16),
        legend.background = element_rect(fill="transparent",colour=NA),
        axis.title.y = element_text (size=18), 
        axis.title.x = element_text (size=18))

#### TIMING GRAPH ####
rbind(autumn_graph, spring_graph, summer_graph) %>% 
  mutate (Season = factor(Season, levels = c( "Autumn", "Spring", "Summer"))) -> graph #reorder levels
#rm(graph_date)
x11()
ggplot()+
  geom_point(data=graph, aes(Season, mean, color=incubator), size =2) +
  facet_grid (.~macroclimate) +
  geom_errorbar(data= graph, aes(Season, mean, ymin = lower, ymax = upper, color = incubator), width = 0.3, size =1) +
  geom_line(data= graph, aes(as.numeric(Season), mean, color = incubator), size =1) +
  scale_color_manual (name= "Incubator", values = c ("Fellfield"= "chocolate2", "Snowbed" ="deepskyblue3")) +
  scale_y_continuous (limits = c(0,1), breaks = seq (0, 1, by= 0.2)) +
  labs(title = "Germination timing", x = "Season", y = "Germination proportion") +
  theme_classic(base_size = 16)+
  theme(plot.title = element_text (hjust = 0.5, size = 30),
        strip.text = element_text(face = "italic", size = 24), 
        legend.position = "right",
        legend.title = element_text(size=20), 
        legend.text = element_text(size=16),
        legend.background = element_rect(fill="transparent",colour=NA),
        axis.title.y = element_text (size=18), 
        axis.title.x = element_text (size=18))

########### GERMINATION X TEMPERATURES #########################################
# NON-GROWING season germination <2ºC ####
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))  %>%
  spread(date, germinated, fill = 0) %>% # wide format for dates, and fill Na with 0
  gather ("date", "germinated", 8: last_col()) %>% # back in long format frrom 8th colum to the last
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
  merge (viables_germ) %>%
  filter(viablePER>25)%>%
  filter(germPER>0)%>%
  select(species, code, incubator, petridish, seeds_germ, viable)%>%
  arrange(species) %>%
  mutate(species= str_replace(species, "Cerastium sp", "Cerastium pumilum"))%>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>% 
  merge(species, by = c("code", "species")) %>%
  group_by (macroclimate, incubator)%>%
  summarise (seeds_germ = sum (seeds_germ),
             viable = sum(viable)) %>%
  mutate (binom.confint(seeds_germ, viable, methods = "wilson")) %>% 
  mutate (Temperature = "< 2ºC")-> nogrow_graph 

ggplot(data= nogrow_graph, aes(incubator, mean, color=macroclimate))+
  geom_point(size =2) +
  geom_signif(comparisons = list(c("Fellfield", "Snowbed")), annotations = "***", y_position = 0.3, tip_length = 0.04, color = "black", size = 1, textsize = 9) +
  #facet_grid (.~macroclimate) +
  geom_errorbar(data=  nogrow_graph, aes(incubator, mean, ymin = lower, ymax = upper, 
                                         color = macroclimate), width = 0.3, size =1.2) +
  #scale_color_manual (name= "Incubator", values = c ("Fellfield"= "chocolate2", "Snowbed" ="deepskyblue3")) +
  scale_color_manual (name = "Macroclimate", values = c ("Mediterranean"= "darkgoldenrod1" , "Temperate" = "forestgreen")) +
  scale_y_continuous (limits = c(0,1), breaks = seq (0, 1, by= 0.2)) +
  labs(title = "Germination at T<2ºC", x = "Incubator", y = "Germination proportion") +
  theme_classic(base_size = 16)+
  theme(plot.title = element_text (hjust = 0.5, size = 30),
        strip.text = element_text(face = "italic", size = 24), 
        legend.position = "right",
        legend.title = element_text(size=20), 
        legend.text = element_text(size=16),
        legend.background = element_rect(fill="transparent",colour=NA),
        axis.title.y = element_text (size=18), 
        axis.title.x = element_text (size=18))

# COLD GERMINATION ####
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  #filter (species == "Carex sempervirens") %>%
  spread(date, germinated, fill = 0) %>% # wide format for dates, and fill Na with 0
  gather ("date", "germinated", 8: last_col()) %>% # back in long format frrom 8th colum to the last
  arrange (species, accession, code, incubator, petridish, date)%>% # sort row observations this way
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>% 
  group_by (species, code, incubator, petridish, time) %>% 
  summarise(seeds_germ = sum(germinated)) %>%
  mutate(cumulative = cumsum(seeds_germ)) %>%
  filter (incubator == "Fellfield") %>% 
  filter (between(time, 249, 299)) %>% ## amount of days calculated from the dates 249 = 5/4 and 299 = 25/05
  group_by (species, code, incubator, petridish) %>%
  summarise (seeds_germ = sum(seeds_germ),
             cumulative = last(cumulative)) %>%
  merge (germ_AW) %>%
  select(species, code, incubator, petridish, seeds_germ, germ_AW, viable, viablePER, germPER)  -> coldgerm_F

read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))  %>%
  spread(date, germinated, fill = 0) %>% # wide format for dates, and fill Na with 0
  gather ("date", "germinated", 8: last_col()) %>% # back in long format frrom 8th colum to the last
  arrange (species, accession, code, incubator, petridish, date)%>% # sort row observations this way
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>% 
  group_by (species, code, incubator, petridish, time) %>% 
  summarise(seeds_germ = sum(germinated)) %>%
  mutate(cumulative = cumsum(seeds_germ)) %>%
  filter (incubator == "Snowbed") %>% 
  filter (between(time, 301, 322)) %>% ## amount of days calculated from the dates 301 = 27/05 and 322 = 17/06
  group_by (species, code, incubator, petridish) %>%
  summarise (seeds_germ = sum(seeds_germ),
             cumulative = last(cumulative)) %>%
  merge (germ_AW) %>%
  select(species, code, incubator, petridish, seeds_germ, germ_AW, viable, viablePER, germPER) -> coldgerm_S

rbind(coldgerm_F, coldgerm_S) %>%
  arrange(species) %>%
  filter(viablePER>25)%>%
  filter(germPER>0)%>%
  mutate(species= str_replace(species, "Cerastium sp", "Cerastium pumilum"))%>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>% 
  merge(species, by = c("code", "species")) %>%
  group_by (macroclimate, incubator)%>%
  summarise (seeds_germ = sum (seeds_germ),
             viable = sum(viable)) %>%
  mutate (binom.confint(seeds_germ, viable, methods = "wilson")) %>% 
  mutate (Temperature = "2-10ºC")-> coldgerm_graph 

ggplot(data= coldgerm_graph , aes(incubator, mean, color=macroclimate))+
  geom_point(size =2) +
  geom_signif(comparisons = list(c("Fellfield", "Snowbed")), annotations = "***", y_position = 0.4, tip_length = 0.04, color = "black", size = 1, textsize = 9) +
  #facet_grid (.~macroclimate) +
  geom_errorbar(data=  coldgerm_graph , aes(incubator, mean, ymin = lower, ymax = upper, 
                                         color = macroclimate), width = 0.3, size =1.2) +
  #scale_color_manual (name= "Incubator", values = c ("Fellfield"= "chocolate2", "Snowbed" ="deepskyblue3")) +
  scale_color_manual (name = "Macroclimate", values = c ("Mediterranean"= "darkgoldenrod1" , "Temperate" = "forestgreen")) +
  scale_y_continuous (limits = c(0,1), breaks = seq (0, 1, by= 0.2)) +
  labs(title = "Germination at T = 2-10 ºC", x = "Incubator", y = "Germination proportion") +
  theme_classic(base_size = 16)+
  theme(plot.title = element_text (hjust = 0.5, size = 30),
        strip.text = element_text(face = "italic", size = 24), 
        legend.position = "right",
        legend.title = element_text(size=20), 
        legend.text = element_text(size=16),
        legend.background = element_rect(fill="transparent",colour=NA),
        axis.title.y = element_text (size=18), 
        axis.title.x = element_text (size=18))
# WARM GERMINATION #### 
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  spread(date, germinated, fill = 0) %>% # wide format for dates, and fill Na with 0
  gather ("date", "germinated", 8: last_col()) %>% # back in long format frrom 8th colum to the last
  arrange (species, accession, code, incubator, petridish, date)%>% # sort row observations this way
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>% 
  group_by (species, code, incubator, petridish, time) %>% 
  summarise(seeds_germ = sum(germinated)) %>%
  mutate(cumulative = cumsum(seeds_germ)) %>%
  filter (incubator == "Fellfield") %>% 
  filter (between(time, 300, 430)) %>% ## amount of days calculated from the dates 300= 26/05 and 430 = 19/09
  summarise (seeds_germ = sum(seeds_germ),
             cumulative = last(cumulative)) %>%
  merge (germ_AW) %>%
  select(species, code, incubator, petridish, seeds_germ, germ_AW, viable, viablePER, germPER)  -> warmgerm_F

read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  spread(date, germinated, fill = 0) %>% # wide format for dates, and fill Na with 0
  gather ("date", "germinated", 8: last_col()) %>% # back in long format frrom 8th colum to the last
  arrange (species, accession, code, incubator, petridish, date)%>% # sort row observations this way
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>% 
  group_by (species, code, incubator, petridish, time) %>% 
  summarise(seeds_germ = sum(germinated)) %>%
  mutate(cumulative = cumsum(seeds_germ)) %>% 
  filter (incubator == "Snowbed") %>% 
  filter (between(time, 323, 430)) %>% ## amount of days calculated from the dates 323 = 18/06 and 430 = 19/09
  group_by (species, code, incubator, petridish) %>%
  summarise (seeds_germ = sum(seeds_germ),
             cumulative = last(cumulative)) %>%
  merge (germ_AW) %>%
  select(species, code, incubator, petridish, seeds_germ, germ_AW, viable, viablePER, germPER) -> warmgerm_S

rbind(warmgerm_F, warmgerm_S) %>%
  arrange(species) %>%
  filter(viablePER>25)%>%
  filter(germPER>0)%>%
  mutate(species= str_replace(species, "Cerastium sp", "Cerastium pumilum"))%>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>% 
  merge(species, by = c("code", "species")) %>%
  group_by (macroclimate, incubator)%>%
  summarise (seeds_germ = sum (seeds_germ),
             viable = sum(viable)) %>%
  mutate (binom.confint(seeds_germ, viable, methods = "wilson")) %>% 
  mutate (Temperature = "> 10ºC")-> warmgerm_graph 

ggplot(data=warmgerm_graph , aes(incubator, mean, color=macroclimate))+
  geom_point(size =2) +
  geom_signif(comparisons = list(c("Fellfield", "Snowbed")),annotations = "***", y_position = 0.45, tip_length = 0.04, color = "black", size = 1, textsize = 9) +
  #geom_signif(comparisons = list(c("Mediterranean", "Temperate")),annotations = "***", y_position = 0.65, tip_length = 0.04, color = "black", size = 1, textsize = 9) +
  #facet_grid (.~macroclimate) +
  geom_errorbar(data= warmgerm_graph , aes(incubator, mean, ymin = lower, ymax = upper, 
                                        color = macroclimate), width = 0.3, size =1.2) +
  #geom_segment (aes(x=1.5, y =0.24, xend=1.5, yend =0.4), color = "black", size = 1)+
  #geom_segment (aes(x=1.45, y = 0.24, xend =1.55, yend =0.24), color = "black", size = 1) +
  #geom_segment (aes(x=1.45, y = 0.4, xend =1.55, yend =0.4), color = "black", size = 1) +
  #annotate (geom ="text", x= 1.62, y = 0.3, label ="**", colour = "black", size = 9) +
  #scale_color_manual (name= "Incubator", values = c ("Fellfield"= "chocolate2", "Snowbed" ="deepskyblue3")) +
  scale_color_manual (name = "Macroclimate", values = c ("Mediterranean"= "darkgoldenrod1" , "Temperate" = "forestgreen")) +
  scale_y_continuous (limits = c(0,1), breaks = seq (0, 1, by= 0.2)) +
  labs(title = "Germination at T>10ºC", x = "Incubator", y = "Germination proportion") +
  theme_classic(base_size = 16)+
  theme(plot.title = element_text (hjust = 0.5, size = 30),
        strip.text = element_text(face = "italic", size = 24), 
        legend.position = "right",
        legend.title = element_text(size=20), 
        legend.text = element_text(size=16),
        legend.background = element_rect(fill="transparent",colour=NA),
        axis.title.y = element_text (size=18), 
        axis.title.x = element_text (size=18)) 


# TEMPERATURE GRAPH ####
rbind(nogrow_graph, coldgerm_graph, warmgerm_graph) %>% 
  mutate (Temperature = factor(Temperature, levels = c( "< 2ºC", "2-10ºC", "> 10ºC"))) -> graph #reorder levels
#rm(graph_date)
x11()
ggplot()+
  geom_point(data=graph, aes(Temperature, mean, color=incubator), size =2) +
  facet_grid (.~macroclimate) +
  geom_errorbar(data= graph, aes(Temperature, mean, ymin = lower, ymax = upper, color = incubator), width = 0.3, size =1) +
  geom_line(data= graph, aes(as.numeric(Temperature), mean, color = incubator), size =1) +
  scale_color_manual (name= "Incubator", values = c ("Fellfield"= "chocolate2", "Snowbed" ="deepskyblue3")) +
  scale_y_continuous (limits = c(0,1), breaks = seq (0, 1, by= 0.2)) +
  labs(title = "Germination response to Temperature", x = "Temperature", y = "Germination proportion") +
  theme_classic(base_size = 16)+
  theme(plot.title = element_text (hjust = 0.5, size = 30),
        strip.text = element_text(face = "italic", size = 24), 
        legend.position = "right",
        legend.title = element_text(size=20), 
        legend.text = element_text(size=16),
        legend.background = element_rect(fill="transparent",colour=NA),
        axis.title.y = element_text (size=18), 
        axis.title.x = element_text (size=18))
# T50 ####
t50model %>%
  merge(viables_germ) %>%
  filter(viablePER>25)%>%
  filter(germPER>0)%>%
  mutate(species= str_replace(species, "Cerastium sp", "Cerastium pumilum"))%>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>% 
  merge(species, by = c("code", "species"))-> t50_graph

ggplot()+
  geom_boxplot(data=t50_graph, aes(incubator, t50lm, fill=macroclimate), size =1) +
  #facet_grid (.~macroclimate) +
  #scale_fill_manual (name= "Incubator", values = c ("Fellfield"= "chocolate2", "Snowbed" ="deepskyblue3")) +
  scale_fill_manual (name = "Macroclimate",values = c ("Mediterranean"= "darkgoldenrod1" , "Temperate" = "forestgreen")) +
  labs(title = "Time to reach 50% germination", x = "Incubator", y = "Time (days)") +
  theme_classic(base_size = 16)+
  theme(plot.title = element_text (hjust = 0.5, size = 30),
        strip.text = element_text(face = "italic", size = 24), 
        legend.position = "right",
        legend.title = element_text(size=20), 
        legend.text = element_text(size=16),
        legend.background = element_rect(fill="transparent",colour=NA),
        axis.title.y = element_text (size=18), 
        axis.title.x = element_text (size=18))

ggplot(data=t50_graph, aes(macroclimate, t50lm, fill=incubator))+
  geom_boxplot() +
  #facet_grid (.~incubator) +
  #geom_signif(comparisons = list(c("Fellfield", "Snowbed")),annotations = "***", y_position = 420, tip_length = 0.04, color = "black", size = 1, textsize = 9) +
  #geom_signif(comparisons = list(c("Mediterranean", "Temperate")),annotations = "***", y_position = 420, tip_length = 0.04, color = "black", size = 1, textsize = 9) +
  geom_segment (aes(x=0.81, y = 430, xend =1.20, yend =430), color = "black", size = 1) +
  annotate (geom ="text", x= 1, y = 435, label ="***", colour = "black", size = 9) +
  geom_segment (aes(x=1.81, y = 430, xend =2.20, yend =430), color = "black", size = 1) +
  annotate (geom ="text", x= 2, y = 435, label ="***", colour = "black", size = 9) +
  geom_segment (aes(x=0.81, y = 330, xend =1.81, yend =330), color = "black", size = 1) +
  annotate (geom ="text", x= 1.4, y = 335, label ="***", colour = "black", size = 9) +
  geom_segment (aes(x=1.20, y = 350, xend =2.2, yend =350), color = "black", size = 1) +
  annotate (geom ="text", x= 1.7, y = 355, label ="***", colour = "black", size = 9) +
  scale_fill_manual (name= "Incubator", values = c ("Fellfield"= "chocolate2", "Snowbed" ="deepskyblue3")) +
  #scale_fill_manual (name = "Macroclimate", values = c ("Mediterranean"= "darkgoldenrod1" , "Temperate" = "forestgreen")) +
  scale_y_continuous (limits = c(0,460), breaks = seq (0, 450, by= 100)) +
  labs(title = "Time to reach 50% germination", x = "Macroclimate", y = "Time (days)") +
  theme_classic(base_size = 16)+
  theme(plot.title = element_text (hjust = 0.5, size = 30),
        strip.text = element_text(face = "italic", size = 24), 
        legend.position = "right",
        legend.title = element_text(size=20), 
        legend.text = element_text(size=16),
        legend.background = element_rect(fill="transparent",colour=NA),
        axis.title.y = element_text (size=18), 
        axis.title.x = element_text (size=18))
# Heat sum #####
t50_dates%>%
  group_by (species, code, incubator, petridish) %>%
  do (heatsum(.)) %>%
  mutate(species= str_replace(species, "Cerastium sp", "Cerastium pumilum"))%>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>% 
  merge(species, by = c("code", "species")) -> heatsum_graph

ggplot(data=heatsum_graph, aes(macroclimate, HS, fill=macroclimate))+
  geom_boxplot() +
  facet_grid (.~incubator) +
  geom_signif(comparisons = list(c("Mediterranean", "Temperate")),annotations = "***", y_position = 3000, tip_length = 0.04, color = "black", size = 1, textsize = 9) +
  #geom_segment (aes(x=0.81, y = 2500, xend =1.20, yend =2500), color = "black", size = 1) +
  #annotate (geom ="text", x= 1, y = 2550, label ="***", colour = "black", size = 9) +
  #geom_segment (aes(x=1.81, y = 2500, xend =2.20, yend =2500), color = "black", size = 1) +
  #annotate (geom ="text", x= 2, y = 2550, label ="***", colour = "black", size = 9) +
  #scale_fill_manual (name= "Incubator", values = c ("Fellfield"= "chocolate2", "Snowbed" ="deepskyblue3")) +
  scale_fill_manual (name = "Macroclimate",values = c ("Mediterranean"= "darkgoldenrod1" , "Temperate" = "forestgreen")) +
  labs(title = "Environmental heat sum to t50", x = "Macroclimate", y = "Time (days)") +
  theme_classic(base_size = 16)+
  theme(plot.title = element_text (hjust = 0.5, size = 30),
        strip.text = element_text(face = "italic", size = 24), 
        legend.position = "right",
        legend.title = element_text(size=20), 
        legend.text = element_text(size=16),
        legend.background = element_rect(fill="transparent",colour=NA),
        axis.title.y = element_text (size=18), 
        axis.title.x = element_text (size=18))
# delay S-F ####
t50model %>%
  merge (viables_germ)%>%
  filter(viablePER>25)%>%
  group_by (species, code, incubator) %>%
  summarise(t50lm = mean(t50lm)) %>%
  spread(incubator, t50lm) %>% 
  mutate(delayS_F = Snowbed - Fellfield) %>%  # NAs appear when species don't reach 50% germination (t50lm_days =Na)
  na.omit () %>%
  mutate(species= str_replace(species, "Cerastium sp", "Cerastium pumilum"))%>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>% 
  merge(species, by = c("code", "species")) ->delay_graph

ggplot(data=delay_graph, aes(macroclimate, delayS_F, fill=macroclimate))+
  geom_boxplot() +
  #facet_grid (.~incubator) +
  #geom_signif(comparisons = list(c("Mediterranean", "Temperate")),annotations = "***", y_position = 3000, tip_length = 0.04, color = "black", size = 1, textsize = 9) +
  #geom_segment (aes(x=0.81, y = 2500, xend =1.20, yend =2500), color = "black", size = 1) +
  #annotate (geom ="text", x= 1, y = 2550, label ="***", colour = "black", size = 9) +
  #geom_segment (aes(x=1.81, y = 2500, xend =2.20, yend =2500), color = "black", size = 1) +
  #annotate (geom ="text", x= 2, y = 2550, label ="***", colour = "black", size = 9) +
  #scale_fill_manual (name= "Incubator", values = c ("Fellfield"= "chocolate2", "Snowbed" ="deepskyblue3")) +
  scale_fill_manual (name = "Macroclimate",values = c ("Mediterranean"= "darkgoldenrod1" , "Temperate" = "forestgreen")) +
  labs(title = "Delay time to reach t50", x = "Macroclimate", y = "Time (days)") +
  theme_classic(base_size = 16)+
  theme(plot.title = element_text (hjust = 0.5, size = 30),
        strip.text = element_text(face = "italic", size = 24), 
        legend.position = "none",
        legend.title = element_text(size=20), 
        legend.text = element_text(size=16),
        legend.background = element_rect(fill="transparent",colour=NA),
        axis.title.y = element_text (size=18), 
        axis.title.x = element_text (size=18))

#### 2nd sow data transformation and graph####
# tidyverse transformation to account for the number of viable seeds per each specie and incubator
# summing up petridishes and accesions/populations of the same species (not taking into account weekly germination)
read.csv("data/Second_sow.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date)))) %>%
  group_by(species, incubator, code, petridish) %>%
  filter(date == max(date)) %>%
  select(species, incubator, code, petridish, viable) %>%
  group_by(species, incubator) %>%
  summarise(viable = sum(viable)) -> secondsowviables
###  germination curves 2nd sow
read.csv("data/Second_sow.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date)))) %>%
  group_by(species, incubator, time) %>%
  summarise(germinated = sum(germinated)) %>%
  mutate(germinated = cumsum(germinated)) %>%
  merge(viables) %>%
  mutate(germination = germinated/viable) %>%
  filter(species == "Thymus praecox") %>%
  ggplot(aes(time, germination, color = incubator, fill = incubator)) +
  geom_line(size = 1) +
  scale_color_manual (name= "incubator", values = c ("Fellfield"= "red", "Snowbed" ="turquoise")) +
  facet_wrap(~ species, scales = "free_x", ncol = 2) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Time (days)", y = "Germination proportion") +
  theme(strip.text = element_text(face = "italic", size = 16), 
        legend.title = element_text(size=16), 
        legend.text = element_text(size=12),
        axis.title.y = element_text (size=14), 
        axis.title.x = element_text (size=14)) + 
  geom_vline(xintercept = 39, linetype = "dashed", size= 1) +
  geom_vline(xintercept = 165, linetype = "dashed", size= 1, color = "red") +
  geom_vline(xintercept = 207, linetype = "dashed", size= 1, color = "turquoise") 





#### SPECIES GERMINATION CURVES ####
# tidyverse transformation to account for the number of viable seeds per each specie and incubator
# summing up petridishes and accesions/populations of the same species (not taking into account weekly germination)
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date)))) %>%
  group_by(species, incubator, code, petridish) %>%
  filter(date == max(date)) %>%
  select(species, incubator, code, petridish, viable) %>%
  group_by(species, incubator) %>%
  summarise(viable = sum(viable)) -> viables_sp
# write.csv (viables,"results/viables.csv", row.names = FALSE )
# tidyverse modification to have the accumulated germination along the whole experiment + ggplot
x11()
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date)))) %>%
  group_by(species, incubator, time) %>%
  summarise(germinated = sum(germinated)) %>%
  mutate(germinated = cumsum(germinated)) %>%
  merge(viables_sp) %>%
  mutate(germination = germinated/viable) %>%
  filter(species == "Thymus praecox") %>% ### change species name
  ggplot(aes(time, germination, color = incubator, fill = incubator)) +
  geom_line(size = 1.5) +
  scale_color_manual (name= "Incubator", values = c ("Fellfield"= "chocolate2", "Snowbed" ="deepskyblue3")) +
  facet_wrap(~ species, scales = "free_x", ncol = 2) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Time (days)", y = "Germination proportion") +
  theme(strip.text = element_text(face = "italic", size = 20), 
        legend.title = element_text(size=18), 
        legend.text = element_text(size=14),
        axis.title.y = element_text (size=16), 
        axis.title.x = element_text (size=16)) + 
  geom_vline(xintercept = 122, linetype = "dashed", size= 1.5) +
  geom_vline(xintercept = 248, linetype = "dashed", size= 1.5, color = "chocolate2") +
  geom_vline(xintercept = 290, linetype = "dashed", size= 1.5, color = "deepskyblue3") 
  
# Save the plots
ggsave(filename = "results/FigAgrostistileni.png", Agrostistileni, path = NULL, 
       scale = 1, width = 180, height = 180, units = "mm", dpi = 600)
