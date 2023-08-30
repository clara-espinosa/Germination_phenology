library(FD); library(vegan); library(FactoMineR); library(emmeans);
library(tidyverse); library(ggrepel); library(cowplot);library(ggpubr);
library (binom);library (ggsignif);library (rstatix); library (stringr);
library(plyr);library(patchwork)
Sys.setlocale("LC_ALL","English")
## 
#### temperature experimental programs VISUALIZATION ####
read.csv("data/date_temp.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y")) %>%
  mutate(date = as.POSIXct(date))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  mutate(incubator = as.factor(incubator)) %>%
  as.data.frame() -> temp 
str(temp)
summary(temp)
# visualization 1
# mulitple lines with shaded areas in between
x11()
Sys.setlocale("LC_ALL", "English")
ggplot(temp, aes ()) +
  geom_line (aes (x=date, y=Tmin, colour =incubator), size =1.25) +
  geom_line (aes (x=date, y=Tmax, colour = incubator), size =1.25) +
  geom_ribbon (data = temp, aes (x=date, ymin =Tmin, ymax=Tmax, colour = incubator, fill = incubator), alpha =0.5) +
  scale_fill_manual (name= "Incubator", values =c("chocolate2", "deepskyblue3")) +
  scale_color_manual (name= "Incubator", values =c("chocolate2", "deepskyblue3")) +
  scale_y_continuous (limits = c(-5,28), breaks = seq (-5, 25, by= 5)) +
  labs (title = "Experimental temperature regimes", y= "Temperature ºC", x = "Date") + 
  theme_classic(base_size = 20) +
  theme (plot.title = element_text ( size = 36), #hjust = 0.5,
        axis.title.y = element_text (size=32), 
         axis.title.x = element_text (size=32), 
         axis.text.x= element_text (size=28),
         legend.title = element_text(size = 28),
         legend.text = element_text (size =26),
        legend.position = "top") +
  geom_vline(xintercept = as.POSIXct(as.Date("2021-11-12")), linetype = "dashed", size =1.25) +
  annotate (geom ="text", x= as.POSIXct(as.Date("2021-10-17")), y = 28, label ="Autumn \n germination", size = 8, fontface ="bold") +
  geom_vline(xintercept = as.POSIXct(as.Date("2022-06-15")), linetype = "dashed", size =1.25) +
  annotate (geom ="text", x= as.POSIXct(as.Date("2022-05-20")), y = 28, label ="Spring\n germination", size = 8, fontface ="bold") +
  geom_vline(xintercept = as.POSIXct(as.Date("2022-09-19")), linetype = "dashed", size =1.25) +
  annotate (geom ="text", x= as.POSIXct(as.Date("2022-08-25")), y = 28, label ="Summer\n germination", size = 8, fontface ="bold") +
  geom_hline(yintercept=0, linetype ="dashed", size =1, colour = "red") +
  geom_segment (aes(x=as.POSIXct(as.Date("2021-11-12")), y =-2.8, xend=as.POSIXct(as.Date("2022-04-04")), yend =-2.8), color = "chocolate2", size = 2) +
  geom_segment (aes(x=as.POSIXct(as.Date("2021-11-12")), y = -3.8, xend =as.POSIXct(as.Date("2022-05-26")), yend =-3.8), color = "deepskyblue3", size = 2) +
  annotate (geom ="text", x= as.POSIXct(as.Date("2022-01-15")), y = -5, label ="Winter conditions", colour = "black", size = 8.5, fontface ="bold") 
  
########### GERMINATION RATE ###################################################
#### germination peaks (not used in the final paper) #####
data.frame(community = c("Mediterranean","Mediterranean","Mediterranean","Temperate", "Temperate"),
           incubator  = c("Fellfield", "Fellfield","Snowbed", "Fellfield", "Fellfield"),
           time = c(105, 247, 105, 105, 247),
           germinated = c(0,0,0,0,0)) -> data_vis
x11()
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date)))) %>%
  merge (species) %>%
  group_by(community, incubator, date) %>%
  summarise(germinated = sum(germinated)) %>%
  bind_rows (data_vis)%>%
  merge(viables_community) %>%
  mutate(germination = germinated/viable) %>%
  ggplot(aes(time, germination, color = incubator, fill = incubator)) +
  geom_line(size = 2) +
  scale_color_manual (name= "Incubator", values = c ("Fellfield"= "chocolate2", "Snowbed" ="deepskyblue3")) +
  facet_wrap(~ community, scales = "free_x", ncol = 2) +
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
detach(package:plyr) # to make summarise work properly
data.frame(community = c("Mediterranean","Mediterranean","Mediterranean","Temperate", "Temperate"),
           incubator  = c("Fellfield", "Fellfield","Snowbed", "Fellfield", "Fellfield"),
           date = c("2021-11-12", "2022-04-03", "2021-11-12", "2021-11-12", "2022-04-03"),
           germinated = c(0,0,0,0,0)) -> data_vis
data_vis %>%
  mutate(date = as.POSIXct(date)) -> data_vis
str(data_vis)
x11()
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(date = as.POSIXct(date))%>%
  merge (species) %>%
  group_by(community, incubator, date) %>%
  bind_rows (data_vis)%>%
  summarise(germinated = sum(germinated)) %>%
  mutate(germinated = cumsum(germinated)) %>%
  merge(viables_community) %>%
  mutate(germination = germinated/viable) %>%
  filter (community == "Mediterranean") %>%
  ggplot(aes(date, germination, color = incubator, fill = incubator)) +
  geom_line(size = 2.5) +
  scale_color_manual (name= "Incubator", values = c ("Fellfield"= "chocolate2", "Snowbed" ="deepskyblue3")) +
  scale_y_continuous (limits = c(0,1), breaks = seq (0, 1, by= 0.25)) +
  #facet_wrap(~ community, scales = "free_x", ncol = 2) +
  labs(title = "Mediterranean", x = "Date", y = "Germination proportion") +
  theme_classic(base_size = 16) +
  theme (plot.title = element_text (hjust = 0.5, size = 32), #
          axis.title.y = element_text (size=28), 
          axis.title.x = element_text (size=28), 
          axis.text.x= element_text (size=24),
          plot.margin = margin(r= 40),
          legend.title = element_text(size = 32),
          legend.text = element_text (size =30),
          legend.position = "none") -> cumulative_Mediterranean
  #geom_vline(xintercept = as.POSIXct(as.Date("2021-11-12")), linetype = "dashed", size =1.25) +
  #annotate (geom ="text", x= as.POSIXct(as.Date("2021-10-12")), y = 0.8, label ="Autumn", size = 5, fontface ="bold") +
  #geom_vline(xintercept = as.POSIXct(as.Date("2022-06-15")), linetype = "dashed", size =1.25) +
  #annotate (geom ="text", x= as.POSIXct(as.Date("2022-05-15")), y = 0.8, label ="Spring", size = 5, fontface ="bold") +
  #geom_vline(xintercept = as.POSIXct(as.Date("2022-09-19")), linetype = "dashed", size =1.25) 
  #annotate (geom ="text", x= as.POSIXct(as.Date("2022-08-17")), y = 0.8, label ="Summer", size = 5, fontface ="bold") +
  #geom_segment (aes(x=as.POSIXct(as.Date("2021-11-25")), y =0.05, xend=as.POSIXct(as.Date("2022-04-04")), yend =0.05), color = "chocolate2", size = 1.25) +
  #geom_segment (aes(x=as.POSIXct(as.Date("2021-11-12")), y = 0.03, xend =as.POSIXct(as.Date("2022-05-26")), yend =0.03), color = "deepskyblue3", size = 1.25) + 
  #annotate (geom ="text", x= as.POSIXct(as.Date("2022-01-15")), y = 0, label ="Winter period", colour = "black", size = 5, fontface ="bold") 

ggarrange(cumulative_Mediterranean, cumulative_Temperate, ncol =2, nrow= 1,common.legend = TRUE, legend = "bottom") 

  
########### TOTAL GERMINATION + TIMING  #######################################
# Total germination (test representation) ####
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

#value graph
x11()
ggplot(data=totalgerm_graph, aes(incubator, mean, fill=incubator))+
  geom_bar(stat ="identity", width = 0.8) +
  facet_grid (.~community) +
  #geom_signif(comparisons = list(c("Fellfield", "Snowbed")),annotations = "***", y_position = 0.95, tip_length = 0.05, color = "black", size = 1, textsize = 9) +
  geom_errorbar(data= totalgerm_graph, aes(incubator, mean, ymin = lower, ymax = upper), color = "black",width = 0.2, size =1) +
  # geom_segment (aes(x=1.5, y =0.72, xend=1.5, yend =0.82), color = "black", size = 1)+
  # geom_segment (aes(x=1.45, y = 0.72, xend =1.55, yend =0.72), color = "black", size = 1) +
 # geom_segment (aes(x=1.45, y = 0.82, xend =1.55, yend =0.82), color = "black", size = 1) +
 # annotate (geom ="text", x= 1.60, y = 0.79, label =".", colour = "black", size = 9) +
  scale_fill_manual (name= "Incubator", values = c ("Fellfield"= "chocolate2", "Snowbed" ="deepskyblue3")) +
  scale_y_continuous (limits = c(0,1), breaks = seq (0, 1, by= 0.2)) +
  labs(title = "Total germination", x = "Incubator", y = "Germination proportion") +
  #ggthemes::theme_tufte(base_size = 16) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text (hjust = 0.5, size = 30),
      strip.text = element_text(face = "italic", size = 24),
      #strip.text = element_blank(),
      legend.position = "none",
      panel.background = element_rect(fill="transparent",colour=NA),
      #panel.background = element_rect(color = "grey96", fill = "grey96"),
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

x11()
ggplot(data=autumn_graph, aes(incubator, mean, color=community))+
  geom_point(size =2) +
  geom_signif(comparisons = list(c("Fellfield", "Snowbed")),annotations = "***", y_position = 0.65, tip_length = 0.04, color = "black", size = 1, textsize = 9) +
  #geom_signif(comparisons = list(c("Mediterranean", "Temperate")),annotations = "***", y_position = 0.65, tip_length = 0.04, color = "black", size = 1, textsize = 9) +
  #facet_grid (.~community) +
  geom_errorbar(data= autumn_graph, aes(incubator, mean, ymin = lower, ymax = upper, 
                                        color = community), width = 0.3, size =1.2) +
  geom_segment (aes(x=1.5, y =0.24, xend=1.5, yend =0.4), color = "black", size = 1)+
  geom_segment (aes(x=1.45, y = 0.24, xend =1.55, yend =0.24), color = "black", size = 1) +
  geom_segment (aes(x=1.45, y = 0.4, xend =1.55, yend =0.4), color = "black", size = 1) +
  annotate (geom ="text", x= 1.62, y = 0.3, label ="**", colour = "black", size = 9) +
  #scale_color_manual (name= "Incubator", values = c ("Fellfield"= "chocolate2", "Snowbed" ="deepskyblue3")) +
  scale_color_manual (name = "community", values = c ("Mediterranean"= "darkgoldenrod1" , "Temperate" = "forestgreen")) +
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
x11()
ggplot(data= spring_graph, aes(incubator, mean, color=community))+
  geom_point(size =2) +
  geom_signif(comparisons = list(c("Fellfield", "Snowbed")), annotations = "***", y_position = 0.5, tip_length = 0.04, color = "black", size = 1, textsize = 9) +
  #facet_grid (.~community) +
  geom_errorbar(data=  spring_graph, aes(incubator, mean, ymin = lower, ymax = upper, 
                                        color = community), width = 0.3, size =1.2) +
  #scale_color_manual (name= "Incubator", values = c ("Fellfield"= "chocolate2", "Snowbed" ="deepskyblue3")) +
  scale_color_manual (name = "community", values = c ("Mediterranean"= "darkgoldenrod1" , "Temperate" = "forestgreen")) +
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

x11()
ggplot(data= summer_graph , aes(incubator, mean, color=community))+
  geom_point(size =2) +
  #geom_signif(comparisons = list(c("Mediterranean", "Temperate")),annotations = "***", y_position = 0.45, tip_length = 0.05, color = "black", size = 1, textsize = 9) +
  geom_signif(comparisons = list(c("Fellfield", "Snowbed")), annotations = "***", y_position = 0.45, tip_length = 0.04, color = "black", size = 1, textsize = 9) +
  #facet_grid (.~community) +
  geom_errorbar(data=  summer_graph , aes(incubator, mean, ymin = lower, ymax = upper, 
                                         color = community), width = 0.3, size =1.2) +
  #geom_segment (aes(x=1.5, y =0.08, xend=1.5, yend =0.13), color = "black", size = 1)+
  #geom_segment (aes(x=1.45, y = 0.08, xend =1.55, yend =0.08), color = "black", size = 1) +
  #geom_segment (aes(x=1.45, y = 0.13, xend =1.55, yend =0.13), color = "black", size = 1) +
  #annotate (geom ="text", x= 1.62, y = 0.1, label ="**", colour = "black", size = 9) +
  #scale_color_manual (name= "Incubator", values = c ("Fellfield"= "chocolate2", "Snowbed" ="deepskyblue3")) +
  scale_color_manual (name = "community", values = c ("Mediterranean"= "darkgoldenrod1" , "Temperate" = "forestgreen")) +
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
  facet_grid (.~community) +
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

# NON-GROWING season germination <2ºC ####
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

ggplot(data= nogrow_graph, aes(incubator, mean, color=community))+
  geom_point(size =2) +
  geom_signif(comparisons = list(c("Fellfield", "Snowbed")), annotations = "***", y_position = 0.3, tip_length = 0.04, color = "black", size = 1, textsize = 9) +
  #facet_grid (.~community) +
  geom_errorbar(data=  nogrow_graph, aes(incubator, mean, ymin = lower, ymax = upper, 
                                         color = community), width = 0.3, size =1.2) +
  #scale_color_manual (name= "Incubator", values = c ("Fellfield"= "chocolate2", "Snowbed" ="deepskyblue3")) +
  scale_color_manual (name = "community", values = c ("Mediterranean"= "darkgoldenrod1" , "Temperate" = "forestgreen")) +
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

# T50 ####
t50model %>%
  merge(viables) %>%
  filter(viablePER>25)%>%
  filter(germPER>0)%>%
  merge(species, by = c("code", "species"))-> t50_graph
t50_graph %>%
  group_by(incubator, community) %>%
  get_summary_stats(t50lm, type="full")

ggplot()+
  geom_boxplot(data=t50_graph, aes(incubator, t50lm, fill=community), size =1) +
  #facet_grid (.~community) +
  #scale_fill_manual (name= "Incubator", values = c ("Fellfield"= "chocolate2", "Snowbed" ="deepskyblue3")) +
  scale_fill_manual (name = "community",values = c ("Mediterranean"= "darkgoldenrod1" , "Temperate" = "forestgreen")) +
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

ggplot(data=t50_graph, aes(incubator, t50lm, fill=incubator))+
  geom_boxplot() +
  facet_grid (~community) +
  #geom_signif(comparisons = list(c("Fellfield", "Snowbed")),annotations = "***", y_position = 420, tip_length = 0.04, color = "black", size = 1, textsize = 9) +
  #geom_signif(comparisons = list(c("Mediterranean", "Temperate")),annotations = "***", y_position = 420, tip_length = 0.04, color = "black", size = 1, textsize = 9) +
  #geom_segment (aes(x=0.81, y = 430, xend =1.20, yend =430), color = "black", size = 1) +
  #annotate (geom ="text", x= 1, y = 435, label ="***", colour = "black", size = 9) +
  #geom_segment (aes(x=1.81, y = 430, xend =2.20, yend =430), color = "black", size = 1) +
  #annotate (geom ="text", x= 2, y = 435, label ="***", colour = "black", size = 9) +
  #geom_segment (aes(x=0.81, y = 330, xend =1.81, yend =330), color = "black", size = 1) +
  #annotate (geom ="text", x= 1.4, y = 335, label ="***", colour = "black", size = 9) +
  #geom_segment (aes(x=1.20, y = 350, xend =2.2, yend =350), color = "black", size = 1) +
  #annotate (geom ="text", x= 1.7, y = 355, label ="***", colour = "black", size = 9) +
  scale_fill_manual (name= "Incubator", values = c ("Fellfield"= "chocolate2", "Snowbed" ="deepskyblue3")) +
  #scale_fill_manual (name = "community", values = c ("Mediterranean"= "darkgoldenrod1" , "Temperate" = "forestgreen")) +
  scale_y_continuous (limits = c(0,460), breaks = seq (0, 450, by= 100)) +
  labs(title = "Time to reach 50% germination", x = "Incubator", y = "Time (days)") +
  theme_classic(base_size = 16)+
  theme(plot.title = element_text (hjust = 0.5, size = 30),
        strip.text = element_text(face = "bold", size = 12), 
        legend.position = "none",
        panel.background = element_rect(color = "grey96", fill = "grey96"),
        axis.title.y = element_text (size=18), 
        axis.title.x = element_text (size=18))-> fig5b
# Heat sum #####
t50_dates%>%
  group_by (species, code, incubator, petridish) %>%
  do (heatsum(.)) %>%
  merge(species, by = c("code", "species"))-> heatsum_graph
heatsum_graph %>%
  group_by(incubator, community) %>%
  get_summary_stats(HS, type= "full")
x11()
ggplot(data=heatsum_graph, aes(incubator, HS, fill=incubator))+
  geom_boxplot() +
  facet_grid (~community) +
  #geom_signif(comparisons = list(c("Mediterranean", "Temperate")),annotations = "*", y_position = 3000, tip_length = 0.04, color = "black", size = 1, textsize = 9) +
  #geom_segment (aes(x=0.81, y = 2500, xend =1.20, yend =2500), color = "black", size = 1) +
  #annotate (geom ="text", x= 1, y = 2550, label ="*", colour = "black", size = 9) +
  #geom_segment (aes(x=1.81, y = 2500, xend =2.20, yend =2500), color = "black", size = 1) +
  #annotate (geom ="text", x= 2, y = 2550, label ="*", colour = "black", size = 9) +
  scale_fill_manual (name= "Incubator", values = c ("Fellfield"= "chocolate2", "Snowbed" ="deepskyblue3")) +
  #scale_fill_manual (name = "community",values = c ("Mediterranean"= "darkgoldenrod1" , "Temperate" = "forestgreen")) +
  labs(title = "Environmental heat sum to t50", x = "Incubator", y = "Degrees ºC") +
  theme_classic(base_size = 16)+
  theme(plot.title = element_text (hjust = 0.5, size = 30),
        strip.text = element_text(face = "bold", size = 12), 
        legend.position = "none",
        panel.background = element_rect(color = "grey96", fill = "grey96"),
        axis.title.y = element_text (size=18), 
        axis.title.x = element_text (size=18)) -> fig5c
##### effect size#####
library(viridis)
x11()
read.csv("data/test_effectsize.csv", sep=";") %>%
  #filter(!Trait== "Heat sum")%>%
  #filter(!Trait == "t50") %>%
  convert_as_factor(Trait, community, terms) %>%
  mutate (Trait = fct_relevel(Trait, "Environmental heat sum", "t50","Germination in winter conditions", 
                              "Summer germination", "Spring germination","Autumn germination", "Total germination" )) -> effect_size

ggplot(effect_size, aes(x= Trait, y =effect_size, ymin = L95, ymax = U95, color = Trait))+
  geom_point( size = 3) +
  geom_errorbar (width = 0.2, size =1.2) +
  facet_wrap (~ community, ncol = 2, nrow =7, scales = "free_x") +
  scale_color_manual (values = c("#AC1926", "#891171", "#33407D", "#077395", "#00BC7F", "#AADB41", "#FDE333")) +
  scale_x_discrete(labels = function(Trait) str_wrap(Trait, width = 13)) +
  scale_y_continuous (limits = c(-5,5), breaks = seq (-4, 4, by = 2)) +
  geom_hline(yintercept = 0, linetype = "dashed", size =1, color = "black") +
  coord_flip() +
  labs(title = "A", y = "Effect size") + #title = "", x =""
  theme_classic(base_size = 14) +
  #ggthemes::theme_tufte(base_size = 16) +
  theme(plot.title = element_text (size = 22),
        strip.text = element_text(face = "bold", size = 22),
        #strip.background = element_rect(fill = )
        
        #strip.text = element_blank(),
        legend.position = "none",
        #legend.title = element_text(size=20), 
        #legend.text = element_text(size=16),
        #legend.background = element_rect(fill="transparent",colour=NA),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.x = element_text (size=18)) +
  annotate("text", x= 0.70, y = -2.5, label = "Higher in Fellfield", size = 3.5, color = "chocolate2" , fontface ="bold" ) +
  #annotate ("rect", xmin = 0.5, xmax= 7.7, ymin = -5, ymax = -0.1, alpha =0.2, color = "chocolate2", fill = "chocolate2") +
  annotate ("segment", x = 0.5, xend = 0.5, y= -0.1, yend=-4.9, colour= "chocolate2", size = 2, arrow = arrow()) +
  annotate("text", x= 0.70, y = 2.6, label = "Higher in Snowbed", size = 3.5, colour= "deepskyblue3", fontface ="bold" ) +
  annotate ("segment", x = 0.5, xend = 0.5, y= 0.1, yend= 4.9, colour= "deepskyblue3", size = 2, arrow = arrow()) ->effect_size_graph

ggsave(filename = "results/Fig4v3(2).png", effect_size_graph, path = NULL, 
       scale = 1, width = 360, height = 360, units = "mm", dpi = 600)
######  mean value graph #####
library("plyr")
x11()
read.csv("data/meanvalues_graph.csv", sep = ";")%>%
  convert_as_factor(trait, community, incubator) %>%
  mutate (trait = fct_relevel(trait, "Total germination", "Autumn germination","Spring germination", "Summer germination", 
                             "Germination in winter conditions",  "t50","Environmental heat sum")) -> mean_values
ggplot(mean_values, aes(incubator, mean, fill=incubator))+
  geom_col(position = position_dodge(0.7), width = 0.75) +
  facet_grid (trait ~ community, scales = "free_y", labeller = label_wrap_gen (width = 13)) +
  #facet_wrap (trait ~ community, ncol = 2, nrow =7, scales = "free_y") +
  #geom_signif(comparisons = list(c("Fellfield", "Snowbed")),annotations = "***", y_position = 0.95, tip_length = 0.05, color = "black", size = 1, textsize = 9) +
  geom_errorbar(aes(incubator, mean, ymin = lower, ymax = upper), color = "black",width = 0.2, size =1) +
  # geom_segment (aes(x=1.5, y =0.72, xend=1.5, yend =0.82), color = "black", size = 1)+
  # geom_segment (aes(x=1.45, y = 0.72, xend =1.55, yend =0.72), color = "black", size = 1) +
  # geom_segment (aes(x=1.45, y = 0.82, xend =1.55, yend =0.82), color = "black", size = 1) +
  # annotate (geom ="text", x= 1.60, y = 0.79, label =".", colour = "black", size = 9) +
  scale_fill_manual (name= "Incubator", values = c ("Fellfield"= "chocolate2", "Snowbed" ="deepskyblue3")) +
  #labs(title = "B", x = "Incubator", y = "Degres (ºC) Time (days)                               Germination proportion                                ") + #title = "",
  
  labs(title = "B", x = "Incubator", y = "Degres (ºC)         Time (days)                                                                           Germination proportion                                ") + #title = "",
  #ggthemes::theme_tufte(base_size = 16) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text (size = 22),
        #strip.text = element_blank(),
        strip.text.x = element_text(face = "bold", size = 22),
        strip.text.y = element_text(size = 14, angle = 360),
        legend.position = "none",
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title.y = element_text (size=16), 
        axis.title.x = element_text (size=16)) -> mean_values_graph

ggsave(filename = "results/Fig5v3(2).png", mean_values_graph, path = NULL, 
       scale = 1, width = 240, height = 360, units = "mm", dpi = 600)

Fig4_5 <- effect_size_graph + mean_values_graph

ggsave(filename = "results/Fig4_5(2).png", Fig4_5, path = NULL, 
       scale = 1, width = 400, height = 360, units = "mm", dpi = 600)


label_wrap_gen <- function (width = 13) { ### function to make labels in several lines with facet grid
  function(variable, value) {
    laply(strwrap (as.character(value), width = width, simplify = FALSE), paste, collapse ="\n")
  }
}
mean_values %>%
  filter (trait == "t50") %>%
  #filter (trait == "Environmental heat sum") %>%
  #filter (!trait == "Total germination") %>%
  #filter (!trait == "Autumn germination") %>%
  #filter (!trait == "Spring germination") %>%
  #filter (!trait == "Summer germination") %>%
  #filter (!trait == "Germination in winter conditions") %>%
ggplot(aes(incubator, mean, fill=incubator))+
  geom_bar(stat ="identity", width = 0.5) +
  facet_grid (trait ~ community, scales = "free_y", labeller = label_wrap_gen (width = 13)) +
  #facet_wrap (trait ~ community, ncol = 2, nrow =7, scales = "free_y") +
  #geom_signif(comparisons = list(c("Fellfield", "Snowbed")),annotations = "***", y_position = 0.95, tip_length = 0.05, color = "black", size = 1, textsize = 9) +
  geom_errorbar(aes(incubator, mean, ymin = lower, ymax = upper), color = "black",width = 0.2, size =1) +
  # geom_segment (aes(x=1.5, y =0.72, xend=1.5, yend =0.82), color = "black", size = 1)+
  # geom_segment (aes(x=1.45, y = 0.72, xend =1.55, yend =0.72), color = "black", size = 1) +
  # geom_segment (aes(x=1.45, y = 0.82, xend =1.55, yend =0.82), color = "black", size = 1) +
  # annotate (geom ="text", x= 1.60, y = 0.79, label =".", colour = "black", size = 9) +
  scale_fill_manual (name= "Incubator", values = c ("Fellfield"= "chocolate2", "Snowbed" ="deepskyblue3")) +
  labs( x = "Incubator", y = "Time (days)") + #title = "",
  #ggthemes::theme_tufte(base_size = 16) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text (hjust = 0.5, size = 30),
        strip.text.x = element_blank(),
        #strip.text.x = element_text(face = "bold", size = 20),
        strip.text.y = element_text(size = 14, angle = 360),
        legend.position = "none",
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title.y = element_text (size=18), 
        axis.title.x = element_text (size=18)) ->fig5b_2

ggarrange(fig5a_2,fig5b_2, fig5c_2, ncol = 1, nrow = 3, common.legend = TRUE)               
      # Label of the line plot
library(patchwork);library(gapminder)
 fig5 <- fig5a + (fig5b/fig5c)
#### SPECIES GERMINATION CURVES ####
# tidyverse transformation to account for the number of viable seeds per each specie and incubator
# summing up petridishes and accesions/populations of the same species (not taking into account weekly germination)
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date)))) %>%
  dplyr::group_by(species, incubator, accession, code, petridish) %>%
  filter(date == max(date)) %>%
  select(species, incubator, accession, code, petridish, viable) %>%
  dplyr::group_by(species, accession, incubator) %>%
  dplyr::summarise(viable = sum(viable)) -> viables_sp
# write.csv (viables,"results/viables.csv", row.names = FALSE )
# tidyverse modification to have the accumulated germination along the whole experiment + ggplot
read.csv("data/extrapoint_visualization.csv", sep = ";")-> data_vis 
x11()

read.csv("data/all_data.csv", sep = ";") %>%
  rbind(data_vis)%>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  spread(date, germinated, fill = 0) %>% # wide format for dates, and fill Na with 0
  gather ("date", "germinated", 8: last_col() )%>% # back in long format frrom 8th colum to the last
  arrange (species, accession, code, incubator, petridish, date)%>% # sort row observations this way
  mutate(date = as.POSIXct(date))%>%
  merge (species) %>%
  dplyr::group_by(community, species, accession,  incubator, date) %>%
  dplyr::summarise(germinated = sum(germinated)) %>%
  dplyr::mutate(germinated = cumsum(germinated)) %>%
  merge(viables_sp) %>%
  mutate(germination = germinated/viable)  %>%
  mutate(ID = paste(community, accession)) %>%
  filter(species == "Cerastium ramosissimum") %>% ### change species name
  ggplot(aes(date, germination, color = incubator, fill = incubator)) +
  geom_line(size = 2) +
  scale_color_manual (name= "Incubator", values = c ("Fellfield"= "chocolate2", "Snowbed" ="deepskyblue3")) +
  facet_wrap(~ accession, scales = "free_x", ncol = 2) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title= "Cerastium ramosissimum", x = "Time ", y = "Germination proportion") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text (size = 30),
        strip.text = element_text (size = 24, face = "italic"),
        axis.title.y = element_text (size=24), 
        axis.title.x = element_text (size=24), 
        axis.text.x= element_text (size=18, angle = 75, vjust = 0.5),
        legend.position = "right",
        plot.margin = margin(t=0.5, l =0.5, b = 0.5, r =0.5, unit ="cm")) #-> C
  
  
# Save the plots
ggsave(filename = "results/Species germination curves/Thymus praecox.png", C, path = NULL, 
       scale = 1, width = 360, height = 360, units = "mm", dpi = 600)

#loop for many graphs 
for (var in unique(germination_curves$species)) {
  curve_plot = ggplot(germination_curves[germination_curves$species==var,], aes(x = date,
                                                      y = germination, ymin = 0, ymax = 1,
                                                      color = incubator, fill = incubator)) +
    geom_line(size = 2) +
    scale_color_manual (name= "Incubator", values = c ("Fellfield"= "chocolate2", "Snowbed" ="deepskyblue3")) +
    labs(title= var, x = "Time ", y = "Germination proportion") +
    scale_y_continuous(limits = c(0, 1),  breaks = c(0, 0.25, 0.50, 0.75, 1)) +
    facet_wrap(~ accession, nrow = 2) +
    theme_classic(base_size = 14) +
    theme(plot.title = element_text (size = 30),
          strip.text = element_text (size = 24, face = "italic"),
          axis.title.y = element_text (size=24), 
          axis.title.x = element_text (size=24), 
          axis.text.x= element_text (size=18, angle = 75, vjust = 0.5),
          legend.position = "right",
          plot.margin = margin(t=0.5, l =0.5, b = 0.5, r =0.5, unit ="cm"))
  ggsave(curve_plot, file = paste0("Germination ", var,".png"),
         path = NULL, scale = 1, width = 360, height = 250, units = "mm", dpi = 600)
}
