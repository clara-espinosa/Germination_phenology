library(FD); library(vegan); library(FactoMineR); library(emmeans);
library(tidyverse); library(ggrepel); library(cowplot);library(ggpubr);
library (binom);library (ggsignif);library (rstatix); library (stringr);
library(plyr);library(patchwork);library(scales)
theme_set(theme_cowplot(font_size = 10)) 
Sys.setlocale("LC_ALL","English")

## FIG 2B bioclim study area ####
# visualization test for figure 1 (Study area localization and description)
read.csv("data/Bioclim_study_area.csv", sep= ";") %>%
  #mutate(community = fct_relevel (community, "Temperate","Mediterranean"))%>%
  ggplot()+
  geom_boxplot (aes(x= bio1, y= community, fill = community), color = "black", position = "dodge2")+
  coord_flip() +
  scale_fill_manual (name= "Community", values = c ("Mediterranean"= "darkgoldenrod1" , "Temperate" ="forestgreen")) +
  labs(title = "Mean annual air temperature", x= "Temperature (ºC)", y= "Community", tag= "(b)")+
  theme_classic(base_size = 20) +
  theme (plot.title = element_text (size = 23), #hjust = 0.5,
         axis.title.y = element_text (size=14),
         axis.text.y = element_text (size = 13),
         plot.tag.position = c(0.015,1),
         axis.title.x = element_blank(), 
         axis.text.x= element_text (size = 14, color = "black"),
         legend.position = "none") -> fig1b_1;fig1b_1

read.csv("data/Bioclim_study_area.csv", sep= ";") %>%
  #mutate(community = fct_relevel (community, "Temperate","Mediterranean"))%>%
  ggplot()+
  geom_boxplot (aes(x= bio17, y= community, fill = community), color = "black", position = "dodge2")+
  coord_flip()+
  scale_fill_manual (name= "Community", values = c ("Mediterranean"= "darkgoldenrod1" , "Temperate" ="forestgreen")) +
  labs(title = "Summer precipitation", x= "Precipitation (kg m-2)", y= "Community")+
  theme_classic(base_size = 16) +
  theme (plot.title = element_text (size = 23), #hjust = 0.5,
         axis.title.y = element_text (size=14),
         axis.text.y = element_text (size = 13),
         axis.title.x = element_blank(), 
         axis.text.x= element_text (size = 14, color = "black"),
         legend.position = "none")-> fig1b_2;fig1b_2

## Fig 2C graph to compare with temperature programs ####
x11()
read.csv("data/weekly_picos_graph.csv", sep = ",") %>%
  #mutate(order2 = paste(month, week))%>%
  ggplot(aes ()) +
  geom_line (aes (x=order, y=N, colour = Site), linewidth =1.25) + #order for sort as experiment
  geom_line (aes (x=order, y=X, colour = Site), linewidth =1.25) +
  geom_ribbon (aes (x=order, ymin =N, ymax=X, colour = Site, fill = Site), alpha =0.5) +
  scale_color_manual (name = "Microsites", values = c ("Fellfield site"= "chocolate2" , "Snowbed site" = "deepskyblue3")) +
  scale_fill_manual (name = "Microsites", values = c ("Fellfield site"= "chocolate2" , "Snowbed site" = "deepskyblue3")) +
  scale_y_continuous (limits = c(-3,25), breaks = seq (0, 25, by= 5)) +
  labs (title = "Soil temperature records 2008 - 2019", y= "Temperature (ºC)", x = "Time (weeks)", tag = "(c)") + #
  theme_classic() +
  theme (plot.title = element_text ( size = 24), #hjust = 0.5,
         axis.title.y = element_text (size=20), 
         axis.text.y = element_text (size = 18),
         axis.title.x = element_blank(), 
         axis.text.x= element_text (size = 18, color = "black"),
         plot.tag.position = c(0.015,1),
         legend.title = element_text (size =14),
         legend.text = element_text (size =13),
         legend.background = element_rect(fill="transparent"),
         legend.position = c(0.9, 0.28)) + # legend.position = "top") + 
  geom_hline(yintercept=0, linetype ="dashed", size =1, colour = "red") -> fig1c; fig1c

## GERMINATION EXPERIMENTS FIG 2D####
#### dataframe  AND VISUALIZATION for temperature programs 
read.csv("data/date_temp.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y")) %>%
  mutate(date = as.POSIXct(date))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  mutate(incubator = as.factor(incubator), 
         germination_period = as.factor (germination_period)) %>%
  as.data.frame() -> temp 
str(temp)
summary(temp)

x11()
Sys.setlocale("LC_ALL", "English")
ggplot(temp, aes ()) +
  geom_line (aes (x=date, y=Tmin, colour =incubator), linewidth =1.25) +
  geom_line (aes (x=date, y=Tmax, colour = incubator), linewidth =1.25) +
  geom_ribbon (data = temp, aes (x=date, ymin =Tmin, ymax=Tmax, colour = incubator, fill = incubator), alpha =0.5) +
  scale_fill_manual (name= "Incubator", values =c("chocolate2", "deepskyblue3")) +
  scale_color_manual (name= "Incubator", values =c("chocolate2", "deepskyblue3")) +
  scale_y_continuous (limits = c(-3,25), breaks = seq (0, 25, by= 5)) +
  scale_x_datetime(date_breaks = "2 month", date_labels = "%b %y")+
  labs (title = "Experimental temperature programs", y= "Temperature (ºC)", x = "Date", tag= "(d)") + #
  theme_classic() +
  theme (plot.title = element_text ( size = 24), #hjust = 0.5,
         axis.title.y = element_text (size=20), 
         axis.title.x = element_blank(), 
         axis.text.x= element_text (size = 18, color = "black"),
         axis.text.y = element_text(size = 18),
         plot.tag.position = c(0.015,1),
         legend.title = element_text (size =14),
         legend.text = element_text (size =13),
         legend.background = element_rect(fill="transparent"),
         legend.position = c(0.9, 0.28)) + 
  geom_hline(yintercept=0, linetype ="dashed", size =1, colour = "red") -> fig1d; fig1d
 
# Ensamble of all panels with an extra map in canva