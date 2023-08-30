library(FD); library(vegan); library(FactoMineR); library(emmeans);
library(tidyverse); library(ggrepel); library(cowplot);library(ggpubr);
library (binom);library (ggsignif);library (rstatix); library (stringr);
library(plyr);library(patchwork)
theme_set(theme_cowplot(font_size = 10)) 
Sys.setlocale("LC_ALL","English")
## GERMINATION EXPERIMENTS FIG 2!!!
#### dataframe  AND VISUALIZATION for temperature programs ####
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
