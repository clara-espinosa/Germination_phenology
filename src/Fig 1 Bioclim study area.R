library(vegan); library(FactoMineR); library(emmeans);
library(tidyverse); library(ggrepel); library(cowplot);library(ggpubr);
library (binom);library (ggsignif);library (rstatix); library (stringr);
library(plyr);library(patchwork);library(scales)
theme_set(theme_cowplot(font_size = 10)) 
Sys.setlocale("LC_ALL","English")

## FIG 2B bioclim study area ####
# visualization test for figure 1 (Study area localization and description)
read.csv("data/Bioclim_study_area.csv", sep= ",") %>%
  #mutate(community = fct_relevel (community, "Temperate","Mediterranean"))%>%
  ggplot()+
  geom_boxplot (aes(x= bio1, y= community, fill = community), color = "black", position = "dodge2")+
  coord_flip() +
  scale_fill_manual (name= "Community", values = c ("Mediterranean"= "darkgoldenrod1" , "Temperate" ="forestgreen")) +
  labs(title = "Mean annual air temperature", x= "Temperature (ºC)", y= "Community", tag= "(b)")+
  theme_classic(base_size = 12) +
  theme (plot.title = element_text (size = 14), #hjust = 0.5,
         axis.title.y = element_text (size=12),
         axis.text.y = element_text (size = 12),
         plot.tag.position = c(0.015,1),
         axis.title.x = element_blank(), 
         axis.text.x= element_text (size = 10, color = "black"),
         legend.position = "none") -> fig1b_1;fig1b_1

ggsave(filename = "fig1b_1.png", plot =fig1b_1, path = "results/Figures/", 
        device = "png", dpi = 600)

read.csv("data/Bioclim_study_area.csv", sep= ",") %>%
  #mutate(community = fct_relevel (community, "Temperate","Mediterranean"))%>%
  ggplot()+
  geom_boxplot (aes(x= bio17, y= community, fill = community), color = "black", position = "dodge2")+
  coord_flip()+
  scale_fill_manual (name= "Community", values = c ("Mediterranean"= "darkgoldenrod1" , "Temperate" ="forestgreen")) +
  labs(title = "Summer precipitation", x= "Precipitation (mm)", y= "Community")+
  theme_classic(base_size = 12) +
  theme (plot.title = element_text (size = 14), #hjust = 0.5,
         axis.title.y = element_text (size=12),
         axis.text.y = element_text (size = 10),
         axis.title.x = element_blank(), 
         axis.text.x= element_text (size = 12, color = "black"),
         legend.position = "none")-> fig1b_2;fig1b_2
ggsave(filename = "fig1b_2.png", plot =fig1b_2, path = "results/Figures/", 
       device = "png", dpi = 600)

# combine plots
library(patchwork)
fig1b_1 +fig1b_2 + plot_layout()-> fig1b;fig1b

ggsave(filename = "fig1b.png", plot =fig1b, path = "results/Figures/", 
       device = "png", dpi = 600)


## Fig 2C graph to compare with temperature programs ####
x11()
read.csv("data/weekly_picos_graph.csv", sep = ",") %>%
  #mutate(order2 = paste(month, week))%>%
  ggplot(aes ()) +
  geom_line (aes (x=order, y=N, colour = Site), linewidth =1.25) + #order for sort as experiment
  geom_line (aes (x=order, y=X, colour = Site), linewidth =1.25) +
  geom_ribbon (aes (x=order, ymin =N, ymax=X, colour = Site, fill = Site), alpha =0.5) +
  scale_color_manual (name = "Sites", values = c ("Fellfield site"= "chocolate2" , "Snowbed site" = "deepskyblue3"),
                      labels = c("Fellfield", "Snowbed")) +
  scale_fill_manual (name = "Sites", values = c ("Fellfield site"= "chocolate2" , "Snowbed site" = "deepskyblue3"),
                     labels = c("Fellfield", "Snowbed")) +
  scale_y_continuous (limits = c(-3,25), breaks = seq (0, 25, by= 5)) +
  scale_x_continuous(breaks = seq(0,60, by=8), labels = c("", "Sep", "Nov", "Jan", "Mar", "May", "July", "Sep"))+
  labs (title = "Average soil temperatures 2008 - 2019", y= "Temperature (ºC)", x = "Time (weeks)", tag = "(c)") + #
  theme_classic(base_size = 8) +
  theme (plot.title = element_text ( size = 12), #hjust = 0.5,
         axis.title.y = element_text (size=10), 
         axis.text.y = element_text (size = 8),
         axis.title.x = element_blank(), 
         axis.text.x= element_text (size = 8, color = "black"),
         plot.tag.position = c(0.015,1),
         legend.key.size = unit(0.45, "cm"),
         legend.title = element_text (size =10),
         legend.text = element_text (size =8),
         legend.background = element_rect(fill="transparent"),
         legend.position = c(0.9, 0.3)) + # legend.position = "top") + 
  geom_hline(yintercept=0, linetype ="dashed", size =1, colour = "red") -> fig1c; fig1c

ggsave(filename = "fig1c", plot =fig1c, path = "results/Figures/", 
       device = "png", dpi = 600)

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
  scale_fill_manual (name= "Climate regime", values =c("chocolate2", "deepskyblue3")) +
  scale_color_manual (name= "Climate regime", values =c("chocolate2", "deepskyblue3")) +
  scale_y_continuous (limits = c(-3,25), breaks = seq (0, 25, by= 5)) +
  scale_x_datetime(date_breaks = "2 month", date_labels = "%b %y")+
  labs (title = "Climate regimes temperature programs", y= "Temperature (ºC)", x = "Date", tag= "(d)") + #
  theme_classic(base_size = 8) +
  theme (plot.title = element_text ( size = 12), #hjust = 0.5,
         axis.title.y = element_text (size=10), 
         axis.title.x = element_blank(), 
         axis.text.x= element_text (size = 8, color = "black"),
         axis.text.y = element_text(size = 8),
         plot.tag.position = c(0.015,1),
         legend.key.size = unit(0.45, "cm"),
         legend.title = element_text (size =10),
         legend.text = element_text (size =8),
         legend.background = element_rect(fill="transparent"),
         legend.position = c(0.89, 0.29)) + 
  geom_hline(yintercept=0, linetype ="dashed", size =1, colour = "red") -> fig1d; fig1d
 

ggsave(filename = "fig1d", plot =fig1d, path = "results/Figures/", 
       device = "png", dpi = 600)

# combine plots
# combine plots
library(patchwork)
x11()
fig1c +fig1d + plot_layout()-> fig1cd;fig1cd

ggsave(filename = "fig1cd.png", plot =fig1cd, path = "results/Figures/", 
       device = "png", dpi = 600)
# Ensamble of all panels with an extra map in powerpoint