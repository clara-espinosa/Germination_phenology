library(FD); library(vegan); library(FactoMineR); library(emmeans);
library(tidyverse); library(ggrepel); library(cowplot);library(ggpubr);
library (binom);library (ggsignif);library (rstatix); library (stringr);
library(plyr);library(patchwork);library(scales)
theme_set(theme_cowplot(font_size = 10)) 
Sys.setlocale("LC_ALL","English")
## GERMINATION EXPERIMENTS FIG 2!!!
#### dataframe  AND VISUALIZATION for temperature programs ####
read.csv("data/date_temp.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y")) %>%
  mutate(date = as.POSIXct(date))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  mutate(incubator = as.factor(incubator), 
         germination_period = as.factor (germination_period)) %>%
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
  scale_y_continuous (limits = c(-3,25), breaks = seq (0, 25, by= 5)) +
  scale_x_datetime(date_breaks = "2 month", date_labels = "%b %y")+
  labs (title = "Experimental temperature programs", y= "Temperature ºC", x = "Date") + #, tag= "B"
  theme_classic(base_size = 20) +
  theme (plot.title = element_text ( size = 24), #hjust = 0.5,
         axis.title.y = element_text (size=20), 
         axis.title.x = element_blank(), 
         axis.text.x= element_text (size = 18, color = "black"),
         axis.text.y = element_text(size = 18),
         plot.tag.position = c(0,1),
         legend.title = element_text (size =14),
         legend.text = element_text (size =13),
         legend.background = element_rect(fill="transparent"),
         legend.position = c(0.9, 0.28)) + 
  geom_hline(yintercept=0, linetype ="dashed", size =1, colour = "red") -> fig2b; fig2b

#stacked bar 
detach(package:plyr)
temp %>%
  group_by(incubator, germination_period) %>%
  summarise(n= length(date))%>%
  mutate(incubator = fct_relevel(incubator, "Snowbed", "Fellfield")) %>%
  mutate(germination_period = fct_relevel(germination_period, "Autumn", "Winter","Spring", "Summer" )) %>%
  ggplot(aes(x= incubator, y = n, fill = germination_period)) +
  geom_bar(stat = "identity", width=0.99, color = "black", position = position_stack(reverse= TRUE)) +
  coord_flip()+
  labs (tag = "C") + 
  scale_fill_manual (name= "Germination phenology periods", values =c("brown", "cadetblue3","chartreuse3", "darkgoldenrod1"),
                     guide = guide_legend (title.position = "top",direction = "horizontal")) +
  theme_minimal(base_size = 20) +
  theme (plot.title = element_text (size = 24), #hjust = 0.5,
         axis.title.y = element_blank (), 
         axis.text.y = element_text(size= 18, color = c( "deepskyblue3","chocolate2")),
         axis.title.x = element_blank (), 
         axis.text.x= element_blank (),
         plot.tag.position = c(0,1),
         legend.position = "bottom",
         legend.margin = margin(0, 0, 0, 0)) -> fig2c;fig2c
         #legend.title = element_text (size =14),
         #legend.text = element_text (size =14),
          
  
         

  #geom_vline(xintercept = as.POSIXct(as.Date("2021-11-12")), linetype = "dashed", size =1.25) +
  #annotate (geom ="text", x= as.POSIXct(as.Date("2021-10-17")), y = 28, label ="Autumn \n germination", size = 8, fontface ="bold") +
  #geom_vline(xintercept = as.POSIXct(as.Date("2022-06-15")), linetype = "dashed", size =1.25) +
  #annotate (geom ="text", x= as.POSIXct(as.Date("2022-05-20")), y = 28, label ="Spring\n germination", size = 8, fontface ="bold") +
  #geom_vline(xintercept = as.POSIXct(as.Date("2022-09-19")), linetype = "dashed", size =1.25) +
  #annotate (geom ="text", x= as.POSIXct(as.Date("2022-08-25")), y = 28, label ="Summer\n germination", size = 8, fontface ="bold") +
  #geom_segment (aes(x=as.POSIXct(as.Date("2021-11-12")), y =-2.8, xend=as.POSIXct(as.Date("2022-04-04")), yend =-2.8), color = "chocolate2", size = 2) +
  #geom_segment (aes(x=as.POSIXct(as.Date("2021-11-12")), y = -3.8, xend =as.POSIXct(as.Date("2022-05-26")), yend =-3.8), color = "deepskyblue3", size = 2) +
  #annotate (geom ="text", x= as.POSIXct(as.Date("2022-01-15")), y = -5, label ="Winter conditions", colour = "black", size = 8.5, fontface ="bold") 
