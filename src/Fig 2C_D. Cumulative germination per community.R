library(FD); library(vegan); library(FactoMineR); library(emmeans);
library(tidyverse); library(ggrepel); library(cowplot);library(ggpubr);
library (binom);library (ggsignif);library (rstatix); library (stringr);
library(plyr);library(patchwork)
detach(package:plyr)
theme_set(theme_cowplot(font_size = 10)) 
Sys.setlocale("LC_ALL","English")

#### cumulative germination ####
# tidyverse modification to have the accumulated germination along the whole experiment + ggplot
data.frame(community = c("Mediterranean","Mediterranean","Mediterranean","Temperate", "Temperate"),
           incubator  = c("Fellfield", "Fellfield","Snowbed", "Fellfield", "Fellfield"),
           date = c("2021-11-12", "2022-04-03", "2021-11-12", "2021-11-12", "2022-04-03"),
           germinated = c(0,0,0,0,0)) -> data_vis
data_vis %>%
  mutate(date = as.POSIXct(date)) -> data_vis
str(data_vis)
x11()
#temperate
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(date = as.POSIXct(date))%>%
  merge (species) %>%
  group_by(community, incubator, date) %>%
  bind_rows (data_vis)%>%
  summarise(germinated = sum(germinated)) %>%
  mutate(germinated = cumsum(germinated)) %>%
  merge(viables1_community) %>%
  mutate(germination = germinated/viable) %>%
  filter (community == "Temperate") %>%
  ggplot(aes(date, germination, color = incubator, fill = incubator)) +
  geom_line(size = 2.5) +
  scale_color_manual (name= "Incubator", values = c ("Fellfield"= "chocolate2", "Snowbed" ="deepskyblue3")) +
  scale_y_continuous (limits = c(0,1), breaks = seq (0, 1, by= 0.25)) +
  scale_x_datetime(date_breaks = "2 month", date_labels = "%b-%Y")+
  labs(title = "D) Temperate (N = 38)", x = "Date", y = "Germination proportion") +
  theme_classic(base_size = 16) +
  theme (plot.title = element_text ( size = 32), #hjust = 0.5,
         axis.title.y = element_text (size=28),
         axis.text.y = element_text (size = 22),
         axis.title.x = element_blank(), 
         axis.text.x= element_blank(),
         plot.margin = margin(r= 40),
         legend.title = element_text (size =24),
         legend.text = element_text (size =24),
         legend.position = c(0.9, 0.15)) -> fig2d;fig2d

#mediterranean
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(date = as.POSIXct(date))%>%
  merge (species) %>%
  group_by(community, incubator, date) %>%
  bind_rows (data_vis)%>%
  summarise(germinated = sum(germinated)) %>%
  mutate(germinated = cumsum(germinated)) %>%
  merge(viables1_community) %>%
  mutate(germination = germinated/viable) %>%
  filter (community == "Mediterranean") %>%
  ggplot(aes(date, germination, color = incubator, fill = incubator)) +
  geom_line(size = 2.5) +
  scale_color_manual (name= "Incubator", values = c ("Fellfield"= "chocolate2", "Snowbed" ="deepskyblue3")) +
  scale_y_continuous (limits = c(0,1), breaks = seq (0, 1, by= 0.25)) +
  scale_x_datetime(date_breaks = "2 month", date_labels = "%b-%Y")+
  labs(title = "E) Mediterranean (N = 21)", x = "Date", y = "Germination proportion") +
  theme_classic(base_size = 16) +
  theme (plot.title = element_text ( size = 32), #hjust = 0.5,
         axis.title.y = element_text (size=28),
         axis.text.y = element_text (size = 22),
         axis.title.x = element_blank(), 
         axis.text.x= element_text (size=24),
         plot.margin = margin(r= 40),
         legend.title = element_text (size =24),
         legend.text = element_text (size =24),
         legend.position = c(0.9, 0.15)) -> fig2e; fig2e




ggsave(ggarrange(fig2a, fig2b, fig2c, fig2d,  ncol =1, nrow= 4,common.legend = FALSE), file="fig2.png", width = 210, height = 297, units = "mm")
