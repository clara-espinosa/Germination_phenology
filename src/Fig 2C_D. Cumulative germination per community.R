library(FD); library(vegan); library(FactoMineR); library(emmeans);
library(tidyverse); library(ggrepel); library(cowplot);library(ggpubr);
library (binom);library (ggsignif);library (rstatix); library (stringr);
library(plyr);library(patchwork)
detach(package:plyr)
theme_set(theme_cowplot(font_size = 10)) 
Sys.setlocale("LC_ALL","English")

#### cumulative germination ####
# tidyverse modification to have the accumulated germination along the whole experiment + ggplot
data.frame(community = c("Mediterranean","Mediterranean","Mediterranean","Temperate", "Temperate" ),
           incubator  = c("Fellfield", "Fellfield","Snowbed", "Fellfield", "Fellfield"),
           date = c("2021-11-12", "2022-04-03", "2021-11-12", "2021-11-12", "2022-04-03"),
           germinated = c(0,0,0,0,0)) -> data_vis
data_vis %>%
  mutate(date = as.POSIXct(date)) -> data_vis
str(data_vis)
x11()

#both communities
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(date = as.POSIXct(date))%>%
  merge (species) %>%
  group_by(community, incubator, date) %>%
  bind_rows (data_vis)%>%
  summarise(germinated = sum(germinated)) %>%
  mutate(germinated = cumsum(germinated)) %>%
  merge(viables_community) %>%
  mutate(germination = germinated/viable) %>%
  mutate(community = fct_relevel (community, "Temperate","Mediterranean"))%>%
  mutate (community = recode (community, "Temperate" = "Temperate (N = 38)", "Mediterranean" = "Mediterranean (N = 21)"))%>%
  #filter (community == "Temperate") %>%
  ggplot(aes(date, germination, color = incubator, fill = incubator)) +
  facet_wrap(~community, ncol = 2) +
  geom_line(size = 2.5) +
  scale_color_manual (name= "Incubator", values = c ("Fellfield"= "chocolate2", "Snowbed" ="deepskyblue3")) +
  scale_y_continuous (limits = c(0,1), breaks = seq (0, 1, by= 0.25)) +
  scale_x_datetime(date_breaks = "2 month", date_labels = "%b %y")+
  labs(x = "Date", y = "Germination proportion") +
  theme_classic(base_size = 16) +
  theme (plot.title = element_text (size = 24), #hjust = 0.5,
         axis.title.y = element_text (size=20),
         axis.text.y = element_text (size = 14),
         axis.title.x = element_blank(), 
         axis.text.x= element_text (size = 13.5, color = "black"),
         strip.text = element_text(face = "bold", size = 20, hjust = 0),
         strip.background = element_blank(), 
         panel.background = element_blank(), #element_rect(color = "black", fill = NULL), 
         panel.grid = element_blank(),
         legend.title = element_text (size =14),
         legend.text = element_text (size =14),
         legend.position = "bottom", # legend.position = c(0.85, 0.5),
         legend.box.background = element_rect(color = "black", size = 2)) -> fig2d;fig2d

#ggarrange (fig2a, fig2b, fig2c, fig2d,ncol =1, nrow= 4,common.legend = FALSE)

#inidividual community ex: mediterranean
read.csv("data/clean data.csv", sep = ";") %>%
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


ggarrange( fig2d, fig2e, ncol =1, nrow= 2,common.legend = TRUE)

ggsave(ggarrange( fig2d, fig2e, ncol =1, nrow= 4,common.legend = TRUE, file="fig2.png", width = 210, height = 297, units = "mm")
