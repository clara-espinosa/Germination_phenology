library(FD); library(vegan); library(FactoMineR); library(emmeans);
library(tidyverse); library(ggrepel); library(cowplot);library(ggpubr);
library (binom);library (ggsignif);library (rstatix); library (stringr);
library(patchwork)
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
  labs(title = "Cumulative germination curves", x = "Date", y = "Germination proportion") +
  theme_classic(base_size = 16) +
  theme (plot.title = element_text (face = "bold",hjust = 0.5,size = 24), #
         axis.title.y = element_text (size=16),
         axis.text.y = element_text (size = 14),
         axis.title.x = element_blank(), 
         axis.text.x= element_text (size = 13, color = "black"),
         strip.text = element_text( size = 20, hjust = 0),
         strip.background = element_blank(), 
         panel.background = element_blank(), #element_rect(color = "black", fill = NULL), 
         #panel.grid = element_blank(),
         legend.title = element_text (size =14),
         legend.text = element_text (size =14),
         legend.position = "bottom", # legend.position = c(0.85, 0.5),
         legend.box.background = element_rect(color = "black", size = 2)) -> fig2d;fig2d

#ggarrange (fig2a, fig2b, fig2c, fig2d,ncol =1, nrow= 4,common.legend = FALSE)

##### individual species curves #####
x11()
read.csv("data/clean data.csv", sep = ";") %>%
  rbind(data_vis)%>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  spread(date, germinated, fill = 0) %>% # wide format for dates, and fill Na with 0
  gather ("date", "germinated", 8: last_col() )%>% # back in long format frrom 8th colum to the last
  arrange (species, code, incubator, petridish, date)%>% # sort row observations this way
  mutate(date = as.POSIXct(date))%>%
  merge (species) %>%
  dplyr::group_by(community, species, incubator, date) %>%
  dplyr::summarise(germinated = sum(germinated)) %>%
  dplyr::mutate(germinated = cumsum(germinated)) %>%
  merge(viables_sp) %>%
  mutate(germination = germinated/viable)  %>%   
  mutate(community = fct_relevel (community, "Temperate","Mediterranean"))%>%
  mutate (community = recode (community, "Temperate" = "Temperate (N = 38)", "Mediterranean" = "Mediterranean (N = 21)"))%>%
  filter(incubator == "Snowbed") %>% ### change species name
  ggplot(aes(date, germination, group = species, color = incubator, fill = incubator)) +
  geom_line(size = 2) +
  scale_color_manual (name= "Incubator", values = c ("Snowbed" ="deepskyblue3")) + #,"Fellfield"= "chocolate2" 
  facet_wrap(~ community, scales = "free_x", ncol = 2) +
  scale_y_continuous (limits = c(0,1), breaks = seq (0, 1, by= 0.25)) +
  scale_x_datetime(date_breaks = "2 month", date_labels = "%b %y")+
  coord_cartesian(ylim = c(0, 1)) +
  labs(title= "Snowbed germination curves", x = "Time ", y = "Germination proportion") +
  theme_classic(base_size = 16) +
  theme (plot.title = element_text (face = "bold",hjust = 0.5,size = 24), #
         axis.title.y = element_text (size=16),
         axis.text.y = element_text (size = 14),
         axis.title.x = element_blank(), 
         axis.text.x= element_text (size = 13, color = "black"),
         strip.text = element_text( size = 20, hjust = 0),
         strip.background = element_blank(), 
         panel.background = element_blank(), #element_rect(color = "black", fill = NULL), 
         #panel.grid = element_blank(),
         legend.title = element_text (size =14),
         legend.text = element_text (size =14),
         legend.position = "none", # legend.position = c(0.85, 0.5),
         legend.box.background = element_rect(color = "black", size = 2))

ggarrange( fig2d, fig2e, ncol =1, nrow= 2,common.legend = TRUE)

ggsave(ggarrange( fig2d, fig2e, ncol =1, nrow= 4,common.legend = TRUE, file="fig2.png", width = 210, height = 297, units = "mm")
