library(FD); library(vegan); library(FactoMineR); library(emmeans);
library(tidyverse); library(ggrepel); library(cowplot);library(ggpubr);
library (binom);library (ggsignif);library (rstatix); library (stringr);
library(patchwork);library (RColorBrewer);library(scales);library(plyr)
theme_set(theme_cowplot(font_size = 10)) 
Sys.setlocale("LC_ALL","English")

#### Fig 2A cumulative germination ####
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
  ggplot(aes(date, germination, color = incubator, fill = incubator)) +
  facet_wrap(~community, ncol = 2) +
  geom_line(size = 2.5) +
  scale_color_manual (name= "Incubator", values = c ("Fellfield"= "chocolate2", "Snowbed" ="deepskyblue3")) +
  scale_y_continuous (limits = c(0,1), breaks = seq (0, 1, by= 0.25)) +
  scale_x_datetime(date_breaks = "2 month", date_labels = "%b %y")+
  labs(title = "Cumulative germination curves (all species)", x = "Date", y = "Germination proportion", tag = "A") +
  theme_classic(base_size = 16) +
  theme (plot.title = element_text (face = "bold",size = 20), #hjust = 0.5,
         plot.tag.position = c(0,1),
         axis.title.y = element_text (size=14),
         axis.text.y = element_text (size = 13),
         axis.title.x = element_blank(), 
         axis.text.x= element_text (size = 12, color = "black"),
         strip.text = element_text( size = 18, hjust = 0),
         strip.background = element_blank(), 
         panel.background = element_blank(), #element_rect(color = "black", fill = NULL), 
         legend.title = element_text (size =14),
         legend.text = element_text (size =14),
         legend.position = "none", # legend.position = c(0.85, 0.5),
         legend.box.background = element_rect(color = "black", size = 2)) -> fig2a;fig2a

#stacked bar bottom fig2a
detach(package:plyr)
temp %>%
  group_by(incubator, germination_period) %>%
  summarise(n= length(date))%>%
  mutate(incubator = fct_relevel(incubator, "Snowbed", "Fellfield")) %>%
  mutate(germination_period = fct_relevel(germination_period, "Autumn", "Winter","Spring", "Summer" )) %>%
  ggplot(aes(x= incubator, y = n, fill = germination_period)) +
  geom_bar(stat = "identity", width=0.99, color = "black", position = position_stack(reverse= TRUE)) +
  coord_flip()+
  #labs (tag = "C") + 
  scale_fill_manual (name= "Germination phenology periods", 
                     values =c("brown", "cadetblue3","chartreuse3", "darkgoldenrod1"),
                     guide = guide_legend (title.position = "top",direction = "horizontal")) +
  theme_minimal(base_size = 20) +
  theme (plot.title = element_text (size = 24), #hjust = 0.5,
         axis.title.y = element_blank (), 
         axis.text.y = element_text(size= 22, color = c( "deepskyblue3","chocolate2")),
         axis.title.x = element_blank (), 
         axis.text.x= element_blank (),
         plot.tag.position = c(0,1),
         legend.title = element_text(size =24),
         legend.text = element_text(size=22),
         legend.position = "bottom",
         legend.margin = margin(0, 0, 0, 0)) -> fig2a_b;fig2a_b

##### individual species curves #####

# create/extent colorpalette (https://www.datanovia.com/en/blog/easy-way-to-expand-color-palettes-in-r/)
nb.cols <- 59
Fcolors <- colorRampPalette(brewer.pal(9, "Oranges"))(nb.cols)
Scolors <- colorRampPalette(brewer.pal(9, "Blues"))(nb.cols)

#fellfield
x11()
read.csv("data/clean data.csv", sep = ";") %>%
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
  filter(incubator == "Fellfield") %>% ### change species name
  ggplot(aes(date, germination, group = species, color = species, fill = species)) +
  geom_line(size = 1.2) +
  facet_wrap(~ community, scales = "free_x", ncol = 2) +
  scale_y_continuous (limits = c(0,1), breaks = seq (0, 1, by= 0.25)) +
  scale_x_datetime(date_breaks = "2 month", date_labels = "%b %y")+
  scale_color_manual(values = Fcolors) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title= "Fellfield cumulative germination curves (individual species)", x = "Time ", y = "Germination proportion", tag="B") +
  theme_classic(base_size = 16) +
  theme (plot.title = element_text (face = "bold",size = 20), #,hjust = 0.5
         plot.tag.position = c(0,1),
         axis.title.y = element_text (size=14),
         axis.text.y = element_text (size = 13),
         axis.title.x = element_blank(), 
         axis.text.x= element_text (size = 12, color = "black"),
         strip.text = element_text( size = 18, hjust = 0),
         strip.background = element_blank(), 
         panel.background = element_blank(), #element_rect(color = "black", fill = NULL), 
         #panel.grid = element_blank(),
         legend.title = element_text (size =14),
         legend.text = element_text (size =14),
         legend.position = "none", # legend.position = c(0.85, 0.5),
         legend.box.background = element_rect(color = "black", size = 2))-> fig3b;fig3b


# snowbed
x11()
read.csv("data/clean data.csv", sep = ";") %>%
  #rbind(read.csv("data/extrapoint_visualization.csv", sep = ";") )%>%
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
  ggplot(aes(date, germination, group = species, color = species, fill = species)) +
  geom_line(size = 1.2) +
  facet_wrap(~ community, scales = "free_x", ncol = 2) +
  scale_y_continuous (limits = c(0,1), breaks = seq (0, 1, by= 0.25)) +
  scale_x_datetime(date_breaks = "2 month", date_labels = "%b %y")+
  scale_color_manual(values = Scolors) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title= "Snowbed cumulative germination curves (individual species)", x = "Time ", y = "Germination proportion", tag = "C") +
  theme_classic(base_size = 16) +
  theme (plot.title = element_text (face = "bold",size = 20), #,hjust = 0.5
         plot.tag.position = c(0,1),
         axis.title.y = element_text (size=14),
         axis.text.y = element_text (size = 13),
         axis.title.x = element_blank(), 
         axis.text.x= element_text (size = 12, color = "black"),
         strip.text = element_text( size = 18, hjust = 0),
         strip.background = element_blank(), 
         panel.background = element_blank(), #element_rect(color = "black", fill = NULL), 
         #panel.grid = element_blank(),
         legend.title = element_text (size =14),
         legend.text = element_text (size =14),
         legend.position = "none", # legend.position = c(0.85, 0.5),
         legend.box.background = element_rect(color = "black", size = 2))-> fig2c;fig2c

ggarrange( fig2a, fig2b, fig2c,ncol =1, nrow= 3,common.legend = FALSE)
