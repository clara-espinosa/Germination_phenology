library(vegan); 
library(tidyverse); library(ggrepel); library(cowplot);library(ggpubr);
library (binom);library (ggsignif);library (rstatix); library (stringr);
library(patchwork);library (RColorBrewer);library(scales)
detach(package:plyr)
theme_set(theme_cowplot(font_size = 10)) 
Sys.setlocale("LC_ALL","English")
##### individual species curves #####

### seeds germ + viables/species and incubator 
read.csv("data/clean data.csv", sep = ";") %>%
  group_by(species, code, incubator, petridish) %>% # grouping variables
  summarise(total_germ = sum(germinated)) %>% # sum total germination
  merge(viables) %>% # merge viables object created above to filter raw data
  merge(species) %>% # join species header data
  group_by (species, incubator) %>% # new grouping variables
  summarise (viable = sum (viable))-> viables_sp # viables per species and incubator

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
  labs(title= "Fellfield cumulative germination curves", x = "Time ", y = "Germination proportion", tag="(a)") +
  theme_classic(base_size = 16) +
  theme(plot.margin = unit(c(0.25, 0,0,0.1), "cm"))+
  theme (plot.title = element_text (face = "bold",size = 16, margin=margin(0,0,0,0)), #,hjust = 0.5
         plot.tag.position = c(0.015,1),
         axis.title.y = element_text (size=14),
         axis.text.y = element_text (size = 12),
         axis.title.x = element_blank(), 
         axis.text.x= element_text (size = 12, color = "black", angle= 20),
         strip.text = element_text( size = 14, hjust = 0),
         strip.background = element_blank(), 
         panel.background = element_blank(), #element_rect(color = "black", fill = NULL), 
         #panel.grid = element_blank(),
         legend.title = element_text (size =14),
         legend.text = element_text (size =14),
         legend.position = "none", # legend.position = c(0.85, 0.5),
         legend.box.background = element_rect(color = "black", size = 2))-> fig2b;fig2b


# snowbed
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
  filter(incubator == "Snowbed") %>% ### change species name
  ggplot(aes(date, germination, group = species, color = species, fill = species)) +
  geom_line(size = 1.2) +
  facet_wrap(~ community, scales = "free_x", ncol = 2) +
  scale_y_continuous (limits = c(0,1), breaks = seq (0, 1, by= 0.25)) +
  scale_x_datetime(date_breaks = "2 month", date_labels = "%b %y")+
  scale_color_manual(values = Scolors) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title= "Snowbed cumulative germination curves", x = "Time ", y = "Germination proportion", tag = "(b)") +
  theme_classic(base_size = 16) +
  theme(plot.margin = unit(c(0.25, 0,0,0.1), "cm"))+
  theme (plot.title = element_text (face = "bold",size = 16, margin=margin(0,0,0,0)), #,hjust = 0.5
         plot.tag.position = c(0.015,1),
         axis.title.y = element_text (size=14),
         axis.text.y = element_text (size = 12),
         axis.title.x = element_blank(), 
         axis.text.x= element_text (size = 12, color = "black", angle= 20),
         strip.text = element_text( size = 14, hjust = 0),
         strip.background = element_blank(), 
         panel.background = element_blank(), #element_rect(color = "black", fill = NULL), 
         #panel.grid = element_blank(),
         legend.title = element_text (size =14),
         legend.text = element_text (size =14),
         legend.position = "none", # legend.position = c(0.85, 0.5),
         legend.box.background = element_rect(color = "black", size = 2))-> fig2c;fig2c

# partially arrenged with cowplot, final formatting in canva to add fig2a bottom
x11() # rearrange size 
fig2b/fig2c + plot_layout(heights = c(1,1))-> FigS2;FigS2
ggsave(filename = "Fig S2.png", plot =FigS2, path = "results/Supporting information/", 
       width = 18, height = 18, unit = "cm", device = "png", dpi = 600)
