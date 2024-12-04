library(vegan); 
library(tidyverse); library(ggrepel); library(cowplot);library(ggpubr);
library (binom);library (ggsignif);library (rstatix); library (stringr);
library(patchwork);library (RColorBrewer);library(scales)
detach(package:plyr)
theme_set(theme_cowplot(font_size = 10)) 
Sys.setlocale("LC_ALL","English")

#### Fig 2A cumulative germination ####
# tidyverse modification to have the accumulated germination along the whole experiment + ggplot
#both communities
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(date = as.POSIXct(date))%>%
  merge (species) %>%
  group_by(community, incubator, date) %>%
  summarise(germinated = sum(germinated)) %>%
  mutate(germinated = cumsum(germinated)) %>%
  merge(viables_community) %>%
  mutate(germination = germinated/viable) %>%
  mutate(community = fct_relevel (community, "Temperate","Mediterranean"))%>%
  mutate (community = recode (community, "Temperate" = "Temperate (N = 38)", "Mediterranean" = "Mediterranean (N = 21)"))%>%
  ggplot(aes(date, germination, color = incubator, fill = incubator)) +
  facet_wrap(~community, ncol = 2) +
  geom_line(linewidth = 2) +
  scale_color_manual (name= "Climate regime", values = c ("Fellfield"= "chocolate2", "Snowbed" ="deepskyblue3")) +
  scale_y_continuous (limits = c(0,1), breaks = seq (0, 1, by= 0.25)) +
  scale_x_datetime(date_breaks = "2 month", date_labels = "%b %y")+
  labs(title = "Cumulative germination curves (all species)", x = "Date", y = "Germination proportion") +
  theme_classic(base_size = 16) +
  theme(plot.margin = unit(c(0.25, 0,0,0.1), "cm"))+
  theme (plot.title = element_text (face = "bold",size = 16, margin=margin(0,0,0,0)), #hjust = 0.5,
         plot.tag.position = c(0.015,1),
         axis.title.y = element_text (size=14),
         axis.text.y = element_text (size = 12),
         axis.title.x = element_blank(), 
         axis.text.x= element_text (size = 12, color = "black", angle= 20),
         strip.text = element_text( size = 14, hjust = 0),
         strip.background = element_blank(), 
         panel.background = element_blank(), #element_rect(color = "black", fill = NULL), 
         legend.title = element_text (size =13),
         legend.text = element_text (size =12),
         legend.position = c(0.9, 0.25))-> fig2a;fig2a  #legend.position = "none", 
         #legend.box.background = element_rect(color = "black", size = 2)) 
ggsave(filename = "fig2.png", plot =fig2a , path = "results/Figures/", width = 18, heigh = 10, units= "cm",
       device = "png", dpi = 600)

#stacked bar bottom fig2
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
  theme_minimal(base_size = 10) +
  theme (plot.title = element_text (size = 14), #hjust = 0.5,
         axis.title.y = element_blank (), 
         axis.text.y = element_text(size= 12, color = c( "deepskyblue3","chocolate2")),
         axis.title.x = element_blank (), 
         axis.text.x= element_blank (),
         plot.tag.position = c(0,1),
         legend.key.size = unit(0.35, "cm"),
         legend.title = element_text(size =12),
         legend.text = element_text(size=10),
         legend.position = "bottom",
         legend.margin = margin(0, 0, 0, 0)) -> fig2a_b1;fig2a_b1

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
  theme_minimal(base_size = 10) +
  theme (plot.title = element_text (size = 14), #hjust = 0.5,
         axis.title.y = element_blank (), 
         axis.text.y = element_blank (),
         axis.title.x = element_blank (), 
         axis.text.x= element_blank (),
         plot.tag.position = c(0,1),
         legend.key.size = unit(0.35, "cm"),
         legend.title = element_text(size =12),
         legend.text = element_text(size=10),
         legend.position = "bottom",
         legend.margin = margin(0, 0, 0, 0)) -> fig2a_b2;fig2a_b2

# combine figures
fig2a_b1+fig2a_b2 + plot_layout(guides= "auto", axis_titles= "collect", widths=c(1, 1)) &theme (legend.position = "bottom")->fig2_bottom;fig2_bottom

# combine figures# combine figuresaxis_titles = 
library(patchwork)
fig2a /fig2_bottom + plot_layout(heights = c(3,0.4))-> fig2;fig2

ggsave(filename = "fig2.png", plot =fig2 , path = "results/Figures/",
       scale=1, width =180, height = 120, units = "mm", device = "png", dpi = 600)
