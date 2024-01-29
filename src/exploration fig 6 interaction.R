library(FD); library(vegan); library(FactoMineR); library(emmeans);
library(tidyverse); library(ggrepel); library(cowplot);library(ggpubr);
library (binom);library (ggsignif);library (rstatix); library (stringr);
library(plyr);library(patchwork)
Sys.setlocale("LC_ALL","English")

# germination timing for explaining interaction, points + lines####
x11()
mean_values%>%
  filter(!trait=="T50") %>%
  filter(!trait=="Environmental heat sum") %>%
  filter(!trait=="Total germination") %>%
  mutate(trait = recode_factor(trait, "Autumn germination" = "Autumn", 
                               "Winter germination"="Winter", 
                               "Spring germination" = "Spring", 
                               "Summer germination" = "Summer"))%>%
  mutate(grouping_lines = paste(community, incubator))%>%
  ggplot()+
  geom_line (aes(x= trait, y= mean, group = grouping_lines, color = incubator), linewidth = 1.3)+
  geom_point(aes(x= trait, y= mean, fill = incubator), shape = 21, size = 3) +
  #geom_jitter(height = 0.1)+
  #geom_errorbar(aes(trait, mean, ymin = lower, ymax = upper),color = "black", width = 0.2, size =1) +
  scale_fill_manual (name= "", values =c("Fellfield" = "chocolate2", "Snowbed"="deepskyblue3")) + #values = c ("Mediterranean"= "darkgoldenrod1" , "Temperate" ="forestgreen")
  scale_color_manual (name= "", values =c("chocolate2", "deepskyblue3")) +
  facet_wrap(~community)+
  #scale_y_continuous (limits = c(0,1), breaks = seq (0, 1, by= 0.25)) +
  labs(title= "Germination phenology",  x = "Germination period", y = "Relative germination") + #tag = "B",
  theme_classic(base_size =12) +
  theme(plot.title = element_text (size = 24),
        strip.text.x = element_text( size = 20),
        legend.position = "bottom",
        legend.margin = margin(0, 0, 0, 0),
        plot.tag.position = c(0,1),
        legend.text = element_text(size=14),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.text.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.y = element_text (size=14), 
        axis.title.x = element_blank())

mean_values%>%
  filter(!trait=="Autumn germination") %>%
  filter(!trait=="Winter germination") %>%
  filter(!trait=="Spring germination") %>%
  filter(!trait=="Summer germination") %>%
  mutate(grouping_bars = paste(community, incubator))%>%
  mutate(grouping_bars = fct_relevel (grouping_bars, "Temperate Fellfield","Temperate Snowbed",
                                      "Mediterranean Fellfield", "Mediterranean Snowbed"))%>%
  ggplot(aes(fill=incubator, y=mean, x=community))+
  geom_bar(position="dodge", stat="identity", color = "black")+
  facet_wrap(~trait, ncol =1, scales = "free")+
  scale_fill_manual (name= "Incubator", values =c("chocolate2", "deepskyblue3")) +
  theme_classic(base_size =12) +
  theme(plot.title = element_text (size = 22),
        strip.text.x = element_text( size = 18),# face = "bold",
        strip.text.y = element_text(size = 14, angle = 360),
        legend.position = "none",
        legend.margin = margin(0, 0, 0, 0),
        plot.tag.position = c(0,1),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.text.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.y = element_text (size=15), 
        axis.title.x = element_blank ()) 

# germination timing for explaining interaction, stacked bar#####
mean_values%>%
  filter(!trait=="T50") %>%
  filter(!trait=="Environmental heat sum") %>%
  filter(!trait=="Total germination") %>%
  mutate(trait = recode_factor(trait, "Autumn germination" = "Autumn", 
                               "Winter germination"="Winter", 
                               "Spring germination" = "Spring", 
                               "Summer germination" = "Summer"))%>%
  mutate(trait = fct_relevel (trait, "Summer","Spring","Winter", "Autumn"))%>%
  mutate(grouping_bars = paste(community, incubator))%>%
  mutate(grouping_bars = fct_relevel (grouping_bars, "Temperate Fellfield","Temperate Snowbed",
                                      "Mediterranean Fellfield", "Mediterranean Snowbed"))%>%
  ggplot(aes(fill=trait, y=mean, x=grouping_bars))+
  geom_bar(position="fill", stat= "identity")+
  scale_y_continuous (name= "Germination proportion", limits = c(0,1), breaks = seq (0, 1, by= 0.25)) +
  scale_x_discrete (labels = c("Temperate\nFellfield","Temperate\nSnowbed",
                               "Mediterranean\nFellfield","Mediterranean\nSnowbed"))+
  scale_fill_manual (name= "Germination per phenology periods", 
                     values =c("darkgoldenrod1","chartreuse3","cadetblue3", "brown"),
                     guide = guide_legend (title.position = "top",direction = "horizontal"))+
  theme_classic(base_size = 16) +
  theme (plot.title = element_text (size = 24), #hjust = 0.5,
         axis.title.y = element_text(size= 18),
         axis.text.y = element_text(size= 16),
         axis.title.x = element_blank (), 
         axis.text.x= element_text(size= 16, angle = -25, vjust= 0.5, hjust = 0.4, color= "black"), #
         plot.tag.position = c(0,1),
         legend.title = element_text(size =20),
         legend.text = element_text(size=16),
         legend.position = "top",
         legend.margin = margin(0, 0, 0, 0))
