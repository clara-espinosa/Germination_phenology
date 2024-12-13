library(tidyverse);library (rstatix);library (stringr); library (vegan);library(cowplot)

# significances differences between germination 

read.csv("results/Field germination segment significances.csv", sep = ";")%>%
  convert_as_factor(species, microhabitat_buried, retrieval_season)%>%
  mutate(x= as.numeric(x),
         xend= as.numeric(xend),
         y= as.numeric(y),
         yend= as.numeric(yend))->sig.segment
str(sig.segment)

read.csv("results/Field germination label significances.csv", sep = ";")%>%
  convert_as_factor(species, microhabitat_buried, retrieval_season)%>%
  mutate(x= as.numeric(x),
         y= as.numeric(y))->sig.label
str(sig.label)

#bar plot for visualization
read.csv ("data/field_germination.csv") %>%
  select(retrieval_season, community, microhabitat_buried, species, seeds_initial, bag1:bag3)%>%
  convert_as_factor(retrieval_season, community, microhabitat_buried, species) %>%
  mutate(species = fct_relevel (species, "Armeria duriaei", "Dianthus langeanus", "Plantago holosteum",
                                "Luzula caespitosa", "Phyteuma hemisphaericum", "Silene ciliata", 
                                "Androsace villosa",  "Carex sempervirens","Gypsophila repens", 
                                "Armeria cantabrica","Festuca glacialis","Jasione cavanillesii"))%>%
  gather(bag, germ, bag1:bag3)%>%
  group_by(retrieval_season, microhabitat_buried, species)%>%
  dplyr::summarise(germ=sum(germ))%>%
  spread(retrieval_season, germ)%>%
  mutate(Autumn_23 = Autumn_23-Spring_23)%>%
  mutate(Autumn_23 = ifelse(Autumn_23>0,Autumn_23, 0))%>%
  gather(retrieval_season, field_germ, Spring_23:Autumn_23)%>%
  mutate(retrieval_season = recode_factor (retrieval_season, "Spring_23" = "Spring","Autumn_23"= "2nd Autumn"))%>%
  #mutate(retrieval_season = fct_relevel (retrieval_season, "Spring","Autumn"))%>%
  #filter(species=="Armeria duriaei")%>%
  ggplot(aes(retrieval_season, field_germ, condition=microhabitat_buried, fill = microhabitat_buried))+
  geom_bar(width = 0.7, position=position_dodge(width = 0.8), stat="identity", color="black")+
  scale_fill_manual (name="Microhabitat",values = c ("chocolate2", "deepskyblue3")) + #
  labs (y= "Number of field germinated seeds", x= "Retrieval season")+ #title= "Field germination",
  scale_y_continuous (limits = c(0,60), breaks = seq (0, 60, by= 15)) +
  facet_wrap(~species, ncol=3)+
  geom_segment (data=sig.segment, aes(x= x,xend =xend,  y = y, yend= yend), color = "black", linewidth = 1, show.legend = F)+
  geom_text (data= sig.label, aes(x= x, y= y, label = label), size= 4.5)+
  labs(subtitle = "")+
  theme_classic(base_size = 16) +
  theme(plot.title = element_text (size = 22),
        plot.subtitle= element_text(size= 16, hjust=0.5),
        plot.margin = margin(20,2, 2, 2),
        strip.text.x = element_text( size = 12, face = "italic"),# face = "bold",
        strip.text.y = element_text(size = 14, angle = 360),
        legend.position = "bottom", #bottom
        legend.box.margin= margin(-10,-10,0,-10),
        plot.tag.position = c(0,1),
        panel.spacing.y=unit(1.5,"lines"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.text.x = element_text(size = 13, color = "black"), #, angle = 10
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.y = element_text (size=15), 
        axis.title.x = element_text (size=15)) -> fig5;fig5

library(cowplot)
ggdraw(fig5) + draw_plot_label(x = c(0.445,0.465), y = c(0.97, 0.55), 
                            label =c("Mediterranean","Temperate"),
                            hjust = 0, size = 15)-> fig5;fig5

ggsave(filename = "fig5.png", plot =fig5, path = "results/Figures/", 
       scale=1, width =180, height = 210, units = "mm",device = "png", dpi = 600)
