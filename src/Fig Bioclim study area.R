library (tidyverse)
# bioclim study area
# visualization test for figure 1 (Study area localization and description)
read.csv("data/Bioclim_study_area.csv", sep= ";") %>%
  #mutate(community = fct_relevel (community, "Temperate","Mediterranean"))%>%
  ggplot()+
  geom_boxplot (aes(x= bio1, y= community, fill = community), color = "black", position = "dodge2")+
  coord_flip() +
  scale_fill_manual (name= "Community", values = c ("Mediterranean"= "firebrick3", "Temperate" ="cyan3")) +
  labs(title = "Mean annual air temperature", x= "Temperature (ºC)", y= "Community")+
  theme_classic(base_size = 16) +
  theme (plot.title = element_text (face = "bold",size = 18), #hjust = 0.5,
         axis.title.y = element_text (size=14),
         axis.text.y = element_text (size = 13),
         axis.title.x = element_blank(), 
         axis.text.x= element_text (size = 12, color = "black"),
         legend.position = "none")

read.csv("data/Bioclim_study_area.csv", sep= ";") %>%
  mutate(community = fct_relevel (community, "Temperate","Mediterranean"))%>%
  ggplot()+
  geom_boxplot (aes(x= bio14, y= community,  fill = community),color = "black", position = "dodge2")+
  coord_flip()+
  labs(title = "Precipitation of the driest month", x= "Precipitation (kg m-2)", y= "Community")+
  theme_classic(base_size = 16) +
  theme (plot.title = element_text (face = "bold",size = 16), #hjust = 0.5,
         axis.title.y = element_text (size=14),
         axis.text.y = element_text (size = 13),
         axis.title.x = element_blank(), 
         axis.text.x= element_text (size = 12, color = "black"),
         legend.position = "none")

read.csv("data/Bioclim_study_area.csv", sep= ";") %>%
  #mutate(community = fct_relevel (community, "Temperate","Mediterranean"))%>%
  ggplot()+
  geom_boxplot (aes(x= bio17, y= community, fill = community), color = "black", position = "dodge2")+
  coord_flip()+
  scale_fill_manual (name= "Community", values = c ("Mediterranean"= "firebrick3", "Temperate" ="cyan3")) +
  labs(title = "Summer precipitation", x= "Precipitation (kg m-2)", y= "Community")+
  theme_classic(base_size = 16) +
  theme (plot.title = element_text (face = "bold",size = 18), #hjust = 0.5,
         axis.title.y = element_text (size=14),
         axis.text.y = element_text (size = 13),
         axis.title.x = element_blank(), 
         axis.text.x= element_text (size = 12, color = "black"),
         legend.position = "none")
