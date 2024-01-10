library (tidyverse)
# bioclim study area
# visualization test for figure 1 (Study area localization and description)
read.csv("data/Bioclim_study_area.csv", sep= ";") %>%
  #mutate(community = fct_relevel (community, "Temperate","Mediterranean"))%>%
  ggplot()+
  geom_boxplot (aes(x= bio1, y= community, fill = community), color = "black", position = "dodge2")+
  coord_flip() +
  scale_fill_manual (name= "Community", values = c ("Mediterranean"= "darkgoldenrod1" , "Temperate" ="forestgreen")) +
  labs(title = "Mean annual air temperature", x= "Temperature (ºC)", y= "Community", tag= "B")+
  theme_classic(base_size = 20) +
  theme (plot.title = element_text (size = 23), #hjust = 0.5,
         axis.title.y = element_text (size=14),
         axis.text.y = element_text (size = 13),
         plot.tag.position = c(0,1),
         axis.title.x = element_blank(), 
         axis.text.x= element_text (size = 14, color = "black"),
         legend.position = "none")

read.csv("data/Bioclim_study_area.csv", sep= ";") %>%
  #mutate(community = fct_relevel (community, "Temperate","Mediterranean"))%>%
  ggplot()+
  geom_boxplot (aes(x= bio17, y= community, fill = community), color = "black", position = "dodge2")+
  coord_flip()+
  scale_fill_manual (name= "Community", values = c ("Mediterranean"= "darkgoldenrod1" , "Temperate" ="forestgreen")) +
  labs(title = "Summer precipitation", x= "Precipitation (kg m-2)", y= "Community")+
  theme_classic(base_size = 16) +
  theme (plot.title = element_text (size = 23), #hjust = 0.5,
         axis.title.y = element_text (size=14),
         axis.text.y = element_text (size = 13),
         axis.title.x = element_blank(), 
         axis.text.x= element_text (size = 14, color = "black"),
         legend.position = "none")
