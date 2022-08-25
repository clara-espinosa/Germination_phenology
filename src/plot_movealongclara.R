library(FD); library(vegan); library(FactoMineR); library(emmeans);
library(tidyverse); library(ggrepel); library(cowplot)
theme_set(theme_cowplot(font_size = 10)) 

## GERMINATION EXPERIMENTS

# Read and prepare the germination data

# dataframe for temperature programs
temp <- read.csv("data/date_temp.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y")) %>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  mutate(condition = as.factor(condition))
temp <- as.data.frame(temp)
str(temp)
summary(temp)
# visualization 1
# mulitple lines with shaded areas in between
ggplot(temp, aes ()) +
  geom_line (aes (x=time, y=Tmin, colour =condition), size =1) +
  geom_line (aes (x=time, y=Tmax, colour = condition), size =1) +
  geom_ribbon (data = temp, aes (x=time, ymin =Tmin, ymax=Tmax, colour = condition, fill = condition), alpha =0.5) +
  labs (title = "Move-along temperature regimes", y= "Temperature ºC", x = "Time (days)") + 
  theme (axis.title.y = element_text (size=12), 
         axis.title.x = element_text (size=12), 
         plot.title = element_text (size =14), 
         legend.title = element_text(size = 12)) +
  geom_vline(xintercept = 122, linetype = "dashed", size =1) +
  annotate (geom ="text", x= 155, y = 22, label ="winter begins", colour = "black", size = 3,  fontface ="bold") +
  geom_vline(xintercept = 241, linetype = "dashed", size =1, color = "red") +
  annotate (geom ="text", x= 265, y = 22, label ="fellfield\nwinter ends", colour = "red", size = 3,  fontface ="bold") +
  geom_vline(xintercept = 290, linetype = "dashed", size =1,color = "turquoise") +
  annotate (geom ="text", x= 315, y = 22, label ="snowbed\nwinter ends", colour = "turquoise", size = 3,  fontface ="bold") 
  
# visualization 2
# all area shaded below x temperature (PROBLEM: don't know why but all values are higher than what appear in data)
ggplot(temp, aes())+
geom_area(aes(time,  Tmean, colour = condition, fill = condition), size = 1, alpha = 0.4) +
  geom_vline(xintercept = 122, linetype = "dashed", size =1) +
  annotate (geom ="text", x= 155, y = 20, label ="winter begins", colour = "black", size = 4,  fontface ="bold") +
  geom_vline(xintercept = 241, linetype = "dashed", size =1, color = "red") +
  annotate (geom ="text", x= 260, y = 20, label ="fellfield", colour = "red", size = 4,  fontface ="bold") +
  geom_vline(xintercept = 290, linetype = "dashed", size =1,color = "turquoise") +
  annotate (geom ="text", x= 312, y = 20, label ="snowbed", colour = "turquoise", size = 4,  fontface ="bold")

# tidyverse transformation to account for the number of viable seeds per each specie and condition
# summing up petridishes and accesions/populations of the same species (not taking into account weekly germination)
read.csv("data/R long data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date)))) %>%
  group_by(species, condition, code, petridish) %>%
  filter(date == max(date)) %>%
  select(species, condition, code, petridish, viable) %>%
  group_by(species, condition) %>%
  summarise(viable = sum(viable)) -> viables
 
write.csv (viables,"results/viables.csv", row.names = FALSE )

# tidyverse modification to have the accumulated germination along the whole experiment + ggplot
read.csv("data/R long data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date)))) %>%
  group_by(species, condition, time) %>%
  summarise(germinated = sum(germinated)) %>%
  mutate(germinated = cumsum(germinated)) %>%
  merge(viables) %>%
  mutate(germination = germinated/viable) %>%
  filter(species == "Agrostis tileni") %>%
  ggplot(aes(time, germination, color = condition, fill = condition)) +
  geom_line(size = 1) +
  scale_color_manual (name= "Condition", values = c ("Fellfield"= "red", "Snowbed" ="turquoise")) +
  facet_wrap( ~ species, scales = "free_x", ncol = 2) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Time (days)", y = "Germination proportion") +
  theme(strip.text = element_text(face = "italic", size = 14), 
        legend.title = element_text(size=14), 
        legend.text = element_text(size=10)) + 
  geom_vline(xintercept = 70, linetype = "dashed", size= 0.75) +
  geom_vline(xintercept = 126, linetype = "dashed", size= 0.75, color = "red") +
  geom_vline(xintercept = 168, linetype = "dashed", size= 0.75, color = "turquoise") 
  
# Save the plots

ggsave(filename = "results/FigAgrostistileni.png", Agrostistileni, path = NULL, 
       scale = 1, width = 180, height = 180, units = "mm", dpi = 600)

# option to aggrupate species per mountain??


# Create the plots, divided for 4 pages

dfg %>%
  filter(species %in% unique(dfg$species)[4:8]) %>% # 12 species
  fg -> FigGerminationA

dfg %>%
  filter(species %in% unique(dfg$species)[13:26]) %>%
  fg -> FigGerminationB

dfg %>%
  filter(species %in% unique(dfg$species)[27:40]) %>%
  fg -> FigGerminationC

dfg %>%
  filter(species %in% unique(dfg$species)[41:55]) %>%
  fg -> FigGerminationD

