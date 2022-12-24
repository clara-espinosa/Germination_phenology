library(FD); library(vegan); library(FactoMineR); library(emmeans);
library(tidyverse); library(ggrepel); library(cowplot)
theme_set(theme_cowplot(font_size = 10)) 

## GERMINATION EXPERIMENTS

# Read and prepare the germination data

#### dataframe  AND VISUALIZATION for temperature programs ####
temp <- read.csv("data/date_temp.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y")) %>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  mutate(incubator = as.factor(incubator))
temp <- as.data.frame(temp)
str(temp)
summary(temp)
# visualization 1
# mulitple lines with shaded areas in between
x11()
ggplot(temp, aes ()) +
  geom_line (aes (x=time, y=Tmin, colour =incubator), size =1) +
  geom_line (aes (x=time, y=Tmax, colour = incubator), size =1) +
  geom_ribbon (data = temp, aes (x=time, ymin =Tmin, ymax=Tmax, colour = incubator, fill = incubator), alpha =0.5) +
  scale_fill_manual (values =c("chocolate2", "deepskyblue3")) +
  scale_color_manual (values =c("chocolate2", "deepskyblue3")) +
  labs (title = "Move-along temperature regimes", y= "Temperature ºC", x = "Time (days)") + 
  theme (axis.title.y = element_text (size=14), 
         axis.title.x = element_text (size=14), 
         plot.title = element_text (size =16), 
         legend.title = element_text(size = 14),
         legend.text = element_text (size =12)) +
  geom_vline(xintercept = 122, linetype = "dashed", size =1) +
  annotate (geom ="text", x= 155, y = 22, label ="winter begins", colour = "black", size = 3.5,  fontface ="bold") +
  geom_vline(xintercept = 240, linetype = "dashed", size =1, color = "chocolate2") +
  annotate (geom ="text", x= 265, y = 22, label ="fellfield\nwinter ends", colour = "chocolate2", size = 3.5,  fontface ="bold") +
  geom_vline(xintercept = 291, linetype = "dashed", size =1,color = "deepskyblue3") +
  annotate (geom ="text", x= 317, y = 22, label ="snowbed\nwinter ends", colour = "deepskyblue3", size = 3.5,  fontface ="bold") +
  geom_hline(yintercept=0, linetype ="dashed", size =0.5, colour = "red")
  
#### main data transformation and visualization ####
#### SPECIES GERMINATION CURVES ####
# tidyverse transformation to account for the number of viable seeds per each specie and incubator
# summing up petridishes and accesions/populations of the same species (not taking into account weekly germination)
read.csv("data/R long data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date)))) %>%
  group_by(species, incubator, code, petridish) %>%
  filter(date == max(date)) %>%
  select(species, incubator, code, petridish, viable) %>%
  group_by(species, incubator) %>%
  summarise(viable = sum(viable)) -> viables
# write.csv (viables,"results/viables.csv", row.names = FALSE )
# tidyverse modification to have the accumulated germination along the whole experiment + ggplot
x11()
read.csv("data/R long data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date)))) %>%
  group_by(species, incubator, time) %>%
  summarise(germinated = sum(germinated)) %>%
  mutate(germinated = cumsum(germinated)) %>%
  merge(viables) %>%
  mutate(germination = germinated/viable) %>%
  filter(species == "Thymus praecox") %>% ### change species name
  ggplot(aes(time, germination, color = incubator, fill = incubator)) +
  geom_line(size = 1.5) +
  scale_color_manual (name= "Incubator", values = c ("Fellfield"= "chocolate2", "Snowbed" ="deepskyblue3")) +
  facet_wrap(~ species, scales = "free_x", ncol = 2) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Time (days)", y = "Germination proportion") +
  theme(strip.text = element_text(face = "italic", size = 20), 
        legend.title = element_text(size=18), 
        legend.text = element_text(size=14),
        axis.title.y = element_text (size=16), 
        axis.title.x = element_text (size=16)) + 
  geom_vline(xintercept = 122, linetype = "dashed", size= 1.5) +
  geom_vline(xintercept = 248, linetype = "dashed", size= 1.5, color = "chocolate2") +
  geom_vline(xintercept = 290, linetype = "dashed", size= 1.5, color = "deepskyblue3") 
  
# Save the plots
ggsave(filename = "results/FigAgrostistileni.png", Agrostistileni, path = NULL, 
       scale = 1, width = 180, height = 180, units = "mm", dpi = 600)

#### germination peaks identification #####
read.csv("data/R long data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date)))) %>%
  group_by(mountain, species, incubator, code, petridish) %>%
  filter(date == max(date)) %>%
  select(mountain, species, incubator, code, petridish, viable) %>%
  group_by(mountain, incubator) %>%
  summarise(viable = sum(viable)) -> viables_peak
x11()
read.csv("data/R long data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date)))) %>%
  group_by(mountain, incubator, time) %>%
  summarise(germinated = sum(germinated)) %>%
  #mutate(germinated = cumsum(germinated)) %>%
  merge(viables_peak) %>%
  mutate(germination = germinated/viable) %>%
  ggplot(aes(time, germination, color = incubator, fill = incubator)) +
  geom_line(size = 1.5) +
  scale_color_manual (name= "Incubator", values = c ("Fellfield"= "chocolate2", "Snowbed" ="deepskyblue3")) +
  facet_wrap(~ mountain, scales = "free_x", ncol = 2) +
  labs(x = "Time (days)", y = "Germination proportion") +
  theme(strip.text = element_text(face = "italic", size = 20), 
        legend.title = element_text(size=18), 
        legend.text = element_text(size=14),
        axis.title.y = element_text (size=16), 
        axis.title.x = element_text (size=16)) + 
  geom_vline(xintercept = 122, linetype = "dashed", size= 1.5) +
  geom_vline(xintercept = 248, linetype = "dashed", size= 1.5, color = "chocolate2") +
  geom_vline(xintercept = 290, linetype = "dashed", size= 1.5, color = "deepskyblue3")


#### cumulative germination all species together compare incubators only ####
read.csv("data/R long data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date)))) %>%
  group_by(species, incubator, code, petridish) %>%
  filter(date == max(date)) %>%
  select(species, mountain, incubator, code, petridish, viable) %>%
  group_by(incubator, mountain) %>%
  summarise(viable = sum(viable)) -> viables
# write.csv (viables,"results/viables.csv", row.names = FALSE )
# tidyverse modification to have the accumulated germination along the whole experiment + ggplot
read.csv("data/R long data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date)))) %>%
  group_by(mountain, incubator, time) %>%
  summarise(germinated = sum(germinated)) %>%
  mutate(germinated = cumsum(germinated)) %>%
  merge(viables) %>%
  mutate(germination = germinated/viable) %>%
  ggplot(aes(time, germination, color = incubator, fill = incubator)) +
  geom_line(size = 1.5) +
  scale_color_manual (name= "Incubator", values = c ("Fellfield"= "chocolate2", "Snowbed" ="deepskyblue3")) +
  facet_wrap(~ mountain, scales = "free_x", ncol = 2) +
  labs(title = "", x = "Time (days)", y = "Germination proportion") +
  theme(strip.text = element_text(face = "italic", size = 20), 
             legend.title = element_text(size=18), 
             legend.text = element_text(size=14),
             axis.title.y = element_text (size=16), 
             axis.title.x = element_text (size=16)) + 
  geom_vline(xintercept = 122, linetype = "dashed", size= 1.5) +
  geom_vline(xintercept = 248, linetype = "dashed", size= 1.5, color = "chocolate2") +
  geom_vline(xintercept = 290, linetype = "dashed", size= 1.5, color = "deepskyblue3") 

#### 2nd sow data transformation and graph####
# tidyverse transformation to account for the number of viable seeds per each specie and incubator
# summing up petridishes and accesions/populations of the same species (not taking into account weekly germination)
read.csv("data/Second_sow.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date)))) %>%
  group_by(species, incubator, code, petridish) %>%
  filter(date == max(date)) %>%
  select(species, incubator, code, petridish, viable) %>%
  group_by(species, incubator) %>%
  summarise(viable = sum(viable)) -> secondsowviables
###  germination curves 2nd sow
read.csv("data/Second_sow.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date)))) %>%
  group_by(species, incubator, time) %>%
  summarise(germinated = sum(germinated)) %>%
  mutate(germinated = cumsum(germinated)) %>%
  merge(viables) %>%
  mutate(germination = germinated/viable) %>%
  filter(species == "Thymus praecox") %>%
  ggplot(aes(time, germination, color = incubator, fill = incubator)) +
  geom_line(size = 1) +
  scale_color_manual (name= "incubator", values = c ("Fellfield"= "red", "Snowbed" ="turquoise")) +
  facet_wrap(~ species, scales = "free_x", ncol = 2) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Time (days)", y = "Germination proportion") +
  theme(strip.text = element_text(face = "italic", size = 16), 
        legend.title = element_text(size=16), 
        legend.text = element_text(size=12),
        axis.title.y = element_text (size=14), 
        axis.title.x = element_text (size=14)) + 
  geom_vline(xintercept = 39, linetype = "dashed", size= 1) +
  geom_vline(xintercept = 165, linetype = "dashed", size= 1, color = "red") +
  geom_vline(xintercept = 207, linetype = "dashed", size= 1, color = "turquoise") 


#### intraespecific variation (filter each species with 2 populations) ####
read.csv("data/R long data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date)))) %>%
  group_by(species, incubator, code, petridish) %>%
  filter(date == max(date)) %>%
  select(species, incubator, code, petridish, viable) %>%
  group_by(species, incubator) %>%
  summarise(viable = sum(viable)) -> viables

# write.csv (viables,"results/viables.csv", row.names = FALSE )

# tidyverse modification to have the accumulated germination along the whole experiment + ggplot
read.csv("data/R long data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date)))) %>%
  group_by(species, incubator, code, time) %>%
  summarise(germinated = sum(germinated)) %>%
  mutate(germinated = cumsum(germinated)) %>%
  merge(viables) %>%
  mutate(germination = germinated/viable) %>%
  filter(species == "Phalacrocarpum oppositifolium") %>%
  ggplot(aes(time, germination, color = incubator, fill = incubator)) +
  geom_line(size = 1.5) +
  scale_color_manual (name= "Incubator", values = c ("Fellfield"= "chocolate2", "Snowbed" ="deepskyblue3")) +
  facet_wrap(~ code, scales = "free_x", ncol = 2) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Phalacrocarpum oppositifolium", x = "Time (days)", y = "Germination proportion") +
  theme(title = element_text(face ="italic", size =16),
        legend.title = element_text(size=18), 
        legend.text = element_text(size=14),
        axis.title.y = element_text (size=16), 
        axis.title.x = element_text (size=16)) + 
  geom_vline(xintercept = 122, linetype = "dashed", size= 1.5) +
  geom_vline(xintercept = 248, linetype = "dashed", size= 1.5, color = "chocolate2") +
  geom_vline(xintercept = 290, linetype = "dashed", size= 1.5, color = "deepskyblue3") 

