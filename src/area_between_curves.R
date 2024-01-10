library (tidyverse); library(scales)

##### Calculate areas between curves (Fellfield-Snowbed) for each species ####
# CALCULATED INDIVIDUALLY PER EACH SPECIES AND ADDED MANUALLY TO all_info.csv
#visualization
x11()
read.csv("data/clean data.csv", sep = ";") %>%
  rbind(data_vis)%>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  spread(date, germinated, fill = 0) %>% # wide format for dates, and fill Na with 0
  gather ("date", "germinated", 8: last_col() )%>% # back in long format frrom 8th colum to the last
  arrange (species, code, incubator, petridish, date)%>% # sort row observations this way
  mutate(date = as.POSIXct(date))%>%
  merge (species) %>%
  dplyr::group_by(community, species, code,  incubator, date) %>%
  dplyr::summarise(germinated = sum(germinated)) %>%
  dplyr::mutate(germinated = cumsum(germinated)) %>%
  merge(viables_pop) %>%
  mutate(germination = germinated/viable)  %>%
  mutate(ID = paste(community, code)) %>%
  filter(species == "Sesleria caerula") %>% ### change species name
  ggplot(aes(date, germination, color = incubator, fill = incubator)) +
  geom_line(size = 2) +
  scale_color_manual (name= "Incubator", values = c ("Fellfield"= "chocolate2", "Snowbed" ="deepskyblue3")) +
  facet_wrap(~ code, scales = "free_x", ncol = 2) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title= "Sesleria caerula", x = "Time ", y = "Germination proportion") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text (size = 30),
        strip.text = element_text (size = 24, face = "italic"),
        axis.title.y = element_text (size=24), 
        axis.title.x = element_text (size=24), 
        axis.text.x= element_text (size=18, angle = 75, vjust = 0.5),
        legend.position = "right",
        plot.margin = margin(t=0.5, l =0.5, b = 0.5, r =0.5, unit ="cm")) #-> C

#example with 1 species (x = time (not dates))
read.csv("data/clean data.csv", sep = ";") %>%
  #rbind(data_vis)%>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  spread(time, germinated, fill = 0) %>% # wide format for dates, and fill Na with 0
  mutate(date = as.POSIXct(date))%>%
  gather ("time", "germinated", 9: last_col()) %>% # back in long format from 9th colum to the last
  mutate (time = as.numeric(time)) %>%
  arrange (species, code, incubator, petridish, time)%>% # sort row observations this way
  merge (species) %>%
  dplyr::group_by(community, species, code,  incubator, time) %>%
  dplyr::summarise(germinated = sum(germinated)) %>%
  dplyr::mutate(germinated = cumsum(germinated)) %>%
  merge(viables_pop) %>%
  mutate(germination = germinated/viable)  %>%
  filter (community == "Temperate") %>%
  filter(species == "Sesleria caerula")%>%
  select( species, code, incubator, community, time, germination)%>%
  #dplyr::group_by(species, incubator, time) %>%   # command line for 2 populations sp
  #dplyr::summarise(germination = mean(germination)) %>%  # command line for 2 populations sp
  spread (incubator, germination) %>%
  data.frame()-> df

#calculate area by trapezoidal rule for 1 pop species
# 1st define each curve coordenates (within each species)
d0 <- df[,c(4,6)] # snowbed data
d1 <- df[,4:5] # snowbed fellfield
f0 <- approxfun(x= d0$time, y= d0$Snowbed)
f1 <- approxfun(x= d1$time, y= d1$Fellfield) 
# define x range of the density overlap
ovrng <- c(min (d0$time), max(d0$time))
# dividing it to sections (for example n=500)
i <- seq(min(ovrng), max(ovrng), length.out=500)
# calculating the distance between the density curves
h1 <- f0(i)-f1(i) # snowbed
h2 <- f1(i)-f0(i) # fellfield
#and using the formula for the area of a trapezoid we add up the areas
area1<-sum( (h1[-1]+h1[-length(h1)]) /2 *diff(i) *(h1[-1]>=0+0))     # for the regions where d1>d0
area2<-sum( (h2[-1]+h2[-length(h2)]) /2 *diff(i) *(h2[-1]>=0+0))     # for the regions where d1<d0
area_total <- area2 - area1 
area_total

#calculate area by trapezoidal rule for 2 pop species 
# https://stackoverflow.com/questions/24742677/how-to-measure-area-between-2-distribution-curves-in-r-ggplot2
# 1st define each curve coordenates (within each species)
d0 <- df[,c(2,4)] # snowbed data
d1 <- df[,2:3] # snowbed fellfield
f0 <- approxfun(x= d0$time, y= d0$Snowbed)
f1 <- approxfun(x= d1$time, y= d1$Fellfield) 
# define x range of the density overlap
ovrng <- c(min (d0$time), max(d0$time))
# dividing it to sections (for example n=500)
i <- seq(min(ovrng), max(ovrng), length.out=500)
# calculating the distance between the density curves
h1 <- f0(i)-f1(i) # snowbed
h2 <- f1(i)-f0(i) # fellfield
#and using the formula for the area of a trapezoid we add up the areas
area1<-sum( (h1[-1]+h1[-length(h1)]) /2 *diff(i) *(h1[-1]>=0+0))     # for the regions where d1>d0
area2<-sum( (h2[-1]+h2[-length(h2)]) /2 *diff(i) *(h2[-1]>=0+0))     # for the regions where d1<d0
area_total <- area2 - area1 
area_total

##### visualization density plots #####
read.csv("data/all_info.csv", sep = ";") %>%
  select(community, species, family, habitat, ABC_clean_data) %>%
  group_by (community, species, family, habitat) %>%
  summarise(area_curves = mean(ABC_clean_data))%>%
  filter (! area_curves==0) %>% # remove species with 0 difference meaning they had 0 germ
  filter (community == "Temperate")%>%
  ggplot(aes(x= area_curves))+  #
  geom_density(alpha =0.2, fill = "azure4", color="azure4" ) + #, position="stack"
  geom_vline(xintercept = 0, linetype = "dashed", size =1) +
  scale_x_continuous( limits = c(-100, 200), breaks = seq (-100, 200, by= 100)) +
  scale_y_continuous(labels = percent,limits = c(0,0.0105)) +
  #facet_wrap(~community)+
  theme_classic(base_size = 16) +
  labs (title= "Germination shift")+
  theme (plot.title = element_text (face = "bold",size = 16), #hjust = 0.5,
         axis.title.y = element_text (size=14),
         axis.text.y = element_text (size = 13),
         axis.title.x = element_blank(), 
         axis.text.x= element_text (size = 13, color = "black"),
         #strip.text = element_text( size = 20, hjust = 0),
         #strip.background = element_blank(), 
         panel.background = element_blank(), #element_rect(color = "black", fill = NULL), 
         #panel.grid = element_blank(),
         legend.title = element_text (size =14),
         legend.text = element_text (size =14),
         legend.position = "bottom", # legend.position = c(0.85, 0.5),
         legend.box.background = element_rect(color = "black", size = 2))


# test if area between curves depends on habitat specialist vs generalist
read.csv("data/all_info.csv", sep =";") %>%
  select(species, family, community, habitat, ABC_clean_data, germ_period_F, germ_period_S) %>%
  group_by (community, habitat) %>%
  count(germ_period_S)%>%
  ggplot( aes(x= germ_period_F))+
  geom_bar(stat= "count", position = "stack") +

#### delay time to reach 50% germ in days between incubators ####
t50model %>%
  ungroup() %>%
  #select(species, code, incubator, t50lm) %>% 
  group_by (species, code, incubator) %>%# only possible if we join data by species and incubator, adding petridish produce an error
  summarise(t50lm = mean(t50lm)) %>%
  spread(incubator, t50lm) %>% 
  mutate(delayS_F = Snowbed - Fellfield) %>% 
  group_by(species) %>%
  summarise(t50_Fellfield = mean(Fellfield),
            t50_Snowbed = mean (Snowbed),
            delayS_F= mean (delayS_F))-> delaytime # NAs appear when species don't reach 50% germination (t50lm_days =Na)
