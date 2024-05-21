library(tidyverse);library(ggpubr);library (lubridate)

## WEEKLY MEANS (2009-2018) FOR PICOS ####
read.csv("data/temp_picos_raw.csv", sep = ";") %>% # raw field temperature data from 2008 to 2019 in four sampling site from temperate system 
  mutate(Time = strptime(as.character(Time), "%d/%m/%Y %H:%M"))%>% #specify format of TIME variable
  mutate(Time = as.POSIXct(Time, tz = "UTC")) %>%
  mutate(Month = lubridate::month(Time)) %>% 
  mutate(Day = lubridate::day(Time)) %>% 
  group_by(Site, Day = lubridate::floor_date(Time, "day")) %>% # group temperature data for each day
  summarise(T = mean(Temperature), X = max(Temperature), N = min(Temperature)) %>% # Daily mean, max, min
  mutate(Snow = ifelse(X < 0.5 & N > -0.5, 1, 0)) %>% # Day with or without snow
  mutate (week = strftime(Day, format = "%V")) %>% # create new variable to mean data by weeks to compare to experimental temperature regimes programmed 
  group_by (Site, week) %>%
  summarise_all(mean, na.rm = TRUE)%>% 
  mutate (Site = as.factor(Site)) %>%
  mutate(week = as.numeric(week)) %>%
  filter(!Site== "Hou Sin Tierri")%>%
  filter(!Site== "Los Cazadores")%>%
  as.data.frame() -> weekly_picos

write.csv(weekly_picos, "data/weekly_picos.csv") 
#modification in excel to add week order to match experimental sequence, referred to as weekly_picos_graph for fig 1