library(tidyverse);library(ggpubr)

## WEEKLY MEANS (2009-2018) FOR PICOS ####
read.csv("data/temp_picos_raw.csv", sep = ";") %>%
  mutate(Time = strptime(as.character(Time), "%d/%m/%Y %H:%M"))%>% #specify format of TIME variable
  mutate(Time = as.POSIXct(Time, tz = "UTC")) %>%
  mutate(Year = lubridate::year(Time)) %>% 
  mutate(Month = lubridate::month(Time)) %>% 
  mutate(Day = lubridate::day(Time)) %>% 
  mutate(Hour = lubridate::hour(Time)) %>% 
  group_by(Site, Day = lubridate::floor_date(Time, "day")) %>%
  summarise(T = mean(Temperature), X = max(Temperature), N = min(Temperature), n = length(Time)) %>% # Daily mean, max, min
  mutate(Snow = ifelse(X < 0.5 & N > -0.5, 1, 0)) %>% # Day is snow day or not
  mutate(FreezeThaw = ifelse(X > 0.5 & N < -0.5, 1, 0)) %>% # Day with freeze-thaw cycles
  mutate(FDD = ifelse(T < 0, T, 0)) %>% # Freezing degrees per day
  mutate(GDD = ifelse(T >= 1, T, 0)) %>% # Growing degrees day per month https://link.springer.com/article/10.1007/s00035-021-00250-1
  mutate (week = strftime(Day, format = "%V")) %>%
  select(Site, T, X, N, Snow, FDD, GDD, week)%>%
  group_by (Site, week) %>%
  summarise_all(mean, na.rm = TRUE)%>% 
  mutate (Site = as.factor(Site)) %>%
  mutate(week = as.numeric(week)) %>%
  as.data.frame() -> weekly_picos

#write.csv(weekly_picos, "data/weekly_picos.csv") 
#modification in excel to add week order and  extend period to match experiment length !!