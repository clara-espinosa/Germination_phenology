library(tidyverse)

### Annual indices PICOS 10 years of data ####
read.csv("data/temperatures-picos.csv", sep = ";") %>%
  mutate(Time = strptime(as.character(Time), "%d/%m/%Y %H:%M"))%>% #specify format of TIME variable
  mutate(Time = as.POSIXct(Time, tz = "UTC")) %>% 
  mutate(Month = lubridate::month(Time)) %>% 
  mutate(Day = lubridate::day(Time)) %>% 
  mutate(Hour = lubridate::hour(Time)) %>% 
  #filter(Hour %in% c(0, 4, 8, 12, 16, 20)) %>% # Keep same recording hours as ibuttons
  filter(! Month %in% 9) %>% # Remove missing September days
  filter(! (Month %in% 8 & Day > 6)) %>% # Remove missing August days
  filter(! (Month %in% 10 & Day < 3)) %>% # Remove missing October days %>%
  dplyr::select(-c(Month, Day, Hour)) %>%
  group_by(Site, Day = lubridate::floor_date(Time, "day")) %>%
  summarise(T = mean(Temperature), X = max(Temperature), N = min(Temperature), n = length(Time)) %>% # Daily mean, max, min
  mutate(Snow = ifelse(X < 0.5 & N > -0.5, 1, 0)) %>% # Day is snow day or not
  # mutate(Snow = ifelse((X - N) <= 0.5, 1, 0)) %>% # Day is snow day or not
  mutate(FreezeThaw = ifelse(X > 0.5 & N < -0.5, 1, 0)) %>% # Day with freeze-thaw cycles
  mutate(FDD = ifelse(T < 0, T, 0)) %>% # Freezing degrees per day
  mutate(GDD = ifelse(T >= 1, T, 0)) %>% # Growing degrees day per month https://link.springer.com/article/10.1007/s00035-021-00250-1
  group_by(Site, Month = lubridate::floor_date(Day, "month")) %>%
  summarise(T = mean(T), X = mean(X), N = mean(N), # Daily mean, max, min
            Snow = sum(Snow), # Snow days per month
            FreezeThaw = sum(FreezeThaw), # Freeze-thaw days per month
            FDD = sum(FDD), # FDD per month
            GDD = sum(GDD)) %>% # GDD per month
  group_by(Site, Year = lubridate::floor_date(Month, "year")) %>%
  summarise(bio1 = mean(T), # Annual Mean Temperature
            bio2 = mean(X - N), # Mean Diurnal Range (Mean of monthly (max temp - min temp))
            bio7 = max(X) - min(N), # Temperature Annual Range (BIO5-BIO6)
            Snw = sum(Snow),
            FDD = abs(sum(FDD)), # FDD per year
            GDD = sum(GDD)) %>%
  mutate(Year = as.character(lubridate::year(Year))) %>% 
  #rename(Plot = Year) %>%
  #mutate(Survey = "Temporal") 
  select(Site, bio1, bio2, bio7, Snw, FDD, GDD)%>% 
  group_by(Site)%>%
  summarise_all(mean, na.rm = TRUE) -> annual_picos

#write.csv("data/annual-indices.csv")

### Annual indices Omaña 1year of data ####
read.csv("data/temperatures-villabandin.csv", sep = ";") %>% # exclude Cañada site bc data only goes for 2 months
  mutate(Time = strptime(as.character(Time), "%d/%m/%Y %H:%M"))%>% #specify format of TIME variable
  mutate(Time = as.POSIXct(Time, tz = "UTC")) %>% 
  mutate(Month = lubridate::month(Time)) %>% 
  mutate(Day = lubridate::day(Time)) %>% 
  mutate(Hour = lubridate::hour(Time)) %>% 
  #filter(! Month %in% 6) %>% # Remove missing June days (loggers buried later in rabinalto)
  dplyr::select(-c(Month, Day, Hour)) %>%
  group_by(Site, Day = lubridate::floor_date(Time, "day")) %>%
  summarise(T = mean(Temperature), X = max(Temperature), N = min(Temperature), n = length(Time)) %>% # Daily mean, max, min
  mutate(Snow = ifelse(X < 0.5 & N > -0.5, 1, 0)) %>% # Day is snow day or not
  mutate(FreezeThaw = ifelse(X > 0.5 & N < -0.5, 1, 0)) %>% # Day with freeze-thaw cycles
  mutate(FDD = ifelse(T < 0, T, 0)) %>% # Freezing degrees per day
  mutate(GDD = ifelse(T >= 1, T, 0)) %>% # Growing degrees day per month https://link.springer.com/article/10.1007/s00035-021-00250-1
  group_by(Site, Month = lubridate::floor_date(Day, "month")) %>%
  summarise(T = mean(T), X = mean(X), N = mean(N), # Daily mean, max, min
            Snow = sum(Snow), # Snow days per month
            FreezeThaw = sum(FreezeThaw), # Freeze-thaw days per month
            FDD = sum(FDD), # FDD per month
            GDD = sum(GDD)) %>% # GDD per month
  group_by(Site) %>% #, Year = lubridate::floor_date(Month, "year")) %>%
  summarise(bio1 = mean(T), # Annual Mean Temperature
            bio2 = mean(X - N), # Mean Diurnal Range (Mean of monthly (max temp - min temp))
            bio7 = max(X) - min(N), # Temperature Annual Range (BIO5-BIO6)
            Snw = sum(Snow),
            FDD = abs(sum(FDD)), # FDD per year
            GDD = sum(GDD)) -> annual_villa
rbind (annual_picos, annual_villa) %>%
write.csv("data/annual_indices.csv")
##graph to compare with temperature programs ####
x11()
read.csv("data/clima_week_20-21.csv", sep =";") %>%
  #filter (!(Site == "Los Boches")) %>%
  group_by(Macroclimate, week, order) %>%
  summarise(T = mean(T), X = mean(X), N = mean(N)) %>% ### max and min for range, Or mean 
  ggplot(aes ()) +
  geom_line (aes (x=order, y=N, colour = Macroclimate), size =1.25) + #order for sort as experiment
  geom_line (aes (x=order, y=X, colour = Macroclimate), size =1.25) +
  geom_ribbon (aes (x=order, ymin =N, ymax=X, colour = Macroclimate, fill = Macroclimate), alpha =0.3) +
  scale_color_manual (name = "Macroclimate", values = c ("Mediterranean"= "darkgoldenrod1" , "Temperate" = "forestgreen")) +
  scale_fill_manual (name = "Macroclimate", values = c ("Mediterranean"= "darkgoldenrod1" , "Temperate" = "forestgreen")) +
  scale_y_continuous (limits = c(-5,35), breaks = seq (-5, 35, by= 5)) +
  labs (title = "Mean Tmax and Tmin 2021-2022", y= "Temperature ºC", x = "Time (weeks)") + 
  theme_classic(base_size = 16) +
  theme (plot.title = element_text ( size = 30), #hjust = 0.5,
         axis.title.y = element_text (size=18), 
         axis.title.x = element_text (size=18), 
         legend.title = element_text(size = 20),
         legend.text = element_text (size =16)) +
  geom_hline(yintercept=0, linetype ="dashed", size =1, colour = "red")
## BOCHES EXTRA DATA ####
read.csv("data/temp_picos.csv", sep = ";") %>%
  mutate(Time = strptime(as.character(Time), "%d/%m/%Y %H:%M"))%>% #specify format of TIME variable
  mutate(Time = as.POSIXct(Time, tz = "UTC")) %>%
  mutate(Year = lubridate::year(Time)) %>% 
  mutate(Month = lubridate::month(Time)) %>% 
  mutate(Day = lubridate::day(Time)) %>% 
  mutate(Hour = lubridate::hour(Time)) %>% 
  filter(Site == "Los Boches") %>%
  #filter (! Year < 2021) %>%
  #filter(! (Year %in% 2021 & Month < 7 )) %>% # Remove all data from June 2021 
  #filter(! (Year %in% 2021 & Month %in% 7 & Day < 6)) %>% # remove data in july 2021 until day 6
  #filter(! (Year %in% 2022 & Month %in% 7 & Day > 5)) %>% # remove data in july 2022 from day 5
  #dplyr::select(-c(Month, Day, Hour)) %>%
  group_by(Site, Day = lubridate::floor_date(Time, "day")) %>%
  summarise(T = mean(Temperature), X = max(Temperature), N = min(Temperature), n = length(Time)) %>% # Daily mean, max, min
  mutate(Snow = ifelse(X < 0.5 & N > -0.5, 1, 0)) %>% # Day is snow day or not
  # mutate(Snow = ifelse((X - N) <= 0.5, 1, 0)) %>% # Day is snow day or not
  mutate(FreezeThaw = ifelse(X > 0.5 & N < -0.5, 1, 0)) %>% # Day with freeze-thaw cycles
  mutate(FDD = ifelse(T < 0, T, 0)) %>% # Freezing degrees per day
  mutate(GDD = ifelse(T >= 1, T, 0)) %>% # Growing degrees day per month https://link.springer.com/article/10.1007/s00035-021-00250-1
  mutate (week = strftime(Day, format = "%V")) %>%
  select(Site, T, X, N, Snow, FDD, GDD, week)%>%
  group_by (Site, week) %>%
  summarise_all(mean, na.rm = TRUE)%>% 
  mutate (Site = as.factor(Site)) %>%
  mutate(week = as.numeric(week)) %>%
  as.data.frame() -> weekly_boches
### Omaña weekly values 1 year of data ####
read.csv("data/temp_villa.csv", sep = ";") %>% # exclude Cañada site bc data only goes for 2 months
  mutate(Time = strptime(as.character(Time), "%d/%m/%Y %H:%M"))%>% #specify format of TIME variable
  mutate(Time = as.POSIXct(Time, tz = "UTC")) %>% 
  mutate(Year = lubridate::year(Time)) %>% 
  mutate(Month = lubridate::month(Time)) %>% 
  mutate(Day = lubridate::day(Time)) %>% 
  mutate(Hour = lubridate::hour(Time)) %>% 
  filter(Site == "Penauta" | Site =="Rabinalto" | Site == "Solana") %>% # | = "or" , avoid PenautaS/PenautaN
  filter(! (Year %in% 2021 & Month %in% 6)) %>% # Remove all data from June 2021 
  filter(! (Year %in% 2021 & Month %in% 7 & Day < 6)) %>% # remove data in july 2021 until day 6
  filter(! (Year %in% 2022 & Month %in% 7 & Day > 5)) %>% # remove data in july 2022 from day 5
  dplyr::select(-c(Year, Month, Day, Hour)) %>%
  group_by(Site, Day = lubridate::floor_date(Time, "day")) %>%
  summarise(T = mean(Temperature), X = max(Temperature), N = min(Temperature), n = length(Time)) %>% # Daily mean, max, min
  mutate(Snow = ifelse(X < 0.5 & N > -0.5, 1, 0)) %>% # Day is snow day or not
  mutate(FreezeThaw = ifelse(X > 0.5 & N < -0.5, 1, 0)) %>% # Day with freeze-thaw cycles
  mutate(FDD = ifelse(T < 0, T, 0)) %>% # Freezing degrees per day
  mutate(GDD = ifelse(T >= 1, T, 0)) %>% # Growing degrees day per month https://link.springer.com/article/10.1007/s00035-021-00250-1
  mutate (week = strftime(Day, format = "%V")) %>%
  #mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  select(Site, week, T, X, N, Snow, FDD, GDD)%>%
  group_by (Site, week) %>%
  summarise_all(mean, na.rm = TRUE) %>%
  mutate (Site = as.factor(Site)) %>%
  mutate(week = as.numeric(week)) %>%
  as.data.frame()-> weekly_villa
 
#write.csv(weekly_villa, "data/weekly-villa.csv")

rbind (picos_vis, villa_vis) -> temp_clima
# write.csv(temp_clima, "data/clima_week.csv")
# visualization 1
# mulitple lines with shaded areas in between
# data reordered in excel so it matches tem_regimes order
read.csv("data/clima_week.csv", sep = ";" ) -> temp_clima
temp_clima %>%
  filter(Site == "Hoyo Sin Tierra" | Site == "Los Boches") ->temp_clima 
  #filter (Site == "Villabandín") -> temp_clima
x11()
ggplot(weekly_villa, aes ()) +
  geom_line (aes (x=week, y=N, colour = Site), size =1.25) + #order for sort as experiment
  geom_line (aes (x=week, y=X, colour = Site), size =1.25) +
  geom_ribbon (data = weekly_villa, aes (x=week, ymin =N, ymax=X, colour = Site, fill = Site), alpha =0.3) +
  #scale_fill_manual (values =c("chocolate2", "deepskyblue3")) + #, "red"
  #scale_color_manual (values =c("chocolate2", "deepskyblue3")) + #, "red"
  scale_y_continuous (limits = c(-10,40), breaks = seq (-10, 40, by= 10)) +
  labs (title = "Mediterranean 2021-2022", y= "Temperature ºC", x = "Time (weeks)") + 
  theme_classic(base_size = 16) +
  theme (plot.title = element_text ( size = 30), #hjust = 0.5,
         axis.title.y = element_text (size=18), 
         axis.title.x = element_text (size=18), 
         legend.title = element_text(size = 20),
         legend.text = element_text (size =16)) +
  #geom_vline(xintercept = 18, linetype = "dashed", size =1.25) + # 16
  #annotate (geom ="text", x= 20, y = 22, label ="Winter begins", colour = "black", size = 4,  fontface ="bold") +
  #geom_vline(xintercept = 35, linetype = "dashed", size =1.25, color = "chocolate2") + # 35
  #annotate (geom ="text", x= 255, y = 22, label ="Fellfield\nwinter ends", colour = "chocolate2", size = 4,  fontface ="bold") +
  #geom_vline(xintercept = 42, linetype = "dashed", size =1.25,color = "deepskyblue3") + # 42
  #annotate (geom ="text", x= 307, y = 22, label ="Snowbed\nwinter ends", colour = "deepskyblue3", size = 4,  fontface ="bold") +
  geom_hline(yintercept=0, linetype ="dashed", size =1, colour = "red")
  



## PICOS WEEKLY values 1 year of data####
read.csv("data/temp_picos.csv", sep = ";") %>%
  mutate(Time = strptime(as.character(Time), "%d/%m/%Y %H:%M"))%>% #specify format of TIME variable
  mutate(Time = as.POSIXct(Time, tz = "UTC")) %>%
  mutate(Year = lubridate::year(Time)) %>% 
  mutate(Month = lubridate::month(Time)) %>% 
  mutate(Day = lubridate::day(Time)) %>% 
  mutate(Hour = lubridate::hour(Time)) %>% 
  filter(Site == "Los Cazadores" | Site =="Hoyo Sin Tierra" | Site == "Hou Sin Tierri") %>%
  filter (! Year < 2021) %>%
  filter(! (Year %in% 2021 & Month < 7 )) %>% # Remove all data from June 2021 
  filter(! (Year %in% 2021 & Month %in% 7 & Day < 6)) %>% # remove data in july 2021 until day 6
  filter(! (Year %in% 2022 & Month %in% 7 & Day > 5)) %>% # remove data in july 2022 from day 5
  dplyr::select(-c(Month, Day, Hour)) %>%
  group_by(Site, Day = lubridate::floor_date(Time, "day")) %>%
  summarise(T = mean(Temperature), X = max(Temperature), N = min(Temperature), n = length(Time)) %>% # Daily mean, max, min
  mutate(Snow = ifelse(X < 0.5 & N > -0.5, 1, 0)) %>% # Day is snow day or not
  # mutate(Snow = ifelse((X - N) <= 0.5, 1, 0)) %>% # Day is snow day or not
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
str(weekly_picos)
summary(weekly_picos)
#write.csv(weekly_picos, "data/weekly-picos.csv")

rbind (weekly_picos, weekly_boches)-> weekly_picosb

ggplot(weekly_picosb, aes ()) +
  geom_line (aes (x=week, y=N, colour = Site), size =1.25) + #order for sort as experiment
  geom_line (aes (x=week, y=X, colour = Site), size =1.25) +
  geom_ribbon (data = weekly_picosb, aes (x=week, ymin =N, ymax=X, colour = Site, fill = Site), alpha =0.3) +
  scale_y_continuous (limits = c(-10,40), breaks = seq (0, 40, by= 10)) +
  labs (title = "Temperate 2021-2022", y= "Temperature ºC", x = "Time (weeks)") + 
  theme_classic(base_size = 16) +
  theme (plot.title = element_text ( size = 30), #hjust = 0.5,
         axis.title.y = element_text (size=18), 
         axis.title.x = element_text (size=18), 
         legend.title = element_text(size = 20),
         legend.text = element_text (size =16)) +
  geom_hline(yintercept=0, linetype ="dashed", size =1, colour = "red")

rbind(weekly_picosb, weekly_villa) %>%
  write.csv( "data/clima_week_20-21.csv")
  
weekly_picosb %>%
  #filter(!(Site =="Los Boches")) %>%
  group_by(week) %>%
  summarise(T = mean(T), X = max(X), N = min(N)) %>%
  mutate (Macroclimate = "Temperate") -> temperate_graph

weekly_villa %>%
  group_by(week) %>%
  summarise(T = mean(T), X = max(X), N = min(N)) %>%
  mutate (Macroclimate = "Mediterranean") -> mediterranean_graph

rbind (temperate_graph, mediterranean_graph) -> clima_graph
x11()
ggplot(clima_graph, aes ()) +
  geom_line (aes (x=week, y=N, colour = Macroclimate), size =1.25) + #order for sort as experiment
  geom_line (aes (x=week, y=X, colour = Macroclimate), size =1.25) +
  geom_ribbon (data = clima_graph, aes (x=week, ymin =N, ymax=X, colour = Macroclimate, fill = Macroclimate), alpha =0.3) +
  scale_color_manual (name = "Macroclimate", values = c ("Mediterranean"= "darkgoldenrod1" , "Temperate" = "forestgreen")) +
  scale_y_continuous (limits = c(-10,40), breaks = seq (0, 40, by= 10)) +
  labs (title = "Temperatures range 2021-2022", y= "Temperature ºC", x = "Time (weeks)") + 
  theme_classic(base_size = 16) +
  theme (plot.title = element_text ( size = 30), #hjust = 0.5,
         axis.title.y = element_text (size=18), 
         axis.title.x = element_text (size=18), 
         legend.title = element_text(size = 20),
         legend.text = element_text (size =16)) +
  geom_hline(yintercept=0, linetype ="dashed", size =1, colour = "red")
