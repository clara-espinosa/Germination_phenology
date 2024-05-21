library (tidyverse);library(dplyr)
library(plyr)
# based on MCMC_traits script, refer there for detailed explanations
# appendix table with sp data + traits
read.csv ("data/all_info.csv", sep = ";") -> species
# viables + final germ + filtering variables ####
detach(package:plyr)
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  group_by(species, code, incubator, petridish) %>%
  summarise(total_germ = sum(germinated),
            viable= last(viable),
            total = last(total)) %>% 
  merge(species)%>% 
  group_by (community, species, code, incubator, family, recolection) %>%
  summarise (total_germ = sum (total_germ), 
             viable = sum (viable), 
             total = sum(total), 
             abundance = mean (abundance),
             area_curves= mean(ABC_clean_data))%>%
  mutate(germPER = (total_germ/viable) *100, # 
         germPER = round (germPER, digit =2)) %>%
  mutate(viablePER=(viable/total) *100, 
         viablePER = round(viablePER, digit =2)) %>%
  arrange(community, species, code)  -> appendix 
#autumn
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  group_by (species, code, incubator, petridish, time) %>%
  summarise(seeds_germ = sum(germinated)) %>%
  filter (between(time, 1, 105)) %>% ## 105 = 12/11 last check before winter
  group_by (species,code, incubator) %>%
  summarise (seeds_germ = sum(seeds_germ)) %>%
  merge(appendix)%>%
  group_by (community, species,  code, incubator)%>%
  mutate(autumn_germ = (seeds_germ/viable)*100, 
         autumn_germ = round (autumn_germ, digit = 2))%>%
  select (community, species, abundance, family, code,  recolection, incubator, 
           total,viable, viablePER, total_germ,germPER, area_curves, autumn_germ) -> appendix
# winter
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))  %>%
  spread(date, germinated, fill = 0) %>% # wide format for dates, and fill Na with 0
  gather ("date", "germinated", 8: last_col()) %>% # back in long format from 8th colum to the last
  arrange (species, code, incubator, petridish, date)%>% # sort row observations this way
  group_by (species, code, incubator, petridish, date)%>% 
  summarise(seeds_germ = sum(germinated)) %>%
  filter (date %in% c("2022-04-04", # # first check after winter before Tmin >2ºC
                      "2022-05-11","2022-05-12", "2022-05-18", "2022-05-19", 
                      "2022-05-25", "2022-05-26")) %>% #undersnow + first checks Tmin>2ºC
  ungroup () %>%
  group_by (species, code, incubator) %>%
  summarise (seeds_germ = sum(seeds_germ)) %>%
  merge(appendix)%>%
  group_by (community, species, code, incubator)%>%
  mutate(winter_germ = (seeds_germ/viable)*100, 
         winter_germ = round (winter_germ, digit = 2))%>%
  select (community, species, abundance, family, code,  recolection, incubator, 
          total,viable, viablePER, total_germ,germPER,  area_curves, autumn_germ, winter_germ) -> appendix
#spring
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))  %>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  spread(time, germinated, fill = 0) %>% # wide format for dates, and fill Na with 0
  dplyr::select (!date) %>%
  gather ("time", "germinated", 7: last_col()) %>% # back in long format from 8th colum to the last
  group_by (species, code, incubator, petridish, time) %>% 
  filter (incubator == "Fellfield") %>% 
  filter (time %in% 254:331) %>% ## amount of days calculated from the dates 254 = 10/4 and 331 = 25/06
  group_by (species, code, incubator, petridish) %>%
  summarise(seeds_germ = sum(germinated)) %>%
  merge (viables)  -> Springgerm_F

read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))  %>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  spread(time, germinated, fill = 0) %>% # wide format for dates, and fill Na with 0
  dplyr::select (!date) %>%
  gather ("time", "germinated", 7: last_col()) %>% # back in long format from 8th colum to the last
  group_by (species, code, incubator, petridish, time) %>% 
  filter (incubator == "Snowbed") %>% 
  filter (time %in% 301:331) %>% ## amount of days calculated from the dates 301 = 27/5 and 331 = 25/06
  group_by (species, code, incubator, petridish) %>%
  summarise(seeds_germ = sum(germinated)) %>%
  merge (viables)  -> Springgerm_S

rbind(Springgerm_F, Springgerm_S) %>%
  select(species, code, incubator, petridish, seeds_germ, viable) %>%
  group_by (species, code, incubator) %>%
  summarise (seeds_germ = sum(seeds_germ)) %>%
  merge(appendix)%>%
  group_by (community, species,  code,incubator)%>%
  mutate(spring_germ = (seeds_germ/viable)*100, 
         spring_germ = round (spring_germ, digit = 2))%>%
  select (community, species, abundance, family,  code, recolection, incubator, 
          total,viable, viablePER, total_germ,germPER,  area_curves, autumn_germ, winter_germ, spring_germ) -> appendix
#summer
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))  %>%
  spread(date, germinated, fill = 0) %>% # wide format for dates, and fill Na with 0
  gather ("date", "germinated", 8: last_col() )%>% # back in long format frrom 8th colum to the last
  arrange (species, code, incubator, petridish, date)%>% # sort row observations this way
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  group_by (species, code, incubator, petridish, time)%>% 
  summarise(seeds_germ = sum(germinated)) %>%
  mutate(cumulative = cumsum(seeds_germ)) %>%
  filter (between(time, 331, 430)) %>% # amount of days calculated from the dates 331 = 25/06 and 430 = 19/09
  ungroup () %>%
  group_by (species,code, incubator) %>%
  summarise (seeds_germ = sum(seeds_germ)) %>%
  merge(appendix)%>%
  group_by (community, species, code, incubator)%>%
  mutate(summer_germ = (seeds_germ/viable)*100, 
         summer_germ = round (summer_germ, digit = 2))%>%
  select (community, species, abundance, family,  code, recolection, incubator, 
          total,viable, viablePER, total_germ,germPER,  area_curves, autumn_germ, winter_germ, spring_germ,
          summer_germ) -> appendix

# t50
f50 <- function(df0) {
  lm(t50g ~ t50times, data = df0) -> mf1 # Linear model between time before and time after
  as.data.frame((0.5 - as.numeric(mf1$`coefficients`[1])) / as.numeric(mf1$`coefficients`[2])) # Use linear model to interpolate t50
}

read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y")) %>%
  group_by(species, code, incubator, petridish) %>%
  mutate(days = difftime(date, min(date), units = "days")) %>% # calculate time from sowing to x date!
  arrange(species, code,  incubator, petridish, days) %>%
  group_by (species, code, incubator, petridish) %>% 
  mutate(cs = cumsum(germinated), # cumulative sum of germinated seeds 
         fg = max(cs) / viable, # max germinated (final germination) divided by viable seeds (proportion)
         g = cs / viable, # proportion of germination at each date
         Timeframe = ifelse(g >= 0.5, "Upper", "Lower")) %>% # # New column that divides dataset in germination < 0.5 or > 0.5
  group_by(species, code, incubator, petridish, Timeframe) %>%
  mutate(t50times = ifelse(Timeframe == "Lower", max(days), min(days)), # keep only the number of days just before and after reaching t50
         t50g = ifelse(days == t50times, g, NA)) %>% # copy of germination proportion the days just before and after t50
  na.omit %>%
  select(species,  code, incubator, petridish, Timeframe, fg, t50times, t50g) %>%
  unique %>%
  group_by(species,  code, incubator, petridish) %>%
  filter(fg >= .50) %>% # filter to keep only the species that reach more than 0.5 in max final germination
  group_by(species,  code, incubator, petridish) %>%
  do(f50(.)) %>%
  rename(#FG = fg, 
    t50lm = `(0.5 - as.numeric(mf1$coefficients[1]))/as.numeric(mf1$coefficients[2])`) ->t50model

t50model %>%
 group_by(species, code, incubator) %>%
  summarise(t50lm = mean (t50lm))%>%
  right_join(appendix, by =c("species",  "code", "incubator")) %>%
  select(community, species, abundance, family,  code, recolection, incubator, 
        total,viable, viablePER, total_germ, germPER,  area_curves, autumn_germ, winter_germ, spring_germ,
        summer_germ,t50lm) %>% 
  arrange (species, code, incubator)-> appendix

#heat sum
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  group_by (species, code, incubator, petridish) %>%
  summarise (datesow = first (date)) %>% 
  merge (viables)%>% 
  filter (germPER>49) %>%
  merge (t50model) %>%
  mutate (date50 = as.Date(datesow) + t50lm)  -> t50_dates 

heatsum <-function (t50_dates) {
  t50_dates %>%
    pull (incubator)%>%
    unique() -> inc
  t50_dates  %>%
    pull(datesow) -> date1 #convierte columna en vector
  t50_dates %>%
    pull(date50) -> date2
  temp %>% 
    mutate (date = as.Date(date)) %>%
    filter(incubator == inc) %>%
    filter(date >= date1 & date <= date2) %>%
    summarise (HS =sum(Tmean)) 
}

t50_dates%>%
  group_by (species,code, incubator, petridish) %>%
  do (heatsum(.)) -> HS

HS %>%
  group_by (species,code, incubator) %>%
  summarise(HS = mean(HS)) %>%
  right_join(appendix, by =c("species",  "code", "incubator")) %>%
  select(community, species, abundance, family,  code, recolection, incubator, 
         total,viable, viablePER, total_germ,germPER,  area_curves, autumn_germ,winter_germ,  spring_germ,
         summer_germ, t50lm, HS) %>% 
  arrange (species,  code, incubator)-> appendix


# final table, only selected variables
appendix%>%
  ungroup()%>%
  select(community, species, family, abundance, incubator, area_curves, 
         autumn_germ,winter_germ,  spring_germ, summer_germ, t50lm, HS) %>% 
  group_by(community, species, family, incubator)%>% 
  summarise(across(everything(), .f= mean , na.rm = TRUE)) %>%
  write.csv("results/Tables/2. Traits per sp summary.csv")

appendix %>% write.csv("results/Supplementary/Species traits summary.csv")
  str(appendix)
appendix %>%
  select(species, incubator, HS)%>%
  spread(incubator, HS)%>%
  print(n=95)
