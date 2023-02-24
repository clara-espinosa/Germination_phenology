library(FD); library(vegan); library(FactoMineR); library(emmeans);
library(tidyverse); library(ggrepel); library(cowplot); library (ggplot2);
library (lubridate); library(binom); library (GerminaR)
theme_set(theme_cowplot(font_size = 10)) 

####dataframe with temperature programs x incubator####
read.csv("data/date_temp.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y")) %>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  mutate(incubator = as.factor(incubator)) %>% 
  as.data.frame()-> temp
str(temp)
# dataframe with species data ####
read.csv ("data/species.csv", sep =";")-> sp_data
############################# GENERAL GERMINATION TRAITS ##########################
# viable seeds calculation ####
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  group_by(species, code, incubator, petridish) %>%
  filter(date == max(date)) %>%
  mutate (viablePER = (viable/total)*100,
          viablePER = round (viablePER, digit = 2)) -> viables

read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  merge(accessions, by =c("code"))%>%
  group_by(region, species, code, incubator, petridish) %>%
  filter(date == max(date)) %>%
  select(region, species, code, incubator, petridish, viable, total) %>%
  group_by (region, incubator) %>%
  summarise(viable = sum(viable),
            total = sum(total)) -> viables_region

# final germination percentage calculation ####
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  group_by(species, code, incubator, petridish) %>%
  summarise(seeds_germ = sum(germinated)) %>% # 
  merge(viables) %>%
  mutate(germPER = (seeds_germ/viable) *100, # 
         germPER = round (germPER, digit =2)) %>%
  select (species, code, incubator, petridish, seeds_germ, germPER)-> finalgerm 
# seed germinated in spring + summer (needed for early sesason germ calculation) ####
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>% 
  group_by (species, code, incubator, petridish, time) %>%
  filter (incubator == "Fellfield") %>%
  filter (between(time, 257, 430)) %>%
  group_by(species, code, incubator, petridish) %>%
  summarise(germ_AW = sum(germinated)) -> germ_AF_F

read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>% 
  group_by (species, code, incubator, petridish, time) %>%
  filter (incubator == "Snowbed") %>%
  filter (between(time, 301, 430)) %>%
  group_by(species, code, incubator, petridish) %>%
  summarise(germ_AW = sum(germinated)) -> germ_AF_S


rbind(germ_AF_F, germ_AF_S)%>%
  arrange(species) -> germ_AF

# t50 according to germ checks (NOT USED)  ####
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  group_by (species, code, incubator, petridish) %>%
  mutate (datesow= first (date)) %>% # a different column for sowing date
  mutate(days = difftime(date, min(date), units = "days")) %>% # time count from sowing of each petridish
  arrange(species, code, incubator, petridish, days) %>% 
  group_by (species, code, incubator, petridish, days) %>% 
  summarise(seeds_germ = sum(germinated)) %>%
  mutate(cum_seeds_germ = cumsum(seeds_germ)) %>%
  merge (viables) %>%
  mutate(germPER = (cum_seeds_germ/viable) *100, # 
         germPER = round (germPER, digit =2)) %>%
  mutate (t50 = ifelse(germPER >= 50, days, 0)) %>% # CAMBIAR 0 POR double TIME MAX???
  group_by (species, code, incubator, petridish) %>%
  mutate (t50check = first(t50[t50>0])) %>%
  summarise(t50check = min (t50check)) %>%
  ungroup()-> time50 # values discrete bc count t50 according to germination check days

### modeling t50 from raw data ####

f50 <- function(df0) {
  lm(t50g ~ t50times, data = df0) -> mf1 # Linear model between time before and time after
  as.data.frame((0.5 - as.numeric(mf1$`coefficients`[1])) / as.numeric(mf1$`coefficients`[2])) # Use linear model to interpolate t50
}

read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y")) %>%
  group_by(species, code, incubator, petridish) %>%
  mutate(days = difftime(date, min(date), units = "days")) %>% # calculate time from sowing to x date!
  arrange(species, code, incubator, petridish, days) %>%
  group_by (species, code, incubator, petridish) %>% 
  mutate(cs = cumsum(germinated), # cumulative sum of germinated seeds 
         fg = max(cs) / viable, # max germinated (final germination) divided by viable seeds (proportion)
         g = cs / viable, # proportion of germination at each date
         Timeframe = ifelse(g >= 0.5, "Upper", "Lower")) %>% # # New column that divides dataset in germination < 0.5 or > 0.5
  group_by(species, code, incubator, petridish, Timeframe) %>%
  mutate(t50times = ifelse(Timeframe == "Lower", max(days), min(days)), # keep only the number of days just before and after reaching t50
         t50g = ifelse(days == t50times, g, NA)) %>% # copy of germination proportion the days just before and after t50
  na.omit %>%
  select(species, code, incubator, petridish, Timeframe, fg, t50times, t50g) %>%
  unique %>%
  group_by(species, code, incubator, petridish) %>%
  filter(fg >= .50) %>% # filter to keep only the species that reach more than 0.5 in max final germination
  group_by(species, code, incubator, petridish) %>%
  do(f50(.)) %>%
  rename(#FG = fg, 
         t50lm = `(0.5 - as.numeric(mf1$coefficients[1]))/as.numeric(mf1$coefficients[2])`) ->t50model

# COMPARISON T50 FROM CHECK DATE AND T50 FROM LIENAR MODEL , NOT NECESARY, WE WILL USE T50 FROM MODEL
t50_trait <- full_join(time50, t50model, by = c("species", "code", "incubator", "petridish"))

# sowing date + t50 date or final date check, NECESARY FOR ENV HEAT FUNCTION ####
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  group_by (species, code, incubator, petridish) %>%
  mutate (datesow= first (date)) %>% 
  mutate (datelast= last (date)) %>% 
  select(species, code, incubator, petridish, viable, germinated, datesow, datelast) %>%
  merge (t50model) %>%
  group_by (species, code, incubator, petridish, datesow, datelast, t50lm)%>% 
  summarise(t50lm = min (t50lm)) %>% 
  mutate (date50 = as.Date(datesow) + t50lm) %>%  
  mutate (date50 = replace_na (date50, as.character(datelast))) -> t50_dates 
# we dont leave NAs in order to apply heat function properly
### Heat_sum function, CAREFUL FOR <50% GERM HAS BEEN CALCULATED FROM LAST EXPERIMENT DATE####
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
  group_by (species, code, incubator, petridish) %>%
  do (heatsum(.))-> heat_sum # punto significa que el objeto al que aplicar la funcion heat sum es el de la linea de arriba
  
  
### modeling t10 from raw data instead of germination onset######
f10 <- function(t10) {
  lm(t10g ~ t10times, data = t10) -> mf2 # Linear model between time before and time after
  as.data.frame((0.1 - as.numeric(mf2$`coefficients`[1])) / as.numeric(mf2$`coefficients`[2])) # Use linear model to interpolate t50
}

read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y")) %>%
  group_by(species, code, incubator, petridish) %>%
  mutate(days = difftime(date, min(date), units = "days")) %>% # calculate time from sowing to x date!
  arrange(species, code, incubator, petridish, days) %>%
  group_by (species, code, incubator, petridish) %>% 
  mutate(cs = cumsum(germinated), # cumulative sum of germinated seeds 
         fg = max(cs) / viable, # max germinated (final germination) divided by viable seeds (proportion)
         g = cs / viable, # proportion of germination at each date
         Timeframe = ifelse(g >= 0.1, "Upper", "Lower")) %>% # # New column that divides dataset in germination < 0.5 or > 0.5
  group_by(species, code, incubator, petridish, Timeframe) %>%
  mutate(t10times = ifelse(Timeframe == "Lower", max(days), min(days)), # keep only the number of days just before and after reaching t50
         t10g = ifelse(days == t10times, g, NA)) %>% # copy of germination proportion the days just before and after t50
  na.omit %>%
  select(species, code, incubator, petridish, petricode, Timeframe, fg, t10times, t10g) %>%
  unique %>%
  group_by(species, code, incubator, petridish, petricode) %>%
  filter(fg >= .10) %>% # filter to keep only the species that reach more than 0.5 in max final germination
  filter(length(petricode) > 1) %>% # not sure why this comand about lenght of petricode variable, to fix some strange error?
  group_by(species, code, incubator, petridish) %>%
  do(f10(.)) %>%
  rename(#FG = fg, 
    t10lm = `(0.1 - as.numeric(mf2$coefficients[1]))/as.numeric(mf2$coefficients[2])`) ->t10model

model_dates <-full_join(t50_dates, t10model, by= c ("species", "code", "incubator", "petridish"))

# t10 dates from sowing date + t10 model ####
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  group_by (species, code, incubator, petridish) %>%
  mutate (datesow= first (date)) %>% 
  select(species, code, incubator, petridish, viable, germinated, datesow) %>%
  merge (t10model) %>%
  group_by (species, code, incubator, petridish, datesow, t10lm)%>% 
  summarise(t10lm = min (t10lm)) %>% 
  mutate (date10 = as.Date(datesow) + t10lm) %>% 
  ungroup () %>% 
  select (species,code, incubator, petridish, t10lm, date10)-> t10_dates


## germination onset, first check date with seeds germinated (NOT USED)####
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  group_by (species, code, incubator, petridish) %>%
  mutate (datesow= first (date)) %>% # a different column for sowing date
  mutate(days = difftime(date, min(date), units = "days")) %>% # time count from sowing of each petridish
  arrange(species, code, incubator, petridish, days) %>% 
  group_by (species, code, incubator, petridish, days) %>% 
  summarise(seeds_germ = sum(germinated)) %>%
  mutate(cum_seeds_germ = cumsum(seeds_germ)) %>%
  merge (viables) %>%
  mutate(germPER = (cum_seeds_germ/viable) *100, # 
         germPER = round (germPER, digit =2)) %>%
  mutate (firstgerm_days = ifelse(cum_seeds_germ > 0, days, 0)) %>% 
  group_by (species, code, incubator, petridish) %>%
  mutate (firstgerm_days= first(firstgerm_days[firstgerm_days>0])) %>%
  summarise(firstgerm_days = min (firstgerm_days)) %>%
  mutate (firstgerm_days = coalesce(firstgerm_days, 430))-> days_first_germ
  
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  group_by (species, code, incubator, petridish) %>%
  mutate (datesow= first (date)) %>% 
  select(species, code, incubator, petridish, datesow) %>%
  merge (days_first_germ) %>%
  group_by (species, code, incubator, petridish, datesow, firstgerm_days)%>% 
  summarise( firstgerm_days = min ( firstgerm_days)) %>% 
  mutate (onset_date = as.Date(datesow) + firstgerm_days) %>% 
  ungroup() %>% 
  select (species, code, incubator, petridish, firstgerm_days, onset_date) -> germ_onset

############################# GERMINATION TIMING TRAITS ########################### 
## MID-NOVEMBER GERM (proxy of Autumn germination = germination before winter) last check 11/11/2021 ####
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  group_by (species, code, incubator, petridish, time)%>%
  filter (between(time, 1, 105)) %>%
  group_by (species, code, incubator, petridish) %>%
  summarise(seeds_germ = sum(germinated))%>%
  merge(viables) %>%
  mutate(Mid_nov= (seeds_germ/viable) *100,  
         Mid_nov= round (Mid_nov, digit =2)) %>%
  select(species, code, incubator, petridish, Mid_nov) -> Mid_november 


## Winter germination  considering until Tmin<3ºC#####
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))  %>%
  group_by (species, code, incubator, petridish, date)%>% 
  filter (date %in% c("2022-04-04","2022-04-11", "2022-04-12", # # first check after winter before Tmin<3ºC
                      "2022-05-11","2022-05-12", "2022-05-18", "2022-05-19", 
                      "2022-05-25", "2022-05-26", "2022-06-01", "2022-06-02" )) %>% #undersnow + first checks Tmin>2ºC
  group_by (species, code, incubator, petridish) %>% 
  summarise(seeds_germ = sum(germinated))%>%
  merge (viables) %>%
  select(species, code, incubator, petridish, seeds_germ, viable) %>%
  mutate (Winter_germ = (seeds_germ/viable) *100,
          Winter_germ = round (Winter_germ, digit=2)) %>%
  select (species, code, incubator, petridish, Winter_germ) -> Wintergerm

# SNOWBED:  undersnow germination (nº of seeds) +  
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))  %>%
  filter (incubator=="Snowbed") %>%# only able to determine in the snowbed incubator (due to env conditions)
  group_by (species, code, incubator, petridish, date)%>% 
  filter (date %in% c("2022-05-11","2022-05-12")) %>% # date of germination check right before end of winter conditions (dark + 0ÂºC)
  mutate (snowgerm = germinated) %>%
  merge (viables) %>%
  select(species, code, incubator, petridish, viable, snowgerm) %>%
  mutate (snowPER = (snowgerm/viable)*100, 
          snowPER= round (snowPER, digit =2)) %>% 
  replace (is.na(.), 0) %>% 
  select (species, code, incubator, petridish, snowPER) -> snow_trait  # percentage of germinated under snow 
  

## MID-JUNE germination #####
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  filter(species == "Gypsophila repens") %>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  group_by (species, code, incubator, petridish, time)%>% 
  filter (between(time, 106, 322)) %>% # amount of days calculated from the dates 106 = 13/11 and 322 = 16/6
  group_by (species, code, incubator, petridish) %>%
  summarise(seeds_germ = sum(germinated)) %>%
  merge (viables) %>%
  select(species, code, incubator, petridish, seeds_germ, viable) %>%
  mutate (WSgerm = (seeds_germ/viable) *100,
          WSgerm = round (WSgerm, digit=2)) %>%
  select (species, code, incubator, petridish, WSgerm) %>%
  merge (Wintergerm)%>%
  mutate (Mid_june = WSgerm-Wintergerm) %>%
  select (species, code, incubator, petridish, Mid_june)-> Mid_june
## SPRING GERMINATION (when Tmean >10 for 1 week) ####
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))  %>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  group_by (species, code, incubator, petridish, time) %>% 
  filter (incubator == "Fellfield") %>% 
  filter (between(time, 262, 305)) %>% ## amount of days calculated from the dates 262 = 18/5 and 305 = 31/05
  group_by (species, code, incubator, petridish) %>%
  summarise(seeds_germ = sum(germinated)) %>%
  merge (viables) %>%
  select(species, code, incubator, petridish, seeds_germ, viable) %>%
  mutate (spring_germ = (seeds_germ/viable) *100,
          spring_germ = round (spring_germ, digit=2)) %>%
  select (species, code, incubator, petridish, spring_germ) -> Springgerm_F

read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))  %>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  group_by (species, code, incubator, petridish, time) %>% 
  filter (incubator == "Snowbed") %>% 
  filter (between(time, 301, 335)) %>% ## amount of days calculated from the dates 301 = 27/5 and 355 = 30/06
  group_by (species, code, incubator, petridish) %>%
  summarise(seeds_germ = sum(germinated)) %>%
  merge (viables) %>%
  select(species, code, incubator, petridish, seeds_germ, viable) %>%
  mutate (spring_germ = (seeds_germ/viable) *100,
          spring_germ = round (spring_germ, digit=2)) %>%
  select (species, code, incubator, petridish, spring_germ) -> Springgerm_S

Springgerm <- full_join(Springgerm_F, Springgerm_S, by = c("species", "code", "incubator", "petridish", "spring_germ"))
## EARLY SEASON GERMINATION (when Tmax < 10 ) ####
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>% 
  group_by (species, code, incubator, petridish, time) %>% 
  filter (incubator == "Fellfield") %>% 
  filter (between(time, 257, 285)) %>% ## amount of days calculated from the dates 257 = 13/4 and 285 = 11/05
  group_by (species, code, incubator, petridish) %>%
  summarise(seeds_germ = sum(germinated)) %>%
  merge (germ_AF_F) %>%
  mutate (spring_germ = (seeds_germ/germ_AW) *100,
          spring_germ = round (spring_germ, digit=2)) %>%
  mutate (spring_germ= replace_na(spring_germ, 0))%>% # not sure about this line... what to do with species with no germ
  select (species, code, incubator, petridish, spring_germ) -> Earlyseason_F

read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))  %>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  group_by (species, code, incubator, petridish, time) %>% 
  filter (incubator == "Snowbed") %>% 
  filter (between(time, 301, 322)) %>% ## amount of days calculated from the dates 301 = 27/5 and 322 = 17/06
  group_by (species, code, incubator, petridish) %>%
  summarise(seeds_germ = sum(germinated)) %>%
  merge (germ_AF_S) %>%
  mutate (spring_germ = (seeds_germ/germ_AW) *100,
          spring_germ = round (spring_germ, digit=2)) %>%
  mutate (spring_germ= replace_na(spring_germ, 0))%>% # not sure about this line... what to do with species with no germ
  select (species, code, incubator, petridish, spring_germ)  -> Earlyseason_S

rbind(Earlyseason_F, Earlyseason_S) %>%
  arrange(species) -> Earlyseasongerm 

## SUMMER GERMINATION (check after 1 week with max T)####
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))  %>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  group_by (species, code, incubator, petridish, time) %>% 
  filter (incubator == "Fellfield") %>% 
  filter (between(time, 305, 395)) %>% ## amount of days calculated from the dates 305 = 31/05 and 395 = 29/08
  group_by (species, code, incubator, petridish) %>%
  summarise(seeds_germ = sum(germinated)) %>%
  merge (viables) %>%
  select(species, code, incubator, petridish, seeds_germ, viable) %>%
  mutate (summer_germ = (seeds_germ/viable) *100,
          summer_germ = round (summer_germ, digit=2)) %>%
  select (species, code, incubator, petridish, summer_germ) -> summergerm_F

read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))  %>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  group_by (species, code, incubator, petridish, time) %>% 
  filter (incubator == "Snowbed") %>% 
  filter (between(time, 335, 395)) %>% ## amount of days calculated from the dates 355 = 30/06 and 395 = 29/08
  group_by (species, code, incubator, petridish) %>%
  summarise(seeds_germ = sum(germinated)) %>%
  merge (viables) %>%
  select(species, code, incubator, petridish, seeds_germ, viable) %>%
  mutate (summer_germ = (seeds_germ/viable) *100,
          summer_germ = round (summer_germ, digit=2)) %>%
  select (species, code, incubator, petridish, summer_germ) -> summergerm_S

Summergerm <- full_join(summergerm_F, summergerm_S, by = c("species", "code", "incubator", "petridish", "summer_germ"))

## MID-SEPTEMBER germination (end of experiment) ####
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))  %>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  group_by (species, code, incubator, petridish, time)%>% 
  filter (between(time, 322, 430)) %>% # amount of days calculated from the dates 322 = 16/6 and 430 = 19/09
  group_by (species, code, incubator, petridish) %>%
  summarise(seeds_germ = sum(germinated)) %>%
  merge (viables) %>%
  select(species, code, incubator, petridish, seeds_germ, viable) %>%
  mutate (Mid_sep = (seeds_germ/viable) *100,
          Mid_sep = round (Mid_sep, digit=2)) %>%
  select (species, code, incubator, petridish, Mid_sep) -> Mid_september

#### seed germinated in spring + summer (needed for cold/hot germ calculation)####
read.csv("data/all_data.csv", sep = ";") %>%
  spread(date, germinated, fill = 0) %>% # wide format for dates, and fill Na with 0
  gather ("date", "germinated", 8: last_col()) %>% # back in long format frrom 8th colum to the last
  arrange (species, accession, code, incubator, petridish, date)%>% # sort row observations this way
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>% 
  group_by (species, code, incubator, petridish, time) %>%
  filter (incubator == "Fellfield") %>%
  filter (between(time, 249, 430)) %>% # dates 257 = 5/4; 430 = 19/09 end of experiment
  group_by(species, code, incubator, petridish) %>%
  summarise(germ_AW = sum(germinated)) -> germ_AW_F

read.csv("data/all_data.csv", sep = ";") %>%
  spread(date, germinated, fill = 0) %>% # wide format for dates, and fill Na with 0
  gather ("date", "germinated", 8: last_col()) %>% # back in long format frrom 8th colum to the last
  arrange (species, accession, code, incubator, petridish, date)%>% # sort row observations this way
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>% 
  group_by (species, code, incubator, petridish, time) %>%
  filter (incubator == "Snowbed") %>%
  filter (between(time, 301, 430)) %>% # dates 301 = 27/06; 430 = 19/09 end of experiment 
  group_by(species, code, incubator, petridish) %>%
  summarise(germ_AW = sum(germinated)) -> germ_AW_S


rbind(germ_AW_F, germ_AW_S)%>%
  merge(viables_germ) %>%
  arrange(species) -> germ_AW

#### COLD GERMINATION (BETWEEN  Tmin>2ºc AND Tmean<10 )####
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  #filter (species == "Carex sempervirens") %>%
  spread(date, germinated, fill = 0) %>% # wide format for dates, and fill Na with 0
  gather ("date", "germinated", 8: last_col()) %>% # back in long format frrom 8th colum to the last
  arrange (species, accession, code, incubator, petridish, date)%>% # sort row observations this way
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>% 
  group_by (species, code, incubator, petridish, time) %>% 
  summarise(seeds_germ = sum(germinated)) %>%
  mutate(cumulative = cumsum(seeds_germ)) %>%
  filter (incubator == "Fellfield") %>% 
  filter (between(time, 249, 299)) %>% ## amount of days calculated from the dates 249 = 5/4 and 299 = 25/05
  group_by (species, code, incubator, petridish) %>%
  summarise (seeds_germ = sum(seeds_germ),
             cumulative = last(cumulative)) %>%
  merge (germ_AW) %>%
  select(species, code, incubator, petridish, seeds_germ, germ_AW, viable, viablePER, germPER)  -> coldgerm_F

read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))  %>%
  spread(date, germinated, fill = 0) %>% # wide format for dates, and fill Na with 0
  gather ("date", "germinated", 8: last_col()) %>% # back in long format frrom 8th colum to the last
  arrange (species, accession, code, incubator, petridish, date)%>% # sort row observations this way
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>% 
  group_by (species, code, incubator, petridish, time) %>% 
  summarise(seeds_germ = sum(germinated)) %>%
  mutate(cumulative = cumsum(seeds_germ)) %>%
  filter (incubator == "Snowbed") %>% 
  filter (between(time, 301, 322)) %>% ## amount of days calculated from the dates 301 = 27/05 and 322 = 17/06
  group_by (species, code, incubator, petridish) %>%
  summarise (seeds_germ = sum(seeds_germ),
             cumulative = last(cumulative)) %>%
  merge (germ_AW) %>%
  select(species, code, incubator, petridish, seeds_germ, germ_AW, viable, viablePER, germPER) -> coldgerm_S

rbind(coldgerm_F, coldgerm_S) %>%
  arrange(species) %>%
  filter(viablePER>25)%>%
  filter(germPER>0)%>%
  mutate(species= str_replace(species, "Cerastium sp", "Cerastium pumilum"))%>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>% 
  merge(species, by = c("code", "species")) %>%
  mutate(code=factor(code)) %>%
  mutate(species=factor(species)) %>%
  mutate(incubator=factor(incubator)) %>%
  mutate(family=factor(family)) %>%
  mutate(macroclimate=factor(macroclimate)) %>%
  mutate(habitat=factor(habitat)) %>%
  mutate(germ_strategy=factor(germ_strategy)) %>%
  arrange(species, code, accession, incubator, petridish)  %>%
  mutate(ID = gsub(" ", "_", species), animal = ID) %>% 
  select(!family) %>%  
  na.omit ()-> df
summary(df)
df %>%
  select (code, species, incubator, petridish, seeds_germ, viable, macroclimate)%>%
  group_by (macroclimate, incubator) %>%
  summarise (germ = (seeds_germ/viable)*100) %>%
  get_summary_stats(type ="full")
#### WARM GERMINATION (BETWEEN Tmean>10 until end of experiment!) ####
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  #filter (species == "Carex sempervirens") %>%
  spread(date, germinated, fill = 0) %>% # wide format for dates, and fill Na with 0
  gather ("date", "germinated", 8: last_col()) %>% # back in long format frrom 8th colum to the last
  arrange (species, accession, code, incubator, petridish, date)%>% # sort row observations this way
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>% 
  group_by (species, code, incubator, petridish, time) %>% 
  summarise(seeds_germ = sum(germinated)) %>%
  mutate(cumulative = cumsum(seeds_germ)) %>%
  filter (incubator == "Fellfield") %>% 
  filter (between(time, 300, 430)) %>% ## amount of days calculated from the dates 300= 26/05 and 430 = 19/09
  summarise (seeds_germ = sum(seeds_germ),
             cumulative = last(cumulative)) %>%
  merge (germ_AW) %>%
  select(species, code, incubator, petridish, seeds_germ, germ_AW, viable, viablePER, germPER)  -> warmgerm_F

read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  #filter (species == "Carex sempervirens") %>%
  spread(date, germinated, fill = 0) %>% # wide format for dates, and fill Na with 0
  gather ("date", "germinated", 8: last_col()) %>% # back in long format frrom 8th colum to the last
  arrange (species, accession, code, incubator, petridish, date)%>% # sort row observations this way
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>% 
  group_by (species, code, incubator, petridish, time) %>% 
  summarise(seeds_germ = sum(germinated)) %>%
  mutate(cumulative = cumsum(seeds_germ)) %>% 
  filter (incubator == "Snowbed") %>% 
  filter (between(time, 323, 430)) %>% ## amount of days calculated from the dates 323 = 18/06 and 430 = 19/09
  group_by (species, code, incubator, petridish) %>%
  summarise (seeds_germ = sum(seeds_germ),
             cumulative = last(cumulative)) %>%
  merge (germ_AW) %>%
  select(species, code, incubator, petridish, seeds_germ, germ_AW, viable, viablePER, germPER) -> warmgerm_S

rbind(warmgerm_F, warmgerm_S) %>%
  arrange(species) %>%
  filter(viablePER>25)%>%
  filter(germPER>0)%>%
  mutate(species= str_replace(species, "Cerastium sp", "Cerastium pumilum"))%>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>% 
  merge(species, by = c("code", "species")) %>%
  mutate(code=factor(code)) %>%
  mutate(species=factor(species)) %>%
  mutate(incubator=factor(incubator)) %>%
  mutate(family=factor(family)) %>%
  mutate(macroclimate=factor(macroclimate)) %>%
  mutate(habitat=factor(habitat)) %>%
  mutate(germ_strategy=factor(germ_strategy)) %>%
  arrange(species, code, accession, incubator, petridish)  %>%
  mutate(ID = gsub(" ", "_", species), animal = ID) %>% 
  select(!family) %>%  
  na.omit ()-> df
summary(df)
df %>%
  select (code, species, incubator, petridish, seeds_germ, viable, macroclimate)%>%
  group_by (macroclimate, incubator) %>%
  summarise (germ = (seeds_germ/viable)*100) %>%
  get_summary_stats(type ="full")
## Germination timing table traits ####
germ_timing_list = list (viables, finalgerm, Mid_november, Wintergerm, Springgerm, Summergerm, Mid_june,  Mid_september)
germ_timing_list %>% 
  reduce(full_join, by = c("species", "code", "incubator", "petridish")) %>% 
#  select(species, code, incubator, petridish, total, viable, viablePER, germPER, 
#         Mid_nov, Wintergerm, Spring_germ, Mid_june, Summer_germ, Mid_sep)%>%
  mutate (Mid_nov = replace_na(Mid_nov, 100)) %>%
  mutate (Winter_germ = replace_na(Winter_germ, 0)) %>% 
  mutate (Mid_june= replace_na(Mid_june, 0)) %>% 
  mutate (Mid_sep = replace_na(Mid_sep, 0)) %>% 
  mutate (spring_germ= replace_na(spring_germ, 0))%>% 
  mutate (summer_germ = replace_na(summer_germ, 0))-> germ_timing 

germ_timing %>% 
  group_by (species, code, incubator) %>%
  summarise (total = sum(total),
             viable = sum (viable), 
             viablePER = mean(viablePER), 
             germPER = mean(germPER),
             Mid_nov = mean (Mid_nov), 
             Winter_germ = mean (Winter_germ), 
             spring_germ = mean (spring_germ), 
             summer_germ  = mean (summer_germ),
             Mid_june = mean (Mid_june),
             Mid_sep = mean (Mid_sep))%>%
  write.csv ("results/germination_timing2.csv")

### VISUALIZATION germination timing following Edu's instructions #####
# mid-november
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  group_by(region, species, code, incubator, petridish, date) %>%
  summarise(seeds_germ = sum(germinated)) %>%
  mutate(cum_seeds_germ = cumsum(seeds_germ)) %>%
  filter (date %in% c("2021-11-05","2021-11-08","2021-11-09","2021-11-10", "2021-11-11" )) %>% # select date with any of these values
  select(region, species, code, incubator, petridish, cum_seeds_germ) %>%
  group_by (region, incubator)%>%
  summarise (cum_seeds_germ = sum(cum_seeds_germ)) %>%
  merge(viables_region) %>%
  mutate (binom.confint(cum_seeds_germ, viable, methods = "wilson")) %>% 
  mutate (season = "Mid-November") %>% 
  select (region, incubator, cum_seeds_germ, viable, mean, lower, upper, season) -> Mid_november_graph
#rm(Mid_november_graph)
# winter
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))  %>%
  group_by (region, species, code, incubator, petridish, date)%>% 
  filter (date %in% c("2022-04-04","2022-04-11", "2022-04-12", # # first check after winter with T>2ºC
                      "2022-05-11","2022-05-12", "2022-05-18", "2022-05-19", "2022-05-25", "2022-05-26" )) %>% #undersnow + first checks T>2ºC
  group_by (region, incubator) %>% 
  summarise(cum_seeds_germ = sum(germinated)) %>%
  merge (viables_region) %>%
  mutate (binom.confint(cum_seeds_germ, viable, methods = "wilson")) %>% 
  mutate (season = "Winter") %>% 
  select (region, incubator, cum_seeds_germ, viable, mean, lower, upper, season) -> Winter_graph
#rm(Winter_graph)
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))  %>%
  group_by (region, species, code, incubator, petridish, date)%>% 
  filter (date %in% c("2022-04-04","2022-04-11", "2022-04-12", # # first check after winter with T>2ºC
                      "2022-05-11","2022-05-12", "2022-05-18", "2022-05-19", "2022-05-25", "2022-05-26" )) %>% #undersnow + first checks T>2ºC
  group_by (region, species, code, incubator, petridish) %>% 
  summarise(winter_seeds = sum(germinated)) -> winter_seeds ### necesary for mid-june graph calculation
#rm(winter_seeds)
# Mid_june
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))  %>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  group_by (region, species, code, incubator, petridish, time)%>% 
  filter (between(time, 106, 322)) %>% # amount of days calculated from the dates 106 = 13/11 and 322 = 16/6
  group_by (region, species, code, incubator, petridish) %>%
  summarise(seeds_germ = sum(germinated)) %>%
  merge(winter_seeds) %>%
  mutate (spring_germ = seeds_germ - winter_seeds) %>%
  group_by (region, incubator)%>%
  summarise (cum_seeds_germ = sum(spring_germ)) %>%
  merge(viables_region)  %>%
  mutate (binom.confint(cum_seeds_germ, viable, methods = "wilson")) %>% 
  mutate (season = "Mid-June") %>% 
  select (region, incubator, cum_seeds_germ, viable, mean, lower, upper, season) -> Mid_june_graph
#rm (Mid_june_graph)

# Mid-september
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))  %>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  group_by (region, species, code, incubator, petridish, time)%>% 
  filter (between(time, 322, 430)) %>% # amount of days calculated from the dates 322 = 16/6 and 430 = 19/09
  group_by (region,incubator) %>%
  summarise(cum_seeds_germ = sum(germinated)) %>%
  merge (viables_region) %>%
  mutate (binom.confint(cum_seeds_germ, viable, methods = "wilson")) %>% 
  mutate (season = "Mid-September") %>% 
  select (region, incubator, cum_seeds_germ, viable, mean, lower, upper, season) -> Mid_september_graph
#rm (Mid_september_graph)
# Spring
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))  %>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  group_by (species, code, incubator, petridish, time) %>% 
  filter (incubator == "Fellfield") %>% 
  filter (between(time, 262, 305)) %>% ## amount of days calculated from the dates 262 = 18/5 and 305 = 31/05
  group_by (region, incubator) %>%
  summarise(cum_seeds_germ = sum(germinated)) %>%
  merge (viables_region) %>%
  mutate (binom.confint(cum_seeds_germ, viable, methods = "wilson")) %>% 
  mutate (season = "Spring") %>% 
  select (region, incubator, cum_seeds_germ, viable, mean, lower, upper, season) -> Springgraph_F

read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))  %>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  group_by (species, code, incubator, petridish, time) %>% 
  filter (incubator == "Snowbed") %>% 
  filter (between(time, 301, 335)) %>% ## amount of days calculated from the dates 301 = 27/5 and 355 = 30/06
  group_by (region, incubator) %>%
  summarise(cum_seeds_germ = sum(germinated)) %>%
  merge (viables_region) %>%
  mutate (binom.confint(cum_seeds_germ, viable, methods = "wilson")) %>% 
  mutate (season = "Spring") %>% 
  select (region, incubator, cum_seeds_germ, viable, mean, lower, upper, season) -> Springgraph_S

Spring_graph <- rbind(Springgraph_F, Springgraph_S)
#rm(Spring_graph,Springgraph_F, Springgraph_S )
# summer
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))  %>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  group_by (species, code, incubator, petridish, time) %>% 
  filter (incubator == "Fellfield") %>% 
  filter (between(time, 305, 395)) %>% ## amount of days calculated from the dates 305 = 31/05 and 395 = 29/08
  group_by (region, incubator) %>%
  summarise(cum_seeds_germ = sum(germinated)) %>%
  merge (viables_region) %>%
  mutate (binom.confint(cum_seeds_germ, viable, methods = "wilson")) %>% 
  mutate (season = "Summer") %>% 
  select (region, incubator, cum_seeds_germ, viable, mean, lower, upper, season) -> Summergraph_F

read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))  %>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  group_by (species, code, incubator, petridish, time) %>% 
  filter (incubator == "Snowbed") %>% 
  filter (between(time, 335, 395)) %>% ## amount of days calculated from the dates 355 = 30/06 and 395 = 29/08
  group_by (region, incubator) %>%
  summarise(cum_seeds_germ = sum(germinated)) %>%
  merge (viables_region) %>%
  mutate (binom.confint(cum_seeds_germ, viable, methods = "wilson")) %>% 
  mutate (season = "Summer") %>% 
  select (region, incubator, cum_seeds_germ, viable, mean, lower, upper, season)  -> Summergraph_S

Summer_graph <- rbind(Summergraph_F, Summergraph_S)
#rm (Summer_graph, Summergraph_F, Summergraph_S)
# date_graph put together
graph_date <- rbind(Mid_november_graph, Winter_graph, Mid_june_graph, Mid_september_graph) %>% 
  mutate (season = factor(season, levels = c( "Mid-November", "Winter", "Mid-June", "Mid-September"))) #reorder levels
#rm(graph_date)
x11()
ggplot()+
  geom_point(data=graph_date, aes(season, mean, color=incubator), size =2) +
  facet_grid (.~region) +
  geom_errorbar(data= graph_date, aes(season, mean, ymin = lower, ymax = upper, color = incubator), size =1) +
  geom_line(data= graph_date, aes(as.numeric(season), mean, color = incubator), size =1) +
  theme_bw()

# season_graph put together
graph_season <- rbind(Mid_november_graph, Winter_graph, Spring_graph, Summer_graph) %>% 
    mutate (season = factor(season, levels = c( "Mid-November", "Winter", "Spring", "Summer"))) #reorder levels
  x11()
ggplot()+
    geom_point(data=graph_season, aes(season, mean, color=incubator), size =2) +
    facet_grid (.~region) +
    geom_errorbar(data= graph_season, aes(season, mean, ymin = lower, ymax = upper, color = incubator), size =1) +
    geom_line(data= graph_season, aes(as.numeric(season), mean, color = incubator), size =1) +
    theme_bw()
#rm (graph_season)
### delay to reach t50 check days between incubators ####
t50model %>%
  merge (viables_germ)%>%
  filter(viablePER>25)%>%
  group_by (species, code, incubator) %>%
  summarise(t50lm = mean(t50lm)) %>%
  spread(incubator, t50lm) %>% 
  mutate(delayS_F = Snowbed - Fellfield) %>%  # NAs appear when species don't reach 50% germination (t50lm_days =Na)
  na.omit () %>%
  mutate(species= str_replace(species, "Cerastium sp", "Cerastium pumilum"))%>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>% 
  merge(species, by = c("code", "species")) %>%
  mutate(code=factor(code)) %>%
  mutate(species=factor(species)) %>%
  mutate(family=factor(family)) %>%
  mutate(macroclimate=factor(macroclimate)) %>%
  mutate(habitat=factor(habitat)) %>%
  mutate(germ_strategy=factor(germ_strategy)) %>%
  mutate(ID = gsub(" ", "_", species), animal = ID) %>% 
  select(!family) %>%  
  na.omit () -> df 
summary(df)
df %>%
  select (code, species,delayS_F, macroclimate)%>%
  group_by (macroclimate) %>%
  summarise (Delay = delayS_F) %>%
  get_summary_stats(type ="full")
# delay to reach t50 check days between incubators ####
t50model %>%
  ungroup() %>%
  #select(species, code, incubator, t50lm) %>% 
  group_by (species, code, incubator) %>%# only possible if we join data by species and incubator, adding petridish produce an error
  summarise(t50lm = mean(t50lm)) %>%
  spread(incubator, t50lm) %>% 
  mutate(delayS_F = Snowbed - Fellfield) %>% 
  select (species, code, delayS_F)-> delaytime # NAs appear when species don't reach 50% germination (t50lm_days =Na)


#############################  DORMANCY ##############################
## Response to cold stratification ##
# compare autumn germination to after winter germination (winter+spring+summer)
germ_timing %>%
  mutate (coldstratgerm = ((germPER - Mid_nov)/germPER)*100, 
          coldstratgerm = round (coldstratgerm, digit=2)) %>%
  mutate (coldstratgerm = replace_na(coldstratgerm, 0)) %>%
  select (species, code, incubator, petridish, coldstratgerm) -> cold_strat
  

######## MOVE-ALONG TRAITS #########
germ_traits_list = list (germ_timing, t50_dates, heat_sum, cold_strat) 
germ_traits_list %>% 
  reduce(full_join, by = c("species", "code", "incubator", "petridish")) %>% 
  write.csv ("results/move_along_traits2.csv")
############################ GERMINAR ###############################
# Germination indices:
  # (1) Germination percentage
  # GRS - number of germinated seeds
  # GRP - germination percentage (from 0 to 100%)
  
  # (2) Germination speed
  # MGT - mean germination time (time units)
  # MGR - mean germination rate (time units)
  # GSP - germination speed (%)
  # NB! MGT, MGR and GSP tightly correlate to each other, since they express the very same process
  
  # (3) Germination synchrony 
  # UNC - uncertanity index (bits)
  # SYN - syncronization index (from 0 to 1)
  # NB! UNC and SYN tightly correlate to each other, since they express the very same process
  
  # 'auxiliary' indices - additional indices of germination process
  # VGT - germination variance
  # SDG - germination standard deviation
  # CVG - coefficient of variation

#  We will focus on uncertainty index, we need to create a function to do the same procedure for each species #  
  
UNC <- function (germ_ind) {
  read.csv("data/clean data.csv", sep = ";") %>%
    mutate(date = strptime(as.character(date), "%d/%m/%Y")) %>%
    pull(species) %>%
    unique() -> sp
  
  read.csv("data/clean data.csv", sep = ";") %>%
    mutate(date = strptime(as.character(date), "%d/%m/%Y")) %>% 
    filter (species == sp) %>% 
    group_by (region, species, code, incubator, petridish) %>%
    mutate(days = difftime(date, min(date), units = "days"), 
           days = round (days,digit =0)) %>% # calculate time from sowing to x date!
    select (region, species, code, incubator, petridish, total, viable, germinated, days) %>%
    rename (seeds = viable) %>%
    mutate (days = paste("D", days, sep = "")) %>%
    spread(days, germinated) %>%
    data.frame () %>%
    replace(is.na(.),0) -> dat

 ger_summary(SeedN = "seeds", evalName = "D", data=dat[,1: ncol(dat)])%>%
      select(region, species, code, incubator, petridish, unc)  
    
}

  read.csv("data/clean data.csv", sep = ";") %>%
    #mutate(date = strptime(as.character(date), "%d/%m/%Y")) %>%
    #group_by (species, code, incubator, petridish) %>%
    do (UNC(.))-> uncertainty
    
read.csv("data/clean data.csv", sep = ";") %>%
    mutate(date = strptime(as.character(date), "%d/%m/%Y")) %>%
    filter (species == "Sedum anglicum") %>% ## error en row id 250 y 254 info repetida ARREGLAR
    mutate(days = difftime(date, min(date), units = "days"), 
           days = round (days,digit =0)) %>% # calculate time from sowing to x date!
    select (species, region, code, incubator, petridish, total, viable, germinated, days) %>%
    group_by (region, species, code, incubator, petridish) %>%
    rename (seeds = viable) %>%
    mutate (days = paste("D", days, sep = "")) %>%
    spread(days, germinated)%>%
    data.frame () %>%
    replace(is.na(.),0) -> dat
ger_summary(SeedN = "seeds", evalName = "D", data=dat[,1: ncol(dat)])%>%
  mutate (unc = round(unc, digit=3))%>%
   select(region, species, code, incubator, petridish, unc) 
    
# calculate seed germination indices, SeedN - total number of VIABLE seeds in a replicate; D - columns with time units ('days')



  