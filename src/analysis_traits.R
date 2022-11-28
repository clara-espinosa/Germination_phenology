library(FD); library(vegan); library(FactoMineR); library(emmeans);
library(tidyverse); library(ggrepel); library(cowplot); library (ggplot2);
library (lubridate)
theme_set(theme_cowplot(font_size = 10)) 

####dataframe with temperature programs x incubator####
temp <- read.csv("data/date_temp.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y")) %>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  mutate(incubator = as.factor(incubator)) 
temp <- as.data.frame(temp)

# viable seeds calculation ####
# summing up petridishes and accesions/populations of the same species (not taking into account weekly germination)
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date)))) %>%
  group_by(species, code, incubator, petridish) %>%
  filter(date == max(date)) %>%
  select(species, code, incubator, petridish, viable) %>%
  summarise(viable = sum(viable)) -> viables

# final germination percentage calculation ####
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date)))) %>%
  group_by(species, code, incubator, petridish) %>%
  summarise(seeds_germ = sum(germinated)) %>% # 
  merge(viables) %>%
  mutate(germPER = (seeds_germ/viable) *100, # 
         germPER = round (germPER, digit =2)) -> finalgerm 

# tidyverse modification to have time to reach 50% germination according to germ checks  ####
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
  mutate (t50 = ifelse(germPER >= 50, days, 0)) %>% # CAMBIAR 0 POR TIME MAX
  group_by (species, code, incubator, petridish) %>%
  mutate (t50check = first(t50[t50>0])) %>%
  summarise(t50check = min (t50check)) %>%
  ungroup()-> time50 # values discrete bc count t50 according to germination check days


### modeling t50 from raw data ####

f1 <- function(df0) {
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
  select(species, code, incubator, petridish, petricode, Timeframe, fg, t50times, t50g) %>%
  unique %>%
  group_by(species, code, incubator, petridish, petricode) %>%
  filter(fg >= .50) %>% # filter to keep only the species that reach more than 0.5 in max final germination
  filter(length(petricode) > 1) %>% # not sure why this comand about lenght of petricode variable, to fix some strange error?
  group_by(species, code, incubator, petridish) %>%
  do(f1(.)) %>%
  rename(#FG = fg, 
         t50lm = `(0.5 - as.numeric(mf1$coefficients[1]))/as.numeric(mf1$coefficients[2])`) ->t50model

t50_trait <- full_join(time50, t50model, by = c("species", "code", "incubator", "petridish")) %>%
  mutate (t50lm1 = t50lm) %>%
  mutate_at(vars(t50lm1), ~ as.integer(round(.x))) %>% # round decimals to integer numbers
  rename (t50lm_days = t50lm1) 
str(t50model)

# sowing date + t50 date or final date check ####
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  group_by (species, code, incubator, petridish) %>%
  mutate (datesow= first (date)) %>% 
  mutate (datelast= last (date)) %>% 
  select(species, code, incubator, petridish, viable, germinated, datesow, datelast) %>%
  merge (t50_trait) %>%
  group_by (species, code, incubator, petridish, datesow, datelast, t50lm_days)%>% 
  summarise(t50check = min (t50check)) %>% 
  mutate (date50 = as.Date(datesow) + t50lm_days) %>% 
  mutate (date50 = replace_na (date50, as.character(datelast))) -> t50_dates


## traits Fellfield-Snowbed dataframe version 1 ####
traits_FS1 <- full_join(t50_dates, finalgerm,  by= c("species", "code", "incubator", "petridish")) 

### Heat:sum function####
heatsum <-function (traits_FS1) {
  traits_FS1 %>%
    pull (incubator)%>%
    unique() -> inc
  traits_FS1 %>%
    pull(datesow) -> date1 #convierte columna en vector
  traits_FS1%>%
    pull(date50) -> date2
  temp %>% 
    mutate (date = as.Date(date)) %>%
    filter(incubator == inc) %>%
    filter(date >= date1 & date <= date2) %>%
    summarise (HS =sum(Tmean)) 
}

traits_FS1%>%
  group_by (species, code, incubator, petridish) %>%
  do (heatsum(.)) %>% # punto significa que el objeto al que aplicar la funcion heat sum es el de la linea de arriba
  data.frame -> heat_sum

## traits Fellfield-Snowbed dataframe version 2 ####
traits_FS2 <- full_join(traits_FS1, heat_sum, by= c("species", "code", "incubator", "petridish")) %>%
  ungroup()

# undersnow germination (nº of seeds) ####
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))  %>%
  filter (incubator=="Snowbed") %>%# only able to determine in the snowbed incubator (due to env conditions)
  group_by (species, code, incubator, petridish, date)%>% 
  filter (date %in% c("2022-05-11","2022-05-12")) %>% # date of germination check right before end of winter conditions (dark + 0ÂºC)
  mutate (snowgerm = germinated) %>%
  merge (finalgerm) %>%
  select(species, code, incubator, petridish, viable, seeds_germ, germPER, snowgerm) %>%
  mutate (snowPER = (snowgerm/seeds_germ)*100) %>% 
  select (species, code, incubator, petridish, snowgerm, snowPER) -> snow_trait  # percentage of germinated under snow 
## traits Fellfield-Snowbed dataframe version 3 ####
traits_FS3 <- full_join(traits_FS2, snow_trait, by= c("species", "code", "incubator", "petridish")) %>%
  select(species, code, incubator, petridish, viable, seeds_germ, germPER, snowgerm, snowPER, t50lm_days,t50check, HS)

# delay to reach t50 check days between incubators ####
traits_FS3 %>%
  ungroup() %>%
  select(species, code, incubator, t50lm_days) %>%# only posible if we join data by species and incubator, addind petridish produce an error
  group_by (species, code, incubator) %>%
  summarise(t50lm = mean(t50lm_days)) %>%
  spread(incubator, t50lm) %>% 
  mutate(delayS_F = Snowbed - Fellfield) %>% 
  select (species, code, delayS_F)-> delaytime 
 
snow_trait %>%
  group_by (species, code)%>%
  summarise(snowgerm = sum(snowgerm), 
            snowPER = mean(snowPER))  %>%
  merge(delaytime) -> traits_sp

finalgerm %>%
  group_by (species, code) %>%
  summarise(viable = sum(viable), 
            seeds_germ= sum(seeds_germ),
            germPER = mean(germPER))%>%
  merge(traits_sp) %>%
  write.csv("results/traits_sintesis_sp.csv")


##### weighted mean of germination timing ####
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y")) %>%
  filter (species == "Minuartia recurva") %>%
  select (species, code, incubator, petridish, viable, germinated, date)%>%
  group_by(species, code, incubator, petridish) %>%
  mutate(days = difftime(date, min(date), units = "days")) %>% # calculate time from sowing to x date!
  mutate(cs = cumsum(germinated), # cumulative sum of germinated seeds 
         fg = max(cs) / viable, # max germinated (final germination) divided by viable seeds (proportion)
         g = cs / viable)%>%
  ungroup()%>%
  group_by (species, code, incubator, petridish) %>%
  mutate (weightmean = weighted.mean (as.numeric (days), germinated)) %>%
  summarise (germinated = sum(germinated),
             cs = max (cs),
             fg = fg,
             weightmean = weightmean) -> weightmean
             
#### sintesis table ####
traits_FS3 %>%
  ungroup () %>%
  group_by(species, code, incubator) %>%
  summarise(viable = sum(viable), 
            seeds_germ= sum(seeds_germ),
            germPER = mean(germPER),
            snowgerm= sum(snowgerm), 
            snowPER= mean(snowPER),
            t50lm = mean(t50lm_days),
            t50check = mean (t50check), 
            Hs = mean(HS)) %>%
write.csv("results/traits_sintesis_incubator.csv")
