library(FD); library(vegan); library(FactoMineR); library(emmeans);
library(tidyverse); library(ggrepel); library(cowplot); library (ggplot2);
library (lubridate)
theme_set(theme_cowplot(font_size = 10)) 

####dataframe with temperature programs x incubator####
temp <- read.csv("data/date_temp.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y")) %>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  mutate(incubator = as.factor(incubator)) 
  #mutate(heatsum = cumsum(Tmean))
temp <- as.data.frame(temp)
str(temp)
summary(temp)

# viable seeds calculation ####
# summing up petridishes and accesions/populations of the same species (not taking into account weekly germination)
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date)))) %>%
  group_by(species, incubator, code, petridish) %>%
  filter(date == max(date)) %>%
  select(species, incubator, code, petridish, viable) %>%
  group_by(species, incubator) %>%
  summarise(viable = sum(viable)) -> viables

ggplot (viables, aes(viable, incubator, fill = incubator)) +
  geom_boxplot () +
  theme_bw()

# final germination percentage calculation ####
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date)))) %>%
  group_by(species, incubator) %>%
  summarise(germinated = sum(germinated)) %>%
  merge(viables) %>%
  mutate(germination = (germinated/viable) *100) -> finalgerm 

ggplot (finalgerm, aes(germination, incubator, fill= incubator))+
  geom_boxplot ()+
  theme_bw ()

# tidyverse modification to have time to reach 50% germination according to germ checks  ####
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  group_by (species) %>%
  mutate (datesow= first (date)) %>% # a different column for sowing date
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date)))) %>% # time count from sowing of each sp
  ungroup () %>%
  group_by (species, incubator, time)%>% 
  summarise(germinated = sum(germinated)) %>%
  mutate(germination = cumsum(germinated)) %>%
  merge (viables) %>%
  mutate(germination = (germination/viable) *100) %>%
  ungroup() %>%
  mutate (t50 = ifelse(germination > 50, time, 0)) %>% # CAMBIAR 0 POR TIME MAX
  group_by (species, incubator) %>%
  mutate (t50check = first(t50[t50>0])) %>%
  summarise(t50check = min (t50check)) %>%
  ungroup()-> time50 # values discrete bc count t50 according to germination check days

ggplot (time50, aes(t50check, incubator, fill= incubator)) +
  geom_boxplot () +
  theme_bw()

# sowing date + t50 date or final date check ####
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  group_by (species, incubator) %>%
  mutate (datesow= first (date)) %>% 
  mutate (datelast= last(date)) %>%
  select(species, incubator, viable, germinated, datesow, datelast) %>%
  merge (time50) %>%
  group_by (species, incubator, datesow, datelast)%>% 
  summarise(t50check = min (t50check)) %>%
  mutate (date50 = as.Date(datesow) + t50check) %>%
  mutate (date50 = replace_na (date50, as.character(datelast))) -> time50dates

## traits Fellfield-Snowbed dataframe version 1 ####
traits_FS <- full_join(finalgerm, time50dates, by= c("species", "incubator")) 

### Heat:sum function####
heatsum <-function (traits_FS) {
  traits_FS %>%
    pull (incubator)%>%
    unique() -> inc
  traits_FS %>%
    pull(datesow) -> date1 #convierte columna en vector
  traits_FS%>%
    pull(date50) -> date2
  temp %>% 
    mutate (date = as.Date(date)) %>%
    filter(incubator == inc) %>%
    filter(date >= date1 & date <= date2) %>%
    summarise (HS =sum(Tmean)) 
}

traits_FS%>%
  group_by (species, incubator) %>%
  do (heatsum(.)) %>% # punto significa que el objeto al que aplicar la funcion heat sum es el de la linea de arriba
  data.frame -> heat_sum

## traits Fellfield-Snowbed dataframe version 2 ####
traits_FS <- full_join(traits_FS, heat_sum, by= c("species", "incubator")) 

ggplot (traits_df, aes(HS, incubator, fill=HS)) +
  geom_boxplot () +
  theme_bw()

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
         Timeframe = ifelse(g > 0.5, "Upper", "Lower")) %>% # # New column that divides dataset in germination < 0.5 or > 0.5
  group_by(species, code, incubator, petridish, Timeframe) %>%
  mutate(t50times = ifelse(Timeframe == "Lower", max(days), min(days)), # keep only the number of days just before and after reaching t50
         t50g = ifelse(days == t50times, g, NA)) %>% # copy of germination proportion the days just before and after t50
  na.omit %>%
  select(species, code, incubator, petridish, petricode, Timeframe, fg, t50times, t50g) %>%
  unique %>%
  group_by(species, code, incubator, petridish, petricode) %>%
  filter(fg >= .50) %>% # filter to keep only the species that reach more than 0.5 in max final germination
  filter(length(petricode) > 1) %>% # not sure why this comand about lenght of petricode variable, to fix some strange error?
  #group_by(species, code, incubator, petridish, petricode, fg) %>%
  group_by (species, incubator) %>%
  do(f1(.)) %>%
  rename(#FG = fg, 
         t50lm = `(0.5 - as.numeric(mf1$coefficients[1]))/as.numeric(mf1$coefficients[2])`) ->t50model

t50comparison <- full_join(time50, t50model, by = c ("species", "incubator"))

# undersnow germination (nÂº of seeds) ####
snow <- read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y")) 
snow <- filter (snow, incubator=="Snowbed") # only able to determine in the snowbed incubator (due to env conditions)
snow %>%
  group_by (species, incubator, date)%>% 
  summarise(germinated = sum(germinated)) %>%
  mutate(germination = cumsum(germinated)) %>%
  inner_join (viables) %>%
  filter (date %in% c("2022-05-11","2022-05-12")) %>% # date of germination check right before end of winter conditions (dark + 0ÂºC)
  mutate (snowgerm = germinated) %>%
  group_by(species) %>%
  select (species, incubator, snowgerm)  %>%
  merge (finalgerm) %>%
  select(species, snowgerm, germinated) %>%
  mutate (snowPER = (snowgerm/germinated)*100)-> snow_trait # percentage of germinated under snow 
snow_trait <- rename (snow_trait, germinated_snowbed = germinated)  
#### traits sp dataframe  (average/combination Fellfield + Snowbed)####
traits_FS %>%
  group_by (species)%>%
  summarise (germinated =sum(germinated),
             viable = sum (viable),
             germination = mean (germination),
             t50check= mean (t50check),
             Heatsum= mean (HS)) -> traits_sp
traits_sp <- full_join(traits_sp, snow_trait, by= "species")

# delay to reach t50 check days between incubators ####
traits_FS %>%
  select(species, incubator, t50check)-> delaytime
spread(delaytime, incubator, t50check)-> delaytime 
delaytime %>% 
  mutate(delayF_S = Snowbed - Fellfield) %>% 
  select (species, delayF_S)-> delaytime

traits_sp <- full_join(traits_sp, delaytime, by= "species") 
# dormancy/cold stratification needed + warm cues?####