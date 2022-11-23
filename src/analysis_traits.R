library(FD); library(vegan); library(FactoMineR); library(emmeans);
library(tidyverse); library(ggrepel); library(cowplot); library (ggplot2);
library (lubridate)
theme_set(theme_cowplot(font_size = 10)) 

## dataframe for temperature programs ##
temp <- read.csv("data/date_temp.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y")) %>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  mutate(incubator = as.factor(incubator)) 
  #mutate(heatsum = cumsum(Tmean))
temp <- as.data.frame(temp)
str(temp)
summary(temp)

## main data transformation and visualization ####
# tidyverse transformation to account for the number of viable seeds per each specie and incubator
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

# tidyverse modification to have the accumulated germination along the whole experiment
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

# tidyverse modification to have time to reach 50% germination
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

# sowing date + t50 date
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

traits_df <- full_join(finalgerm, time50dates, by= c("species", "incubator")) 

### FUNCION FUNCIONA ####
# FUNCION BUENA
heatsum_bis <-function (traits_df) {
  traits_df %>%
    pull (incubator)%>%
    unique() -> inc
  
  traits_df %>%
    pull(datesow) -> date1 #convierte columna en vector
  
  traits_df%>%
    pull(date50) -> date2
  
  temp %>% 
    mutate (date = as.Date(date)) %>%
    filter(incubator == inc) %>%
    filter(date >= date1 & date <= date2) %>%
    summarise (HS =sum(Tmean)) 
}


traits_df %>%
  group_by (species, incubator) %>%
  do (heatsum_bis(.)) %>% # punto significa que el objeto al que aplicar la funcion heat sum es el de la linea de arriba
  data.frame -> heat_sum

traits_df <- full_join(traits_df, heat_sum, by= c("species", "incubator")) 

ggplot (traits_df, aes(HS, incubator, fill=HS)) +
  geom_boxplot () +
  theme_bw()
### how to model exact p50 from our data points?? lm? non linear model? 

### t50

f1 <- function(df0) {
  lm(t50g ~ t50times, data = df0) -> mf1 # Linear model between time before and time after
  as.data.frame((0.5 - as.numeric(mf1$`coefficients`[1])) / as.numeric(mf1$`coefficients`[2])) # Use linear model to interpolate t50
}

read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y")) %>%
  group_by(species, code, incubator, petridish) %>%
  mutate(days = difftime(date, min(date), units = "days")) %>%
  arrange(species, code, incubator, petridish, days) %>%
  group_by (species, code, incubator, petridish) %>% 
  mutate(cs = cumsum(germinated), 
         fg = max(cs) / viable,
         g = cs / viable,
         Timeframe = ifelse(g > 0.5, "Upper", "Lower")) %>%
  group_by(species, code, incubator, petridish, Timeframe) %>%
  mutate(t50times = ifelse(Timeframe == "Lower", max(days), min(days)),
         t50g = ifelse(days == t50times, g, NA)) %>%
  na.omit %>%
  select(species, code, incubator, petridish, petricode, Timeframe, fg, t50times, t50g) %>%
  unique %>%
  group_by(species, code, incubator, petridish, petricode) %>%
  filter(fg >= .50) %>%
  filter(length(petricode) > 1) %>%
  group_by(species, code, incubator, petridish, petricode, fg) %>%
  do(f1(.)) %>%
  rename(FG = fg, 
         t50 = `(0.5 - as.numeric(mf1$coefficients[1]))/as.numeric(mf1$coefficients[2])`)