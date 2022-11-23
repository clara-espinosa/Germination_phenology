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
#FUNCION de EDUARDO

f1 <- function(df0) {
  df0 %>% # Data frame with germination scoring
    mutate(t = T / PG) %>% # Calculate cumulative germination proportion
    mutate(Timeframe = ifelse(t > 0.5, "Upper", "Lower")) %>% # New column that divides dataset in germination < 0.5 or > 0.5
    group_by(Accession, Treatment, Timeframe) -> dff1 # Group by accession, treatment, timeframe
  rbind(
    filter(dff1, Timeframe == "Lower") %>% do(tail(., n = 1)), # Scoring time just before t50
    filter(dff1, Timeframe == "Upper") %>% do(head(., n = 1)) # Scoring time just after t50
  ) -> dff2
  lm(t ~ Time, data = dff2) -> mf1 # Linear model between time before and time after
  (0.5 - as.numeric(mf1$`coefficients`[1])) / as.numeric(mf1$`coefficients`[2]) # Use linear model to interpolate t50
}

f1 <- function(df0) {
    read.csv("data/clean data.csv", sep = ";") %>% # Data frame with germination scoring
    mutate(t = germinated / viable) %>% # Calculate cumulative germination proportion
    mutate(Timeframe = ifelse(t > 0.5, "Upper", "Lower")) %>% # New column that divides dataset in germination < 0.5 or > 0.5
    group_by(species, incubator, Timeframe)-> dff1 # Group by accession, treatment, timeframe
  rbind(
    filter(dff1, Timeframe == "Lower") %>% do(tail(., n = 1)), # Scoring time just before t50
    filter(dff1, Timeframe == "Upper") %>% do(head(., n = 1)) # Scoring time just after t50
  ) -> dff2
  lm(t ~ date, data = dff2) -> mf1 # Linear model between time before and time after
  as.data.frame((0.5 - as.numeric(mf1$`coefficients`[1])) / as.numeric(mf1$`coefficients`[2])) # Use linear model to interpolate t50
}

read.csv("data/clean data.csv", sep = ";") %>% # Data frame with germination scoring
  mutate(t = germinated / viable) %>% # Calculate cumulative germination proportion
  mutate(Timeframe = ifelse(t > 0.5, "Upper", "Lower")) %>% # New column that divides dataset in germination < 0.5 or > 0.5
  group_by(species, incubator, Timeframe)-> dff1 # Group by accession, treatment, timeframe
rbind(
  filter(dff1, Timeframe == "Lower") %>% do(tail(., n = 1)), # Scoring time just before t50
  filter(dff1, Timeframe == "Upper") %>% do(head(., n = 1)) # Scoring time just after t50
) -> dff2
dff2 %>%
  arrange (species, incubator, Timeframe)
lm(t ~ date, data = dff2) -> mf1 # Linear model between time before and time after
(0.5 - as.numeric(mf1$`coefficients`[1])) / as.numeric(mf1$`coefficients`[2]) # Use linear model to interpolate t50

read.csv("data/clean data.csv", sep = ";") %>%
  group_by (species, incubator) %>%
  do (f1(.)) %>% # punto significa que el objeto al que aplicar la funcion heat sum es el de la linea de arriba
  data.frame -> time50lm
# seguimos teniendo el problema de aquellas especies que no han llegado al 50% de germinación. 
f1 <- function(df0) {
  lm(t50g ~ t50times, data = df0) -> mf1 # Linear model between time before and time after
  as.data.frame((0.5 - as.numeric(mf1$`coefficients`[1])) / as.numeric(mf1$`coefficients`[2])) # Use linear model to interpolate t50
}

