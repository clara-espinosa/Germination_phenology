library(FD); library(vegan); library(FactoMineR); library(emmeans);
library(tidyverse); library(ggrepel); library(cowplot); library (ggplot2);
library (lubridate); library(binom); library (GerminaR)
theme_set(theme_cowplot(font_size = 10)) 

#### dataframe with temperature programs x incubator####
read.csv("data/date_temp.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y")) %>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  filter (between(time, 83, 430)) %>% #keep data only between start of 2n sow
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  mutate(incubator = as.factor(incubator)) %>% 
  as.data.frame()-> temp2
str()
# dataframe with species and accessions data ####
#read.csv ("data/species.csv", sep =";")-> sp_data
# rm(sp_data)
read.csv ("data/accessions 2nd sow.csv", sep =";") -> accessions2

# viable seeds calculation ####
read.csv("data/clean data 2n sow.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  group_by(species, code, incubator, petridish) %>%
  filter(date == max(date)) %>%
  summarise(viable = sum(viable),
            total = sum(total)) %>%
  mutate (viablePER = (viable/total)*100,
          viablePER = round (viablePER, digit = 2)) -> viables2

read.csv("data/clean data 2n sow.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  merge(accessions2, by =c("code", "accession"))%>%
  group_by(region, species, code, incubator, petridish) %>%
  filter(date == max(date)) %>%
  select(region, species, code, incubator, petridish, viable, total) %>%
  group_by (region, incubator) %>%
  summarise(viable = sum(viable),
            total = sum(total)) -> viables_region2
# final germination percentage calculation ####
read.csv("data/clean data 2n sow.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  group_by(species, code, incubator, petridish) %>%
  summarise(seeds_germ = sum(germinated)) %>% # 
  merge(viables2) %>%
  mutate(germPER = (seeds_germ/viable) *100, # 
         germPER = round (germPER, digit =2)) %>%
  select (species, code, incubator, petridish, seeds_germ, germPER)-> finalgerm2 

# seed germinated in spring + summer (needed for early season germ calculation) ####
read.csv("data/clean data 2n sow.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>% 
  group_by (species, code, incubator, petridish, time) %>%
  filter (incubator == "Fellfield") %>%
  filter (between(time, 174, 347)) %>%
  group_by(species, code, incubator, petridish) %>%
  summarise(germ_AW = sum(germinated)) -> germ2_AF_F

read.csv("data/clean data 2n sow.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>% 
  group_by (species, code, incubator, petridish, time) %>%
  filter (incubator == "Snowbed") %>%
  filter (between(time, 218, 347)) %>%
  group_by(species, code, incubator, petridish) %>%
  summarise(germ_AW = sum(germinated)) -> germ2_AF_S

rbind(germ2_AF_F, germ2_AF_S)%>%
  arrange(species) -> germ2_AF
# MID-NOVEMBER GERM (proxy of Autumn germination = before winter germination) ####
# last check 5-11/11/2021, Tmean < 3ºC  
read.csv("data/clean data 2n sow.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  group_by (species, code, incubator, petridish, time)%>%
  filter (between(time, 1, 21)) %>%
  group_by (species, code, incubator, petridish) %>%
  summarise(seeds_germ = sum(germinated))%>%
  merge(viables2) %>%
  mutate(Mid_nov= (seeds_germ/viable) *100,  
         Mid_nov= round (Mid_nov, digit =2)) %>%
  select(species, code, incubator, petridish, Mid_nov) -> Mid_november2

# EARLY SEASON GERMINATION (when Tmax < 10 ) ####
read.csv("data/clean data 2n sow.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>% 
  group_by (species, code, incubator, petridish, time) %>% 
  filter (incubator == "Fellfield") %>% 
  filter (between(time, 174, 202)) %>% ## amount of days calculated from the dates 174 = 13/4 and 202 = 11/05
  group_by (species, code, incubator, petridish) %>%
  summarise(seeds_germ = sum(germinated)) %>%
  merge (germ2_AF_F) %>%
  mutate (spring_germ = (seeds_germ/germ_AW) *100,
          spring_germ = round (spring_germ, digit=2)) %>%
  mutate (spring_germ= replace_na(spring_germ, 0))%>% # not sure about this line... what to do with species with no germ
  select (species, code, incubator, petridish, spring_germ) -> Earlyseason2_F

read.csv("data/clean data 2n sow.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))  %>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  group_by (species, code, incubator, petridish, time) %>% 
  filter (incubator == "Snowbed") %>% 
  filter (between(time, 218, 239)) %>% ## amount of days calculated from the dates 218 = 27/5 and 239 = 17/06
  group_by (species, code, incubator, petridish) %>%
  summarise(seeds_germ = sum(germinated)) %>%
  merge (germ2_AF_S) %>%
  mutate (spring_germ = (seeds_germ/germ_AW) *100,
          spring_germ = round (spring_germ, digit=2)) %>%
  mutate (spring_germ= replace_na(spring_germ, 0))%>% # not sure about this line... what to do with species with no germ
  select (species, code, incubator, petridish, spring_germ)  -> Earlyseason2_S

rbind(Earlyseason2_F, Earlyseason2_S) %>%
  arrange(species) -> Earlyseasongerm2 
