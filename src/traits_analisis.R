library(FD); library(vegan); library(FactoMineR); library(emmeans);
library(tidyverse); library(ggrepel); library(cowplot); library (ggplot2);
library (lubridate); library(binom); library (GerminaR)
theme_set(theme_cowplot(font_size = 10)) 

#### dataframe with temperature programs x incubator####
read.csv("data/date_temp.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y")) %>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  mutate(incubator = as.factor(incubator)) %>% 
  as.data.frame()-> temp
str()
############################### THIS VERSION WITH 1ST + 2ND SOW DATA (ALL_DATA.CSV) #############################
# dataframe with species and accessions data ####
#read.csv ("data/species.csv", sep =";")-> sp_data
# rm(sp_data)
#read.csv ("data/accessions.csv", sep =";") -> accessions # for only 1st sow
read.csv ("data/all_code.csv", sep =";") -> codes # 1st + 2nd sow to calculate all traits in analysis? 
                                                  # only after winter / early season?
# viable seeds calculation ####
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  group_by(species, code, incubator, petridish) %>%
  filter(date == max(date)) %>%
  summarise(viable = sum(viable),
            total = sum(total)) %>%
  mutate (viablePER = (viable/total)*100,
          viablePER = round (viablePER, digit = 2)) -> viables

read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  merge(codes, by =c("code"))%>%
  group_by(region, species, code, incubator, petridish) %>%
  filter(date == max(date)) %>%
  select(region, species, code, incubator, petridish, viable, total) %>%
  group_by (region, incubator) %>%
  summarise(viable = sum(viable),
            total = sum(total)) -> viables_region
# final germination percentage calculation ####
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  group_by(species, code, incubator, petridish) %>%
  summarise(seeds_germ = sum(germinated)) %>% # 
  merge(viables) %>%
  mutate(germPER = (seeds_germ/viable) *100, # 
         germPER = round (germPER, digit =2)) %>%
  select (species, code, incubator, petridish, seeds_germ, germPER)-> finalgerm 

# seed germinated in spring + summer (needed for early season germ calculation) ####
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>% 
  group_by (species, code, incubator, petridish, time) %>%
  filter (incubator == "Fellfield") %>%
  filter (between(time, 257, 430)) %>%
  group_by(species, code, incubator, petridish) %>%
  summarise(germ_AW = sum(germinated)) -> germ_AF_F

read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>% 
  group_by (species, code, incubator, petridish, time) %>%
  filter (incubator == "Snowbed") %>%
  filter (between(time, 301, 430)) %>%
  group_by(species, code, incubator, petridish) %>%
  summarise(germ_AW = sum(germinated)) -> germ_AF_S


rbind(germ_AF_F, germ_AF_S)%>%
  arrange(species) -> germ_AF
# MID-NOVEMBER GERM (proxy of Autumn germination = before winter germination) ####
# last check 5-11/11/2021, Tmean < 3ºC  
read.csv("data/all_data.csv", sep = ";") %>%
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

# EARLY SEASON GERMINATION (when Tmax < 10 ) ####
read.csv("data/all_data.csv", sep = ";") %>%
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

read.csv("data/all_data.csv", sep = ";") %>%
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

# table for analysis and models ####
table_analisis_list = list (viables, finalgerm, Mid_november, Earlyseasongerm)
table_analisis_list %>% 
  reduce(full_join, by = c("species", "code", "incubator", "petridish"))%>% 
  select(species, code, incubator, petridish, viable, viablePER, germPER, Mid_nov, spring_germ) %>% 
  mutate (spring_germ= replace_na(spring_germ, 0)) -> dataframe # not sure about this line... what to do with species with no germ
