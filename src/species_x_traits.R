library(tidyverse)

#species data
read.csv("data/all_info.csv", sep =";") %>% # modify column values of 3 sp names:
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>% # Minuartia CF as Minuartia arctica;
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album"))-> species  #Sedum sp (album cf) as Sedum album

  

# final germination /petridish calculation ####
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  merge(read.csv("data/all_info.csv", sep =";"))%>%
  group_by(species, macroclimate, code, incubator, petridish) %>%
  summarise(seeds_germ = sum(germinated)) %>% # 
  merge(viables) %>%
  mutate(germPER = (seeds_germ/viable) *100, # 
         germPER = round (germPER, digit =2)) %>%
  group_by(species, macroclimate, incubator) %>%
  summarise(seeds_germ = sum(seeds_germ),
            viable = sum(viable),
            viablePER = mean(viablePER),
            germPER = mean (germPER))-> vgerm
         
# AUTUMN (Mid November) ####
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  merge(read.csv("data/all_info.csv", sep =";"))%>%
  group_by (species, macroclimate, code, incubator, petridish, time) %>%
  summarise(seeds_germ = sum(germinated)) %>%
  mutate(cumulative = cumsum(seeds_germ)) %>%
  filter (between(time, 1, 105)) %>% ## 105 = 12/11 last check before winter
  group_by(species, macroclimate, incubator) %>%
  summarise (autumn_germ = sum(seeds_germ)) -> autumn


# SPRING (Mid June)  ####
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  spread(date, germinated, fill = 0) %>% # wide format for dates, and fill Na with 0
  gather ("date", "germinated", 8: last_col() )%>% # back in long format frrom 8th colum to the last
  arrange (species, accession, code, incubator, petridish, date)%>% # sort row observations this way
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  merge(read.csv("data/all_info.csv", sep =";"))%>%
  group_by (species, macroclimate, code, incubator, petridish, time)%>% 
  summarise(seeds_germ = sum(germinated)) %>%
  mutate(cumulative = cumsum(seeds_germ)) %>%
  filter (between(time, 106, 322)) %>% # amount of days calculated from the dates 106 = 13/11 and 322 = 16/6
  ungroup () %>%
  group_by(species, macroclimate, incubator) %>%
  summarise (spring_germ = sum(seeds_germ)) -> spring

# SUMMER (Mid September) ####
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))  %>%
  spread(date, germinated, fill = 0) %>% # wide format for dates, and fill Na with 0
  gather ("date", "germinated", 8: last_col() )%>% # back in long format frrom 8th colum to the last
  arrange (species, accession, code, incubator, petridish, date)%>% # sort row observations this way
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  merge(read.csv("data/all_info.csv", sep =";"))%>%
  group_by (species, macroclimate, code, incubator, petridish, time)%>% 
  summarise(seeds_germ = sum(germinated)) %>%
  mutate(cumulative = cumsum(seeds_germ)) %>%
  filter (between(time, 323, 430)) %>% # amount of days calculated from the dates 323 = 17/6 and 430 = 19/09
  ungroup () %>%
  group_by(species, macroclimate, incubator) %>%
  summarise (summer_germ = sum(seeds_germ)) -> summer

# winter germ ####
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))  %>%
  spread(date, germinated, fill = 0) %>% # wide format for dates, and fill Na with 0
  gather ("date", "germinated", 8: last_col()) %>% # back in long format frrom 8th colum to the last
  arrange (species, accession, code, incubator, petridish, date)%>% # sort row observations this way
  merge(read.csv("data/all_info.csv", sep =";"))%>%
  group_by (species, macroclimate, code, incubator, petridish, date)%>% 
  summarise(seeds_germ = sum(germinated)) %>%
  mutate(cumulative = cumsum(seeds_germ)) %>%
  filter (date %in% c("2022-04-04", # # first check after winter before Tmin 2ºC
                      "2022-05-11","2022-05-12", "2022-05-18", "2022-05-19", 
                      "2022-05-25", "2022-05-26")) %>% #undersnow + first checks Tmin 2ºC
  ungroup () %>%
  group_by(species, macroclimate, incubator) %>%
  summarise (winter_germ = sum(seeds_germ)) -> winter

# t50 model
t50model %>%
  merge(read.csv("data/all_info.csv", sep =";"))%>%
  group_by(species, macroclimate, incubator) %>%
  summarise (t50 = mean(t50lm)) -> t50

# heat sum
t50_dates%>%
  merge(read.csv("data/all_info.csv", sep =";"))%>%
  group_by (species, macroclimate, code, incubator, petridish) %>%
  do (heatsum(.)) %>%
  group_by(species, macroclimate, incubator) %>%
  summarise(heat_sum = mean(HS)) -> heat_sum

# seed mass
read.csv("data/seed_mass.csv", sep = ";") %>%
  group_by(species, macroclimate) %>%
  summarise(mass_50 = mean(mass_50)) %>%
  filter (!(species == "Avenella flexuosa")) %>%  # filter species that are not in the move-along
  filter (!(species == "Festuca glacialis")) %>%
  filter (!(species == "Galium pyrenaicum")) %>%
  filter (!(species == "Jasione laevis")) %>%
  filter (!(species == "Linaria alpina")) %>%
  filter (!(species == "Saxifraga oppositifolia")) %>%
  filter (!(species == "Sedum atratum")) %>%
  filter (!(species == "Sempervivum vicentei")) %>%
  filter (!(species == "Teesdalia conferta")) %>%
  filter (!(species == "Veronica nummularia")) -> seed_mass

setdiff(vgerm$species, seed_mass$species)  # 7 sp from move along so not have seed mass own data (TRY database?)
setdiff(vgerm$species, t50$species)  # 7 sp from move along so not have seed mass own data (TRY database?)
# traits table x PCA WORK IN PROGRESS, CONTINUE HERE ON MONDAY!!!!!!
germ_traits_list = list (vgerm, autumn, spring, summer, winter, t50, heat_sum) #, delay

germ_traits_list %>% 
  reduce(full_join, by = c("species", "macroclimate", "incubator")) %>% 
  mutate (autumn_germ = autumn_germ/viable,
          spring_germ = spring_germ/viable,
          summer_germ = spring_germ/viable, 
          winter_germ = winter_germ/viable) %>% 
  full_join(seed_mass, by= c("species", "macroclimate")) %>% 
  merge(read.csv("data/species.csv", sep =";")) %>% 
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album"))%>%
  write.csv("data/traits_sp_incubator.csv")



######### EXTRA traits ###############
# cold germ 
rbind(coldgerm_F, coldgerm_S) %>%
  group_by (species) %>%
  summarise (cold_germ = sum(seeds_germ)) -> cold

# warm germ
rbind(warmgerm_F, warmgerm_S) %>%
  group_by (species) %>%
  summarise (warm_germ = sum(seeds_germ)) -> warm 

  # delay 
t50model %>%
  group_by (species, code, incubator) %>%
  summarise(t50lm = mean(t50lm)) %>%
  spread(incubator, t50lm) %>% 
  mutate (delayS_F = Snowbed - Fellfield) %>%
  na.omit() %>%
  group_by(species) %>%
  summarise (delayS_F = mean(delayS_F))-> delay   # NAs appear when species don't reach 50% germination (t50lm_days =Na)



  
