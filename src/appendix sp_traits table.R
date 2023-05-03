library (tidyverse)

# appendic table with sp data + traits
read.csv ("data/all_info.csv", sep = ";") -> sp
# viables + final germ + filtering variables ####
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  group_by(species, accession, code, incubator, petridish) %>%
  summarise(total_germ = sum(germinated),
            viable= last(viable),
            total = last(total)) %>% 
  merge(species)%>% 
  group_by (community, species, accession, code, incubator, family, recolection) %>%
  summarise (total_germ = sum (total_germ), 
             viable = sum (viable), 
             total = sum(total), 
             abundance = mean (abundance))%>%
  mutate(germPER = (total_germ/viable) *100, # 
         germPER = round (germPER, digit =2)) %>%
  mutate(viablePER=(viable/total) *100, 
         viablePER = round(viablePER, digit =2)) %>%
  arrange(community, species, accession, code) -> appendix 
#autumn
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  group_by (species, accession, code, incubator, petridish, time) %>%
  summarise(seeds_germ = sum(germinated)) %>%
  mutate(cumulative = cumsum(seeds_germ)) %>%
  filter (between(time, 1, 105)) %>% ## 105 = 12/11 last check before winter
  group_by (species, accession, code, incubator) %>%
  summarise (seeds_germ = sum(seeds_germ)) %>%
  merge(appendix)%>%
  group_by (community, species, accession, code, incubator)%>%
  mutate(autumn_germ = (seeds_germ/viable)*100, 
         autumn_germ = round (autumn_germ, digit = 2))%>%
  select (community, species, abundance, family, accession, code,  recolection, incubator, 
           total,viable, viablePER, total_germ,germPER,  autumn_germ) -> appendix
#spring
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  spread(date, germinated, fill = 0) %>% # wide format for dates, and fill Na with 0
  gather ("date", "germinated", 8: last_col() )%>% # back in long format frrom 8th colum to the last
  arrange (species, accession, code, incubator, petridish, date)%>% # sort row observations this way
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  group_by (species, accession, code, incubator, petridish, time)%>% 
  summarise(seeds_germ = sum(germinated)) %>%
  mutate(cumulative = cumsum(seeds_germ)) %>%
  filter (between(time, 106, 322)) %>% # amount of days calculated from the dates 106 = 13/11 and 322 = 16/6
  ungroup () %>%
  group_by (species, accession, code, incubator) %>%
  summarise (seeds_germ = sum(seeds_germ)) %>%
  merge(appendix)%>%
  group_by (community, species, accession, code,incubator)%>%
  mutate(spring_germ = (seeds_germ/viable)*100, 
         spring_germ = round (spring_germ, digit = 2))%>%
  select (community, species, abundance, family, accession, code, recolection, incubator, 
          total,viable, viablePER, total_germ,germPER,  autumn_germ, spring_germ) -> appendix
#summer
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))  %>%
  spread(date, germinated, fill = 0) %>% # wide format for dates, and fill Na with 0
  gather ("date", "germinated", 8: last_col() )%>% # back in long format frrom 8th colum to the last
  arrange (species, accession, code, incubator, petridish, date)%>% # sort row observations this way
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  group_by (species, accession, code, incubator, petridish, time)%>% 
  summarise(seeds_germ = sum(germinated)) %>%
  mutate(cumulative = cumsum(seeds_germ)) %>%
  filter (between(time, 323, 430)) %>% # amount of days calculated from the dates 323 = 17/6 and 430 = 19/09
  ungroup () %>%
  group_by (species, accession, code, incubator) %>%
  summarise (seeds_germ = sum(seeds_germ)) %>%
  merge(appendix)%>%
  group_by (community, species, accession, code, incubator)%>%
  mutate(summer_germ = (seeds_germ/viable)*100, 
         summer_germ = round (summer_germ, digit = 2))%>%
  select (community, species, abundance, family, accession, code, recolection, incubator, 
          total,viable, viablePER, total_germ,germPER,  autumn_germ, spring_germ,
          summer_germ) -> appendix
#winter
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))  %>%
  spread(date, germinated, fill = 0) %>% # wide format for dates, and fill Na with 0
  gather ("date", "germinated", 8: last_col()) %>% # back in long format frrom 8th colum to the last
  arrange (species, accession, code, incubator, petridish, date)%>% # sort row observations this way
  group_by (species, accession, code,incubator, petridish, date)%>% 
  summarise(seeds_germ = sum(germinated)) %>%
  mutate(cumulative = cumsum(seeds_germ)) %>%
  filter (date %in% c("2022-04-04", # # first check after winter before Tmin 2ºC
                      "2022-05-11","2022-05-12", "2022-05-18", "2022-05-19", 
                      "2022-05-25", "2022-05-26")) %>% #undersnow + first checks Tmin 2ºC
  ungroup () %>%
  group_by (species, accession, code, incubator) %>%
  summarise (seeds_germ = sum(seeds_germ)) %>%
  merge(appendix)%>%
  group_by (community, species, accession, code,)%>%
  mutate(winter_germ = (seeds_germ/viable)*100, 
         winter_germ = round (winter_germ, digit = 2))%>%
  select (community, species, abundance, family, accession, code, recolection, incubator, 
          total,viable, viablePER, total_germ,germPER,  autumn_germ, spring_germ,
          summer_germ, winter_germ) -> appendix
# t50
f50 <- function(df0) {
  lm(t50g ~ t50times, data = df0) -> mf1 # Linear model between time before and time after
  as.data.frame((0.5 - as.numeric(mf1$`coefficients`[1])) / as.numeric(mf1$`coefficients`[2])) # Use linear model to interpolate t50
}

read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y")) %>%
  group_by(species, accession, code, incubator, petridish) %>%
  mutate(days = difftime(date, min(date), units = "days")) %>% # calculate time from sowing to x date!
  arrange(species,accession, code,  incubator, petridish, days) %>%
  group_by (species, accession, code, incubator, petridish) %>% 
  mutate(cs = cumsum(germinated), # cumulative sum of germinated seeds 
         fg = max(cs) / viable, # max germinated (final germination) divided by viable seeds (proportion)
         g = cs / viable, # proportion of germination at each date
         Timeframe = ifelse(g >= 0.5, "Upper", "Lower")) %>% # # New column that divides dataset in germination < 0.5 or > 0.5
  group_by(species,accession, code, incubator, petridish, Timeframe) %>%
  mutate(t50times = ifelse(Timeframe == "Lower", max(days), min(days)), # keep only the number of days just before and after reaching t50
         t50g = ifelse(days == t50times, g, NA)) %>% # copy of germination proportion the days just before and after t50
  na.omit %>%
  select(species, accession, code, incubator, petridish, Timeframe, fg, t50times, t50g) %>%
  unique %>%
  group_by(species, accession, code, incubator, petridish) %>%
  filter(fg >= .50) %>% # filter to keep only the species that reach more than 0.5 in max final germination
  group_by(species, accession, code, incubator, petridish) %>%
  do(f50(.)) %>%
  rename(#FG = fg, 
    t50lm = `(0.5 - as.numeric(mf1$coefficients[1]))/as.numeric(mf1$coefficients[2])`) ->t50model

t50model %>%
 group_by(species, accession, code, incubator) %>%
  summarise(t50lm = mean (t50lm))%>%
  right_join(appendix, by =c("species", "accession", "code", "incubator")) %>%
  select(community, species, abundance, family, accession, code, recolection, incubator, 
        total,viable, viablePER, total_germ,germPER,  autumn_germ, spring_germ,
        summer_germ, winter_germ, t50lm) %>% 
  arrange (species, accession, code, incubator)-> appendix

#heat sum
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  group_by (species, accession, code, incubator, petridish) %>%
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
  group_by (species, accession, code, incubator, petridish) %>%
  do (heatsum(.)) -> HS

HS %>%
  group_by (species, accession, code, incubator) %>%
  summarise(HS = mean(HS)) %>%
  right_join(appendix, by =c("species", "accession", "code", "incubator")) %>%
  select(community, species, abundance, family, accession, code, recolection, incubator, 
         total,viable, viablePER, total_germ,germPER,  autumn_germ, spring_germ,
         summer_germ, winter_germ, t50lm, HS) %>% 
  arrange (species, accession, code, incubator)-> appendix
write.csv(appendix, "doc/table_sp.csv")

# mean germination rate (MGR) + synchrony (SYN) from GerminaR package
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
library (GerminaR)
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y")) %>%
  filter (species == "Thymus praecox") %>% ## 
  mutate(days = difftime(date, min(date), units = "days"), 
         days = round (days,digit =0)) %>% # 
  select (species, accession, code, incubator, petridish, total, viable, germinated, days) %>%
  group_by ( species, accession, code, incubator, petridish) %>%
  rename (seeds = viable) %>%
  mutate(days = factor(days)) %>%
  mutate(days = recode_factor(days, "3"="003", "6"="006", "35"="035", "48"="048", "63"="063", "68"="068","7"="007",
                              "36"="036", "49"="049", "17"="017", "16"="016", "50"="050", "46"="046", "14"="014",
                              "27"="027", "12"="012", "37"="037", "13"="013","38"="038", "39"="039", "41"="041", 
                              "56"="056", "59"="059", "61"="061", "19"="019", "4"="004", "20"="020", "34"="034", 
                              "70"="070", "83"="083", "98"="098", "45"="045", "33"="033", "30"="030", "42"="042", 
                              "28"="028", "15"="015", "51" ="051", "47"="047", "54"="054", "2"="002", "9"="009",
                              "43"="043", "29"="029", "53"="053", "40"="040", "25"="025", "11"="011", "65"="065", 
                              "60"="060", "31"="031", "64"="064", "55"="055", "1"="001", "5"="005", "21"="021", 
                              "18"="018", "62"="062", "44"="044", "10"="010"))%>% 
  mutate (days = paste("D", days, sep = "")) %>%
  spread(days, germinated)%>%
  data.frame () %>%
  replace(is.na(.),0) -> dat
ger_summary(SeedN = "seeds", evalName = "D", data=dat[,1: ncol(dat)])%>%
  select(species, accession, code, incubator, petridish, mgr, syn)%>%
  group_by (species, accession, code, incubator) %>%
  na.omit()%>%
  summarise (mgr = mean(mgr), syn = mean(syn)) %>%
  arrange (species, accession, code, incubator)


# final table x accession (no code = 1st + 2ndsow together)
read.csv("doc/table_sp_codes.csv", sep= ";") -> code_df
str(code_df)
code_df$Syn <- as.numeric(code_df$Syn)
str(code_df)
summary(code_df)
code_df %>%
  group_by(community, species, family, accession, recolection, incubator) %>%
  summarise(abundance = mean(abundance),
           total_seeds = sum(total_seeds), 
            viable_seeds = sum(viable_seeds), 
            viablePER = mean(viablePER), 
            seeds_germ = sum(seeds_germ), 
            germPER = mean(germPER), 
            autumn_germ = mean(autumn_germ), 
            spring_germ = mean (spring_germ), 
            summer_germ = mean(summer_germ),
            winter_germ = mean (winter_germ), 
            t50 = mean(t50), 
            HS= mean(HS), 
            MGT = mean (MGT), 
            Syn = mean (Syn)) %>%
  arrange(species)%>%
  write.csv("doc/table_population.csv")
