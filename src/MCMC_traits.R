library(vegan); library(tidyverse); library(ggrepel); 
library (lubridate); library(binom); library(rstatix);
library(MCMCglmm); library (ggplot2);
library(dplyr); library(plyr)

####dataframe with temperature programs x incubator####
read.csv("data/date_temp.csv", sep = ";") %>% # read raw data for chambers temperature programs
  mutate(date = strptime(as.character(date), "%d/%m/%Y")) %>% # modify date format
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>% # calculate the time from the start of the experiment
  mutate(incubator = as.factor(incubator)) %>% 
  as.data.frame()-> temp
str(temp)

### Read species data ####
read.csv("data/all_info.csv", sep =";") -> species # species header information
# viables + final germ + filtering variables ####
detach(package:plyr)
read.csv("data/clean data.csv", sep = ";") %>% # read raw germination data
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>% # modify date format
  group_by(species, code, incubator, petridish) %>% # grouping variables
  summarise(total_germ = sum(germinated), 
            viable= last(viable),
            total = last(total)) %>% 
  mutate(germPER = (total_germ/viable) *100, # calculate germination proportion per petri dish
         germPER = round (germPER, digit =2)) %>%
  mutate(viablePER=(viable/total) *100, # calculate percentage of viable seeds per petri dish
         viablePER = round(viablePER, digit =2))  -> viables # viables per petri dish used to filter 
### seeds germ + viables/community filtered ####
read.csv("data/clean data.csv", sep = ";") %>%
  #mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  group_by(species, code, incubator, petridish) %>% # grouping variables
  summarise(total_germ = sum(germinated)) %>% # sum total germination
  merge(viables) %>% # merge viables object created above to filter raw data
  filter(viablePER>25)%>% # do no consider petridish with less than 25% of viable seeds at the ends of the experiment
  filter(germPER>0)%>% # do not consider petridish with 0 germination at the end of the experiment
  merge(species) %>% # join species header data
  group_by (community, incubator) %>% # new grouping variables
  summarise (total_germ = sum(total_germ),
             viable = sum (viable))-> viables_community # viables per community and incubator

# ANALISIS from raw data instead of percentages, response variables should be a proportion
###################################### GERMINATION RESPONSE TO MICROHABIAT ###################################
# start by calculating each dataframe pecific for each germination trait used in the study
#### AUTUMN (Mid November) ####
read.csv("data/clean data.csv", sep = ";") %>% # read raw germination data
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>% # change date format
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>% # calculate amount of days since start of the experiment
  group_by (species, code, incubator, petridish, time) %>% # grouping variables
  summarise(seeds_germ = sum(germinated)) %>% # sum germinated seeds by petri dish
  filter (between(time, 1, 105)) %>% ## filter time to consider only autumn period 105 days = 12/11 last check before winter
  group_by (species, code, incubator, petridish) %>% # new grouping variables
  summarise (seeds_germ = sum(seeds_germ)) %>% # 
  merge(viables, by = c("code", "species", "incubator", "petridish"))%>%
  filter(viablePER>25)%>% # do no consider petridish with less than 25% of viable seeds at the ends of the experiment 
  filter(germPER>0)%>% # do not consider petridish with 0 germination at the end of the experiment
  select(species, code, incubator, petridish, seeds_germ, viable) %>%
  merge(species, by = c("code", "species")) %>%
  convert_as_factor(code, species, incubator, family, community) %>% # explanatory variables to factors
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>% # change species name of a species we later confirmed the correct identification
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>% # change species name of a species we later confirmed the correct identification
  arrange(species, code, site, incubator, petridish)  %>%
  mutate(ID = gsub(" ", "_", species), animal = ID) %>% # necessary to include phylogeny as a random factor in our MCMC-GLMM models
  select(!family) %>%  
  #filter(community == "Temperate") %>% # activate or deactivate to condider community together or separately
  na.omit ()-> df # data frame correctly formated for MCMC-GLMM model analysis
summary (df)

#### WINTER germination  considering first checks after winter until Tmean>2ºC#####
# basic steps as above 
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))  %>%
  spread(date, germinated, fill = 0) %>% # wide format for dates, and fill Na with 0 for those species with already 100% germination
  gather ("date", "germinated", 8: last_col()) %>% # back in long format from 8th column to the last
  arrange (species, code, incubator, petridish, date)%>% # sort row observations this way
  group_by (species, code, incubator, petridish, date)%>% 
  summarise(seeds_germ = sum(germinated)) %>%
  filter (date %in% c("2022-04-04", # # first check after winter before T min >2ºC For FELLFIELD incubator
                      "2022-05-11","2022-05-12", "2022-05-18", "2022-05-19", 
                      "2022-05-25", "2022-05-26")) %>% #undersnow + first checks Tmin >2ºC For SNOWBED incubator
  ungroup () %>%
  group_by (species, code, incubator, petridish) %>%
  summarise (seeds_germ = sum(seeds_germ)) %>%
  merge(viables)%>%
  filter(viablePER>25)%>%
  filter(germPER>0)%>%
  select(species, code, incubator, petridish, seeds_germ, viable) %>%
  merge(species, by = c("code", "species")) %>%
  convert_as_factor(code, species, incubator, family, community) %>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>% 
  arrange(species, code, site, incubator, petridish)  %>%
  mutate(ID = gsub(" ", "_", species), animal = ID) %>% 
  select(!family) %>%  
  #filter(community == "Mediterranean") %>%
  na.omit ()-> df

summary(df)
df
#### SPRING (Tmean >2 to Mid-June = solstice) ####
detach(package:plyr)
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))  %>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  spread(time, germinated, fill = 0) %>% # wide format for dates, and fill Na with 0
  dplyr::select (!date) %>%
  gather ("time", "germinated", 7: last_col()) %>% # back in long format from 8th column to the last
  group_by (species, code, incubator, petridish, time) %>% 
  filter (incubator == "Fellfield") %>% 
  filter (time %in% 254:331) %>% ## amount of days corresponding to spring period calculated from the dates 254 = 10/4 and 331 = 25/06 in fellfield incubator
  group_by (species, code, incubator, petridish) %>%
  summarise(seeds_germ = sum(germinated)) %>%
  merge (viables)  -> Springgerm_F # spring germination only for fellfield incubator

read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))  %>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  spread(time, germinated, fill = 0) %>% # wide format for dates, and fill Na with 0
  dplyr::select (!date) %>%
  gather ("time", "germinated", 7: last_col()) %>% # back in long format from 8th colum to the last
  group_by (species, code, incubator, petridish, time) %>% 
  filter (incubator == "Snowbed") %>% 
  filter (time %in% 301:331) %>% ## amount of days corresponding to spring period from the dates 301 = 27/5 and 331 = 25/06 in snowbed incubator
  group_by (species, code, incubator, petridish) %>%
  summarise(seeds_germ = sum(germinated)) %>%
  merge (viables)  -> Springgerm_S # spring germination insnowbed incubator

rbind(Springgerm_F, Springgerm_S) %>% # merge both object with spring germination
  filter(viablePER>25)%>% # filtering same as above
  filter(germPER>0)%>%# filtering same as above
  select(species, code, incubator, petridish, seeds_germ, viable) %>%
  merge(species, by = c("code", "species")) %>%
  convert_as_factor(code, species, incubator, family, community) %>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>% # check explanation above
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>% # check explanation above
  arrange(species, code,site, incubator, petridish)  %>%
  mutate(ID = gsub(" ", "_", species), animal = ID) %>% 
  select(!family) %>%  
  #filter(community == "Mediterranean") %>%
  na.omit ()-> df
summary(df)

#### SUMMER (Mid- September) ####
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))  %>%
  spread(date, germinated, fill = 0) %>% # wide format for dates, and fill Na with 0
  gather ("date", "germinated", 8: last_col() )%>% # back in long format frrom 8th colum to the last
  arrange (species, code, incubator, petridish, date)%>% # sort row observations this way
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  group_by (species, code, incubator, petridish, time)%>% 
  summarise(seeds_germ = sum(germinated)) %>%
  filter (between(time, 331, 430)) %>% # amount of days considered wihtin summer period calculated from the dates 331 = 25/06 and 430 = 19/09
  ungroup () %>%
  group_by (species, code, incubator, petridish) %>%
  summarise (seeds_germ = sum(seeds_germ)) %>%
  merge(viables)%>%
  filter(viablePER>25)%>% # filtering as above
  filter(germPER>0)%>% # filtering as above 
  select(species, code, incubator, petridish, seeds_germ, viable) %>%
  merge(species, by = c("code", "species")) %>%
  convert_as_factor(code, species, incubator, family, community) %>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>% # see explanation above
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>%  # see explanation above
  arrange(species, code, site, incubator, petridish)  %>%
  mutate(ID = gsub(" ", "_", species), animal = ID) %>% 
  select(!family) %>%  
  filter(community == "Temperate") %>%
  na.omit ()-> df
summary(df)

#### TOTAL GERMINATION  ####
viables %>% 
  filter(viablePER>25)%>%
  filter(germPER>0)%>%
  select (species, code, incubator, petridish, viable, total_germ) %>%
  merge(species, by = c("code", "species")) %>%
  convert_as_factor(code, species, incubator, family, community) %>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>% 
  arrange(species, code, site, incubator, petridish)  %>%
  mutate(ID = gsub(" ", "_", species), animal = ID) %>% 
  select(!family) %>%  
  #filter(community == "Temperate") %>%
  na.omit ()-> df
summary(df)

#### modeling t50 from raw data ####

f50 <- function(df0) {
  lm(t50g ~ t50times, data = df0) -> mf1 # Linear model between time before and time after reaching 50% germination
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
  na.omit() %>%
  select(species, code, incubator, petridish, Timeframe, fg, t50times, t50g) %>%
  unique %>%
  group_by(species, code, incubator, petridish) %>%
  filter(fg >= .50) %>% # filter to keep only the species that reach more than 0.5 in max final germination
  group_by(species, code, incubator, petridish) %>%
  do(f50(.)) %>%
  rename(t50lm = `(0.5 - as.numeric(mf1$coefficients[1]))/as.numeric(mf1$coefficients[2])`) -> t50model # object with t50 dates from model 

t50model %>%
  merge(viables)%>%
  filter(viablePER>25)%>%
  filter(germPER>0)%>%
  select(species, code, incubator, petridish, t50lm) %>%
  merge(species, by = c("code", "species")) %>%
  convert_as_factor(code, species, incubator, family, community) %>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>% 
  arrange(species, code, incubator, petridish)  %>%
  mutate(ID = gsub(" ", "_", species), animal = ID) %>% 
  select(!family) %>%  
  #filter(community == "Mediterranean") %>%
  na.omit ()-> df
summary(df)

#### Heat_sum function ##########  
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  group_by (species, code, incubator, petridish) %>%
  summarise (datesow = first (date)) %>% 
  merge (viables)%>% 
  filter (germPER>49) %>%
  merge (t50model) %>%
  mutate (date50 = as.Date(datesow) + t50lm) %>% 
  filter(viablePER>25) -> t50_dates 

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
  do (heatsum(.)) %>%
  merge(species, by = c("code", "species")) %>%
  convert_as_factor(code, species, incubator, family, community) %>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>% 
  arrange(species, code, incubator, petridish)  %>%
  mutate(ID = gsub(" ", "_", species), animal = ID) %>% 
  select(!family) %>%  
  filter (community == "Mediterranean") %>%
  na.omit () -> df # punto significa que el objeto al que aplicar la funcion heat sum es el de la linea de arriba
summary(df)


#### PHYLO TREE AND MODEL SPECIFICATION FOR MULTINOMIAL####
 ### Read tree in resutls created following the phylo tree scripts in scr
phangorn::nnls.tree(cophenetic(ape::read.tree("results/tree.tree")), 
                      ape::read.tree("results/tree.tree"), method = "ultrametric") -> 
    nnls_orig
  
  nnls_orig$node.label <- NULL
  
### Set number of iterations
  nite = 1000000
  nthi = 100
  nbur = 100000
  
# shorter iterations
  # nite = 10000
   #nthi = 10
   #nbur = 100

### Set priors for germination models (as many prior as random factors)
priors <- list(R = list(V = 1, nu = 50), 
                  G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500), 
                           G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500)))  
### TEST
MCMCglmm::MCMCglmm(cbind(seeds_germ, viable - seeds_germ) ~ incubator, # community term will be added to answer the second question
                   random = ~ animal + code:ID, # animal reference to plant family and code:Id refence to collection site nested within species
                   family = "multinomial2", pedigree = nnls_orig, prior = priors, data = df,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> m1

#save(m1, file = "results/mcmc.Rdata")
x11()
plot(m1) # to check diagnostic plots


# load("results/mcmc.Rdata")
summary(m1)

### Random and phylo
# Calculate lambda http://www.mpcm-evolution.com/practice/online-practical-material-chapter-11/chapter-11-1-simple-model-mcmcglmm

  lambda <- m1$VCV[,"animal"]/(m1$VCV[,"animal"] + m1$VCV[,"units"]) 

mean(m1$VCV[,"animal"]/(m1$VCV[,"animal"] + m1$VCV[,"units"])) %>% round(2)
coda::HPDinterval(m1$VCV[,"animal"]/(m1$VCV[,"animal"] + m1$VCV[,"units"]))[, 1] %>% round(2)
coda::HPDinterval(m1$VCV[,"animal"]/(m1$VCV[,"animal"] + m1$VCV[,"units"]))[, 2] %>% round(2)

# Random effects animal
summary(m1)$Gcovariances[1, 1] %>% round(2) 
summary(m1)$Gcovariances[1, 2] %>% round(2) 
summary(m1)$Gcovariances[1, 3] %>% round(2)

# Random effects code:ID
summary(m1)$Gcovariances[2, 1] %>% round(2)
summary(m1)$Gcovariances[2, 2] %>% round(2) 
summary(m1)$Gcovariances[2, 3] %>% round(2) 


#### PHYLO AND MODEL SPECIFICATION FOR GAUSSIAN data (t50 and EHS)####
### Read tree from in resutls created following the phylo tree scripts in scr
phangorn::nnls.tree(cophenetic(ape::read.tree("results/tree.tree")), 
                    ape::read.tree("results/tree.tree"), method = "ultrametric") -> 
  nnls_orig

nnls_orig$node.label <- NULL

### Set number of iterations
nite = 1000000
nthi = 100
nbur = 100000

# shorter iterations
# nite = 10000
#nthi = 10
#nbur = 100
### Gaussian priors
priors <- list(R = list(V = 1, nu = 0.2),
               G = list(G1 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        G2 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3)))
                        #G3 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3)))


# Gaussian model
MCMCglmm::MCMCglmm(scale(HS) ~ incubator , #*community
                   random = ~animal + code:ID,
                   family = "gaussian", pedigree = nnls_orig, prior = priors, data = df,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> g1
# save(m1, file = "results/mcmc.Rdata")
x11()
plot(g1)

# load("results/mcmc.Rdata")
summary(g1)
### Random and phylo
# Calculate lambda http://www.mpcm-evolution.com/practice/online-practical-material-chapter-11/chapter-11-1-simple-model-mcmcglmm

lambda <- g1$VCV[,"animal"]/(g1$VCV[,"animal"] + g1$VCV[,"units"]) 

mean(g1$VCV[,"animal"]/(g1$VCV[,"animal"] + g1$VCV[,"units"])) %>% round(2)
coda::HPDinterval(g1$VCV[,"animal"]/(g1$VCV[,"animal"] + g1$VCV[,"units"]))[, 1] %>% round(2)
coda::HPDinterval(g1$VCV[,"animal"]/(g1$VCV[,"animal"] + g1$VCV[,"units"]))[, 2] %>% round(2)

# Random effects animal
summary(g1)$Gcovariances[1, 1] %>% round(2) 
summary(g1)$Gcovariances[1, 2] %>% round(2) 
summary(g1)$Gcovariances[1, 3] %>% round(2)

# Random effects code:ID
summary(g1)$Gcovariances[2, 1] %>% round(2)
summary(g1)$Gcovariances[2, 2] %>% round(2) 
summary(g1)$Gcovariances[2, 3] %>% round(2) 
