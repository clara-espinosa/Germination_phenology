library(FD); library(vegan); library(FactoMineR); library(emmeans);
library(tidyverse); library(ggrepel); library(cowplot); library (ggplot2);
library (lubridate); library(binom); library (GerminaR); library(rstatix);
library(MCMCglmm); 
library(dplyr)
library(plyr)
theme_set(theme_cowplot(font_size = 10)) 

#setdiff () interesting function to id missing value between 2 variables
# unique(df$species) # to check the use of the same names
# unique(species$species)
# unique(nnls_orig$tip.label)
# setdiff(nnls_orig$tip.label,df$animal)
# df %>% 
#  filter(! animal %in% nnls_orig$tip.label) %>% 
    #  select(animal) %>% 
#  unique

####dataframe with temperature programs x incubator####
read.csv("data/date_temp.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y")) %>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  mutate(incubator = as.factor(incubator)) %>% 
  as.data.frame()-> temp
str(temp)
### Read species data ####
read.csv("data/all_info.csv", sep =";") -> species
# viables + final germ + filtering variables ####
# 1st sow only
detach(package:plyr)
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  group_by(species, code, incubator, petridish) %>%
  summarise(total_germ = sum(germinated),
            viable= last(viable),
            total = last(total)) %>% 
  mutate(germPER = (total_germ/viable) *100, # 
         germPER = round (germPER, digit =2)) %>%
  mutate(viablePER=(viable/total) *100, 
         viablePER = round(viablePER, digit =2))  -> viables 
### seeds germ + viables/community filtered ####
# 1st sow only
read.csv("data/clean data.csv", sep = ";") %>%
  #mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  group_by(species, code, incubator, petridish) %>%
  summarise(total_germ = sum(germinated)) %>% # 
  merge(viables) %>%
  filter(viablePER>25)%>%
  filter(germPER>0)%>%
  merge(species) %>% 
  group_by (community, incubator) %>%
  summarise (total_germ = sum(total_germ),
             viable = sum (viable))-> viables_community

# ANALISIS from raw data instead of percentages, response variables should be a proportion
###################################### GERMINATION RESPONSE TO MICROHABIAT ###################################
#### AUTUMN (Mid November) ####
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  group_by (species, code, incubator, petridish, time) %>%
  summarise(seeds_germ = sum(germinated)) %>%
  mutate(cumulative = cumsum(seeds_germ)) %>%
  filter (between(time, 1, 105)) %>% ## 105 = 12/11 last check before winter
  group_by (species, code, incubator, petridish) %>%
  summarise (seeds_germ = sum(seeds_germ)) %>%
  merge(viables, by = c("code", "species", "incubator", "petridish"))%>%
  filter(viablePER>25)%>%
  filter(germPER>0)%>%
  select(species, code, incubator, petridish, seeds_germ, viable) %>%
  merge(species, by = c("code", "species")) %>%
  convert_as_factor(code, species, incubator, family, community, habitat, germ_strategy) %>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>% 
  arrange(species, code, site, accession, incubator, petridish)  %>%
  mutate(ID = gsub(" ", "_", species), animal = ID) %>% 
  select(!family) %>%  
  #filter(community == "Temperate") %>%
  na.omit ()-> df 
summary (df)

#### WINTER germination  considering first checks after winter until Tmin>1ºC#####
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))  %>%
  spread(date, germinated, fill = 0) %>% # wide format for dates, and fill Na with 0
  gather ("date", "germinated", 8: last_col()) %>% # back in long format from 8th colum to the last
  arrange (species, code, incubator, petridish, date)%>% # sort row observations this way
  group_by (species, code, incubator, petridish, date)%>% 
  summarise(seeds_germ = sum(germinated)) %>%
  filter (date %in% c("2022-04-04", # # first check after winter before Tmin >2ºC
                      "2022-05-11","2022-05-12", "2022-05-18", "2022-05-19", 
                      "2022-05-25", "2022-05-26")) %>% #undersnow + first checks Tmin>2ºC
  ungroup () %>%
  group_by (species, code, incubator, petridish) %>%
  summarise (seeds_germ = sum(seeds_germ)) %>%
  merge(viables)%>%
  filter(viablePER>25)%>%
  filter(germPER>0)%>%
  select(species, code, incubator, petridish, seeds_germ, viable) %>%
  merge(species, by = c("code", "species")) %>%
  convert_as_factor(code, species, incubator, family, community, habitat, germ_strategy) %>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>% 
  arrange(species, code, site, accession, incubator, petridish)  %>%
  mutate(ID = gsub(" ", "_", species), animal = ID) %>% 
  select(!family) %>%  
  #filter(community == "Mediterranean") %>%
  na.omit ()-> df

summary(df)
df
#### SPRING (Tmin >1 to Mid-June = solstice) ####
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))  %>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  spread(time, germinated, fill = 0) %>% # wide format for dates, and fill Na with 0
  dplyr::select (!date) %>%
  gather ("time", "germinated", 7: last_col()) %>% # back in long format from 8th colum to the last
  group_by (species, code, incubator, petridish, time) %>% 
  filter (incubator == "Fellfield") %>% 
  filter (time %in% 254:331) %>% ## amount of days calculated from the dates 254 = 10/4 and 331 = 25/06
  group_by (species, code, incubator, petridish) %>%
  summarise(seeds_germ = sum(germinated)) %>%
  merge (viables)  -> Springgerm_F

read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))  %>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  spread(time, germinated, fill = 0) %>% # wide format for dates, and fill Na with 0
  dplyr::select (!date) %>%
  gather ("time", "germinated", 7: last_col()) %>% # back in long format from 8th colum to the last
  group_by (species, code, incubator, petridish, time) %>% 
  filter (incubator == "Snowbed") %>% 
  filter (time %in% 301:331) %>% ## amount of days calculated from the dates 301 = 27/5 and 331 = 25/06
  group_by (species, code, incubator, petridish) %>%
  summarise(seeds_germ = sum(germinated)) %>%
  merge (viables)  -> Springgerm_S

rbind(Springgerm_F, Springgerm_S) %>%
  filter(viablePER>25)%>%
  filter(germPER>0)%>%
  select(species, code, incubator, petridish, seeds_germ, viable) %>%
  merge(species, by = c("code", "species")) %>%
  convert_as_factor(code, species, incubator, family, community, habitat, germ_strategy) %>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>% 
  arrange(species, code,site,  accession, incubator, petridish)  %>%
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
  filter (between(time, 331, 430)) %>% # amount of days calculated from the dates 331 = 25/06 and 430 = 19/09
  ungroup () %>%
  group_by (species, code, incubator, petridish) %>%
  summarise (seeds_germ = sum(seeds_germ)) %>%
  merge(viables)%>%
  filter(viablePER>25)%>%
  filter(germPER>0)%>%
  select(species, code, incubator, petridish, seeds_germ, viable) %>%
  merge(species, by = c("code", "species")) %>%
  convert_as_factor(code, species, incubator, family, community, habitat, germ_strategy) %>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>% 
  arrange(species, code, site, accession, incubator, petridish)  %>%
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
  #filter(str_length (code)<10)%>% ## filter accesions codes shorter than 10 character (all 2nd sow has an extra "b")
  merge(species, by = c("code", "species")) %>%
  convert_as_factor(code, species, incubator, family, community, habitat, germ_strategy) %>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>% 
  arrange(species, code, site, accession, incubator, petridish)  %>%
  mutate(ID = gsub(" ", "_", species), animal = ID) %>% 
  select(!family) %>%  
  filter(community == "Temperate") %>%
  na.omit ()-> df
summary(df)

#### modeling t50 from raw data ####

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

t50model %>%
  merge(viables)%>%
  filter(viablePER>25)%>%
  filter(germPER>0)%>%
  select(species, code, incubator, petridish, t50lm) %>%
  merge(species, by = c("code", "species")) %>%
  convert_as_factor(code, species, incubator, family, community, habitat, germ_strategy) %>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>% 
  arrange(species, code, accession, incubator, petridish)  %>%
  mutate(ID = gsub(" ", "_", species), animal = ID) %>% 
  select(!family) %>%  
  #filter(community == "Mediterranean") %>%
  na.omit ()-> df
summary(df)

#### Heat_sum function ##########    ####
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
  convert_as_factor(code, species, incubator, family, community, habitat, germ_strategy) %>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>% 
  arrange(species, code, accession, incubator, petridish)  %>%
  mutate(ID = gsub(" ", "_", species), animal = ID) %>% 
  select(!family) %>%  
  filter (community == "Mediterranean") %>%
  na.omit () -> df # punto significa que el objeto al que aplicar la funcion heat sum es el de la linea de arriba
summary(df)


#### PHYLO TREE AND MODEL SPECIFICATION FOR MULTINOMIAL####
 ### Read tree
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
# change nu = 2 and alpha.V = 1000 for germination rate (originally nu = 1 and alpha.v = 500)
priors <- list(R = list(V = 1, nu = 50), 
                  G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500), 
                           G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500)))
                           #G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500)))   
### TEST
MCMCglmm::MCMCglmm(cbind(total_germ, viable - total_germ) ~ incubator, # *community
                   random = ~ animal + code:ID,
                   family = "multinomial2", pedigree = nnls_orig, prior = priors, data = df,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> m1

#save(m1, file = "results/mcmc.Rdata")
x11()
plot(m1)


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


#### PHYLO AND MODEL SPECIFICATION FOR GAUSSIAN (t50 and env heat)####
### Read tree
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
