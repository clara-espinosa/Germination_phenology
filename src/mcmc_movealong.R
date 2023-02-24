library(FD); library(vegan); library(FactoMineR); library(emmeans);
library(tidyverse); library(ggrepel); library(cowplot); library (ggplot2);
library (lubridate); library(binom); library (GerminaR); library(rstatix);
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
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  group_by(species, code, incubator, petridish) %>%
  summarise(total_germ = sum(germinated),
            viable= last(viable),
            total = last(total)) %>% # 
  mutate(germPER = (total_germ/viable) *100, # 
         germPER = round (germPER, digit =2)) %>%
  mutate(viablePER=(viable/total) *100, 
         viablePER = round(viablePER, digit =2)) -> viables 
### seeds germ + viables/community filtered ####
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  group_by(species, code, incubator, petridish) %>%
  summarise(total_germ = sum(germinated)) %>% # 
  merge(viables) %>%
  filter(viablePER>25)%>%
  filter(germPER>0)%>%
  merge(read.csv("data/all_info.csv", sep =";")) %>% 
  group_by (community, incubator) %>%
  summarise (total_germ = sum(total_germ),
             viable = sum (viable))-> viables_community

# ANALISIS from raw data instead of percentages, response variables should be a proportion
###################################### GERMINATION RESPONSE TO MICROHABIAT ###################################
#### GERMINATION RATE  ####
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  group_by (species, code, incubator, petridish) %>%
  mutate (cumulative = cumsum(germinated))%>%
  merge(viables)%>%
  arrange(species, code, accession, incubator, petridish, time) %>%
  filter(viablePER>25)%>%
  filter(germPER>0)%>%
  select(species, code, incubator, petridish, germinated, cumulative, viable, time)%>%
  merge(species, by = c("code", "species")) %>%
  convert_as_factor(code, species, incubator, family, community, habitat, germ_strategy) %>%
  mutate(species= str_replace(species, "Cerastium sp", "Cerastium pumilum"))%>% # Cerastium sp as Cerastium pumilum;
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>% # Minuartia CF as Minuartia arctica;
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>%
  mutate(ID = gsub(" ", "_", species), animal = ID) %>% 
  select(!family) %>%
  filter (community == "Temperate") %>%
  na.omit () -> df
summary(df)
#### TOTAL GERMINATION  ####
viables_germ %>% 
  filter(viablePER>25)%>%
  filter(germPER>0)%>%
  select (species, code, incubator, petridish, viable, total_germ) %>%
  mutate(species= str_replace(species, "Cerastium sp", "Cerastium pumilum"))%>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>% 
  #filter(str_length (code)<10)%>% ## filter accesions codes shorter than 10 character (all 2nd sow has an extra "b")
  merge(species, by = c("code", "species")) %>%
  mutate(code=factor(code)) %>%
  mutate(species=factor(species)) %>%
  mutate(incubator=factor(incubator)) %>%
  mutate(family=factor(family)) %>%
  mutate(community=factor(community)) %>%
  mutate(habitat=factor(habitat)) %>%
  mutate(germ_strategy=factor(germ_strategy)) %>%
  arrange(species, code, accession, incubator, petridish)  %>%
  mutate(ID = gsub(" ", "_", species), animal = ID) %>% 
  select(!family) %>%  
  na.omit ()-> df
summary(df)
df %>%
  select (code, species, incubator, petridish, total_germ, viable, community)%>%
  group_by (community, incubator) %>%
  summarise (germ = (total_germ/viable)*100) %>%
  get_summary_stats(type ="full")
#### AUTUMN (Mid November) ####
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  #filter (species == "Minuartia verna")%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  group_by (species, code, incubator, petridish, time) %>%
  summarise(seeds_germ = sum(germinated)) %>%
  mutate(cumulative = cumsum(seeds_germ)) %>%
  filter (between(time, 1, 105)) %>% ## 105 = 12/11 last check before winter
  group_by (species, code, incubator, petridish) %>%
  summarise (seeds_germ = sum(seeds_germ)) %>%
  merge(viables_germ, by = c("code", "species", "incubator", "petridish"))%>%
  filter(viablePER>25)%>%
  filter(germPER>0)%>%
  select(species, code, incubator, petridish, seeds_germ, viable) %>%
  mutate(species= str_replace(species, "Cerastium sp", "Cerastium pumilum"))%>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>% 
  merge(species, by = c("code", "species")) %>%
  mutate(code=factor(code)) %>%
  mutate(species=factor(species)) %>%
  mutate(incubator=factor(incubator)) %>%
  mutate(family=factor(family)) %>%
  mutate(community=factor(community)) %>%
  mutate(habitat=factor(habitat)) %>%
  mutate(germ_strategy=factor(germ_strategy)) %>%
  arrange(species, code, accession, incubator, petridish)  %>%
  mutate(ID = gsub(" ", "_", species), animal = ID) %>% 
  select(!family) %>%  
  na.omit ()-> df 
summary (df)
df %>%
  select (code, species, incubator, petridish, seeds_germ, viable, community)%>%
  group_by (community, incubator) %>%
  summarise (germ = (seeds_germ/viable)*100) %>%
  get_summary_stats(type ="full")
#### SPRING (Mid-June) ####
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  #filter(species == "Gypsophila repens") %>%
  spread(date, germinated, fill = 0) %>% # wide format for dates, and fill Na with 0
  gather ("date", "germinated", 8: last_col() )%>% # back in long format frrom 8th colum to the last
  arrange (species, accession, code, incubator, petridish, date)%>% # sort row observations this way
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  group_by (species, code, incubator, petridish, time)%>% 
  summarise(seeds_germ = sum(germinated)) %>%
  mutate(cumulative = cumsum(seeds_germ)) %>%
  filter (between(time, 106, 322)) %>% # amount of days calculated from the dates 106 = 13/11 and 322 = 16/6
  ungroup () %>%
  group_by (species, code, incubator, petridish) %>%
  summarise (seeds_germ = sum(seeds_germ),
             cumulative = last(cumulative)) %>%
  merge(viables_germ)%>%
  filter(viablePER>25)%>%
  filter(germPER>0)%>%
  select(species, code, incubator, petridish, seeds_germ, cumulative, viable) %>%
  mutate(species= str_replace(species, "Cerastium sp", "Cerastium pumilum"))%>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>% 
  merge(species, by = c("code", "species")) %>%
  mutate(code=factor(code)) %>%
  mutate(species=factor(species)) %>%
  mutate(incubator=factor(incubator)) %>%
  mutate(family=factor(family)) %>%
  mutate(community=factor(community)) %>%
  mutate(habitat=factor(habitat)) %>%
  mutate(germ_strategy=factor(germ_strategy)) %>%
  arrange(species, code, accession, incubator, petridish)  %>%
  mutate(ID = gsub(" ", "_", species), animal = ID) %>% 
  select(!family) %>%  
  na.omit ()-> df   
summary (df)
df %>%
  select (code, species, incubator, petridish, seeds_germ, viable, community)%>%
  group_by (community, incubator) %>%
  summarise (germ = (seeds_germ/viable)*100) %>%
  get_summary_stats(type ="full")
#### END-SUMMER (Mid- September) ####
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))  %>%
  #filter (species == "Armeria cantabrica") %>%
  spread(date, germinated, fill = 0) %>% # wide format for dates, and fill Na with 0
  gather ("date", "germinated", 8: last_col() )%>% # back in long format frrom 8th colum to the last
  arrange (species, accession, code, incubator, petridish, date)%>% # sort row observations this way
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  group_by (species, code, incubator, petridish, time)%>% 
  summarise(seeds_germ = sum(germinated)) %>%
  mutate(cumulative = cumsum(seeds_germ)) %>%
  filter (between(time, 323, 430)) %>% # amount of days calculated from the dates 323 = 17/6 and 430 = 19/09
  ungroup () %>%
  group_by (species, code, incubator, petridish) %>%
  summarise (seeds_germ = sum(seeds_germ),
             cumulative = last(cumulative)) %>%
  merge(viables_germ)%>%
  filter(viablePER>25)%>%
  filter(germPER>0)%>%
  select(species, code, incubator, petridish, seeds_germ, cumulative, viable) %>%
  mutate(species= str_replace(species, "Cerastium sp", "Cerastium pumilum"))%>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>% 
  merge(species, by = c("code", "species")) %>%
  mutate(code=factor(code)) %>%
  mutate(species=factor(species)) %>%
  mutate(incubator=factor(incubator)) %>%
  mutate(family=factor(family)) %>%
  mutate(community=factor(community)) %>%
  mutate(habitat=factor(habitat)) %>%
  mutate(germ_strategy=factor(germ_strategy)) %>%
  arrange(species, code, accession, incubator, petridish)  %>%
  mutate(ID = gsub(" ", "_", species), animal = ID) %>% 
  select(!family) %>%  
  na.omit ()-> df
df %>%
  select (code, species, incubator, petridish, seeds_germ, viable, community)%>%
  group_by (community, incubator) %>%
  summarise (germ = (seeds_germ/viable)*100) %>%
  get_summary_stats(type ="full")
## Winter germination  considering first checks after winter until Tmin>2ºC#####
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))  %>%
  #filter (species == "Conopodium majus") %>%
  spread(date, germinated, fill = 0) %>% # wide format for dates, and fill Na with 0
  gather ("date", "germinated", 8: last_col()) %>% # back in long format frrom 8th colum to the last
  arrange (species, accession, code, incubator, petridish, date)%>% # sort row observations this way
  #mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  group_by (species, code, incubator, petridish, date)%>% 
  summarise(seeds_germ = sum(germinated)) %>%
  mutate(cumulative = cumsum(seeds_germ)) %>%
  filter (date %in% c("2022-04-04", # # first check after winter before Tmin 2ºC
                      "2022-05-11","2022-05-12", "2022-05-18", "2022-05-19", 
                      "2022-05-25", "2022-05-26")) %>% #undersnow + first checks Tmin 2ºC
  ungroup () %>%
  group_by (species, code, incubator, petridish) %>%
  summarise (seeds_germ = sum(seeds_germ),
             cumulative = last(cumulative)) %>%
  merge (viables_germ) %>%
  filter(viablePER>25)%>%
  filter(germPER>0)%>%
  select(species, code, incubator, petridish, seeds_germ, viable)%>%
  arrange(species) %>%
  mutate(species= str_replace(species, "Cerastium sp", "Cerastium pumilum"))%>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>% 
  merge(species, by = c("code", "species")) %>%
  mutate(code=factor(code)) %>%
  mutate(species=factor(species)) %>%
  mutate(incubator=factor(incubator)) %>%
  mutate(family=factor(family)) %>%
  mutate(community=factor(community)) %>%
  mutate(habitat=factor(habitat)) %>%
  mutate(germ_strategy=factor(germ_strategy)) %>%
  arrange(species, code, accession, incubator, petridish)  %>%
  mutate(ID = gsub(" ", "_", species), animal = ID) %>% 
  select(!family) %>%  
  na.omit ()-> df
summary(df)
df %>%
  select (code, species, incubator, petridish, seeds_germ, viable, community)%>%
  group_by (community, incubator) %>%
  summarise (germ = (seeds_germ/viable)*100) %>%
  get_summary_stats(type ="full")
### modeling t50 from raw data ####

f50 <- function(df0) {
  lm(t50g ~ t50times, data = df0) -> mf1 # Linear model between time before and time after
  as.data.frame((0.5 - as.numeric(mf1$`coefficients`[1])) / as.numeric(mf1$`coefficients`[2])) # Use linear model to interpolate t50
}

read.csv("data/all_data.csv", sep = ";") %>%
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
  merge(viables_germ) %>%
  filter(viablePER>25)%>%
  filter(germPER>0)%>%
  mutate(species= str_replace(species, "Cerastium sp", "Cerastium pumilum"))%>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>% 
  merge(species, by = c("code", "species")) %>%
  mutate(code=factor(code)) %>%
  mutate(species=factor(species)) %>%
  mutate(incubator=factor(incubator)) %>%
  mutate(family=factor(family)) %>%
  mutate(community=factor(community)) %>%
  mutate(habitat=factor(habitat)) %>%
  mutate(germ_strategy=factor(germ_strategy)) %>%
  arrange(species, code, accession, incubator, petridish)  %>%
  mutate(ID = gsub(" ", "_", species), animal = ID) %>% 
  select(!family) %>%  
  na.omit () -> df
summary(df)
df %>%
  select (code, species, incubator, petridish, t50lm, community)%>%
  group_by (community, incubator) %>%
  summarise (t50 = t50lm) %>%
  get_summary_stats(type ="full")
### Heat_sum function, ##########    ####
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  group_by (species, code, incubator, petridish) %>%
  summarise (datesow = first (date)) %>% 
  merge (viables_germ)%>% 
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
  mutate(species= str_replace(species, "Cerastium sp", "Cerastium pumilum"))%>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>% 
  merge(species, by = c("code", "species")) %>%
  mutate(code=factor(code)) %>%
  mutate(species=factor(species)) %>%
  mutate(incubator=factor(incubator)) %>%
  mutate(family=factor(family)) %>%
  mutate(community=factor(community)) %>%
  mutate(habitat=factor(habitat)) %>%
  mutate(germ_strategy=factor(germ_strategy)) %>%
  arrange(species, code, accession, incubator, petridish)  %>%
  mutate(ID = gsub(" ", "_", species), animal = ID) %>% 
  select(!family) %>%  
  na.omit () -> df # punto significa que el objeto al que aplicar la funcion heat sum es el de la linea de arriba
summary(df)
df %>%
  select (code, species, incubator, petridish, HS, community)%>%
  group_by (community, incubator) %>%
  summarise (Env_heat_sum = HS) %>%
  get_summary_stats(type ="full")

#### PHYLO TREE AND MODEL SPECIFICATION FOR MULTINOMIAL####
 ### Read tree
phangorn::nnls.tree(cophenetic(ape::read.tree("results/tree.tree")), 
                      ape::read.tree("results/tree.tree"), rooted = TRUE) -> 
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
                           G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500),
                           G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500)))   
### TEST
MCMCglmm::MCMCglmm(cbind(cumulative, viable - cumulative) ~ incubator*time,
                   random = ~ animal + ID +  code:ID,
                   family = "multinomial2", pedigree = nnls_orig, prior = priors, data = df,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> m1


# save(m1, file = "results/mcmc.Rdata")
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

# Random effects species ID
summary(m1)$Gcovariances[2, 1] %>% round(2)
summary(m1)$Gcovariances[2, 2] %>% round(2) 
summary(m1)$Gcovariances[2, 3] %>% round(2) 

# Random effects code:ID
summary(m1)$Gcovariances[3, 1] %>% round(2)
summary(m1)$Gcovariances[3, 2] %>% round(2) 
summary(m1)$Gcovariances[3, 3] %>% round(2)

# Random effects code:ID
summary(m1)$Gcovariances[4, 1] %>% round(2)
summary(m1)$Gcovariances[4, 2] %>% round(2) 
summary(m1)$Gcovariances[4, 3] %>% round(2)

#### model specification for  gaussian (t50 + env heat)

#### PHYLO AND MODEL SPECIFICATION FOR GAUSSIAN ####
### Read tree
phangorn::nnls.tree(cophenetic(ape::read.tree("results/tree.tree")), 
                    ape::read.tree("results/tree.tree"), rooted = TRUE) -> 
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
                        G2 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        G3 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3)))
                        #G4 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3)))

# Gaussian model
MCMCglmm::MCMCglmm(delayS_F ~ habitat,
                   random = ~animal + ID + code:ID,
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

# Random effects species ID
summary(g1)$Gcovariances[2, 1] %>% round(2)
summary(g1)$Gcovariances[2, 2] %>% round(2) 
summary(g1)$Gcovariances[2, 3] %>% round(2) 

# Random effects code:ID
summary(g1)$Gcovariances[3, 1] %>% round(2)
summary(g1)$Gcovariances[3, 2] %>% round(2) 
summary(g1)$Gcovariances[3, 3] %>% round(2)

# Random effects code:ID
summary(g1)$Gcovariances[4, 1] %>% round(2)
summary(g1)$Gcovariances[4, 2] %>% round(2) 
summary(g1)$Gcovariances[4, 3] %>% round(2)

