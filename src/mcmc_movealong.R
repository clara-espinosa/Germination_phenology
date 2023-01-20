library(tidyverse)

#setdiff () interesting function to id missing value between 2 variables
# unique(df$species) # to check the use of the same names
# unique(species$species)
# unique(nnls_orig$tip.label)
# setdiff(nnls_orig$tip.label,df$animal)
# df %>% 
#  filter(! animal %in% nnls_orig$tip.label) %>% 
    #  select(animal) %>% 
#  unique
### Read data ####
read.csv("data/all_info.csv", sep =";") -> species
str (species)
# modify column values of 3 sp names: 
# Cerastium sp as Cerastium pumilum; Minuartia CF as Minuartia arctica; Sedum sp (album cf) as Sedum album
species %>%
  mutate(species= str_replace(species, "Cerastium sp", "Cerastium pumilum"))%>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) -> species

viables %>% ### necessary for filtering!!
  merge(finalgerm) %>%
  select (species, code, incubator, petridish, viable, viablePER, germPER) -> viables_germ

# seed germinated in spring + summer (needed for early season germ calculation)
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
  merge(viables_germ) %>%
  arrange(species) -> germ_AF

# from raw data instead of percentages, response variables should be a proportion
#### Mid-November  all info ### change input files all_data.csv or clean data.csv #####
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  group_by (species, code, incubator, petridish, time)%>%
  filter (between(time, 1, 105)) %>%
  group_by (species, code, incubator, petridish) %>%
  summarise(seeds_germ = sum(germinated))%>%
  merge(viables) %>%
  select(species, code, incubator, petridish, seeds_germ, viable) %>%
  mutate(species= str_replace(species, "Cerastium sp", "Cerastium pumilum"))%>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>% 
  merge(species, by = c("code", "species")) %>%
  mutate(code=factor(code)) %>%
  mutate(species=factor(species)) %>%
  mutate(incubator=factor(incubator)) %>%
  mutate(family=factor(family)) %>%
  mutate(region=factor(region)) %>%
  mutate(habitat=factor(habitat)) %>%
  mutate(germ_strategy=factor(germ_strategy)) %>%
  arrange(species, code, accession, incubator, petridish)  %>%
  mutate(ID = gsub(" ", "_", species), animal = ID) %>% 
  select(!family) %>%  
  na.omit ()-> df

#### Mid November filter viablePER>25 and total germination >0 ####
##change input files all_data.csv or clean data.csv 
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  group_by (species, code, incubator, petridish, time)%>%
  filter (between(time, 1, 105)) %>%
  group_by (species, code, incubator, petridish) %>%
  summarise(seeds_germ = sum(germinated))%>%
  merge(viables_germ)%>%
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
  mutate(region=factor(region)) %>%
  mutate(habitat=factor(habitat)) %>%
  mutate(germ_strategy=factor(germ_strategy)) %>%
  arrange(species, code, accession, incubator, petridish)  %>%
  mutate(ID = gsub(" ", "_", species), animal = ID) %>% 
  select(!family) %>%  
  na.omit ()-> df 

#### EARLY SEASON GERMINATION (BETWEEN when Tmin>2ºc AND Tmax < 10 ) ALL INFO####
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>% 
  group_by (species, code, incubator, petridish, time) %>% 
  filter (incubator == "Fellfield") %>% 
  filter (between(time, 257, 285)) %>% ## amount of days calculated from the dates 257 = 13/4 and 285 = 11/05
  group_by (species, code, incubator, petridish) %>%
  summarise(seeds_germ = sum(germinated)) %>%
  merge (germ_AF) %>%
  select(species, code, incubator, petridish, seeds_germ, germ_AW, viablePER, germPER)  -> Earlyseason_F

read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))  %>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  group_by (species, code, incubator, petridish, time) %>% 
  filter (incubator == "Snowbed") %>% 
  filter (between(time, 301, 322)) %>% ## amount of days calculated from the dates 301 = 27/5 and 322 = 17/06
  group_by (species, code, incubator, petridish) %>%
  summarise(seeds_germ = sum(germinated)) %>%
  merge (germ_AF) %>%
  select(species, code, incubator, petridish, seeds_germ, germ_AW, viablePER, germPER) -> Earlyseason_S

rbind(Earlyseason_F, Earlyseason_S) %>%
  arrange(species) %>%
  mutate(species= str_replace(species, "Cerastium sp", "Cerastium pumilum"))%>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>% 
  merge(species, by = c("code", "species")) %>%
  mutate(code=factor(code)) %>%
  mutate(species=factor(species)) %>%
  mutate(incubator=factor(incubator)) %>%
  mutate(family=factor(family)) %>%
  mutate(region=factor(region)) %>%
  mutate(habitat=factor(habitat)) %>%
  mutate(germ_strategy=factor(germ_strategy)) %>%
  arrange(species, code, accession, incubator, petridish)  %>%
  mutate(ID = gsub(" ", "_", species), animal = ID) %>% 
  select(!family) %>%  
  na.omit ()-> df

#### EARLY SEASON GERMINATION (BETWEEN when Tmin>2ºc AND Tmax < 10 ) filtered####
rbind(Earlyseason_F, Earlyseason_S) %>%
  arrange(species) %>%
  filter(viablePER>25)%>%
  filter(germPER>0)%>%
  select(species, code, incubator, petridish, seeds_germ, germ_AW) %>%
  mutate(species= str_replace(species, "Cerastium sp", "Cerastium pumilum"))%>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>% 
  merge(species, by = c("code", "species")) %>%
  mutate(code=factor(code)) %>%
  mutate(species=factor(species)) %>%
  mutate(incubator=factor(incubator)) %>%
  mutate(family=factor(family)) %>%
  mutate(region=factor(region)) %>%
  mutate(habitat=factor(habitat)) %>%
  mutate(germ_strategy=factor(germ_strategy)) %>%
  arrange(species, code, accession, incubator, petridish)  %>%
  mutate(ID = gsub(" ", "_", species), animal = ID) %>% 
  select(!family) %>%  
  na.omit ()-> df

#### TOTAL GERMINATION (all data)####
viables %>% 
  merge(finalgerm) %>% 
  select (species, code, incubator, petridish, viable, seeds_germ) %>%
  mutate(species= str_replace(species, "Cerastium sp", "Cerastium pumilum"))%>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>% 
  filter(str_length (code)<10)%>% ## filter accesions codes shorter  than 10 character (all 2nd sow has an extra "b")
  merge(species, by = c("code", "species")) %>%
  mutate(code=factor(code)) %>%
  mutate(species=factor(species)) %>%
  mutate(incubator=factor(incubator)) %>%
  mutate(family=factor(family)) %>%
  mutate(region=factor(region)) %>%
  mutate(habitat=factor(habitat)) %>%
  mutate(germ_strategy=factor(germ_strategy)) %>%
  arrange(species, code, accession, incubator, petridish)  %>%
  mutate(ID = gsub(" ", "_", species), animal = ID) %>% 
  select(!family) %>%  
  na.omit ()-> df
### TOTAL GERMINATION (filtered) ####
viables %>% 
  merge(finalgerm) %>% 
  filter(viablePER>25)%>%
  filter(germPER>0)%>%
  select (species, code, incubator, petridish, viable, seeds_germ) %>%
  mutate(species= str_replace(species, "Cerastium sp", "Cerastium pumilum"))%>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>% 
  filter(str_length (code)<10)%>% ## filter accesions codes shorter  than 10 character (all 2nd sow has an extra "b")
  merge(species, by = c("code", "species")) %>%
  mutate(code=factor(code)) %>%
  mutate(species=factor(species)) %>%
  mutate(incubator=factor(incubator)) %>%
  mutate(family=factor(family)) %>%
  mutate(region=factor(region)) %>%
  mutate(habitat=factor(habitat)) %>%
  mutate(germ_strategy=factor(germ_strategy)) %>%
  arrange(species, code, accession, incubator, petridish)  %>%
  mutate(ID = gsub(" ", "_", species), animal = ID) %>% 
  select(!family) %>%  
  na.omit ()-> df

#### PHYLO TREE AND MODEL SPECIFICATION ####
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
  #  nite = 10000
  #  nthi = 10
  #  nbur = 100

### Set priors for germination models (as many prior as random factors)
   
priors <- list(R = list(V = 1, nu = 50), 
                  G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500), 
                           G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500),
                           G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500), 
                           G4 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500)))   

### TEST: 
MCMCglmm::MCMCglmm(cbind(seeds_germ, viable - seeds_germ) ~ incubator,
                   random = ~ animal + ID + region + code:ID,
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

# Random effects incubator/region
summary(m1)$Gcovariances[3, 1] %>% round(2)
summary(m1)$Gcovariances[3, 2] %>% round(2) 
summary(m1)$Gcovariances[3, 3] %>% round(2)

# Random effects code:ID
summary(m1)$Gcovariances[4, 1] %>% round(2)
summary(m1)$Gcovariances[4, 2] %>% round(2) 
summary(m1)$Gcovariances[4, 3] %>% round(2)


