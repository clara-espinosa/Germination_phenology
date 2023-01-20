library(tidyverse)

#setdiff () interesting function to id missing value between 2 variables
### Read data
read.csv("data/all_info.csv", sep =";") -> species
str (species)
str(dataframe)

# modify column values of 3 sp names: 
# Cerastium sp as Cerastium pumilum; Minuartia CF as Minuartia arctica; Sedum sp (album cf) as Sedum album

species %>%
  mutate(species= str_replace(species, "Cerastium sp", "Cerastium pumilum"))%>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) -> species

# from raw data instead of percentatges, response variables should be a proportion
read.csv("data/all_data.csv", sep = ";") %>%
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
  na.omit ()%>% 
  filter(viable >1)%>% 
  filter (seeds_germ>1)-> df
summary(df)

str (df)
df <- as.data.frame(df)
df
# unique(df$species) # to check the use of the same names
# unique(species$species)
# unique(nnls_orig$tip.label)
# setdiff(nnls_orig$tip.label,df$animal)

df %>% 
  filter(! animal %in% nnls_orig$tip.label) %>% 
  select(animal) %>% 
  unique


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
  nite = 10000
  nthi = 10
  nbur = 100

### Set priors for germination models (as many prior as random factors)
   
priors <- list(R = list(V = 1, nu = 50), 
                  G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500), 
                           G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500),
                           G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500), 
                           G4 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500)))   

### TEST 1: compare germination before winter between regions
MCMCglmm::MCMCglmm(cbind(seeds_germ, viable - seeds_germ) ~ region,
                   random = ~ animal + ID + incubator + species,
                   family = "multinomial", pedigree = nnls_orig, prior = priors, data = df,
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

# Random effects bedrock
summary(m1)$Gcovariances[3, 1] %>% round(2)
summary(m1)$Gcovariances[3, 2] %>% round(2) 
summary(m1)$Gcovariances[3, 3] %>% round(2)

# Random effects site:bedrock
summary(m1)$Gcovariances[4, 1] %>% round(2)
summary(m1)$Gcovariances[4, 2] %>% round(2) 
summary(m1)$Gcovariances[4, 3] %>% round(2)


# data as traits in percentatges %
dataframe %>%
  mutate(species= str_replace(species, "Cerastium sp", "Cerastium pumilum"))%>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>%
  mutate (spring_germ= replace_na(spring_germ, 0)) %>%
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
  na.omit ()-> df ## Need to add filters for low viable and 0 germ


MCMCglmm::MCMCglmm(Mid_nov ~ region,
                   random = ~ animal + ID + incubator + species + code:species,
                   family = "multinomial", pedigree = nnls_orig, prior = priors, data = df,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> m1