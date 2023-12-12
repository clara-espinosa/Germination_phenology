library (tidyverse); library (dplyr)
library(MCMCglmm); library(DHARMa);
# Species levels: ecological preferences and germination responses
# Species preferences calculated from ibuttons data and inventories
# original scripts and data in Omaña80 and Picos Github repository (respectively)

# FROM ibuttons data: in picos from 02/10/2018 to 07/08/2019; in villabandin from 12/07/2021 to 29/05/2022
    # FDD: sum of degrees (daily mean) when daily Tmean is below 0 ºC (in absolute values)
    # GDD: sum of degrees (daily mean) when daily Tmean is above 5 ºC 
    # Snowdays = sum of days with snow cover, calculated as days when Tmax <0.5 and Tmin > -0.5 ºC
# These variables (altogether with bio1, bio2 and bio7) calculated mean x day, then month and finally year (whole datset)

# From inventory data: 80 plots per each community, at each plot species coverage estimated in %
# To calculate species preferences several filters applied:
    # consider only plots where each species has at least 10% coverage
    # climatic variables weighted per coverage 


####load data and select variables ####
read.csv("data/species_preferences/sp_pref_villa.csv", sep = ";") %>%
  select(species, community, bio1, bio2, bio7, FDD, GDD, Snw) -> sp_pref_villa

read.csv("data/species_preferences/sp_pref_picos.csv", sep = ";") %>%
  select(species, community, bio1, bio2, bio7, FDD, GDD, Snw) -> sp_pref_picos

read.csv("data/all_info.csv", sep =";") -> species
species %>% 
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) -> species
# merge datasets + fitler species used in move-along
# setdiff () interesting function to id missing value between 2 variables
# unique(df$species) # to check the use of the same names
# unique(species$species)
# unique(nnls_orig$tip.label)
# setdiff(nnls_orig$tip.label,df$animal)
# df %>% 
#  filter(! animal %in% nnls_orig$tip.label) %>% 
#  select(animal) %>% 
#  unique
#setdiff(sp_pref_villa$species, species$species)

#### PCA and correlations ####
sp_pref_villa%>%
  merge(species)%>%
  dplyr::select(community, species, accession, code,site, family, abundance, bio1:Snw) %>%
  group_by (species)%>%
  summarise(across (abundance:Snw, ~ mean(.x, na.rm = TRUE))) -> sp_pref_villa

sp_pref_villa[, 3:8] %>%
  FactoMineR::PCA() -> pca1

pca1$var$contrib

sp_pref_villa[, 3:8] %>%
  cor()

cbind((sp_pref_villa %>%  dplyr::select(species)), data.frame(pca1$ind$coord[, 1:2])) %>%
  mutate(species = factor(species)) -> pcaInds

pca1$var$coord[, 1:2] %>%
  data.frame %>%
  rownames_to_column(var = "Variable") %>%
  mutate(Variable = fct_recode(Variable, "Snow" = "Snw"))-> pcaVars

### Plot PCA
ggplot(pcaInds, aes(x = Dim.1, y = Dim.2)) +
  coord_fixed() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_segment(data = pcaVars, aes(x = 0, y = 0, xend = 3*Dim.1, yend = 3*Dim.2)) +
  geom_point(aes(fill = species, color = species), size = 4) +
  geom_label(data = pcaVars, aes(x = 3*Dim.1, y = 3*Dim.2, label = Variable),  show.legend = FALSE, size = 4) +
  ggthemes::theme_tufte() + 
  theme(text = element_text(family = "sans"),
        legend.position = "right", 
        legend.title = element_blank(),
        legend.text = element_text(size = 8, color = "black"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12, color = "black"),
        plot.margin = unit(c(0, 0.1, 0, 0), "cm")) +
  guides(fill = guide_legend(override.aes = list(shape = 22))) +
  scale_x_continuous(name = paste("Axis 1 (", round(pca1$eig[1, 2], 0),
                                  "% variance explained)", sep = "")) + 
  scale_y_continuous(name = paste("Axis 2 (", round(pca1$eig[2, 2], 0), 
                                  "% variance explained)", sep = "")) 

#setdiff(sp_pref_picos$species, species$species)
#setdiff(species$species, sp_pref_picos$species)
sp_pref_picos%>%
  merge(species)%>%
  dplyr::select(community, species, accession, code,site, family, abundance, bio1:Snw) %>%
  group_by (species)%>%
  summarise(across (abundance:Snw, ~ mean(.x, na.rm = TRUE))) -> sp_pref_picos

sp_pref_picos[, 3:8] %>%
  FactoMineR::PCA() -> pca1

pca1$var$contrib

sp_pref_picos[, 3:8] %>%
  cor()

cbind((sp_pref_picos %>%  dplyr::select(species)), data.frame(pca1$ind$coord[, 1:2])) %>%
  mutate(species = factor(species)) -> pcaInds

pca1$var$coord[, 1:2] %>%
  data.frame %>%
  rownames_to_column(var = "Variable") %>%
  mutate(Variable = fct_recode(Variable, "Snow" = "Snw"))-> pcaVars

### Plot PCA
ggplot(pcaInds, aes(x = Dim.1, y = Dim.2)) +
  coord_fixed() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_segment(data = pcaVars, aes(x = 0, y = 0, xend = 3*Dim.1, yend = 3*Dim.2)) +
  geom_point(aes(fill = species, color = species), size = 4) +
  geom_label(data = pcaVars, aes(x = 3*Dim.1, y = 3*Dim.2, label = Variable),  show.legend = FALSE, size = 4) +
  ggthemes::theme_tufte() + 
  theme(text = element_text(family = "sans"),
        legend.position = "right", 
        legend.title = element_blank(),
        legend.text = element_text(size = 8, color = "black"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12, color = "black"),
        plot.margin = unit(c(0, 0.1, 0, 0), "cm")) +
  guides(fill = guide_legend(override.aes = list(shape = 22))) +
  scale_x_continuous(name = paste("Axis 1 (", round(pca1$eig[1, 2], 0),
                                  "% variance explained)", sep = "")) + 
  scale_y_continuous(name = paste("Axis 2 (", round(pca1$eig[2, 2], 0), 
                                  "% variance explained)", sep = "")) 

#### GLMs/LM #### ####
# area between curves depending realized niche (FDD, GDD, SNW)
species %>%
  merge(sp_pref_villa) %>%
  group_by(species) %>%
  summarise(across(ABC_clean_data:Snw, ~ mean(.x, na.rm = TRUE))) %>%
  mutate(species= factor(species),
         Area_curves=ABC_clean_data) %>%
  filter(!Area_curves==0)%>%
  data.frame() %>%
  mutate(species = factor(species)) %>%
  mutate(ID = gsub(" ", "_", species), animal = ID) %>% 
  na.omit () -> sp_villa

species %>%
  merge(sp_pref_picos) %>%
  group_by(species) %>%
  summarise(across(ABC_clean_data:Snw, ~ mean(.x, na.rm = TRUE))) %>%
  mutate(species= factor(species),
         Area_curves=ABC_clean_data) %>%
  filter(!Area_curves==0)%>%
  data.frame() %>%
  mutate(species = factor(species)) %>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>% 
  mutate(ID = gsub(" ", "_", species), animal = ID) %>% 
  na.omit ()-> sp_picos

hist(sp_picos$Area_curves)
hist(sp_villa$Area_curves)

glmm(Area_curves ~ FDD + GDD + Snw, family = "gaussian", data= sp_picos) -> m1
summary(m1)

par(mfrow = c(2, 2))
lm(Area_curves ~ GDD , data= sp_villa) -> m1
lm(Area_curves ~ FDD , data= sp_villa) -> m1
lm(Area_curves ~ Snw , data= sp_villa) -> m1
lm(Area_curves ~ FDD + GDD + Snw, data=sp_villa) -> m1
summary(m1)
plot(m1)
summary(m1)$adj.r.squared

#Linear regression makes several assumptions about the data, such as :
# 1- Linearity of the data. The relationship between the predictor (x) and the outcome (y) is assumed to be linear.
# 2- Normality of residuals. The residual errors are assumed to be normally distributed.
# 3- Homogeneity of residuals variance. The residuals are assumed to have a constant variance (homoscedasticity)
# 4- Independence of residuals error terms.
# to check the assmptions we can create diagnostic plots, with the 2 following code lines
# par(mfrow = c(2, 2))
# plot(model)
# The diagnostic plots show residuals in four different ways:
# 1- Residuals vs Fitted. Used to check the linear relationship assumptions. 
#A horizontal line, without distinct patterns is an indication for a linear relationship, what is good.
# 2- Normal Q-Q. Used to examine whether the residuals are normally distributed. 
#It's good if residuals points follow the straight dashed line.
# 3- Scale-Location (or Spread-Location). Used to check the homogeneity of variance of the residuals (homoscedasticity). 
# Horizontal line with equally spread points is a good indication of homoscedasticity. 
# 4- Residuals vs Leverage. Used to identify influential cases, that is extreme values that might influence the 
# regression results when included or excluded from the analysis.

plot(x= bio1, y= GDD, data= sp)

#MCMC GLMM for AREA between curves ####
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
nite = 10000
nthi = 10
nbur = 100
### Gaussian priors
priors <- list(R = list(V = 1, nu = 0.2),
               G = list(G1 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3),
                        G2 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3))) 
#G3 = list(V = 1, nu = 0.2, alpha.mu = 0, alpha.V = 1e3)))


# Gaussian model
MCMCglmm::MCMCglmm( Area_curves ~ FDD + GDD + Snw,
                   random = ~ animal + ID,
                   family = "gaussian", pedigree = nnls_orig, prior = priors, data = sp_villa,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> g1
# save(m1, file = "results/mcmc.Rdata")
x11()
plot(g1)

# load("results/mcmc.Rdata")
summary(g1)

# exploratory plots ####
par(mfrow = c(1, 1))
library(ggrepel)
str(sp_picos)
ggplot (sp_villa, aes(x=Snw, y=Area_curves)) + #, color = species
  geom_point() +
  geom_smooth(method = "lm") +
  geom_text_repel (aes(label= species))+
  #geom_text(x = 25, y = 300, label = lm_eqn(df), parse = TRUE)+
  #theme_bw()+
  theme(legend.position = "none")

# GLMs traits depending realized niche (FDD, GDD, SNW) ####
read.csv("doc/table_population.csv", sep = ";") %>%
  merge(sp_pref_picos) %>%
  merge(species) %>%
  select(community, species, accession, incubator, final_germ, autumn_germ, spring_germ, summer_germ,
          winter_germ, t50, Env_heat_sum, MGT, Syn, FDD, GDD, Snw, Area_curves) -> traits_picos

glm(Syn ~ FDD + GDD + Snw + incubator, family = "gaussian", data= traits_picos) -> m1
summary(m1)

