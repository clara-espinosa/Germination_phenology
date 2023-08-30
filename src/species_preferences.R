library (tidyverse); library (dplyr)
library(MCMCglmm); 
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

##### Calculate areas between curves (Fellfield-Snowbed) for each species ####
# CALCULATED INDIVIDUALLY PER EACH SPECIES AND ADDED MANUALLY TO all_info.csv
#visualization
read.csv("data/all_data.csv", sep = ";") %>%
  rbind(data_vis)%>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  spread(date, germinated, fill = 0) %>% # wide format for dates, and fill Na with 0
  gather ("date", "germinated", 8: last_col() )%>% # back in long format frrom 8th colum to the last
  arrange (species, accession, code, incubator, petridish, date)%>% # sort row observations this way
  mutate(date = as.POSIXct(date))%>%
  merge (species) %>%
  dplyr::group_by(community, species, accession,  incubator, date) %>%
  dplyr::summarise(germinated = sum(germinated)) %>%
  dplyr::mutate(germinated = cumsum(germinated)) %>%
  merge(viables_sp) %>%
  mutate(germination = germinated/viable)  %>%
  mutate(ID = paste(community, accession)) %>%
  filter(species == "Iberis carnosa") %>% ### change species name
  ggplot(aes(date, germination, color = incubator, fill = incubator)) +
  geom_line(size = 2) +
  scale_color_manual (name= "Incubator", values = c ("Fellfield"= "chocolate2", "Snowbed" ="deepskyblue3")) +
  facet_wrap(~ accession, scales = "free_x", ncol = 2) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title= "Iberis carnosa", x = "Time ", y = "Germination proportion") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text (size = 30),
        strip.text = element_text (size = 24, face = "italic"),
        axis.title.y = element_text (size=24), 
        axis.title.x = element_text (size=24), 
        axis.text.x= element_text (size=18, angle = 75, vjust = 0.5),
        legend.position = "right",
        plot.margin = margin(t=0.5, l =0.5, b = 0.5, r =0.5, unit ="cm")) #-> C

#example with 1 species (x = time (not dates))
read.csv("data/all_data.csv", sep = ";") %>%
  #rbind(data_vis)%>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  spread(time, germinated, fill = 0) %>% # wide format for dates, and fill Na with 0
  mutate(date = as.POSIXct(date))%>%
  gather ("time", "germinated", 9: last_col()) %>% # back in long format frrom 9th colum to the last
  mutate (time = as.numeric(time)) %>%
  arrange (species, accession, code, incubator, petridish, time)%>% # sort row observations this way
  merge (species) %>%
  dplyr::group_by(community, species, accession,  incubator, time) %>%
  dplyr::summarise(germinated = sum(germinated)) %>%
  dplyr::mutate(germinated = cumsum(germinated)) %>%
  merge(viables_sp) %>%
  mutate(germination = germinated/viable)  %>%
  filter (community == "Temperate") %>%
  filter(species == "Silene ciliata")%>%
  select( species, accession, incubator, community, time, germination)%>%
  dplyr::group_by(species, incubator, time) %>%   # command line for 2 populations sp
  dplyr::summarise(germination = mean(germination)) %>%  # command line for 2 populations sp
  spread (incubator, germination) %>%
  data.frame()-> df

#calculate area by trapezoidal rule for 2 pop species 
# https://stackoverflow.com/questions/24742677/how-to-measure-area-between-2-distribution-curves-in-r-ggplot2
# 1st define each curve coordenates (within each species)
d0 <- df[,c(2,4)] # snowbed data
d1 <- df[,2:3] # snowbed fellfield
f0 <- approxfun(x= d0$time, y= d0$Snowbed)
f1 <- approxfun(x= d1$time, y= d1$Fellfield) 
# define x range of the density overlap
ovrng <- c(min (d0$time), max(d0$time))
# dividing it to sections (for example n=500)
i <- seq(min(ovrng), max(ovrng), length.out=500)
# calculating the distance between the density curves
h1 <- f0(i)-f1(i) # snowbed
h2 <- f1(i)-f0(i) # fellfield
#and using the formula for the area of a trapezoid we add up the areas
area1<-sum( (h1[-1]+h1[-length(h1)]) /2 *diff(i) *(h1[-1]>=0+0))     # for the regions where d1>d0
area2<-sum( (h2[-1]+h2[-length(h2)]) /2 *diff(i) *(h2[-1]>=0+0))     # for the regions where d1<d0
area_total <- area2 - area1 
area_total

#calculate area by trapezoidal rule for 1 pop species
# 1st define each curve coordenates (within each species)
d0 <- df[,c(4,6)] # snowbed data
d1 <- df[,4:5] # snowbed fellfield
f0 <- approxfun(x= d0$time, y= d0$Snowbed)
f1 <- approxfun(x= d1$time, y= d1$Fellfield) 
# define x range of the density overlap
 ovrng <- c(min (d0$time), max(d0$time))
# dividing it to sections (for example n=500)
i <- seq(min(ovrng), max(ovrng), length.out=500)
# calculating the distance between the density curves
h1 <- f0(i)-f1(i) # snowbed
h2 <- f1(i)-f0(i) # fellfield
#and using the formula for the area of a trapezoid we add up the areas
area1<-sum( (h1[-1]+h1[-length(h1)]) /2 *diff(i) *(h1[-1]>=0+0))     # for the regions where d1>d0
area2<-sum( (h2[-1]+h2[-length(h2)]) /2 *diff(i) *(h2[-1]>=0+0))     # for the regions where d1<d0
area_total <- area2 - area1 
area_total

#### GLMs ####
# area between curves depending realized niche (FDD, GDD, SNW)
species %>%
  merge(sp_pref_villa) %>%
  group_by(species) %>%
  summarise(across(Area_curves:Snw, ~ mean(.x, na.rm = TRUE))) %>%
  mutate(species= factor(species)) %>%
  data.frame() %>%
  mutate(species = factor(species)) %>%
  mutate(ID = gsub(" ", "_", species), animal = ID) %>% 
  na.omit () -> sp_villa

species %>%
  merge(sp_pref_picos) %>%
  group_by(species) %>%
  summarise(across(Area_curves:Snw, ~ mean(.x, na.rm = TRUE))) %>%
  mutate(species= factor(species)) %>%
  data.frame() %>%
  mutate(species = factor(species)) %>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>% 
  mutate(ID = gsub(" ", "_", species), animal = ID) %>% 
  na.omit ()-> sp_picos

rbind (sp_picos, sp_villa) -> sp
hist(sp$Area_curves)

glm(Area_curves ~ FDD + GDD + Snw, family = "gaussian", data= sp) -> m1
summary(m1)

lm(Area_curves ~ FDD + GDD , data= sp_picos) -> m1
summary(m1)

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
library(ggrepel)
str(sp_villa)
ggplot (sp_villa, aes(x=GDD, y= bio1)) + #, color = species
  geom_point() +
  geom_smooth(method = "lm") +
  geom_text_repel (aes(label= species))+
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

