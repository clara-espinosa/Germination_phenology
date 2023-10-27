#PCA
library(tidyverse);library(FactoMineR);library(factoextra)
library(vegan);library(ggrepel)

# seed mass data ######
read.csv("data/seed_mass.csv", sep = ";") %>%
  group_by(community, species)%>%
  summarise(seed_mass = mean(mass_50)) -> seed_mass
seed_mass$species

ggplot()+
  geom_boxplot(data=seed_mass, aes(community, seed_mass, fill=community), size =1) +
  #facet_grid (.~macroclimate) +
  #scale_fill_manual (name= "Incubator", values = c ("Fellfield"= "chocolate2", "Snowbed" ="deepskyblue3")) +
  scale_fill_manual (name = "Macroclimate",values = c ("Mediterranean"= "darkgoldenrod1" , "Temperate" = "forestgreen")) +
  labs(title = "Seed mass", x = "Macroclimate", y = "Mass (mg)") +
  theme_classic(base_size = 16)+
  theme(plot.title = element_text (hjust = 0.5, size = 30),
        strip.text = element_text(face = "italic", size = 24), 
        legend.position = "right",
        legend.title = element_text(size=20), 
        legend.text = element_text(size=16),
        legend.background = element_rect(fill="transparent",colour=NA),
        axis.title.y = element_text (size=18), 
        axis.title.x = element_text (size=18))

t.test(seed_mass~community, data=seed_mass)  

####traits correlation ####
read.csv("data/traits_inc_sp.csv", sep = ";") %>%
  merge(seed_mass, by= c("community", "species"))%>%
  select(total_germ, autumn_germ, spring_germ, summer_germ, winter_germ, t50, EHS, delay, area_curves, seed_mass) %>%
  na.omit()-> traits_cor
           
traits_cor %>% cor() -> biocor_all
biocor_all
# correlation significance test
cor.mtest <- function(biocor_all, ...) {
  biocor_all <- as.matrix(biocor_all)
  n <- ncol(biocor_all)
  p.biocor_all<- matrix(NA, n, n)
  diag(p.biocor_all) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(biocor_all[, i], biocor_all[, j], ...)
      p.biocor_all[i, j] <- p.biocor_all[j, i] <- tmp$p.value
    }
  }
  colnames(p.biocor_all) <- rownames(p.biocor_all) <- colnames(biocor_all)
  p.biocor_all
}
# matrix of the p-value of the correlation
p.biocor_all <- cor.mtest(traits_cor)
# no correlation significance

library (corrplot)
x11()
corrplot (biocor_all, method="number", type = "upper", order = "hclust", 
          p.mat = p.biocor_all, sig.level = 0.05)

########## mean value x species and incubator ####################
setdiff(sp_traits$species, seed_mass$species)

read.csv("data/traits_inc_sp.csv", sep = ";") -> sp_traits   #%>%
  #merge(seed_mass, by= c("community", "species")) 

summary(sp_traits)
str(sp_traits)

# problems with NA (mainly from t50/EHS/Delay), 
# SUBSET DATA x COMMUNITY and INCUBATOR
# first ALL TRAITS  (some NA will be removed)
# second germination traits, area curves and seed mass
# first subset data according to  germ traits and t 50 traits
#germination traits + t50 traits (remove species with NAs, from 55 to 42)
read.csv("data/traits_inc_sp.csv", sep = ";") %>%
  mutate(across(c(community, species, family, incubator), as.factor))%>%
  merge(seed_mass, by= c("community", "species"))%>%
  filter (community == "Mediterranean") %>%
  filter (incubator == "Snowbed") %>%
  select(community, species, family, incubator, total_germ, autumn_germ, spring_germ, 
         summer_germ, winter_germ, EHS, area_curves, seed_mass ) %>% #    all traits
  #select(total_germ, autumn_germ, spring_germ, summer_germ, winter_germ, area_curves, seed_mass) %>% 
  na.omit() -> traits

# USEFUL!!Shorten sp names with 4 letters from genus name and 4 letter from species name 
# Sp$Species <- make.cepnames(Sp$Species, seconditem = FALSE) works if Sp are in columns 
# Sp <- t(Sp) # transpose the matrix so species are in columns instead of in rows.

traits$species <-make.cepnames(traits$species)

### PCA

traits[, 5:12] %>%
  FactoMineR::PCA() -> pca1

cbind((traits %>%  dplyr::select(species, family)), data.frame(pca1$ind$coord[, 1:4]))-> pcaInds

pca1$var$coord[, 1:2] %>%
  data.frame %>%
  rownames_to_column(var = "Variable") -> pcaVars

### Plot PCA
#x11()
ggplot(pcaInds, aes(x = Dim.1, y = Dim.2)) +
  coord_fixed() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_segment(data = pcaVars, aes(x = 0, y = 0, xend = 3*Dim.1, yend = 3*Dim.2)) +
  geom_point(aes(fill = species, color= family), size = 4) +
  geom_label(data= pcaVars, aes(x = 3*Dim.1, y = 3*Dim.2, label = Variable),  show.legend = FALSE, size = 4.5) +
  geom_text_repel (data =pcaInds, aes(x=Dim.1, y = Dim.2, label = species, color = family ), show.legend = FALSE, size = 4.5)+
  labs (title = "Mediterranean community", subtitle = "Snowbed incubator") +
  ggthemes::theme_tufte() + 
  theme(title = element_text(family= "sans", face = "bold", size= 24),
        text = element_text(family = "sans"),
        legend.position = "none", 
        legend.title = element_blank(),
        legend.text = element_text(size = 18, color = "black"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 18, color = "black"),
        plot.margin = unit(c(0, 0.1, 0, 0), "cm")) +
  #guides(fill = guide_legend(override.aes = list(shape = 22))) +
  scale_x_continuous(name = paste("Axis 1 (", round(pca1$eig[1, 2], 0),
                                  "% variance explained)", sep = ""),limits = c(-4, 4) ) + 
  scale_y_continuous(name = paste("Axis 2 (", round(pca1$eig[2, 2], 0), 
                                  "% variance explained)", sep = ""), limits = c(-2, 3)) 
pca1$eig
pca1$var




PCA_traits<- rda(traits[,5:12], scale = TRUE) # scale = TRUE to standarize
head(summary (PCA_traits))

cbind((traits %>%  dplyr::select(species, family, incubator)), data.frame(PCA_traits$ind$coord[, 1:2])) %>%
spe.scores <- scores(PCA_traits, display= "species", scaling = 0)
spe.scores
sort(abs(spe.scores[,1]), decreasing = TRUE) # sort variables with highest absolute correlation with 1st axis
sort(abs(spe.scores[,2]), decreasing = TRUE)# sort variables with highest absolute correlation with 2nd axis
site.scores <- scores(PCA_traits, display= "sites", scaling = 0)
site.scores

biplot (PCA_traits, display = c("species", "sites"), scaling = "species")

smry <- summary(PCA_traits)
speciescores  <- data.frame(smry$sites[,1:2])       # species scores
traitscores  <- data.frame(smry$species[,1:2])     # traits scores

# all traits x species
read.csv ("data/traits_sp.csv", header = TRUE, sep = ";") %>%
  na.omit() -> sp_traits
# germ + t50 + heat x species
read.csv ("data/traits_sp.csv", header = TRUE, sep = ";") %>%
  select (!(mass_50))%>%
  na.omit() -> sp_traits
# germ traits x species
read.csv ("data/traits_sp.csv", header = TRUE, sep = ";") %>%
  select (!(mass_50))%>%
  select (!(t50))%>%
  select (!(heat_sum))%>%
  na.omit()-> sp_traits
# plot
x11()
ggplot() + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point (data= speciescores,aes(x=PC1, y=PC2, color = sp_traits$macroclimate), size = 0.1) + 
  geom_text(data = speciescores, 
                   aes(x=PC1,y=PC2,label=sp_traits$species, color = sp_traits$macroclimate), size=4, show.legend = FALSE) +
#scale_color_manual(name = "Habitat", values = c ("Specialist" = "blue", "Not_specialist" = "darkred")) +
  scale_color_manual (name = "Macroclimate",values = c ("Mediterranean"= "darkgoldenrod1" , "Temperate" = "forestgreen")) +
  labs(title="Germination + t50 traits", subtitle = "Mean values per species") +
  scale_x_continuous (limits = c(-2,2), breaks = seq (-1, 1, by= 0.5)) +
  theme_classic(base_size = 16)+
  theme(plot.title = element_text (hjust = 0.02, size = 30),
        plot.subtitle = element_text (hjust = 0.02, size = 22),
        strip.text = element_text(face = "italic", size = 24), 
        legend.position = "right",
        legend.title = element_text(size=20), 
        legend.text = element_text(size=16),
        legend.background = element_rect(fill="transparent",colour=NA),
        axis.title.y = element_text (size=18), 
        axis.title.x = element_text (size=18)) +
  guides(color = guide_legend(override.aes = list(linetype = 0, size=5))) + 
  geom_segment(data=traitscores, aes(x=0, xend=PC1, y=0, yend=PC2), size = 1, 
               color="red", arrow=arrow(length=unit(0.02,"npc"))) +
  geom_text(data=traitscores, 
            aes(x=PC1,y=PC2,label=rownames(traitscores),
                hjust=0.5*(1-sign(PC1)),vjust=0.5*(1-sign(PC2))), 
            color="black", size=5)
