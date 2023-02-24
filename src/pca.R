#PCA
library(tidyverse);library(FactoMineR);library(factoextra)
library(vegan);library(ggrepel)

#traits correlation
read.csv("data/traits_sp_incubator.csv", sep = ";") %>%
  rename(total_germ = germPER)%>%
  select(total_germ, autumn_germ, spring_germ, summer_germ, winter_germ, t50, heat_sum, mass_50) %>%
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

########## mean value x species ####################
sp_traits <-read.csv ("data/traits_sp.csv", header = TRUE, sep = ";")
summary(sp_traits)
str(sp_traits)
# problem with NA, first subset data according to  germ traits and t 50 traits
#germination traits + t50 traits (remove species with NAs, from 55 to 42)
read.csv ("data/traits_sp.csv", header = TRUE, sep = ";") %>%
  rename(total_germ = germPER)%>%
  select(total_germ, autumn_germ, spring_germ, summer_germ, winter_germ, t50, heat_sum, mass_50) %>% # all traits
  #select(total_germ, autumn_germ, spring_germ, summer_germ, winter_germ, t50, heat_sum) %>% #germ + t50/heat
  #select(total_germ, autumn_germ, spring_germ, summer_germ, winter_germ) %>% #germ traits
  na.omit() -> traits

PCA_traits<- rda(traits, scale = TRUE) # scale = TRUE to standarize
head(summary (PCA_traits))

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

###################### Mean value x incubator #########################################  
sp_traits <-read.csv ("data/traits_sp_incubator.csv", header = TRUE, sep = ";")
summary(sp_traits)
str(sp_traits)
# problem with NA, first subset data according to  germ traits and t 50 traits
#germination traits + t50 traits (remove species with NAs, from 55 to 42)
read.csv ("data/traits_sp_incubator.csv", header = TRUE, sep = ";") %>%
  rename(total_germ = germPER)%>%
  #filter (incubator == "Fellfield") %>% # choose incubator
  filter (incubator == "Snowbed") %>% # choose incubator
  #select(total_germ, autumn_germ, spring_germ, summer_germ, winter_germ, t50, heat_sum, mass_50) %>% # all traits
  select(total_germ, autumn_germ, spring_germ, summer_germ, winter_germ, t50, heat_sum) %>% #germ + t50/heat
  #select(total_germ, autumn_germ, spring_germ, summer_germ, winter_germ) %>% #germ traits
  na.omit() -> traits

PCA_traits<- rda(traits, scale = TRUE) # scale = TRUE to standarize
head(summary (PCA_traits))

spe.scores <- scores(PCA_traits, display= "species", scaling = 0)
spe.scores
sort(abs(spe.scores[,1]), decreasing = TRUE) # sort variables with highest absolute correlation with 1st axis
sort(abs(spe.scores[,2]), decreasing = TRUE)# sort variables with highest absolute correlation with 2nd axis

smry <- summary(PCA_traits)
speciescores  <- data.frame(smry$sites[,1:2])       # species scores
traitscores  <- data.frame(smry$species[,1:2])     # traits scores

# all traits x species
read.csv ("data/traits_sp_incubator.csv", header = TRUE, sep = ";") %>%
  filter(incubator == "Snowbed")%>%
  na.omit() -> sp_traits
# germ + t50 + heat x species
read.csv ("data/traits_sp_incubator.csv", header = TRUE, sep = ";") %>%
  filter(incubator == "Snowbed")%>%
  select (!(mass_50))%>%
  na.omit() -> sp_traits
# germ traits x species
read.csv ("data/traits_sp_incubator.csv", header = TRUE, sep = ";") %>%
  filter(incubator == "Snowbed")%>%
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
  labs(title="Germination + t50 traits", subtitle = "Mean values per incubator: SNOWBED") +
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

###### seed mass #####
read.csv("data/seed_mass.csv", sep = ";") %>%
  group_by(macroclimate) %>%
  summarise (mass_50 = mean(mass_50))

ggplot()+
  geom_boxplot(data=seed_mass, aes(macroclimate, mass_50, fill=macroclimate), size =1) +
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
t.test(mass_50~macroclimate, data=seed_mass)         

# seed mass values are notably higher in EMditerraenan but not significantly different

#### species centroid ####
