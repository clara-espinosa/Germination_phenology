library(tidyverse);library(V.PhyloMaker)
library(phylosignal);library(phylobase);library(ape);library(tidytree)

### Phylo tree both communitites #####
# always check family names with http://www.mobot.org/MOBOT/research/APweb/
read.csv("data/species.csv", sep =";") %>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>%
  select (species, family) %>%
  unique %>%
  separate(species, into = c("genus", "species"), sep = " ") %>%
  mutate(species = paste(genus, species),
         genus = genus,
         family = family) %>%
  #mutate(family = fct_recode(family, 
   #                           "Fabaceae" = "Leguminosae",
    #                          "Apiaceae" = "Umbelliferae",
     #                         "Asteraceae"= "Compositae",
      #                        "Orobanchaceae"= "Scrophulariaceae",
       #                       "Brassicaceae"= "Cruciferae",
        #                      "Asparagaceae" = "Liliaceae", 
         #                     "Lamiaceae" = "Labiatae")) %>%
  arrange(species) %>%
  na.omit %>%
  select(species, genus, family)-> 
  ranks1
#devtools::install_github("jinyizju/V.PhyloMaker")
phylo.maker(sp.list = ranks1,
            tree = GBOTB.extended, 
            nodes = nodes.info.1, 
            scenarios = "S3") ->
  tree

rm(ranks1)
x11()
plot(tree$scenario.3)

write.tree(tree$scenario.3, file = "results/tree.tree")

### Phylo tree temperate community #####
# always check family names with http://www.mobot.org/MOBOT/research/APweb/
read.csv("data/species.csv", sep =";") %>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>%
  filter(community == "Temperate") %>%
  select (species, family) %>%
  unique %>%
  separate(species, into = c("genus", "species"), sep = " ") %>%
  mutate(species = paste(genus, species),
         genus = genus,
         family = family) %>%
  #mutate(family = fct_recode(family, 
  #                           "Fabaceae" = "Leguminosae",
  #                          "Apiaceae" = "Umbelliferae",
  #                         "Asteraceae"= "Compositae",
  #                        "Orobanchaceae"= "Scrophulariaceae",
  #                       "Brassicaceae"= "Cruciferae",
  #                      "Asparagaceae" = "Liliaceae", 
  #                     "Lamiaceae" = "Labiatae")) %>%
  arrange(species) %>%
  na.omit %>%
  select(species, genus, family)-> 
  ranks1
#devtools::install_github("jinyizju/V.PhyloMaker")
phylo.maker(sp.list = ranks1,
            tree = GBOTB.extended, 
            nodes = nodes.info.1, 
            scenarios = "S3") ->
  tree_tem

rm(ranks1)
x11()
plot(tree_tem$scenario.3)

write.tree(tree_tem$scenario.3, file = "results/temperate_tree.tree")

### Phylo tree Mediterranean community #####
# always check family names with http://www.mobot.org/MOBOT/research/APweb/
read.csv("data/species.csv", sep =";") %>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>%
  filter(community == "Mediterranean") %>%
  select (species, family) %>%
  unique %>%
  separate(species, into = c("genus", "species"), sep = " ") %>%
  mutate(species = paste(genus, species),
         genus = genus,
         family = family) %>%
  #mutate(family = fct_recode(family, 
  #                           "Fabaceae" = "Leguminosae",
  #                          "Apiaceae" = "Umbelliferae",
  #                         "Asteraceae"= "Compositae",
  #                        "Orobanchaceae"= "Scrophulariaceae",
  #                       "Brassicaceae"= "Cruciferae",
  #                      "Asparagaceae" = "Liliaceae", 
  #                     "Lamiaceae" = "Labiatae")) %>%
  arrange(species) %>%
  na.omit %>%
  select(species, genus, family)-> 
  ranks1
#devtools::install_github("jinyizju/V.PhyloMaker")
phylo.maker(sp.list = ranks1,
            tree = GBOTB.extended, 
            nodes = nodes.info.1, 
            scenarios = "S3") ->
  tree_med

rm(ranks1)
x11()
plot(tree_med$scenario.3)

write.tree(tree_med$scenario.3, file = "results/mediterranean_tree.tree")

#### traits Area between curves /delay to 50% plotted according to phylogeny#####
# temperate #
ape::read.tree("results/temperate_tree.tree")-> tree
x11()
plot(tree)
appendix %>%
  select(community, species, area_curves, delayS_F)%>%
  group_by(community, species) %>%
  summarise(area_curves= mean(area_curves), 
            delay = mean(delayS_F))%>%
  filter (community == "Temperate") %>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>%
  select(species,area_curves, delay )%>%
  mutate(label= gsub(" ","_", species))%>%
  remove_rownames %>% 
  column_to_rownames(var="label")  %>% 
  select(area_curves)-> tem_phylo #,delay


obj_tem <- phylo4d(as(tree, "phylo4"), data.frame(tem_phylo), match.data=TRUE)

mat.col <- ifelse(tdata(obj_tem , "tip") < 0, "red", "grey35")
barplot.phylo4d(obj_tem, tree.ratio = 0.5,bar.col = mat.col, center=TRUE) #rainbow(37)
barplot(obj_tem,tree.type = "fan", tip.cex = 0.6, tree.open.angle = 160, trait.cex = 0.6)
dotplot.phylo4d (obj_tem, tree.type= "cladogram")
gridplot.phylo4d (obj_tem, tree.type = "fan", tip.cex = 0.6, show.trait = FALSE)
phyloSignal(obj_tem)
phyloCorrelogram (obj_tem)
# mediterranean#
ape::read.tree("results/mediterranean_tree.tree")-> tree
plot(tree)
appendix %>%
  select(community, species, area_curves, delayS_F)%>%
  group_by(community, species) %>%
  summarise(area_curves= mean(area_curves), 
            delay = mean(delayS_F))%>%
  filter (community == "Mediterranean") %>%
  select(species,area_curves, delay )%>%
  mutate(label= gsub(" ","_", species))%>%
  remove_rownames %>% 
  column_to_rownames(var="label")  %>% 
  select(area_curves)-> med_phylo #, delay


obj_med <- phylo4d(as(tree, "phylo4"), data.frame(med_phylo), match.data=TRUE)

barplot.phylo4d(obj_med) 
dotplot.phylo4d (obj_med)
gridplot.phylo4d (obj_med)
phyloSignal(obj_med)