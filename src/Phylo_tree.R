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
  filter (! species == "Kobresia myosuroides")%>% # germ = 0
  filter (! species == "Salix breviserrata")%>% # germ = 0
  filter (! species == "Sedum album")%>% # germ = 0
  filter (! species == "Gentianella campestris")%>% # germ = 0
  filter (! species == "Festuca rubra")%>% #only sowed in fellfield
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
  tree_tem2

rm(ranks1)
x11()
plot(tree_tem2$scenario.3)

write.tree(tree_tem2$scenario.3, file = "results/temperate2_tree.tree")

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

#### traits Area between curves plotted according to phylogeny#####
# temperate #
ape::read.tree("results/temperate2_tree.tree")-> tree
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
  filter (! species == "Kobresia myosuroides")%>%
  filter (! species == "Salix breviserrata")%>%
  filter (! species == "Sedum album")%>%
  filter (! species == "Gentianella campestris")%>%
  filter (! species == "Festuca rubra")%>% #only sowed in fellfield
  select(species,area_curves, delay )%>%
  mutate(label= gsub(" ","_", species))%>%
  remove_rownames %>% 
  column_to_rownames(var="label")  %>% 
  select(area_curves)%>% 
  rename(germination_shift =area_curves)%>% 
  mutate(germination_shift = -germination_shift)-> tem_phylo #,delay


obj_tem <- phylo4d(as(tree, "phylo4"), data.frame(tem_phylo), match.data=TRUE)
mat.col <- ifelse(tdata(obj_tem , "tip") < 0,  "chocolate2","deepskyblue3")
barplot.phylo4d(obj_tem, tree.ratio = 0.2,bar.col = mat.col, center=F,
                trait.bg.col = "white", show.box = F, grid.vertical = F,
                grid.horizontal = F, tip.labels = NULL, tree.ladderize =T) #rainbow(37)

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
  select(area_curves)%>% 
  rename(germination_shift =area_curves)%>% 
  mutate(germination_shift = -germination_shift)-> med_phylo #, delay


obj_med <- phylo4d(as(tree, "phylo4"), data.frame(med_phylo), match.data=TRUE)
mat.col <- ifelse(tdata(obj_med , "tip") < 0,  "chocolate2","deepskyblue3")
barplot.phylo4d(obj_med, tree.ratio = 0.2,bar.col = mat.col, center=F,
                trait.bg.col = "white", show.box = F, grid.vertical = F,
                grid.horizontal = F, tip.labels = NULL, tree.ladderize =T)

barplot.phylo4d(obj_med) 
dotplot.phylo4d (obj_med)
gridplot.phylo4d (obj_med)
phyloSignal(obj_med)