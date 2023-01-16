library(tidyverse)

### Phylo tree
# always check family names with http://www.mobot.org/MOBOT/research/APweb/
read.csv("data/species.csv", sep =";") %>%
  mutate(species= str_replace(species, "Cerastium sp", "Cerastium pumilum"))%>%
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
library(V.PhyloMaker)
phylo.maker(sp.list = ranks1,
            tree = GBOTB.extended, 
            nodes = nodes.info.1, 
            scenarios = "S3") ->
  tree

rm(ranks1)
x11()
plot(tree$scenario.3)

write.tree(tree$scenario.3, file = "results/tree.tree")
