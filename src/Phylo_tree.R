library(tidyverse);library(V.PhyloMaker); library(scales)
library(phylosignal);library(phylobase);library(ape);library(tidytree)


### Phylo tree both communitites #####
# always check family names with http://www.mobot.org/MOBOT/research/APweb/
read.csv("data/all_info.csv", sep =";") %>% # read species header data
  select(species, family)%>%
  unique()%>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>% # change species name of a species we later confirmed the correct identification
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>%# change species name of a species we later confirmed the correct identification
  separate(species, into = c("genus", "species"), sep = " ") %>%
  mutate(species = paste(genus, species),
         genus = genus,
         family = family) %>%
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
