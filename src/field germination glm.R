library(tidyverse);library (rstatix);library (stringr); library (vegan)

##### Field germination for graph + glms ####
#dataframe for testing field germination differences
read.csv("data/field_germination.csv") %>% # read raw field germination data
  select(retrieval_season, community, site_buried,microhabitat_buried, species, seeds_initial, bag1:bag3)%>%
  convert_as_factor(retrieval_season, community, site_buried, microhabitat_buried, species) %>%
  gather(bag, germ, bag1:bag3)%>% # transform data from wide to long format
  mutate(seeds_initial = 10) %>%
  filter(retrieval_season == "Spring_23" ) -> earlyseason_df # create dataframe for spring removal

read.csv("data/field_germination.csv") %>% # read raw field germination data
  select(retrieval_season, community, site_buried,microhabitat_buried, species, seeds_initial, bag1:bag3)%>%
  convert_as_factor(retrieval_season, community, site_buried, microhabitat_buried, species) %>%
  gather(bag, germ, bag1:bag3)%>% # transform data from wide to long format
  mutate(seeds_initial = 10) %>%  
  filter(retrieval_season == "Autumn_23")-> lateseason_df # create one dataframe for autumn removal


### Get GLM coefficients with a function for each species 

glms <- function(x) {
  glm(cbind(germ, seeds_initial - germ) ~ microhabitat_buried, family = "binomial", data = x) -> m1  
  broom::tidy(m1)
}

earlyseason_df%>%
  group_by(species) %>%
  do(glms(.))%>%
  print(n=36)

lateseason_df%>%
  group_by(species) %>%
  do(glms(.))%>%
  print(n=36)
