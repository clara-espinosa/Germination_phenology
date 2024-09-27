library (tidyverse); library(scales)
library(tidyverse);library(V.PhyloMaker)
library(phylosignal);library(phylobase);library(ape);library(tidytree)
library (tidyverse); library(scales)

##### Calculate GERMINATION SHIFT (areas between curves Fellfield-Snowbed) for each species population ####
### seeds germ + viables/species and incubator 
read.csv("data/clean data.csv", sep = ";") %>%
  group_by(species, code, incubator, petridish) %>% # grouping variables
  summarise(total_germ = sum(germinated)) %>% # sum total germination
  merge(viables) %>% # merge viables object created above to filter raw data
  merge(species) %>% # join species header data
  group_by (species, code, incubator) %>% # new grouping variables
  summarise (viable = sum (viable))-> viables_pop # viables per cspeceis and population

# CALCULATED INDIVIDUALLY PER EACH SPECIES AND ADDED MANUALLY TO all_info.csv as ABC_clean_data.csv
read.csv("data/clean data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  spread(time, germinated, fill = 0) %>% # wide format for dates, and fill Na with 0
  mutate(date = as.POSIXct(date))%>%
  gather ("time", "germinated", 9: last_col()) %>% # back in long format from 9th colum to the last
  mutate (time = as.numeric(time)) %>%
  arrange (species, code, incubator, petridish, time)%>% # sort row observations this way
  merge (species) %>%
  dplyr::group_by(community, species, code,  incubator, time) %>%
  dplyr::summarise(germinated = sum(germinated)) %>%
  dplyr::mutate(germinated = cumsum(germinated)) %>%
  merge(viables_pop) %>%
  mutate(germination = germinated/viable)  %>%
  filter (community == "Temperate") %>%
  filter(species == "Sesleria caerula")%>% 
  select( species, code, incubator, community, time, germination)%>%
  #dplyr::group_by(species, incubator, time) %>%   # command line for 2 populations sp
  #dplyr::summarise(germination = mean(germination)) %>%  # command line for 2 populations sp
  spread (incubator, germination) %>%
  data.frame()-> df

#calculate area by trapezoidal rule for 1 pop species
# https://stackoverflow.com/questions/24742677/how-to-measure-area-between-2-distribution-curves-in-r-ggplot2
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

##### Fig 3A visualization density plots #####

read.csv("data/all_info.csv", sep = ",") %>%
  select(community, species, family, ABC_clean_data) %>%
  group_by (community, species, family) %>%
  summarise(area_curves = mean(ABC_clean_data))%>%
  filter (! area_curves==0) %>% # remove species with 0 difference meaning they had 0 germ
  #filter (community == "Mediterranean")%>%
  mutate(community = as.factor (community))%>%
  mutate(community = fct_relevel (community, "Temperate","Mediterranean"))%>%
  mutate (community = recode (community, "Temperate" = "Temperate (N = 38)", "Mediterranean" = "Mediterranean (N = 21)"))%>%
  ggplot(aes(x= -area_curves))+  #
  geom_density(alpha =0.2, fill = "azure4", color="azure4" ) + #, position="stack"
  geom_vline(xintercept = 0, linetype = "dashed", size =1) +
  scale_x_continuous( limits = c(-200, 100), breaks = seq (-200, 100, by= 100)) +
  scale_y_continuous(labels = percent,limits = c(0,0.0105)) +
  facet_wrap(~community)+
  theme_classic(base_size = 16) +
  labs ( x = "Germination shift", y = "Values density")+
  theme (plot.title = element_text (hjust = 0.5,face = "bold",size = 24), #hjust = 0.5,
         plot.subtitle = element_text(face = "bold",size = 16),
         axis.title.y = element_text (size=14),
         axis.text.y = element_text (size = 13),
         axis.title.x = element_text (size=14),
         axis.text.x= element_text (size = 13, color = "black"),
         strip.text = element_text( size = 18, hjust = 0),
         strip.background = element_blank(), 
         panel.background = element_blank(), #element_rect(color = "black", fill = NULL), 
         #panel.grid = element_blank(),
         legend.title = element_text (size =14),
         legend.text = element_text (size =14),
         legend.position = "bottom", # legend.position = c(0.85, 0.5),
         legend.box.background = element_rect(color = "black", size = 2))-> fig3a;fig3a

ggsave(filename = "fig3a", plot =fig3a, path = "results/Figures/", 
       device = "png", dpi = 600)
### Phylo tree temperate community #####
# always check family names with http://www.mobot.org/MOBOT/research/APweb/
read.csv("data/all_info.csv", sep =";") %>% # read species header data
  select(community, species, family)%>%
  unique()%>%
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
read.csv("data/all_info.csv", sep =",") %>% # read species header data
  select(community, species, family)%>%
  unique()%>%
  mutate(species= str_replace(species, "Minuartia CF", "Minuartia arctica"))%>%
  mutate(species= str_replace(species, "Sedum album cf", "Sedum album")) %>%
  filter(community == "Mediterranean") %>%
  select (species, family) %>%
  unique %>%
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
  tree_med

rm(ranks1)
x11()
plot(tree_med$scenario.3)

write.tree(tree_med$scenario.3, file = "results/mediterranean_tree.tree")

#### Fig 3B traits Area between curves plotted according to phylogeny#####
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
  filter (! species == "Kobresia myosuroides")%>% # species with 0 germination
  filter (! species == "Salix breviserrata")%>%# species with 0 germination
  filter (! species == "Sedum album")%>%# species with 0 germination
  filter (! species == "Gentianella campestris")%>%# species with 0 germination
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
barplot.phylo4d(obj_tem, tree.ratio = 0.2,bar.col = mat.col, center=F,tip.cex = 1.5,
                trait.bg.col = "white", show.box = F, grid.vertical = F,
                grid.horizontal = F, tip.labels = NULL, tree.ladderize =T) #rainbow(37)

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

# plot panel arranged in canva 