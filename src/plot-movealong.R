library(FD); library(vegan); library(FactoMineR); library(emmeans);
library(tidyverse); library(ggrepel); library(cowplot)
theme_set(theme_cowplot(font_size = 10)) 

## GERMINATION EXPERIMENTS

# Read and prepare the germination data

read.csv("data/example/Germination.csv") %>%
  group_by(species, country, date) %>%
  summarise(germination = (sum(germinated)/sum(viable)) * 100,
            establishment = (sum(seedling)/sum(viable)) * 100) %>%
  merge(read.csv("data/example/Germination.csv"), by = c("species", "country", "date")) %>% 
  group_by(species, country, dish) %>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))),
         germination = cumsum(germination),
         germinated = cumsum(germinated), 
         establishment = cumsum(establishment),
         seedling = cumsum(germinated)) %>%
  merge(read.csv("data/example/Incubators.csv"), by = c("date", "country"), all.x = T) %>%
  arrange(species, country, dish, date) ->
  dfg

# Set a function to create the germination plots

fg <- function(df) {
  ggplot(data = df, aes(x = time, y = germinated/viable * 100)) +
    facet_wrap(species ~ country, scales = "free_x", ncol = 5) +
    geom_point(size = 1, alpha = 0.2, color = "turquoise4") +
    geom_line(data = filter(df, dish == "A"), 
              aes(y = germination), color = "turquoise4", size = 1) +
    coord_cartesian(ylim = c(0, 100)) +
    labs(x = "Time (days)",y = "Germination (%)") +
    geom_area(data = filter(df, dish == "A"), 
              aes(y = temperature/max(temperature)*100), alpha = 0.4, fill = "red") +
    scale_y_continuous(sec.axis = sec_axis(~.*max(df$temperature)/100, 
                                           name = "Temperature (?C)")) +
    theme(strip.text = element_text(face = "italic", size = 7))}

# Create the plots, divided for 2 pages

dfg %>%
  filter(species %in% unique(dfg$species)[1:19]) %>%
  fg -> FigGerminationA

dfg %>%
  filter(species %in% unique(dfg$species)[20:43]) %>%
  fg -> FigGerminationB

# Save the plots

ggsave(filename = "results/FigGerminationA.png", FigGerminationA, path = NULL, 
       scale = 1, width = 180, height = 180, units = "mm", dpi = 600)

ggsave(filename = "results/FigGerminationB.png", FigGerminationB, path = NULL, 
       scale = 1, width = 180, height = 180, units = "mm", dpi = 600)
