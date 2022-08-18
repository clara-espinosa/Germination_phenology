library(FD); library(vegan); library(FactoMineR); library(emmeans);
library(tidyverse); library(ggrepel); library(cowplot)
theme_set(theme_cowplot(font_size = 10)) 

## GERMINATION EXPERIMENTS

# Read and prepare the germination data
temp <- read.csv("data/date_temp.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))

dfg <- read.csv("data/R long data.csv", sep = ";") %>%
  group_by(species, condition, code, date) %>%
  summarise(germination = (sum(germinated)/sum(viable)) * 100) %>%
  merge(read.csv("data/R long data.csv", sep = ";"), by = c("species", "condition", "code", "date")) %>% 
  group_by(species, condition, code) %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  mutate(germination = cumsum(germination),
         germinated = cumsum(germinated)) %>%
  merge(temp, by = c("date", "condition"), all.x = T) %>%
  arrange(species, condition, code, time) 
str (dfg)

# Set a function to create the germination plots

fg <- function(df) {
  ggplot(data = df, aes(x = time, y = germinated/viable * 100)) +
    facet_wrap(species ~ condition, scales = "free_x", ncol = 2) +
    geom_point(size = 1, alpha = 0.2, color = "turquoise4") +
    geom_line(data = filter(df, petridish == "2"), 
              aes(y = germination), color = "turquoise4", size = 1) +
    coord_cartesian(ylim = c(0, 100)) +
    labs(x = "Time (days)",y = "Germination (%)") +
    geom_area(data = filter(df, petridish == "1"), 
              aes(y = Tmean), alpha = 0.4, fill = "red") +
    scale_y_continuous(sec.axis = sec_axis(~.*max(df$Tmean)/100, 
                                           name = "Temperature (ºC)")) +
    theme(strip.text = element_text(face = "italic", size = 7))}


# Create the plots, divided for 4 pages

dfg %>%
  filter(species %in% unique(dfg$species)[4:8]) %>% # 12 species
  fg -> FigGerminationA

dfg %>%
  filter(species %in% unique(dfg$species)[13:26]) %>%
  fg -> FigGerminationB

dfg %>%
  filter(species %in% unique(dfg$species)[27:40]) %>%
  fg -> FigGerminationC

dfg %>%
  filter(species %in% unique(dfg$species)[41:55]) %>%
  fg -> FigGerminationD
# Save the plots

ggsave(filename = "results/FigGerminationA.png", FigGerminationA, path = NULL, 
       scale = 1, width = 180, height = 180, units = "mm", dpi = 600)

ggsave(filename = "results/FigGerminationB.png", FigGerminationB, path = NULL, 
       scale = 1, width = 180, height = 180, units = "mm", dpi = 600)
