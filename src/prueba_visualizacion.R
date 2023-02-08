library(FD); library(vegan); library(FactoMineR); library(emmeans);
library(tidyverse); library(ggrepel); library(cowplot); library (ggplot2);
library (lubridate); library(binom); library (GerminaR)

# AUTUMN (test representation) ####
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  group_by (species, code, incubator, petridish, time)%>%
  filter (between(time, 1, 105)) %>%## 105 = 12/11 last check before winter
  group_by (species, code, incubator, petridish) %>%
  summarise(cumulative = sum(germinated), 
            viable = last(viable))%>%
  merge(viables_germ)%>%
  filter(viablePER>25)%>%
  filter(germPER>0)%>%
  merge(species, by = c("code", "species")) %>%
  select (code, species, region, incubator, petridish, cumulative, viable, family, habitat)%>%
  group_by (region, incubator)%>%
  summarise (cumulative = sum (cumulative),
             viable = sum(viable)) %>%
  mutate (binom.confint(cumulative, viable, methods = "wilson")) %>% 
  mutate (Season = "Autumn") -> autumn_graph

# WINTER (test representation) ####
read.csv("data/all_data.csv", sep = ";") %>%
  #filter (species == "Luzula caespitosa") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))  %>%
  group_by (species, code, incubator, petridish, date)%>% 
  filter (date %in% c("2022-04-04","2022-04-11", "2022-04-12", # # first check after winter before Tmin 2ºC
                      "2022-05-11","2022-05-12", "2022-05-18", "2022-05-19", 
                      "2022-05-25", "2022-05-26", "2022-06-01", "2022-06-02" )) %>% #undersnow + first checks Tmin 2ºC
  group_by (species, code, incubator, petridish) %>% 
  summarise(cumulative = sum(germinated), 
            viable = last(viable))%>%
  merge (viables_germ) %>%
  filter(viablePER>25)%>%
  filter(germPER>0)%>%
  merge(species, by = c("code", "species")) %>%
  select(code, species, region, incubator, petridish, cumulative, viable, family, habitat)%>%
  group_by ( incubator)%>%
  summarise (cumulative = sum (cumulative),
             viable = sum(viable)) %>%
  mutate (binom.confint(cumulative, viable, methods = "wilson")) %>% 
  mutate (Season = "Winter")-> winter_graph  

# COLD GERM  (test representation) ####
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>% 
  group_by (species, code, incubator, petridish, time) %>% 
  filter (incubator == "Fellfield") %>% 
  filter (between(time, 257, 285)) %>% ## amount of days calculated from the dates 257 = 13/4 and 285 = 11/05
  group_by (species, code, incubator, petridish) %>%
  summarise(seeds_germ = sum(germinated)) %>%
  merge (germ_AF) %>%
  select(species, code, incubator, petridish, seeds_germ, germ_AW, viable, viablePER, germPER)  -> Earlyseason_F

read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))  %>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  group_by (species, code, incubator, petridish, time) %>% 
  filter (incubator == "Snowbed") %>% 
  filter (between(time, 308, 322)) %>% ## amount of days calculated from the dates 308 = 03/06 and 322 = 17/06
  group_by (species, code, incubator, petridish) %>%
  summarise(seeds_germ = sum(germinated)) %>%
  merge (germ_AF) %>%
  select(species, code, incubator, petridish, seeds_germ, germ_AW, viable, viablePER, germPER) -> Earlyseason_S

rbind(Earlyseason_F, Earlyseason_S) %>%
  mutate (cumulative = seeds_germ)%>%
  arrange(species) %>%
  filter(viablePER>25)%>%
  filter(germPER>0) %>%
  merge(species, by = c("code", "species")) %>%
  select(code, species, region, incubator, petridish, cumulative, viable, family, habitat)%>%
  group_by (incubator)%>%
  summarise (cumulative = sum (cumulative),
             viable = sum(viable)) %>%
  mutate (binom.confint(cumulative, viable, methods = "wilson")) %>% 
  mutate (Season = "Early season")-> Earlyseason_graph  

# WARM germ (test representation) ####
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>% 
  group_by (species, code, incubator, petridish, time) %>% 
  filter (incubator == "Fellfield") %>% 
  filter (between(time, 286, 430)) %>% ## amount of days calculated from the dates 286 = 12/05 and 430 = 19/09
  group_by (species, code, incubator, petridish) %>%
  summarise(seeds_germ = sum(germinated)) %>%
  merge (germ_AF) %>%
  select(species, code, incubator, petridish, seeds_germ, germ_AW, viable, viablePER, germPER)  -> Lateseason_F

read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))  %>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  group_by (species, code, incubator, petridish, time) %>% 
  filter (incubator == "Snowbed") %>% 
  filter (between(time, 323, 430)) %>% ## amount of days calculated from the dates 323 = 18/06 and 430 = 19/09
  group_by (species, code, incubator, petridish) %>%
  summarise(seeds_germ = sum(germinated)) %>%
  merge (germ_AF) %>%
  select(species, code, incubator, petridish, seeds_germ, germ_AW, viable, viablePER, germPER) -> Lateseason_S

rbind(Lateseason_F, Lateseason_S) %>%
  mutate (cumulative = seeds_germ)%>%
  arrange(species) %>%
  filter(viablePER>25)%>%
  filter(germPER>0)%>%
  merge(species, by = c("code", "species")) %>%
  select(code, species, region, incubator, petridish, cumulative, viable, family, habitat)%>%
  group_by (region)%>%
  summarise (cumulative = sum (cumulative),
             viable = sum(viable)) %>%
  mutate (binom.confint(cumulative, viable, methods = "wilson")) %>% 
  mutate (Season = "Late season")-> Lateseason_graph  

# individual season graphs ####
x11()
ggplot()+
  geom_point(data=autumn_graph, aes(incubator, mean, color=incubator), size =2) +
  facet_grid (.~region) +
  geom_errorbar(data= autumn_graph, aes(incubator, mean, ymin = lower, ymax = upper, color = incubator), size =1.2) +
  #scale_y_continuous (limits = c(0,1), breaks = seq (0, 1, by= 0.2)) +
  theme_bw(base_size = 16)+
  theme(legend.position = "right", 
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10),
        legend.background = element_rect(fill="transparent",colour=NA)) #change legend text font size theme_bw()
##### ALL seasons, total amount of viables ##########
read.csv("data/all_data.csv", sep = ";") %>%
  group_by(species, code, incubator, petridish) %>%
  summarise(germinated = sum(germinated),
            viable = last(viable)) %>%     
  merge(species, by = c("code", "species")) %>%
  merge(viables_germ)%>%
  filter(viablePER>25)%>%
  filter(germPER>0)%>%
  group_by(region, incubator)%>%
  summarise(viable = sum(viable)) -> viable_graph
# autumn
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  group_by (species, code, incubator, petridish, time)%>%
  filter (between(time, 1, 105)) %>%
  group_by (species, code, incubator, petridish) %>%
  summarise(cumulative = sum(germinated))%>%
  merge(viables_germ)%>%
  filter(viablePER>25)%>%
  filter(germPER>0)%>%
  merge(species, by = c("code", "species")) %>%
  select (code, species, region, incubator, petridish, cumulative, family, habitat)%>%
  group_by (region, incubator)%>%
  summarise (cumulative = sum (cumulative)) %>%
  merge(viable_graph) %>%
  mutate (binom.confint(cumulative, viable, methods = "wilson")) %>% 
  mutate (Season = "Autumn") -> autumn_graph
#winter
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))  %>%
  group_by (species, code, incubator, petridish, date)%>% 
  filter (date %in% c("2022-04-04","2022-04-11", "2022-04-12", # # first check after winter before Tmin 2ºC
                      "2022-05-11","2022-05-12", "2022-05-18", "2022-05-19", 
                      "2022-05-25", "2022-05-26", "2022-06-01", "2022-06-02" )) %>% #undersnow + first checks Tmin 2ºC
  group_by (species, code, incubator, petridish) %>% 
  summarise(cumulative = sum(germinated))%>%
  merge (viables_germ) %>%
  filter(viablePER>25)%>%
  filter(germPER>0)%>%
  merge(species, by = c("code", "species")) %>%
  select(code, species, region, incubator, petridish, cumulative, family, habitat)%>%
  group_by (region, incubator)%>%
  summarise (cumulative = sum (cumulative)) %>%
  merge(viable_graph) %>%
  mutate (binom.confint(cumulative, viable, methods = "wilson")) %>% 
  mutate (Season = "Winter")-> winter_graph  

# early season
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>% 
  group_by (species, code, incubator, petridish, time) %>% 
  filter (incubator == "Fellfield") %>% 
  filter (between(time, 257, 285)) %>% ## amount of days calculated from the dates 257 = 13/4 and 285 = 11/05
  group_by (species, code, incubator, petridish) %>%
  summarise(seeds_germ = sum(germinated)) %>%
  merge (germ_AF) %>%
  select(species, code, incubator, petridish, seeds_germ, germ_AW, viable, viablePER, germPER)  -> Earlyseason_F

read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))  %>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  group_by (species, code, incubator, petridish, time) %>% 
  filter (incubator == "Snowbed") %>% 
  filter (between(time, 308, 322)) %>% ## amount of days calculated from the dates 308 = 03/06 and 322 = 17/06
  group_by (species, code, incubator, petridish) %>%
  summarise(seeds_germ = sum(germinated)) %>%
  merge (germ_AF) %>%
  select(species, code, incubator, petridish, seeds_germ, germ_AW, viable, viablePER, germPER) -> Earlyseason_S
rbind(Earlyseason_F, Earlyseason_S) %>%
  mutate (cumulative = seeds_germ)%>%
  arrange(species) %>%
  filter(viablePER>25)%>%
  filter(germPER>0) %>%
  merge(species, by = c("code", "species")) %>%
  select(code, species, region, incubator, petridish, cumulative, family, habitat)%>%
  group_by (region, incubator)%>%
  summarise (cumulative = sum (cumulative)) %>%
  merge(viable_graph)%>%
  mutate (binom.confint(cumulative, viable, methods = "wilson")) %>% 
  mutate (Season = "Spring")-> Earlyseason_graph 

# late season
read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))%>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>% 
  group_by (species, code, incubator, petridish, time) %>% 
  filter (incubator == "Fellfield") %>% 
  filter (between(time, 286, 430)) %>% ## amount of days calculated from the dates 286 = 12/05 and 430 = 19/09
  group_by (species, code, incubator, petridish) %>%
  summarise(seeds_germ = sum(germinated)) %>%
  merge (germ_AF) %>%
  select(species, code, incubator, petridish, seeds_germ, germ_AW, viable, viablePER, germPER)  -> Lateseason_F

read.csv("data/all_data.csv", sep = ";") %>%
  mutate(date = strptime(as.character(date), "%d/%m/%Y"))  %>%
  mutate(time = as.numeric(as.Date(date)) - min(as.numeric(as.Date(date))))%>%
  group_by (species, code, incubator, petridish, time) %>% 
  filter (incubator == "Snowbed") %>% 
  filter (between(time, 323, 430)) %>% ## amount of days calculated from the dates 323 = 18/06 and 430 = 19/09
  group_by (species, code, incubator, petridish) %>%
  summarise(seeds_germ = sum(germinated)) %>%
  merge (germ_AF) %>%
  select(species, code, incubator, petridish, seeds_germ, germ_AW, viable, viablePER, germPER) -> Lateseason_S
rbind(Lateseason_F, Lateseason_S) %>%
  mutate (cumulative = seeds_germ)%>%
  arrange(species) %>%
  filter(viablePER>25)%>%
  filter(germPER>0)%>%
  merge(species, by = c("code", "species")) %>%
  select(code, species, region, incubator, petridish, cumulative, family, habitat)%>%
  group_by (region, incubator)%>%
  summarise (cumulative = sum (cumulative)) %>%
  merge(viable_graph)%>%
  mutate (binom.confint(cumulative, viable, methods = "wilson")) %>% 
  mutate (Season = "Summer")-> Lateseason_graph  


graph <- rbind(autumn_graph, winter_graph, Earlyseason_graph, Lateseason_graph) %>% 
  mutate (Season = factor(Season, levels = c( "Autumn", "Winter", "Spring", "Summer"))) #reorder levels
#rm(graph_date)
x11()
ggplot()+
  geom_point(data=graph, aes(Season, mean, color=incubator), size =2) +
  facet_grid (.~region) +
  geom_errorbar(data= graph, aes(Season, mean, ymin = lower, ymax = upper, color = incubator), size =1) +
  geom_line(data= graph, aes(as.numeric(Season), mean, color = incubator), size =1) +
  scale_y_continuous (limits = c(0,1), breaks = seq (0, 1, by= 0.2)) +
  theme_bw(base_size = 16)+
  theme(legend.position = "right", 
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10),
        legend.background = element_rect(fill="transparent",colour=NA)) #change legend text font size theme_bw()
## Total germination (test representation) ####
viables %>% 
  merge(finalgerm) %>% 
  filter(viablePER>25)%>%
  filter(germPER>0)%>%
  select (species, code, incubator, petridish, viable, seeds_germ) %>%
  merge(species, by = c("code", "species")) %>%
  mutate (cumulative = seeds_germ) %>%
  select (code, species, region, incubator, petridish, cumulative, viable, family, habitat)%>%
  group_by (incubator)%>%
  summarise (cumulative = sum (cumulative),
             viable = sum(viable)) %>%
  mutate (binom.confint(cumulative, viable, methods = "wilson")) -> total_germ
# graph #### #### WORK IN PROGRESS, FULL Y-AXIS LENGTH (0-1)###
x11()
ggplot()+
  geom_point(data=total_germ, aes(incubator, mean, color=incubator), size =2) +
 # facet_grid (.~region) +
  geom_errorbar(data= total_germ, aes(incubator, mean, ymin = lower, ymax = upper, color = incubator), size =1.2) +
  scale_y_continuous (limits = c(0,1), breaks = seq (0, 1, by= 0.2)) +
  #geom_line(data= graph_date, aes(as.numeric(season), mean, color = incubator), size =1) +
  theme_bw(base_size = 16)+
  theme(legend.position = "right", 
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10),
        legend.background = element_rect(fill="transparent",colour=NA)) #change legend text font size theme_bw()

# date_graph put together
graph_date <- rbind(Mid_november_graph, Winter_graph, Mid_june_graph, Mid_september_graph) %>% 
  mutate (season = factor(season, levels = c( "Mid-November", "Winter", "Mid-June", "Mid-September"))) #reorder levels
#rm(graph_date)
x11()
ggplot()+
  geom_point(data=graph_date, aes(season, mean, color=incubator), size =2) +
  facet_grid (.~region) +
  geom_errorbar(data= graph_date, aes(season, mean, ymin = lower, ymax = upper, color = incubator), size =1) +
  geom_line(data= graph_date, aes(as.numeric(season), mean, color = incubator), size =1) +
  theme_bw()