library(tidyverse);library(ggpubr);

## WEEKLY MEANS (2009-2018) FOR PICOS ####
read.csv("data/temp_picos_raw.csv", sep = ";") %>%
  mutate(Time = strptime(as.character(Time), "%d/%m/%Y %H:%M"))%>% #specify format of TIME variable
  mutate(Time = as.POSIXct(Time, tz = "UTC")) %>%
  mutate(Year = lubridate::year(Time)) %>% 
  mutate(Month = lubridate::month(Time)) %>% 
  mutate(Day = lubridate::day(Time)) %>% 
  mutate(Hour = lubridate::hour(Time)) %>% 
  group_by(Site, Day = lubridate::floor_date(Time, "day")) %>%
  summarise(T = mean(Temperature), X = max(Temperature), N = min(Temperature), n = length(Time)) %>% # Daily mean, max, min
  mutate(Snow = ifelse(X < 0.5 & N > -0.5, 1, 0)) %>% # Day is snow day or not
  mutate(FreezeThaw = ifelse(X > 0.5 & N < -0.5, 1, 0)) %>% # Day with freeze-thaw cycles
  mutate(FDD = ifelse(T < 0, T, 0)) %>% # Freezing degrees per day
  mutate(GDD = ifelse(T >= 1, T, 0)) %>% # Growing degrees day per month https://link.springer.com/article/10.1007/s00035-021-00250-1
  mutate (week = strftime(Day, format = "%V")) %>%
  select(Site, T, X, N, Snow, FDD, GDD, week)%>%
  group_by (Site, week) %>%
  summarise_all(mean, na.rm = TRUE)%>% 
  mutate (Site = as.factor(Site)) %>%
  mutate(week = as.numeric(week)) %>%
  as.data.frame() -> weekly_picos

#write.csv(weekly_picos, "data/weekly_picos.csv") 
# modification in excel to add week order and  extend period to match experiment length !!

##graph to compare with temperature programs ####
x11()
read.csv("data/weekly_picos_graph.csv", sep = ";") %>%
  ggplot(aes ()) +
  geom_line (aes (x=order, y=N, colour = Site), size =1.25) + #order for sort as experiment
  geom_line (aes (x=order, y=X, colour = Site), size =1.25) +
  geom_ribbon (aes (x=order, ymin =N, ymax=X, colour = Site, fill = Site), alpha =0.3) +
  scale_color_manual (name = "Field sites", values = c ("Fellfield site"= "darkgoldenrod1" , "Snowbed site" = "forestgreen")) +
  scale_fill_manual (name = "Field sites", values = c ("Fellfield site"= "darkgoldenrod1" , "Snowbed site" = "forestgreen")) +
  scale_y_continuous (limits = c(-3,25), breaks = seq (0, 25, by= 5)) +
  labs (title = "Soil temperature records 2008 - 2019", y= "Temperature ºC", x = "Time (weeks)", tag = "A") + #
  theme_classic(base_size = 16) +
  theme (plot.title = element_text ( size = 24), #hjust = 0.5,
         axis.title.y = element_text (size=20), 
         axis.text.y = element_text (size = 18),
         axis.title.x = element_blank(), 
         axis.text.x= element_blank(),
         axis.ticks.x = element_blank(),
         plot.tag.position = c(0,1),
         legend.title = element_text (size =14),
         legend.text = element_text (size =13),
         legend.background = element_rect(fill="transparent"),
         legend.position = c(0.9, 0.28)) + # legend.position = "top") + 
  geom_hline(yintercept=0, linetype ="dashed", size =1, colour = "red") -> fig2a; fig2a

ggarrange(fig2a, fig2b, ncol =1, nrow= 2,common.legend = FALSE)
#ggsave(ggarrange(fig2a, fig2b,  ncol =1, nrow= 2,common.legend = FALSE), file="fig2.png", width = 210, height = 297, units = "mm")
