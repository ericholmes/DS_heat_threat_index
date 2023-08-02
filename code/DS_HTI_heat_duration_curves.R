## Delta Smelt heat threat index
## Heat duration curves

library(tidyverse)

load("data/DS_HTI_Wtemp_data.Rdata")

# Prepare data for analysis --------------------------------------------------------

wtempdata$jday <- as.numeric(format(wtempdata$Datetime, format = "%j"))
wtempdata$Year <- as.numeric(format(wtempdata$Datetime, format = "%Y"))
wtempdata$Month <- as.numeric(format(wtempdata$Datetime, format = "%m"))

##Analysis goals:
## 1) Define recurrence interval of "fish killing" heat waves at LIS (longest record). Apply elsewhere.
##    Plan to use flood typing methodology. First discretize heat events above some threshold. Bin by month?
## 2) Use code for flow duration curve to apply to "heat duration curve". Apply to daily, min, max, and means.

wtempply <- wtempdata[is.na(wtempdata$Param_val) == FALSE & wtempdata$Month %in% 6:10,] %>% 
  group_by(Year, jday, Month, Site_no) %>% 
  summarize(meantemp = mean(Param_val), mintemp = min(Param_val),
            maxtemp = max(Param_val), rangetemp = maxtemp - mintemp) %>% data.frame()

ggplot(wtempdata, aes(x = Datetime, Param_val, color = Site_no)) + geom_line()

ggplot(wtempply[wtempply$Month %in% 7:9,], aes(maxtemp, fill = Site_no)) +
  stat_density() + geom_vline(xintercept = 27, color = "red") + theme_bw() +
  labs(x = "Daily maximum water temperature (C)")

ggplot(wtempply[wtempply$Month %in% 7:9,], aes(mintemp, fill = Site_no)) +
  stat_density() + geom_vline(xintercept = 24, color = "red") + theme_bw() +
  labs(x = "Daily minimum water temperature (C)")

wtempply <- wtempply %>% group_by(Site_no) %>% 
  mutate(maxrank = rank(-maxtemp),
         minrank = rank(-mintemp)) %>%
  mutate(maxP = 100 * (maxrank / (length(maxtemp) + 1)),
         minP = 100 * (minrank / (length(mintemp) + 1)))

png("output/DS_HTI_heat_duraton_curve.png", height = 4, width = 6, unit = "in", res = 1000)
lisintmax <- max(c(0, unlist(wtempply[wtempply$Site_no == "LIS" & wtempply$maxtemp > 27, "maxP"])))
rvbintmax <- max(c(0, unlist(wtempply[wtempply$Site_no == "RVB" & wtempply$maxtemp > 27, "maxP"])))

ggplot(wtempply, aes(x = maxP, y = maxtemp, color = Site_no))+
  geom_line() + labs(title = "Heat duration curve - Daily max water temp") +
  xlab("% Time temp equalled or exceeded") + theme_bw() +
  ylab("Temp (C)") + coord_cartesian(ylim = c(10,32)) +
  scale_color_manual(values = c("LIS" = "firebrick", "RVB" = "darkblue")) +
  geom_hline(yintercept = 27, color = "black", linetype = 2) +
  geom_segment(aes(x = lisintmax, xend = lisintmax, y = 0, yend = 27), color = "firebrick")+
  geom_segment(aes(x = rvbintmax, xend = rvbintmax, y = 0, yend = 27))

lisintmin <- max(c(0, unlist(wtempply[wtempply$Site_no == "LIS" & wtempply$mintemp > 24, "minP"])))
rvbintmin <- max(c(0, unlist(wtempply[wtempply$Site_no == "RVB" & wtempply$mintemp > 24, "minP"])))

ggplot(wtempply, aes(x = minP, y = mintemp, color = Site_no))+
  geom_line() + labs(title = "Heat duration curve - Daily min water temp") +
  xlab("% Time temp equalled or exceeded") + theme_bw() +
  ylab("Temp (C)") + coord_cartesian(ylim = c(10,32)) +
  scale_color_manual(values = c("LIS" = "firebrick", "RVB" = "darkblue")) +
  geom_hline(yintercept = 24, color = "black", linetype = 2) +
  geom_segment(aes(x = lisintmin, xend = lisintmin, y = 0, yend = 24), color = "firebrick")+
  geom_segment(aes(x = rvbintmin, xend = rvbintmin, y = 0, yend = 24))

dev.off()