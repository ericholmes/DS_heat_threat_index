# DS_HTI along NDFS axis from LIS to RVB
library(tidyverse)

source("code/useful_functions.R")

NDFSwqcont <- read.csv("C:/Users/eholmes/OneDrive - California Department of Water Resources/Documents/NDFS/data/NDFA_synthesis_wqcont_combined_all.csv")
stations <- data.frame(StationCode = c("RCS", "WWT", "RD22", "I80", "LIS", "STTD","BL5", "PRS", "LIB", "RYI", "RVB"),
                       Dist = c(0, 7.7, 10.1, 17.39, 24.31, 33.28, 38.31, 40.17, 41.47, 43.58, 47.62))

NDFSwtemp <- NDFSwqcont[is.na(NDFSwqcont$WaterTemp) == F, c("StationCode", "DateTime","WaterTemp")]

dput(unique(NDFSwtemp$StationCode))

NDFSwtemp <- NDFSwtemp[NDFSwtemp$StationCode %in% c("LIBCUT", "LIB", "LIS","RVB", 
                                                    "RYI", "SDI", "SRH", "STTD"),]

##Combine with shag slough,USGS DWSC 11455142
usgssites <- data.frame(StationCode = c("Cache", "DWSC"),
           IDs = c( 11455280, 11455142))
##janky for-loop, replace with vectorized apply function?
datparams <- data.frame()

for(i in usgssites$IDs){
  print(i)
  tempdat <- downloadNWIS(i, "00010", "2010-01-01", "2023-10-01")
  tempdat$StationCode = i
  datparams <- rbind(datparams, tempdat)
}

datparamsub <- datparams[, c(2,9,6)]
colnames(datparamsub) <- c("StationCode", "DateTime","WaterTemp")

# Prepare data fields -----------------------------------------------------
NDFSwtemp$DateTime <- as.POSIXct(NDFSwtemp$DateTime, format = "%Y/%m/%d %H:%M:%S")

NDFSwtemp <- rbind(NDFSwtemp, datparamsub)
NDFSwtemp$Year <- format(NDFSwtemp$DateTime, format = "%Y")
NDFSwtemp$Month <- as.numeric(format(NDFSwtemp$DateTime, format = "%m"))
NDFSwtemp$jday <- as.numeric(format(NDFSwtemp$DateTime, format = "%j"))
NDFSwtemp$Date <- as.Date(NDFSwtemp$DateTime)

ggplot(NDFSwtemp, aes(x = DateTime, y = WaterTemp)) + geom_line() + facet_wrap(StationCode ~ .)

# Create heat duration curves and percent exceedance values ---------------

wtempply <- NDFSwtemp[is.na(NDFSwtemp$WaterTemp) == FALSE & NDFSwtemp$Month %in% 6:10,] %>% 
  group_by(Year, jday, Month, StationCode) %>% 
  summarize(meantemp = mean(WaterTemp, na.rm = T), mintemp = min(WaterTemp),
            maxtemp = max(WaterTemp), rangetemp = maxtemp - mintemp) %>% data.frame()

ggplot(wtempply[wtempply$Month %in% 7:9,], aes(maxtemp, fill = StationCode)) +
  stat_density(alpha = .9, color = "white") + geom_vline(xintercept = 27, color = "red") + theme_bw() +
  labs(x = "Daily maximum water temperature (C)")

ggplot(wtempply[wtempply$Month %in% 7:9,], aes(mintemp, fill = StationCode)) +
  stat_density() + geom_vline(xintercept = 24, color = "red") + theme_bw() +
  labs(x = "Daily minimum water temperature (C)")

wtempply <- wtempply %>% group_by(StationCode) %>% 
  mutate(maxrank = rank(-maxtemp),
         minrank = rank(-mintemp)) %>%
  mutate(maxP = 100 * (maxrank / (length(maxtemp) + 1)),
         minP = 100 * (minrank / (length(mintemp) + 1)))

# Plot heat duration curves for daily mins and maxes ----------------------

xint <- data.frame()
for(i in unique(wtempply$StationCode)){
  print(i)
  xint <- rbind(xint, data.frame("StationCode" = i, 
                                 "intmax" = max(c(0, unlist(wtempply[wtempply$StationCode == i & 
                                                                       wtempply$maxtemp > 27, "maxP"]))),
                                 "intmin" = max(c(0, unlist(wtempply[wtempply$StationCode == i & 
                                                                       wtempply$mintemp > 24, "minP"])))))
}

xint$Sitelabel <- c("CACHE", "RVB", "RYI", "SDI", "SRH", "LIS", "LIB", "STTD", 
                     "DWSC", "LIBCUT")

png("output/DS_HTI_heat_duration_curve_NDFSaxis_%02d.png", height = 4, width = 6, unit = "in", res = 1000)

ggplot(wtempply, aes(x = maxP, y = maxtemp, color = StationCode))+
  geom_line() + labs(title = "Heat duration curve - Daily max water temp") +
  xlab("% Time temp equalled or exceeded") + theme_bw() +
  ylab("Temp (C)") + coord_cartesian(ylim = c(10,32)) +
  geom_hline(yintercept = 27, color = "black", linetype = 2) +
  scale_color_viridis_d(option = "D") +
  geom_segment(data = xint, aes(x = intmax, xend = intmax, y = 0, yend = 27))
   
ggplot(wtempply, aes(x = minP, y = mintemp, color = StationCode))+
  geom_line() + labs(title = "Heat duration curve - Daily min water temp") +
  xlab("% Time temp equalled or exceeded") + theme_bw() +
  ylab("Temp (C)") + coord_cartesian(ylim = c(10,32)) +
  geom_hline(yintercept = 24, color = "black", linetype = 2) +
  scale_color_viridis_d(option = "D") +
  geom_segment(data = xint, aes(x = intmin, xend = intmin, y = 0, yend = 24))

cowplot::plot_grid(ggplot(xint, aes(x = reorder(Sitelabel, -intmax), y = intmax)) + geom_bar(stat = "identity") +
                     theme_bw() + labs(x = NULL, y = "daily max > 27°C (%)"),
                   
                   ggplot(xint, aes(x = reorder(Sitelabel, -intmax), y = intmin)) + geom_bar(stat = "identity") +
                     theme_bw() + labs(x = NULL, y = "daily min > 24°C (%)"),
                   nrow = 2)

dev.off()

# Spatial -----------------------------------------------------------------

# Leaflet
library(sp)
library(leaflet)
sites <- read.csv("data/NDFS_longitudinal_dist.csv")
sites <- rbind(sites[,c(1:3)], data.frame("Sitename" = c("LIBCUT", "SDI", 11455280, 11455142),
                          "Lat" = c(38.3289, 38.0934, 38.2755, 38.341666),
                          "Lon" = c(-121.6675, -121.736, -121.7092, -121.64388)))
sites <- merge(sites, xint, by.x = "Sitename", by.y = "StationCode", all.y = T)
sites <- sites[is.na(sites$Lat) == F, ]


pal <- colorNumeric(c("blue", "red", "black"), c(0,20),
                    na.color = "transparent")

leaflet(sites) %>% addTiles() %>%
  addCircleMarkers(~Lon, ~Lat,
                   color = ~pal(intmax),
                   label = ~paste0(signif(intmax, 3), "%")) %>% 
  addLegend(pal = pal, values = c(0,20),
            title = "daily max > 27°C<br>% exceedance")

leaflet(sites) %>% addTiles() %>%
  addProviderTiles("Esri.NatGeoWorldMap", group = "OSM") %>%
  addProviderTiles("Esri.WorldImagery", group = "Imagery") %>% 
  addCircleMarkers(~Lon, ~Lat,
                   color = ~pal(intmin),
                   label = ~paste0(signif(intmin, 3), "%")) %>% 
  addLegend(pal = pal, values = c(0,20),
            title = "daily min > 24°C<br>% exceedance") %>% 
  addLayersControl(
    baseGroups = c("Imagery", "OSM"),
    options = layersControlOptions(collapsed = FALSE))

# static maps
library(sf)
library(ggrepel)
library(ggthemes)

ndelta <- read_sf("data/ndelta.geojson")
sites_sf <- st_as_sf(sites, coords = c("Lon", "Lat"), crs = "+proj=longlat +datum=WGS84 +no_defs +type=crs")
set.seed(1)
(imaxplot <- ggplot(sites_sf) + geom_sf(data = ndelta, fill = "aquamarine3", color = "aquamarine3") + 
    geom_sf(size = 5.5, color = "white") + theme_bw() +
    theme_map() + labs(title = "daily max > 27°C % exceedance") +
    geom_sf(aes(color = intmax), size = 4) + theme(legend.position = "none") +
    ylim(min(sites$Lat), max(sites$Lat)) +
    geom_label_repel(data = sites_sf,
      aes(label = Sitelabel, geometry = geometry),
      stat = "sf_coordinates",
      min.segment.length = 0) +
    scale_color_viridis_c(option = "C", limits = c(0,max(c(sites$intmax, sites$intmin)))))
set.seed(1)
(iminplot <- ggplot(sites_sf) + geom_sf(data = ndelta, fill = "aquamarine3", color = "aquamarine3") + 
    geom_sf(size = 5.5, color = "white") +
    theme_map() + labs(title = "daily min > 24°C % exceedance") +
    geom_sf(aes(color = intmin), size = 4) + theme(legend.position = "none") +
    ylim(min(sites$Lat), max(sites$Lat)) +
    geom_label_repel(data = sites_sf,
      aes(label = Sitelabel, geometry = geometry),
      stat = "sf_coordinates",
      min.segment.length = 0) +
    scale_color_viridis_c(option = "C", limits = c(0,max(c(sites$intmax, sites$intmin)))))

legendary = cowplot::get_legend(iminplot + labs(color = "% exceedance")) + 
                                  theme(legend.position = "bottom",legend.justification="center") 
legendairy = cowplot::get_legend(imaxplot + labs(color = "% exceedance")+
                                  theme(legend.position = "right"))

png("output/DS_HTI_heat_exceedance_map_%02d.png", 
    height = 6, width = 7, unit = "in", res = 1000, family = "serif")

cowplot::plot_grid(
  cowplot::plot_grid(imaxplot, iminplot, nrow = 1), 
  legendary, nrow = 2, rel_heights = c(10,2))


cowplot::plot_grid(imaxplot, iminplot, legendairy, nrow = 1, rel_widths = c(10,10,2.5), rel_heights = c(1,1,5))

dev.off()


