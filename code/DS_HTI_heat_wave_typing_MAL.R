## Delta Smelt heat threat index
## Heat wave typing and recurrence intervals

library(tidyverse)
library(lubridate)
library(ggrepel)
library(vegan)
library(ggrepel)
# library(streamgraph)

##Water year julian day
hydro.day.new = function(x, start.month = 10L){
  start.yr = year(x) - (month(x) < start.month)
  start.date = make_date(start.yr, start.month, 1L)
  as.integer(x - start.date + 1L)
}

# Data download -----------------------------------------------------------

## Pull in Lisbon weir temp data from water data library
## B91560Q - Yolo Bypass near Lisbon Weir Lat:38.4747608 Long:-121.5885961 Elev:1ft
## 363.00 - ADCP Water Temperature (degrees C)
# Qualities:
# 1 - Good data
# 2 - Good quality edited data
# 151 - Data Missing
# 170 - Unreliable Data
# 255 - No data exists

# listemp <- read.csv("C:/Users/eholmes/Documents/R/Projects/Sentinel_CSC/data/B91560Q_Water_Temperature_ADCP_Raw.csv")

load("data/DS_HTI_Wtemp_RVB+MAL.Rdata")

listemp <- wtempdata[wtempdata$Site_no == "MAL" & is.na(wtempdata$Param_val) == F & wtempdata$Param_val < 30,]
listemp$Date <- as.Date(listemp$Datetime)
listemp$Year <- format(listemp$Datetime, format = "%Y")
listemp <- listemp[is.na(listemp$Datetime) == F,]
listemp$jdayfrac <- sapply(listemp$Datetime, function(x) julian(x, origin = as.POSIXct(paste0(format(x, "%Y"),'-01-01'), tz = 'GMT')))

# Data processing ---------------------------------------------------------

#Plot full time series
ggplot(listemp, aes(x = Date, y = Param_val)) + geom_line() + theme_bw()

lisply <- listemp %>% group_by(Date) %>% summarize(mintemp = min(Param_val),
                                                   meantemp = mean(Param_val),
                                                   maxtemp = max(Param_val))

lis <- lisply

#Set heating threshold (Maybe set above 25)
lis$heat <- ifelse(lis$mintemp > 21 | lis$maxtemp > 23, "over", "under")

#Discretize contiguous heat events
lis$runs <- rep(seq_along(rle(lis$heat)$lengths), rle(lis$heat)$lengths)

#Add water year
lis$Year <- as.numeric(format(lis$Date, format = "%Y"))
lis$Month <- as.numeric(format(lis$Date, format = "%m"))
# lis$WY <-  ifelse(lis$Month %in% c(10:12), lis$Year + 1, lis$Year)
lis$jday <- as.integer(format(lis$Date, format = "%j"))

#Subset to Summer/fall period
lis <- lis[lis$Month %in% 6:11,]

#Provide each heat event with a unique ID
lis$heat_ID <- paste(lis$Year, sprintf("%06d", lis$runs), sep = "_")

# Calculate metrics for each heat event ----------------------------------

##Calculate rate of change
lis$ratechange <- c(NA, diff(lis$meantemp))

#Identify rising or descending limbs of hydrograph
lis$limb <- ifelse(lis$ratechange > 0, "rising", "descending")

#Identify heat peaks and troughs
lis <- merge(lis, data.frame(Date = lis[cumsum(rle(lis$limb)$lengths), "Date"], 
                             Peak = ifelse(rle(lis$limb)$values == "descending", "trough", "peak")), 
             by = "Date", all.x = T)

##Add peak location and number of peaks per heat
lispeaks <- lis %>% filter(heat == "over" & Peak == "peak") %>% group_by(heat_ID, Peak) %>% 
  add_tally() %>% group_by(heat_ID) %>% summarise(Peaks = mean(n))

# ##Add trough location and number of troughs per heat wave
# listroughs <- lis %>% filter(heat == "over" & Peak == "peak") %>% group_by(heat_ID, Peak) %>% 
#   add_tally() %>% group_by(heat_ID) %>% summarise(Peaks = mean(n))

#Calculate heat event summary metrics
lissum <- lis %>% filter(heat == "over") %>% group_by(heat_ID) %>%
  summarise(startdate = min(Date), enddate = max(Date), 
            maxtemp = max(maxtemp), meantemp = mean(meantemp), sumtemp = sum(meantemp),
            riserate = max(ratechange), fallrate = min(ratechange))

##Merge heat summary metrics with heat peak metrics
lissum <- merge(lissum, lispeaks, by = "heat_ID", all.x = T)

##Calculate heat centroid date
lissum$centdate <- lissum$startdate + floor((lissum$enddate-lissum$startdate)/2)

##Calculate heat duration
lissum$duration <- as.numeric(lissum$enddate - lissum$startdate) + 1

##Convert dates to julian water days
lissum$startday <- hydro.day.new(lissum$startdate)
lissum$endday <- hydro.day.new(lissum$enddate)
lissum$centday <- hydro.day.new(lissum$enddate)

##Center of heat timing prep
lis$cttop <- lis$jday*lis$meantemp

##Calculate heat shape and timing metrics
heatshape <- data.frame()
for(heat in lissum$heat_ID){
  print(heat)
  heatshape <- rbind(heatshape, data.frame(heat_ID = heat,
                                           cumtemp = sum(lis[lis$Year %in% as.integer(substr(heat, 1, 4)) & 
                                                               lis$jday < lissum[lissum$heat_ID == heat, "startday"], "meantemp"]),
                                           peakday = lis[lis$heat_ID %in% heat & lis$maxtemp == lissum[lissum$heat_ID %in% heat, "maxtemp"], "jday"][1],
                                           heatct = sum(lis[lis$heat_ID %in% heat, "cttop"])/sum(lis[lis$heat_ID %in% heat, "maxtemp"])))
}

##Merge summary metrics and heat shape metrics
lissum <- merge(lissum, heatshape, by = "heat_ID", all.x = T)

##Calculate peak location
lissum$peakloc <- (lissum$peakday - lissum$startday)/(lissum$endday - lissum$startday + .000001)

##Calculate centroid volume location
lissum$centvolloc <- (lissum$heatct - lissum$startday)/(lissum$endday - lissum$startday + .000001)

# K-means clustering ------------------------------------------------------

row.names(lissum) <- lissum$heat_ID

##Set minimum heat duration threshold
minduration <- 2

##Set metrics to include in the scaling, K-means and PCA analyses
# "endday","startday", "sumflow", "peakloc", "centvolloc", "meanflow", "cumflow" ##unused metrics , "fallrate"
metrics <- c("duration", "centday", "maxtemp")

##Scale heat summary data
lisscale <- (scale(lissum[lissum$duration > minduration, metrics]))

set.seed(25842)                             
##K-means clustering analysis
lismeans <- kmeans(lisscale, centers = 3)

##append K-means clusters to scaled heat metrics
lispairs <- data.frame(cbind(lissum[lissum$duration > minduration, c("heat_ID")], lisscale, lismeans[["cluster"]]))

##Convert clusters to factors
lispairs$clust <- as.factor(lismeans[["cluster"]])

lispairs$duration <- as.numeric(lispairs$duration)
# lispairs$mintemp <- as.numeric(lispairs$mintemp)
lispairs$centday <- as.numeric(lispairs$centday)
lispairs$maxtemp <- as.numeric(lispairs$maxtemp)
# lispairs$peakday <- as.numeric(lispairs$peakday)

##Correlation plot
pairs(lispairs[,metrics], col = lismeans[["cluster"]])
lispairs <- lispairs[lispairs$V1 != "2021_000781",]

ggplot(lispairs, aes(x = maxtemp, y = duration, color = clust)) + geom_point() + theme_bw()
ggplot(lispairs, aes(x = maxtemp, y = centday, color = clust)) + geom_point() + theme_bw()
ggplot(lispairs, aes(x = duration, y = centday, color = clust)) + geom_point() + theme_bw()
ggplot(lispairs, aes(x = centday, y = maxtemp, color = clust)) + geom_point() + theme_bw()

# PCA ---------------------------------------------------------------------

#Compute principle component analysis
lispca <- prcomp(lisscale)

#stock PCA plots
plot(lispca)
biplot(lispca)

#Prepare PCA data for ggplot
lispcadf <- data.frame(lispca[["x"]])
lispcadf$heat_ID <- row.names(lispcadf)
lispcadf$Cluster <- factor(lismeans$cluster)
lispcadf$Year <- substr(lispcadf$heat_ID, 1, 4)

#Prepare PCA rotations for ggplot
lisrot <- data.frame(lispca$rotation)
lisrot$var <- row.names(lisrot)

centroids <- lispcadf %>% group_by(Cluster) %>% summarise(PC1mean = mean(PC1), PC2mean = mean(PC2))

centroids$PC1offset <- c(0.4, .5, .15)
centroids$PC2offset <- c(0, .6, -.4)
clustername <- data.frame(Cluster = 1:3, 
                          label = c("Intense", "Late-short", "Early-short"),
                          
                          Key = c("Intense", "Late-short", "Early-short"))

clustername$Key <- factor(clustername$Key, levels = c("Intense", "Late-short", "Early-short"))

centroids <- merge(centroids, clustername, by = "Cluster", all.x = T)
lispcadf <- merge(lispcadf, clustername, by = "Cluster", all.x = T)

##Convex hull computation for each cluster
hull_clust <- lispcadf %>% group_by(Key) %>% slice(chull(PC1, PC2))
##Summary pre and post dam
##Calculate annual frequency

lispcadfannual <- lispcadf %>% group_by(Key, label) %>% summarise(count = length(PC1))

lispcadfannual$label <- factor(lispcadfannual$label, levels = c("Early\nsmall", "Late\nsmall", "Intermediate", "Long\nduration", "Ravaging"))
lispcadfannual$totyears = length(unique(lis[lis$Year < 8,"Year"]))

##Convert frequency to observed recurrence frequency
lispcadfannual$Recurrence <- (lispcadfannual$count)/lispcadfannual$totyears
lispcadfannual$Recurrence_int <- (lispcadfannual$totyears +1)/lispcadfannual$count

lisclust <- merge(lis, lispcadf[, c("heat_ID", "Cluster", "Key")], by = "heat_ID", all.x = T)

medheat <- lisclust %>% group_by(Key, jday) %>% 
  summarise(mediantemp = median(meantemp), meantemp = mean(meantemp, na.rm = TRUE),
            sdtemp = sd(meantemp, na.rm = TRUE),
            ntemp = n()) %>%
  mutate(setemp = sdtemp / sqrt(ntemp),
         lower.citemp = meantemp - qt(1 - (0.05 / 2), ntemp - 1) * setemp,
         upper.citemp = meantemp + qt(1 - (0.05 / 2), ntemp - 1) * setemp)

heatseas <- lisclust %>% filter(is.na(Cluster) == F) %>% group_by(Year) %>% 
  summarise(heatstart = min(jday), heatend= max(jday))

# if(saveplots == T){png("C:/Users/ejholmes/Box/Holmes/CV_heattyping/Output/Verona/Verona_heat_typing_Kmeans%03d.png", 
#     family = "serif", width = 6.5, height= 5, units = "in", res = 500)}

##Plot PCA with K-means convex hulls
kmeanspca <- ggplot(lispcadf, aes(x = PC1, y = PC2, color = Key)) + 
  geom_point(alpha = .25, show.legend = F) +
  geom_polygon(data = hull_clust, aes(fill = Key), alpha = .07, show.legend = F) + 
  geom_text(data = lisrot, aes(x = PC1*2, y = PC2*2, label = var),color = "black") + 
  geom_segment(data = lisrot, aes(x = 0, y = 0, xend = PC1*2, yend = PC2*2), color = "black") +
  scale_fill_brewer(palette = "Set1") + scale_color_brewer(palette = "Set1") +
  scale_shape_manual(values = c("Pre" = 1, "Post" = 3)) + 
  theme_bw() +
  geom_text(data = centroids, aes(x = PC1mean + PC1offset, y = PC2mean + PC2offset, label = label), 
            fontface = "bold", show.legend = F)

kmeanspca

##Plot PCA with representative water years
kmeanspca + geom_point(data = lispcadf[lispcadf$Year %in% c("2017"),], shape = "7", color = "black", size = 3) + 
  geom_point(data = lispcadf[lispcadf$Year %in% c("2019"),], shape = "9", color = "black", size = 3) +
  geom_point(data = lispcadf[lispcadf$Year %in% c("2018"),], shape = "8", color = "black", size = 3) +
  geom_point(data = lispcadf[lispcadf$Year %in% c("2024"),], shape = "4", color = "red", size = 3) 

##Plot heat recurrence timing
# recur_rate <- ggplot(lispcadfannual, aes(x = Year)) + geom_bar(aes(y = Recurrence), stat = "identity") + 
#   facet_grid(. ~ Key) + theme_bw() + labs(y = "Recurrence rate")
# recur_rate

recur_int <- ggplot(lispcadfannual, aes(x = Key)) + geom_bar(aes(y = Recurrence_int), stat = "identity") + 
  theme_bw() + labs(y = "Recurrence interval", x = "Shasta dam construction")
recur_int

recur_inthoriz <- ggplot(lispcadfannual, aes(x = Key)) + 
  geom_bar(aes(y = Recurrence_int, fill = Key), stat = "identity", show.legend = F) + 
  theme_bw() + labs(y = "Recurrence interval", x = "Shasta dam construction") + coord_flip() +
  scale_fill_brewer(palette = "Set1")
recur_inthoriz

##Plot heats by cluster with median day temp
(heattypes <- ggplot() + geom_line(data = lisclust[is.na(lisclust$Cluster) == F,], 
                                   aes(x = jday, y = meantemp, group = heat_ID, color = Key), alpha = .2, show.legend = F) +
    geom_ribbon(data =lisclust[is.na(lisclust$Cluster) == F,], 
                aes(x = jday, ymin = 22, ymax = maxtemp, group = heat_ID, fill = Key), alpha = .1, show.legend = F) + 
    geom_line(data = medheat[is.na(medheat$Key) == F,], aes(x = jday, y = meantemp, color = Key), show.legend = F) + 
    facet_grid(Key ~ .) + theme_bw() + labs(x = "Day of year", "temp (C)") +
    scale_color_brewer(palette = "Set1") + scale_fill_brewer(palette = "Set1"))

heattypes

heat_timeline <- ggplot() + 
  geom_line(data = lisclust[is.na(lisclust$Cluster) == F,], show.legend = F,
            aes(x = jday, y = Year, group = heat_ID, color = Key), linewidth = 5) +
  scale_color_brewer(palette = "Set1") + theme_bw() + labs(y = NULL, x = "Day of water year") +
  scale_y_continuous(breaks = seq(1980, 2024, 1), limits = c(1987,2024))

heat_timeline 

gg1 <- ggplot_build(ggplot() + stat_smooth(data = heatseas, aes(x = Year, y = heatstart), span = .15, size = .5,se = F) +
                      stat_smooth(data = heatseas, aes(x = Year, y = heatend), span = .15, size = .5,se = F))

df2 <- data.frame(x = gg1$data[[1]]$x,
                  ymin = gg1$data[[1]]$y,
                  ymax = gg1$data[[2]]$y) 

if(saveplots == T){png("output/DS_HTI_heat_wave_timeline_MAL_%02d.png", 
                       family = "serif", width = 6.5, height= 7.5, units = "in", res = 500)}


left <- cowplot::plot_grid(kmeanspca, heattypes, labels = c("A", "B"),nrow = 2, greedy = T)
leftrecur <- cowplot::plot_grid(kmeanspca, recur_inthoriz, labels = c("A", "B"), 
                                nrow = 2, rel_heights = c(5,6.5), greedy = T)

cowplot::plot_grid(left, heat_timeline, labels = c("", "C"), ncol = 2)

# cowplot::plot_grid(leftrecur, heat_timeline_flipped, labels = c("", "C"), ncol = 2)

dev.off()

# Water year typing -----------------------------------
# Center of water year timing (DONE) correlation with heat types defined above
# Summarize water year by heat types

##Center of water year timing
lis$cttop <- lis$jday*lis$meantemp

t <- aggregate(lis$cttop, by = list(lis$Year), FUN = "sum")
b <- aggregate(lis$meantemp, by = list(lis$Year), FUN = "sum")
CT <- t/b
CT$yr <- unique(lis$Year)

plot(x~yr, data = CT[CT$yr %in% 2014:2022,], ylab="Days Since Start of Water Year", xlab="Year")
# abline(lm(x ~ yr, data = CT[CT$yr %in% 2014:2022,]))
# abline(lm(x ~ yr, data = CT[CT$yr>=8,]))
abline(lm(CT$x~CT$yr))
lines(loess.smooth(CT$yr,CT$x, span = .55), col="red", lty=2)
grid()
legend("topleft", c("Linear trend", "Loess Line"),lty=c(1,2),col=c("black","red"),bty="n", cex=.9)

if(saveplots == T){png("C:/Users/ejholmes/Box/Holmes/CV_heattyping/Output/Verona/Verona_heat_typing_CT_trends%03d.png", 
                       family = "serif", width = 6.5, height= 4.5, units = "in", res = 500, Param_valsize = 7.5)}

(CTdam <- ggplot(CT[CT$yr < 2024,], aes(y = x, x = yr)) + geom_point() +
    labs(y = "Water year center of mass", x = "year") +
    stat_smooth(data = CT[CT$yr < 2024,], method = "lm", se = F) + theme_bw())

dev.off()

##Tally heattypes per year
Yearlong <- lispcadf %>% group_by(Year, Key) %>% tally() 
Year <- Yearlong%>% spread(Key, n, fill = 0)
unique(lispcadf$Key)

Yearclass <- data.frame()
for(wyear in as.character(2014:2022)){
  if(wyear %in% Year$Year){print(wyear); Yearclass <- rbind(Yearclass, Year[Year$Year == wyear,])}
  else{print("NOT in dataset")
    Yearclass <- rbind(Yearclass, data.frame("Year" = wyear, "Earlysmall" = 0, "Intermediate" = 0, 
                                         "Latesmall" = 0, "Longduration" = 0, "Ravaging" = 0))}
}

Yearclass$events <- rowSums(Yearclass[,c(2:4)])

Yearnoheat <- Yearclass[Yearclass$events == 0,]
Yearnoheat$Yearnum <- as.numeric(Yearnoheat$Year)

ggplot(Yearnoheat, aes(x = Yearnum)) + geom_dotplot(binwidth = 1) +
  theme_bw()

Yearlong$Yearnum <- as.numeric(Yearlong$Year)
Yearlong$cuts <- cut(Yearlong$Yearnum, breaks = seq(1890, 2020, 10))

Yearclasslong <- reshape2::melt(Yearclass)
Yearclasslong$Year <- as.integer(Yearclasslong$Year)
# Yearlong %>% group_by(Key)
# 
# unique(Yearlong$Key)


2/60; 15/77

if(saveplots == T){png("C:/Users/ejholmes/Box/Holmes/CV_heattyping/Output/Verona/Verona_heat_typing_violins%03d.png", 
                       family = "serif", width = 6.5, height= 4.5, units = "in", res = 500, Param_valsize = 7.5)}
(pyplot <- ggplot(Yearlong, aes(x = Key, y = Yearnum)) + geom_violin(aes(fill = Key), alpha = .5, bw = 5) + 
    geom_dotplot(aes(y = Yearnum), binaxis = "y", binwidth = 5, stackdir = "center", dotsize = .5) + 
    theme_minimal() + #scale_y_continuous(breaks = seq(1890, 2020, 10)) + scale_fill_brewer(palette = "Set1") +
    theme(legend.position = "none")) + labs(x = NULL, y = NULL)

(hyplot <- ggplot(Yearlong, aes(x = Yearnum)) + geom_histogram(aes(fill = Key), alpha = .5) + 
    theme_minimal() + scale_fill_brewer(palette = "Set1") +
    theme(legend.position = "none")) + labs(x = NULL, y = NULL) + facet_grid(Key ~ .)

(wyplot <- ggplot(Yearclasslong[Yearclasslong$variable != "events",], aes(x = Year, y = value)) + geom_point(aes(color = variable), alpha = .5) + 
    stat_smooth(aes(color = variable), se = F) +
    theme_minimal() + scale_color_brewer(palette = "Set1") +
    theme(legend.position = "none") + labs(x = NULL, y = NULL))

dev.off()

# heat frequency analysis ------------------------------------------------
# flow duration curve pre and post dam by water year type DONE
# Median flow hydrographs pre and post with 95% confidence intervals Done
listemp$jday <- as.numeric(format(listemp$Datetime, format = "%j"))
listemp$Date <- as.Date(listemp$Datetime)
listemp <- listemp %>% group_by(Date) %>% mutate(mintemp = min(Param_val),
                                      meantemp = mean(Param_val),
                                      maxtemp = max(Param_val))

medtemp <- listemp %>% group_by(jday) %>%
  summarise(mediantemp = median(meantemp), meantemp = mean(meantemp, na.rm = TRUE),
            sdtemp = sd(meantemp, na.rm = TRUE),
            ntemp = n()) %>%
  mutate(setemp = sdtemp / sqrt(ntemp),
         lower.ciflow = meantemp - qt(1 - (0.05 / 2), ntemp - 1) * setemp,
         upper.ciflow = meantemp + qt(1 - (0.05 / 2), ntemp - 1) * setemp)
# 
# lis$period <- cut(lis$Year, seq(1890, 2030, 10), labels = paste(seq(1890, 2020, 10), "s", sep = ""))
# 
# periodflow <- lis %>% filter(Year > 1930 & Year < 2021 & jday < 365) %>% group_by(period, jday) %>% 
#   summarise(medianflow = median(Flow), meanflow = mean(Flow, na.rm = TRUE),
#             sdflow = sd(Flow, na.rm = TRUE),
#             nflow = n()) %>%
#   mutate(seflow = sdflow / sqrt(nflow),
#          lower.ciflow = meanflow - qt(1 - (0.05 / 2), nflow - 1) * seflow,
#          upper.ciflow = meanflow + qt(1 - (0.05 / 2), nflow - 1) * seflow)
# 
# meanflow <- data.frame(lis %>% filter(Month %in% c(12,1:4)) %>% group_by(dam, Year, period) %>% 
#                          summarise(meanflow = mean(Flow, na.rm = T), medianflow = median(Flow, na.rm = T)))

summerflow <- data.frame(lis) %>% group_by(Year) %>% 
                           summarise(meantemp = mean(meantemp, na.rm = T), mediantemp = median(meantemp, na.rm = T))
# meantemp <- lis %>% filter(Month %in% c(12,1:4)) %>% group_by(dam, Year, period) %>% summarise(meantemp = sum(temp, na.rm = T))
# wytype <- read.csv("C:/Users/ejholmes/Box/Holmes/CV_heattyping/Data/CDEC_Water_year_class.csv")
# wytype$dacr <- wytype$Dec + wytype$Jan + wytype$Feb + wytype$Mar + wytype$Apr
# 
# meantemp <- merge(meantemp, wytype, by = "Year", all.x = T)
# meantemp$Type <- ifelse(meantemp$Year %in% 2021, "C", meantemp$Type)
# meantemp$Typefac <- factor(meantemp$Type, levels = c("W", "AN", "BN",   "D",  "C"))
# levels(meantemp$Typefac) <- c("Wet", "Above normal", "Below normal", "Dry", "Critical")

png("output/DS_HTI_heat_wave_types_MAL_%02d.png", 
    family = "serif", width = 6.5, height= 7.5, units = "in", res = 1000, pointsize = 7.5)

listemp2 <- merge(listemp, lisclust[, c("Date", "Cluster", "Key", "heat_ID")], by = "Date", all.x = T)

ggplot(lis[lis$Year == 2019,], aes(x = jday, y = meantemp)) + theme_bw() +
  geom_line(data = listemp[listemp$Year == 2019,], aes(x = jdayfrac, y = Param_val)) +
  geom_ribbon(data = lisclust[lisclust$Year == 2019 & is.na(lisclust$Cluster) != T,], alpha = .5, 
              aes(ymin = 18, ymax = meantemp, group = heat_ID, fill = Key)) + scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "bottom")  + labs(x = NULL) + ylim(18,26)+
  geom_hline(yintercept = c(24,27), linetype = 2) +
  scale_x_continuous(limits = c(150, 270),
                     breaks = c(0, 31, 59, 90, 120, 151, 181, 212, 242, 272, 303, 333, 364),
                     labels = c(month.abb, month.abb[1]))

ggplot(lis[lis$Year == c(2019, 2021, 2023, 2024),], aes(x = jday, y = meantemp)) + theme_bw() +
  geom_line(data = listemp[listemp$Year == c(2019, 2021, 2023, 2024),], aes(x = jdayfrac, y = Param_val)) +
  geom_ribbon(data = listemp2[listemp2$Year == c(2019, 2021, 2023, 2024),], alpha = .5, 
              aes(ymin = 18, ymax = meantemp, group = heat_ID, fill = Key)) + scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "bottom")  + labs(x = NULL) + ylim(18,26)+
  geom_hline(yintercept = c(22,24), linetype = 2) +
  scale_x_continuous(limits = c(150, 270),
                     breaks = c(0, 31, 59, 90, 120, 151, 181, 212, 242, 272, 303, 333, 364),
                     labels = c(month.abb, month.abb[1]))+ facet_grid(Year ~ .)

ggplot(lis[lis$Year == c(2005:2024),], aes(x = jday, y = meantemp)) + theme_bw() +
  geom_line(data = listemp[listemp$Year == c(2005:2024),], aes(x = jdayfrac, y = Param_val)) +
  geom_ribbon(data = listemp2[listemp2$Year == c(2005:2024),], alpha = .5, 
              aes(ymin = 18, ymax = meantemp, group = heat_ID, fill = Key)) + scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "bottom")  + labs(x = NULL) + ylim(18,26)+
  geom_hline(yintercept = c(22,24), linetype = 2) +
  scale_x_continuous(limits = c(150, 270),
                     breaks = c(0, 31, 59, 90, 120, 151, 181, 212, 242, 272, 303, 333, 364),
                     labels = c(month.abb, month.abb[1]))+ facet_grid(Year ~ .)

ggplot(lis, aes(x = jday, y = meantemp)) + theme_bw() +
  geom_line(data = listemp, aes(x = jdayfrac, y = Param_val)) +
  geom_ribbon(data = listemp2, alpha = .5, 
              aes(ymin = 18, ymax = meantemp, group = heat_ID, fill = Key)) + scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "bottom")  + labs(x = NULL) + ylim(18,26)+
  geom_hline(yintercept = c(24,27), linetype = 2, color = "grey30") +
  scale_x_continuous(limits = c(150, 270),
                     breaks = c(0, 31, 59, 90, 120, 151, 181, 212, 242, 272, 303, 333, 364),
                     labels = c(month.abb, month.abb[1])) + facet_grid(Year ~ .)
dev.off()

