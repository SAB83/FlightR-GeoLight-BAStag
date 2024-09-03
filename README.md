# FlightR-GeoLight-BAStag
#To analyze purple martin migration and its stopover

library(maps)
library(GeoLight)
library(raster)
library(ks)
library(RCurl) #this package needed for getURL function

sourceLigURL <- getURL('https://raw.githubusercontent.com/SCBI-MigBirds/MigBirds/master/source/sourceRead_lig.R')  # function written by Simon Wotherspoon accessed from GitHUB
sourceLightBugURL <- getURL('https://raw.githubusercontent.com/SCBI-MigBirds/MigBirds/master/source/sourceRead_LightBug.R')  # read_lig function modified by M.T.Hallworth to read in LightBug data

#Code to load function that reads .lig file
#From Github https://github.com/SWotherspoon/BAStag/blob/master/R/BAStag.R
readLig <- function(file,skip=0) {
  ## Read csv file and add column names
  d <- read.csv(file,header=FALSE,skip=skip,
                col.names=c("Valid","Date","Julian","Light"),
                colClasses=c("character","character","numeric","integer"))
  ## Parse date
  d$Date <- as.POSIXct(strptime(d$Date,"%d/%m/%y %H:%M:%S",tz="GMT"))
  d
}

#or
#Import .lig file
devtools::install_github("SWotherspoon/BAStag")

d.lux <- readLig(file.choose())
head(d.lux)
attach(d.lux)

GL_008 <- d.lux

#--------------------------------------------------------------------------
#Determining light levels
#LightThreshold - setting threshold level (threshold 32 for martins in BasTrack)
#ask=TRUE - can go through every twilight
GL_008_transitions <- twilightCalc(datetime = GL_008[,2],
                                   light= GL_008[,4],
                                   LightThreshold=32,
                                   ask=TRUE)

head(GL_008_transitions)

#--------------------------------------------------------------------------
#Using a Loess Filter to remove outliers in defined twilight times based on 
#smoother function
loess_GL_008<-loessFilter(tFirst=GL_008_transitions[,1], 
                            tSecond=GL_008_transitions[,2], type=2, 
                            twl=GL_008_transitions, k = 2, plot = TRUE)
head(loess_GL_008)

#Loess filter is the automatic equlivanent in BasTrack of us manually removing 
#outliers
#Used this when combining time and location to delete outlier values

#--------------------------------------------------------------------------
#Sun elevation angle
SunElev<-getElevation(tFirst= GL_008_transitions[33:63,1],     #ED:last week of July and first week of August
                      tSecond= GL_008_transitions[33:63,2],
                      type=GL_008_transitions[33:63,3],      
                      known.coord=c(-97.130305,49.734761),     #ED:Known coordinates at the breeding grounds, changed for MB Clifton's for GEO 008
                      plot=TRUE)

GL_008Locations<-coord(tFirst= GL_008_transitions[,1],
                       tSecond= GL_008_transitions[,2],
                       type=GL_008_transitions[,3], 
                       degElevation=SunElev, tol=0.13)  #ED:changed to 0.13 from 0                         
#Setting tolerance (positions discarded around equinox) helps remove some outliers
#Tol puts NAs around equinox lats and longs, not 30 days around equinox though
#ED:Investigate what tol=0 is doing, is this removing no positions around equinox?
head(GL_008Locations)  

#-------------------------------------------------------------------------------
#Plot the location data
plot(GL_008Locations, 
     pch="*", 
     col="red",
     xlab="Longitude",
     ylab="Latitude")
map("world",add=TRUE)
#Had to alter from original code, it was not working

#-------------------------------------------------------------------------------
#Determine residency with changelight
stop<-changeLight(tFirst=GL_008_transitions[,1], tSecond=GL_008_transitions[,2], 
                  type=2, twl=GL_008_transitions, quantile = 0.9, rise.prob = NA,
                  set.prob = NA, days = 2, plot = TRUE, summary = TRUE)

#ED:This is the output:::
Probability threshold(s):
  
  Sunrise:  0.04165	Sunset:  0.07321


Migration schedule table:
  
  Site             Arrival           Departure  Days     P.start      P.end Days.1   P.start.1
1     a 2016-07-14 06:27:08 2016-07-29 06:34:12  15.0 0.000000000 0.04496828   15.0 0.000000000
2     b 2016-07-30 06:30:05 2016-08-11 06:42:16  12.0 0.000000000 0.00000000   12.0 0.000000000
3     c 2016-08-12 06:48:16 2016-08-15 18:45:19   3.5 0.000000000 0.00000000    3.5 0.000000000
4     d 2016-08-17 18:41:30 2016-08-20 18:42:27   3.0 0.023865878 0.00500000    3.0 0.023865878
5     e 2016-08-21 18:40:39 2016-08-25 06:39:46   3.5 0.000000000 0.00000000    3.5 0.000000000
6     f 2016-08-26 06:44:01 2016-08-30 18:23:16   4.5 0.004166667 0.02217416    4.5 0.004166667
7     g 2016-08-31 18:22:08 2016-09-03 18:22:04   3.0 0.000000000 0.01574074    3.0 0.000000000
8     h 2016-09-04 18:06:54 2016-09-11 05:34:47   6.5 0.000000000 0.00000000    6.5 0.000000000
9     i 2016-09-15 04:47:35 2016-09-19 04:26:30   4.0 0.000000000 0.00000000    4.0 0.000000000
10    j 2016-09-20 04:26:51 2017-04-11 04:02:01 203.0 0.003703704 0.00000000  203.0 0.003703704
11    k 2017-04-17 05:39:36 2017-04-21 17:58:16   4.5 0.000000000 0.00000000    4.5 0.000000000
12    l 2017-04-22 18:10:49 2017-04-28 06:18:05   5.5 0.037500000 0.00000000    5.5 0.037500000
13    m 2017-05-02 06:21:25 2017-05-05 18:24:34   3.5 0.000000000 0.00000000    3.5 0.000000000
14    n 2017-05-11 03:19:53                <NA>    NA 0.000000000 0.00000000     NA 0.000000000


#5 days is probably too long. Not able to pick up spring migration stopovers
#2 days seems good
#Changelight output is what I need to look at for residency periods
#--------------------------------------------------------------------------
  
#dataframe summarizing residency and movement
schedule(tFirst=GL_008_transitions[,1], tSecond=GL_008_transitions[,2], stop$site)

#schedule is weird for residency periods - they overlap with each other
#Use chart from changeLigh instead
#--------------------------------------------------------------------------
#Draw positions and trip on map
tripMap(GL_008Locations, equinox = TRUE, map.range = c("America"), legend = TRUE)

#--------------------------------------------------------------------------
#Draw sites of residency - sort of gives me what I want with locations combined with dates
siteMap(GL_008Locations, stop$site, type = "points", quantiles = c(0.25, 0.75), hull = T,
        map.range = c("America"))

#--------------------------------------------------------------------------
#Convert to dataframes

#loess filter
frameloess<-data.frame(loess_GL_008)
write.csv(frameloess, "frameloess.csv")

#location file
framelocation<-data.frame(GL_008Locations)
write.csv(framelocation, "framelocation.csv")

#timing file
frametime<-data.frame(GL_008_transitions)
write.csv(frametime, "frametime.csv")

#--------------------------------------------------------------------------

#How to pull out the lats/longs associated with residency periods?
#No way to do with automatically in GeoLight
#Use the three .csv files and combine by hand. Then highlight the residency
#periods given ny the changeLight table. Use date and location to determine
#the different timing periods

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#SECOND HALF: *******FLightR**********

#trying to input things into flightR. Using FLightR package vingettes and also using info from "FLIGHTR: an R package for reconstructing animal paths from solar geolocation loggers" by Rakhimberfiev et al, 2017, Methods in Ecology and Evolution

#From GeoLight: just using lig file and transitions defined from twilightCalc function

library(FLightR)
#taking only "Date" and "Light" column from raw GL_008 file
GL_008b <- subset(GL_008, select=c("Date", "Light"))

#changing data into TAGS format 
TAGS_twilights <- GeoLight2TAGS(GL_008b, GL_008_transitions, threshold=32)

#removing NA
TAGS_twilights <- na.omit(TAGS_twilights)

#changing data into .csv format
write.csv(TAGS_twilights, "TAGS_twilights")

#defining start and end date for migration period (can change later, testing for now), need to use start and end for get.tags.data function
start_date <- as.POSIXct("2016-07-01 00:01:00", "%Y-%m-%d %H:%M:%S", tz = "GMT")
end_date <- as.POSIXct("2017-05-01 00:01:57", "%Y-%m-%d %H:%M:%S", tz = "GMT")

#reads csv file in TAGS format and returns it as a FLIGHTR data object
Proc.data <- get.tags.data("TAGS_twilights", start_date, end_date, log.light.borders = "auto", log.irrad.borders = "auto", saves = c("auto", "mean", "max"), measurement.period = 120, impute.on.boundaries = FALSE)

#calibration period (last week of July & first week of August)
calibration.start<- as.POSIXct("2016-07-25 00:00:00", "%Y-%m-%d %H:%M:%S", tz = "GMT")
calibration.stop<- as.POSIXct("2016-08-07 00:00:00", "%Y-%m-%d %H:%M:%S", tz = "GMT")

Calibration.periods <- data.frame(calibration.start, calibration.stop, lon=-97.130305,lat=49.734761)
print(Calibration_periods)

MB_calib <- make.calibration(Proc.data, Calibration.periods, model.ageing = FALSE, plot.each = FALSE, plot.final = FALSE, likelihood.correction = "auto", fixed.logSlope = c(NA, NA), suggest.irrad.borders = FALSE, return.slopes = FALSE)

Grid <-make.grid(left = -110, bottom = -40, right = -50, top = 80,
                        distance.from.land.allowed.to.use = c(-Inf, Inf),
                        distance.from.land.allowed.to.stay = c(-Inf, Inf), plot = TRUE,
                        return.distances = FALSE, probability.of.staying = 0.5)

#combining data, calibration, spatial extent and movement priors and estimates spatial likelihoods
Combined_things <- make.prerun.object(Proc.data, Grid, start=c(-97.130305,49.734761), end = start, MB_calib, threads = -1, Decision = 0.05, Direction = 0, Kappa = 0, M.mean = 300, M.sd = 500, likelihood.correction = TRUE)

library(BAStag)
library(GeoLight)
library(geosphere)
library(FLightR)
library(ggplot2)
library(ggmap)
library(maptools)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Converting Existing GL data to Tags format##
## Pre-Processing light (Twilight)##
tags.raw <- readLig(file.choose())
tags.raw <- subset(tags.raw, select=c("Date","Light")) 
tags.raw=tags.raw[-1,]


lightImage(tags.raw,offset=19, zlim = c(0, 12))

twl= preprocessLight(tags.raw, offset= offset, threshold=32, zlim = c(0, 12))  ## QQ to exit ##


twl <- preprocessLight(tagdata=tags.raw, 
                       threshold=32, 
                       offset=19, 
                       zlim=c(0,12), 
                       dark.min=240)
threshold=32

TAGS.twilight <- BAStag2TAGS(tags.raw, twl,threshold) ## this converts the lig to BASTAGs##

TAGS.twilight$datetime<-format(TAGS.twilight$datetime, 
                               format="%Y-%m-%dT%T.000Z")


### Saving the GL data as a TAGS format file for use in next section##

write.csv(TAGS.twilight, file="3652VA2012re.csv", quote=FALSE, row.names=FALSE) ## Remember to change .csv file name each run##

#-------------------------------------------------------------------------------

# What if you are not sure what data are good for calibrating?
# Here's a process that uses FlightR to delineate good calibration periods

Proc.data<-get.tags.data(file.choose()) # opens and formats data straight from TAGS 
# formatted csv file
# indicates kind of tag and mode set in

Lat.calib<- (38.61301)
Lon.calib <- (-77.262662) ### from word file or excel###


start <- c(Lon.calib, Lat.calib)  # tracking orgin, start location 

## finding claibration period from twilight data ## ## aim for ~3-4 weeks but use what is available from plot##
plot_slopes_by_location(Proc.data=Proc.data, location=start)

# Black and red points should be together during the appropriate calibrating 
# period

# Using locator to determine the start of calibrating period
start.calib <- as.POSIXct(locator(n=1)$x, origin = "1970-01-01")
start.calib
# Use locator function to determine the end of the appropriate calibrating 
# period
end.calib <- as.POSIXct(locator(n=1)$x, origin = "1970-01-01")
end.calib

#or define period
start.calib <- as.POSIXct("2012-07-15 18:23:24 CDT")
end.calib <- as.POSIXct("2012-07-30 02:04:29 CDT")


# Next create a data.frame with a separate line is designated for each 
# calibration period. 
# The columns are: 
#     (1) start of the calibration period, 
#     (2) end of the calibration period, 
#     (3) longitude of the calibration location and, 
#     (4) latitude of the calibration location.

Calibration.periods<-data.frame(   #This will create two lines of data
  calibration.start=as.POSIXct(start.calib),   
  calibration.stop=as.POSIXct(end.calib),
  lon=Lon.calib, lat=Lat.calib) # use c() also for the coordinates, 
# if you have more than one calibration location 
# (e. g.,  lon=c(5.43, 6.00), lat=c(52.93,52.94))



print(Calibration.periods) # viewing calibration periods
PMCalibration<- make.calibration(Proc.data, Calibration.periods)

#-------------------------------------------------------------------------------
### Establishing a spatial grid and rules for migration paths ###

# Default grid resolution is is 50 X 50 km.
# Terms "left," "right," "bottom," and "top" define your bounding box. 
# 
# distance.from.land.allowed.to.use should be vector with length of two, 
#   first number is negative distance allowed to use while over land 
#   (restricts birds to flying only over coastlines and water) 
#   and second is distance from land allowed to use while over water 
#   (restricts birds to flying only over coastlines and land). 
#   
# distance.from.land.allowed.to.stay should be vector with length of two, 
#   first number is negative distance where bird is allowed to be stationary, 
#   (restricts birds to landing only on coastlines and land)
#   and second is distance from land allowed to fly over during twilight while 
#   over water. (restricts birds to landing only on coastlines and water)


### 50 km grid, Inf means infinity, if antelope -Inf, 0

Grid<- make.grid(left=-130, bottom=-40, right=-34, top=70,
                 distance.from.land.allowed.to.use=c(-Inf, Inf),  
                 distance.from.land.allowed.to.stay=c(-Inf, Inf))

#-------------------------------------------------------------------------------
# Create a proposal
# Here we create an array of settings and data that incorporates all the objects
# created at earlier steps: 
#    - the light data with the detected twilight events (Proc.data), 
#    - the spatial parameters (Grid), 
#    - geographic coordinates of the starting location (start)
#    - and the calibration parameters (Calibration).
# This can take a while.

a= Sys.time()

all.inGEO<-make.prerun.object(Proc.data, Grid, start=c(Lon.calib, Lat.calib),end = NA, Calibration=PMCalibration)


b= Sys.time()
b-a

#-------------------------------------------------------------------------------
# Run the particle filter

# Here is where the results are are calculated 
# (coordinates, behavior, stationarity).
# Within the function run.particle.filter, the following parameters can be 
# preset:
#   -number of particles (1e4 is recommended for test and 1e6 for the analysis) 
#   -known.last = TRUE if you know the track ends where it began  (FALSE is the 
#    default) 
#   -check.outliers = TRUE, for the "on a fly" discard of outliers (only 
#    recommended to make pretty maps).

#### NOTE: First include outliers, then you can exclude outliers if you decide

nParticles=1e6  # Change to a thousand (1e4) for trial (to get an initial idea) 
# and 1 million (1e6) for real analysis
a= Sys.time()   # Measure the analysis time

Result3652VA2012re<-run.particle.filter(all.inGEO, threads=-1,
                                        nParticles=nParticles, known.last=FALSE,
                                        precision.sd=25, check.outliers=T, b=1500)


b= Sys.time()
b-a                 # How long did it take? 19 min for 1e4 

# Now save your results are as an RData object.

save(Result3652VA2012re, file="Result3652VA2012re.RData")

## Mapping Results##

#Get the google key

ggmap:: register_google('???')

?register_google

map.FLightR.ggmap(Result3652VA2012re, seasonal.donut.location = 'bottomleft' )+  
  
  ?map.FLightR.ggmap

#____________________________________________________________________________________
###Stopover Summary######
stationary.migration.summary(Result3652VA2012, prob.cutoff= 0.1, min.stay= 1)


#-------------------------------------------------------------------------------
#map by leaflet

library(leaflet)
#Test widget
m = leaflet()
m = addTiles(m)
m = m %>% addProviderTiles(providers$Esri.WorldImagery)

fall_path = data.frame(DayDat8898[c(1:181),c(2,3)]) #from fall departure to fall arrival
spring_path = data.frame(DayDat8898[c(540:597),c(2,3)]) #from spring departure to spring arrival
#winter = data.frame(DayDat5520[c(184:539),c(2,3)])

falldep = DayDat8898[42,c(2,3)]
fallarr = DayDat8898[181,c(2,3)]
sprdep = DayDat8898[540,c(2,3)]
sprarr = DayDat8898[1,c(2,3)]

map8898 = addCircleMarkers(map_leaflet, lng = Lon, lat = Lat, radius = 2, color = "yellow", opacity = 100,
                           fillColor = "yellow", fillOpacity = 100, label = Date)%>%
  addPolylines(lng = fall_path$Lon, lat = fall_path$Lat, group = "Fall", weight = 4, color="red", opacity = 100)%>%
  addPolylines(lng = spring_path$Lon, lat = spring_path$Lat, group = "Spring", weight = 4, color="dodgerblue", opacity = 100)%>%
  addPolylines(lng=DayDat8898[c(597,1),3], lat=DayDat8898[c(597,1),2], weight = 4, color="dodgerblue", opacity = 100,
               dashArray = 12,12) %>%
  addLegend(values = 1, group = "Fall", position = "bottomright", labels = "Fall Migration", colors = "red")%>%
  addLegend(values = 2, group = "Spring", position = "bottomright", labels = "Spring Migration", colors = "dodgerblue")%>%
  #  addCircleMarkers(lng = -101.933, lat = 35.04, color="green",
  #                  label = "Breeding colony origin", labelOptions = labelOptions(
  #                   noHide=T,textOnly=T,direction="bottom-left",style = list("color"="white")))%>%
  addCircleMarkers(lng = falldep$Lon, lat = falldep$Lat, radius=1,color="purple",opacity=100,fillColor="purple",fillOpacity=100,
                   label="Jul 6, 2018", labelOptions=labelOptions(noHide=T,direction="right",textsize="15px"))%>%
  addCircleMarkers(lng = fallarr$Lon, lat = fallarr$Lat, radius=1,color="purple",opacity=100,fillColor="purple",fillOpacity=100,
                   label="Sep 16, 2018", labelOptions=labelOptions(noHide=T,direction="right",textsize="15px"))%>%
  addCircleMarkers(lng = sprdep$Lon, lat = sprdep$Lat, radius=1,color="purple",opacity=100,fillColor="purple",fillOpacity=100,
                   label="Mar 15, 2019", labelOptions=labelOptions(noHide=T,direction="right",textsize="15px"))%>%
  addCircleMarkers(lng = sprarr$Lon, lat = sprarr$Lat, radius=1,color="purple",opacity=100,fillColor="purple",fillOpacity=100,
                   label="Apr 12, 2019", labelOptions=labelOptions(noHide=T,direction="left",textsize="15px"))%>%
  addScaleBar(position="topright")
print(map8898) #save it as 1000 x 797 png or jpeg


mapshot(map5520, url = paste0(getwd(),"/5520Map_updated.html"),remove_controls = NULL)



