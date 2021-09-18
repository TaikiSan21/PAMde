# all the way through
source('./devel/gridFunctions.R')
gpsAll <- read.csv('./devel/straights/straightPathWeffort_1641.csv', stringsAsFactors = FALSE)
gpsAll$effort <- gpsAll$straight & gpsAll$aeffort == 'on'
gpsAll$Longitude <- ifelse(gpsAll$Longitude <= 0, gpsAll$Longitude + 360, gpsAll$Longitude)
pmDets <- read.csv('./devel/yberddap/Spermies_20200310.csv', stringsAsFactors = FALSE)
pmDets <- rename(pmDets, Longitude = lon, Latitude=lat, distance = pdist)
pmDets$Longitude <- ifelse(pmDets$Longitude <= 0, pmDets$Longitude + 360, pmDets$Longitude)
gpsAll$UTC <- lubridate::ymd_hms(gpsAll$UTC)
pmDets$UTC <- lubridate::mdy_hm(pmDets$UTC)
pmDets$distance <- abs(pmDets$distance)
dsm <- Distance::ds(pmDets, key='hr')

library(dplyr)
detTest <- pmDets %>%
    # arrange(Longitude, Latitude) %>%
    # arrange(UTC) %>%
    # head(14)
    filter(survey == 1641)
gpsTest <- gpsAll %>%
    # filter(Longitude < max(detTest$Longitude),
    #        Longitude > min(detTest$Longitude),
    #        Latitude < max(detTest$Latitude),
    #        Latitude > min(detTest$Latitude))
    filter(UTC < max(detTest$UTC),
           UTC > min(detTest$UTC))
# rename all to Longitude Latitude UTC, convert itmes to POSIXct, to make detfun need distance column
allTest <- doAllGrid(gps = gpsTest,
                     bounds = NULL,
                     dets = detTest,
                     trunc_m = 20e3, #METERS
                     dsmodel = dsm,
                     pixel = NULL,
                     grid = NULL,
                     plot = FALSE)

plotGridResult(allTest)
detGridIx <- sapply(1:nrow(allTest$detections), function(x) {
    searchPoint(c(allTest$detections$Longitude[x], allTest$detections$Latitude[x]), grid = allTest$grid, gridMap = attr(allTest$grid, 'map'))
})
allTest$effort[detGridIx]
plot(allTest$grid)
lines(x=allTest$gps$Longitude, y=allTest$gps$Latitude, col='blue')
points(x=detTest$Longitude, y= detTest$Latitude, col='red')

detDistance <- 5
boundary <- dfToBounds(pmDets) # need all to be same range, say 0-360
pix <- .09
grid <- makeGrid(boundary, pixel = pix, buffer_km = detDistance, plot = T) # did with .09
grid <- connectFast(grid) # this takes long, add prog bar?
ends <- getEndPoints(gps)
plot(grid)
areaList <- lapply(ends, function(x) {
    if(is.null(x)) return(0)
    doAreaShit2(x[1,], x[2,], grid=grid, plot=TRUE)
})

effArea <- purrr::reduce(areaList, `+`)
realArea <- as.numeric(st_area(grid))
coveragePct <- round(effArea/realArea, 3) * 100

# erddap stuff
library(rerddap)
library(rerddapXtracto)

getEnv <- function(data, dataset, variable, pixel=c(0,0)) {
    rerddap::cache_delete_all(force = TRUE)
    dataInfo <- rerddap::info(datasetid = dataset)
    if(!all(variable %in% dataInfo$variables$variable_name)) {
        stop('Yo dawg we aint find those variables, best check yo spellin.')
    }
    rxtracto(dataInfo,
             parameter = variable,
             xcoord = data$Longitude,
             ycoord = data$Latitude,
             tcoord = data$UTC,
             xlen = pixel[1],
             ylen = pixel[2],
             progress_bar = TRUE)
}

dataset <- myDatasets$MURSSTmday$dataset
parameter <- c('sst')
# Some datasets have an altitude dimension. If so, then zcood must be included in the rxtracto call.
# If the dataInfo shows an altitude dimension, uncomment "zcoord <- 0" and include tcoord=tcoord in the rxtracto call.
# zcoord <- 0.
# CAN ONLY DO ONE VAR T A TIME

pmDets$UTC <- lubridate::mdy_hm(pmDets$UTC)
testSst <- getEnv(pmDets[1:10,], dataset, parameter, c(0, .09))

# okay so next thing - matching evniro data to your grid cells. Need to turn grid to a df, but then we dont really
# want to just match by centroid? like may need to average all within dakine like rxtracto do. so this would mean
# figuring out buffer in index units for each one based on grid size? which is reasonable.
# rxtracto does it by making a separate ncvar_get command for each one and just sets the box based on those bounded
# coords of x - xlen, x + xlen so each one is a new request. this seems reasonable.

# for each var:
# data -> if need ncdf, get bounds and download. once have ncdf:
# data -> start and count for ncvar_get call -> ncvar_get. We do this for each datapoint so that we can do box around??
# or should we read all at once
testData <- pmDets
testData$UTC <- lubridate::mdy_hm(testData$UTC)

testData[1, ]
buffer <- c(0, 0.09, 0)/2 # XYT VALUES
ncfile <- downloadEnv(data = testData[1:10, ], env = myDatasets$MURSSTmday, fileName = 'NcAllTest.nc', varPick = TRUE, buffer=buffer)
myEnv <- ncToData(testData[1:10,], ncFile = ncfile, buffer = buffer)


ncToData(testData[10,], ncFile = suppressWarnings(downloadEnv(testData[10,], env=myDatasets$MURSSTmday, varPick = TRUE)), quiet=TRUE)
# IT WORKS THO??? NEED TESTING AND CHECK AGIANST RXTRACTO VERSION AND ADD SMART GRID BUFFERING BUT KINDA WORKS

# downloads file, returns invisibly the fileName if it succeeded or FALSE if it failed


maybeFile <- downloadEnv(testData[1,], env = myDatasets$MURSSTmday)
nc_open(maybeFile)

# NOTEEEE
# HAVE A KEEPCOORDS OPTION TO STORE LAT LONG TIME WHATEVER WITH THE OUTPUT

# STILL DOING SHIT
# CHECK RANGE FROM DATALIST AGAINST RANGE OF YOUR DATA
# SPLIT CROSS DATELINE CALLS LOL FUCK YOU

# NOTES SEARCHPOINT WILL ONLY WORK FOR NICE CONCAVE SHAPES
