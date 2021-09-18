# all the way through
source('gridFunctions.R')
# EXAMPLE WORK FLOW
# 1) First read in gps, then you need to convert longitudes to 0-360,
# convert UTC to POSIXct format, and create a column called "effort" that is TRUE or FALSE
# Column names needed: Longitude, Latitude, UTC, effort
gpsAll <- read.csv('./devel/straights/straightPathWeffort_1641.csv', stringsAsFactors = FALSE)
gpsAll$effort <- gpsAll$straight & gpsAll$aeffort == 'on'
gpsAll$Longitude <- ifelse(gpsAll$Longitude <= 0, gpsAll$Longitude + 360, gpsAll$Longitude)
gpsAll$UTC <- lubridate::ymd_hms(gpsAll$UTC)

# 2) Read in detections, similar to above convert column names to Longitude, Latitude, UTC,
# and "distance" if you need to create a detection function using Distance package
# Longitude must be matching 0-360 like above
pmDets <- read.csv('./devel/yberddap/Spermies_20200310.csv', stringsAsFactors = FALSE)
pmDets <- rename(pmDets, Longitude = lon, Latitude=lat, distance = pdist)
pmDets$Longitude <- ifelse(pmDets$Longitude <= 0, pmDets$Longitude + 360, pmDets$Longitude)
pmDets$UTC <- lubridate::mdy_hm(pmDets$UTC)
pmDets$distance <- abs(pmDets$distance)

# 3) fit a detection function if you need one
dsm <- Distance::ds(pmDets, key='hr')

# 4) Run grid stuff. This will take a while, and should have a couple progress bars showing you how long the slow parts
# are taking. Set "trunc_m" to whatever truncation distance you might want to use
# It will ask you what pixel size you want  (this is the grid size), in the future you can specify this when you call it
debug(doAllGrid)
gridEffort <- doAllGrid(gps = gpsAll,
                     bounds = NULL,
                     dets = pmDets,
                     trunc_m = 20e3, #METERS
                     dsmodel = dsm,
                     pixel = NULL,
                     grid = NULL,
                     plot = FALSE)

# Result of this is a list with the gps you provided, the grid it created, the effort values for each grid, and
# the detections you provided (these have the effortArea attached to them, as well as the actualArea of that grid
# cell - actualArea is mostly there as a reference point for sanity check, prob not used for much)

# Theres also a plotty function, but it will look pretty messy for all your data. Trackline is (roughly) in blue,
# detections are red dots, amount of effort is shaded in gray, but for a big grid you wont really be able to see the shading
plotGridResult(gridEffort)