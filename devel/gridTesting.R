# test the grid maker process
list.vertices <- list(
    c(243, 32),
    c(238, 33),
    c(234, 38),
    c(234, 48.12),
    c(235.210364011, 48.5375173285),
    c(235.494102374, 48.3877105965),
    c(237, 48),
    c(237, 40),
    c(240, 36),
    c(243, 36),
    c(243, 32)
)

poly.df <- data.frame(do.call(rbind, list.vertices))

# pixel <- .225          # 25km
# pixel <- .090          # 10km
# pixel <- 0.10          # 0.1 degree for SeaGrant modeling project
# pixel <- 0.045         # 5km
pixel <- 0.027         # 3km
# pixel <- 0.018         # 2km

poly.bound <- st_sfc(st_polygon(list(as.matrix(poly.df))), crs = 4326)
plot(poly.bound, axes = TRUE, border = "red")



### Create grid (for areas) and get *centroids* of polygons.
# I dont need areas because I will be getting effort, right?
# but maybe get it anyway. need centroid to match up to enviro
# most likely
grid <- st_make_grid(
    poly.bound, cellsize = pixel,
    offset = st_bbox(poly.bound)[c("xmin", "ymin")] - pixel / 2
)
realGrid <- grid[st_intersects(poly.bound, grid)[[1]]]

grid.cent <- st_set_precision(st_centroid(realGrid), 1e+10)

### Get the centroids that are within the study area
grid.cent.which <- unlist(st_intersects(poly.bound, grid.cent))

# ask for x coords, or get bounding box first to supply?
makeGrid <- function(x, pixel=NULL, buffer_km=0, plot=FALSE) {
    # pixel <- .225          # 25km
    # pixel <- .090          # 10km
    # pixel <- 0.10          # 0.1 degree for SeaGrant modeling project
    # pixel <- 0.045         # 5km
    # pixel <- 0.027         # 3km
    # pixel <- 0.018         # 2km
    # always .009 pixel degrees / km
    # browser()
    pixOpts <- c(.225, .09, .1, .045, .027, .018)
    pixText <- paste0(' (', c('~25km', '~10km', '.1 degrees', '~5km', '~3km', '~2km'), ')')
    if(is.null(pixel)) {
        pChoice <- menu(title = 'Choose a pixel size (degrees) for creating your grid:',
                        choices = paste0(pixOpts, pixText))
        if(pChoice==0) stop('Must supply pixel size.')
        pixel <- pixOpts[pChoice]
    }
    # list of coord vertices to polgon
    polyBound <- listToPoly(x)

    # make new so can plot old later
    if(buffer_km > 0) {
        gridBound <- suppressWarnings(st_buffer(polyBound, units::as_units(buffer_km*.009, 'degrees')))
    } else {
        gridBound <- polyBound
    }
    grid <- st_make_grid(gridBound, cellsize = pixel,
                         offset = st_bbox(gridBound)[c('xmin', 'ymin')] - pixel / 2)
    # gridCent <- st_set_precision(st_centroid(grid), 1e+10)
    # grid makes box around, limit to only study area - do we need buffer?
    # browser()
    grid <- suppressMessages(grid[st_intersects(gridBound, grid)[[1]]])
    if(plot) {
        plot(grid)
        if(buffer_km > 0) {
            plot(gridBound, add=TRUE, border='darkgreen', lwd=2)
        }
        plot(polyBound, add=TRUE, border='blue', lwd=2)

    }
    grid
}

listToPoly <- function(x) {
    x <- t(matrix(unlist(x), nrow = 2))
    if(!identical(x[1, ], x[nrow(x), ])) {
        x <- rbind(x, x[1, ])
    }
    st_sfc(st_polygon(list(x)), crs = 4326)
}
testGrid <- makeGrid(list.vertices)

cg <- connectFast(smallGrid)

library(microbenchmark)

matVerts <- t(matrix(c(c(243, 32),
                     c(238, 33),
                     c(234, 38),
                     c(234, 48.12),
                     c(235.210364011, 48.5375173285),
                     c(235.494102374, 48.3877105965),
                     c(237, 48),
                     c(237, 40),
                     c(240, 36),
                     c(243, 36),
                     c(243, 32)), nrow=2))

identical(st_sfc(st_polygon(list(as.matrix(poly.df))), crs = 4326), st_sfc(st_polygon(list(matVerts)), crs=4326))
microbenchmark(
    toDf = {
        as.matrix(data.frame(do.call(rbind, list.vertices)))
    },
    fromList = {
        t(matrix(unlist(list.vertices), nrow=2))
    },
    times=1e3)

dfToBounds <- function(x, buffer_km = 0) {
    buffer <- buffer_km * .009
    # .009 pix degrees/ km
    if(is.data.frame(x) &&
       all(c('Longitude', 'Latitude') %in% colnames(x))) {
        x <- x[, c('Longitude', 'Latitude')]
    }
    list(c(min(x[,1]) - buffer, min(x[,2]) - buffer),
         c(min(x[,1]) - buffer, max(x[,2]) + buffer),
         c(max(x[,1]) + buffer, max(x[,2]) + buffer),
         c(max(x[,1]) + buffer, min(x[,2]) - buffer),
         c(min(x[,1]) - buffer, min(x[,2]) - buffer)
    )
}

gps <- read.csv('./devel/TestPath_1303AllSpeeds.csv', stringsAsFactors = FALSE)
gpsCoords <- gps[, c('Longitude', 'Latitude')]
myBounds <- dfToBounds(gps[, c('Longitude', 'Latitude')])
testGps <- makeGrid(myBounds)
plot(testGps)
bbPoly <- listToPoly(myBounds)
plot(bbPoly, add=TRUE)
points(x=gpsCoords$Longitude, y=gpsCoords$Latitude)

noBuff <- listToPoly(dfToBounds(gps, 0))
buff <- st_buffer(noBuff, units::as_units(5*.009, 'degrees'))
plot(buff)

# dis da order
myBounds <- dfToBounds(gps, buffer_km = 0)
# this a rectangle, can also be a list of vertices
detDistance <- 25
myGrid <- makeGrid(list.vertices, pixel=.225, buffer_km=detDistance, plot=TRUE) # change to size and units?
myConnect <- connectFast(myGrid)
# myBounds can also be a polygon, but make_grid makes a rectangle around it first
# then we intersect later. We need to extend this for detection distance with
# buffer

gps <- read.csv('./devel/straights/straightPathWeffort_1303.csv', stringsAsFactors = FALSE)
gps$aeffort <- ifelse(gps$aeffort == 'off', FALSE, TRUE)
gps$overallEffort <- gps$straight & gps$aeffort
plot(x=gps$Longitude, y=gps$Latitude, type='l')
bounds <- dfToBounds(gps)
grid <- makeGrid(bounds, pixel = .09, buffer_km = 15, plot=T)

doAllGrid <- function(gps, bound=NULL, pixel=NULL, buffer_km=15, plot=FALSE) {
    if(is.null(bound)) {
        bound <- dfToBounds(gps)
    }
    grid <- makeGrid(bound, pixel=pixel, buffer_km=buffer_km, plot=plot)
    if(plot) {
        lines(x=gps$Longitude, y=gps$Latitude, col='purple', lwd=2)
    }
    grid
}
gps$alt <- gps$overallEffort[c(1, 1:(nrow(gps)-1))] != gps$overallEffort
gps$effGroup <- cumsum(gps$alt)
hmm <- doAllGrid(gps, pixel=.09, buffer_km=15, plot=T)
# possible i should intersect grid with gps points.
library(ggplot2)
library(dplyr)
gps$approxDist <- gps$timeDiff * gps$Speed * .51444
gps <- select(gps, UTC, Latitude, Longitude, Speed, Heading, timeDiff, timeGroup, straight, aeffort, overallEffort, alt, effGroup, approxDist)
gps %>%
    filter(overallEffort, effGroup %in% 1:1000) %>%
    ggplot(aes(x=Longitude, y=Latitude, color=as.factor(effGroup))) +
               geom_path(show.legend = FALSE)

gps <- gps %>%
    filter(overallEffort) %>%
    group_by(effGroup) %>%
    mutate(distGroup = ceiling(cumsum(approxDist) / 1e3)) %>%
    ungroup() %>%
    mutate(distGroup = paste0(effGroup, '_', distGroup))
length(unique(gps$distGroup))

splitGps <- split(gps, gps$distGroup)
endPoints <- lapply(splitGps, function(x) {
    if(nrow(x) == 1) return(NULL)
    tmp <- as.matrix(x[c(1, nrow(x)), c('Longitude', 'Latitude')])
    attr(tmp, 'dimnames') <- NULL
    tmp
})
getEndPoints <- function(gps, length=1e3) {
    gps$aeffort <- ifelse(gps$aeffort == 'off', FALSE, TRUE)
    gps$overallEffort <- gps$straight & gps$aeffort
    gps$alt <- gps$overallEffort[c(1, 1:(nrow(gps)-1))] != gps$overallEffort
    gps$effGroup <- cumsum(gps$alt)
    # gps$approxDist <- gps$timeDiff * gps$Speed * .51444
    gps$approxDist <- geosphere::distGeo(gps[c(1, 1:(nrow(gps)-1)), c('Longitude', 'Latitude')],
                                         gps[,c('Longitude', 'Latitude')])
    gps <- gps %>%
        filter(overallEffort) %>%
        group_by(effGroup) %>%
        mutate(distGroup = ceiling(cumsum(approxDist) / length)) %>%
        ungroup() %>%
        mutate(distGroup = paste0(effGroup, '_', distGroup))
    lapply(split(gps, gps$distGroup), function(x) {
        if(nrow(x) == 1) return(NULL)
        tmp <- as.matrix(x[c(1, nrow(x)), c('Longitude', 'Latitude')])
        attr(tmp, 'dimnames') <- NULL
        tmp
    })
}

grid <- connectFast(grid)
miniGps <- filter(gps, effGroup %in% 1:3)
miniGrid <- makeGrid(dfToBounds(miniGps), .09, buffer_km=5)
miniGrid <- connectFast(miniGrid)
miniEnds <- getEndPoints(miniGps)
plot(miniGrid)
n <- 1
hmm <- lapply(miniEnds, function(x) {
    cat(n)
    n <<- n+1
    if(is.null(x)) return(NULL)
    doAreaShit2(x[1,], x[2,], grid = miniGrid, plot = T)
})
purrr::reduce(hmm ,`+`)
lines(x=miniGps$Longitude, y=miniGps$Latitude, col='red')

library(geosphere)

testDoAll <- function(gps, pixel=NULL, buffer_km) {
    grid <- makeGrid(dfToBounds(gps), pixel, buffer_km, plot = FALSE)
    grid <- connectFast(grid)
    ends <- getEndPoints(gps)
    plot(grid)
    areaList <- lapply(ends, function(x) {
        if(is.null(x)) return(0)
        doAreaShit2(x[1,], x[2,], grid=grid, plot=TRUE)
    })

    effArea <- purrr::reduce(areaList, `+`)
    realArea <- as.numeric(st_area(grid))
    coveragePct <- round(effArea/realArea, 3) * 100
    # plot(grid[coveragePct > 100], col='red', add=TRUE)
    # browser()
    # text(x = sapply(grid, function(x) {
    #   mean(x[[1]][1:2,1])
    # }),
    # y = sapply(grid, function(y) {
    #   mean(y[[1]][2:3,2])
    # }),
    # label = coveragePct
    # )
    coveragePct
}

ntest <- 10e3
testAll <- testDoAll(gps=gps[1:ntest,], buffer_km=5)

debugGrid <- makeGrid(dfToBounds(gps[1:ntest,]), .09, buffer_km=5, plot = F)
plot(debugGrid, col = gray(1-testAll/max(testAll), alpha=.9))
lines(x=gps[1:ntest, 'Longitude'], y=gps[1:ntest, 'Latitude'], col='purple')
