library(dplyr)
library(sf)
library(Matrix)
library(geosphere)
library(Distance)
library(mrds)

# updated 07-02-2020 to add times to effort values as $first $last
calcEffort <- function(points, grid, n=10, d=4e3, dsmodel, plot=FALSE) {
    if(is.null(attr(grid, 'map'))) {
        grid <- connectGrid(grid)
    }
    times <- points$times
    points <- points$points
    gridMap <- attr(grid, 'map')
    area <- rep(0, length(grid))
    first <- rep(NA, length(grid))
    datesAdded <- vector('list', length = length(grid))
    last <- rep(NA, length(grid))
    cat('Calculating effort...\n')
    pb <- txtProgressBar(min=0, max=length(points), style=3)
    # browser()
    # detFunList <- for each unique combo of covars, assign id and eval the detection probabilites
    # and for each point store that detFun number. If it adds a lot of time to access this from each
    # point cuz is big, create a separate detFun map here and use that
    # browser()
    detFun <- makeDetFun(data.frame(x=1), dsmodel) # just need a dummy df here if no covars, eventaully covars from grid
    probs <- detFun(seq(from=0, to=d, length.out=n+1))
    # seq(from=0, to=d, length.out=n+1)

    for(i in seq_along(points)) {
        cat('\rNum', i)
        setTxtProgressBar(pb, value=i)
        # mybe check here for any added areaz.
        # browser()
        newArea <- calcArea(points[[i]][1, ], points[[i]][2, ], grid=grid, detProb = probs, gridMap=gridMap, n=n, d=d, plot=plot)
        whichNew <- which(newArea > 0)
        if(length(whichNew) > 0) {
            for(j in whichNew) {
                if(is.na(first[j])) {
                    first[j] <- times[i]
                    datesAdded[[j]] <- times[i]
                }
                last[j] <- times[i]
                datesAdded[[j]] <- c(datesAdded[[j]], times[i])
            }
        }
        area <- area + newArea
    }
    cat('\n')
    first <- as.POSIXct(first, tz='UTC', origin='1970-01-01')
    last <- as.POSIXct(last, tz='UTC', origin='1970-01-01')
    list(area=area, first=first, last=last, datesAdded=datesAdded)
}

connectGrid <- function(grid, pbshow=FALSE) {
    # gridMat <- matrix(unlist(grid), nrow = 5)
    conMat <- matrix(NA, ncol = 4, nrow = length(grid))
    cat('Connecting grid...\n')
    if(pbshow) {
        pb <- txtProgressBar(min=0, max=length(grid), style=3)
    }
    for(i in seq_along(grid)) {
        for(side in 1:4) {
            if(!is.na(conMat[i, side])) {
                next
            }
            matched <- matchEdge(i, side, grid, start=i)
            conMat[i, side] <- matched
            if(matched == -1) {
                next
            }
            conMat[matched, (side+1) %% 4 + 1] <- i
        }
        if(pbshow) setTxtProgressBar(pb, value=i)
    }
    cat('\n')
    class(grid) <- c(class(grid), 'connectedGrid')
    attr(grid, 'map') <- conMat
    grid
}

matchEdge <- function(id, edgeNo, grid, start=2) {
    # we grab the opposing edge for each grid to check, and we only need to check if a single
    # point is equal since all straight lines
    thisEdge <- .subset2(.subset2(grid, id), 1)[edgeNo, ]
    nextEdgeNo <- (edgeNo + 1) %% 4 + 1
    for(i in start:length(grid)) {
        nextEdge <- .subset2(.subset2(grid, i), 1)[(nextEdgeNo+1), ]
        if(thisEdge[1] == nextEdge[1] &&
           thisEdge[2] == nextEdge[2]) {
            return(i)
        }
    }
    -1
}

calcArea <- function(point1, point2,grid, gridMap, detProb=NULL, n=10, d=4e3, plot = FALSE) {
    # browser()
    length <- geosphere::distGeo(point1, point2)
    # cat('Length :', length, '\n')
    points <- rbind(makeDetLines(point1, point2, n, d, direction = 90)[(n+1):1,],
                    makeDetLines(point1, point2, n, d, direction = -90)[-1,])
    gridLength <- apply(points, 1, function(x) gridPropCalc(x, grid, gridMap))
    gridLength <- as(gridLength, 'dgCMatrix') * length
    if(is.null(detProb)) {
        detFun <- rep(1, length.out=n+1)
    }
    # probs <- detFun(seq(from=0, to=d, length.out=n+1))

    # probs <- c(rev(probs), probs[-1])
    probs <- c(rev(detProb), detProb[-1])
    # all but last, h2 all but first. same for probs w1,w2
    h1 <- gridLength[, 1:(2*n)]
    # h2 <- gridLength[, 2:(2*n + 1)]
    h2 <- gridLength[, -1]
    w1 <- probs[1:(2*n)]
    w2 <- probs[-1]
    # browser()
    areas <- .5*(h1*w1 + h2*w2) - (h1-h2)*(w1-w2)/6
    gridArea <- rowSums(areas) * d / n
    if(plot) {
        # plot(grid)
        for(i in 1:nrow(points)) {
            lines(x=points[i, c(1, 3)], y=points[i, c(2, 4)], col = 'blue')
        }
        points(x=c(point1[1], point2[1]),
               y=c(point1[2], point2[2]),
               col = 'black', pch=16)
        lines(x=c(point1[1], point2[1]),
              y=c(point1[2], point2[2]),
              col = 'black', lwd=2)
    }
    browser()
    gridArea
}

makeDetLines <- function(point1, point2, n=10, d = 4e3, direction = 90) {
    suppressWarnings({
        pathBearing <- geosphere::bearing(point1, point2)
        end1 <- geosphere::destPoint(point1, pathBearing + direction, d) %% 360
        end2 <- geosphere::destPoint(point2, pathBearing + direction, d) %% 360
    })
    lats1 <- seq(from = point1[2], to = end1[2], length.out = n + 1)
    longs1 <- seq(from = point1[1], to = end1[1], length.out = n + 1)
    lats2 <- seq(from = point2[2], to = end2[2], length.out = n + 1)
    longs2 <- seq(from = point2[1], to = end2[1], length.out = n + 1)
    # list(start = matrix(c(longs1, lats1), ncol = 2),
    #      end = matrix(c(longs2, lats2), ncol = 2))
    matrix(c(longs1, lats1, longs2, lats2), ncol=4)
}

gridPropCalc <- function(points, grid,  gridMap) {
    # If we do a search thing, need to check if one/both points are off grid
    # broken if both points are off grid - this doesnt mean nothing exists,
    # line can still intersect
    # browser()
    thisGrid <- searchPoint(c(points[1], points[2]), grid, gridMap)
    if(thisGrid == -1) {
        points <- points[c(3,4,1,2)]
        thisGrid <- searchPoint(c(points[1], points[2]), grid, gridMap)
    }
    props <- vector('integer', length(grid))
    if(thisGrid == -1) {
        # if both start and end are off the grid, cant do search method
        # browser()
        cat('\nYOURE OFF THE GRID ITS SCARY OUT HERE. AND ALSO SLOW. TELL TAIKI THINGS BROKE.\n')
        return(gridPropMat(points, grid))
    }
    toGo <- 1
    thisSide <- 5
    lastProp <- 0
    for(i in 1:1e3) {
        thisBox <- .subset2(.subset2(grid, thisGrid), 1)
        thisProp <- 0
        for(side in 1:4) {
            if(side == thisSide) {
                next
            }
            thisProp <- getIntersection(points[1], points[2], points[3], points[4],
                                        thisBox[side, 1], thisBox[side, 2], thisBox[side+1, 1], thisBox[side+1, 2])
            if(thisProp > 0) {
                break
            }
        }
        if(thisProp == 0) {
            props[thisGrid] <- toGo
            break
        }
        props[thisGrid] <- thisProp - lastProp
        toGo <- toGo - thisProp + lastProp
        lastProp <- thisProp
        thisSide <- (side + 1) %% 4 + 1
        thisGrid <- gridMap[thisGrid, side]
        # cat('grid ', thisGrid,'point ', thisSide, '\n')
        if(thisGrid == -1 ||
           toGo < 1e-8) {
            break
        }
        # i <- i+1
    }
    props
}

searchPoint <- function(point, grid, gridMap) {
    i <- 0
    thisGrid <- 1
    for(i in  1:nrow(gridMap)) {
        if(thisGrid == -1) {
            break
        }
        # thisBox <- gridMat[, (2*thisGrid - 1):(2*thisGrid)]
        thisBox <- .subset2(.subset2(grid, thisGrid), 1)
        # is to east?
        if(point[1] > thisBox[2, 1]) {
            nextGrid <- gridMap[thisGrid, 2]
            if(nextGrid != -1) {
                thisGrid <- nextGrid
                next
            }
        }
        # is to west?
        if(point[1] < thisBox[1,1]) {
            nextGrid <- gridMap[thisGrid, 4]
            if(nextGrid != -1) {
                thisGrid <- nextGrid
                next
            }
        }
        # is to north?
        if(point[2] > thisBox[3, 2]) {
            nextGrid <- gridMap[thisGrid, 3]
            if(nextGrid != -1) {
                thisGrid <- nextGrid
                next
            }
        }
        # is to south?
        if(point[2] < thisBox[1, 2]) {
            nextGrid <- gridMap[thisGrid, 1]
            if(nextGrid != -1) {
                thisGrid <- nextGrid
                next
            }
        }
        return(thisGrid)
    }
    -1
}

getIntersection <- function(t1x, t1y, t2x, t2y, g1x, g1y, g2x, g2y) {
    # t  transect g  grid
    s1x <- t2x - t1x
    s1y <- t2y - t1y
    s2x <- g2x - g1x
    s2y <- g2y - g1y
    # browser()
    denom <- (-s2x * s1y + s1x * s2y)
    if(abs(denom) < 1e-8) {
        return(FALSE)
    } else {
        denom <- 1/denom
    }
    s <- (-s1y * (t1x - g1x) + s1x * (t1y - g1y)) * denom
    if(s >= 0 && s <= 1) {
        t <- (s2x * (t1y - g1y) - s2y * (t1x - g1x)) * denom
        if(t >=0 && t <= 1) {
            return(t)
        }
    }
    FALSE
}

getEndPoints <- function(gps, length=1e3) {
    if(!('effort' %in% colnames(gps))) {
        stop('Gps data must have an "effort" column that is TRUE/FALSE')
    }
    # gps$aeffort <- ifelse(gps$aeffort == 'off', FALSE, TRUE)
    # gps$overallEffort <- gps$straight & gps$aeffort
    # gps$alt <- gps$overallEffort[c(1, 1:(nrow(gps)-1))] != gps$overallEffort
    gps$alt <- gps$effort[c(1, 1:(nrow(gps)-1))] != gps$effort
    gps$effGroup <- cumsum(gps$alt)
    # gps$approxDist <- gps$timeDiff * gps$Speed * .51444
    gps$approxDist <- geosphere::distGeo(gps[c(1, 1:(nrow(gps)-1)), c('Longitude', 'Latitude')],
                                         gps[,c('Longitude', 'Latitude')])
    gps <- gps %>%
        filter(effort) %>%
        group_by(effGroup) %>%
        mutate(distGroup = ceiling(cumsum(approxDist) / length)) %>%
        ungroup() %>%
        mutate(distGroup = paste0(effGroup, '_', distGroup))
    # browser()
    result <- lapply(split(gps, gps$distGroup), function(x) {
        if(nrow(x) == 1) return(NULL)
        tmp <- as.matrix(x[c(1, nrow(x)), c('Longitude', 'Latitude')])
        attr(tmp, 'dimnames') <- NULL
        tmp
    })
    times <- group_by(gps, distGroup) %>%
        summarise(time = min(UTC))
    whichNotNull <- sapply(result, function(x) !is.null(x))
    times <- times$time[whichNotNull]
    result <- result[whichNotNull]
    timeOrder <- order(times)
    # result[sapply(result, function(x) !is.null(x))]
    list(points=result[timeOrder], times=times[timeOrder])
}

makeGrid <- function(x, pixel=NULL, buffer_m=0, plot=FALSE) {
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
    if(buffer_m > 0) {
        gridBound <- suppressWarnings(st_buffer(polyBound, units::as_units(buffer_m * .009 / 1000, 'degrees')))
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
        if(buffer_m > 0) {
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

dfToBounds <- function(x, buffer_m = 0) {
    buffer <- buffer_m * .009 / 1000
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

makeDetFun <- function(covarDf, dsModel=NULL) {
    # base test case perfect detection
    if(is.null(dsModel) ||
       is.numeric(dsModel)) {
        if(is.null(dsModel)) {
            probs <- 1
        } else {
            probs <- dsModel
        }
        return(
            function(distance) {
                rep(probs, length.out = length(distance))
            }
        )
    }
    UNITDENOM <- 1
    if(dsModel$ddf$meta.data$width < 100) {
        warning('Looks like your detection was fit using kilometers, expected meters. Converting units.')
        UNITDENOM <- 1e3
    }
    function(distance) {
        distance <- distance / UNITDENOM
        lessLeft <- which(distance < dsModel$ddf$meta.data$left)

        if(nrow(covarDf) > 1) {
            warning('Cannot create detection function for more than 1 set of covariates at a time.')
            covarDf <- covarDf[1, ]
        }
        result <- eval_with_covars(distance, newdata = covarDf, model=dsModel)
        if(length(lessLeft) > 0) {
            warning('Some distances less than left truncation distance, setting to 0 probability.')
            result[lessLeft] <- 0
        }
        result
    }
}

eval_with_covars <- function(distance, newdata, model){
    if(inherits(model, 'dsmodel')) {
        model <- model$ddf
    }
    ddfobj <- model$ds$aux$ddfobj
    left <- model$meta.data$left
    width <- model$meta.data$width
    fpar <- model$par
    ddfobj <- assign.par(ddfobj, fpar)
    # Get integration ranges either from specified argument or from
    # values stored in the model.
    if(is.data.frame(newdata)){
        nr <- nrow(newdata)
    }else{
        nr <- 1
    }

    if(is.null(model$ds$aux$int.range)){
        int.range <- cbind(rep(0, nr), rep(width, nr))
    }else{
        int.range <- model$ds$aux$int.range
        if(is.vector(int.range)){
            int.range <- cbind(rep(int.range[1], nr),
                               rep(int.range[2], nr))
            #}else if(nrow(int.range) == (nrow(x)+1)){
            #int.range <- int.range[2:nrow(int.range), , drop=FALSE]
        }
    }

    # set the distance column to be the left truncation distance
    # this gets around an issue that Nat Kelly found where later process.data
    # will remove entires with distance < left truncation
    # BUT save the NAs!
    nas <- is.na(newdata$distance)
    newdata$distance <- left
    newdata$distance[nas] <- NA

    newdata_save <- newdata

    # get the data in the model
    model_dat <- model$data

    # counter for NAs...
    naind <- rep(FALSE, nrow(newdata))

    # do this for both scale and shape parameters
    for(df_par in c("scale", "shape")){
        # if that parameter exists...
        if(!is.null(ddfobj[[df_par]])){
            # save the column names from the design matrix
            znames <- colnames(ddfobj[[df_par]]$dm)

            # pull out the columns in the formula and the distances column
            fvars <- all.vars(as.formula(model$ds$aux$ddfobj[[df_par]]$formula))

            if(!all(fvars %in% colnames(newdata))){
                stop("columns in `newdata` do not match those in fitted model\n")
            }

            model_dat <- model_dat[, c("distance", fvars), drop=FALSE]

            if(df_par=="scale"){
                # which rows have NAs?
                naind <- naind | apply(newdata_save[, c("distance", fvars), drop=FALSE],
                                       1, function(x) any(is.na(x)))
            }

            # setup the covariate matrix, using the model data to ensure that
            # the levels are right
            newdata <- rbind(model_dat,
                             newdata_save[, c("distance", fvars), drop=FALSE])
            dm <- mrds:::setcov(newdata, as.formula(ddfobj[[df_par]]$formula))

            # now check that the column names are the same for the model
            # and prediction data matrices
            if(!identical(colnames(dm), znames)){
                stop("fields or factor levels in `newdata` do not match data used in fitted model\n")
            }

            # get only the new rows for prediction
            dm <- dm[(nrow(model_dat)+1):nrow(dm), , drop=FALSE]
            # assign that!
            ddfobj[[df_par]]$dm <- dm

        }
    }

    # handle data setup for uniform key case
    if(ddfobj$type == "unif"){
        model_dat <- model_dat[, "distance", drop=FALSE]
        # which rows have NAs?
        naind <- is.na(newdata_save$distance)

        newdata <- rbind(model_dat,
                         newdata_save[, "distance", drop=FALSE])
        dm <- setcov(newdata, ~1)
        dm <- dm[(nrow(model_dat)+1):nrow(dm), , drop=FALSE]
    }

    # get the bins when you have binned data
    # use the breaks specified in the model!
    if(model$meta.data$binned){
        nanana <- apply(newdata[, c("distance", fvars), drop=FALSE],
                        1, function(x) any(is.na(x)))
        newdata_b <- create.bins(newdata[!nanana, , drop=FALSE], model$meta.data$breaks)
        newdata$distbegin <- NA
        newdata$distend <- NA
        newdata[!nanana, ] <- newdata_b
    }

    # update xmat too
    datalist <- mrds:::process.data(newdata, model$meta.data, check=FALSE)
    ddfobj$xmat <- datalist$xmat[(nrow(model_dat)+1):nrow(datalist$xmat),,drop=FALSE]
    ddfobj$xmat <- ddfobj$xmat[!naind, , drop=FALSE]
    int.range <- int.range[!naind, , drop=FALSE]
    # reset newdata to be the right thing
    newdata <- newdata[(nrow(model_dat)+1):nrow(newdata), , drop=FALSE]


    detfct(distance, ddfobj, select=NULL, index=NULL, width=width,
           standardize=TRUE, stdint=FALSE, left=left)
}

# if bounds, chekc all gps and dets within bounds or we suffer the slow fate of matrix search
doAllGrid <- function(gps, dets, bounds=NULL, trunc_m, dsmodel=NULL, pixel=NULL, grid=NULL, plot=FALSE) {
    # browser()
    # if(trunc > 100) {
    #     cont <- menu(title = paste0('Truncation distance should be in kilometers, are you sure ',
    #                                 'you want to use ', trunc,'?'),
    #                  choices = c('Yes', 'No, stop and set a new value.'))
    #     if(cont == 2) return(NULL)
    # }
    #
    if(is.character(gps)) {
        gps <- read.csv(gps, stringsAsFactors = FALSE)
    }
    if(is.character(dets)) {
        dets <- read.csv(dets, stringsAsFactors = FALSE)
    }
    # if(any(dets$distance) > 100) {
    #     warning('Detection ')
    # }
    # should I check names and format all to 0-360 in a helper?

    if(is.null(bounds)) {
        boundary <- dfToBounds(rbind(gps[, c('Longitude', 'Latitude')],
                                     dets[, c('Longitude', 'Latitude')])) # need all to be same range, say 0-
    } else {
        boundary <- dfToBounds(bounds)
    }
    # pixel <- .225          # 25km
    # pixel <- .090          # 10km
    # pixel <- 0.10          # 0.1 degree for SeaGrant modeling project
    # pixel <- 0.045         # 5km
    # pixel <- 0.027         # 3km
    # pixel <- 0.018         # 2km
    # always .009 pixel degrees / km
    pixOpts <- c(.225, .09, .1, .045, .027, .018)
    pixText <- paste0(' (', c('~25km', '~10km', '.1 degrees', '~5km', '~3km', '~2km'), ')')
    if(is.null(pixel)) {
        pChoice <- menu(title = 'Choose a pixel size (degrees) for creating your grid:',
                        choices = paste0(pixOpts, pixText))
        if(pChoice==0) stop('Must supply pixel size.')
        pixel <- pixOpts[pChoice]
    }
    if(is.null(grid)) {
        grid <- makeGrid(boundary, pixel = pixel, buffer_m = trunc_m * 1.05, plot = plot) # did with .09
    }
    if(inherits(grid, 'connectedGrid')) {
        conGrid <- grid
    } else {
        conGrid <- connectGrid(grid, pbshow = TRUE) # this takes long, add prog bar?
    }
    # this needs "effort" column, so we do this first
    ends <- getEndPoints(gps, length = 1e3) # length of segs to break into, less should more accurato
    # browser()
    effort <- calcEffort(points=ends, grid=conGrid, n=10, d=trunc_m, dsmodel=dsmodel)
    detGridIx <- sapply(1:nrow(dets), function(x) {
        searchPoint(c(dets$Longitude[x], dets$Latitude[x]), grid = conGrid, gridMap = attr(conGrid, 'map'))
    })
    dets$gridIndex <- detGridIx
    dets$effortArea <- effort$area[detGridIx]
    dets$actualArea <- as.numeric(st_area(conGrid[detGridIx]))
    result <- list(gps=gps, grid=conGrid, effort=effort, detections=dets)
    # browser()

    if(plot) {
        plotGridResult(result)
    }
    result
}

plotGridResult <- function(x) {
    actualArea <- as.numeric(st_area(x$grid))
    coveragePct <- round(x$effort$area / actualArea, 3)
    coveragePct <- ifelse(coveragePct > 1, 1, coveragePct)
    plot(x$grid, col = gray(1 - coveragePct, alpha = .9))
    lines(x=x$gps$Longitude, y=x$gps$Latitude, col='blue')
    points(x=x$detections$Longitude, y=x$detections$Latitude, col='red')
}

# effArea <- purrr::reduce(areaList, `+`)
# realArea <- as.numeric(st_area(grid))
# coveragePct <- round(effArea/realArea, 3) * 100
# plot(debugGrid, col = gray(1-testAll/max(testAll), alpha=.9))