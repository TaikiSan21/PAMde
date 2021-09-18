# testing on grid from gridmaker_SMW
library(geosphere)
library(sf)
library(dplyr)
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

makeDetLines <- function(point1, point2, n=10, d = 4e3, direction = 90) {
  suppressWarnings({
    pathBearing <- bearing(point1, point2)
    end1 <- destPoint(point1, pathBearing + direction, d) %% 360
    end2 <- destPoint(point2, pathBearing + direction, d) %% 360
  })
  lats1 <- seq(from = point1[2], to = end1[2], length.out = n + 1)
  longs1 <- seq(from = point1[1], to = end1[1], length.out = n + 1)
  lats2 <- seq(from = point2[2], to = end2[2], length.out = n + 1)
  longs2 <- seq(from = point2[1], to = end2[1], length.out = n + 1)
  # list(start = matrix(c(longs1, lats1), ncol = 2),
  #      end = matrix(c(longs2, lats2), ncol = 2))
  matrix(c(longs1, lats1, longs2, lats2), ncol=4)
}

doAreaShit2 <- function(point1, point2, grid, n=10, d=4e3, plot = FALSE) {
  length <- distGeo(point1, point2)
  # cat('Length :', length, '\n')
  points <- rbind(makeDetLines(point1, point2, n, d, direction = 90)[(n+1):1,],
                  makeDetLines(point1, point2, n, d, direction = -90)[-1,])
  gridMat <- matrix(unlist(grid), nrow=5)
  if(is.null(attr(grid, 'map'))) {
    grid <- connectFast(grid)
  }
  gridMap <- attr(grid, 'map')
  gridLength <- apply(points, 1, function(x) gridPropSearch2(x, gridMat, gridMap)) * length
  probs <- rep(1, n+1) # function of d, 0 to n. H
  # probs <- seq(from=1, to=.5, length.out = n+1)
  probs <- c(rev(probs), probs[-1])

  h1 <- gridLength[, 1:(2*n)]
  # h2 <- gridLength[, 2:(2*n + 1)]
  h2 <- gridLength[, -1]
  w1 <- probs[1:(2*n)]
  w2 <- probs[-1]

  areas <- .5*(h1*w1 + h2*w2) - (h1-h2)*(w1-w2)/6
  gridArea <- apply(areas, 1, sum) * d / n
  if(plot) {
    # plot(grid)

    for(i in 1:nrow(points)) {
      lines(x=points[i, c(1, 3)], y=points[i, c(2, 4)], col = 'blue')
    }
    # text(x = sapply(grid, function(x) {
    #   mean(x[[1]][1:2,1])
    # }),
    # y = sapply(grid, function(y) {
    #   mean(y[[1]][2:3,2])
    # }),
    # label = round(gridArea / 1e6, 3)
    # )
    points(x=c(point1[1], point2[1]),
           y=c(point1[2], point2[2]),
           col = 'black', pch=16)
    lines(x=c(point1[1], point2[1]),
          y=c(point1[2], point2[2]),
          col = 'black', lwd=2)
  }
  gridArea
}

# need this as backup to search
gridPropMat <- function(points, gridMat) {
  # If we do a search thing, need to check if one/both points are off grid
  startGrid <- fastPoint(c(points[1], points[2]), gridMat)
  endGrid <- fastPoint(c(points[3], points[4]), gridMat)
  props <- rep(0, ncol(gridMat)/2)
  toGo <- 1
  for(g in 1:(ncol(gridMat)/2)) {
    thisBox <- gridMat[, (2*g-1):(2*g)]
    thisProp <- rep(0, 4)
    for(side in 1:4) {
      thisProp[side] <- getIntersection(points[1], points[2], points[3], points[4],
                                        thisBox[side, 1], thisBox[side, 2], thisBox[side+1, 1], thisBox[side+1, 2])
    }
    thisProp <- thisProp[thisProp > 0]
    if(length(thisProp) == 0) {
      next
    }
    if(length(thisProp) == 2) {
      thisProp <- abs(thisProp[1]-thisProp[2])
    }
    if(length(thisProp) == 1 &&
       g == endGrid) {
      thisProp <- 1 - thisProp
    }
    props[g] <- thisProp
    toGo <- toGo - thisProp
    if(toGo < 1e-8) {
      break
    }
  }
  props
}

# this can fail if both points of grid, could still have intersect in between
# so have backup of gridPropMat
gridPropSearch2 <- function(points, gridMat, gridMap) {
  # If we do a search thing, need to check if one/both points are off grid
  # broken if both points are off grid - this doesnt mean nothing exists,
  # line can still intersect
  # browser()
  thisGrid <- searchPoint(c(points[1], points[2]), gridMat, gridMap)
  if(thisGrid == -1) {
    points <- points[c(3,4,1,2)]
    thisGrid <- searchPoint(c(points[1], points[2]), gridMat, gridMap)
  }
  # gridMat is lat/long so num grids * 2 is ncol
  props <- vector('integer', ncol(gridMat)/2)
  if(thisGrid == -1) {
    # if both start and end are off the grid, cant do search method
    return(gridPropMat(points, gridMat))
  }
  toGo <- 1
  i <- 1
  thisSide <- 5
  lastProp <- 0
  for(i in 1:1e3) {
    thisBox <- gridMat[, (2*thisGrid-1):(2*thisGrid)]
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

connectFast <- function(grid, pbshow=FALSE) {
  gridMat <- matrix(unlist(grid), nrow = 5)
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
      if(i==1) {
        matched <- matchEdge(i, side, gridMat, start=i)
      } else{
        matched <- matchEdge(i, side, gridMat[, -c(1:(2*(i-1)))], start=i)
      }
      conMat[i, side] <- matched
      if(matched == -1) {
        next
      }
      conMat[matched, (side+1) %% 4 + 1] <- i
    }
    if(pbshow) setTxtProgressBar(pb, value=i)
  }
  class(grid) <- c(class(grid), 'connectedGrid')
  attr(grid, 'map') <- conMat
  grid
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
      matched <- matchList(i, side, grid, start=i)
      conMat[i, side] <- matched
      if(matched == -1) {
        next
      }
      conMat[matched, (side+1) %% 4 + 1] <- i
    }
    if(pbshow) setTxtProgressBar(pb, value=i)
  }
  class(grid) <- c(class(grid), 'connectedGrid')
  attr(grid, 'map') <- conMat
  grid
}

fastPoint <- function(point, gridMat) {
  for(i in 1:(ncol(gridMat)/2)) {
    thisBox <- gridMat[1:3, (2*i - 1):(2*i)]
    if(point[1] >= thisBox[1, 1] &&
       point[1] <= thisBox[2, 1] &&
       point[2] >= thisBox[1, 2] &&
       point[2] <= thisBox[3, 2]) {
      return(i)
    }
  }
  -1
}
searchPoint <- function(point, gridMat, gridMap) {
  i <- 0
  thisGrid <- 1
  for(i in  1:nrow(gridMap)) {
    if(thisGrid == -1) {
      break
    }
    thisBox <- gridMat[, (2*thisGrid - 1):(2*thisGrid)]
    if(point[1] > thisBox[2, 1]) {
      thisGrid <- gridMap[thisGrid, 2]
      next
    }
    if(point[1] < thisBox[1,1]) {
      thisGrid <- gridMap[thisGrid, 4]
      next
    }
    if(point[2] > thisBox[3, 2]) {
      thisGrid <- gridMap[thisGrid, 3]
      next
    }
    if(point[2] < thisBox[1, 2]) {
      thisGrid <- gridMap[thisGrid, 1]
      next
    }
    return(thisGrid)
  }
  -1
}

matchEdge <- function(id, edgeNo, gridMat, start=2) {
  id <- id - start + 1
  thisEdge <- gridMat[edgeNo:(edgeNo+1), (2*id-1):(2*id)]
  nextEdgeNo <- (edgeNo + 1) %% 4 + 1
  for(i in 1:(ncol(gridMat)/2)) {
    nextEdge <- gridMat[nextEdgeNo:(nextEdgeNo+1), (2*i-1):(2*i)]
    if(thisEdge[1, 1] == nextEdge[2, 1] &&
       thisEdge[1, 2] == nextEdge[2, 2]) {
      return(i+start-1)
    }
  }
  -1
}

matchList <- function(id, edgeNo, grid, start=2) {
  # id <- id - start + 1
  # thisEdge <- gridMat[edgeNo:(edgeNo+1), (2*id-1):(2*id)]
  # thisEdge <- .subset2(.subset2(grid, id), 1)[edgeNo:(edgeNo+1), ]
  thisEdge <- .subset2(.subset2(grid, id), 1)[edgeNo, ]
  nextEdgeNo <- (edgeNo + 1) %% 4 + 1
  for(i in start:length(grid)) {
    # nextEdge <- gridMat[nextEdgeNo:(nextEdgeNo+1), (2*i-1):(2*i)]
    # nextEdge <- .subset2(.subset2(grid, i), 1)[nextEdgeNo:(nextEdgeNo+1), ]
    nextEdge <- .subset2(.subset2(grid, i), 1)[(nextEdgeNo+1), ]
    # if(thisEdge[1, 1] == nextEdge[2, 1] &&
    #    thisEdge[1, 2] == nextEdge[2, 2]) {
    if(thisEdge[1] == nextEdge[1] &&
       thisEdge[2] == nextEdge[2]) {
      return(i)
    }
  }
  -1
}

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

  result <- lapply(split(gps, gps$distGroup), function(x) {
    if(nrow(x) == 1) return(NULL)
    tmp <- as.matrix(x[c(1, nrow(x)), c('Longitude', 'Latitude')])
    attr(tmp, 'dimnames') <- NULL
    tmp
  })
  result[sapply(result, function(x) !is.null(x))]
}

areaFast <- function(point1, point2, grid, n=10, d=4e3, plot = FALSE) {
  length <- distGeo(point1, point2)
  # cat('Length :', length, '\n')
  points <- rbind(makeDetLines(point1, point2, n, d, direction = 90)[(n+1):1,],
                  makeDetLines(point1, point2, n, d, direction = -90)[-1,])
  gridMat <- matrix(unlist(grid), nrow=5)
  if(is.null(attr(grid, 'map'))) {
    grid <- connectFast(grid)
  }
  gridMap <- attr(grid, 'map')
  gridLength <- apply(points, 1, function(x) gridPropSearch2(x, gridMat, gridMap)) * length
  probs <- rep(1, n+1) # function of d, 0 to n. H
  # probs <- seq(from=1, to=.5, length.out = n+1)
  probs <- c(rev(probs), probs[-1])

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
    # text(x = sapply(grid, function(x) {
    #   mean(x[[1]][1:2,1])
    # }),
    # y = sapply(grid, function(y) {
    #   mean(y[[1]][2:3,2])
    # }),
    # label = round(gridArea / 1e6, 3)
    # )
    points(x=c(point1[1], point2[1]),
           y=c(point1[2], point2[2]),
           col = 'black', pch=16)
    lines(x=c(point1[1], point2[1]),
          y=c(point1[2], point2[2]),
          col = 'black', lwd=2)
  }
  gridArea
}

areaSparse <- function(point1, point2, grid, n=10, d=4e3, plot = FALSE) {
  length <- distGeo(point1, point2)
  # cat('Length :', length, '\n')
  points <- rbind(makeDetLines(point1, point2, n, d, direction = 90)[(n+1):1,],
                  makeDetLines(point1, point2, n, d, direction = -90)[-1,])
  gridMat <- matrix(unlist(grid), nrow=5)
  if(is.null(attr(grid, 'map'))) {
    grid <- connectFast(grid)
  }
  gridMap <- attr(grid, 'map')
  # gridLength <- apply(points, 1, function(x) gridPropSearch2(x, gridMat, gridMap)) * length
  # browser()
  gridLength <- apply(points, 1, function(x) gridPropSearch2(x, gridMat, gridMap))
  gridLength <- as(gridLength, 'dgCMatrix') * length
  probs <- rep(1, n+1) # function of d, 0 to n. H
  # probs <- seq(from=1, to=.5, length.out = n+1)
  probs <- c(rev(probs), probs[-1])

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
  gridArea
}

ptsToArea <- function(points, grid, n=10, d=4e3, plot=FALSE) {
  gridMat <- matrix(unlist(grid), nrow=5)
  if(is.null(attr(grid, 'map'))) {
    grid <- connectFast(grid)
  }
  gridMap <- attr(grid, 'map')
  result <- sapply(points, function(x) {
    justArea(x[1, ], x[2, ], grid=grid, gridMap=gridMap, n=n, d=d, plot=plot)
  })
  rowSums(result)
}
ptsFor <- function(points, grid, n=10, d=4e3, plot=FALSE, pb=TRUE) {
  # gridMat <- matrix(unlist(grid), nrow=5)
  if(is.null(attr(grid, 'map'))) {
    grid <- connectFast(grid)
  }
  gridMap <- attr(grid, 'map')
  area <- rep(0, length(grid))
  pb <- txtProgressBar(min=1, max=length(points), style=3)
  for(i in seq_along(points)) {
    # cat('\rNum', i)
    setTxtProgressBar(pb, value=i)
    area <- area + justArea(points[[i]][1, ], points[[i]][2, ], grid=grid, gridMap=gridMap, n=n, d=d, plot=plot)
  }
  area
}

justArea <- function(point1, point2,grid, gridMap, n=10, d=4e3, plot = FALSE) {
  length <- distGeo(point1, point2)
  # cat('Length :', length, '\n')
  points <- rbind(makeDetLines(point1, point2, n, d, direction = 90)[(n+1):1,],
                  makeDetLines(point1, point2, n, d, direction = -90)[-1,])
  gridLength <- apply(points, 1, function(x) gridPropSearch3(x, grid, gridMap))
  gridLength <- as(gridLength, 'dgCMatrix') * length
  probs <- rep(1, n+1) # function of d, 0 to n. H
  # probs <- seq(from=1, to=.5, length.out = n+1)
  probs <- c(rev(probs), probs[-1])
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
  gridArea
}

searchPoint2 <- function(point, grid, gridMap) {
  i <- 0
  thisGrid <- 1
  for(i in  1:nrow(gridMap)) {
    if(thisGrid == -1) {
      break
    }
    # thisBox <- gridMat[, (2*thisGrid - 1):(2*thisGrid)]
    thisBox <- .subset2(.subset2(grid, thisGrid), 1)
    if(point[1] > thisBox[2, 1]) {
      thisGrid <- gridMap[thisGrid, 2]
      next
    }
    if(point[1] < thisBox[1,1]) {
      thisGrid <- gridMap[thisGrid, 4]
      next
    }
    if(point[2] > thisBox[3, 2]) {
      thisGrid <- gridMap[thisGrid, 3]
      next
    }
    if(point[2] < thisBox[1, 2]) {
      thisGrid <- gridMap[thisGrid, 1]
      next
    }
    return(thisGrid)
  }
  -1
}

gridPropSearch3 <- function(points, grid,  gridMap) {
  # If we do a search thing, need to check if one/both points are off grid
  # broken if both points are off grid - this doesnt mean nothing exists,
  # line can still intersect
  # browser()
  thisGrid <- searchPoint2(c(points[1], points[2]), grid, gridMap)
  if(thisGrid == -1) {
    points <- points[c(3,4,1,2)]
    thisGrid <- searchPoint2(c(points[1], points[2]), grid, gridMap)
  }
  # gridMat is lat/long so num grids * 2 is ncol
  props <- vector('integer', length(grid))
  if(thisGrid == -1) {
    # if both start and end are off the grid, cant do search method
    cat('\nYOURE OFF THE GRID ITS SCARY OUT HERE. AND ALSO SLOW.\n')
    return(gridPropMat(points, grid))
  }
  toGo <- 1
  i <- 1
  thisSide <- 5
  lastProp <- 0
  for(i in 1:1e3) {
    # thisBox <- gridMat[, (2*thisGrid-1):(2*thisGrid)]
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
