#' @importFrom geosphere distGeo destPoint bearing
#'
doGridEffort <- function(point1, point2, grid, n=10, d=4e3, detFun=NULL, plot = FALSE) {
    length <- distGeo(point1, point2)
    points <- rbind(makeDetLines(point1, point2, n, d, direction = 90)[(n+1):1,],
                    makeDetLines(point1, point2, n, d, direction = -90)[-1,])
    gridMat <- matrix(unlist(grid), nrow=5)
    if(is.null(attr(grid, 'map'))) {
        grid <- connectFast(grid)
    }
    gridMap <- attr(grid, 'map')
    gridLength <- apply(points, 1, function(x) gridPropSearch(x, gridMat, gridMap)) * length
    if(is.null(detFun)) {
        detFun <- 1
    }
    if(is.numeric(detFun)) {
        if(detFun > 1) detFun <- 1
        if(detFun < 0) detFun <- 0
        detFun <- function(d) rep(detFun, length(d))
    }
    if(inherits(detFun, 'ds')) {
        # model from Distance, nonsense doesnt work yet wtffff
    }
    dist <- seq(from=0, to=d, length.out = n+1)
    # probs <- rep(1, n+1) # function of d, 0 to n. H
    probs <- detFun(dist)
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
        plot(grid)
        for(i in 1:nrow(points)) {
            lines(x=points[i, c(1, 3)], y=points[i, c(2, 4)], col = 'blue')
        }
        text(x = sapply(grid, function(x) {
            mean(x[[1]][1:2,1])
        }),
        y = sapply(grid, function(y) {
            mean(y[[1]][2:3,2])
        }),
        label = round(gridArea / 1e6, 3)
        )
        points(x=c(point1[1], point2[1]),
               y=c(point1[2], point2[2]),
               col = 'black', pch=16)
        lines(x=c(point1[1], point2[1]),
              y=c(point1[2], point2[2]),
              col = 'black', lwd=2)
    }
    gridArea
}

gridPropSearch <- function(points, gridMat, gridMap) {
    # If we do a search thing, need to check if one/both points are off grid
    # broken if both points are off grid - this doesnt mean nothing exists,
    # line can still intersect
    thisGrid <- searchPoint(c(points[1], points[2]), gridMat, gridMap)
    if(thisGrid == -1) {
        points <- points[c(3,4,1,2)]
        thisGrid <- searchPoint(c(points[1], points[2]), gridMat, gridMap)
    }
    props <- rep(0, ncol(gridMat)/2)
    if(thisGrid == -1) {
        # if both start and end are off the grid, cant do search method
        return(gridPropMat(points, gridMat))
    }
    toGo <- 1
    i <- 1
    thisSide <- 5
    lastProp <- 0
    while(i < 1e3) {
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
        i <- i+1
    }
    props
}

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

connectFast <- function(grid) {
    gridMat <- matrix(unlist(grid), nrow = 5)
    conMat <- matrix(NA, ncol = 4, nrow = length(grid))
    for(i in seq_along(grid)) {
        for(side in 1:4) {
            if(!is.na(conMat[i, side])) {
                next
            }
            if(i==1) {
                matched <- edgeFast(i, side, gridMat, start=i)
            } else{
                matched <- edgeFast(i, side, gridMat[, -c(1:(2*(i-1)))], start=i)
            }
            conMat[i, side] <- matched
            if(matched == -1) {
                next
            }
            conMat[matched, (side+1) %% 4 + 1] <- i
        }
    }
    class(grid) <- c(class(grid), 'connectedGrid')
    attr(grid, 'map') <- conMat
    grid
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

searchPoint <- function(point, gridMat, gridMap) {
    i <- 0
    thisGrid <- 1
    while(i < nrow(gridMap)) {
        if(thisGrid == -1) {
            break
        }
        thisBox <- gridMat[1:3, (2*thisGrid - 1):(2*thisGrid)]
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

edgeFast <- function(id, edgeNo, gridMat, start=2) {
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