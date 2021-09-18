# tmep not delete
doAreaFor <- function(point1, point2, grid, n=10, d=4e3, plot = FALSE) {
    length <- distGeo(point1, point2)
    points <- rbind(makeDetLines(point1, point2, n, d, direction = 90)[(n+1):1,],
                    makeDetLines(point1, point2, n, d, direction = -90)[-1,])
    gridMat <- matrix(unlist(grid), nrow=5)
    gridLength <- apply(points, 1, function(x) gridPropMat(x, gridMat)) * length
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

makeConnectedGrid <- function(grid) {
    # these match edge number
    dirs <- c('south', 'east', 'north', 'west')
    for(i in seq_along(grid)) {
        for(side in 1:4) {
            # already filled, skip
            if(!is.null(attr(grid[[i]], dirs[side]))) {
                next
            }
            matched <- matchEdge(i, side, grid)
            attr(grid[[i]], dirs[side]) <- matched
            if(matched == -1) {
                next
            }
            # if match was found, set the opposite of that one too
            attr(grid[[matched]], dirs[(side + 1) %% 4 + 1]) <- i
        }
    }
    class(grid) <- c('connectedGrid', class(grid))
    grid
}

matchEdge <- function(id, edgeNo, gridData) {
    # x1, y1, x2, y2
    thisEdge <- .subset2(gridData, id)[[1]][edgeNo:(edgeNo+1), ]
    nextEdgeNo <- (edgeNo + 1) %% 4 + 1
    for(i in seq_along(gridData)) {
        nextEdge <- .subset2(.subset2(gridData, i), 1)[nextEdgeNo:(nextEdgeNo+1), 1:2]
        # if(thisEdge[1, 1] == nextEdge[2, 1] &&
        #    thisEdge[1, 2] == nextEdge[2, 2] &&
        #    thisEdge[2, 1] == nextEdge[1, 1] &&
        #    thisEdge[2, 2] == nextEdge[1, 2]) {
        #   return(i)
        # }
        if(thisEdge[1, 1] == nextEdge[2, 1] &&
           thisEdge[1, 2] == nextEdge[2, 2]) {
            return(i)
        }
    }
    -1
}

matchPoint <- function(point, grid) {
    for(i in seq_along(grid)) {
        thisBox <- .subset2(.subset2(grid, i), 1)
        if(point[1] >= thisBox[1, 1] &&
           point[1] <= thisBox[2, 1] &&
           point[2] >= thisBox[1, 2] &&
           point[2] <= thisBox[3, 2]) {
            return(i)
        }
    }
    -1
}

gridProportion <- function(startPoint, endPoint, grid) {
    # If we do a search thing, need to check if one/both points are off grid
    startGrid <- matchPoint(startPoint, grid)
    endGrid <- matchPoint(endPoint, grid)
    props <- sapply(seq_along(grid), function(g) {
        thisBox <- grid[[g]][[1]]
        thisProp <- sapply(1:4, function(side) {
            getIntersection(startPoint[1], startPoint[2], endPoint[1], endPoint[2],
                            thisBox[side, 1], thisBox[side, 2], thisBox[side+1, 1], thisBox[side+1, 2])
        })
        thisProp <- thisProp[thisProp > 0]
        if(length(thisProp) == 2) {
            return(abs(thisProp[1]-thisProp[2]))
        }
        if(length(thisProp) == 1) {
            if(g == startGrid) {
                return(thisProp)
            }
            if(g == endGrid) {
                return(1 - thisProp)
            }
        }
        0
    })
    props
}

gridProportionFor <- function(startPoint, endPoint, grid) {
    # If we do a search thing, need to check if one/both points are off grid
    startGrid <- matchPoint(startPoint, grid)
    endGrid <- matchPoint(endPoint, grid)
    props <- rep(0, length(grid))
    for(g in seq_along(grid)) {
        thisBox <- grid[[g]][[1]]
        thisProp <- rep(0, 4)
        for(side in 1:4) {
            thisProp[side] <- getIntersection(startPoint[1], startPoint[2], endPoint[1], endPoint[2],
                                              thisBox[side, 1], thisBox[side, 2], thisBox[side+1, 1], thisBox[side+1, 2])
        }
        thisProp <- thisProp[thisProp > 0]
        if(length(thisProp) == 2) {
            props[g] <- abs(thisProp[1]-thisProp[2])
            next
        }
        if(length(thisProp) == 1) {
            if(g == startGrid) {
                props[g] <- thisProp
                next
            }
            if(g == endGrid) {
                props[g] <- 1 - thisProp
                next
            }
        }
    }
    props
}

getIntFast <- function(t1x, t1y, t2x, t2y, g1x, g1y, g2x, g2y) {
    # t  transect g  grid
    s1x <- t2x - t1x
    s1y <- t1y - t2y # this is always negative later
    s2x <- g2x - g1x
    s2y <- g2y - g1y
    # browser()
    denom <- (s2x * s1y + s1x * s2y)
    if(denom < 1e-8 &&
       denom > -1e-8) {
        return(FALSE)
    } else {
        denom <- 1/denom
    }
    tgx <- t1x - g1x
    tgy <- t1y - g1y
    s <- (s1y * tgx + s1x * tgy) * denom
    if(s >= 0 && s <= 1) {
        t <- (s2x * tgy - s2y * tgx) * denom
        if(t >=0 && t <= 1) {
            return(t)
        }
    }
    FALSE
}

getIntHoriz <- function(t1x, t1y, t2x, t2y, g1x, g1y, g2x, g2y) {
    # t  transect g  grid
    s1x <- t2x - t1x
    s1y <- t2y - t1y
    s2x <- g2x - g1x
    # s2y <- g2y - g1y  these are 0
    # browser()
    denom <- -s2x * s1y
    if(abs(denom) < 1e-8) {
        return(FALSE)
    } else {
        denom <- 1/denom
    }
    t <- s2x * (t1y - g1y) * denom
    if(t >=0 && t <= 1) {
        s <- (-s1y * (t1x - g1x) + s1x * (t1y - g1y)) * denom
        if(s >= 0 && s <= 1) {
            return(t)
        }
    }
    FALSE
}

getIntVert <- function(t1x, t1y, t2x, t2y, g1x, g1y, g2x, g2y) {
    # t  transect g  grid
    s1x <- t2x - t1x
    s1y <- t2y - t1y
    # s2x <- g2x - g1x these are 0
    s2y <- g2y - g1y
    # browser()
    denom <- (s1x * s2y)
    if(abs(denom) < 1e-8) {
        return(FALSE)
    } else {
        denom <- 1/denom
    }
    t <- - s2y * (t1x - g1x) * denom
    if(t >=0 && t <= 1) {
        s <- (-s1y * (t1x - g1x) + s1x * (t1y - g1y)) * denom
        if(s >= 0 && s <= 1) {
            return(t)
        }
    }
    FALSE
}