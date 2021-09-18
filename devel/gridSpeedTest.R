# need for run fast
# all the way through
source('./devel/taikiTest.R')
gps <- read.csv('./devel/straights/straightPathWeffort_1303.csv', stringsAsFactors = FALSE)
pmDets <- read.csv('./devel/yberddap/Spermies_20200310.csv', stringsAsFactors = FALSE)
pmDets <- rename(pmDets, Longitude = lon, Latitude=lat)
pmDets$Longitude <- ifelse(pmDets$Longitude <= 0, pmDets$Longitude + 360, pmDets$Longitude)


detDistance <- 20
boundary <- dfToBounds(gps) # need all to be same range, say 0-360
pix <- .09
grid <- makeGrid(boundary, pixel = pix, buffer_km = detDistance, plot = T) # did with .09
conGrid <- connectFast(grid) # this takes long, add prog bar?

profvis(connectFast(grid))
profvis(connectGrid(grid))
microbenchmark(
    old = connectGrid(grid, pbshow = TRUE),
    new = connectGrid(grid),
    times=5)
ends <- getEndPoints(gps)
plot(conGrid)


library(profvis)
profvis({for(i in 1:1000) {
    areaFast(ends[[i]][1,], ends[[i]][2,], grid=conGrid, plot=FALSE)
}}) #10160 for 100 4/24 OLD, 6290 rowSums. 4330 sparse
profvis(ptsFor(ends[1:100], grid = conGrid)) # 1880 do matrix part once

library(microbenchmark)

gn <- 500
gridDt <- data.table(matrix(unlist(grid), nrow=5))
microbenchmark(
    # rep = rep(0L, 5e5),
    dt = gridDt[, (2*gn-1):(2*gn)],
    norm = gridMat[, (2*gn-1):(2*gn)],
    listy = conGrid[[gn]][[1]],
    subs = .subset2(.subset2(conGrid, gn), 1),
    times = 5000)

identical(gridMat[, (2*gn-1):(2*gn)],.subset2(.subset2(conGrid, gn), 1))
microbenchmark(
    oldA = doAreaShit2(ends[[1]][1,], ends[[1]][2,], grid=conGrid, plot=FALSE), # 110ms 4/24
    newA = areaFast(ends[[1]][1,], ends[[1]][2,], grid=conGrid, plot=FALSE), # 53ms 4/24
    sparseA = areaSparse(ends[[1]][1,], ends[[1]][2,], grid=conGrid, plot=FALSE), # 44ms 4/27
    sparseFast = ptsFor(ends[1], conGrid, plot=FALSE), # 15ms 4/28
    times=100) #

myDsPredict <- function(model, newdata) {
    # detfct will use the design matrix in scale and shape to get covariates
    # so build a new one using your new data, all factors and shit must be
    # code the same or you can go fys
    if(inherits(model, 'dsmodel')) {
        model <- model$ddf
    }
    if(inherits(model, 'ds')) {
        model <- model$ds$aux$ddfobj
    }
    if(!all(c('scale', 'shape', 'xmat') %in% names(model))) {
        stop('model input does not appear to be in the right format, should be of class dsmodel, ds, or ddf')
    }
    browser()
    model$scale$dm <- model.matrix(formula(model$scale$formula), data=newdata)
    model$shape$dm <- model.matrix(formula(model$shape$formula), data=newdata)
    model$xmat <- model$xmat[rep(1, nrow(newdata)), , drop=FALSE]
    detfct(newdata$distance, ddfobj = model)
}

library(Distance)
pmDets$distance <- abs(pmDets$pdist)
withCov <- Distance::ds(pmDets, key='hr', formula = ~ acid + peak, truncation = list(left=1000, right=30e3))
dsm <- Distance::ds(pmDets, key='hr')
plot(withCov)
newDat <- data.frame(distance = c(15000, 15000, 20000, 20000),
                     peak = c('A', 'B', 'A', 'B'),
                     Latitude = c(20, 21, 22, 23),
                     acid = c(150, 150, 250, 250))
myDsPredict(withCov$ddf, newdata=newDat)
add_df_covar_line(dsm, newDat)

hm <- myDsPredict(withCov$ddf, newdata = data.frame(distance = c(15000, 15000, 20000, 20000),
                                                    peak = c('A', 'B', 'A', 'B')))

eval_with_covars(distance = newDat$distance, newdata = newDat, model = dsm)
add_df_covar_line(withCov, newDat)

# if shits less than left, make 0
makeDetFun <- function(covarDf, dsModel) {
    function(distance) {
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

# okay so make all the different detction functions you might have and number them, on the track line
# just store the number of the detfun and we go grab it from a list. these have been formatted so that
# its just p = f(d) to avoid worrying about inputting covariates in the "calcArea" step
gps$aeffort <- ifelse(gps$aeffort == 'off', FALSE, TRUE)
gps$overallEffort <- gps$straight & gps$aeffort
gps$effort <- gps$overallEffort
library(sf)
library(Matrix)
system.time(eff <- doAllGrid(gps=gps, dets = pmDets, trunc=20e3,dsmodel = dsm, pixel = .225))



boundary <- dfToBounds(rbind(gps[, c('Longitude', 'Latitude')],
                             pmDets[, c('Longitude', 'Latitude')])) # need all to be same range, say 0-


grid <- makeGrid(boundary, buffer_km = 20, plot = FALSE) # did with .09


conGrid <- connectGrid(grid) # this takes long, add prog bar?

# this needs "effort" column, so we do this first
ends <- getEndPoints(gps, length = 1e3) # length of segs to break into, less should more accurato
effort <- calcEffort(points=ends, grid=conGrid, n=10, d=trunc, dsmodel=dsmodel)