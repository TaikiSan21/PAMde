

downloadEnv <- function(data, env, fileName = NULL, buffer = c(0, 0, 0)) {
    dataBounds <- dataToBounds(data, buffer)
    if(is.list(env)) {
        dataBounds$Longitude <- to180(dataBounds$Longitude, inverse = !env$is180)
        url <- datalistToURL(env, ranges=dataBounds)
    } else if(is.character(env)) {
        # do some erddap info checking shit and make a URL for it
        # list above needs base, vars, dataset, fileType, source
        info <- suppressMessages(try(rerddap::info(env)))
        if(inherits(info, 'try-error')) {
            stop('Not a valid erddap dataset')
        }
        env <- erddapToDatalist(info)
        return(downloadEnv(data, env, fileName, buffer))
    }
    # hm shit need some tempdir stuff, either make one in wd with a weird name, or tempdir(),
    # or steal rerddap:::rrcache$cache_path_get() lol
    if(is.null(fileName)) {
        tempDir <- rerddap:::rrcache$cache_path_get()
        if(!dir.exists(tempDir)) {
            dir.create(tempDir, recursive = TRUE)
        }
        fileName <- paste0(tempDir, '/TEMPFILE.nc')
    }

    # FOR FILE NAMES LETS DO "DATASET NAME_NUMBER"
    maxTries <- 3
    nTry <- 1
    while(nTry <= maxTries) {
        envData <- try(suppressMessages(httr::GET(url, verbose(), progress(), write_disk(fileName, overwrite = TRUE))))
        if(inherits(envData, 'try-error')) {
            nTry <- nTry + 1
            next
        } else {
            return(invisible(fileName))
        }
    }
    FALSE
}

ncToData <- function(data, ncFile, buffer = c(0,0,0), quiet=FALSE) {
    nc <- nc_open(ncFile)
    on.exit(nc_close(nc))
    # this finds the name of the longitude coordinate from my list that matches lon/long/whatever to Longitude
    nc180 <- ncIs180(nc)
    data180 <- dataIs180(data)
    if(nc180 != data180) {
        data <- to180(data, inverse = !nc180)
    }
    # for each variable, make the ncvar_get call which needs start and count#
    # these are XYZT if Z is present, -1 count means get all
    # OPTION TO ONLY GET A CERTAIN Z VALUE LIKE DEPTH FIRST VALUE SOMEHOW
    varNames <- names(nc$var)
    allVar <- vector('list', length = length(varNames))
    names(allVar) <- varNames
    for(v in seq_along(allVar)) {
        allVar[[v]] <- vector('list', length = nrow(data))
    }
    for(i in 1:nrow(data)) {
        for(v in varNames) {
            allVar[[v]][[i]] <- getVarData(data[i,], nc=nc, var=v, buffer = buffer, quiet=quiet)
        }
    }
    allVar <- bind_rows(
        lapply(purrr::transpose(allVar), function(i) {
            result <- vector('list', length = 3 * length(i))
            names(result) <- paste0(rep(names(i), each = 3),
                                    c('_mean', '_median', '_stdev'))
            for(v in seq_along(i)) {
                result[(v-1)*3 + 1] <- mean(i[[v]])
                result[(v-1)*3 + 2] <- median(i[[v]])
                result[(v-1)*3 + 3] <- sd(i[[v]])
            }
            bind_cols(result)
        })
    )
    cbind(data, allVar)
}



getVarData <- function(data, nc, var, buffer, quiet=FALSE) {
    names(nc$dim) <- standardCoordNames(names(nc$dim))
    thisVar <- nc[['var']][[var]]
    names(thisVar$dim) <- names(nc$dim)[thisVar$dimids + 1]
    xIx <- dimToIx(data$Longitude, thisVar$dim$Longitude, buffer[1], quiet)
    yIx <- dimToIx(data$Latitude, thisVar$dim$Latitude, buffer[2], quiet)
    tIx <- dimToIx(data$UTC, thisVar$dim$UTC, buffer[3], quiet)
    hasZ <- any(!(names(thisVar$dim) %in% c('Longitude', 'Latitude', 'UTC')))
    if(hasZ) {
        start <- c(xIx$start, yIx$start, 1, tIx$start)
        count <- c(xIx$count, yIx$count, -1, tIx$count)
    } else {
        start <- c(xIx$start, yIx$start, tIx$start)
        count <- c(xIx$count, yIx$count, tIx$count)
    }
    ncvar_get(nc, varid = var, start = start, count = count)
}
