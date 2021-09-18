# Add detection function lines to a plot.ds() plot!
#
# args:
#  ddf   - a fitted detection function
#  data  - a data.frame with the covariate combination you want to plot
#  ...   - extra arguments to give to line() (like lty or lwd, or col)

# return
#  inivisbly, the values of p over the truncation range

add_df_covar_line <- function(ddf, data, ...){
  # browhmser()
  df <- ddf$ddf

  left <- df$meta.data$left
  width <- df$meta.data$width

  xx <- seq(left, width, length.out=250)
  # xx <- data$distance

  #data <- model.frame(df$ds$aux$ddfobj$scale$formula, data)
  xm <- df$ds$aux$ddfobj$xmat

  data$object <- 1:nrow(data)
  data$distance <- rep(0, nrow(data))
  data$observer <- rep(0, nrow(data))
  data$detected <- rep(0, nrow(data))
  data$binned <- rep(df$ds$aux$ddfobj$xmat$binned[1], nrow(data))

  eval_with_covars <- function(distance, newdata, model){
    ddfobj <- model$ds$aux$ddfobj

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

  linedat <- eval_with_covars(xx, data, df)

  lines(xx, linedat, ...)

  # linedat
  invisible(linedat)
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
