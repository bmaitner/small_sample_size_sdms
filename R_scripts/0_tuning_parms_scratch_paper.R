#Tuning parms


library(disdat)
library(tsutils)
  data_i <- disData(region = "NSW")
  data_i$env <- data_i$env[which(!colnames(data_i$env) %in% c("disturb","soilfert","vegsys"))] 
  data_i$bg <- data_i$bg[which(!colnames(data_i$bg) %in% c("disturb","soilfert","vegsys"))]
  data_i$po <- data_i$po[which(!colnames(data_i$po) %in% c("disturb","soilfert","vegsys"))]
  epsg <- 4326
  

  species_s <- unique(data_i$po$spid)[9]
  group_s <- unique(data_i$po$group[which(data_i$po$spid == species_s)])    
  presence_s <- data_i$po[which(data_i$po$spid == species_s),]
  background_s <- data_i$bg
  
  p.f <- c(rep(1,nrow(presence_s)),rep(0,nrow(background_s)))
  data.f <- rbind(presence_s, background_s)[7:ncol(presence_s)]
  
  data.f <- data.f[1:1000,]
  p.f <- p.f[1:1000]  
  
  cv.
  test <- .cv.maxnet.cm(p.f = p.f,
                data.f = data.f)
  test
  
  
  pbsdm:::dr_maxnet()
  library(maxnet)
  
  test
  test$lambda #199 lambdas?
  test$cvm #mean of cross validated ...something?
  test$cvsd #sd of cross validated ...something?
  test$cvup #upper of cross validated ...something?
  test$cvlo #lower of cross validated ...something?
  test$glmnet.fit
  test$betas #something relating to hinge features?
  test$alpha #na?
  ?maxnet
  
  
  
##################################################


# Modified version of Cory's code below

# should add an option to force a particular lambda rule

#' @param p.f a vector of 1s (for presence points) and 0s (for background points)
#' @param data.f a matrix or data.frame of predictor variables.  Same length at p.f
#' @param form.f = a formula to determine the features to be used. Defaults to maxnet default
#' @param regmult = regularization multiplier. A constant to adjust regularization
#' @param regfun = Regularization function. A function to computer= regularization constant for each feature
#' @param myweights =  
#' @param offset.f =
#' @param maxit Maximum iterations before giving up. Default is 2e5
#' @param Standardize =
#' @param nfolds = 
#' @param foldid =
#' @param parallel =
#' @param verbose Should messages and whatnot be returned? Default is TRUE
#' @param retunPreds =
#' @param quadraticUpperLimit0 =
#' @param lambdaSeq =
#' @export
.cv.maxnet.cm=function(p.f,
                       data.f,
                       form.f = maxnet::maxnet.formula(p.f, data.f),
                       regmult = 1,
                       regfun = maxnet::maxnet.default.regularization,
                       myweights = p.f + (1 - p.f) * 100, #The weights are set this way to make this equivalent to a PPM, per Renner 2015
                       offset.f = NULL, # CM added offset as an arg
                       maxit = 2e5, # CM added larger maxit
                       standardize = TRUE, # added arg for maxent default
                       nfolds = 10,
                       foldid = NULL,
                       parallel = FALSE,
                       verbose = TRUE,
                       returnPreds = FALSE,
                       quadraticUpperLimit0 = FALSE,
                       lambda = NULL, #Lambda Sequence.  If NULL, will use the 
                       ...) {
  #  for testing
  # p.f=p; data.f=m.data; form.f=form; regmult = 1; regfun =  maxnet::maxnet.default.regularization; myweights=p.f + (1 - p.f) * 100; offset.f=full.off; nfolds=10; maxit=1e6; verbose=F; parallel=FALSE; standardize=!toDo$domain$standardizePredictors; returnPreds=F; lambdaSeq=toDo$fitting$lambdaSeq; quadraticUpperLimit0=T
  
  #This bit tosses any points with NA data
  keep <- complete.cases(data.f)
  data.f <- data.f[keep,]
  p.f <- p.f[keep]
  fuck <- which(keep) #apparently this is how Cory keeps track of complete cases
  myweights <- myweights[fuck]
  offset.f <- offset.f[fuck]
  
  mm <- model.matrix(form.f, data.f) #creates a matrix of where rows correspond to env data (data.f) and columns correspond to feature classes
  reg <- regfun(p = p.f, mm) * regmult #creates a regularization constant for each feature, and multiplies them by the reg constant
  glmnet::glmnet.control(pmin = 1e-08, fdev = 1e-5) #sets internal glmnet parameters for a few features
  # set upper.limits on quadratic coeffs to force them to be negative
  
  if(quadraticUpperLimit0){
    cmquadraticUpperLimit0=function(m){
      up=rep(Inf,ncol(m))
      isquadratic <- function(x) grepl("^I\\(.*\\^2\\)", x)
      to0=isquadratic(colnames(m))
      up[to0]=0
      up
    }
    upper.limits=cmquadraticUpperLimit0(mm)
  } else {
    upper.limits=rep(Inf,ncol(mm))
  }
  
  
  if(verbose) print('starting cv model')
  model.cv <- cv.glmnet.cm(x = mm, 
                        y = p.f, 
                        family = "binomial",
                        standardize = standardize, 
                        penalty.factor = reg, 
                        nlambda = 200,
                        lambda = lambda,
                        weights = myweights, 
                        offset=offset.f,
                        maxit=maxit, 
                        nfolds=nfolds,
                        foldRule='PPM', # cm add
                        parallel=FALSE,
                        upper.limits=upper.limits,
                        ...)

    ##Code below is for testing so I don't have to delete the ... above and remember to put it back 
  # 
  # model.cv <- cv.glmnet.cm(x = mm,
  #                          y = p.f,
  #                          family = "binomial",
  #                          standardize = standardize,
  #                          penalty.factor = reg,
  #                          nlambda = 200,
  #                          lambda = lambda,
  #                          weights = myweights,
  #                          offset=offset.f,
  #                          maxit=maxit,
  #                          nfolds=nfolds,
  #                          foldRule='PPM', # cm add
  #                          parallel=FALSE,
  #                          upper.limits=upper.limits)
  # 
  
  print(model.cv)
  class(model.cv) <- c(class(model.cv),"maxnet")
  #if (length(model$beta) < 200 & !length(lambda)==1) # should this be ncol(model$beta)?
  #stop("Error: glmnet failed to complete regularization path")
  
  if(!model.cv$lambda.1se==model.cv$lambda[1])	{			
    lambda=model.cv$lambda.1se
    model.cv$lambdaRule='1se'
  } else if (!model.cv$lambda.min==model.cv$lambda[1]){
    lambda=model.cv$lambda.min
    model.cv$lambdaRule='min'
  } else {
    lambda=tail(model.cv$lambda,1)
    model.cv$lambdaRule='last'
  }
  model.cv$lambdaKeepIndex=which(model.cv$lambda==lambda)
  print('5.4')
  print(' ')
  print(' ')
  #CM: commented on 6/14/18 because removing zeros moved to prediction to avoid losing variable names in the model object
  # CM: 7/24/18: uncommented because we need the betas to be defined 
  
  #-- figure out which lambda to use 
  #if(length(lambda)==1){
  #	bb <- model$beta
  #} else { bb <- model$beta[, keep] }
  bb <- model.cv$glmnet$beta[, model.cv$lambdaKeepIndex]
  model.cv$betas <- bb[bb != 0]   # this lead to some errors by duplicated row names (as of 7/24/18 we haven't gone back to resolve this) # 3/22/21 although its only an issue when hinges are used, and they aren't now
  #cat('model$betas ', length(model$betas),'\n')
  model.cv$alpha <- 0
  # CM: 7/24/28: Annoyingly, you need to predict the whole distribution to get the normalization constant and it seems (but shouldn't be) that estimating the constant in the projection stage deosn't work
  if(returnPreds){
    rr=.predict.maxnetCM2(object=model,
                          newdata=data.frame(data[p.f == 0, , drop = FALSE]),
                          type = "exponent",
                          clamp = F, offset=offset[p.f == 0, drop = FALSE]) # changed to my general predict function with offset
    raw <- rr/sum(rr)
    model.cv$entropy <- -sum(raw * log(raw))
    model.cv$alpha <- -log(sum(rr))
  } else {
    model.cv$entropy=NA
    model.cv$alpha=NA
    raw=NULL
  }
  print('5.5')
  
  
  model.cv$penalty.factor <- reg
  model.cv$featuremins <- apply(mm, 2, min)
  model.cv$featuremaxs <- apply(mm, 2, max)
  vv <- (sapply(data.f, class) != "factor")
  model.cv$varmin <- apply(data.f[, vv, drop = FALSE], 2, min)
  model.cv$varmax <- apply(data.f[, vv, drop = FALSE], 2, max)
  means <- apply(data.f[p.f == 1, vv, drop = FALSE], 2, mean)
  majorities <- sapply(names(data.f)[!vv], function(n) which.max(table(data.f[p.f ==
                                                                                1, n, drop = FALSE])))
  names(majorities) <- names(data.f)[!vv]
  model.cv$samplemeans <- c(means, majorities)
  model.cv$levels <- lapply(data.f, levels)
  print('5.6')
  
  
  model.cv
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.cv.maxnet.cm

#Possibly unneeded helper function to set lambaseq to null  
lambdaSeq <- function(reg,p.f,weights){NULL}  



#' @param x input matrix of covariates
#' @param y response variable (0 and 1 for presence-background modelins)
#' @param weights observation weights (same length as y)
#' @param offset vector of same length as x.  Default is NULL
#' @param lambda lambda sequence
#' @param type.measure 
#' @param nfolds
#' @param foldid
#' @param grouped
#' @param keep
#' @param parallel
#' @param foldrule
#' @export
cv.glmnet.cm <- function (x,
                          y,
                          weights = NULL,
                          offset = NULL,
                          lambda = NULL,
                          type.measure = c("mse", "deviance", "class", "auc", "mae"),
                          nfolds = 10,
                          foldid = NULL,
                          grouped = TRUE,
                          keep = FALSE,
                          parallel = FALSE,
                          foldRule = 'random',...) {
  
  #  For testing
  #  x = mm; y = p; family = "binomial"; penalty.factor = reg; weights = myweights; lambda = 10^(seq(4,0, length.out = 200)) * sum(reg)/length(reg) * sum(p)/sum(weights); foldRule='PPM'; type.measure = "default";   nfolds = 10;  grouped = TRUE;  keep = FALSE;  parallel = FALSE
  
  type.measure = match.arg(type.measure) ##
  
  # #     if (missing(type.measure)) {
  # #       type.measure = "default"
  # #     } else {
  # #     	type.measure = match.arg(type.measure)
  # #     }
  #9/22/21 this might be needed...
  #type.measure = match.arg(arg = type.measure,choices=type.measure,several.ok = TRUE)
  
  if (!is.null(lambda) && length(lambda) < 2)
    stop("Need more than one value of lambda for cv.glmnet")
  N = nrow(x)
  
  if (is.null(weights)) {
    weights = rep(1, N)
  } else {weights = as.double(weights)}
  y = drop(y)
  
  cv.call = glmnet.call = match.call(expand.dots = TRUE)
  which.toss = match(c("type.measure", "nfolds", "foldid", "grouped",
                       "keep"), names(glmnet.call), FALSE) ##
  which.toss[is.na(which.toss)]=0
  ## which.toss = match(c("type.measure", "nfolds", "foldid", "grouped","keep"), names(glmnet.call), F)
  print('5.1.5')
  print(which.toss)
  print(glmnet.call)
  if (any(which.toss) & length(which.toss)>0)  glmnet.call = glmnet.call[-which.toss]
  glmnet.call[[1]] = as.name("glmnet")
  print('5.1.6')
  
  glmnet.object = glmnet(x, y, weights = weights, offset = offset, lambda = lambda, ...)
  print('5.1.7')
  
  glmnet.object$call = glmnet.call
  subclass = class(glmnet.object)[[1]]
  type.measure = glmnet:::cvtype(type.measure, subclass)
  is.offset = glmnet.object$offset
  print('5.1.8')
  
  if (inherits(glmnet.object, "multnet") && !glmnet.object$grouped) {
    nz = predict(glmnet.object, type = "nonzero")
    nz = sapply(nz, function(x) sapply(x, length))
    nz = ceiling(apply(nz, 1, median))
  } else nz = sapply(predict(glmnet.object, type = "nonzero"),
                     length)
  print('5.1.9')
  
  if (missing(foldid)) {
    if(foldRule=='PPM'){
      x.pres=which(y==1)
      # because the code tosses a fold, this i the only thig you need to specify (and maybe could just specify the folds of the bg as NA
      foldid=rep(-1000,length(y)) # dummy value that get's used for background
      # just assign folds to the presences
      foldid[y==1]=sample(rep(seq(nfolds), length = length(x.pres)))
    } else { foldid = sample(rep(seq(nfolds), length = N))}
  } else {nfolds = max(foldid)}
  print('5.1.10')
  
  if (nfolds < 3)
    stop("nfolds must be bigger than 3; nfolds=10 recommended")
  outlist = as.list(seq(nfolds))
  if (parallel) {
    outlist = foreach(i = seq(nfolds), .packages = c("glmnet")) %dopar%{
      which = foldid == i
      if (length(dim(y)) > 1) {
        y_sub = y[!which, ]
      } else y_sub = y[!which]
      if (is.offset) {
        offset_sub = as.matrix(offset)[!which, ]
      } else offset_sub = NULL
      x_sub=x[!which, , drop = FALSE]
      glmnet(x_sub, y_sub, lambda = lambda,
             offset = offset_sub, weights = weights[!which],
             ...)
    }
  } else {
    for (i in seq(nfolds)) {
      which = foldid == i
      if (is.matrix(y)) y_sub = y[!which, ] else y_sub = y[!which]
      if (is.offset)  offset_sub = as.matrix(offset)[!which, ] else offset_sub = NULL
      x_sub=x[!which, , drop = FALSE]
      outlist[[i]] = glmnet::glmnet(x_sub, y_sub, lambda = lambda, offset = offset_sub,weights = weights[!which], ...)
    }
  }
  fun = paste("cv", subclass,'cm', sep = ".")
  lambda = glmnet.object$lambda
  
  cvstuff = do.call(fun, list(outlist, lambda, x, y, weights[y==1], offset[y==1], foldid[y==1], type.measure, grouped, keep))
  
  cvm = cvstuff$cvm
  cvsd = cvstuff$cvsd
  nas = is.na(cvsd)
  if (any(nas)) {
    lambda = lambda[!nas]
    cvm = cvm[!nas]
    cvsd = cvsd[!nas]
    nz = nz[!nas]
  }
  cvname = names(cvstuff$type.measure)
  names(cvname) = cvstuff$type.measure
  out = list(lambda = lambda, cvm = cvm, cvsd = cvsd, cvup = cvm +
               cvsd, cvlo = cvm - cvsd, nzero = nz, name = cvname, glmnet.fit = glmnet.object)
  if (keep)
    out = c(out, list(fit.preval = cvstuff$fit.preval, foldid = foldid))
  lamin = if (cvname == "AUC") { getmin(lambda, -cvm, cvsd) } else getmin(lambda, cvm, cvsd)
  obj = c(out, as.list(lamin))
  class(obj) = "cv.glmnet"
  obj
}


##################################

#' @param outlist
#' @param lambda
#' @param x
#' @param y
#' @param weights
#' @param offset
#' @param foldid
#' @param type.measure
#' @param grouped
#' @param keep
#' @export
cv.lognet.cm <- function (outlist,
                          lambda,
                          x,
                          y,
                          weights,
                          offset,
                          foldid,
                          type.measure,
                          grouped,
                          keep = FALSE) {
  
  prob_min = 1e-05
  prob_max = 1 - prob_min
  nc = dim(y)
  if (is.null(nc)) {
    y = as.factor(y)
    ntab = table(y)
    nc = as.integer(length(ntab))
    y = diag(nc)[as.numeric(y), ]
  }
  N = nrow(y)
  nfolds = max(foldid)
  if ((N/nfolds < 10) && type.measure == "auc") {
    warning("Too few (< 10) observations per fold for type.measure='auc' in cv.lognet; changed to type.measure='deviance'. Alternatively, use smaller value for nfolds",
            call. = FALSE)
    type.measure = cvtype("deviance", "lognet")
  }
  if ((N/nfolds < 3) && grouped) {
    warning("Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold",
            call. = FALSE)
    grouped = FALSE
  }
  if (!is.null(offset)) {
    is.offset = TRUE
    offset = drop(offset)
  }
  else is.offset = FALSE
  mlami = max(sapply(outlist, function(obj) min(obj$lambda)))
  which_lam = lambda >= mlami
  predmat = matrix(NA, nrow(y), length(lambda))
  nlams = double(nfolds)
  for (i in seq(nfolds)) {
    which = foldid == i
    fitobj = outlist[[i]]
    if (is.offset)
      off_sub = offset[which]
    preds = predict(fitobj, x[which, , drop = FALSE], s = lambda[which_lam],
                    newoffset = off_sub, type = "response")
    nlami = sum(which_lam)
    predmat[which, seq(nlami)] = preds
    nlams[i] = nlami
  }
  if (type.measure == "auc") {
    cvraw = matrix(NA, nfolds, length(lambda))
    good = matrix(0, nfolds, length(lambda))
    for (i in seq(nfolds)) {
      good[i, seq(nlams[i])] = 1
      which = foldid == i
      for (j in seq(nlams[i])) {
        cvraw[i, j] = auc.mat(y[which, ], predmat[which,
                                                  j], weights[which])
      }
    }
    N = apply(good, 2, sum)
    weights = tapply(weights, foldid, sum)
  }
  else {
    ywt = apply(y, 1, sum)
    y = y/ywt
    #weights = weights * ywt
    N = nrow(y) - apply(is.na(predmat), 2, sum)
    cvraw = switch(type.measure, mse = (y[, 1] - (1 - predmat))^2 +
                     (y[, 2] - predmat)^2, mae = abs(y[, 1] - (1 - predmat)) +
                     abs(y[, 2] - predmat), deviance = {
                       predmat = pmin(pmax(predmat, prob_min), prob_max)
                       lp = y[, 1] * log(1 - predmat) + y[, 2] * log(predmat)
                       ly = log(y)
                       ly[y == 0] = 0
                       ly = drop((y * ly) %*% c(1, 1))
                       2 * (ly - lp)
                     }, class = y[, 1] * (predmat > 0.5) + y[, 2] * (predmat <=
                                                                       0.5))
    if (grouped) {
      cvob = cvcompute.cm(cvraw, weights, foldid, nlams)
      cvraw = cvob$cvraw
      #weights = cvob$weights
      N = cvob$N
    }
  }
  #cvm = apply(cvraw, 2, weighted.mean, w = weights, na.rm = TRUE)
  cvm = apply(cvraw, 2, mean, na.rm = TRUE)
  
  #cvsd = sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, weighted.mean,
  #    w = weights, na.rm = TRUE)/(N - 1))
  cvsd = sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean,  na.rm = TRUE)/(N - 1))
  out = list(cvm = cvm, cvsd = cvsd, type.measure = type.measure)
  if (keep)
    out$fit.preval = predmat
  out
}

###############################

#' @param mat
#' @param weights
#' @param foldid
#' @param nlams
#' @export
# had to turn off weighted means
cvcompute.cm=function(mat, weights, foldid, nlams) {
  #wisum = tapply(weights, foldid, sum)
  nfolds = max(foldid)
  outmat = matrix(NA, nfolds, ncol(mat))
  good = matrix(0, nfolds, ncol(mat))
  mat[is.infinite(mat)] = NA
  for (i in seq(nfolds)) {
    mati = mat[foldid == i, , drop = FALSE]
    #wi = weights[foldid == i]
    #outmat[i, ] = apply(mati, 2, weighted.mean, w = wi, na.rm = TRUE)
    outmat[i, ] = apply(mati, 2, mean, na.rm = TRUE)
    good[i, seq(nlams[i])] = 1
  }
  N = apply(good, 2, sum)
  list(cvraw = outmat, N = N) #weights = wisum, N = N)
}


################

#' @param lambda
#' @param cvm
#' @param cvsd
#' @export
getmin=function(lambda, cvm, cvsd){
  cvmin=min(cvm,na.rm=TRUE)
  idmin=cvm<=cvmin
  lambda.min=max(lambda[idmin],na.rm=TRUE)
  idmin=match(lambda.min,lambda)
  semin=(cvm+cvsd)[idmin]
  idmin=cvm<=semin
  lambda.1se=max(lambda[idmin],na.rm=TRUE)
  list(lambda.min=lambda.min,lambda.1se=lambda.1se)
}

