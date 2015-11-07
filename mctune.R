# modified version of "tune" from package "e1071"
# licence: GPL-2
# details: http://cran.r-project.org/web/packages/e1071/index.html
mctune <- function(method, train.x, train.y = NULL, data = list(),
                 validation.x = NULL, validation.y = NULL,
                 ranges = NULL, predict.func = predict,
                 tunecontrol = tune.control(),
                 mc.control=NULL,
                 confusionmatrizes=F,
                 ...
) {
  call <- match.call()
 
  require('plyr')
  require('parallel')
 
  ## internal helper functions
  resp <- function(formula, data) {
    model.response(model.frame(formula, data))
  }
 
  classAgreement2 <- function (tab) {
    n <- sum(tab)
    # correct classification rate
    if (!is.null(dimnames(tab))) {
      lev <- intersect(colnames(tab), rownames(tab))
      d <- diag(tab[lev, lev])
      p0 <- sum(d) / n
    } else {
      m <- min(dim(tab))
      d <- diag(tab[1:m, 1:m])
      p0 <- sum(d) / n
    }
    # confusion matrizes
    if(!confusionmatrizes) {
      list(p0=p0)
    } else if(is.null(dimnames(tab))) {
      stop('tables without dimension names are not allowed when generating confusionmatrizes.')
    }
    else {
      # generate confusion matrix for each class
      classnames <- unique(unlist(dimnames(tab)))
      truepositives <-  unlist(lapply(classnames, function(positiveclassname) { sum(d[positiveclassname])}))
      falsepositives <- unlist(lapply(classnames, function(positiveclassname) { sum(tab[positiveclassname,])-tab[positiveclassname,positiveclassname]}))
      falsenegatives <- unlist(lapply(classnames, function(positiveclassname) { sum(tab[,positiveclassname])-tab[positiveclassname,positiveclassname]}))
      truenegatives <- mapply(FUN=function(tp,fn,fp){ sum(tab)-tp-fn-fp }, truepositives, falsenegatives, falsepositives)
      confusions <- data.frame(classnames, truepositives, truenegatives, falsepositives, falsenegatives, row.names=NULL)
      colnames(confusions) <- c('class', 'tp', 'tn', 'fp', 'fn')
      list(p0=p0, confusions=confusions)
    }
  }
 
  ## parameter handling
  if (tunecontrol$sampling == "cross")
    validation.x <- validation.y <- NULL
  useFormula <- is.null(train.y)
  if (useFormula && (is.null(data) || length(data) == 0))
    data <- model.frame(train.x)
  if (is.vector(train.x)) train.x <- t(t(train.x))
  if (is.data.frame(train.y))
    train.y <- as.matrix(train.y)
 
  ## prepare training indices
  if (!is.null(validation.x)) tunecontrol$fix <- 1
  n <- (nrow(if (useFormula) data
            else train.x))
  perm.ind <- sample(n)
  if (tunecontrol$sampling == "cross") {
    if (tunecontrol$cross > n)
      stop(sQuote("cross"), " must not exceed sampling size!")
    if (tunecontrol$cross == 1)
      stop(sQuote("cross"), " must be greater than 1!")
  }
  train.ind <- (if (tunecontrol$sampling == "cross")
                  tapply(1:n, cut(1:n, breaks = tunecontrol$cross), function(x) perm.ind[-x])
                else if (tunecontrol$sampling == "fix")
                  list(perm.ind[1:trunc(n * tunecontrol$fix)])
                else
                ## bootstrap
                  lapply(1:tunecontrol$nboot,
                         function(x) sample(n, n * tunecontrol$boot.size, replace = TRUE))
  )
 
  ## find best model
  parameters <- (if(is.null(ranges))
                  data.frame(dummyparameter = 0)
                else
                  expand.grid(ranges))
  p <- nrow(parameters)
  if (!is.logical(tunecontrol$random)) {
    if (tunecontrol$random < 1)
      stop("random must be a strictly positive integer")
    if (tunecontrol$random > p) tunecontrol$random <- p
    parameters <- parameters[sample(1:p, tunecontrol$random),]
  }
 
  ## - loop over all models
  # concatenate arbitrary mc-arguments with explicit X and FUN arguments
  train_results<-do.call(what="mclapply", args=c(mc.control, list(X=1:p, ..., FUN=function(para.set) {
    sampling.errors <- c()
    sampling.confusions <- c()
 
    ## - loop over all training samples
    for (sample in 1:length(train.ind)) {
      repeat.errors <- c()
      repeat.confusions <- c()
 
      ## - repeat training `nrepeat' times
      for (reps in 1:tunecontrol$nrepeat) {
 
        ## train one model
        pars <- if (is.null(ranges))
          NULL
        else
          lapply(parameters[para.set,,drop = FALSE], unlist)
 
        model <- if (useFormula)
          do.call(method, c(list(train.x,
                                 data = data,
                                 subset = train.ind[[sample]]),
                            pars, list(...)
          )
          )
        else
          do.call(method, c(list(train.x[train.ind[[sample]],],
                                 y = train.y[train.ind[[sample]]]),
                            pars, list(...)
          )
          )
 
        ## predict validation set
        pred <- predict.func(model,
                             if (!is.null(validation.x))
                               validation.x
                             else if (useFormula)
                               data[-train.ind[[sample]],,drop = FALSE]
                             else if (inherits(train.x, "matrix.csr"))
                               train.x[-train.ind[[sample]],]
                             else
                               train.x[-train.ind[[sample]],,drop = FALSE]
        )
 
        ## compute performance measure
        true.y <- if (!is.null(validation.y))
          validation.y
        else if (useFormula) {
          if (!is.null(validation.x))
            resp(train.x, validation.x)
          else
            resp(train.x, data[-train.ind[[sample]],])
        } else
          train.y[-train.ind[[sample]]]
 
        if (is.null(true.y)) true.y <- rep(TRUE, length(pred))
 
        if (!is.null(tunecontrol$error.fun))
          repeat.errors[reps] <- tunecontrol$error.fun(true.y, pred)
        else if ((is.logical(true.y) || is.factor(true.y)) && (is.logical(pred) || is.factor(pred) || is.character(pred))) { ## classification error
          l <- classAgreement2(table(pred, true.y))
          repeat.errors[reps] <- (1 - l$p0) # wrong classification rate
          if(confusionmatrizes) {
            repeat.confusions[[reps]] <- l$confusions
          }
        } else if (is.numeric(true.y) && is.numeric(pred)) ## mean squared error
          repeat.errors[reps] <- crossprod(pred - true.y) / length(pred)
        else
          stop("Dependent variable has wrong type!")
      }
      sampling.errors[sample] <- tunecontrol$repeat.aggregate(repeat.errors)
      # TODO potentially implement separate aggregation of tp tn fp fn values. currently those are taken with correlate to the least error.
      if(confusionmatrizes) {
        sampling.confusions[[sample]] <- repeat.confusions[repeat.errors == sampling.errors[sample]][[1]]
      }
    }
    # TODO potentially implement separate aggregation of tp tn fp fn values. currently uses the same as for error / variance aggregation
    if(!confusionmatrizes) {
      list(
        model.error=tunecontrol$sampling.aggregate(sampling.errors),
        model.variance=tunecontrol$sampling.dispersion(sampling.errors))
    } else {
      # create one confusion data frame
      confusions <- ldply(sampling.confusions, data.frame)
      # calculate aggregate / disperse values per class
      confusions <- ldply(lapply(X=split(confusions, confusions$class), FUN=function(classdf) {
        class=unique(classdf$class)
        # only take numeric values
        classdf[,c('tp','tn','fp','fn')]
        # calculate aggregate / disperse values for this class
        aggregated <- apply(X=classdf[,c('tp','tn','fp','fn')], MAR=2, FUN=tunecontrol$sampling.aggregate)
        dispersions <- apply(X=classdf[,c('tp','tn','fp','fn')], MAR=2, FUN=tunecontrol$sampling.dispersion)
        # make 1 row dataframe out of it (combine rows later with outer ldply)
        t(data.frame(c(value=aggregated,dispersion=dispersions)))
      }), data.frame)
      colnames(confusions) <- c('class', 'tp.value', 'tn.value', 'fp.value', 'fn.value', 'tp.dispersion', 'tn.dispersion', 'fp.dispersion', 'fn.dispersion')
      # calculate mean confusion matrix values (mean of all classes) for best model
      confusions.mean <- data.frame(t(apply(X=confusions[,c('tp.value','tn.value','fp.value','fn.value','tp.dispersion','tn.dispersion','fp.dispersion','fn.dispersion')], MAR=2, FUN=mean)))
      colnames(confusions.mean) <- c('tp.value', 'tn.value', 'fp.value', 'fn.value', 'tp.dispersion', 'tn.dispersion', 'fp.dispersion', 'fn.dispersion')
      list(
        model.error=tunecontrol$sampling.aggregate(sampling.errors),
        model.variance=tunecontrol$sampling.dispersion(sampling.errors),
        model.confusions=confusions,
        model.confusions.mean=confusions.mean
      )
    }
  })))
#   print('mctune: mclapply done.')
#   print(train_results)
  model.errors <- unlist(lapply(train_results,function(x)x$model.error))
  model.variances <- unlist(lapply(train_results,function(x)x$model.variance))
  if(confusionmatrizes){
    model.confusions <- lapply(train_results,function(x)x$model.confusions)
    model.confusions.mean <- ldply(lapply(train_results,function(x)x$model.confusions.mean), data.frame)
  }
 
  ## return results
  best <- which.min(model.errors)
  pars <- if (is.null(ranges))
    NULL
  else
    lapply(parameters[best,,drop = FALSE], unlist)
  structure(list(best.parameters  = parameters[best,,drop = FALSE],
                 best.performance = model.errors[best],
                 method           = if (!is.character(method))
                   deparse(substitute(method)) else method,
                 nparcomb         = nrow(parameters),
                 train.ind        = train.ind,
                 sampling         = switch(tunecontrol$sampling,
                                           fix = "fixed training/validation set",
                                           bootstrap = "bootstrapping",
                                           cross = if (tunecontrol$cross == n) "leave-one-out" else
                                             paste(tunecontrol$cross,"-fold cross validation", sep="")
                 ),
                 performances     = if (tunecontrol$performances) cbind(parameters, error = model.errors, dispersion = model.variances),
                 confusionmatrizes = if (confusionmatrizes) model.confusions,
                 confusionmatrizes.mean = if(confusionmatrizes) model.confusions.mean,
                 best.confusionmatrizes = if(confusionmatrizes) model.confusions[[best]],
                 best.confusionmatrizes.mean = if(confusionmatrizes) model.confusions.mean[best,],
                 best.model       = if (tunecontrol$best.model) {
                   modeltmp <- if (useFormula)
                     do.call(method, c(list(train.x, data = data),
                                       pars, list(...)))
                   else
                     do.call(method, c(list(x = train.x,
                                            y = train.y),
                                       pars, list(...)))
                   call[[1]] <- as.symbol("best.tune")
                   modeltmp$call <- call
                   modeltmp
                 }
  ),
            class = "tune"
  )
}
