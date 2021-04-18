#' @title  Multiple Thresholds detection in mean function
#'
#' @description  Multiple Thresholds detection in mean function for dense or sparse functional data.
#'
#'
#' @param yin A vector containing the observed values.
#' @param xin A vector containing the observation time points corresponding to yin.
#' @param xout A vector containing the observation time points for mean estimation.
#' @param weight A vector of weight for each observations.
#' @param zeta Cut-off threshold for change point detection.
#' @param npoly The degree of polynomial. Default is 1 for local linear smoothing.
#' @param nder The order of derivative. Default is 0 for local linear smoothing, and should smaller than npoly.
#' @param NbGrid The number of support points in the jump points detection; numeric - default: 101.
#' @param kernel Smoothing kernel choice, common for mu and covariance; "rect", "epan", "quar" - default: "epan".
#' @param h_tau The selected (or user specified) bandwidth for jump locations estimate.
#' @param h_d The selected (or user specified) bandwidth for jump sizes estiamte.
#' @param refined A refined stage for jump size estimation after the jump locations are detected; logical - default: TRUE
#'
#'
#'
#' @return A list containing the following fields: If the change points are detected, the function
#' output the followings
#' \item{mu_jumptime}{The detected change point locations. If no change point detected, \code{mu_jumptime=0}}
#' \item{mu_jumpsize}{The estimated jump sizes. If no change point detected, \code{mu_jumpsize=0}}
#' \item{mu_jumpsize_h_tau}{The detected change point locations. If no change point detected, \code{mu_jumptime=0}}
#' \item{mu}{For no change point case. A vector of length obsGrid containing the mean function estimate.}
#' \item{muout}{For change point case. A vector of length xout containing the mean function estimate.}
#' \item{xout}{A vector containing the observation time points for mean estimation.}
#' \item{timings}{A vector with execution times for the basic parts of the FMbreaks call.}
#' \item{h_tau}{The selected (or user specified) bandwidth for jump locations estimate.}
#' \item{h_d}{The selected (or user specified) bandwidth for jump sizes estiamte.}
#'
#'
#'
#' @references
#' \cite{Li, J., Li, Y., and Hsing, T. (2021). " On Functional Processes with Multiple Discontinuities".}
#'
#' @export
#'





######### Get Mean function, Covariance function and sigma2
### mean estimation with multiple breaks
MeanBreaksFP <- function(yin, xin, weight, xout, NbGrid, h_tau, h_d, zeta,
                         npoly = 1, nder = 0, kernel, refined = FALSE) {

  firsttsMeanBreaksFP <- Sys.time() #First time-stamp for MeanBreaksFP


  ####### change-points detection procedure ############
  firsttsJump <- Sys.time()

  # Generate basic grids for jump detect:
  obsGrid = sort(unique(xin));
  if(is.null(NbGrid)){
    jumpGrid = obsGrid;
  } else {
    jumpGrid = seq(min(obsGrid), max(obsGrid),length.out = NbGrid);
  }
  # jumpGrid = seq(0, 1, length.out = NbGrid);
  ## Local linear estimate based on one-sided kernel
  D_h <- jumpGrid[h_tau <= jumpGrid & jumpGrid <= 1- h_tau];
  ##
  mu_est_lr = CPPlwls1d_LR(bw = h_tau, kernel_type = kernel, win = weight,
                           xin = xin, yin = yin, xout = D_h, npoly = npoly)
  mu_diff = mu_est_lr[,2]-mu_est_lr[,1]
  ##
  # mu_est_lr = LL_lr(data = cbind(xin, yin), time = D_h, weight = weight, h = h_tau)
  # mu_diff = mu_est_lr[,1]-mu_est_lr[,2]
  ##
  mu_diff_abs = abs(mu_diff)
  mu_diff_time = D_h[order(-mu_diff_abs)]
  mu_diff_size = mu_diff[order(-mu_diff_abs)]


  ##  only focus on jumps size that are bigger than zeta
  timeindex = which(abs(mu_diff_size) > zeta)
  ll=length(timeindex)

  if(ll == 0) {
    cat('no change point detected.', '\n')

    lasttsJump <- Sys.time()

    ## Mean estimation
    firsttsMu <- Sys.time()

    # mu: the smoothed mean curve evaluated at times 'obsGrid'
    mu = Lwls1D(h_tau, kernel_type = kernel, npoly = npoly, nder = nder,
                xin = xin, yin= yin, xout = obsGrid, win = weight)
    # convert mu to truncated xout
    # muout <- ConvertSupport(obsGrid, toGrid = xout, mu=mu)
    muout = Lwls1D(h_tau, kernel_type = kernel, npoly = npoly,
                   nder = nder, xin = xin, yin= yin, xout = xout,
                   win = weight)
    lasttsMu <- Sys.time()

    timestamps = c(lasttsJump, lasttsMu, firsttsMeanBreaksFP,  firsttsJump, firsttsMu)


    if(is.null(timestamps)) {
      timings = NULL;
    } else {
      timestamps = c(Sys.time(), timestamps)
      timings = round(digits=3, timestamps[1:3]-timestamps[4:6]);
      names(timings) <- c('total', 'jump', 'mu')
    }


    return(list(mu_jumptime = 0,
                mu_jumpsize = 0,
                mu_jumptime_refine = 0,
                mu_jumpsize_refine = 0,
                obsGrid = obsGrid, xout = xout,
                mu = mu, muout = muout,  weight = weight,
                h_tau = h_tau, timings = timings))


  } else {

    mu_diff_time=mu_diff_time[timeindex]
    mu_diff_size=mu_diff_size[timeindex]

    # find cluster centers
    mu_jumptime=mu_diff_time[1]
    mu_jumpsize=mu_diff_size[1]
    for (i in 2:ll){
      if (sum(abs(mu_diff_time-mu_jumptime[i-1]) > 2*h_tau) >0){
        index=which(abs(mu_diff_time-mu_jumptime[i-1]) > 2*h_tau)
        mu_diff_time = mu_diff_time[index]
        mu_diff_size = mu_diff_size[index]
        mu_jumptime=append(mu_jumptime, mu_diff_time[1])
        mu_jumpsize=append(mu_jumpsize, mu_diff_size[1])
      }
      else{
        break
      }
    }
    mu_jumpsize_h_tau = mu_jumpsize[order(mu_jumptime)]
    mu_jumptime = sort(mu_jumptime)

    # cat('change point at:', mu_jumptime, '\n')
    lasttsJump <- Sys.time()


    ## refine stage to estimate jump size
    rho_d = 1.1*h_tau^2
    names(rho_d) = c("rho")
    jumpset_l = c(mu_jumptime - rho_d)
    jumpset_r = c(mu_jumptime + rho_d)

    ## jump size based on the bandwidth h_d
    mu_est_lr_l = CPPlwls1d_LR(bw = h_d, kernel_type = kernel, win = weight,
                               xin = xin, yin = yin, xout = jumpset_l,
                               npoly = npoly)
    mu_est_lr_r = CPPlwls1d_LR(bw = h_d, kernel_type = kernel, win = weight,
                               xin = xin, yin = yin, xout = jumpset_r,
                               npoly = npoly)
    mu_jumpsize_h_d = mu_est_lr_r[,2] - mu_est_lr_l[,1]

    ## refine the jump location
    if(refined){
      mu_jumptime = mu_jumptime[abs(mu_jumpsize_h_d) > zeta]
      mu_jumpsize_h_d = mu_jumpsize_h_d[abs(mu_jumpsize_h_d) > zeta]
    }

    if(length(mu_jumptime)==0){
      mu_jumptime = 0
      mu_jumpsize_h_d = 0

      cat('no change point after refine.', '\n')
    } else {
      cat('change point at:', mu_jumptime, '\n')
    }


    lasttsJump <- Sys.time()

    ## Mean estimation
    firsttsMu <- Sys.time()

    yyin <- yin - sapply(xin, function(z) sum(mu_jumpsize_h_d*(z >= mu_jumptime)))
    # nu: the smoothed mean curve evaluated at times 'xout' with same bandwidth of h_tau
    nu = Lwls1D(h_tau, kernel_type = kernel, npoly = npoly, nder = nder,
                xin = xin, yin= yyin, xout = obsGrid, win = weight)
    mu <- nu + sapply(obsGrid, function(z)sum(mu_jumpsize_h_d *(z >= mu_jumptime)))
    nuout = Lwls1D(bw = h_tau, kernel_type = kernel, npoly = npoly, nder = nder,
                   xin = xin, yin= yyin, xout = xout, win = weight)
    muout <- nuout + sapply(xout, function(z)sum(mu_jumpsize_h_d *(z >= mu_jumptime)))

    lasttsMu <- Sys.time()

    timestamps = c(lasttsJump, lasttsMu, firsttsMeanBreaksFP,  firsttsJump, firsttsMu)


    if(is.null(timestamps)) {
      timings = NULL;
    } else {
      timestamps = c(Sys.time(), timestamps)
      timings = round(digits=3, timestamps[1:3]-timestamps[4:6]);
      names(timings) <- c('total', 'jump', 'mu')
    }



    return(list(mu_jumptime = mu_jumptime,
                mu_jumpsize_h_tau = mu_jumpsize_h_tau,
                mu_jumpsize = mu_jumpsize_h_d,
                mu = mu, muout = muout,
                obsGrid = obsGrid, xout = xout, weight = weight,
                h_tau = h_tau, h_d = h_d, rho_d = rho_d, timings = timings))
  }

}


### covariance estimation: this function is faster than fdapace package
GetRawCov <- function(y,t,obsGridnew, mu, dataType, error){
  #  obtain raw covariance
  #  Input y :       1*n cell array of the observed repeated measurements from n subjects
  #  Input t :       1*n cell array of the observed time points from n subjects
  #  Input obsGridnew:  1*m vector of time points correspond to mu
  #  Input mu:       1*m vector of fitted mean functions from Step I, corresponding to
  #                 pooled unique time points from t
  #  Input dataType: Output of IsRegular()
  #  Input error:    TRUE with measurement error assumption
  #                  FALSE without measurement error assumption
  #
  #  Output res: a list that contains tPairs, cxxn, indx,win and cyy
  #     tPairs:  N * 2  vector denotes the  pairs of time points for subject
  #                 concatenating as two vectors
  #                if error = 1, all (t_ij, t_ij) will be removed
  #       cxxn:    1 * N vector of raw covariance corresponding to tPairs
  #       indx:    1 * N vector of indices for each subject
  #        win:    1 * N weight matrix for the 2-D smoother for covariance function
  #        cyy:    1 * M vector of raw covariance corresponding to all pairs of time points,
  #                i.e., it is the same as cxxn if error = 0
  #       diag:    if error == TRUE: 2-column matrix recording raw covariance along the diagonal direction (col 2)
  #                and the corresponding observed time points (col 1)
  #                if error == FALSE: NULL

  ncohort <- length(y);
  obsGrid <- sort(unique(unlist(t)))
  mu <- MapX1D(x = obsGridnew, y = mu, newx = obsGrid);
  count <- NULL
  indx = NULL
  diag = NULL

  if(dataType %in% c('Sparse', 'DenseWithMV')){

    Ys = lapply(X = y, FUN=pracma::meshgrid) #pracma
    Xs = lapply(X = t, FUN=pracma::meshgrid) #pracma

    # vectorise the grids for y & t
    xx1 = unlist(do.call(rbind, lapply(Xs, '[', 'X')) )
    xx2 = unlist(do.call(rbind, lapply(Xs, '[', 'Y')) )
    yy2 = unlist(do.call(rbind, lapply(Ys, '[', 'Y')) )
    yy1 = unlist(do.call(rbind, lapply(Ys, '[', 'X')) )

    # get id1/2 such that xx1/2 = q(id1/2), where q = unique(xx1/2)
    # id1 = apply(X= sapply(X=xx1, FUN='==',  ...=sort(unique(xx1)) ),MARGIN=2, FUN=which)
    # id2 = apply(X= sapply(X=xx2, FUN='==',  ...=sort(unique(xx2)) ),MARGIN=2, FUN=which)
    # This is more stable and faster than the fdapace package
    id1 = match(xx1, sort(unique(xx1)))
    id2 = match(xx2, sort(unique(xx2)))
    cyy = ( yy1 - mu[ id1]) * (yy2 - mu[id2] )

    # index for subject i
    # indx = unlist(sapply( 1:length(y), function(x) rep(x,  (unlist(lapply(length, X= y))[x])^2) ))
    # This is more stable and faster.
    indx = rep( 1:length(y), times =  unlist(lapply(y, length))^2)

    tPairs = matrix( c(xx1, xx2), nrow=length(xx1), ncol=2);

    if(error){
      tneq = which(xx1 != xx2)
      teq = which(xx1 == xx2)
      indx = indx[tneq];
      diag = matrix(c(tPairs[teq,1], cyy[teq]), ncol = 2)
      tPairs = tPairs[tneq,];
      cxxn = cyy[tneq];
    }else{
      cxxn = cyy;
    }

    # win = pracma::ones(1, length(cxxn));
    # count = GetCount(tPairs)...

  }else if(dataType == 'Dense'){

    yy = t(matrix(unlist(y), length(y[[1]]), ncohort))
    MU = t(matrix( rep(mu, times=length(y)), ncol=length(y)))
    t1 = t[[1]]
    yy = yy - MU;
    cyy = t(yy) %*% yy / ncohort
    cyy = as.vector(t(cyy))
    cxxn = cyy;
    xxyy = pracma::meshgrid(t1); # pracma

    tPairs =  (matrix( c(c(xxyy$X), c(xxyy$Y)), ncol = 2))

    if(error){
      tneq = which(tPairs[,1] != tPairs[,2])
      teq = which(tPairs[,1] == tPairs[,2])
      diag = matrix(c(tPairs[teq,1], cyy[teq]), ncol = 2)
      tPairs = tPairs[tneq,];
      cxxn = cyy[tneq];
    }else{
      cxxn = cyy;
    }

    # win = pracma::ones(1, length(cxxn));
  }else if(dataType == 'RegularWithMV'){
    stop("This is not implemented yet. Contact Pantelis!")
  }else {
    stop("Invalid 'dataType' argument type")
  }

  result <- list( 'tPairs' = tPairs, 'cxxn' = cxxn, 'indx' = indx, # 'win' = win,
                  'cyy' = cyy, 'diag' = diag, 'count' = count, 'error' = error, 'dataType' = dataType);

  class(result) <- "RawCov"
  return(result)
}


# The output outGrid of this function is a (potentially) truncated grid.
GetSmoothedCovarSurface <- function(y, t, mu, obsGrid, regGrid, optns, useBinnedCov=FALSE) {

  dataType <- optns$dataType
  error <- optns$error
  kern <- optns$kernel
  userBwCov <- optns$userBwCov
  methodBwCov <- optns$methodBwCov
  verbose <- optns$verbose
  rotationCut <- optns$rotationCut

  # get the truncation of the output grids.
  outPercent <- optns$outPercent
  buff <- .Machine$double.eps * max(abs(obsGrid)) * 10
  rangeGrid <- range(regGrid)
  minGrid <- rangeGrid[1]
  maxGrid <- rangeGrid[2]
  cutRegGrid <- regGrid[regGrid > minGrid + diff(rangeGrid) * outPercent[1] -
                          buff &
                          regGrid < minGrid + diff(rangeGrid) * outPercent[2] +
                          buff]

  # Get raw covariance, unless user covariance/sigma2 are specified.
  if (is.null(optns[['userCov']]) ||
      (is.null(optns[['userSigma2']]) && error)) {

    rcov <- GetRawCov(y, t, obsGrid, mu, dataType, error)
    if (useBinnedCov && methodBwCov == 'CV') {
      stop('If methodBwCov == \'CV\' then we must use the unbinned rcov.')
    }

    if (useBinnedCov) {
      rcov <- BinRawCov(rcov)
    }
  } else {
    rcov <- NULL
  }

  # Obtain smoothed covariance.
  if( !is.null(optns$userCov)) { # If covariance function is provided
    optns$userCov$t <- as.numeric(optns$userCov$t)
    optns$userCov$cov <- as.numeric(optns$userCov$cov)
    rangeUser <- range(optns$userCov$t)
    rangeCut <- range(cutRegGrid)
    if( rangeUser[1] > rangeCut[1] + buff ||
        rangeUser[2] < rangeCut[2] - buff   ) {
      stop('The range defined by the user provided covariance does not cover the support of the data.')
    }

    bwCov  = NULL
    smoothCov = ConvertSupport(fromGrid = optns$userCov$t, cutRegGrid, Cov =  optns$userCov$cov)

  } else { # estimate the smoothed covariance

    if (userBwCov == 0) { # bandwidth selection
      if (methodBwCov %in% c('GCV', 'GMeanAndGCV')) { # GCV
        gcvObj <- GCVLwls2DV2(obsGrid, regGrid, kern=kern, rcov=rcov, verbose=verbose, t=t)
        bwCov <- gcvObj$h
        if (methodBwCov == 'GMeanAndGCV') {
          bwCov <- sqrt(bwCov * gcvObj$minBW)
        }
      } else if (methodBwCov == 'CV') { # CV 10 fold
        gcvObj <- GCVLwls2DV2(obsGrid, regGrid, kern=kern, rcov=rcov, t=t,
                              verbose=optns$verbose,
                              CV=optns[['kFoldMuCov']], useBW1SE = optns$useBW1SE)
        bwCov <- gcvObj$h
      }
    } else if (userBwCov != 0) {
      bwCov <- userBwCov
    }

    if (!useBinnedCov) {
      smoothCov <- Lwls2D(bwCov, kern, xin=rcov$tPairs, yin=rcov$cxxn,
                          xout1=cutRegGrid, xout2=cutRegGrid)
    } else {
      smoothCov <- Lwls2D(bwCov, kern, xin=rcov$tPairs, yin=rcov$meanVals,
                          win=rcov$count, xout1=cutRegGrid, xout2=cutRegGrid)
    }
  }

  # Obtain the error sigma2.
  if (error) {
    if (!is.null(optns[['userSigma2']])) {
      sigma2 <- optns[['userSigma2']]
    } else if (!is.null(optns[['userCov']])) {
      a0 = min(regGrid)
      b0 = max(regGrid)
      lint = b0 - a0
      middleCutRegGrid <- cutRegGrid > a0 + lint * rotationCut[1] - buff &
        cutRegGrid < a0 + lint * rotationCut[2] + buff
      if (useBinnedCov) {
        diagT <- rcov[['tDiag']]
        diagVal <- rcov[['diagMeans']]
      } else {
        diagTV <- aggregate(rcov[['diag']][, 2], list(rcov[['diag']][, 1]), mean)
        diagT <- diagTV[, 1]
        diagVal <- diagTV[, 2]
      }
      diagEst <- approx(diagT, diagVal, cutRegGrid[middleCutRegGrid])[['y']]
      sigma2 <- mean(diagEst - diag(smoothCov)[middleCutRegGrid])

    } else { # has to estimate sigma2 from scratch
      sigma2 <- PC_CovE(obsGrid, regGrid, bwCov, rotationCut=rotationCut, kernel=kern, rcov=rcov)$sigma2
    }

    if(sigma2 < 0) {
      if(verbose){
        warning("Estimated sigma2 is negative and thus is reset to 1e-6.")
      }
      sigma2 <- 1e-6
    }

  } else { # error=FALSE
    sigma2 <- NULL
  }

  res <- list(rawCov = rcov,
              smoothCov = (smoothCov + t(smoothCov)) / 2,
              bwCov = bwCov,
              sigma2 = sigma2,
              outGrid = cutRegGrid)
  class(res) <- "SmoothCov"
  return(res)
}

