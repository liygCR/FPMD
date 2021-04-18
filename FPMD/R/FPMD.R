#' @title  Functional Processes with Multiple Discontinuities
#'
#' @description  Multiple Thresholds detection in mean function for dense or sparse functional data.
#' The bandwidths are selected by a k-folds cross validation.
#'
#'
#' @param Ly A list of \emph{n} vectors containing the observed values for each individual. Missing values specified by \code{NA}s are supported for dense case (\code{dataType='dense'}).
#' @param Lt A list of \emph{n} vectors containing the observation time points for each individual corresponding to y. Each vector should be sorted in ascending order.
#' @param wi A vector of weight for each observation.
#' @param Wtype Weight structure specified for estimation. It works if \code{wi = NULL}.
#' @param zeta Cut-off threshold for change point detection. The default is \code{NULL}.
#' @param cutoff The method for construct the cut-off statistics. Default is \code{max}. For the case
#'                       observations are much larger than subjections, we suggest to use the \code{min}.
#' @param alpha The confidence level for the cut-off threshold. The default is
#'              \code{alpha = 0.05}, and for heavy tailed process, we suggest a higher level \code{0.1}.
#' @param npoly The degree of polynomial. Default is 1 for local linear smoothing.
#' @param nder The order of derivative. Default is 0 for local linear smoothing, and should smaller than npoly.
#' @param bw.seq The selected (or user specified) bandwidth for cross validation. The default is \code{NULL}.
#' @param nRegGrid The number of support points in each direction of covariance surface; numeric - default: 101.
#' @param NbGrid The number of support points in the jump points detection; numeric - default: 101.
#' @param kernel Smoothing kernel choice, common for mu and covariance; "rect", "epan", "quar" - default: "epan".
#' @param kFolds Number of folds for cross validation - default: 5.
#' @param refined A refined stage for jump size estimation after the jump locations are detected; logical - default: TRUE
#' @param individualCurveJump Jump detection for each curve by the BIC proposed in Xia and Qiu, 2016.
#' @param M_max The maximum number of jumps for individual curve when using the BIC.
#'
#'
#'
#' @return A list containing the following fields: If the change points are detected, the function
#' output the followings
#' \item{mu_jumptime}{The detected change point locations. If no change point detected, \code{mu_jumptime=0}}
#' \item{mu_jumpsize}{The estimated jump sizes. If no change point detected, \code{mu_jumpsize=0}}
#' \item{mu_jumpsize_h_tau}{The detected change point locations. If no change point detected, \code{mu_jumptime=0}}
#' \item{zeta}{The cut-off threshold for change point detection.}
#' \item{mu}{For no change point case. A vector of length obsGrid containing the mean function estimate.}
#' \item{muWork}{For change point case. A vector of length obsGrid containing the mean function estimate.}
#' \item{wi}{The weight for each observation.}
#' \item{mi}{Number of observations for each subject.}
#' \item{timings}{A vector with execution times for the basic parts of the FMbreaks call.}
#' \item{timings_jump}{A vector with execution times for the basic parts of jump detection.}
#' \item{workGrid}{A vector of length nWorkGrid.}
#' \item{regGrid}{A vector of length regGrid.}
#' \item{sigma2}{Variance for measure error.}
#' \item{smoothCov}{A nWorkGrid by nWorkGrid matrix of the smoothed covariance surface.}
#' \item{fittedCov}{A nWorkGrid by nWorkGrid matrix of the fitted covariance surface, which is guaranteed to be non-negative definite.}
#' \item{optns}{A list of actually used options.}
#' \item{h_tau}{The selected (or user specified) bandwidth for jump locations estimate.}
#' \item{h_d}{The selected (or user specified) bandwidth for jump sizes estiamte.}
#' \item{LX}{The fitted individual curves.}
#' \item{indJump}{The jump detection for each individual curve.}
#'
#'
#' @import fdapace
#' @import zoo
#' @import RcppEigen
#' @import Rcpp
#' @import locpol
#' @importFrom extraDistr rtpois
#' @importFrom pracma meshgrid isempty
#' @importFrom stats integrate quantile rnorm rpois rt runif sd aggregate approx density
#' @importFrom graphics legend matplot polygon
#' @importFrom grDevices rgb
#' @useDynLib FPMD, .registration = TRUE
#'
#' @examples
#'
#' \dontrun{
#' # setting 1
#' ## 3 change points
#' tau = c(.25,.5,.75) # jump locations
#' disc = c(.5,-.4,.4) # jump sizes
#' smoothmean=function(x) x^2+sin(2*pi*x) + cos(2*pi*x)
#' jump = function(x) sum((x>=tau)*disc) # jump function
#' mu_fun= function(t){smoothmean(t)+ sapply(t, jump)}
#' n = 400
#' Q = 20
#'
#' # set.seed(123)
#' ## generate gaussian process
#' ## "cos", "sin", "fourier", "legendre01", "poly"
#' data = MakeFPMBData(n = n, Q = Q, proc= "gaussian", muFun = mu_fun, K = 3,
#'                     lambda = 1/(seq_len(3)+1)^2, sigma = 0.2, basisType = "fourier")
#' Lt = data$Lt
#' Ly = data$Ly
#' bw.seq = seq(0.01, 0.13, by = 0.02)
#' system.time(resCP <- FPMD(Ly = Ly, Lt = Lt, wi = NULL, Wtype = "MIX", zeta = NULL,
#'                               bw.seq = 0.06, NbGrid = 101, kFolds = 5, refined = TRUE,
#'                               individualCurveJump = FALSE, nRegGrid = 101,
#'                               npoly = 1, nder = 0, alpha = 0.05, cutoff = max))
#' h_tau = resCP$h_tau
#' h_d = resCP$h_d
#' zeta = resCP$zeta
#' mu_jumptime = resCP$mu_jumptime
#' mu_jumpsize = resCP$mu_jumpsize
#'
#' ## pointwise confidence band
#' PCBplot(res = resCP)
#' lines(resCP$obsGrid, mu_fun(resCP$obsGrid), type="l", col = 'blue' )
#'
#' }
#'
#'
#' @references
#' \cite{Li, J., Li, Y., and Hsing, T. (2021). " On Functional Processes with Multiple Discontinuities".}
#'
#' @export
#'



## main function
FPMD <- function(Ly, Lt, wi = NULL, Wtype = c("OBS", "SUBJ", "MIX", "OPT"),
                 zeta = NULL, cutoff = max, alpha = 0.05, bw.seq = NULL,
                 nRegGrid = 101, NbGrid = 101, kernel='epan',
                 npoly = 1, nder = 0, kFolds = 5, refined = FALSE,
                 individualCurveJump = FALSE, M_max){

  firsttsbreaks <- Sys.time()
  # Check the data validity for further analysis
  CheckData(Ly,Lt)

  inputData <-  HandleNumericsAndNAN(Ly,Lt);
  # inputData <- HandleNumericsAndNAN(Ly,Lt);
  Ly <-  inputData$Ly;
  Lt <-  inputData$Lt;

  # Set the options structure members that are still NULL
  optns <- list(methodBwMu = 'CV',
                kFoldMuCov = kFolds,
                nRegGrid = nRegGrid,
                error = TRUE,
                kernel= kernel,
                verbose=FALSE)
  optns <- SetOptions(Ly, Lt, optns);
  optns$useBW1SE = FALSE;
  # optns$dataType = 'Sparse'

  ## weight type
  Wtype = match.arg(Wtype)
  ### Check the options validity for the PCA function.
  numOfCurves = length(Ly);
  CheckOptions(Lt, optns, numOfCurves)


  ###
  firsttsCVmu <- Sys.time()
  ## mean function
  if(length(bw.seq) == 1) {
    h_tau = h_d = bw.seq
    zeta = zetaFun(bw = bw.seq, alpha = alpha, Lt = Lt, Ly = Ly, optns = optns,
                   cutoff = max, Wtype = Wtype)
  } else {
    ## select by CV procedure
    tunings = CVbandwidth(bw.seq = bw.seq, zeta = zeta, Ly = Ly, Lt = Lt,
                          npoly = npoly, nder = nder, alpha = alpha,
                          optns = optns, NbGrid = NbGrid, Wtype = Wtype,
                          refined = refined, kFolds = kFolds, cutoff = cutoff)
    h_tau = tunings$bopt[1]
    h_d = tunings$bopt[2]
    zeta = tunings$zeta
  }

  ##
  cat(' ', '\n')
  cat('Final selected: ', 'h_tau = ', h_tau, ';', 'h_d = ', h_d, '\n')
  cat('and selected: ', 'zeta = ', zeta, '\n')
  # cat('change point at:', smcObj$mu_jumptime, '\n')
  ###
  obsGrid = sort(unique( c(unlist(Lt))));
  ## cut in the region [h, 1-h]
  regGrid = seq( max(min(obsGrid), h_tau), min(max(obsGrid), 1- h_tau), length.out = optns$nRegGrid);
  # regGrid = seq(min(obsGrid), max(obsGrid),length.out = optns$nRegGrid);
  mi = sapply(Lt, length)
  xin = unlist(Lt);
  yin = unlist(Ly)[order(xin)];
  xin = sort(xin);
  # different type of weight
  win = weightFun(h = h_tau, wi = wi, mi = mi, n = numOfCurves, Wtype = Wtype)
  # Get the mean function using the bandwith estimated above:
  smcObj = MeanBreaksFP(xin = xin, yin= yin, weight = win$weight, xout = regGrid,
                        h_tau = h_tau, h_d = h_d, zeta = zeta, refined = refined,
                        npoly = npoly, nder = nder, kernel = optns$kernel, NbGrid = NbGrid)
  mu = smcObj$mu
  muWork = smcObj$muout
  mu_jumpsize = smcObj$mu_jumpsize
  rho = smcObj$rho_d
  lasttsCVmu <- Sys.time()



  firsttsCov <- Sys.time() #First time-stamp for calculation of the covariance
  ## Covariance function and sigma2
  scsObj = GetSmoothedCovarSurface(y = Ly, t = Lt, mu = mu, obsGrid = obsGrid,
                                   regGrid = regGrid, optns = optns)
  sigma2 <- scsObj[['sigma2']]

  # Get the results for the eigen-analysis
  eigObj = GetEigenAnalysisResults(smoothCov = scsObj$smoothCov, regGrid, optns, muWork = muWork)
  fittedCov = eigObj$fittedCov

  lasttsCov <- Sys.time()


  #
  timestamps = c(lasttsCVmu, lasttsCov, firsttsbreaks, firsttsCVmu, firsttsCov)


  if(is.null(timestamps)) {
    timings = NULL;
  } else {
    timestamps = c(Sys.time(), timestamps)
    timings = round(digits=3, timestamps[1:3]-timestamps[4:6]);
    names(timings) <- c('total', 'mu', 'cov')
  }


  ### individual curve estimate and jump detection
  LX = indJump = NULL
  if(individualCurveJump){

    LX = indJump = vector("list", numOfCurves)
    allObj = MeanBreaksFP(xin = xin, yin= yin, weight = win$weight, xout = sort(unlist(Lt)),
                          h_tau = h_tau, h_d = h_d, zeta = zeta, refined = refined,
                          npoly = npoly, nder = nder, kernel = optns$kernel, NbGrid = NbGrid)
    muall = allObj$muout
    ### predict the individual curve X_i
    Lmu = split(muall, rep(1:numOfCurves, mi))
    rm(muall)
    Leps = mapply(function(y, x) y - x, y = Ly, x = Lmu, SIMPLIFY = FALSE)

    LU = mapply(function(t, y){
      bw = regCVBwSelC(x = t, y = y, deg = npoly, interval = c(0, 0.2))
      # bw = length(t)^(-1/5)
      # print(bw)
      Lwls1D(bw = bw, kernel_type = optns$kernel, npoly = npoly,
             nder = nder, xin = t, yin= y, xout = t,
             win = rep(1, length(t)))
    }, t = Lt, y = Leps, SIMPLIFY = FALSE)
    LX = mapply(function(x,y)x+y, Lmu, LU, SIMPLIFY = FALSE)

    names(indJump) = paste('Curve', 1:numOfCurves)
    for (i in 1:numOfCurves) {
      indJump[[i]] = indMeanbreak(y = Ly[[i]], t = Lt[[i]], M_max = M_max, NbGrid = NbGrid,
                                  kernel = optns$kernel, npoly = npoly, nder = nder)
      # print(i)
    }

  }


  ##
  ret <- list(sigma2 = sigma2,
              obsGrid = obsGrid,
              workGrid = regGrid,
              muWork = muWork,
              mu = mu,
              smoothedCov = scsObj$smoothCov,
              fittedCov = fittedCov,
              optns = optns,
              h_tau = smcObj$h_tau,
              h_d = smcObj$h_d,
              mu_jumptime = smcObj$mu_jumptime,
              mu_jumpsize = smcObj$mu_jumpsize,
              mu_jumpsize_h_tau = smcObj$mu_jumpsize_h_tau,
              timings = timings,
              timings_jump = smcObj$timings,
              wi = win$wi,
              mi = mi,
              zeta = zeta,
              rho = rho,
              LX = LX,
              indJump = indJump)

  return(ret)


}

## different weight
weightFun <- function(h, wi, mi, n, Wtype = c("OBS", "SUBJ", "MIX", "OPT")){

  Wtype = match.arg(Wtype)
  if (is.null(wi)) {

    if(Wtype == "OBS"){
      wi = rep(1/sum(mi), n)
    } else if(Wtype == "SUBJ"){
      wi = 1/(n*mi)
    } else if(Wtype == "OPT"){
      wi = (1/h + mi - 1)^(-1)/sum(mi/(1/h + mi - 1))
    } else if(Wtype == "MIX") {
      c1 <- ( 1/(h*mean(mi)) + mean(mi^2)/mean(mi)^2 )/n
      c2 <- (mean(1/mi)/h + 1)/n
      alp <- c2/(c1+c2)
      wi <- alp/sum(mi) + (1-alp)/(n*mi)

    }

  }

  weight = rep(wi, times = mi)



  return(list(weight = weight, wi = wi))
}

## get the cut-off value zeta
zetaFun <- function(bw, alpha = 0.05, Lt, Ly, optns, cutoff = max,
                     Wtype = c("OBS", "SUBJ", "MIX", "OPT")){


  Wtype = match.arg(Wtype)
  # optns <- list(methodBwMu = 'CV',
  #               nRegGrid = nRegGrid,
  #               error = TRUE,
  #               kernel='epan',
  #               verbose=FALSE)
  # optns <- SetOptions(Ly, Lt, optns);

  ## Mean with one-sided smooth
  mi = sapply(Lt, length)
  numOfCurves = length(Ly);
  ##
  xin = unlist(Lt);
  yin = unlist(Ly)[order(xin)];
  xin = sort(xin);
  win = weightFun(h = bw, wi = NULL, mi = mi, n = numOfCurves, Wtype = Wtype)   # rep(1, length(xin));
  ## two one-side estimate on a cut region
  obsGrid = xin[xin >= bw & xin <= (1-bw)]
  mu_lr = CPPlwls1d_LR(bw = bw, kernel_type = optns$kernel, win = win$weight,
                       xin = xin, yin = yin, xout = obsGrid, npoly = 1)
  ## choose the one with smaller residual mean sequare
  id = which.min(apply(mu_lr, 2, function(x) mean((yin[xin >= bw & xin <= (1-bw)] - x)^2)))
  obsGrid = sort(unique(xin));
  mu = CPPlwls1d_LR(bw = bw, kernel_type = optns$kernel, win = win$weight,
                    xin = xin, yin = yin, xout = obsGrid, npoly = 1)[,id]
  ## Covariance function and sigma2
  # smooth cov and/or sigma2
  # regGrid = seq(bw, 1 - bw, length.out = optns$nRegGrid);
  regGrid = seq( max(min(obsGrid), bw), min(max(obsGrid), 1- bw), length.out = optns$nRegGrid);
  scsObj <- GetSmoothedCovarSurface(Ly, Lt, mu = mu, obsGrid, regGrid, optns)
  sigma2 <- scsObj[['sigma2']]
  smoothedCov <- scsObj$smoothCov
  # workGrid: possibly truncated version of the regGrid
  workGrid <- scsObj$outGrid
  # convert mu to truncated workGrid
  muWork <- ConvertSupport(fromGrid = obsGrid, toGrid = workGrid, mu=mu)
  # muWork = approx(obsGrid, mu, xout= workGrid)$y
  # Get the results for the eigen-analysis
  eigObj <- GetEigenAnalysisResults(smoothCov = smoothedCov, workGrid, optns, muWork = muWork)
  fittedCov <- eigObj$fittedCov


  ### epanechnikov kernel
  K_Epa = function(u) 0.75 * (1 - u^2) * (abs(u) <= 1);
  K_Epa_r = function(u, r) 0.75 * (1 - u^2) *u^r * (abs(u) <= 1);
  nu0 = integrate(K_Epa_r, 0, 1, r=0)$value
  nu1 = integrate(K_Epa_r, 0, 1, r=1)$value
  nu2 = integrate(K_Epa_r, 0, 1, r=2)$value
  phi_fun <- function(u) {
    res = K_Epa(u) * (nu2 - nu1*u)/(nu0*nu2 - nu1^2)
    return(res^2)
  }
  ## one-side kernel norm
  normK_star <- integrate(phi_fun, 0, 1)$value
  C_K = (K_Epa(1)^2*(nu2 - nu1)^2 + 2*nu2^2*K_Epa(0)^2)/(nu0*nu2 - nu1^2)
  B_K = (2*log(1/bw))^(1/2) + (2*log(1/bw))^(-1/2)*(log(log(1/bw))/2 + log(C_K/(2*sqrt(pi)*normK_star)))
  ## density f_T is uniform
  df <- density(runif(10000), kernel = "rectangular")
  f.T <- approx(df$x, df$y, xout=regGrid)$y
  ###
  if (det(smoothedCov) < 1e-10) {
    Omega <- 2*((sum(mi*win$wi^2)/bw)*normK_star * (diag(fittedCov) + sigma2)/f.T +
                  sum(mi*(mi-1)*win$wi^2)*diag(fittedCov))
  } else {
    Omega <- 2*((sum(mi*win$wi^2)/bw)*normK_star * (diag(smoothedCov) + sigma2)/f.T +
                  sum(mi*(mi-1)*win$wi^2)*diag(smoothedCov))
  }
  ## normal distribution
  zeta = stats::qnorm(1-alpha / 2) *cutoff(sqrt(abs(Omega)))
  ## extreme value distribution
  # zeta = (B_K - log(log((1-alpha)^(-1/2)))/(2*log(1/bw))^(1/2))*cutoff(sqrt(Omega))

  return(zeta)
}



######################## k-fold cross validation to select h_tau and h_d
CVbandwidth <- function(bw.seq = NULL, zeta = NULL, Ly, Lt, npoly, nder, optns,
                        NbGrid, Wtype, alpha, refined, kFolds, cutoff){

  dataType = optns$dataType
  # If 'Ly' and 'Lt' are vectors "cheat" and break them in a list of 10 elements
  if ( is.vector(Ly) && is.vector(Lt) && !is.list(Lt) && !is.list(Ly) ){
    if (length(Lt) < 21) {
      stop("You are trying to use a local linear weight smoother in a vector with less than 21 values.\n")
    }
    myPartition =   c(1:10, sample(10, length(t)-10, replace=TRUE));
    Ly = split(Ly, myPartition)
    Lt = split(Lt, myPartition)
    dataType = 'Sparse';
  }

  kernel = optns$kernel

  # Make everything into vectors
  ncohort = length(Lt);
  tt  = unlist(Lt);
  yy  = unlist(Ly);
  ind = unlist(lapply( 1:ncohort, function(j) rep(j, times=length(Lt[[j]]))));
  yyn = yy[order(tt)];
  ind = ind[order(tt)];
  ttn = sort(tt);



  if (is.null(bw.seq)) {

    # Get minimum reasonable bandwidth
    a0=ttn[1];
    b0=ttn[length(ttn)];
    rang = b0-a0;
    m_max = max(sapply(Lt, length))
    # nn = ifelse(m_max/ncohort > 2, npoly + m_max/ncohort - 1, npoly + m_max/ncohort)
    nn = npoly + m_max/ncohort
    dstar = Minb(tt, nn); #
    ## data driven bandwidth selection, dense case intend to have smaller minimal bandwidth
    if (dataType != 'Dense'){
      h0 = 2.5*dstar;
    } else {
      h0 = dstar;
    }
    if (h0 > rang/4){
      h0 = h0*.75;
      warning(sprintf("Warning: the min bandwith choice is too big, reduce to %f !", (h0)  ))
    }


    # Get the candidate bandwidths
    nbw = 7;
    bw = rep(0,nbw-1);
    for (i in 1:(nbw-1)){
      bw[i]= 2.5*rang/ncohort*(ncohort/10)^((i-1)/(nbw-1));
    }
    bw.seq = bw-min(bw)+h0;


    # ###
    # # r = diff(range(t))
    # m_max = max(sapply(Lt, length))
    # N = length(tt);
    # r = tt[N] - tt[1];
    # # Specify the starting bandwidth candidates
    # if ( dataType == "Sparse") {
    #   dstar = fdapace:::Minb(tt, npoly+2);
    #   if ( dstar > r*0.25){
    #     dstar = dstar * 0.75;
    #     warning( c( "The min bandwidth choice is too big, reduce to ", dstar, "!\n"))
    #   }
    #   h0 = 4.5 * dstar;
    # }else if(dataType == "DenseWithMV"){
    #   h0 = 4 * fdapace:::Minb(tt, npoly + 1);
    # } else {
    #   h0 = 3.5 * fdapace:::Minb(tt,npoly + 1);
    # }
    # if ( is.nan(h0) ){
    #   if ( kernel == "gauss" ){
    #     h0 = 0.2 * r;
    #   }else{
    #     stop("The data is too sparse, no suitable bandwidth can be found! Try Gaussian kernel instead!\n")
    #   }
    # }
    # h0 <- min(h0,r)
    # q = (r/(4*h0))^(1/9);
    # bw.seq = sort(q^(0:9)*h0);
    # nbw = length(bw.seq)


  } else {

    nbw = length(bw.seq)
  }



  cv = array(Inf, dim = c(length(bw.seq), length(bw.seq), kFolds));
  #count = c();
  theFolds =  CreateFolds(unique(ind), k= kFolds)
  # theFolds =  caret::createFolds(unique(ind), k= kFolds)
  zeta.seq = rep(0, nbw-1)

  for (j in 1:(nbw-1)){

    if (is.null(zeta)) {
      ## get cut-off value from the hypothesis test
      zeta.seq[j] = zetaFun(bw = bw.seq[j], alpha = alpha, Lt = Lt, Ly = Ly, optns = optns,
                            cutoff = cutoff, Wtype = Wtype)
    } else {
      zeta.seq[j] = zeta
    }


    cat('CV procedure: h_tau = ', bw.seq[j], ';', 'zeta = ', zeta.seq[j], '\n')
    # for (k in 1:(nbw-1))
    for (k in (j):(nbw-1)) {

      cat('and h_d = ', bw.seq[k], '\n')

      for (i in 1:kFolds){

        xout= ttn[ ind %in% theFolds[[i]]];
        obs = yyn[ ind %in% theFolds[[i]]];
        xin = ttn[!ind %in% theFolds[[i]]];
        yin = yyn[!ind %in% theFolds[[i]]];
        win=rep(1/length(yin),length(yin));

        muout = tryCatch(
          MeanBreaksFP(xin = xin, yin= yin, weight = win, xout = xout, NbGrid = NbGrid,
                       h_tau = bw.seq[j], h_d = bw.seq[k], zeta = zeta.seq[j],
                       npoly=npoly, nder= nder, kernel = kernel, refined = refined)$muout,
          error=function(err) {
            warning('Invalid bandwidth during CV. Try enlarging the window size.')
            return(Inf)
          })

        cv[j,k,i] = sum((obs - muout)^2)
        # cv[j,k,i] = trapzRcpp(xout, (obs - muout)^2)
        # print(cv)
        if(is.na(cv[j,k,i])){
          cv[j,k,i] = Inf;
        }
      }
    }
  }
  #cv = cv[(count/ncohort>0.90)];
  #bw = bw[(count/ncohort>0.90)];
  if(min(cv) == Inf){
    stop("All bandwidths resulted in infinite CV costs.")
  }

  cvMean = apply(cv, c(1,2), mean)
  cvMeanid = which(cvMean == min(cvMean), arr.ind=TRUE)
  if(optns$useBW1SE){
    cvMeanid = which(cvMean < min(cvMean) +
                       apply(cv, c(1,2), sd)[cvMeanid]/sqrt(kFolds), arr.ind=TRUE)
    cvMeanid = cvMeanid[which.max(cvMeanid[,1]),]
  }
  bopt = bw.seq[cvMeanid];
  names(bopt) = c('h_tau', 'h_d')
  zeta = zeta.seq[cvMeanid[1]]

  boptList <- list('bopt' = bopt, 'zeta' = zeta)

  return(boptList)

}






## Individual Mean function by JIC in Xia and Qiu 2015
indMeanbreak <- function(y, t, M_max, NbGrid = NbGrid, npoly, nder, kernel) {



  n = length(y)
  BIC = sigma2 = rep(Inf, M_max)
  bw = 0.1*n^(-1/5)
  Obj = list()
  for (i in 1:length(BIC)) {

    # Generate basic grids for jump detect:
    obsGrid = sort(unique(t));
    # obsGrid = sort(t);
    if(is.null(NbGrid)){
      jumpGrid = obsGrid;
    } else {
      jumpGrid = seq(min(obsGrid), max(obsGrid),length.out = NbGrid);
    }
    ## Local linear estimate based on one-sided kernel
    D_h <- jumpGrid[bw <= jumpGrid & jumpGrid <= 1- bw];
    weight = rep(1/length(t), length(t))
    mu_est_lr = CPPlwls1d_LR(bw = bw, kernel_type = kernel, win = weight,
                             xin = t, yin = y, xout = D_h, npoly = npoly)
    mu_diff = mu_est_lr[,2]-mu_est_lr[,1]
    mu_diff_abs = abs(mu_diff)
    mu_diff_time = D_h[order(-mu_diff_abs)]
    mu_diff_size = mu_diff[order(-mu_diff_abs)]

    ##
    mu_jumptime=mu_diff_time[1]
    mu_jumpsize=mu_diff_size[1]
    if(i > 1){
      for (k in 2:i) {
        if (sum(abs(mu_diff_time-mu_jumptime[k-1]) > 2*bw) >0){
          index = which(abs(mu_diff_time - mu_jumptime[k-1]) > 2*bw)
          mu_diff_time = mu_diff_time[index]
          mu_diff_size = mu_diff_size[index]
          mu_jumptime=append(mu_jumptime, mu_diff_time[1])
          mu_jumpsize=append(mu_jumpsize, mu_diff_size[1])
        }
        else{
          break
        }
      }
    }


    yy <- y - sapply(t, function(z) sum(mu_jumpsize*(z >= mu_jumptime)))
    # nu: the smoothed mean curve evaluated at times 'xout' with same bandwidth of h_tau
    nu = Lwls1D(bw, kernel_type = kernel, npoly = npoly, nder = nder,
                xin = t, yin= yy, xout = obsGrid, win = weight)
    mu <- nu + sapply(obsGrid, function(z)sum(mu_jumpsize *(z >= mu_jumptime)))

    #
    sigma2[i] = mean((y - mu)^2)

    ## BIC
    gamma = 1
    pn = ifelse( any(mu_jumpsize == 0), 0,  sum(1/abs(mu_jumpsize)^gamma)*(log(n)*bw/n)^(1/2) )
    BIC[i] = log(sigma2[i]) + pn
    Obj[[i]] = list(mu = mu,
                    obsGrid = obsGrid,
                    bw = bw,
                    mu_jumptime = mu_jumptime,
                    mu_jumpsize = mu_jumpsize)

  }
  ret = Obj[[which.min(BIC)]]


  return(ret);

}


