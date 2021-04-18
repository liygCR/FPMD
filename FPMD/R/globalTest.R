#'  A global test for jump points in mean function
#'
#' @param Ly A list of \emph{n} vectors containing the observed values for each individual. Missing values specified by \code{NA}s are supported for dense case (\code{dataType='dense'}).
#' @param Lt A list of \emph{n} vectors containing the observation time points for each individual corresponding to y. Each vector should be sorted in ascending order.
#' @param Wtype Weight structure specified for estimation. It works if \code{wi = NULL}.
#' @param bw The selected (or user specified) bandwidth.
#' @param nRegGrid The number of support points in each direction of covariance surface; numeric - default: 101.
#' @param alpha The confidence level for the cut-off threshold. The default is
#'              \code{alpha = 0.05}, and for heavy tailed process, we suggest a higher level \code{0.1}.
#'
#'
#' @return
#' \item{reject}{1 indicates that there exists discontinuities.}
#'
#' @export


## global test for jump points
MDTest <- function(bw, alpha = 0.05, Lt, Ly, nRegGrid = 101,
                     Wtype = c("OBS", "SUBJ", "MIX", "OPT")){


  Wtype = match.arg(Wtype)
  optns <- list(methodBwMu = 'CV',
                nRegGrid = nRegGrid,
                error = TRUE,
                kernel='epan',
                verbose=FALSE)
  optns <- SetOptions(Ly, Lt, optns);

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
  ## test statistics
  mu_lr = CPPlwls1d_LR(bw = bw, kernel_type = optns$kernel, win = win$weight,
                       xin = xin, yin = yin, xout = workGrid, npoly = 1)
  mu_diff = mu_lr[,2]-mu_lr[,1]
  mu_diff_abs = abs(mu_diff)
  ##
  reject = ifelse(max(mu_diff_abs/sqrt(Omega)) >
                    B_K - log(log((1-alpha)^(-1/2)))/(2*log(1/bw))^(1/2), 1 , 0)
  cat('Cut = ', (B_K - log(log((1-alpha)^(-1/2)))/(2*log(1/bw))^(1/2))*min(sqrt(Omega)), '\n')
  return(reject)
}



