#' Create a Dense Functional Data sample for a Gaussian process
#'
#' For a Gaussian or t process, create a dense or sparse functional data sample of size n over a [0,1] support.
#'
#'
#'
#' @param n number of samples to generate.
#' @param Q an integral. The observations for per sample are following the poisson distribution with mean Q.
#' @param rdist a sampler for generating the random design time points within [0, 1].
#' @param muFun a function that takes a vector input and output a vector of the corresponding mean (default: zero function).
#' @param K scalar specifying the number of basis to be used (default: 3).
#' @param lambda vector of size K specifying the variance of each components (default: rep(1,K)).
#' @param sigma The standard deviation of the Gaussian noise added to each observation points.
#' @param basisType string specifying the basis type used; possible options are: 'sin', 'cos' and 'fourier' (default: 'cos') (See code of 'CreateBasis' for implementation details.)
#' @param proc The type of random process, 'gaussian' or 't' process.
#'
#'
#' @return A list containing the generated data with \code{List} type:
#' \item{Ly}{A list of \emph{n} vectors containing the observed values for each individual. Missing values specified by \code{NA}s are supported for dense case (\code{dataType='dense'}).}
#' \item{Lt}{A list of \emph{n} vectors containing the observation time points for each individual corresponding to y. Each vector should be sorted in ascending order.}
#' \item{mi}{Number of observations for each subject.}
#'
#'
#' @export
#'

## data generation
MakeFPMBData <- function(n, Q, muFun, rdist=runif, K, lambda, sigma,
                         basisType = "cos", proc = c("gaussian", "t")){
  if(n < 2){
    stop("Samples of size 1 are irrelevant.")
  }
  if(!is.function(rdist)){
    stop("'rdist' needs to be a function.")
  }
  if(!is.function(muFun)){
    stop("'muFun' needs to be a function.")
  }
  if (!is.numeric(sigma) || sigma < 0) {
    stop("'sigma' needs to be a nonnegative number")
  }
  if( !(basisType %in% c('cos','sin','fourier'))) {
    stop("Make sure you provide a valid parametric basis.")
  }

  proc = match.arg(proc)

  if (length(Q) == 1) {
    mi = rtpois(n, Q, a = 2) #rpois(n,Q)
  } else {
    mi = sample(Q, n, replace=TRUE)
  }

  if (proc == "gaussian") {
    ## lambda = 1/(seq_len(K) + 1)^2
    Ksi <- apply(matrix(rnorm(n*K), ncol=K), 2, scale) %*%
      diag(sqrt(lambda), length(lambda))
  } else {
    ## lambda = 1/(seq_len(K) + 1)
    Ksi <- apply(matrix(rt(n*K, df = K), ncol=K), 2, scale) %*%
      diag(sqrt(lambda), length(lambda))

  }


  samp <- lapply(seq_len(n), function(i) {
    ni <- mi[i]
    ti <- sort(rdist(ni))
    Phii <- CreateBasis(K, ti, basisType)

    # Phii <- sapply(seq_len(K), function(k) if (k = 1) {
    #   rep(1, length(ti))
    # }
    # else if (k%%2 == 0) {
    #   sqrt(2) * sin(k * pi * ti)
    # }
    # else {
    #   sqrt(2) * cos((k - 1) * pi * ti)
    # })
    # Phii <- matrix(Phii, ncol = 4)

    yi <- muFun(ti) + as.numeric(tcrossprod(Ksi[i, ], Phii))

    list(ti = ti, yi=yi)
  })

  Lt <- lapply(samp, `[[`, 'ti')
  Ly <- lapply(samp, `[[`, 'yi')


  if (sigma > 0) {
    LyTrue <- Ly
    if (proc == "gaussian") {
      Ly <- lapply(LyTrue, function(x) x + rnorm(length(x), sd=sigma))
    } else {
      Ly <- lapply(LyTrue, function(x) x + 0.1*rt(length(x), df = K))
    }

  }

  ## weight
  # SUBJ
  # wi = 1/(n*mi)
  ## OBJ
  # wi <- rep(1/sum(mi), n)
  ## mix weight
  # h = 0.1 #ifelse(m > n/log(n), 0.075, 0.1)
  # c1 <- ( 1/(h*mean(mi)) + mean(mi^2)/mean(mi)^2 )/n
  # c2 <- (mean(1/mi)/h + 1)/n
  # alp <- c2/(c1+c2)
  # wi <- alp/sum(mi) + (1-alp)/(n*mi)

  res <- list(Ly=Ly, Lt=Lt, xi=Ksi, mi=mi)

  if (sigma > 0) {
    res <- append(res, list(LyTrue=LyTrue))
  }



  return(res)

}



# ## plot mean function and confidence band
# ## true and estimated mean curve are based on workGrid points
# plot.fmb <- function(res, tau, mu_fun){
#
#   ## confidence band
#   cbandMu <- pwCBFun(res)
#   workGrid <- res$workGrid
#   muWork <- res$muWork
#   h_tau <- res$h_tau
#   rho = res$rho
#   # tau_est = c(h_tau, res$mu_jumptime_refine, 1-h_tau)
#   # tau_true = c(h_tau, tau, 1-h_tau)
#   ##
#   tau_est = c(0, res$mu_jumptime, 1)
#   tau_true = c(0, tau, 1)
#
#
#   mudata = t(rbind(cbandMu, muWork))
#
#   plot(NA, ylab = expression(paste(mu,"(t)")), ylim = range(-1.5,3), xlim = c(0, 1), type = "n")
#   sapply(1:(length(tau_true)-1),
#          function(i)lines(workGrid[workGrid < (tau_est[i + 1] - rho) & workGrid >= (tau_est[i] + rho)],
#                           sapply(workGrid, mu_fun)[workGrid < (tau_est[i + 1] - rho) & workGrid >= (tau_est[i] + rho)],
#                           lty = 1, col = 4, lwd = 1.5))
#   ###
#   bandGrid =  t(rbind(cbandMu, workGrid))
#   lband = lapply(1:(length(tau_est) - 1), function(i)
#     as.data.frame(bandGrid[workGrid < (tau_est[i + 1] - rho) &
#                              workGrid >= (tau_est[i] + rho),]))
#   lapply(lband, function(x) polygon(c(x$workGrid, rev(x$workGrid)),c(x$lwCI,rev(x$upCI)),
#                                     col = rgb(0.7,0.7,0.7,0.4) , border = NA))
#   for (i in 1:(length(tau_est)-1)) {
#
#     matplot(workGrid[workGrid < (tau_est[i+1]- rho) & workGrid >= (tau_est[i]+ rho) ],
#             mudata[workGrid < (tau_est[i+1] - rho) & workGrid >= (tau_est[i]+rho), ],
#             type = "l", add = TRUE, lty = c(3, 3, 1), col = c(3,3,2), lwd = 2)
#
#   }
#   legend('topleft', legend = c('True mean', 'Estimated mean', 'Pointwise confidence interval'),
#          cex = .8, lwd = 2, col = c(4,2,3), lty = c(1,1,3), bty = "n")
#   points(x = tau, y = rep(-1.7, length(tau)), pch= rep("|", 3), col = c(4,4,4),  xpd = TRUE)
#   text(x = tau, y = rep(-1.7, length(tau)), labels = round(tau,2), xpd = TRUE, pos = 3, cex = 0.8, col = 4)
#
# }

# # setting 1
# ## 3 change points
# tau = c(.25,.5,.75) # jump locations
# disc = c(.5,-.4,.3) # jump sizes
# smoothmean=function(x) x^2+sin(2*pi*x) + cos(2*pi*x)
# jump = function(x) sum((x>=tau)*disc) # jump function
#
#
# # setting 2
# ## 5 change point
# tau= seq(0,1, length.out = 7)[2:6] # jump locations
# disc = c(.5,-.4,.3, -.4, .5) # jump sizes
# smoothmean=function (x) sin(2*pi*x) # smooth function
# jump = function(x) sum((x>=tau)*disc) # jump function


