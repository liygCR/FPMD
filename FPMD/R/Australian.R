#' Sydney Temperature Data.
#'
#' Australian daily minimum temperature climate data form year 1855 to 2012 for
#' Sydney (Observatory Hill) station.
#'
#' @docType data
#'
#' @usage data(Australian)
#'
#' @format
#' \describe{
#'   \item{dataSS}{Data matrix.}
#'   \item{Lt}{Observation time points.}
#'   \item{Ly}{Daily minimum temperature.}
#'   ...
#' }
#' @source \url{http://www.bom.gov.au/climate/data}
#'
#' @examples
#'
#' \dontrun{
#' data(Australian)
#' dataSS = Australian$dataSS
#' Ly = Australian$Ly
#' Lt = Australian$Lt
#'
#' ## global test
#' bw.seq = seq(0.01, 0.11, by = 0.02)
#' reject = rep(NA, length(bw.seq))
#' ## reject H0: no jump points
#' for (i in 1:length(bw.seq)) {
#'   reject[i] = MDTest(bw = bw.seq[i], alpha = 0.05, Lt = Lt, Ly = Ly, Wtype = 'MIX')
#'   cat('reject:', reject[i], 'for h = ', bw.seq[i], '\n')
#' }
#'
#' ###
#' res = FMbreaks(Ly = Ly, Lt = Lt, wi = NULL, Wtype = "MIX", zeta = NULL,
#'                bw.seq = bw.seq, NbGrid = 101, kFolds = 5, refined = TRUE,
#'                individualCurveJump = TRUE, useBW1SE = FALSE, nRegGrid = 51,
#'                npoly = 1, nder = 0, alpha = 0.05, cutoff = max, M_max = 5)
#'
#' ## with change points
#' mu_jumptime <- res$mu_jumptime
#' mu_jumpsize <- res$mu_jumpsize
#' # mu_jumpsize <- res$mu_jumpsize_h_tau
#' h_tau <- res$h_tau
#' h_d <- res$h_d
#' zeta <- res$zeta
#' wi <- res$wi
#' mu <- res$mu
#' muWork <- res$muWork
#' obsGrid <- res$obsGrid
#' workGrid <- res$workGrid
#'
#'
#' ###
#' par(mfrow = c(1, 2))
#' matplot(dataSS[,100:154], type = "p", pch = 1, col = 'gray', xlab = "Day",
#'         ylab = "Minimum Temperature", main= "", xaxt="n")
#' ind_x = seq.int(1, nrow(dataSS), length.out = 6)
#' axis(1, at=ind_x, labels=dataSS$ymd[ind_x])
#' lines(obsGrid*366, mu, col = "red", lwd = 1.5)
#' abline( v = res$mu_jumptime*366, lty = "dashed", col = "red" )
#' dataSS$ymd[round(res$mu_jumptime*366)]
#'
#'
#' }
#'
#'
"Australian"