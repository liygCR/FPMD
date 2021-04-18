#' Monthly max average value weighted returns data.
#'
#' The original data set consists of the daily simple returns of n = 49
#' industry portfolios from 1927 to 2020.
#'
#' @docType data
#'
#' @usage data(MMavwreturn)
#'
#' @format
#' \describe{
#'   \item{data_SS}{Data matrix.}
#'   \item{Lt}{Observation time points.}
#'   \item{Ly}{Daily average value weighted returns.}
#'   ...
#' }
#' @source \url{http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html}
#'
#' @examples
#'
#' \dontrun{
#'
#'
#'
#' library(reshape2)
#' library(zoo)
#' library(dplyr)
#'
#' data(MMavwreturn)
#' data_SS = MMavwreturn$data_SS
#' Ly = MMavwreturn$Ly
#' Lt = MMavwreturn$Lt
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
#'
#' ##
#' res = FPMD(Ly = Ly, Lt = Lt, wi = NULL, Wtype = "MIX", zeta = NULL,
#'            bw.seq = NULL, NbGrid = 101, kFolds = 5, refined = TRUE,
#'            individualCurveJump = TRUE, nRegGrid = 501,
#'            npoly = 1, nder = 0, alpha = 0.05, cutoff = max, M_max = 15)
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
#' par(mfrow = c(1, 2))
#' matplot(data_SS[,-50], type = "p", pch = 1, xlab = "Dates",
#'         ylab = "Average Value Weighted Returns",
#'         main = '',
#'         col = 'gray', xaxt="n")
#' ind_x = seq(1, nrow(data_SS), length.out = 10)
#' axis(1, at=ind_x, labels=data_SS[ind_x, 50])
#' abline( v = round(mu_jumptime *nrow(data_SS)), lty = "dashed", col = "red")
#' lines(obsGrid*nrow(data_SS), mu, col = "red", lwd = 1.5)
#' data_SS$ym[round(mu_jumptime *nrow(data_SS))]
#'
#'
#'
#' ## plot mean function and confidence band
#' ## true and estimated mean curve are based on workGrid points
#' plot.fmb <- function(res){
#'
#'   ## confidence band
#'   cbandMu <- FPMD:::pwCBFun(res)
#'   workGrid <- res$workGrid
#'   muWork <- res$muWork
#'   rho = res$rho
#'   tau_est = c(0, res$mu_jumptime, 1)
#'   mudata = t(rbind(cbandMu, muWork))
#'
#'
#'   plot(x = NULL, y = NULL, ylim = c(0, 9), xlim = range(unlist(Lt)),
#'        xlab = "Dates",
#'        ylab = "Average Value Weighted Returns",
#'        main = '',
#'        col = rgb(0.7,0.7,0.7,0.4), xaxt="n")
#'   ind_x = seq(min(Lt[[1]]), max(Lt[[1]]), length.out = 10)
#'   axis(1, at=ind_x, labels= data_SS[ind_x*nrow(data_SS), 50])
#'   ###
#'
#'   bandGrid =  t(rbind(cbandMu, workGrid))
#'   lband = lapply(1:(length(tau_est) - 1), function(i)
#'     as.data.frame(bandGrid[workGrid < (tau_est[i + 1] - rho) &
#'                              workGrid >= (tau_est[i] + rho),]))
#'   lapply(lband, function(x) polygon(c(x$workGrid, rev(x$workGrid)),c(x$lwCI,rev(x$upCI)),
#'                                     col = rgb(0.7,0.7,0.7,0.4) , border = NA))
#'   for (i in 1:(length(tau_est)-1)) {
#'
#'     matplot(workGrid[workGrid < (tau_est[i+1]- rho) & workGrid >= (tau_est[i]+ rho) ],
#'             mudata[workGrid < (tau_est[i+1] - rho) & workGrid >= (tau_est[i]+rho), ],
#'             type = "l", add = TRUE, lty = c(3, 3, 1), col = c(4,4,2), lwd = 2)
#'
#'   }
#'   legend('topleft', legend = c( 'Estimated mean', 'Pointwise confidence interval'),
#'          lwd = 2, col = c(2,4), lty = c(1,3), bty = "n")
#'   points(x = res$mu_jumptime, y = rep(-.35, length(res$mu_jumptime)),
#'          pch= rep("*", length(res$mu_jumptime)), cex = 2, col = 4,  xpd = TRUE)
#'   # lines(res$obsGrid, res$mu, col = 4)
#'   text(x = res$mu_jumptime, y = rep(-.3, length(res$mu_jumptime)),
#'        labels = data_SS$ym[round(res$mu_jumptime*nrow(data_SS))],
#'        xpd = TRUE, pos = 3, cex = 0.5, col = 4)
#'
#' }
#'
#' plot.fmb(res)
#'
#' ## individual
#' ####
#' industry = colnames(data_SS)[-50]
#' par(mfrow=c(2,2))
#' for (i in c(5,16,35, 36)) {
#'   plot(data_SS[, i],
#'        xlab = "Dates",
#'        ylab = "Average Value Weighted Returns",
#'        main = industry[i],
#'        col = 'gray', #rgb(0.7,0.7,0.7,0.4),
#'        xaxt="n", ylim = c(0, max(data_SS[, i])))
#'   ind_x = seq(1, nrow(data_SS), length.out = 10)
#'   axis(1, at=ind_x, labels=data_SS[ind_x, 50])
#'   ## indivudual
#'   lines(y = res$indJump[[i]]$mu, x = Lt[[i]]*nrow(data_SS),
#'         col = "blue", lty = 1, lwd = 1.5 )
#'   points(x = res$indJump[[i]]$mu_jumptime*nrow(data_SS),
#'          y = rep(-0.55, length(res$indJump[[i]]$mu_jumptime)),
#'          pch= rep("*", length(res$indJump[[i]]$mu_jumptime)),
#'          cex = 2, col = "blue",  xpd = TRUE)
#'   text(x = res$indJump[[i]]$mu_jumptime*nrow(data_SS),
#'        y = rep(-0.55, length(res$indJump[[i]]$mu_jumptime)),
#'        labels = data_SS$ym[res$indJump[[i]]$mu_jumptime*nrow(data_SS)],
#'        xpd = TRUE, pos = 3, cex = 0.5, col = "blue")
#'   ## our
#'   lines(obsGrid*nrow(data_SS), mu, col = "red", lty = 1, lwd = 1.5)
#'   abline( v = round(mu_jumptime *nrow(data_SS)), col = "red", lty = "dashed")
#' }
#'
#'
#' }
#'
#'
"MMavwreturn"
