#' Average value weighted returns data.
#'
#' The original data set consists of the daily simple returns of n = 49
#' industry portfolios from 1927 to 2020.
#'
#' @docType data
#'
#' @usage data(avwreturn)
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
#' data(avwreturn)
#' data_SS = avwreturn$data_SS
#' Ly = avwreturn$Ly
#' Lt = avwreturn$Lt
#'
#'
#' ## global test
#' bw.seq = seq(0.01, 0.11, by = 0.01)
#' reject = rep(NA, length(bw.seq))
#' ## reject H0: no jump points
#' for (i in 1:length(bw.seq)) {
#'   reject[i] = MDTest(bw = bw.seq[i], alpha = 0.05, nRegGrid = 201, Lt = Lt, Ly = Ly, Wtype = 'MIX')
#'   cat('reject:', reject[i], 'for h = ', bw.seq[i], '\n')
#' }
#'
#' ###
#' res = FPMD(Ly = Ly, Lt = Lt, wi = NULL, Wtype = "MIX", zeta = NULL,
#'            bw.seq = NULL, NbGrid = 101, kFolds = 5, refined = TRUE,
#'            individualCurveJump = TRUE, nRegGrid = 301,
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
#' par(mfrow=c(1,1))
#' matplot(data_SS[,-1], type = "p", pch = 1, xlab = "Dates",
#'         ylab = "Average Value Weighted Returns",
#'         main = '',
#'         col = rgb(0.7,0.7,0.7,0.4), xaxt="n")
#' ind_x = seq(1, nrow(data_SS), length.out = 18)
#' axis(1, at=ind_x, labels=data_SS[ind_x, 1])
#' abline( v = round(mu_jumptime *nrow(data_SS)), lty = 2)
#' lines((obsGrid*nrow(data_SS)), mu, col = 4)
#' data_SS$Date[round(mu_jumptime*nrow(data_SS))]
#'
#'
#' ## combine all as iid, Xia and Qiu 2016
#' t = unlist(Lt);
#' y = unlist(Ly)[order(t)];
#' t = sort(t);
#' indJumpall = FPMD:::indMeanbreak(y = y, t = t, M_max = 15, NbGrid = 101,
#'                           kernel = res$optns$kernel, npoly = 1, nder = 0)
#' ##  hdbinseg package from H. Cho and P. Fryzlewicz (2014) JRSSB
#' library(hdbinseg)
#' dd = t(data_SS[,-1])
#' ecp_CHO = dcbs.alg(dd, cp.type=1, phi=0.5, temporal=FALSE, do.parallel=0)$ecp
#'
#'
#'
#'
#'
#' ###
#' par(mfrow = c(2, 2))
#' ind_x = seq(1, nrow(data_SS), length.out = 18)
#' ## 1
#' matplot(data_SS[1:100,-1], type = "p", pch = 1, xlab = "Dates",
#'         ylab = "Average Value Weighted Returns",
#'         main = '', col = 'gray', xaxt="n")
#' axis(1, at=ind_x[1:5], labels=data_SS[ind_x[1:5], 1])
#' abline( v = round(mu_jumptime *nrow(data_SS)), col = 'red', lty = 2)
#' lines((obsGrid*nrow(data_SS))[1:100], mu[1:100], col = 'red', lwd = 1.5)
#' abline( v = round(indJumpall$mu_jumptime*nrow(data_SS)),
#'         lty = 'dotted', col = "green")
#' lines(y = indJumpall$mu[1:100], x = (obsGrid*nrow(data_SS))[1:100],
#'       col = "green", lty = 1, lwd = 1.5)
#' points(x = round(indJumpall$mu_jumptime*nrow(data_SS))[1:2],
#'        y = rep(-8.3, length(indJumpall$mu_jumptime))[1:2],
#'        pch= rep("*", length(indJumpall$mu_jumptime)),
#'        cex = 2, col = "green",  xpd = TRUE)
#' text(x = round(indJumpall$mu_jumptime*nrow(data_SS))[1:2],
#'      y = rep(-8.3, length(indJumpall$mu_jumptime))[1:2],
#'      labels = data_SS$Date[round(indJumpall$mu_jumptime*nrow(data_SS))][1:2],
#'      xpd = TRUE, pos = 1, cex = 0.5, col = "green")
#' ## cho
#' points(x = ecp_CHO,
#'        y = rep(-8.3, length(ecp_CHO)),
#'        pch= rep("*", length(ecp_CHO)),
#'        cex = 2, col = "blue",  xpd = TRUE)
#' text(x = ecp_CHO,
#'      y = rep(-8.3, length(ecp_CHO)),
#'      labels = data_SS$Date[ecp_CHO],
#'      xpd = TRUE, pos = 3, cex = 0.5, col = "blue")
#' ## 2
#' matplot(data_SS[101:160,-1], type = "p", pch = 1, xlab = "Dates",
#'         ylab = "Average Value Weighted Returns",
#'         main = '', col = 'gray', xaxt="n")
#' axis(1, at=ind_x[6:8]-100, labels=data_SS[ind_x[6:8], 1])
#' abline( v = round(mu_jumptime *nrow(data_SS))[-c(1:3)]-100, col = 'red', lty = 2)
#' lines((obsGrid*nrow(data_SS))[1:60], mu[101:160], col = 'red', lwd = 1.5)
#' lines(y = indJumpall$mu[101:160], x = (unique(t)*nrow(data_SS))[1:60],
#'       col = "green", lty = 1, lwd = 1.5)
#' points(x = round(indJumpall$mu_jumptime*nrow(data_SS))-100,
#'        y = rep(-6, length(indJumpall$mu_jumptime)),
#'        pch= rep("*", length(indJumpall$mu_jumptime)),
#'        cex = 2, col = "green",  xpd = TRUE)
#' text(x = round(indJumpall$mu_jumptime*nrow(data_SS))-100,
#'      y = rep(-6, length(indJumpall$mu_jumptime)),
#'      labels = data_SS$Date[round(indJumpall$mu_jumptime*nrow(data_SS))],
#'      xpd = TRUE, pos = 1, cex = 0.5, col = "green")
#'
#' ## 3
#' matplot(data_SS[161:280,-1], type = "p", pch = 1, xlab = "Dates",
#'         ylab = "Average Value Weighted Returns", ylim = c(-8, 8),
#'         main = '', col = 'gray', xaxt="n")
#' axis(1, at=ind_x[9:13]-160, labels=data_SS[ind_x[9:13], 1])
#' abline( v = round(mu_jumptime *nrow(data_SS))[-c(1:4)]-160, col = 'red', lty = 2)
#' lines((obsGrid*nrow(data_SS))[1:120], mu[161:280], col = 'red', lwd = 1.5)
#' lines(y = indJumpall$mu[161:280], x = (unique(t)*nrow(data_SS))[1:120],
#'       col = "green", lty = 1, lwd = 1.5)
#' points(x = round(indJumpall$mu_jumptime*nrow(data_SS))-160,
#'        y = rep(-8.6, length(indJumpall$mu_jumptime)),
#'        pch= rep("*", length(indJumpall$mu_jumptime)),
#'        cex = 2, col = "green",  xpd = TRUE)
#' text(x = round(indJumpall$mu_jumptime*nrow(data_SS))-160,
#'      y = rep(-8.6, length(indJumpall$mu_jumptime)),
#'      labels = data_SS$Date[round(indJumpall$mu_jumptime*nrow(data_SS))],
#'      xpd = TRUE, pos = 1, cex = 0.5, col = "green")
#'
#' ## 4
#' matplot(data_SS[281:354,-1], type = "p", pch = 1, xlab = "Dates",
#'         ylab = "Average Value Weighted Returns",
#'         main = '', col = 'gray', xaxt="n")
#' axis(1, at=ind_x[14:18]-280, labels=data_SS[ind_x[14:18], 1])
#' abline( v = round(mu_jumptime *nrow(data_SS))[-c(1:7)]-280, col = 'red', lty = 2)
#' lines((obsGrid*nrow(data_SS))[1:74], mu[281:354], col = 'red', lwd = 1.5)
#' lines(y = indJumpall$mu[281:354], x = (unique(t)*nrow(data_SS))[1:74],
#'       col = "green", lty = 1, lwd = 1.5)
#' abline( v = round(indJumpall$mu_jumptime*nrow(data_SS))-280,
#'         lty = 'dotted', col = "green")
#' points(x = indJumpall$mu_jumptime*nrow(data_SS) - 280,
#'        y = rep(-24.4, length(indJumpall$mu_jumptime)),
#'        pch= rep("*", length(indJumpall$mu_jumptime)),
#'        cex = 2, col = "green",  xpd = TRUE)
#' text(x = indJumpall$mu_jumptime*nrow(data_SS) - 280,
#'      y = rep(-24.4, length(indJumpall$mu_jumptime)),
#'      labels = data_SS$Date[indJumpall$mu_jumptime*nrow(data_SS)],
#'      xpd = TRUE, pos = 1, cex = 0.5, col = "green")
#' ## cho
#' points(x = ecp_CHO - 280,
#'        y = rep(-24.4, length(ecp_CHO)),
#'        pch= rep("*", length(ecp_CHO)),
#'        cex = 2, col = "blue",  xpd = TRUE)
#' text(x = ecp_CHO - 280,
#'      y = rep(-24.4, length(ecp_CHO)),
#'      labels = data_SS$Date[ecp_CHO],
#'      xpd = TRUE, pos = 3, cex = 0.5, col = "blue")
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
#'   ###
#'   bandGrid =  t(rbind(cbandMu, workGrid))
#'   lband = lapply(1:(length(tau_est) - 1), function(i)
#'     as.data.frame(bandGrid[workGrid < (tau_est[i + 1] - rho) &
#'                              workGrid >= (tau_est[i] + rho),]))
#'
#'   par(mfrow = c(2, 2))
#'   ind_x = seq(1, nrow(data_SS), length.out = 18)
#'   ## 1
#'   plot(x = NULL, y = NULL, ylim = c(-3, 3), xlim = c(0, 0.26),
#'        xlab = "Dates",
#'        ylab = "Average Value Weighted Returns",
#'        main = '',
#'        col = rgb(0.7,0.7,0.7,0.4), xaxt="n")
#'   axis(1, at= (ind_x/nrow(data_SS))[1:5], labels=data_SS[ind_x[1:5], 1])
#'
#'   lapply(lband[1:5], function(x) polygon(c(x$workGrid, rev(x$workGrid)),c(x$lwCI,rev(x$upCI)),
#'                                     col = rgb(0.7,0.7,0.7,0.4) , border = NA))
#'   for (i in 1:5) {
#'     matplot(workGrid[workGrid < (tau_est[i+1]- rho) & workGrid >= (tau_est[i]+ rho) ],
#'             mudata[workGrid < (tau_est[i+1] - rho) & workGrid >= (tau_est[i]+rho), ],
#'             type = "l", add = TRUE, lty = c(3, 3, 1), col = c(4,4,2), lwd = 2)
#'
#'   }
#'   legend('top', legend = c( 'Estimated mean', 'Pointwise confidence interval'),
#'          cex = 1, lwd = 2, col = c(2,4), lty = c(1,3), bty = "n")
#'   points(x = res$mu_jumptime[1:4], y = rep(-3.2, length(res$mu_jumptime[1:4])),
#'          pch= rep("*", length(res$mu_jumptime[1:4])), cex = 2, col = 4,  xpd = TRUE)
#'   # lines(res$obsGrid, res$mu, col = 4)
#'   text(x = res$mu_jumptime[1:4], y = rep(-3.2, length(res$mu_jumptime[1:4])),
#'        labels = data_SS$Date[round(res$mu_jumptime[1:4]*nrow(data_SS))],
#'        xpd = TRUE, pos = 3, cex = 0.8, col = 4)
#'
#'   ## 2
#'   plot(x = NULL, y = NULL, ylim = c(-1.5, 1.5), xlim = c(0.26, 0.5),
#'        xlab = "Dates",
#'        ylab = "Average Value Weighted Returns",
#'        main = '',
#'        col = rgb(0.7,0.7,0.7,0.4), xaxt="n")
#'   axis(1, at= (ind_x/nrow(data_SS))[5:8], labels=data_SS[ind_x[5:8], 1])
#'   ###
#'   lapply(lband[5:8], function(x) polygon(c(x$workGrid, rev(x$workGrid)),c(x$lwCI,rev(x$upCI)),
#'                                     col = rgb(0.7,0.7,0.7,0.4) , border = NA))
#'   for (i in 5:8) {
#'     matplot(workGrid[workGrid < (tau_est[i+1]- rho) & workGrid >= (tau_est[i]+ rho) ],
#'             mudata[workGrid < (tau_est[i+1] - rho) & workGrid >= (tau_est[i]+rho), ],
#'             type = "l", add = TRUE, lty = c(3, 3, 1), col = c(4,4,2), lwd = 2)
#'   }
#'   legend('top', legend = c( 'Estimated mean', 'Pointwise confidence interval'),
#'          cex = 1, lwd = 2, col = c(2,4), lty = c(1,3), bty = "n")
#'   points(x = res$mu_jumptime[5:7], y = rep(-1.6, length(res$mu_jumptime[5:7])),
#'          pch= rep("*", length(res$mu_jumptime[5:7])), cex = 2, col = 4,  xpd = TRUE)
#'   # lines(res$obsGrid, res$mu, col = 4)
#'   text(x = res$mu_jumptime[5:7], y = rep(-1.6, length(res$mu_jumptime[5:7])),
#'        labels = data_SS$Date[round(res$mu_jumptime[5:7]*nrow(data_SS))],
#'        xpd = TRUE, pos = 3, cex = 0.8, col = 4)
#'
#'
#'   ## 3
#'   plot(x = NULL, y = NULL, ylim = c(-1.5, 1.5), xlim = c(0.5, 0.8),
#'        xlab = "Dates",
#'        ylab = "Average Value Weighted Returns",
#'        main = '',
#'        col = rgb(0.7,0.7,0.7,0.4), xaxt="n")
#'   axis(1, at= (ind_x/nrow(data_SS))[9:13], labels=data_SS[ind_x[9:13], 1])
#'   ###
#'   ###
#'   lapply(lband[8:12], function(x) polygon(c(x$workGrid, rev(x$workGrid)),c(x$lwCI,rev(x$upCI)),
#'                                          col = rgb(0.7,0.7,0.7,0.4) , border = NA))
#'   for (i in 8:12) {
#'     matplot(workGrid[workGrid < (tau_est[i+1]- rho) & workGrid >= (tau_est[i]+ rho) ],
#'             mudata[workGrid < (tau_est[i+1] - rho) & workGrid >= (tau_est[i]+rho), ],
#'             type = "l", add = TRUE, lty = c(3, 3, 1), col = c(4,4,2), lwd = 2)
#'   }
#'   legend('top', legend = c( 'Estimated mean', 'Pointwise confidence interval'),
#'          cex = 1, lwd = 2, col = c(2,4), lty = c(1,3), bty = "n")
#'   points(x = res$mu_jumptime[8:11], y = rep(-1.6, length(res$mu_jumptime[8:11])),
#'          pch= rep("*", length(res$mu_jumptime[8:11])), cex = 2, col = 4,  xpd = TRUE)
#'   # lines(res$obsGrid, res$mu, col = 4)
#'   text(x = res$mu_jumptime[8:11], y = rep(-1.6, length(res$mu_jumptime[8:11])),
#'        labels = data_SS$Date[round(res$mu_jumptime[8:11]*nrow(data_SS))],
#'        xpd = TRUE, pos = 3, cex = 0.7, col = 4)
#'
#'   ## 4
#'   plot(x = NULL, y = NULL, ylim = c(-5, 5), xlim = c(0.8, 1),
#'        xlab = "Dates",
#'        ylab = "Average Value Weighted Returns",
#'        main = '',
#'        col = rgb(0.7,0.7,0.7,0.4), xaxt="n")
#'   axis(1, at= (ind_x/nrow(data_SS))[14:18], labels=data_SS[ind_x[14:18], 1])
#'   lapply(lband[12:14], function(x) polygon(c(x$workGrid, rev(x$workGrid)),c(x$lwCI,rev(x$upCI)),
#'                                           col = rgb(0.7,0.7,0.7,0.4) , border = NA))
#'   for (i in 12:14) {
#'     matplot(workGrid[workGrid < (tau_est[i+1]- rho) & workGrid >= (tau_est[i]+ rho) ],
#'             mudata[workGrid < (tau_est[i+1] - rho) & workGrid >= (tau_est[i]+rho), ],
#'             type = "l", add = TRUE, lty = c(3, 3, 1), col = c(4,4,2), lwd = 2)
#'   }
#'   legend('topleft', legend = c( 'Estimated mean', 'Pointwise confidence interval'),
#'          cex = 1, lwd = 2, col = c(2,4), lty = c(1,3), bty = "n")
#'   points(x = res$mu_jumptime[12:13], y = rep(-5.3, length(res$mu_jumptime[12:13])),
#'          pch= rep("*", length(res$mu_jumptime[12:13])), cex = 2, col = 4,  xpd = TRUE)
#'   # lines(res$obsGrid, res$mu, col = 4)
#'   text(x = res$mu_jumptime[12:13], y = rep(-5.3, length(res$mu_jumptime[12:13])),
#'        labels = data_SS$Date[round(res$mu_jumptime[12:13]*nrow(data_SS))],
#'        xpd = TRUE, pos = 3, cex = 0.8, col = 4)
#'
#'
#' }
#'
#' plot.fmb(res)
#'
#'
#' ## individual, Xia and Qiu 2016
#' #### for each curve 1,5,6,16,35, 36
#' industry = colnames(data_SS)[-50]
#' par(mfrow=c(2,2))
#' for (i in c(5,16,35, 36)) {
#'   plot(data_SS[, i+1],
#'        xlab = "Dates",
#'        ylab = "Average Value Weighted Returns",
#'        main = industry[i+1],
#'        col = 'gray', #rgb(0.7,0.7,0.7,0.4),
#'        xaxt="n",
#'        ylim = c(-6,6))
#'   ind_x = seq(1, nrow(data_SS), length.out = 10)
#'   axis(1, at= ind_x, labels = data_SS[ind_x, 1])
#'   ## indivudual
#'   lines(y = res$indJump[[i]]$mu, x = Lt[[i]]*nrow(data_SS),
#'        col = "blue", lty = 1, lwd = 1.5)
#'   points(x = res$indJump[[i]]$mu_jumptime*nrow(data_SS),
#'          y = rep(-6.4, length(res$indJump[[i]]$mu_jumptime)),
#'          pch= rep("*", length(res$indJump[[i]]$mu_jumptime)),
#'          cex = 2, col = "blue",  xpd = TRUE)
#'   text(x = res$indJump[[i]]$mu_jumptime*nrow(data_SS),
#'        y = rep(-6.4, length(res$indJump[[i]]$mu_jumptime)),
#'        labels = data_SS$Date[res$indJump[[i]]$mu_jumptime*nrow(data_SS)],
#'        xpd = TRUE, pos = 3, cex = 0.5, col = "blue")
#'   ## our mean
#'   lines(obsGrid*nrow(data_SS), mu, col = "red", lty = 1, lwd = 1.5)
#'   abline( v = round(mu_jumptime *nrow(data_SS)), col = "red", lty = "dashed")
#'
#'
#' }
#' par(mfrow=c(1,1))
#' }
#'
#'
"avwreturn"
