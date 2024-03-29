% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Maturities.R
\docType{data}
\name{Maturities}
\alias{Maturities}
\title{Monthly U.S. Treasuries Data.}
\format{
\describe{
  \item{dmat}{Data matrix.}
  \item{Lt}{Observation time points.}
  \item{Ly}{Monthly maturity.}
  ...
}
}
\usage{
data(Maturities)
}
\description{
The monthly U.S. treasuries from January 1983 to September 2010.
In total, there are n = 15 interest rates at maturities of
3, 6, 9, 12, 18, 24, 30, 36, 48, 60, 72, 84, 96, 108, and 120 months.
}
\examples{

\dontrun{
data(Maturities)
dmat = Maturities$dmat
Ly = Maturities$Ly
Lt = Maturities$Lt


## global test
bw.seq = seq(0.01, 0.11, by = 0.02)
reject = rep(NA, length(bw.seq))
## reject H0: no jump points
for (i in 1:length(bw.seq)) {
  reject[i] = MDTest(bw = bw.seq[i], alpha = 0.05, Lt = Lt, Ly = Ly, Wtype = 'MIX')
  cat('reject:', reject[i], 'for h = ', bw.seq[i], '\n')
}

##
res = FPMD(Ly = Ly, Lt = Lt, wi = NULL, Wtype = "MIX", zeta = NULL,
               bw.seq = NULL, NbGrid = 101, kFolds = 5, refined = TRUE,
               individualCurveJump = TRUE, nRegGrid = 101,
               npoly = 1, nder = 0, alpha = 0.05, cutoff = max, M_max = 10)


## with change points
mu_jumptime <- res$mu_jumptime
mu_jumpsize <- res$mu_jumpsize
# mu_jumpsize <- res$mu_jumpsize_h_tau
h_tau <- res$h_tau
h_d <- res$h_d
zeta <- res$zeta
wi <- res$wi
mu <- res$mu
muWork <- res$muWork
obsGrid <- res$obsGrid
workGrid <- res$workGrid


library(zoo)
ym <- seq(as.yearmon("1983-01"), as.yearmon("2010-09"), 1/12)
ym = format(as.Date(ym), "\%Y-\%m")


par(mfrow=c(1,2))
matplot(dmat, type = "p", pch = 1,
        xlab = "Dates", ylab = "Maturities(Months)",
        #main = 'Monthly U.S. Treasuries',
        col = 'gray', xaxt="n")
ind_x = c(1, round(mu_jumptime*333), 333)
axis(1, at=ind_x, labels= ym[ind_x])
abline( v = round(mu_jumptime*333), lty = "dashed", col = "red")
lines(obsGrid*333, mu, col = "red", lwd = 1.5)


## combine all as iid, Xia and Qiu 2016
t = unlist(Lt);
y = unlist(Ly)[order(t)];
t = sort(t);
indJumpall = FPMD:::indMeanbreak(y = y, t = t, M_max = 15, NbGrid = 101,
                                 kernel = res$optns$kernel, npoly = 1, nder = 0)
lines(y = indJumpall$mu, x = (unique(t)*nrow(dmat)),
      col = "green", lty = 1, lwd = 1.5)
abline(v = round(indJumpall$mu_jumptime*nrow(dmat)),
       lty = 'dotted', col = "green")
points(x = indJumpall$mu_jumptime*nrow(dmat),
       y = rep(-.55, length(indJumpall$mu_jumptime)),
       pch= rep("*", length(indJumpall$mu_jumptime)),
       cex = 2, col = "green",  xpd = TRUE)
text(x = indJumpall$mu_jumptime*nrow(dmat),
     y = rep(-.55, length(indJumpall$mu_jumptime)),
     labels = ym[indJumpall$mu_jumptime*nrow(dmat)],
     xpd = TRUE, pos = 1, cex = 0.5, col = "green")
## hdbinseg package from H. Cho and P. Fryzlewicz (2014) JRSSB
library(hdbinseg)
dd = t(dmat)
ecp_CHO = dcbs.alg(dd, cp.type=1, phi= 1, temporal = TRUE, do.parallel=0)$ecp
# ecp_CHO = sbs.alg(dd, cp.type=1, temporal = TRUE, do.parallel=0)$ecp
## cho
points(x = ecp_CHO,
       y = rep(-0.55, length(ecp_CHO)),
       pch= rep("*", length(ecp_CHO)),
       cex = 2, col = "blue",  xpd = TRUE)
text(x = ecp_CHO,
     y = rep(-0.55, length(ecp_CHO)),
     labels = ym[ecp_CHO],
     xpd = TRUE, pos = 3, cex = 0.5, col = "blue")



## plot mean function and confidence band
## true and estimated mean curve are based on workGrid points
plot.fmb <- function(res){

  ## confidence band
  cbandMu <- FPMD:::pwCBFun(res)
  workGrid <- res$workGrid
  muWork <- res$muWork
  rho = res$rho
  tau_est = c(0, res$mu_jumptime, 1)
  mudata = t(rbind(cbandMu, muWork))


  plot(x = NULL, y = NULL, ylim = c(0, 11), xlim = range(unlist(Lt)),
          xlab = "Dates",
          ylab = "Maturities(Months)",
          main = '',
          col = rgb(0.7,0.7,0.7,0.4), xaxt="n")
  ind_x = c(min(Lt[[1]]), res$mu_jumptime, max(Lt[[1]]))
  axis(1, at=ind_x, labels= ym[round(ind_x*333)])
  ###

  bandGrid =  t(rbind(cbandMu, workGrid))
  lband = lapply(1:(length(tau_est) - 1), function(i)
    as.data.frame(bandGrid[workGrid < (tau_est[i + 1] - rho) &
                             workGrid >= (tau_est[i] + rho),]))
  lapply(lband, function(x) polygon(c(x$workGrid, rev(x$workGrid)),c(x$lwCI,rev(x$upCI)),
                                    col = rgb(0.7,0.7,0.7,0.4) , border = NA))
  for (i in 1:(length(tau_est)-1)) {

    matplot(workGrid[workGrid < (tau_est[i+1]- rho) & workGrid >= (tau_est[i]+ rho) ],
            mudata[workGrid < (tau_est[i+1] - rho) & workGrid >= (tau_est[i]+rho), ],
            type = "l", add = TRUE, lty = c(3, 3, 1), col = c(3,3,2), lwd = 2)

  }
  legend('top', legend = c( 'Estimated mean', 'Pointwise confidence interval'),
         lwd = 2, col = c(2,3), lty = c(1,3), bty = "n")
  points(x = res$mu_jumptime, y = rep(-.4, length(res$mu_jumptime)),
         pch= rep("*", length(res$mu_jumptime)), col = 4, cex  =2, xpd = TRUE)
  # # lines(res$obsGrid, res$mu, col = 2)
  text(x = res$mu_jumptime, y = rep(-.4, length(res$mu_jumptime)),
       labels= ym[round(res$mu_jumptime*nrow(dmat))],
       xpd = TRUE, pos = 3, cex = 0.8, col = 4)

}

plot.fmb(res)

## individual
####
par(mfrow=c(3,5))
ind = c(3, 6, 9, 12, 18, 24, 30, 36, 48, 60, 72, 84, 96, 108, 120)
for (i in 1:length(ind)) {

  plot(dmat[, i],
       xlab = "Dates",
       ylab = "Maturities",
       main = paste("maturities of", ind[i], "months"),
       col = 'gray', #rgb(0.7,0.7,0.7,0.4),
       xaxt="n")
  ind_x = c(1, round(mu_jumptime*333), 333)
  axis(1, at=ind_x, labels= ym[ind_x])
  lines(y = res$indJump[[i]]$mu, x = Lt[[i]]*333,
        col = "blue", lty = 1, lwd = 1.5 )
  points(x = round(res$indJump[[i]]$mu_jumptime*333),
         y = rep(min(dmat[, i])-0.45, length(res$indJump[[i]]$mu_jumptime)),
         pch= rep("*", length(res$indJump[[i]]$mu_jumptime)),
         cex = 2, col = "blue", xpd = TRUE)
  text(x = round(res$indJump[[i]]$mu_jumptime*333),
       y = rep(min(dmat[, i])-0.45, length(res$mu_jumptime)),
       labels= ym[res$indJump[[i]]$mu_jumptime*nrow(dmat)],
       xpd = TRUE, pos = 3, cex = 0.5, col = "blue")
  ## our
  lines(obsGrid*333, mu, col = "red", lty = 1, lwd = 1.5)
  abline( v = round(mu_jumptime *nrow(dmat)), col = "red", lty = "dashed")

}
par(mfrow=c(1,1))



}


}
\keyword{datasets}
