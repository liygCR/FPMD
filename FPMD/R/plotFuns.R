#' Plot the estimated mean function and the corresponding pointwise confidence interval
#'
#' @param res an FPMD class object returned by \code{FPMD()}.
#'
#' @export

## plot mean function and confidence band
## true and estimated mean curve are based on workGrid points
PCBplot <- function(res){

  ## confidence band
  cbandMu <- pwCBFun(res)
  workGrid <- res$workGrid
  muWork <- res$muWork
  h_tau <- res$h_tau
  rho = res$rho

  # tau_est = c(h_tau, res$mu_jumptime_refine, 1-h_tau)
  tau_est = c(0, res$mu_jumptime, 1)


  mudata = t(rbind(cbandMu, muWork))
  plot(NA, ylab = expression(paste(mu,"(t)")), ylim = range(-1.5,3), xlim = c(0, 1), type = "n")
  ###
  bandGrid =  t(rbind(cbandMu, workGrid))
  lband = lapply(1:(length(tau_est) - 1), function(i)
    as.data.frame(bandGrid[workGrid < (tau_est[i + 1] - rho) & workGrid >= (tau_est[i] + rho),]))
  lapply(lband, function(x) polygon(c(x$workGrid, rev(x$workGrid)),c(x$lwCI,rev(x$upCI)),
                                    col = rgb(0.7,0.7,0.7,0.4) , border = NA))
  for (i in 1:(length(tau_est)-1)) {

    matplot(workGrid[workGrid < (tau_est[i+1]- rho) & workGrid >= (tau_est[i]+ rho) ],
            mudata[workGrid < (tau_est[i+1] - rho) & workGrid >= (tau_est[i]+rho), ],
            type = "l", add = TRUE, lty = c(3, 3, 1), col = c(3,3,2), lwd = 2)

  }
  legend('topleft', legend = c('Estimated mean', 'Pointwise confidence interval'),
         cex = .8, lwd = 2, col = c(2,3), lty = c(1,3), bty = "n")
}


###### pointwise confidence band #######
pwCBFun <- function(par){
  obsGrid <- par$obsGrid
  workGrid <- par$workGrid
  bw <- par$h_tau
  muWork <- par$muWork
  covR.hat <- par$smoothedCov
  if (det(covR.hat) < 1e-16) {
    covR.hat <- par$fittedCov
  }
  sigma2 <- par$sigma2
  mi <- par$mi
  wi <- par$wi
  ## density f_T is uniform
  df <- density(runif(10000), kernel = "rectangular")
  f.T <- approx(df$x,df$y,xout=workGrid)$y
  K_Epa = function(u) 0.75 * (1 - u^2) * (abs(u) <= 1); # epanechnikov kernel
  # f.T <- sapply(workGrid, function(t){ mean(K_Epa( (obsGrid - t)/bw)/bw) })
  K_Epa2 = function(u) (0.75 * (1 - u^2) * (abs(u) <= 1))^2
  normK <- integrate(K_Epa2, -1, 1)$value
  ##
  pwVar <- (sum(mi*wi^2)/bw)*normK * (diag(covR.hat) + sigma2)/f.T +
    sum(mi*(mi-1)*wi^2)*diag(covR.hat)
  ##
  alpha = 0.05
  upCI <- t(muWork + stats::qnorm(1-alpha / 2) * sqrt(pwVar))
  lwCI <- t(muWork - stats::qnorm(1-alpha / 2) * sqrt(pwVar))
  cbandMu <- rbind(upCI, lwCI) ## on workGrid
  rownames(cbandMu) = c('upCI', 'lwCI')

  return(cbandMu)
}

## confidence band cover
coverFun <- function(res, mu_fun) {

  cbandMu <- pwCBFun(res)
  workGrid <- res$workGrid
  rho = res$rho
  h_tau <- res$h_tau
  tau_est = c(h_tau, res$mu_jumptime, 1-h_tau)

  if (res$mu_jumptime[1] > 0) {
    # tau_true = c(h_tau, tau, 1-h_tau)
    # LworkGrid =  lapply(1:(length(tau_est) - 1), function(i)
    #   workGrid[workGrid < min(tau_est[i + 1] - rho, tau_true[i+1]) &
    #              workGrid >= max(tau_est[i] + rho, tau_true[i])])
    LworkGrid =  lapply(1:(length(tau_est) - 1), function(i)
      workGrid[workGrid < tau_est[i + 1] - rho &
                 workGrid >= tau_est[i] + rho])
    LworkGrid = LworkGrid[!sapply(LworkGrid, isempty)]
    Lmutrue = lapply(LworkGrid, mu_fun)
    # LcbandMu = lapply(1:(length(tau_est) - 1), function(i)
    #   cbandMu[,workGrid < min(tau_est[i + 1] - rho, tau_true[i+1]) &
    #             workGrid >= max(tau_est[i] + rho, tau_true[i]), drop = FALSE])
    LcbandMu = lapply(1:(length(tau_est) - 1), function(i)
      cbandMu[,workGrid < tau_est[i + 1] - rho &
                workGrid >= tau_est[i] + rho, drop = FALSE])
    LcbandMu = LcbandMu[!sapply(LcbandMu, isempty)]

    ## simutanious confidence band
    # ind = all(mapply(function(x,y)all(x >= y[2,]) & all(x <= y[1,]), x = Lmutrue, y = LcbandMu))
    ## pointwise confidence band
    pcb = mean(do.call(c,mapply(function(x,y) x >= y[2,] & x <= y[1,],
                                x = Lmutrue, y = LcbandMu, SIMPLIFY = FALSE)))
  } else {

    mutrue = sapply(workGrid, mu_fun)
    ## simutanious confidence band
    # scb = all(mutrue >= cbandMu[2,]) & all(mutrue <= cbandMu[1,])
    ## pointwise confidence band
    pcb = mean(mutrue >= cbandMu[2,] & mutrue <= cbandMu[1,])
  }



  return(pcb)
}




