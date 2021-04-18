library(devtools)
devtools::install_github("liygCR/FPMD/FPMD")
library(FPMD)
library(fdapace)
library(extraDistr)
library(pracma)


##################### 200 replication simulation #################
# setting 1
## 3 change points
tau = c(.25,.5,.75) # jump locations
disc = c(.5,-.4,.4) # jump sizes
smoothmean=function(x) x^2+sin(2*pi*x) + cos(2*pi*x)
jump = function(x) sum((x>=tau)*disc) # jump function
mu_fun= function(t){smoothmean(t)+ sapply(t, jump)}

# setting 2
## 5 change point
tau= round(seq(0,1, length.out = 7)[2:6], 2) # jump locations
disc = c(.45, -.5, .45, -.5, .45) # jump sizes
smoothmean=function(x) .5*sin(3*pi*x) # smooth function
jump = function(x) sum((x>=tau)*disc) # jump function
mu_fun <- function(t){smoothmean(t)+ sapply(t, jump)}



###
simFun <- function(n, Q, mu_fun){

  ## generate gaussian process
  # "cos", "sin", "fourier", "legendre01", "poly"
  data = MakeFPMBData(n = n, Q = Q, proc= "gaussian", muFun = mu_fun, K = 3,
                      lambda = 1/(seq_len(3)+1)^2, sigma = 0.2, basisType = "fourier")
  # data = MakeFPMBData(n = n, Q = Q, proc= "t", muFun = mu_fun, K = 3,
  #                     lambda = 1/(seq_len(3)+1)^2, sigma = 0.1, basisType = "fourier")
  Lt = data$Lt
  Ly = data$Ly
  bw.seq = seq(0.01, 0.13, by = 0.02)
  resCP = FPMD(Ly = Ly, Lt = Lt, wi = NULL, Wtype = "SUBJ", zeta = NULL,
               bw.seq = bw.seq, NbGrid = 101, kFolds = 5, refined = TRUE,
               individualCurveJump = FALSE, nRegGrid = 101,
               npoly = 1, nder = 0, alpha = 0.05, cutoff = max)

  ## fdapace
  optns = resCP$optns
  # optns$kFoldMuCov = 5
  optns$useBW1SE = FALSE;
  obsGrid = sort(unique( c(unlist(Lt))));
  resNCP = SmoothedMean(y = Ly, t = Lt, weight = NULL, obsGrid = obsGrid,
                        regGrid = NULL, optns = optns, Wtype = "MIX")

  # resNCP = NULL

  return(list(resCP = resCP, resNCP= resNCP))

}

### computation time
system.time(res <- simFun(n = 400, Q = 20, mu_fun = mu_fun))


## parallel computing for simulations
## only for Linux\Unix\Mac server
library(doParallel)
library(foreach)
library(doMC)
library(doRNG)

registerDoMC(10)
n_run = 200 # simulation run
# set.seed(11)
ptm <- proc.time()
Results <- foreach(1:n_run) %dorng% {
  simFun(n = 400, Q = 20, mu_fun = mu_fun)
}
proc.time() - ptm



# compare with no change points -------------------------------------------
### compare with no change points
resCP = lapply(Results, function(res) res$resCP)
## Calculate the proportion of intervals that cover.
mean(sapply(resCP, FPMD:::coverFun, mu_fun))*100
##
resNCP = lapply(Results, function(res) res$resNCP)

##
plot(resCP[[1]]$obsGrid, mu_fun(resCP[[1]]$obsGrid), type="l")
lines(resCP[[1]]$obsGrid, resNCP[[1]]$mu, col=2)
lines(resCP[[1]]$obsGrid, resCP[[1]]$mu, col=4)
#
ISE_fun <- function(resNCP, resCP) {

  ISE.ncp <- trapz(resCP$obsGrid, (resNCP$mu - mu_fun(resCP$obsGrid))^2)
  ISE.cp <- trapz(resCP$obsGrid, (resCP$mu - mu_fun(resCP$obsGrid))^2)

  ISE = cbind(ISE.ncp, ISE.cp)

  return(ISE)
}

ISE = mapply(ISE_fun, resNCP, resCP)
apply(ISE, 1, mean)
apply(ISE, 1, sd)


##############
# resCP = lapply(Results, function(x)x$resCP)
resCP <- resCP[!sapply(resCP, function(x)is.null(x$h_d) )]


# bandwidth selection -----------------------------------------------------
## bandwidth
h_tau <- sapply(resCP, function(a){a$h_tau})
mean(h_tau)
table(h_tau)
h_d <-  sapply(resCP, function(a){a$h_d})
mean(h_d)
zeta <-  sapply(resCP, function(a){a$zeta})
mean(zeta)

## without cp from fdapace
bw <- sapply(resNCP, function(a){a$bw_mu})
mean(bw)
table(bw)


# change point locations ---------------------------------------------------
##
Locofcp <- sapply(resCP, function(a){a$mu_jumptime}, simplify = FALSE)
discp <- sapply(resCP, function(a){a$mu_jumpsize_h_tau}, simplify = FALSE)
## if 3 or 5 changes
ncp = length(tau)
mean(sapply(Locofcp, length)) - ncp
median(sapply(Locofcp, length)) - ncp
table(sapply(Locofcp, length))
## which detect true change points
Locofcp2 <- do.call(rbind, Locofcp[sapply(Locofcp, function(x){length(x)==ncp})])
nrow(Locofcp2)/length(resCP)

## jump location consistency
## hausdorff distance
hds <- function (P, Q) {
  stopifnot(is.numeric(P), is.numeric(Q))
  if (is.vector(P))
    P <- matrix(P, ncol = 1)
  if (is.vector(Q))
    Q <- matrix(Q, ncol = 1)
  if (ncol(P) != ncol(Q))
    stop("'P' and 'Q' must have the same number of columns.")
  D <- pracma::distmat(P, Q)
  dhd_PQ <- max(apply(D, 1, min))
  dhd_QP <- max(apply(D, 2, min))
  return(cbind(dhd_PQ, dhd_QP))
}

##
loc_dis <- sapply(Locofcp, hds, Q = tau)
apply(loc_dis,1,mean)
apply(loc_dis,1,sd)

## jump size consistency
id = sapply(Locofcp, function(x){length(x)==ncp})
jump_size_est <- sapply(resCP[id], function(a){a$mu_jumpsize}, simplify = FALSE)
max_jump_size_est <- sapply(jump_size_est, function(x){max(abs(x-disc))})
mean(max_jump_size_est)
sd(max_jump_size_est)


# mean and covariance estimate after jump points -----------------------------------------
## plot the mean
# which.median <- function(x) which.min(abs(x - median(x)))
median.id <-  which.min(abs(ISE[2,] - median(ISE[2,])))
# median.id <-  which.min(abs(performance[2,] - median(performance[2,])))
res <- resCP[[median.id]]
PCBplot(res)
## model error
err <- sapply(resCP, function(a){a$sigma2}, simplify = FALSE)
boxplot(unlist(err[sapply(Locofcp, function(x){length(x)==ncp})]))
