SmoothedMean <- function (y, t, weight, obsGrid, regGrid, optns, 
                          Wtype = c("OBS", "SUBJ", "MIX", "OPT")){
  
  
  ## weight type
  Wtype = match.arg(Wtype)
  mi = sapply(t, length)
  numOfCurves = length(y);
  
  userMu = optns$userMu;
  methodBwMu = optns$methodBwMu;
  npoly = 1
  nder = 0 
  userBwMu = optns$userBwMu; 
  kernel = optns$kernel
  
  # If the user provided a mean function use it
  if ( is.list(userMu) && (length(userMu$mu) == length(userMu$t))){
    
    buff <- .Machine$double.eps * max(abs(obsGrid)) * 10
    rangeUser <- range(optns$userMu$t)
    rangeObs <- range(obsGrid)
    if( rangeUser[1] > rangeObs[1] + buff || 
        rangeUser[2] < rangeObs[2] - buff   ) {
      stop('The range defined by the user provided mean does not cover the support of the data.')
    }
    
    mu = spline(userMu$t, userMu$mu, xout= obsGrid)$y;
    muDense = spline(obsGrid,mu, xout=regGrid)$y;
    bw_mu = NULL;
    
    # otherwise if the user provided a mean bandwidth use it to estimate the mean function (below)
  } else {
    if (userBwMu > 0){
      bw_mu = userBwMu;
      #otherwise estimate the mean bandwith via the method selected to estimnate the mean function (below)
    } else {
      if( any(methodBwMu == c('GCV','GMeanAndGCV') )){
        # get the bandwidth using GCV
        bw_mu =  unlist(GCVLwls1D1(yy = y, tt = t, kernel = kernel, npoly = npoly, nder = nder, dataType = optns$dataType) )[1]    
        if ( 0 == length(bw_mu)){ 
          stop('The data is too sparse to estimate a mean function. Get more data!\n')
        }
        # Uncomment to ensure MATLAB compatibility (AdjustBW1 is removed (3-Jun-16); check older versions.)
        # bw_mu = AdjustBW1(kernel=kernel,bopt=bw_mu,npoly=npoly,dataType=optns$dataType,nder=nder)
        # get the geometric mean between the minimum bandwidth and GCV bandwidth to estimnate the mean function (below)         
        if ( methodBwMu == 'GMeanAndGCV') {
          minbw = Minb( unlist(t),2)
          bw_mu = sqrt(minbw*bw_mu);
        } 
      } else {
        # get the bandwidth using CV to estimnate the mean function (below)
        bw_mu = CVLwls1D(y, t, kernel= kernel, npoly=npoly, nder=nder, dataType= optns$dataType, kFolds = optns$kFoldMuCov, 
                         useBW1SE = optns$useBW1SE); 
      }
    }
    # Get the mean function using the bandwith estimated above:
    xin = unlist(t);    
    yin = unlist(y)[order(xin)];
    xin = sort(xin);   
    if(is.null(weight)){
      win = weightFun(h = bw_mu, wi = weight, mi = mi, n = numOfCurves, Wtype = Wtype)$weight
    } else {
      win = weight #rep(1, length(xin));
    }
    if(is.null(regGrid)){
      regGrid = seq( max(min(obsGrid), bw_mu), min(max(obsGrid), 1- bw_mu), length.out = optns$nRegGrid);
    }
    mu = Lwls1D(bw_mu, kernel_type = kernel, npoly = npoly, nder = nder, xin = xin, yin= yin, xout = obsGrid, win = win)
    muWork = Lwls1D(bw_mu, kernel_type = kernel, npoly = npoly, nder = nder, xin = xin, yin= yin, xout = regGrid, win = win)
  }  
  
  
  ## Covariance function and sigma2
  scsObj = FPMD:::GetSmoothedCovarSurface(y = y, t = t, mu = mu, obsGrid = obsGrid, 
                                   regGrid = regGrid, optns = optns) 
  sigma2 <- scsObj[['sigma2']]
  
  # Get the results for the eigen-analysis
  eigObj = fdapace:::GetEigenAnalysisResults(smoothCov = scsObj$smoothCov, regGrid, optns, muWork = muWork)
  fittedCov = eigObj$fittedCov
  
  
  
  result <- list( 'mu' = mu, 
                  'muWork'= muWork, 
                  'obsGrid' = obsGrid, 
                  'workGrid' = regGrid,
                  'bw_mu' = bw_mu, 
                  'smoothedCov' = scsObj$smoothCov, 
                  'fittedCov' = fittedCov, 
                  'sigma2' = sigma2,
                  'wi' = win,
                  'mi' = mi);
  class(result) = 'resNCP'
  return(result)
}





CVLwls1D <- function(y, t, kernel, npoly, nder, dataType, kFolds = 5, useBW1SE = FALSE ){
  
  # If 'y' and 't' are vectors "cheat" and break them in a list of 10 elements
  if ( is.vector(y) && is.vector(t) && !is.list(t) && !is.list(y) ){
    if (length(t) < 21) {
      stop("You are trying to use a local linear weight smoother in a vector with less than 21 values.\n")
    }
    myPartition =   c(1:10, sample(10, length(t)-10, replace=TRUE));
    y = split(y, myPartition)
    t = split(t, myPartition)
    dataType = 'Sparse';
  } 
  
  # Make everything into vectors
  ncohort = length(t);
  tt  = unlist(t);
  yy  = unlist(y);
  ind = unlist(lapply( 1:ncohort, function(j) rep(j, times=length(t[[j]]))));
  yyn = yy[order(tt)];
  ind = ind[order(tt)];
  ttn = sort(tt);
  
  # Get minimum reasonable bandwidth
  bw = seq(0.01, 0.11, by = 0.02)
  nbw = length(bw)
  

  cv = matrix(Inf, ncol = length(bw), nrow = kFolds);
  #count = c();
  theFolds =  fdapace:::CreateFolds(unique(ind), k= kFolds)
  
  for (j in 1:(nbw-1)){
    # cv[j]=0;
    # count[j]=0;
    #for (i in 1:ncohort){
    for (i in 1:kFolds){
      
      xout= ttn[ ind %in% theFolds[[i]]];
      obs = yyn[ ind %in% theFolds[[i]]];
      xin = ttn[!ind %in% theFolds[[i]]];
      yin = yyn[!ind %in% theFolds[[i]]];
      
      win=rep(1,length(yin));
      #win[ind==i] = NA;        
      #if(dataType=='Dense') {
      #  yyn=(ave*ncohort-t[[i]])/(ncohort-1);
      #  ttn=t[[1]];
      #  win=pracma::ones(1,length(t[[1]]));    
      #  yyn = yyn[order(ttn)]
      #  ttn = sort(ttn)           
      #}  
      
      mu = tryCatch(
        Lwls1D(bw= bw[j], kernel_type = kernel, npoly=npoly, nder= nder, xin = xin, yin= yin, xout=xout, win = win), 
        error=function(err) {
          warning('Invalid bandwidth during CV. Try enlarging the window size.')
          return(Inf)
        })
      
      cv[i,j] = sum((obs-mu)^2)
      # print(cv)
      if(is.na(cv[i,j])){
        cv[i,j] = Inf;
      }
      #count[j] = count[j]+1;
    }
  }
  #cv = cv[(count/ncohort>0.90)];
  #bw = bw[(count/ncohort>0.90)];
  if(min(cv) == Inf){
    stop("All bandwidths resulted in infinite CV costs.")
  }
  if( useBW1SE ){
    # This will pick the bandwidth that is the max but it's average cost is at most
    # 1 standard error of the minimum cost /  I use means so it is more straighforward what the SE is.
    bopt = bw[max(which( 
      colMeans(cv) < min(colMeans(cv)) + apply(cv,2, sd)[which.min(colMeans(cv))]/sqrt(kFolds)))]
  } else {
    bopt = bw[which.min( colMeans(cv))];
  }
  
  return(bopt)
  
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





