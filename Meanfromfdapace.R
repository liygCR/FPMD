SmoothedMean <- function (y, t, weight, obsGrid, regGrid, optns){
  
  # For the case of binned data one may use a weighted mean response for each time-point.
  # This is not currently implemented. \hat{y}_i = \sum_i w_i y_i where w_i are the
  # same points for common t_is so we have: \hat{y}_i = n_t w_i \bar{y}
  
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
    win = weight #rep(1, length(xin));
    mu = Lwls1D(bw_mu, kernel_type = kernel, npoly = npoly, nder = nder, xin = xin, yin= yin, xout = obsGrid, win = win)
    muDense = Lwls1D(bw_mu, kernel_type = kernel, npoly = npoly, nder = nder, xin = xin, yin= yin, xout = regGrid, win = win)
  }  
  
  result <- list( 'mu' = mu, 'muDense'= muDense, 'bw_mu' = bw_mu);
  class(result) <- "SMC"
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





