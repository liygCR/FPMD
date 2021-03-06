#include <RcppEigen.h>
#include <map>          // to map kernels to integers for the switch
#include <string>       // to read in the kernel name
#include <vector>       // to use vectors
#include <algorithm>    // to get the intersect and sort

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]

Eigen::MatrixXd CPPlwls1d_LR( const double & bw, const std::string kernel_type, const Eigen::Map<Eigen::VectorXd> & win, 
                              const Eigen::Map<Eigen::VectorXd> & xin, const Eigen::Map<Eigen::VectorXd> & yin, 
                              const Eigen::Map<Eigen::VectorXd> & xout, const unsigned int & npoly = 1, 
                              const unsigned int & nder = 0){
  
  
  // Convenient constants
  const double invSqrt2pi=  1./(sqrt(2.*M_PI));
  const double factorials[] = {1,1,2,6,24,120,720,5040,40320,362880,3628800};
  
  const unsigned int nXGrid = xin.size();
  const unsigned int nUnknownPoints = xout.size();
  Eigen::MatrixXd result(nUnknownPoints, 2);
  
  // ========================
  // The checks start here:
  
  if(nXGrid == 0) {
    Rcpp::stop("The input X-grid has length zero.");
  }
  
  // Check that we have equal number of readings
  if( nXGrid != yin.size()){
    Rcpp::stop("The input Y-grid does not have the same number of points as input X-grid.");
  }
  
  if( nXGrid != win.size()){
    Rcpp::stop("The input weight vector does not have the same number of points as input X-grid.");
  }
  
  // Check that bandwidth is greater than zero
  if( bw <= 0.){
    Rcpp::stop("The bandwidth supplied for 1-D smoothing is not positive.");
  }
  
  // Check that the degreee of polynomial used and the order of the derivative are reasonable
  if (npoly < nder){
    Rcpp::stop("The degree of polynomial supplied for 1-D smoothing is less than the order of derivative");
  }
  
  // Map the kernel name so we can use switches  
  std::map<std::string,int> possibleKernels;
  possibleKernels["epan"]    = 1;   possibleKernels["rect"]    = 2;
  possibleKernels["gauss"]   = 3;   possibleKernels["gausvar"] = 4;
  possibleKernels["quar"]    = 5;
  
  // If the kernel_type key exists set KernelName appropriately
  int KernelName = 0;
  if ( possibleKernels.count( kernel_type ) != 0){
    KernelName = possibleKernels.find( kernel_type )->second; //Set kernel choice
  } else {
    // otherwise use "epan"as the kernel_type 
    Rcpp::warning("Kernel_type argument was not set correctly; Epanechnikov kernel used.");
    KernelName = possibleKernels.find( "epan" )->second;;
  }
  
  // Check that we do not have zero weights // Should do a try-catch here
  if ( !(win.all()) ){  // 
    Rcpp::warning("Cases with zero-valued windows maybe not be too safe.");
  }
  // Check if the first 5 elements are sorted // Very rough check // Will use issorted in the future
  if ( (xin[0]> xin[1]) ||  (xin[1]> xin[2]) ||  (xin[2]> xin[3]) ||  (xin[3]> xin[4]) ||  (xin[4]> xin[5]) ){
    Rcpp::stop("The X-grid used is not sorted. (or you have less than 6 points)");
  }
  
  // The checks end here.
  // ===================
  
  for (unsigned int i = 0; i != nUnknownPoints; ++i){
    //locating local window (LOL) (bad joke)
    std::vector <unsigned int> indx;
    const double* lower ;
    const double* upper ;
    
    //if the kernel is not Gaussian
    if ( KernelName != 3 && KernelName != 4) {
      //construct listX as vectors / size is unknown originally
      // for (unsigned int y = 0; y != nXGrid; ++y){ if ( std::fabs( xout(i) - xin(y) ) <= bw ) { indx.push_back(y); }  }
      // Get iterator pointing to the first element which is not less than xou(u)
      lower = std::lower_bound(xin.data(), xin.data() + nXGrid, xout(i) - bw);
      upper = std::lower_bound(xin.data(), xin.data() + nXGrid, xout(i) + bw);
      //  const unsigned int firstElement = lower - &xin[0];
      //  for (unsigned int xx1 =0; xx1 != upper-lower; ++xx1){
      //   indx.push_back(  firstElement+ xx1 );
      //  }
    } else {
      lower = xin.data();
      upper = xin.data() + nXGrid;
    }
    
    const unsigned int firstElement = lower - &xin[0];
    for (unsigned int xx1 =0; xx1 != upper-lower; ++xx1){
      indx.push_back(  firstElement+ xx1 );
    }
    
    // for(unsigned int r4=0; r4 != indx.size(); r4++){ Rcpp::Rcout << indx.at(r4) << ", "; } Rcpp::Rcout << "\n";
    
    unsigned int indxSize = indx.size();
    Eigen::VectorXd temp_l(indxSize);
    Eigen::VectorXd temp_r(indxSize);
    Eigen::VectorXd lw(indxSize);
    Eigen::VectorXd ly(indxSize);
    Eigen::VectorXd lx(indxSize);
    
    for (unsigned int y = 0; y != indxSize; ++y){
      lx(y) = xin(indx[y]);
      lw(y) = win(indx[y]);
      ly(y) = yin(indx[y]);
    }
    
    
    Eigen::VectorXd llx = (lx.array() -xout(i)) * (1./bw) ;
    
    //define the kernel used 
    switch (KernelName){
    case 1: // Epan
      for (unsigned int y = 0; y != indxSize; ++y){
        if(-1 <= llx(y) && llx(y) < 0) temp_l(y) = (1- pow(llx(y), 2))*1.5*(lw(y));
        else temp_l(y) = 0;
      }
      for (unsigned int y = 0; y != indxSize; ++y){
        if(0 <= llx(y) && llx(y) <= 1) temp_r(y) = (1- pow(llx(y), 2))*1.5*(lw(y));
        else temp_r(y) = 0;
      }
      break;
    case 2 : // Rect
      for (unsigned int y = 0; y != indxSize; ++y){
        if(-1 <= llx(y) && llx(y) < 0) temp_l(y) = lw(y);
        else temp_l(y) = 0;
      }
      for (unsigned int y = 0; y != indxSize; ++y){
        if(0 <= llx(y) && llx(y) <= 1) temp_r(y) = lw(y);
        else temp_r(y) = 0;
      }
      break;
    case 3 : // Gauss
      for (unsigned int y = 0; y != indxSize; ++y){
        if(-1 <= llx(y) && llx(y) < 0) temp_l(y) = exp(-.5*pow(llx(y), 2)) * invSqrt2pi *lw(y);
        else temp_l(y) = 0;
      }
      for (unsigned int y = 0; y != indxSize; ++y){
        if(0 <= llx(y) && llx(y) <= 1) temp_r(y) = exp(-.5*pow(llx(y), 2)) * invSqrt2pi *lw(y);
        else temp_r(y) = 0;
      }
      break;
    case 4 : // GausVar
      for (unsigned int y = 0; y != indxSize; ++y){
        if(-1 <= llx(y) && llx(y) < 0) temp_l(y) = exp(-.5*pow(llx(y), 2)) * invSqrt2pi *lw(y) *
          (1.25 - 0.25 * pow(llx(y), 2));
        else temp_l(y) = 0;
      }
      for (unsigned int y = 0; y != indxSize; ++y){
        if(-1 <= llx(y) && llx(y) < 0) temp_r(y) = exp(-.5*pow(llx(y), 2)) * invSqrt2pi *lw(y) *
          (1.25 - 0.25 * pow(llx(y), 2));
        else temp_r(y) = 0;
      }
      break;
    case 5 :  // Quar
      for (unsigned int y = 0; y != indxSize; ++y){
        if(-1 <= llx(y) && llx(y) < 0) temp_l(y) = pow((1.- pow(llx(y), 2)), 2) * (15./16.);
        else temp_l(y) = 0;
      }
      for (unsigned int y = 0; y != indxSize; ++y){
        if(-1 <= llx(y) && llx(y) < 0) temp_r(y) = pow((1.- pow(llx(y), 2)), 2) * (15./16.);
        else temp_r(y) = 0;
      }
      break;
    }
    
    if(nder >= indxSize){
      Rcpp::warning("Cannot estimate derivatives of order p with less than p+1 points.");
      result(i) = std::numeric_limits<double>::quiet_NaN();
    } else {
      // make the design matrix
      Eigen::MatrixXd X(indxSize, npoly+1);
      X.setOnes();
      for (unsigned int y = 1; y <= npoly; ++y){
        X.col(y) = (xout(i) - lx.array()).array().pow(y);
      }
      
      Eigen::LDLT<Eigen::MatrixXd> ldlt_XTWlX(X.transpose() * temp_l.asDiagonal() *X);
      Eigen::LDLT<Eigen::MatrixXd> ldlt_XTWrX(X.transpose() * temp_r.asDiagonal() *X);
      // The solver should stop if the value is NaN. See the HOLE example in gcvlwls2dV2.
      Eigen::VectorXd betal = ldlt_XTWlX.solve(X.transpose() * temp_l.asDiagonal() * ly);
      Eigen::VectorXd betar = ldlt_XTWrX.solve(X.transpose() * temp_r.asDiagonal() * ly);
      
      //  Rcpp::Rcout << "lx: " << lx.transpose() << std::endl;
      //  Rcpp::Rcout << "ly: " << ly.transpose() << std::endl;
      //  Rcpp::Rcout << "temp: " << temp.transpose() << std::endl;
      //  Rcpp::Rcout << "llx: " << llx.transpose() << std::endl;
      //  Rcpp::Rcout << "xin: " << xin.transpose() << std::endl;
      //  Rcpp::Rcout << "yin: " << yin.transpose() << std::endl;
      //  Rcpp::Rcout << "xout: " << xout.transpose() << std::endl;
      //  Rcpp::Rcout << "X: " << X.transpose() << std::endl;
      //  Rcpp::Rcout << "beta: " << beta.transpose() << std::endl;
      //  Rcpp::Rcout << "factorials[nder]: " << factorials[nder]  << std::endl;
      //  Rcpp::Rcout << "pow (-1.0, nder): " << pow (-1.0, nder) << std::endl;
      
      result(i,0) = betal(nder+0) * factorials[nder] *  std::pow(-1.0, int(nder));
      result(i,1) = betar(nder+0) * factorials[nder] *  std::pow(-1.0, int(nder));
    }
  }
  return result;
}
