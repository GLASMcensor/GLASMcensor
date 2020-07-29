#include <Rcpp.h>
#include <RcppEigen.h>
#include <omp.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

NumericVector ccumprod(NumericVector x);
//' The proposed marginal screening method based on the BKR independence test
//' @name cMBKRall
//' @param Ti is the vector of the ordered censored Time (ascent)
//' @param censor is the censor indicator for the Ti. 1 means no censor and 0 means censoring.
//' @param X is the complete covariates for the samples.
//' @param weight_C the inverse weight in the survival analysis. it is the output of the function KMestimate_CTime.
//' @export
//' @return the cBKRmeasure

// [[Rcpp::export]]
NumericVector cMBKRall(NumericVector Ti,NumericVector censor,NumericMatrix X,NumericVector weight_C) {
  // omp_set_num_threads(4);
  int N=Ti.length();
  int P=X.ncol();
  NumericVector uti(P);
  NumericVector xk(N);
  NumericVector Fx(N);
  NumericVector Fy(N);  
  NumericVector Iyi(N);
  NumericVector Ixi(N);
  NumericVector weight_fxy(N);
  NumericVector ccoef_fxy(N);
  NumericMatrix Fxy(N,N);
  NumericMatrix w(N,N);
  // int i;
  // int j;
  // int ii;
  int k;
  // double xi;
  // double tj;
  // double result;
  double Ix;
  double It;
  //double w = 0;
  // double weach=0;
  // double numerator=0;
  // double dominator=0;
  
  weight_fxy[0]=1;
  //#pragma omp parallel for 
  for(int i=1;i<N;i++)
  {
    weight_fxy[i]=pow(((double(N-i))/(double(N-i+1))),censor[i]);
  }

  ccoef_fxy=ccumprod(weight_fxy);


  for(k =0;k<P;k++){
    xk = X(_,k);
    xk=(xk-mean(xk))/sd(xk);

 #pragma omp parallel for 
  for(int i=0;i<N;i++)
  {
    double xi=xk[i];
    for(int j=0;j<N;j++)
    {
    double tj=Ti[j];
    double result=0;
        for(int ii=0;ii<N;ii++)
        {
          if(xk[ii]<=xi){Ix=1.0;}else{Ix=0.0;}
          if(Ti[ii]<=tj){It=1.0;}else{It=0.0;}
          result=result + censor[ii]*ccoef_fxy[ii]*Ix*It/(double(N-ii));
        }
        Fxy(i,j)=result;
    }
  }
  //#pragma omp parallel for
  for(int i =0 ;i<N;i++){
    Iyi   = ifelse( Ti<=Ti[i], 1.0, 0.0);
    Ixi   = ifelse( xk<=xk[i], 1.0, 0.0);
    Fx[i]=sum(Ixi)/N;
    Fy[i]=mean(Iyi* censor * weight_C);

  }
  
  
  #pragma omp parallel for
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      double weach=0;
    double  numerator = pow(Fxy(i,j)-Fx[i]*Fy[j],2);
    double dominator = Fx[i]*(1-Fx[i])*Fy[j]*(1-Fy[j]);
      if(dominator<0.000001){
        weach=0;
      }else{
      weach=numerator/dominator;
      }
      w(i,j) =  weach;
    }
  }
  uti[k]=0;
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      uti[k]=uti[k]+w(i,j);
      }
    }
  uti[k]=uti[k]/pow(N,2);
}
  return(uti);
}



NumericVector ccumprod(NumericVector x) {
  return(cumprod(x));
}
