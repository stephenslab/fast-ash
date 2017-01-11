#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix mmult1(NumericVector & pi, const NumericMatrix & matrix_lik){
  int n=matrix_lik.nrow(), k=matrix_lik.ncol();
  NumericMatrix m(n,k);
  for (int i=0;i<k;i++){
    m.column(i)=pi[i]*matrix_lik.column(i);
  }
  return(m);
}

// This code may be relevant
// http://stackoverflow.com/questions/18349053/fastest-way-for-multiplying-a-matrix-to-a-vector
// [[Rcpp::export]]
NumericMatrix mmult( NumericMatrix m , NumericVector v , bool byrow = true ){
  if(byrow){
    if( ! (m.ncol() == v.size()) ) stop("Non-conformable arrays") ;
  }
  if( ! byrow ){
    if( ! (m.nrow() == v.size()) ) stop("Non-conformable arrays") ;
  }
  
  NumericMatrix out(m) ;
  
  if( byrow ){
    for (int j = 0; j < m.ncol(); j++) {
      for (int i = 0; i < m.nrow(); i++) {
        out(i,j) = m(i,j) * v[j];
      }
    }
  }
  if( ! byrow ){
    for (int i = 0; i < m.nrow(); i++) {
      for (int j = 0; j < m.ncol(); j++) {
        out(i,j) = m(i,j) * v[i];
      }
    }
  }
  return out ;
}

// this multiplies columns of matrix m, between rows a and b, by vector v
// [[Rcpp::export]]
NumericMatrix my_mmult(const NumericMatrix & m ,const NumericVector & v , int a, int b){

  if( ! (m.ncol() == v.size()) ) stop("Non-conformable arrays") ;
    
  NumericMatrix out(b-a+1 , m.ncol()) ;
  
  for (int j = 0; j < m.ncol(); j++) {
    for (int i = 0; i < b-a+1; i++) {
      out(i,j) = m(i+a,j) * v[j];
    }
  }
  
  return out ;
}
