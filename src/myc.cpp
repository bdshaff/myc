#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

//' myc_proj
//'
//' Compute vector projection.
//'
//' @param u NumericVector
//' @param a NumericVector
// [[Rcpp::export]]
NumericVector myc_proj(NumericVector u, NumericVector a) {
  
  int length = u.size();
  double s1 = 0;
  double s2 = 0;
  NumericVector proj(length);
  
  for(int i = 0; i < length; i++){
    s1 = s1 + u[i]*u[i];
    s2 = s2 + u[i]*a[i];
  }
  
  double sc = s2/s1;
  
  for(int i = 0; i < length; i++){
    proj[i] = sc*u[i];
  }
  
  return proj;
}

//' myc_matmult
//'
//' Matrix Multiplication
//'
//' @param A NumericMatrix
//' @param B NumericMatrix
// [[Rcpp::export]]
NumericMatrix myc_matmult(NumericMatrix A, NumericMatrix B) {
        
        int n = A.nrow();
        int m = B.ncol();
        int p = B.nrow();
        
        NumericMatrix res(n, m);
        for(int i = 0; i < n; i++){
                for(int j = 0; j < m; j++){
                        double v = 0;
                        int k = 0;
                        while(k < p){
                                double a = A(i,k);
                                double b = B(k,j);
                                v = v + a*b;
                                k = k + 1;
                        }
                        res(i,j) = v;
                }
        }
        
        return res;
}

//' myc_diag
//'
//' Compute vector projection.
//'
//' @param A NumericMatrix
// [[Rcpp::export]]
NumericVector myc_diag(NumericMatrix A){
  int n = A.nrow();
  NumericVector vec(n);
  
  for(int i = 0; i < n; i++){
    vec[i] = A(i,i);
  }
  
  return vec;
}

//' myc_qr
//'
//' Compute QR Factorization
//'
//' @param A NumericMatrix
// [[Rcpp::export]]
List myc_qr(NumericMatrix A) {
  
  int a = A.rows();
  int b = A.cols();
  
  NumericMatrix U(a,b);
  NumericMatrix Q(a,b);
  NumericMatrix::Column Ucol1 = U(_ , 0);
  NumericMatrix::Column Acol1 = A(_ , 0);
  
  Ucol1 = Acol1;
  
  for(int i = 1; i < b; i++){
    NumericMatrix::Column Ucol = U(_ , i);
    NumericMatrix::Column Acol = A(_ , i);
    NumericVector subt(a);
    int j = 0;
    while(j < i){
      NumericVector uj = U(_ , j);
      NumericVector ai = A(_ , i);
      subt = subt + myc_proj(uj, ai);
      j++;
    }
    Ucol = Acol - subt;
  }
  
  for(int i = 0; i < b; i++){
    NumericMatrix::Column ui = U(_ , i);
    NumericMatrix::Column qi = Q(_ , i);
    
    double sum2_ui = 0;
    for(int j = 0; j < a; j++){
      sum2_ui = sum2_ui + ui[j]*ui[j];
    }
    
    qi = ui/sqrt(sum2_ui);
    
  }
  
  NumericMatrix R = myc_matmult(transpose(Q), A);
  
  List L = List::create(Named("Q") = Q,
                        Named("R") = R);
  
  return L;
  
}

//' myc_eigen
//'
//' Compute Eigen
//'
//' @param A NumericMatrix
// [[Rcpp::export]]
List myc_eigen(NumericMatrix A, double margin = 1e-20){
  List QR = myc_qr(A);
  NumericMatrix Q = QR["Q"];
  NumericMatrix U = QR["Q"];
  
  NumericMatrix V = myc_matmult(A,Q);
  NumericMatrix E = myc_matmult(transpose(Q),V);
  
  NumericVector d1 = myc_diag(A);
  NumericVector d2 = myc_diag(E);
  while(sum(pow(d1 - d2, 2.0)) > 1e-20){
    d1 = d2;
    
    List QR = myc_qr(E);
    NumericMatrix Q = QR["Q"];
    
    V = myc_matmult(E,Q);
    E = myc_matmult(transpose(Q),V);

    U = myc_matmult(U,Q);
    
    d2 = myc_diag(E);
  }
  
  List L = List::create(Named("values") = d2,
                        Named("vectors") = U);
  
  return L;

}

//' myc_svd
//'
//' Compute SVD
//'
//' @param A NumericMatrix
// [[Rcpp::export]]
List myc_svd(NumericMatrix A, double margin = 1e-20){
  int ncol = A.cols();
  int nrow = A.rows();
  if(ncol >= nrow){
    NumericMatrix X = myc_matmult(A, transpose(A));
    List E = myc_eigen(X);
    NumericVector d = E["values"];
    d = sqrt(d);
    NumericMatrix U = E["vectors"];
    NumericMatrix S(d.length(), d.length());
    for(int i = 0; i < d.length(); i++){
      for(int j = 0; j < d.length(); j++)
        if(i == j){
          S(i,j) = 1.0/d(i);
        }
    }
    
    NumericMatrix O = myc_matmult(transpose(U),A);
    NumericMatrix V = myc_matmult(S, O);
    V = transpose(V);
    
    List L = List::create(Named("d") = d,
                          Named("U") = U,
                          Named("V") = V,
                          Named("S") = S);
    
    return L;
    
  }else{
    NumericMatrix X = myc_matmult(transpose(A), A);
    List E = myc_eigen(X);
    NumericVector d = E["values"];
    d = sqrt(d);
    NumericMatrix V = E["vectors"];
    NumericMatrix S(d.length(), d.length());
    for(int i = 0; i < d.length(); i++){
      for(int j = 0; j < d.length(); j++)
        if(i == j){
          S(i,j) = 1.0/d(i);
        }
    }
    
    NumericMatrix O = myc_matmult(A,V);
    NumericMatrix U = myc_matmult(O,S);
    
    List L = List::create(Named("d") = d,
                          Named("U") = U,
                          Named("V") = V,
                          Named("S") = S);
    
    return L;
  }
  
  
}

//' myc_dist
//'
//' Compute Distance
//'
//' @param A NumericMatrix
// [[Rcpp::export]]
List myc_dist(NumericMatrix A){
  int nn = A.rows();
  NumericMatrix DD(nn,nn);
  for(int i = 0; i < nn; i++){
    for(int j = nn-1; j >= 0; j--){
      if(i < j){
        DD(i,j) = sqrt(sum(pow(A(i,_) - A(j,_),2.0)));
        DD(j,i) = DD(i,j);
      }
    }
  }
  
  List L = List::create(Named("D") = DD);
  return L;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
#projC(c(1,2,3), c(2,3,4))
*/

/*** R
# A = matrix(1:10, ncol = 2)
# B = matrix(1:10, nrow = 2)
# 
# my_matmult(A,B)
*/

/*** R
# myc_diag(matrix(c(3, 2, 3, 
#                   2, 5, 1, 
#                   5, 4, 9), ncol = 3))
*/

/*** R
# myc_qr(matrix(c(3, 2, 3, 
#                 2, 5, 1, 
#                 5, 4, 9), ncol = 3))
*/

/*** R
# myc_eigen(matrix(c(26, 40, 41,
#                    40, 67, 62,
#                    41, 62, 95), ncol = 3))
*/

/*** R
# myc_dist(matrix(c(1, 2, 3,
#                   2, 3, 4,
#                   4, 5, 6), ncol = 3))
*/


