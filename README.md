# myc

`myc` is an R package with several functions wiritten in C++. The purpous was for me to 
a. learn C++ and Rcpp 
b. revise some linear algebra and impliment matrix factorizations commonly used in data science.

# Download if you want

You can download this little package if you would like using devtools. Obviously this will not provide with any added functionality for R, and for the most part you are better of using the built-in functions that that use LINPACK for linear algebra algorithms.

```r
devtools::install_github("bdshaff/myc")

```

# See how I use it

On my website bdshaff.github.io I authored a couple blog posts that talk about linear algebra and discuss the computation of matrix decompositions along with a couple applications. I use the `myc` package there so you can take a look at that if you want.

## Functions

The main objective for me was to impliment the SVD matrix decomposition that is so central in Data Science. Along the way the QR decompopistion, and the Eigenvalue decompositions naturally happened too. Because I was not using anything beyond Rcpp like RcppArmadillo or RcppEigenI had to write a couple other functions along the way to do matrix multipliaction and extracting the diagonal of a matrix along with a couple other functions.

Here is the list of functions and their code. You can find the actual C++ code in the src/ directory.

`myc_proj`

```c++
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
```

`myc_matmult`

```c++
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
```

`myc_diag`

```c++
NumericVector myc_diag(NumericMatrix A){
  int n = A.nrow();
  NumericVector vec(n);
  
  for(int i = 0; i < n; i++){
    vec[i] = A(i,i);
  }
  
  return vec;
}
```

`myc_dist`
```c++
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
```

`myc_qr`
```c++
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
```

`myc_eigen`
```c++
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
```

`myc_svd`
```c++
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
```
