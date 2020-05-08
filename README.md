# myc

`myc` is an R package with several functions wiritten in C++. The purpous was for me to a. learn C++ and Rcpp b. revise some linear algebra and impliment matrix factorizations commonly used in data science.

## Try this out!

* point 1
* point 2
* point 3

ok..


## Functions

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

`myc_diag`

`myc_dist`

`myc_qr`

`myc_eigen`

`myc_svd`
