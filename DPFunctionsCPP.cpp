#include <iostream>
#include <cstdlib>
#include <string>
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h>
#include <RcppArmadilloExtensions/sample.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma; 
using namespace Rcpp;
using namespace std; 

// [[Rcpp::export]]
vec like(vec Z, vec thetah, vec Sticks, vec beta, mat data){
  int n = thetah.size(); // this code creates the multinomial probability for cluster assignment
  vec store(n); double q; 
  for (int i = 0; i < n; i++){
    q = accu(pow((Z - data*beta - thetah(i)), 2));
    store(i) = std::exp(-.5*q);
  }
  return(Sticks%store);
}

// [[Rcpp::export]]
vec like0(vec Z, vec thetah, vec Sticks){
  int n = thetah.size(); // this code creates the multinomial probability for cluster assignment
  vec store(n); double q; 
  for (int i = 0; i < n; i++){
    q = accu(pow((Z - thetah(i)), 2));
    store(i) = std::exp(-.5*q);
  }
  return(Sticks%store);
}

// This function creates the multinomial likelihood numerator that includes an option for additional covariates
// [[Rcpp::export]]
mat multinomprobCPP(vec Z, vec Gene, int colf, vec thetah, vec Sticks, vec beta, mat data, int G, int H){
  mat S_i(H, G); S_i.fill(0); 
  vec betaCov(colf);
  for (int q = 0; q < colf; q++){ // this gets the betas for fixed effects
    betaCov(q) = beta(q);
  }
  
  for(int k = 0; k < G; k++){
    uvec indic = find(Gene == k);
    vec Zcol = Z.elem(indic); 
    mat dataG = data.rows(indic);
    S_i.col(k) = like(Zcol,thetah,Sticks,betaCov,dataG);
  }
  return(S_i);
}

// This function creates the multinomial likelihood numerator that includes an option for additional covariates
//When we don't have covariates
// [[Rcpp::export]]
mat multinomprob0CPP(vec Z, vec Gene, vec thetah, vec Sticks, int G, int H){
  mat S_i(H, G); S_i.fill(0); 
  for(int k = 0; k < G; k++){
    uvec indic = find(Gene == k);
    vec Zcol = Z.elem(indic); 

    S_i.col(k) = like0(Zcol,thetah,Sticks);
  }
  return(S_i);
}
