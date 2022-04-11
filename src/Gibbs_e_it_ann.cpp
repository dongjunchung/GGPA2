// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "RcppArmadillo.h"
#include <Rmath.h>

// [[Rcpp::export]]
arma::mat Gibbs_e_it_ann(arma::mat beta_mat, arma::mat G_mat, int n_SNP, int n_GibbsStep, arma::mat true_Gamma, arma::mat A_mat){
  
  int n_pheno = beta_mat.n_rows ; // i=1,...,n_pheno,  t=1,...,n_SNP
  int n_annot = true_Gamma.n_cols ; // i=1,...,n_pheno,  t=1,...,n_SNP
	
  arma::mat e_it(n_pheno,n_SNP) ; // filled with 0
  
  for (int i_iter=1; i_iter<=n_GibbsStep; i_iter++){
    int no_eit_1 = 0 ; 
    for (int t=0; t<n_SNP; t++){
      for (int i=0; i<n_pheno; i++){
        double unnorm_logprob1 = beta_mat(i,i) ;
				for (int m=0; m<n_annot; m++){
					unnorm_logprob1 = unnorm_logprob1 + true_Gamma(i,m) * A_mat(m,t) ; 
				} 
        double unnorm_logprob0 = 0.0 ; 
        for (int j=0; j<n_pheno; j++){
          if (G_mat(i,j)==1){
            unnorm_logprob1 = unnorm_logprob1 + beta_mat(i,j) * 1.0 * e_it(j,t) ; 
          }
        }
        double prob_e_it = 1.0 / ( 1.0 + exp(unnorm_logprob0-unnorm_logprob1) ) ; 
        if ( prob_e_it >= R::runif(0,1) ){
          e_it(i,t) = 1 ; no_eit_1 ++ ; 
        } else {
          e_it(i,t) = 0 ; 
        }
      } // i=0~(n_pheno-1)
    } // t=0~(n_SNP-1)
    if (i_iter%100==0){
      Rprintf("iter: %d", i_iter);
      Rprintf(", total no. of e_it=1: %d \n", no_eit_1); 
    }
  } // i_iter 
    
  return e_it ; 
  
}
