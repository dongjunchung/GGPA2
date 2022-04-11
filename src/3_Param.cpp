#include "3_Param.h" 

#define LOG_2_PI 1.83787706640935

CData::CData(){ Debug = false ; }
CData::~CData(){ } //Destructor

CParam::CParam(){
  accept_prob_vec = arma::zeros<arma::vec>(9) ; 
  is_accept_vec = arma::zeros<arma::vec>(9) ; 
}
CParam::~CParam(){ }

void CParam::Initialize(CData &Data){
  n_pheno = Data.Y.n_rows ; n_SNP = Data.Y.n_cols ; // constants but stored in CParam for convienence
  Data.logY = arma::zeros<arma::mat>(n_pheno,n_SNP) ; 
		// Note 1. This is used to calculate sum of log y_it when e_it=1. 
		//      2. (Weak signal) Because e_it=0 with y_it <= 0 by definition (model assumption),
		//          we will not use log y_it for e_it=0 and store zero here. 
		//			3. (Strong signal) Also, we put e_it=1 when y_it > threshold    // V 1.3.1 
  for (int i_pheno=0; i_pheno<n_pheno; i_pheno++ ){
    for (int i_SNP=0; i_SNP<n_SNP; i_SNP++ ){
      if ( Data.Y(i_pheno,i_SNP) > 0){ // V 2.0.2
        Data.logY(i_pheno,i_SNP) = log(Data.Y(i_pheno,i_SNP)) ;
				if ( Data.Y(i_pheno,i_SNP) > Data.threshold_on ) E_mat(i_pheno,i_SNP) = 1 ; 
      } else {
      	E_mat(i_pheno,i_SNP) = 0 ;
      }
    }
  } 
	
	n_annot = gamma_mat.n_cols ; 
  normC_vec = arma::zeros<arma::vec>(pow(2,n_annot)) ; 
	for (int temp_type=0; temp_type<pow(2,n_annot); temp_type++){
		// std::cout << temp_type << std::endl ; 
		arma::vec temp_a_vec = arma::zeros<arma::vec>(n_annot) ; 
		int resid = temp_type ; 
		for (int i_m=0; i_m<n_annot; i_m++){
		  int temp_pow = n_annot-1-i_m ; 
		  if (resid>=pow(2,temp_pow)){
		    temp_a_vec(temp_pow) = 1 ; 
		    resid = resid - pow(2,temp_pow) ; 
		  }
		} // for 
		// std::cout << temp_a_vec << std::endl ; 
		normC_vec(temp_type) = normC_fn(Beta, Data, gamma_mat, temp_a_vec) ;
	  	// Ver_1_4_1
	  if ( normC_vec(temp_type) < 0 ){
	    Rcpp::stop("The initialized normC_vec(temp_type) has a negative value.") ;
	  }
	  Data.msg_level = 0; //0 errors only; 1: errors and warnings; 2: errors, warnings and information
		
	} // 
	// std::cout << normC_vec << std::endl ; 
	// std::cout << gamma_mat << std::endl ; 
	
  sum_E_ijt = arma::zeros<arma::cube>(n_pheno,n_pheno,n_SNP) ; 
	
	p_u = 0.5 ; 
  
  is_initialized = 1 ; 
}

void CParam::check_random_generate(CData &Data){
  RandVec = Rcpp::rnorm(2,2,1000) ;
  std::cout << RandVec(0) << "  " << RandVec(1) << std::endl ; 
  // std::cout << pow(10,-3) << std::endl ; 
}

void CParam::iterate(int iter, CData &Data, int n_simul) {
  if ( is_initialized == 1 ){
    logPost = 0.0 ; loglikelihood = 0.0 ; 
    S1_e_it(Data) ;
    S2_mu_i(Data) ;
    S3_sig2_i(Data) ;
    S4_alpha_i(Data) ;
  	S5_beta_ij(Data) ;
  	S6_G_beta_ij(Data) ;
  	S7_gamma_i(Data) ;
  	S7b_u_gamma(Data) ;
  	S8_p_u(Data) ;
    store_Eit(Data) ;
  } else {
    Rcpp::stop("Need To Run model$Initialize()") ; 
  }
} 

void CParam::clearE_ijt( ){
  sum_E_ijt = arma::zeros<arma::cube>(n_pheno,n_pheno,n_SNP) ; 
}

///////////
void CParam::S1_e_it(CData &Data) {
	
  for (int i_SNP=0; i_SNP<n_SNP; i_SNP++){
    for (int i_pheno=0; i_pheno<n_pheno; i_pheno++){
      if ( ( Data.Y(i_pheno,i_SNP) > 0 ) & ( Data.Y(i_pheno,i_SNP) <= Data.threshold_on ) ){ 
				// V 1.3.1 -> other case: e_it is fixed in CParam::Initialize, 
				//             i.e., e_it=0 if y_it<=0 and 1 if y_it>=Data.threshold_on
        double unnorm_logprob0 = 0.0 ;
        double unnorm_logprob1 = Beta(i_pheno,i_pheno) ;
        for (int m=0; m<n_annot; m++){
          unnorm_logprob1 = unnorm_logprob1 + gamma_mat(i_pheno,m) * Data.annotMat(m,i_SNP) ;
        }
        for (int j_pheno=0; j_pheno<n_pheno; j_pheno++){
          if (G_mat(i_pheno,j_pheno)==1) unnorm_logprob1 = unnorm_logprob1 + Beta(i_pheno,j_pheno) * E_mat(j_pheno,i_SNP) ; 
        }
				  // Note: Not count G_mat(i,j)=0 or G_mat(i,j)=9,i.e., diagonal
        unnorm_logprob0 = unnorm_logprob0 + R::dnorm(Data.Y(i_pheno,i_SNP),0,1,1) ; // log = T
        unnorm_logprob1 = unnorm_logprob1 + R::dlnorm(Data.Y(i_pheno,i_SNP),mu_vec(i_pheno),sqrt(sig2_vec(i_pheno)),1) ; // log = T
        double prob_e_it = 1.0 / ( 1.0 + exp(unnorm_logprob0-unnorm_logprob1) ) ; 
        	// Note: p1 = f1/(f1+f0) = 1 / (1+f0/f1) = 1 / (1+exp(log f0 - log f1)) 
				RandVec = Rcpp::runif(1,0,1) ;
        if ( prob_e_it >= RandVec(0) ){
          E_mat(i_pheno,i_SNP) = 1 ; 
        } else {
          E_mat(i_pheno,i_SNP) = 0 ; 
        }
      }         
      // To check loglikelihood
      if (E_mat(i_pheno,i_SNP)==1){
        loglikelihood = loglikelihood + R::dlnorm(Data.Y(i_pheno,i_SNP),mu_vec(i_pheno),sqrt(sig2_vec(i_pheno)),1)  ;
      } else {
        loglikelihood = loglikelihood + R::dnorm(Data.Y(i_pheno,i_SNP),0,1,1) ;
      }
		} 
  }
	
  is_accept_vec(0) = 1 ; 
} 

//////////
void CParam::S2_mu_i(CData &Data) {

	for (int i=0; i<n_pheno; i++){
		double sig2_i = sig2_vec(i) ;
		arma::vec e_i = E_mat.row(i).t() ; 
		double n_i = sum(e_i) ; 
		arma::vec logy_i = Data.logY.row(i).t() ;  
		arma::vec e1_logy = logy_i % e_i ; 
		double sum_e1_logy = sum(e1_logy) ; 
			// Note: % Schur product: element-wise multiplication of two objects
			//         sum of log_y_it with e_it=1
		double mean_star = (sig2_i * Data.theta_mu + Data.tau2_mu * sum_e1_logy) / (sig2_i + Data.tau2_mu * n_i) ; 
		double var_star = (sig2_i * Data.tau2_mu)/(sig2_i + Data.tau2_mu * n_i) ; 
		RandVec = Rcpp::rnorm(1, mean_star, sqrt(var_star)) ; 
		mu_vec(i) = RandVec(0) ;
  } 
	
	is_accept_vec(1) = 1 ; 
} 

///////////
void CParam::S3_sig2_i(CData &Data) {

	for (int i=0; i<n_pheno; i++){
		arma::vec mu_i_onevec(n_SNP) ; mu_i_onevec.fill(mu_vec(i)) ; 
		arma::vec e_i = E_mat.row(i).t() ; 
		double n_i = sum(e_i) ; 
		arma::vec logy_i = Data.logY.row(i).t() ;  
		arma::vec logy_minus_mu = logy_i - mu_i_onevec ;  
		arma::vec e1_logy_minus_mu = logy_minus_mu % e_i ;
		arma::vec sum_e1_squares = e1_logy_minus_mu.t() * e1_logy_minus_mu ; 
		double a_star = Data.a_sigma + 0.5 * n_i ; 
		double b_star = Data.b_sigma + 0.5 * sum_e1_squares(0) ; 
		sig2_vec(i) = rinvgamma(a_star, b_star) ; 
  } 

	is_accept_vec(2) = 1 ;
}

///////////
void CParam::S4_alpha_i(CData &Data) {

  is_accept_vec(3) = 0 ;
  for (int i_pheno=0; i_pheno<n_pheno; i_pheno++){  // V_1_3_5 ; Update each alpha_i
    arma::mat Beta_q = Beta ;
    RandVec = Rcpp::rnorm(1, Beta(i_pheno,i_pheno), Data.stepsize_alpha ) ;
    Beta_q(i_pheno,i_pheno) = RandVec(0) ;
    double logP_numer = R::dnorm(Beta_q(i_pheno,i_pheno), Data.theta_alpha, sqrt(Data.tau2_alpha),1) ; // log=TRUE
    double logP_denom = R::dnorm(Beta(i_pheno,i_pheno), Data.theta_alpha, sqrt(Data.tau2_alpha),1) ; // log=TRUE
    
    arma::vec normC_vec_q = normC_vec ; 
    is_normC_q_type_updated = arma::zeros<arma::vec>(pow(2,n_annot)) ; 
    for (int i_SNP=0; i_SNP<n_SNP; i_SNP++){
      arma::vec a_vec = Data.annotMat.col(i_SNP) ;
      int temp_type = 0 ; 
      for (int i_m=0; i_m<n_annot; i_m++){
        if (a_vec(i_m)==1) temp_type = temp_type + pow(2,i_m) ; 
      } // for
      if (is_normC_q_type_updated(temp_type)==0){
        normC_vec_q(temp_type) = normC_fn(Beta_q, Data, gamma_mat, a_vec) ;
        is_normC_q_type_updated(temp_type) = 1 ;
      } // if
      double e_it = E_mat(i_pheno,i_SNP) ;
      logP_numer = logP_numer + log( normC_vec(temp_type) ) + Beta_q(i_pheno,i_pheno) * e_it  ;
      logP_denom = logP_denom + log( normC_vec_q(temp_type) ) + Beta(i_pheno,i_pheno) * e_it ;
    } // for (i_SNP)
    double accept_prob = exp( logP_numer - logP_denom ) ;
    accept_prob_vec(3) = accept_prob ;
    RandVec = Rcpp::runif(1, 0, 1) ;
    if ( accept_prob >= RandVec(0) ){
      Beta = Beta_q ; normC_vec = normC_vec_q ;
      is_accept_vec(3) = is_accept_vec(3) + 1.0 / n_pheno ;
    }
  } // for (int i=0; i<n_pheno; i++){
  
}

///////////
void CParam::S5_beta_ij(CData &Data) {

	is_accept_vec(4) = 0 ;
	int w_cur = 0 ; // just for calculation of is_accept_vec( )
	for (int i=0; i<(n_pheno-1); i++){
		for (int j=(i+1); j<n_pheno; j++){
			if ( G_mat(i,j)==1 ){
				w_cur ++ ;
			  arma::mat Beta_q = Beta ;

				double temp_Beta_q = 0.0 ;
				while ( temp_Beta_q <= 0.0 ){
					temp_Beta_q = rtruncNorm_uppertail_fn(Beta(i,j), Data.stepsize_beta, 0) ;
				}
				Beta_q(i,j) = temp_Beta_q ;
				Beta_q(j,i) = temp_Beta_q ;
				double logP_numer = R::dgamma(Beta_q(i,j),Data.a_beta,(1.0/Data.b_beta),1) ;   // (shape,scale,log), i.e., b_beta = rate // V_1_3_5
			  double logP_denom = R::dgamma(Beta(i,j),Data.a_beta,(1.0/Data.b_beta),1) ;         //
        
        // double normC_prop = normC_fn(Beta_q, Data) ;
        arma::vec normC_vec_q = normC_vec ; 
        is_normC_q_type_updated = arma::zeros<arma::vec>(pow(2,n_annot)) ; 
			  for (int i_SNP=0; i_SNP<n_SNP; i_SNP++){
			    arma::vec a_vec = Data.annotMat.col(i_SNP) ;
			    int temp_type = 0 ; 
			    for (int i_m=0; i_m<n_annot; i_m++){
			      if (a_vec(i_m)==1) temp_type = temp_type + pow(2,i_m) ; 
			    } // for
			    if (is_normC_q_type_updated(temp_type)==0){
			      normC_vec_q(temp_type) = normC_fn(Beta_q, Data, gamma_mat, a_vec) ;
			      is_normC_q_type_updated(temp_type) = 1 ;
			    } // if
			    double e_it = E_mat(i,i_SNP) ;
			    double e_jt = E_mat(j,i_SNP) ;
			    logP_numer = logP_numer + log(normC_vec(temp_type)) + Beta_q(i,j) * e_it * e_jt ;
			    logP_denom = logP_denom + log(normC_vec_q(temp_type)) + Beta(i,j) * e_it * e_jt ;
			  } // for (i_SNP)
			  
			  double logQ_numer = log( dtruncnorm_uppertail_fn(Beta(i,j), 0, Beta_q(i,j), Data.stepsize_beta) ) ; // CHECK
			  double logQ_denom = log( dtruncnorm_uppertail_fn(Beta_q(i,j), 0, Beta(i,j), Data.stepsize_beta) ) ;
				double accept_prob = exp( logP_numer - logP_denom + logQ_numer - logQ_denom ) ;
			  accept_prob_vec(4) = accept_prob ;
        RandVec = Rcpp::runif(1, 0, 1) ;
			  if ( accept_prob >= RandVec(0) ){
			    Beta = Beta_q ; normC_vec = normC_vec_q ;
					is_accept_vec(4) = is_accept_vec(4) + 1 ;
			  }
			}
		}
	}

	if ( w_cur == 0 ){
	  is_accept_vec(4) = 0.2 ;
	} else {
	  is_accept_vec(4) = 1.0 / w_cur * is_accept_vec(4) ;
	}
}

///////////
void CParam::S6_G_beta_ij(CData &Data) {

  arma::mat G_mat_q = G_mat ; arma::mat Beta_q = Beta ;
  double logQ_numer, logQ_denom ;
  double logP_numer, logP_denom ;
  // double normC_prop ;

  bool is_forcein_edge_selected = false ;

	// double P_G_q, P_G ; // For updated E(i,j)

  // Step 1
  int w_cur = 0 ; int w_max = 0 ; int w_prop ;
  for (int i=0; i<(n_pheno-1); i++){
    for (int j=(i+1); j<n_pheno; j++){
      w_max ++ ; w_cur = w_cur + G_mat(i,j) ;
    }
  }
  if ( w_cur==0 ){ w_prop = 1 ; logQ_numer = log(0.5) ; logQ_denom = log(1.0) ; } // q(w|w^q) // q(w^q|w)
  if ( w_cur==w_max ){ w_prop = w_max - 1 ; logQ_numer = log(0.5) ;  logQ_denom = log(1.0) ; } // q(w|w^q) // q(w^q|w)
  if ( (w_cur > 0) && (w_cur < w_max) ){
    RandVec = Rcpp::runif(1, 0, 1) ;
    if ( RandVec(0) < 0.5 ){ w_prop = w_cur + 1 ; } else { w_prop = w_cur - 1 ; }
    logQ_denom = log(0.5) ; // q(w^q|w)
    if ( (w_prop==0) || (w_prop==w_max) ){ logQ_numer = log(1.0) ; } else { logQ_numer = log(0.5) ; } // q(w|w^q)
  }

  // Step 2
  arma::vec normC_vec_q = normC_vec ; 
  
  if ( w_prop > w_cur  ){ // step 2-a

		int id_added = rDiscrete(w_max-w_cur) + 1 ; // 1 ~ total. of empty edges
    int count_empty_edge = 0 ;

    for (int i=0; i<(n_pheno-1); i++){
      for (int j=(i+1); j<n_pheno; j++){
        if ( G_mat(i,j)==0 ){
          count_empty_edge++;
          if ( id_added==count_empty_edge ){
            G_mat_q(i,j) = 1 ; G_mat_q(j,i) = G_mat_q(i,j) ;
            RandVec = Rcpp::rgamma(1, Data.a_betaG, 1.0/Data.b_betaG) ; // q(beta_ij^q)
            Beta_q(i,j) = RandVec(0) ;
            Beta_q(j,i) = Beta_q(i,j) ;
            id_added = -9 ;

						// if ( Data.PriorSetting==2 ){
						//	P_G = 1.0 - Data.priorprob_G(i,j) ; // Bernoulli(E_ij=0; p_ij)
						//	P_G_q = Data.priorprob_G(i,j) ; // Bernoulli(E_ij^q=1; p_ij)
						// }

            logP_numer = R::dgamma(Beta_q(i,j),Data.a_beta,(1.0/Data.b_beta),1) ; // f(beta_ij^q|G_ij^q)
            logP_denom = 0 ; // f(beta_ij|G_ij) cancelled with q(beta_ij)
            
            // normC_prop = normC_fn(Beta_q, Data) ;
            // arma::vec normC_vec_q = normC_vec ; 
            is_normC_q_type_updated = arma::zeros<arma::vec>(pow(2,n_annot)) ; 
            for (int i_SNP=0; i_SNP<n_SNP; i_SNP++){
              arma::vec a_vec = Data.annotMat.col(i_SNP) ;
              int temp_type = 0 ; 
              for (int i_m=0; i_m<n_annot; i_m++){
                if (a_vec(i_m)==1) temp_type = temp_type + pow(2,i_m) ; 
              } // for
              if (is_normC_q_type_updated(temp_type)==0){
                normC_vec_q(temp_type) = normC_fn(Beta_q, Data, gamma_mat, a_vec) ;
                is_normC_q_type_updated(temp_type) = 1 ;
              } // if
              double e_it = E_mat(i,i_SNP) ;
              double e_jt = E_mat(j,i_SNP) ;
              logP_numer = logP_numer + log(normC_vec(temp_type)) + Beta_q(i,j) * e_it * e_jt ; // f(e_t|alpha,beta^q,G^q)
              logP_denom = logP_denom + log(normC_vec_q(temp_type)) + Beta(i,j) * e_it * e_jt ; // f(e_t|alpha,beta,G)
            } // for (i_SNP) 
            
            logQ_numer = logQ_numer + 0 ; // q(beta_ij) cancelled with f(beta_ij|G_ij)
            logQ_denom = logQ_denom + R::dgamma(Beta_q(i,j),Data.a_betaG,(1.0/Data.b_betaG),1) ; // q(beta_ij^q)
          } // if ( id_added==count_empty_edge )
        } // for ( G_mat(i,j)==0 )
      } // for (j)
    } // for (i)

    logQ_numer = logQ_numer - log(w_cur) ; // q(G|G^q,w)
    logQ_denom = logQ_denom - log(w_max-w_cur)  ; // q(G^q|G,w^q)

	} else { // step 2-b

    int id_deleted = rDiscrete(w_cur)+1 ; // 1 ~ total no. of connected edges
    int count_connected_edge = 0 ;
    for (int i=0; i<(n_pheno-1); i++){
      for (int j=(i+1); j<n_pheno; j++){
        if ( G_mat(i,j)==1 ){
          count_connected_edge++;
          if ( id_deleted==count_connected_edge ){
            G_mat_q(i,j) = 0 ; G_mat_q(j,i) = G_mat_q(i,j) ;
            Beta_q(i,j) = 0 ; // q(beta_ij^q)
            Beta_q(j,i) = Beta_q(i,j) ;
            id_deleted = -9 ;

            if (Data.isforcein==true){
              if (Data.E_forcein_mat(i,j)==1) is_forcein_edge_selected = true ;
            }

						// if ( Data.PriorSetting==2 ){
						//	P_G = Data.priorprob_G(i,j) ; // Bernoulli(E_ij=1; p_ij)
						//	P_G_q = 1.0 - Data.priorprob_G(i,j) ; // Bernoulli(E_ij^q=0; p_ij)
						// }

            logP_numer = 0 ; // f(beta_ij^q|G_ij^q) cancelled with q(beta_ij^q)
            logP_denom = R::dgamma(Beta(i,j),Data.a_beta,(1.0/Data.b_beta),1) ; // f(beta_ij|G_ij)
            // normC_prop = normC_fn(Beta_q, Data) ;
            
            is_normC_q_type_updated = arma::zeros<arma::vec>(pow(2,n_annot)) ; 
            for (int i_SNP=0; i_SNP<n_SNP; i_SNP++){
              arma::vec a_vec = Data.annotMat.col(i_SNP) ;
              int temp_type = 0 ; 
              for (int i_m=0; i_m<n_annot; i_m++){
                if (a_vec(i_m)==1) temp_type = temp_type + pow(2,i_m) ; 
              } // for
              if (is_normC_q_type_updated(temp_type)==0){
                normC_vec_q(temp_type) = normC_fn(Beta_q, Data, gamma_mat, a_vec) ;
                is_normC_q_type_updated(temp_type) = 1 ;
              } // if
              double e_it = E_mat(i,i_SNP) ;
              double e_jt = E_mat(j,i_SNP) ;
              logP_numer = logP_numer + log(normC_vec(temp_type)) + Beta_q(i,j) * e_it * e_jt ; // f(e_t|alpha,beta^q,G^q)
              logP_denom = logP_denom + log(normC_vec_q(temp_type)) + Beta(i,j) * e_it * e_jt ; // f(e_t|alpha,beta,G)
            }
            logQ_numer = logQ_numer + R::dgamma(Beta(i,j),Data.a_betaG,(1.0/Data.b_betaG),1) ; // q(beta_ij)
            logQ_denom = logQ_denom + 0 ; // q(beta_ij^q) cancelled with f(beta_ij^q|G_ij^q)
          }
        }
      }
    }
    logQ_numer = logQ_numer - log(w_max-w_cur) ; // q(G|G^q,w)
    logQ_denom = logQ_denom - log(w_cur) ; // q(G^q|G,w^q)

  }

	// if ( Data.PriorSetting==1 ){
	  if ( w_prop > 0 ) logP_numer = logP_numer - log(w_prop) ; // f(G^q) // if w_prop=0, let 1/w_prop = 1, so that log(w_prop) = 0
	  if ( w_cur > 0 ) logP_denom = logP_denom - log(w_cur) ; // f(G) // if w_cur=0, let 1/w_cur = 1, so that log(w_cur) = 0
	// }
	// if ( Data.PriorSetting==2 ){
	//  if ( w_prop > 0 ) logP_numer = logP_numer + log(P_G_q) ;
	//  if ( w_cur > 0 ) logP_denom = logP_denom + log(P_G) ;
	// }

  // Step 3
  double accept_prob = exp( logP_numer - logP_denom + logQ_numer - logQ_denom ) ;

  if (is_forcein_edge_selected==true) accept_prob = 0 ;

  accept_prob_vec(5) = accept_prob ;
  RandVec = Rcpp::runif(1, 0, 1) ;
  if ( accept_prob >= RandVec(0) ){
    Beta = Beta_q ; G_mat = G_mat_q ; normC_vec = normC_vec_q ;
    is_accept_vec(5) = 1 ;
  } else {
    is_accept_vec(5) = 0 ;
  }
  
} // void CParam::S6_G_beta_ij(CData &Data)

// ///////////
// void CParam::S7_gamma_i(CData &Data) {
//
//   is_accept_vec(6) = 0 ;
//
//   for (int i_pheno=0; i_pheno<n_pheno; i_pheno++){
//
//     arma::mat gamma_mat_q = gamma_mat ;
//
//     double logP_numer = 0 ; double logP_denom = 0 ;
//     for (int i_m=0; i_m<n_annot; i_m++){
//       gamma_mat_q(i_pheno,i_m) = rtruncNorm_uppertail_fn(gamma_mat(i_pheno,i_m), sqrt(Data.stepsize_gamma), 0) ; // mu, sigma, lowerlimit
//       logP_numer = logP_numer + log( dtruncnorm_uppertail_fn(gamma_mat_q(i_pheno,i_m), Data.theta_gamma, sqrt(Data.tau2_gamma), 0) ) ;
//       logP_denom = logP_denom + log( dtruncnorm_uppertail_fn(gamma_mat(i_pheno,i_m), Data.theta_gamma, sqrt(Data.tau2_gamma), 0) ) ;
//       logP_numer = logP_numer + log( dtruncnorm_uppertail_fn(gamma_mat(i_pheno,i_m), gamma_mat_q(i_pheno,i_m), sqrt(Data.stepsize_gamma), 0) ) ;
//       logP_denom = logP_denom + log( dtruncnorm_uppertail_fn(gamma_mat_q(i_pheno,i_m), gamma_mat(i_pheno,i_m), sqrt(Data.stepsize_gamma), 0) ) ;
//     } // i_m
//
//     arma::vec normC_vec_q = normC_vec ;
//     is_normC_q_type_updated = arma::zeros<arma::vec>(pow(2,n_annot)) ;
//     for (int i_SNP=0; i_SNP<n_SNP; i_SNP++){
//       arma::vec a_vec = Data.annotMat.col(i_SNP) ;
//       // log C - log C_q
//       int temp_type = 0 ;
//       for (int i_m=0; i_m<n_annot; i_m++){
//         if (a_vec(i_m)==1) temp_type = temp_type + pow(2,i_m) ;
//       } // for
//       if (is_normC_q_type_updated(temp_type)==0){
//         normC_vec_q(temp_type) = normC_fn(Beta, Data, gamma_mat_q, a_vec) ;
//         is_normC_q_type_updated(temp_type) = 1 ;
//       } // if
//       logP_numer = logP_numer + log( normC_vec(temp_type) )  ;
//       logP_denom = logP_denom + log( normC_vec_q(temp_type) ) ;
//       // sum gamma_q a e - sum gamma a e
//       double e_it = E_mat(i_pheno,i_SNP) ;
//       if (e_it==1){
//         for (int m=0; m<n_annot; m++){
//           logP_numer = logP_numer + gamma_mat_q(i_pheno,m) * Data.annotMat(m,i_SNP) ;
//           logP_denom = logP_denom + gamma_mat(i_pheno,m) * Data.annotMat(m,i_SNP) ;
//         }
//       } // if (e_it)
//     } // for (i_SNP)
//
//     double accept_prob = exp( logP_numer - logP_denom ) ;
//     accept_prob_vec(6) = accept_prob ;
//     RandVec = Rcpp::runif(1, 0, 1) ;
//     if ( accept_prob >= RandVec(0) ){
//       gamma_mat = gamma_mat_q ; normC_vec = normC_vec_q ;
//       is_accept_vec(6) = is_accept_vec(6) + 1.0 / n_pheno ;
//     }
//   } // for (int i=0; i<n_pheno; i++){
//
// } // void CParam::S7_gamma_i(CData &Data)

///////////
void CParam::S7_gamma_i(CData &Data) {

  is_accept_vec(6) = 0 ;

  int count_attempt = 0 ; int count_success = 0 ; 
  
  for (int i_pheno=0; i_pheno<n_pheno; i_pheno++){
		for (int i_m=0; i_m<n_annot; i_m++){
		  
			if ( u_mat(i_pheno,i_m)==1 ){
			  
			  count_attempt = count_attempt + 1 ; 
			  
			  arma::mat gamma_mat_q = gamma_mat ;
			  double logP_numer = 0 ; double logP_denom = 0 ;
			  
			  gamma_mat_q(i_pheno,i_m) = rtruncNorm_uppertail_fn(gamma_mat(i_pheno,i_m), sqrt(Data.stepsize_gamma), 0) ; // mu, sigma, lowerlimit
			  
			  logP_numer = logP_numer + R::dgamma(gamma_mat_q(i_pheno,i_m),Data.a_gamma,(1.0/Data.b_gamma),1) ;   // (shape,scale,log), i.e., b_beta = rate // V_1_3_5
			  logP_denom = logP_denom + R::dgamma(gamma_mat(i_pheno,i_m),Data.a_gamma,(1.0/Data.b_gamma),1) ;   // (shape,scale,log), i.e., b_beta = rate // V_1_3_5
			  
			  logP_numer = logP_numer + log( dtruncnorm_uppertail_fn(gamma_mat(i_pheno,i_m), gamma_mat_q(i_pheno,i_m), sqrt(Data.stepsize_gamma), 0) ) ;
			  logP_denom = logP_denom + log( dtruncnorm_uppertail_fn(gamma_mat_q(i_pheno,i_m), gamma_mat(i_pheno,i_m), sqrt(Data.stepsize_gamma), 0) ) ;
			  
			  arma::vec normC_vec_q = normC_vec ;
			  is_normC_q_type_updated = arma::zeros<arma::vec>(pow(2,n_annot)) ;
			  for (int i_SNP=0; i_SNP<n_SNP; i_SNP++){
			    arma::vec a_vec = Data.annotMat.col(i_SNP) ;
			    // log C - log C_q 
			    int temp_type = 0 ;
			    for (int i_m=0; i_m<n_annot; i_m++){
			      if (a_vec(i_m)==1) temp_type = temp_type + pow(2,i_m) ;
			    } // for
			    if (is_normC_q_type_updated(temp_type)==0){
			      normC_vec_q(temp_type) = normC_fn(Beta, Data, gamma_mat_q, a_vec) ;
			      is_normC_q_type_updated(temp_type) = 1 ;
			    } // if
			    logP_numer = logP_numer + log( normC_vec(temp_type) )  ;
			    logP_denom = logP_denom + log( normC_vec_q(temp_type) ) ;
			    logP_numer = logP_numer + gamma_mat_q(i_pheno,i_m) * Data.annotMat(i_m,i_SNP) * E_mat(i_pheno,i_SNP) ;
			    logP_denom = logP_denom + gamma_mat(i_pheno,i_m) * Data.annotMat(i_m,i_SNP) * E_mat(i_pheno,i_SNP) ;
			  } // for (i_SNP)
			  
			  double accept_prob = exp( logP_numer - logP_denom ) ;
			  accept_prob_vec(6) = accept_prob ;
			  RandVec = Rcpp::runif(1, 0, 1) ;
			  if ( accept_prob >= RandVec(0) ){
			    gamma_mat = gamma_mat_q ; normC_vec = normC_vec_q ;
			    count_success = count_success + 1 ; 
			  }
			  
			} // if ( u_mat(i_pheno,i_m)==1 )
		} // for (int i_m) 
  } // for (int i=0; i<n_pheno; i++){

  if (count_attempt > 0){
    is_accept_vec(6) = count_success / count_attempt ;
  }
  
} // void CParam::S7_gamma_i(CData &Data)


///////////
void CParam::S7b_u_gamma(CData &Data) {
  
  is_accept_vec(7) = 0 ;
  
  for (int i_pheno=0; i_pheno<n_pheno; i_pheno++){
    for (int i_m=0; i_m<n_annot; i_m++){
      
      arma::mat gamma_mat_q = gamma_mat ;
      arma::mat u_mat_q = u_mat ;
      double logP_numer = 0 ; double logP_denom = 0 ; 
      u_mat_q(i_pheno,i_m) = 1 - u_mat(i_pheno,i_m) ; 
      
      if ( u_mat_q(i_pheno,i_m) == 1 ){
        RandVec = Rcpp::rgamma(1, Data.a_gamma,(1.0/Data.b_gamma)) ; // q(beta_ij^q)
        gamma_mat_q(i_pheno,i_m) = RandVec(0) ;
        // logP_numer = logP_numer + R::dgamma(gamma_mat_q(i_pheno,i_m),Data.a_gamma,(1.0/Data.b_gamma),1) ; // logP_denom = ++ log(1) i.e., add nothing 
        // logP_denom = logP_denom + R::dgamma(gamma_mat_q(i_pheno,i_m),1.0,1.0,1) ; 
      } else { 
        gamma_mat_q(i_pheno,i_m) = 0 ;
        // logP_denom = logP_denom + R::dgamma(gamma_mat(i_pheno,i_m),Data.a_gamma,(1.0/Data.b_gamma),1) ; // logP_numer = ++ log(1) i.e., add nothing 
        // logP_numer = logP_numer + R::dgamma(gamma_mat_q(i_pheno,i_m),1.0,1.0,1) ; 
      }

      logP_numer = logP_numer + u_mat_q(i_pheno,i_m) * log(p_u) + (1-u_mat_q(i_pheno,i_m)) * log(1-p_u) ;
      logP_denom = logP_denom + u_mat(i_pheno,i_m) * log(p_u) + (1-u_mat(i_pheno,i_m)) * log(1-p_u) ;
      
      arma::vec normC_vec_q = normC_vec ;
      is_normC_q_type_updated = arma::zeros<arma::vec>(pow(2,n_annot)) ;
      for (int i_SNP=0; i_SNP<n_SNP; i_SNP++){
        arma::vec a_vec = Data.annotMat.col(i_SNP) ;
        // log C - log C_q 
        int temp_type = 0 ;
        for (int i_m=0; i_m<n_annot; i_m++){
          if (a_vec(i_m)==1) temp_type = temp_type + pow(2,i_m) ;
        } // for
        if (is_normC_q_type_updated(temp_type)==0){
          normC_vec_q(temp_type) = normC_fn(Beta, Data, gamma_mat_q, a_vec) ;
          is_normC_q_type_updated(temp_type) = 1 ;
        } // if
        logP_numer = logP_numer + log( normC_vec(temp_type) )  ;
        logP_denom = logP_denom + log( normC_vec_q(temp_type) ) ;
        // sum gamma_q a e - sum gamma a e
        logP_numer = logP_numer + gamma_mat_q(i_pheno,i_m) * Data.annotMat(i_m,i_SNP) * E_mat(i_pheno,i_SNP) ;
        logP_denom = logP_denom + gamma_mat(i_pheno,i_m) * Data.annotMat(i_m,i_SNP) * E_mat(i_pheno,i_SNP) ;
      } // for (i_SNP)
      
      double accept_prob = exp( logP_numer - logP_denom ) ;
      accept_prob_vec(7) = accept_prob ;
      RandVec = Rcpp::runif(1, 0, 1) ;
      if ( accept_prob >= RandVec(0) ){
        gamma_mat = gamma_mat_q ; u_mat = u_mat_q ; normC_vec = normC_vec_q ;
        is_accept_vec(7) = is_accept_vec(7) + 1.0 / ( n_pheno * n_annot ) ;
      }
      
    } // i_m
  } // for (int i=0; i<n_pheno; i++){
  
} // void CParam::S7b_u_gamma(CData &Data)


///////////
void CParam::S8_p_u(CData &Data) {
  
  double sum_u_mat = accu(u_mat) ;
  RandVec = Rcpp::rbeta(1, 1+sum_u_mat, 1+(n_pheno*n_annot-sum_u_mat)) ;
  p_u = RandVec(0) ; 
  
} // 

// ///////////
void CParam::store_Eit(CData &Data) {
  for (int i_pheno=0; i_pheno<n_pheno; i_pheno++){
    for (int j_pheno=0; j_pheno<n_pheno; j_pheno++){
      for (int i_SNP=0; i_SNP<n_SNP; i_SNP++){
        sum_E_ijt(i_pheno,j_pheno,i_SNP) = sum_E_ijt(i_pheno,j_pheno,i_SNP) + E_mat(i_pheno,i_SNP) * E_mat(j_pheno,i_SNP) ; 
      }   
    }   
  }
}

// void CParam::store_Eit(CData &Data) {
//   
//   FILE *FILE_E_it; 
//   
//   std::string Name1 = "E_it_1" ; 
//   Name1.append(".dat") ; 
//   
//   FILE_E_it = fopen(Name1.c_str(),"a");
//   for (int i_SNP=0; i_SNP<n_SNP; i_SNP++){
//     fprintf(FILE_E_it, "%.0f ",E_mat(0,i_SNP));
//   }
//   fprintf(FILE_E_it, "\n");
//   fclose(FILE_E_it);
//   
//   std::string Name2 = "E_it_2" ; 
//   Name2.append(".dat") ; 
//   
//   FILE_E_it = fopen(Name2.c_str(),"a");
//   for (int i_SNP=0; i_SNP<n_SNP; i_SNP++){
//     fprintf(FILE_E_it, "%.0f ",E_mat(0,i_SNP));
//   }
//   fprintf(FILE_E_it, "\n");
//   fclose(FILE_E_it);
//   
// }
  
//////////////////////////////////////
// Function

// normC_vec(temp_type) = normC_fn(Beta, Data, gamma_mat, temp_a_vec) ;

double CParam::normC_fn(arma::mat Beta_input, CData &Data, arma::mat gamma_mat, arma::vec temp_a_vec){
  
    double normC_temp = 1.0 ; // exp(0) 
    arma::vec e_temp(n_pheno) ; 
    int n_pow_pheno = pow(2,n_pheno) ;
    
    for (int i_pow=1; i_pow<n_pow_pheno; i_pow++){ // exclude (0 0 0)
      
      int temp_no = i_pow ; 
      for (int i_pheno=0; i_pheno<n_pheno; i_pheno++){
        int temp_denominator = pow(2,(n_pheno-i_pheno-1)) ; 
        if ( floor(temp_no/temp_denominator) == 1 ){
          e_temp(i_pheno) = 1 ; 
          temp_no = temp_no - temp_denominator ; 
        } else {
          e_temp(i_pheno) = 0 ; 
        }
      } // for i_pheno
      
      double temp_beta_sum = 0.0 ; 
      for (int i_pheno=0; i_pheno<n_pheno; i_pheno++){
        if (e_temp(i_pheno)==1){
          arma::vec temp_vec = gamma_mat.row(i_pheno) * temp_a_vec ; 
          temp_beta_sum = temp_beta_sum + Beta_input(i_pheno,i_pheno) + temp_vec(0) ; 
          if (i_pheno<(n_pheno-1)){
            for (int j_pheno=(i_pheno+1); j_pheno<n_pheno; j_pheno++){
              if (e_temp(j_pheno)==1){
                temp_beta_sum = temp_beta_sum + Beta_input(i_pheno,j_pheno) ; 
              } // if (e_temp(j_pheno)==1)
            } // for j_pheno 
          } // if (i_pheno<(n_pheno-1))
        } // if (e_temp(i_pheno)==1)
      } // for i_pheno
      
      normC_temp = normC_temp + exp(temp_beta_sum) ; 
      
    } // for i_pow
    
    return normC_temp ;
} 

//////////////////////////////////////
// Distribution

double CParam::rinvgamma(double alpha, double beta){
  RandVec = Rcpp::rgamma(1, alpha, (1.0/beta) ) ; 
  return ( 1.0 / RandVec(0) ) ; 
  // Note that b in Rcpp::rgamma(a,b) is scale, 
  // i.e., its mean is ab, NOT a/b  
  // If X ~ Gamma(a,b), then 1/X ~ IG(a,1/b)
}

double CParam::rtruncNorm_lowertail_fn(double mu, double sigma, double upperlimit){
  RandVec = Rcpp::runif(1,0,1) ; 
  double std_upperlimit = (upperlimit-mu) / sigma ; 
  double temp_pnorm = R::pnorm(std_upperlimit,0,1,1,0) * RandVec(0) ; // x,mu,sigma,lt,lg
  double std_x = R::qnorm(temp_pnorm,0,1,1,0) ;  // p,mu,sigma,lt,lg
  return( mu + std_x * sigma ) ;
}

double CParam::rtruncNorm_uppertail_fn(double mu, double sigma, double lowerlimit){
  RandVec = Rcpp::runif(1,0,1) ; 
  double std_lowerlimit = (lowerlimit-mu) / sigma ; 
  double temp_pnorm = (1.0-R::pnorm(std_lowerlimit,0,1,1,0)) * RandVec(0) ; // x,mu,sigma,lt,lg
  double std_x = -1.0 * R::qnorm(temp_pnorm,0,1,1,0) ;  // p,mu,sigma,lt,lg
  return( mu + std_x * sigma ) ;
}

double CParam::dtruncnorm_uppertail_fn(double x, double mu, double sigma, double lowerlimit){
	double alpha = (lowerlimit-mu) / sigma ; 
	double z = (x-mu) / sigma ;
	double res = R::dnorm(z,0,1,0) / (sigma * ( 1.0-R::pnorm(alpha,0,1,1,0) ) ) ; 	
  return ( res ) ;
	
}

int CParam::rDiscrete(int max_no){
  // generate an integer from 0 to (max_no-1)
  double max_no_ = 1.0 * max_no ; 
  RandVec = Rcpp::runif(1, 0.0, max_no_) ; 
  return ( floor(RandVec(0)) ) ; 
}
