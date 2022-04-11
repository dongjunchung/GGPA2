#if !defined(_Param_H)
#define _Param_H

#include <RcppArmadillo.h>
#include <sstream>      // For convert int to string

class CData {
  
public:
  
  CData() ; //constructor
  ~CData() ; //destructor
  bool Debug ;

  arma::mat Y, logY ;
  double theta_mu, tau2_mu ; 
  double a_sigma, b_sigma ;  
  double theta_alpha, tau2_alpha, stepsize_alpha ; 
  double a_beta, b_beta, stepsize_beta, stepsize_gamma ;
  double a_betaG, b_betaG ; 
  double threshold_on ; 	
	arma::mat E_forcein_mat ; bool isforcein ;
	
	bool have_annot ; 
	double a_gamma, b_gamma ; 
	arma::mat annotMat ;
		
	int msg_level ; // 0: errors only; 1: errors and warnings; 2: errors, warnings and information
	// bool OptionEit ; // TRUE: store E_it directly from c++
// private:
	  
};

class CParam {
	
	public:
  
	CParam(); 
	virtual ~CParam(); //Destructor
	
  int n_pheno, n_SNP, n_annot ; // constants but stored in CParam for convienence 
  
	arma::mat E_mat ;		// Naming - If the name is too short, attach _double / _vec / _mat 
	arma::vec mu_vec, sig2_vec ;
  arma::mat Beta, G_mat, gamma_mat, u_mat ;
  double p_u ; 
	
	// double normC ;
	arma::vec normC_vec, is_normC_q_type_updated ;
	
	// Result or log
  arma::vec accept_prob_vec, is_accept_vec ;
  double logPost, loglikelihood ; 
	arma::cube sum_E_ijt ; 
	
  // Initialize
  int is_initialized ; 
  void Initialize(CData &Data) ;
  void check_random_generate(CData &Data) ; 
  void clearE_ijt( ) ; 
    
  // MCMC
  void iterate(int iter, CData &Data, int n_simul);
  
  // Function
  double normC_fn(arma::mat Beta_input, CData &Data, arma::mat, arma::vec) ; 
  
private:
  
  // For random number 
  Rcpp::NumericVector RandVec ; 
  
  // MCMC steps  
  void S1_e_it(CData &Data);
	void S2_mu_i(CData &Data);
	void S3_sig2_i(CData &Data);
  void S4_alpha_i(CData &Data);
  void S5_beta_ij(CData &Data);
	void S6_G_beta_ij(CData &Data);
	void S7_gamma_i(CData &Data); 
	void S7b_u_gamma(CData &Data);
	void S8_p_u(CData &Data);
	
	void store_Eit(CData &Data); 
	
  // Distribution
  double rinvgamma(double, double);
  double rtruncNorm_lowertail_fn(double, double, double) ;
  double rtruncNorm_uppertail_fn(double, double, double) ;
	double dtruncnorm_uppertail_fn(double, double, double, double);
  int rDiscrete(int max_no); // generate an integer from 0 to (max_no-1)  
  
};

#endif
