#include "1_Main.h"

RCPP_MODULE(cGGPA2module){
  
  using namespace R;
  using namespace Rcpp;
  
  class_<CMain>( "cGGPA2" )
  
  .constructor<arma::mat>()     // expose the default constructor
  .method("Initialize", &CMain::Initialize, "Initialize" )
  .method("check_random_generate", &CMain::check_random_generate, "check_random_generate" )
  .method("Iterate", &CMain::Iterate, "Run one iteration of MCMC algorithm" )
  .method("Run", &CMain::Run, "Run MCMC algorithm for given times" )
  .property("msg.level", &CMain::GetMsgLevel, &CMain::SetMsgLevel, "Message Level")
  .method("clear.sum_E_ijt", &CMain::clearE_ijt, "clear.sum_E_ijt" )
  
	// Initial values
  .property("E_mat", &CMain::GetE_mat, &CMain::SetE_mat, "Latent Biological Signal")
	.property("mu_vec", &CMain::Getmu_vec, &CMain::Setmu_vec, "mu_vec, parameter for N(y_it;mu_i,sig2_i)")
	.property("sig2_vec", &CMain::Getsig2_vec, &CMain::Setsig2_vec, "sig2_vec, parameter for N(y_it;mu_i,sig2_i)")
  .property("Beta", &CMain::GetBeta, &CMain::SetBeta, "Autologistic Model Coefficients")
  .property("G_mat", &CMain::GetG_mat, &CMain::SetG_mat, "MRF Graph")
  
	// Hyperparameter 
	.property("theta_mu", &CMain::Gettheta_mu,  &CMain::Settheta_mu, "Get theta_mu")
	.property("tau2_mu", &CMain::Gettau2_mu,  &CMain::Settau2_mu, "Get tau2_mu")
	.property("a_sigma", &CMain::Geta_sigma, &CMain::Seta_sigma, "Get a_sigma")
	.property("b_sigma", &CMain::Getb_sigma, &CMain::Setb_sigma, "Get b_sigma")
	.property("theta_alpha", &CMain::Gettheta_alpha,  &CMain::Settheta_alpha, "Get theta_alpha")
	.property("tau2_alpha", &CMain::Gettau2_alpha,  &CMain::Settau2_alpha, "Get tau2_alpha")
	.property("stepsize_alpha", &CMain::Getstepsize_alpha,  &CMain::Setstepsize_alpha, "Get stepsize_alpha")
	.property("a_beta", &CMain::Geta_beta,  &CMain::Seta_beta, "Get a_beta")
	.property("b_beta", &CMain::Getb_beta,  &CMain::Setb_beta, "Get b_beta")
	.property("stepsize_beta", &CMain::Getstepsize_beta,  &CMain::Setstepsize_beta, "Get stepsize_beta")
	.property("a_betaG", &CMain::Geta_betaG,  &CMain::Seta_betaG, "Get a_betaG")
	.property("b_betaG", &CMain::Getb_betaG,  &CMain::Setb_betaG, "Get b_betaG")
  .property("threshold_on", &CMain::Getthreshold_on,  &CMain::Setthreshold_on, "Get threshold_on")
	.property("E_forcein_mat", &CMain::GetE_forcein_mat, &CMain::SetE_forcein_mat, "Forced-in graph structure")
  .property("stepsize_gamma", &CMain::Getstepsize_gamma,  &CMain::Setstepsize_gamma, "Get stepsize_gamma")	
	
	.property("a_gamma", &CMain::Geta_gamma,  &CMain::Seta_gamma, "Get a_gamma")
	.property("b_gamma", &CMain::Getb_gamma,  &CMain::Setb_gamma, "Get b_gamma")
	.property("annotMat", &CMain::GetannotMat,  &CMain::SetannotMat, "Get annotMat")
	.property("have_annot", &CMain::Gethave_annot,  &CMain::Sethave_annot, "Get have_annot")
	.property("gamma_mat", &CMain::Getgamma_mat,  &CMain::Setgamma_mat, "Get gamma_mat")
  .property("u_mat", &CMain::Getu_mat,  &CMain::Setu_mat, "Get u_mat")
	
	// Print data or result 
	.property("Y.input", &CMain::GetY, "Input Data")
	.property("Accept", &CMain::GetAccept, "Are Proposals Accepted")
	.property("AccProb", &CMain::GetAccProb, "Acceptence Prob. for Proposals")	
	.property("normC", &CMain::GetnormC, "normC")
	.property("loglikelihood", &CMain::Getloglikelihood, "loglikelihood")	
  .property("sum_E_ijt", &CMain::Getsum_E_ijt, "sum_E_ijt")
	
	.property("isforcein", &CMain::Getisforcein, "isforcein")
	  
  ;
}       

CMain::CMain(arma::mat Y_) {
  Data.Y = Y_ ;	
  IterCount = 0 ;
  
  Data.isforcein = false ;
}

CMain::~CMain(){ } //Destructor

void CMain::Initialize() { Param.Initialize(Data) ; }

void CMain::check_random_generate() { Param.check_random_generate(Data) ; }
void CMain::Iterate() {
  IterCount++;
  Param.iterate(IterCount, Data, 100);
}
void CMain::Run(int iter) {
  for (int i = 1; i <= iter; i++) {
    IterCount++;
    Param.iterate(IterCount, Data, 100);
    if (IterCount%100==0){
      Rprintf("Iter: %d \n",IterCount) ; 
    }
  }
}
void CMain::SetMsgLevel(int level) { Data.msg_level = level ; }
int CMain::GetMsgLevel() { return Data.msg_level ; }
void CMain::clearE_ijt() { Param.clearE_ijt() ; }

void CMain::SetE_mat(arma::mat E_mat_) { Param.E_mat = E_mat_ ; }
arma::mat CMain::GetE_mat() { return Param.E_mat ; }
void CMain::Setmu_vec(arma::vec mu_vec_) { Param.mu_vec = mu_vec_ ; }
arma::vec CMain::Getmu_vec() { return Param.mu_vec ; }
void CMain::Setsig2_vec(arma::vec sig2_vec_) { Param.sig2_vec = sig2_vec_ ; }
arma::vec CMain::Getsig2_vec() { return Param.sig2_vec ; }
void CMain::SetBeta(arma::mat Beta_) {
  Param.Beta = Beta_;
}
arma::mat CMain::GetBeta() { return Param.Beta ; }
void CMain::SetG_mat(arma::mat G_mat_) { Param.G_mat = G_mat_ ; }
arma::mat CMain::GetG_mat() { return Param.G_mat ; }

// Hyperparameters	

void CMain::Settheta_mu(double theta_mu_) { Data.theta_mu = theta_mu_ ; }
double CMain::Gettheta_mu() { return Data.theta_mu ; }
void CMain::Settau2_mu(double tau2_mu_) { Data.tau2_mu = tau2_mu_ ; }
double CMain::Gettau2_mu() { return Data.tau2_mu ; }
void CMain::Seta_sigma(double a_sigma_) { Data.a_sigma = a_sigma_ ; }
double CMain::Geta_sigma() { return Data.a_sigma ; }
void CMain::Setb_sigma(double b_sigma_) { Data.b_sigma = b_sigma_ ; }
double CMain::Getb_sigma() { return Data.b_sigma ; }
void CMain::Settheta_alpha(double theta_alpha_) { Data.theta_alpha = theta_alpha_ ; }
double CMain::Gettheta_alpha() { return Data.theta_alpha ; }
void CMain::Settau2_alpha(double tau2_alpha_) { Data.tau2_alpha = tau2_alpha_ ; }
double CMain::Gettau2_alpha() { return Data.tau2_alpha ; }
void CMain::Setstepsize_alpha(double stepsize_alpha_) { Data.stepsize_alpha = stepsize_alpha_ ; }
double CMain::Getstepsize_alpha() { return Data.stepsize_alpha ; }
void CMain::Seta_beta(double a_beta_) { Data.a_beta = a_beta_ ; }
double CMain::Geta_beta() { return Data.a_beta ; }
void CMain::Setb_beta(double b_beta_) { Data.b_beta = b_beta_ ; }
double CMain::Getb_beta() { return Data.b_beta ; }
void CMain::Setstepsize_beta(double stepsize_beta_) { Data.stepsize_beta = stepsize_beta_ ; }
double CMain::Getstepsize_beta() { return Data.stepsize_beta ; }
void CMain::Seta_betaG(double a_betaG_) { Data.a_betaG = a_betaG_ ; }
double CMain::Geta_betaG() { return Data.a_betaG ; }
void CMain::Setb_betaG(double b_betaG_) { Data.b_betaG = b_betaG_ ; }
double CMain::Getb_betaG() { return Data.b_betaG ; }
void CMain::Setthreshold_on(double threshold_on_) { Data.threshold_on = threshold_on_ ; }
double CMain::Getthreshold_on() { return Data.threshold_on ; }
void CMain::SetE_forcein_mat(arma::mat E_forcein_mat_) { 
  Data.E_forcein_mat = E_forcein_mat_ ; 
  Data.isforcein = true ; 
}
arma::mat CMain::GetE_forcein_mat() { return Data.E_forcein_mat ; }
bool CMain::Getisforcein() { return Data.isforcein ; }
void CMain::Setstepsize_gamma(double stepsize_gamma_) { Data.stepsize_gamma = stepsize_gamma_ ; }
double CMain::Getstepsize_gamma() { return Data.stepsize_gamma ; }

void CMain::Seta_gamma(double a_gamma_) { Data.a_gamma = a_gamma_ ; }
double CMain::Geta_gamma() { return Data.a_gamma ; }
void CMain::Setb_gamma(double b_gamma_) { Data.b_gamma = b_gamma_ ; }
double CMain::Getb_gamma() { return Data.b_gamma ; }
void CMain::SetannotMat(arma::mat annotMat_) { Data.annotMat = annotMat_ ; }
arma::mat CMain::GetannotMat() { return Data.annotMat ; }
void CMain::Sethave_annot(bool have_annot_) { Data.have_annot = have_annot_ ; }
bool CMain::Gethave_annot() { return Data.have_annot ; }
void CMain::Setgamma_mat(arma::mat gamma_mat_) { Param.gamma_mat = gamma_mat_ ; }
arma::mat CMain::Getgamma_mat() { return Param.gamma_mat ; }
void CMain::Setu_mat(arma::mat u_mat_) { Param.u_mat = u_mat_ ; }
arma::mat CMain::Getu_mat() { return Param.u_mat ; }

// Print data or result

arma::mat CMain::GetY() { return Data.Y ; }
arma::vec CMain::GetAccept() { return Param.is_accept_vec ; }
arma::vec CMain::GetAccProb() { return Param.accept_prob_vec ; }
arma::vec CMain::GetnormC() { return Param.normC_vec ; }
double CMain::GetlogPost() { return Param.logPost ; }
double CMain::Getloglikelihood() { return Param.loglikelihood ; }
arma::cube CMain::Getsum_E_ijt() { return Param.sum_E_ijt ; }
