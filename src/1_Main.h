#if !defined(_Main_H)
#define _Main_H

#include <RcppArmadillo.h>
#include "3_Param.h"

class CMain {
	
  public:
		
    CMain(arma::mat Y_);
	  ~CMain(); //destructor

    void Initialize() ; void check_random_generate() ; 
    void Iterate(); void Run(int);
    void SetMsgLevel(int); int GetMsgLevel();
    void clearE_ijt() ; 
    // void SetOptionEit(bool OptionEit_); bool GetOptionEit();
    
		// Initial values // 
    void SetE_mat(arma::mat); arma::mat GetE_mat();
		void Setmu_vec(arma::vec); arma::vec Getmu_vec();
		void Setsig2_vec(arma::vec); arma::vec Getsig2_vec();
    void SetBeta(arma::mat); arma::mat GetBeta();
    void SetG_mat(arma::mat); arma::mat GetG_mat();
		
		void SetE_forcein_mat(arma::mat); arma::mat GetE_forcein_mat();
		
		// Hyperparameters	
		void Settheta_mu(double) ; double Gettheta_mu() ;  
		void Settau2_mu(double) ; double Gettau2_mu() ;
		void Seta_sigma(double) ; double Geta_sigma() ; 
		void Setb_sigma(double) ; double Getb_sigma() ; 
		void Settheta_alpha(double) ; double Gettheta_alpha(); 
		void Settau2_alpha(double) ; double Gettau2_alpha() ;
		void Setstepsize_alpha(double) ; double Getstepsize_alpha() ;
		void Seta_beta(double) ; double Geta_beta() ;
		void Setb_beta(double) ; double Getb_beta() ;
		void Setstepsize_beta(double) ; double Getstepsize_beta() ;
		void Seta_betaG(double) ; double Geta_betaG() ;
		void Setb_betaG(double) ; double Getb_betaG() ;
		void Setthreshold_on(double) ; double Getthreshold_on() ;
		void SetPriorSetting(int) ; int GetPriorSetting() ;
		void Setpriorprob_G(arma::mat) ; arma::mat Getpriorprob_G() ;
		void Setstepsize_gamma(double) ; double Getstepsize_gamma() ;
		
		void Seta_gamma(double) ; double Geta_gamma() ; 
		void Setb_gamma(double) ; double Getb_gamma() ; 
		void SetannotMat(arma::mat) ; arma::mat GetannotMat() ;
		void Sethave_annot(bool) ; bool Gethave_annot() ;
		void Setgamma_mat(arma::mat) ; arma::mat Getgamma_mat() ;
		void Setu_mat(arma::mat) ; arma::mat Getu_mat() ;
		
		// Print data or result
		arma::mat GetY();
		arma::vec GetAccept() ; arma::vec GetAccProb();
		arma::vec GetnormC() ;
		double GetlogPost() ;
		double Getloglikelihood() ;
		arma::cube Getsum_E_ijt() ;
		
		bool Getisforcein() ;
		
  private:
    
    CData Data;
    CParam Param;
    int IterCount;

};

#endif  //_CMain_H
