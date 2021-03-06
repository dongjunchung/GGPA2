// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// Gibbs_e_it_ann
arma::mat Gibbs_e_it_ann(arma::mat beta_mat, arma::mat G_mat, int n_SNP, int n_GibbsStep, arma::mat true_Gamma, arma::mat A_mat);
RcppExport SEXP _GGPA2_Gibbs_e_it_ann(SEXP beta_matSEXP, SEXP G_matSEXP, SEXP n_SNPSEXP, SEXP n_GibbsStepSEXP, SEXP true_GammaSEXP, SEXP A_matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type beta_mat(beta_matSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type G_mat(G_matSEXP);
    Rcpp::traits::input_parameter< int >::type n_SNP(n_SNPSEXP);
    Rcpp::traits::input_parameter< int >::type n_GibbsStep(n_GibbsStepSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type true_Gamma(true_GammaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A_mat(A_matSEXP);
    rcpp_result_gen = Rcpp::wrap(Gibbs_e_it_ann(beta_mat, G_mat, n_SNP, n_GibbsStep, true_Gamma, A_mat));
    return rcpp_result_gen;
END_RCPP
}
// Gibbs_e_it_no_ann
arma::mat Gibbs_e_it_no_ann(arma::mat beta_mat, arma::mat G_mat, int n_SNP, int n_GibbsStep);
RcppExport SEXP _GGPA2_Gibbs_e_it_no_ann(SEXP beta_matSEXP, SEXP G_matSEXP, SEXP n_SNPSEXP, SEXP n_GibbsStepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type beta_mat(beta_matSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type G_mat(G_matSEXP);
    Rcpp::traits::input_parameter< int >::type n_SNP(n_SNPSEXP);
    Rcpp::traits::input_parameter< int >::type n_GibbsStep(n_GibbsStepSEXP);
    rcpp_result_gen = Rcpp::wrap(Gibbs_e_it_no_ann(beta_mat, G_mat, n_SNP, n_GibbsStep));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP _rcpp_module_boot_cGGPA2module();

static const R_CallMethodDef CallEntries[] = {
    {"_GGPA2_Gibbs_e_it_ann", (DL_FUNC) &_GGPA2_Gibbs_e_it_ann, 6},
    {"_GGPA2_Gibbs_e_it_no_ann", (DL_FUNC) &_GGPA2_Gibbs_e_it_no_ann, 4},
    {"_rcpp_module_boot_cGGPA2module", (DL_FUNC) &_rcpp_module_boot_cGGPA2module, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_GGPA2(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
