#include <RcppDist.h> 
#include <RcppArmadilloExtensions/sample.h> 
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
using namespace Rcpp;
using namespace RcppArmadillo;
//------------------------------------------------------------------------------
// HELPER FUNCTIONS
//------------------------------------------------------------------------------
double LogSumExp_cpp(arma::rowvec logX){
  double a = max(logX);
  return(  a + log(accu( exp( logX-a ) )));
}

//------------------------------------------------------------------------------
// [[Rcpp::export]]
arma::colvec reverse_cumsum_cpp(arma::colvec X){
  return( accu(X) - arma::cumsum(X));
}

//------------------------------------------------------------------------------
double log_Const_prod_gamma(int D, double nu){
  double Q = 0.0;
  for( int d=1; d<(D+1); d++){
    Q += lgamma( ( nu + 1.0 - (d) ) * .5);
  }
  return(Q);
}

//----------------------------------------------------------------------------
// [[Rcpp::export]]
double Const_sum_digamma(int D, double nu){
  double Q = 0.0;
  for( int d=1; d<(D+1); d++){
    Q += R::digamma( (nu + 1.0 - (d)) * .5);
  }
  return(Q);
}

//------------------------------------------------------------------------- B.79
double logB_wish(arma::mat W, double nu, double log_Const_prod_gamma){
  
  int D = W.n_cols;
  double p1 = - 0.5 * nu * log(arma::det(W));
  double p2 = (nu * D * .5) * log(2.0) + ( D * (D-1.0) / 4.0 ) * log(arma::datum::pi);
  return( p1 - p2 - log_Const_prod_gamma );
}

//------------------------------------------------------------------------- B.81
// [[Rcpp::export]]
double ElogDetLambda(arma::mat W, double Const_sum_digamma){
  
  int D = W.n_cols;
  return( Const_sum_digamma + log(arma::det(W)) + D*log(2.0) );
  
}

//------------------------------------------------------------------------- B.82
double  H_Lambda(arma::mat W, double nu, 
                 double Const_sum_digamma, 
                 double log_Const_prod_gamma){
  int D = W.n_cols;
  double ElDL = ElogDetLambda(W, Const_sum_digamma);
  double lnB = logB_wish(W, nu, log_Const_prod_gamma);
  return( - lnB - ( nu - D - 1.0) * 0.5 * ElDL + nu * D * 0.5);
  
} 
// [[Rcpp::export]]
arma::colvec E_log_beta(arma::colvec a,
                        arma::colvec b){
  
  int n = a.n_elem;
  arma::colvec res(n);
  for(int i=0; i<n; i++){
    res[i] = R::digamma(a[i]) - R::digamma(a[i] + b[i]);
  }
  return(res);
}

// -----------------------------------------------------------------------------
// [[Rcpp::export]]
arma::colvec E_log_p_Y_Mtheta_cpp_mvt(arma::mat Y,
                                      arma::colvec ml,
                                      double tl,
                                      double cl,
                                      arma::mat Dl){
  
  int p = Y.n_cols;
  int N = Y.n_rows;
  
  double CSDg = Const_sum_digamma(p, cl);
  double ell1 = ElogDetLambda(Dl, CSDg);                   // \ell1 in algorithm
  
  arma::colvec fq(N);
  arma::mat DIFF = Y - arma::repelem(ml.t(), N, 1); // Nj x D
  arma::mat DIFFt = DIFF.t();
  for(int i=0; i<N; i++){
    arma::vec temp = (DIFF.row(i) * Dl * (DIFFt.col(i)));
    fq(i) = temp[0];
  }
  
  arma::colvec ell2 = ( - p * 1.0/(tl) - cl * fq ); // Nj x 1
  
  return(.5 * (ell1 + ell2));
}

// -----------------------------------------------------------------------------
// [[Rcpp::export]]
arma::colvec E_log_p_X_cpp_mvt(arma::mat X,
                               arma::colvec ml,
                               double tl,
                               double cl,
                               arma::mat Dl){
  
  int q = X.n_cols;
  int J = X.n_rows;
  
  double CSDg = Const_sum_digamma(q, cl);
  double ell1 = ElogDetLambda(Dl, CSDg);                   // \ell1 in algorithm
  
  arma::colvec fq(J);
  arma::mat DIFF = X - arma::repelem(ml.t(), J, 1); // J x q
  arma::mat DIFFt = DIFF.t();
  for(int i=0; i<J; i++){
    arma::vec temp = (DIFF.row(i) * Dl * (DIFFt.col(i)));
    fq(i) = temp[0];
  }
  
  arma::colvec ell2 = ( - q * 1.0/(tl) - cl * fq ); // J x 1
  
  return(.5 * (ell1 + ell2));
}

//------------------------------------------------------------------------------
// MAIN FUNCTIONS
//------------------------------------------------------------------------------
// --------------------------------------------------------------------------  1
// [[Rcpp::export]]
arma::mat Update_Vk_DP_conc_cpp(double b_bar,
                                arma::mat RHO_jk){
  //b_bar is s1/s2, which needs to be updated every time s1 and s2 are updated
  int K = RHO_jk.n_cols;

  arma::colvec mk = (arma::sum(RHO_jk,0)).t(); // This calculates \sum_{j=1}^{J}RHO_jk for k = 1,...,H. mk is Hx1 vector
  arma::colvec rev_cs_mk = reverse_cumsum_cpp(mk); // This calculates \sum_{l=k+1}^{H}\sum_{j=1}^{J}RHO_jl, 

  arma::colvec a_tilde_vk      = mk        + 1.0; //This is \bar{a}_k = 1 + \sum_{j=1}^{J}RHO_jk
  a_tilde_vk[K-1] = 1.0; //Since \bar{a}_k update goes from 1,...,H, We fix \bar{a}_H = 1
  arma::colvec b_tilde_vk      = rev_cs_mk + b_bar; //This is \bar{b}_k = s1/s2 + \sum_{l=k+1}\sum_{j=1}^{J}RHO_jl
  b_tilde_vk[K-1] = 1e-10; //Since \bar{b}_k update goes from 1,...,H, We fix \bar{b}_H to a ver small number

  arma::colvec E_ln_Vk    = E_log_beta(a_tilde_vk, b_tilde_vk); //This gives g(\bar{a}_k, \bar{b}_k) for k = 1,...,H as an Hx1 vector
  arma::colvec E_ln_1mVk  = E_log_beta(b_tilde_vk, a_tilde_vk); //This gives g(\bar{b}_k, \bar{a}_k) for k = 1,...,H as an Hx1 vector
  arma::colvec sE_ln_1mVk = shift(E_ln_1mVk, +1); // This shifts g(\bar{b}_k, \bar{a}_k) for k = 1,...,H by 1 cyclically, which yields the vector (g(\bar{b}_H, \bar{a}_H), g(\bar{b}_1, \bar{a}_1), ... , g(\bar{b}_H-1, \bar{a}_H-1)) 

  sE_ln_1mVk[0] = 0; // Set g(\bar{b}_H, \bar{a}_H) = 0 in the vector (g(\bar{b}_H, \bar{a}_H), g(\bar{b}_1, \bar{a}_1), ... , g(\bar{b}_H-1, \bar{a}_H-1)) 
  // This is because, we never really consider  g(\bar{b}_H, \bar{a}_H) when calculating \sum{r=1}^{k-1} g(\bar{b}_r, \bar{a}_r) for k = 1,...,H. when r = H, we consider g(\bar{b}_H-1, \bar{a}_H-1)
  arma::colvec CS_E_ln_1mVk = arma::cumsum(sE_ln_1mVk); // This gives\sum{r=1}^{k-1} g(\bar{b}_r, \bar{a}_r) for k = 1,...,H as an Hx1 vector
  arma::mat results(K,3);

  results.col(0) = a_tilde_vk;
  results.col(1) = b_tilde_vk;
  results.col(2) = E_ln_Vk + CS_E_ln_1mVk ; // This gives g(\bar{a}_k, \bar{b}_k) + \sum{r=1}^{k-1} g(\bar{b}_r, \bar{a}_r) for k = 1,...,H as an Hx1 vector, which is needed for ELBO calculation

  return(results);
}
// --------------------------------------------------------------------------- 2
// [[Rcpp::export]]
arma::cube Update_Ulk_cpp(arma::field<arma::mat> XI_jil,
                          arma::mat RHO_jk,
                          double const b_bar,
                          int const L,
                          int const J,
                          int const H){
  
  //b_bar is r1/r2, which needs to be updated every time r1 and r2 are updated
  arma::mat N_jl(J,L);
  
  for(int j=0; j<J; j++){
    N_jl.row(j) = arma::sum(XI_jil(j),0); //This calculates \sum_{i=1}^{n_j} XI_jil. So N_jl is a JxL matrix
  }
  
  arma::mat Q_lk = N_jl.t() * RHO_jk ; // This calculates \sum_{j=1}^{J}\sum_{i=1}^{n_j} XI_jil RHO_jk. So this gives a LxH matrix
  arma::mat a_bar_Ulk(L,H);
  arma::mat b_bar_Ulk(L,H);
  arma::mat E_lnOmega_lk(L,H);
  
  
  for(int k = 0; k<H; k++){
    
    arma::colvec rQl_k  = reverse_cumsum_cpp(Q_lk.col(k)); // This calculates \sum_{t=l+1}^{L}\sum_{j=1}^{J}\sum_{i=1}^{n_j} XI_jit RHO_jk
    arma::colvec G1 = 1.0 + Q_lk.col(k);  //This is \bar{a}_lk = 1 + \sum_{j=1}^{J}\sum_{i=1}^{n_j} XI_jil RHO_jk, for fixed k
    arma::colvec G2 = b_bar + rQl_k;      //This is \bar{b}_lk = r1/r2 + \sum_{t=l+1}^{L}\sum_{j=1}^{J}\sum_{i=1}^{n_j} XI_jit RHO_jk, for fixed k
    
    G1[L-1] = 1;  //Since \bar{a}_lk update goes from 1,...,L-1, We fix \bar{a}_{L-1 k} = 1 for all k
    G2[L-1] = 1e-10; //Since \bar{b}_lk update goes from 1,...,L-1, We fix \bar{b}_{L-1 k} to a very small number for all k
    
    a_bar_Ulk.col(k) =  G1;
    b_bar_Ulk.col(k) =  G2;
    
    arma::colvec E_ln_Ul_k    = E_log_beta(G1, G2); //This gives g(\bar{a}_lk, \bar{b}_lk) for l = 1,...,L as an Lx1 vector for each k = 1,...,H
    arma::colvec E_ln_1mUl_k  = E_log_beta(G2, G1); //This gives g(\bar{b}_lk, \bar{a}_lk) for l = 1,...,L as an Lx1 vector for each k = 1,...,H
    arma::colvec sE_ln_1mUl_k = shift(E_ln_1mUl_k, +1);  //This shifts the vector(g(\bar{b}_lk, \bar{a}_lk)), l=1,...,L by 1 position cyclically
    sE_ln_1mUl_k[0] = 0; //This sets (g(\bar{b}_Lk, \bar{a}_Lk)) = 0
    
    arma::colvec CS_E_ln_1mUL_k = arma::cumsum(sE_ln_1mUl_k); // This gives \sum{r=1}^{l-1} g(\bar{b}_rk, \bar{a}_rk)
    E_lnOmega_lk.col(k) = CS_E_ln_1mUL_k + E_ln_Ul_k; // This gives g(\bar{a}_lk, \bar{b}_lk) + \sum{r=1}^{l-1} g(\bar{b}_rk, \bar{a}_rk) for l = 1,...,L as an Lx1 vector, which is needed for ELBO calculation
  }//Running over 
  
  
  arma::cube results(L,H,3);
  results.slice(0) = a_bar_Ulk;
  results.slice(1) = b_bar_Ulk;
  results.slice(2) = E_lnOmega_lk;
  
  return(results);
}

// --------------------------------------------------------------------------- 3

// [[Rcpp::export]]
arma::colvec Update_alpha_concentration_par(arma::colvec a_tilde_Vk,
                                              arma::colvec b_tilde_Vk,
                                              arma::colvec conc_hyper){

  int H = b_tilde_Vk.n_rows;
  arma::colvec upd_par(2);
  a_tilde_Vk.shed_row(H-1); // Removing the Hth row of a_tilde_Vk
  b_tilde_Vk.shed_row(H-1); // Removing the Hth row of b_tilde_Vk
  //This is needed because we need \sum_{k=1}^{H-1}g(\bar{b}_k, \bar(a)_k)

  upd_par[0] = conc_hyper[0] + H - 1.0 ; //This is s1 = a_{\alpha} + (H-1)
  upd_par[1] = conc_hyper[1] - arma::accu(E_log_beta(b_tilde_Vk,a_tilde_Vk)); //This is s2 = b_{\alpha} - \sum_{k=1}^{H-1}g(\bar{b}_k, \bar(a)_k)

  return(upd_par);
}

// --------------------------------------------------------------------------- 4
// [[Rcpp::export]]
arma::colvec Update_beta_concentration_par(arma::mat a_bar_Ulk,
                                           arma::mat b_bar_Ulk,
                                           arma::colvec conc_hyper,
                                           int L,
                                           int H){
  
  arma::colvec upd_par(2);
  a_bar_Ulk.shed_row(L-1); // Removing the Lth row of a_bar_Ulk. a_bar_Ulk is an LxH matrix originally
  b_bar_Ulk.shed_row(L-1); // Removing the Lth row of b_bar_Ulk. b_bar_Ulk is an LxH matrix originally
  
  arma::colvec R(H);
  
  for(int k = 0; k < H; k ++){
    
    R[k] =  arma::accu( E_log_beta(b_bar_Ulk.col(k),a_bar_Ulk.col(k)) ); // This calculates \sum_{k=1}^{H}g(\bar{b}_{lk}, \bar(a)_{lk})
    
  }
  

  upd_par[0] = conc_hyper[0] + H * (L - 1.0); //This is r1 = a_{\beta} + H(L-1)
  upd_par[1] = conc_hyper[1] - arma::accu(R); //This is r2 = b_{\beta} - \sum_{l=1}^{L-1}\sum_{k=1}^{H}g(\bar{b}_{lk}, \bar(a)_{lk})
  //This is why we removed the row L-1 a_bar_Ulk and b_bar_Ulk
  
  
  return(upd_par);
}

// --------------------------------------------------------------------------- 5
// [[Rcpp::export]]
Rcpp::List Update_THETAl_Y_cpp_mvt(arma::field<arma::mat> Y_grouped,
                                 arma::field<arma::mat> XI_jil,
                                 arma::colvec m0,
                                 double lambda0,
                                 double nu0,
                                 arma::mat W0,
                                 arma::mat iW0){

  int D = W0.n_rows;
  int J = XI_jil.n_elem;
  int L = XI_jil(0).n_cols;


  arma::rowvec Nl(L, arma::fill::zeros);
  arma::mat SumY_l(D, L, arma::fill::zeros);
  arma::mat Ybarl(D, L, arma::fill::zeros);
  arma::cube Sl(D,D, L, arma::fill::zeros);
  arma::cube Sl_over_Nl(D,D, L, arma::fill::zeros);

  arma::mat ml(D,L);
  arma::rowvec lambdal(L);
  arma::rowvec nul(L);
  arma::cube Wl(D,D,L);


  for(int j=0; j <J; j++){
    // For a fixed j, XI_jil is an n_j x L matrix
    arma::mat subY = Y_grouped[j];       // This gives Y_j, which is a n_j x p matrix
    Nl      += arma::sum( XI_jil(j), 0); // This calculates [\sum_{i=1}^{n_j} XI_jil] for a fixed j for l = 1,...,L, and running the loop over j gives \sum_{j=1}^{J}\sum_{i=1}^{n_j} XI_jil. So Nl is 1xL vector.
    SumY_l  += subY.t() * XI_jil(j);  // p x L // This gives Y_j^T  XI_jil, which is \sum_{j=1}^{J}\sum_{i=1}^{n_j} XI_jil y_{ji}, which is p x L
  }

  arma::mat B0M0  = arma::repelem(lambda0 * m0,1,L); // p x L. B0M0 is repeating \lambda_0\mu_0, which is a px1 vector L times along a column to get a pxL matrix. This is to facilitate matrix addition 
  lambdal         = Nl + lambda0; // This is t_l = N_l + \lambda0
  nul             = Nl + nu0 ;    // This is c_l = N_l + \nu_0
  arma::mat iBL   = arma::repelem(1.0 / lambdal, D, 1); // p x L. iBL is repeating 1/t_l, which is 1xL vector p times along the rows to get a pxL matrix
  ml              = ( B0M0 + SumY_l ) % ( iBL ); // Element-wise multiplying ( B0M0 + SumY_l ) with iBL gives m_l


  for(int l=0; l<L; l++){

    if(Nl(l)>0){
      Ybarl.col(l) = SumY_l.col(l)/Nl(l); // This calculates \bar{y}_l
    }

    for(int j=0; j <J; j++){
      arma::mat subY = Y_grouped[j].t(); // p x Nj
      for(int ii=0; ii<subY.n_cols; ii++){
        arma::mat temp = ( subY.col(ii) - Ybarl.col(l)) *
          (subY.col(ii) - Ybarl.col(l)).t(); //This calculating (y_ji - \bar{y}_l)(y_ji - \bar{y}_l)^T for a fixed j, i, and l
        double xi = (XI_jil(j))(ii,l); // This is XI_{jil} for a fixed j, i, and l
        Sl.slice(l) += xi * temp; //This calculating XI_{jil}(y_ji - \bar{y}_l)(y_ji - \bar{y}_l)^T for a fixed j, i, and l
      }// Finally this sums over i, for a fixed j, and fixed l
    }//This sums over j, for a fixed l. Returning S_l as Sl.slice(l)
    //Rcpp::Rcout<< Sl_Nl.slice(l) << ".....\nNl = " << Nl(l) << "\n.....\n";

    if(Nl(l)>0){
      Sl_over_Nl.slice(l) = Sl.slice(l)/Nl(l); // This calculates S_l/N_l, to be used in the calculation of the ELBO
    }

    Wl.slice(l) = arma::inv_sympd(
      (iW0)+
        Sl.slice(l) +
        lambda0*Nl(l)/(lambda0+Nl(l)) *
        (Ybarl.col(l)-m0) * (Ybarl.col(l)-m0).t()
    ); // This gives D_l^{-1}

  }// This completes the loop over l
  //Rcpp::Rcout << "Sl = " << Sl;

  Rcpp::List results = Rcpp::List::create(
    Rcpp::_["ml"]  = ml,
    Rcpp::_["lambdal"] = lambdal.t(),
    Rcpp::_["nul"] = nul.t(),
    Rcpp::_["Wl"] = Wl,
    Rcpp::_["Ybar"] = Ybarl,
    Rcpp::_["Nl"] = Nl,
    Rcpp::_["Sl_over_Nl"] = Sl_over_Nl

  );
  return(results);

}

// --------------------------------------------------------------------------- 6
// [[Rcpp::export]]
Rcpp::List Update_THETAl_X_cpp_mvt(arma::mat X,
                                   arma::mat RHO_jk,
                                   arma::colvec m0,
                                   double lambda0,
                                   double nu0,
                                   arma::mat W0,
                                   arma::mat iW0){
  
  int D = W0.n_rows;
  int J = RHO_jk.n_rows;
  int H = RHO_jk.n_cols;
  
  
  arma::rowvec Nk(H, arma::fill::zeros);
  arma::mat SumX_k(D, H, arma::fill::zeros);
  arma::mat Xbark(D, H, arma::fill::zeros);
  arma::cube Sk(D,D, H, arma::fill::zeros);
  arma::cube Sk_over_Nk(D,D, H, arma::fill::zeros);
  
  arma::mat mk(D,H);
  arma::rowvec lambdak(H);
  arma::rowvec nuk(H);
  arma::cube Wk(D,D,H);
  
 // RHO_jk is  J x H matrix
    Nk  = arma::sum(RHO_jk, 0); // This calculates \sum_{j=1}^{J}RHO_jk  for k = 1,...,H. So Nk is 1xH vector.
    SumX_k  = X.t() * RHO_jk;  // q x H. X is a J x q matrix. This gives \sum_{j=1}^{J}RHO_jk x_j
  
  
  arma::mat B0M0  = arma::repelem(lambda0 * m0,1,H); // q x H. B0M0 is repeating \lambda_0\mu_0, which is a qx1 vector H times along a column to get a qxH matrix. This is to facilitate matrix addition
  lambdak         = Nk + lambda0; // This is t_k = N_k + \lambda0
  nuk             = Nk + nu0 ;    // This is c_k = N_k + \nu_0
  arma::mat iBL   = arma::repelem(1.0 / lambdak, D, 1); // q x H. iBL is repeating 1/t_k, which is 1xH vector q times along the rows to get a qxH matrix
  mk              = ( B0M0 + SumX_k ) % ( iBL ); // Element-wise multiplying ( B0M0 + SumX_k ) with iBL gives m_k
  
  
  for(int k=0; k<H; k++){
    
    if(Nk(k)>0){
      Xbark.col(k) = SumX_k.col(k)/Nk(k);  // This calculates \bar{x}_k
    }
    
    
      arma::mat subX = X.t(); // D x J
      for(int j=0; j<subX.n_cols; j++){
        arma::mat temp = ( subX.col(j) - Xbark.col(k)) *
          (subX.col(j) - Xbark.col(k)).t(); //This calculating (x_j - \bar{x}_k)(x_j - \bar{x}_k)^T for a fixed j and k
        double xi = RHO_jk(j,k);  //This is RHO_jk for a fixed j and k
        Sk.slice(k) += xi * temp; //This calculating RHO_jk(x_j - \bar{x}_k)(x_j - \bar{x}_k)^T for a fixed j and k
      }// Finally this sums over j, for a fixed k
    
    
    //Rcpp::Rcout<< Sl_Nl.slice(l) << ".....\nNl = " << Nl(l) << "\n.....\n";
    
    if(Nk(k)>0){
      Sk_over_Nk.slice(k) = Sk.slice(k)/Nk(k);  // This calculates S_k/N_k, to be used in the calculation of the ELBO
    }
    
    Wk.slice(k) = arma::inv_sympd(
      (iW0)+
        Sk.slice(k) +
        lambda0*Nk(k)/(lambda0+Nk(k)) *
        (Xbark.col(k)-m0) * (Xbark.col(k)-m0).t()
    ); // This gives D_k^{-1}
    
  }
  //Rcpp::Rcout << "Sl = " << Sl;
  
  Rcpp::List results = Rcpp::List::create(
    Rcpp::_["mk"]  = mk,
    Rcpp::_["lambdak"] = lambdak.t(),
    Rcpp::_["nuk"] = nuk.t(),
    Rcpp::_["Wk"] = Wk,
    Rcpp::_["Xbar"] = Xbark,
    Rcpp::_["Nk"] = Nk,
    Rcpp::_["Sk_over_Nk"] = Sk_over_Nk
  
  );
  return(results);
  
}

// --------------------------------------------------------------------------- 7
// [[Rcpp::export]]
arma::field<arma::mat> Update_XIjil_cpp(arma::field<arma::mat> Y_grouped,
                                            arma::mat RHO_jk,
                                            arma::mat ElnOM_lk,

                                            arma::colvec Nj,

                                            arma::mat ml,
                                            arma::colvec tl,
                                            arma::colvec cl,
                                            arma::cube Dl,

                                            int const L,
                                            int const J,
                                            int const K){


  arma::field<arma::mat> XI_jil(J); // J different nj x L matrices
  // ElnOM_lk contains the matrix [g(\bar{a}_lh, \bar{b}_lh) + \sum_{t=1}^{l-1}g(\bar{b}_th, \bar{a}_th)]
  // For l = 1,...,L; h = 1,...,H
  // So ElnOM_lk is an LxH matrix
  // RHO_jk is a JxH matrix 
  arma::mat M_LJ  =  ElnOM_lk * RHO_jk.t(); // This gives \sum_{h=1}^{H} RHO_jh g(\bar{a}_lh, \bar{b}_lh) + \sum_{t=1}^{l-1}g(\bar{b}_th, \bar{a}_th)
  // M_LJ is an LxJ matrix
  arma::mat M_JL  = M_LJ.t(); // Transposing this gives M_JL is a JxL matrix

  arma::rowvec logunn(L);

  for(int j= 0; j<J; j++){
    int Nj = Y_grouped[j].n_rows;
    arma::mat TT(Nj,L);

    for(int l= 0; l<L; l++){
    // TT is an n_j x L matrix
      TT.col(l) = E_log_p_Y_Mtheta_cpp_mvt(Y_grouped[j],
             ml.col(l), tl(l), cl(l), Dl.slice(l));

    }
    // This calculates 0.5[\sum_{d=1}^{p}\Psi(c_l^y - d + 1)/2 + plog2 + log|D_l^y| -
    //  p/t_l^y + c_l^y (y_ji - m_l^y)^T D_l^y (y_ji - m_l^y)]; for i = 1,...,n_j; l=1,..,L
    arma::mat temp1 = TT;  // Nj x L
    arma::mat tempres(Nj,L);

    for(int i=0; i<Nj; i++){
    // Adding jth row of M_JL and ith row of temp1 
      logunn = M_JL.row(j) + temp1.row(i);
    // This gives an n_j x L matrix
      tempres.row(i) =  logunn - LogSumExp_cpp(logunn); // Normalizing by LogSumExp trick
    }
    // For a fixed j, returning the exp of normalized XI_jil to gives the probabilities XI_jil
    XI_jil(j) = exp(tempres);

  }// Running over the loop of j
  return(XI_jil); // This is essentially a list of J probabilities with jth one being of dimension n_jxL
}





// -----------------------------------------------------------------------------
// FOR VI WITH REGRESSION OF Y ON X
// -----------------------------------------------------------------------------
// HELPER FUNCTION
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------1
// [[Rcpp::export]]
arma::colvec E_log_p_Y_X_Mtheta_cpp_mvt(arma::mat Y,
                                        arma::mat Ml,
                                        arma::colvec X,
                                        arma::mat Ul,
                                        double cl,
                                        arma::mat Dl){
  
  int p = Y.n_cols;
  int q = X.size();
  int N = Y.n_rows;
  
  // Ml is qxp dimensional 
  double CSDg = Const_sum_digamma(p, cl);
  double ell1 = ElogDetLambda(Dl, CSDg);                   // \ell1 in algorithm
  
  arma::colvec Ml_t_X = Ml.t() * X;
  
  arma::colvec fq(N);
  arma::mat DIFF = Y - arma::repelem(Ml_t_X.t(), N, 1); // Nj x D
  arma::mat DIFFt = DIFF.t();
  for(int i=0; i<N; i++){
    arma::vec temp = (DIFF.row(i) * Dl * (DIFFt.col(i)));
    fq(i) = temp[0];
  }
  double tl = arma::trace(X.t() * Ul * X);
  arma::colvec ell2 = ( - p * (tl) - cl * fq ); // Nj x 1
  
  return(.5 * (ell1 + ell2));
}

// -----------------------------------------------------------------------------
// MAIN FUNCTIONS
// --------------------------------------------------------------------------- 1
// [[Rcpp::export]]
arma::field<arma::mat> Update_XIjil_X_cpp(arma::field<arma::mat> Y_grouped,
                                          arma::mat RHO_jk,
                                          arma::mat ElnOM_lk,
                                          
                                          arma::colvec Nj,
                                          
                                          arma::mat X,
                                          arma::cube Ml,
                                          arma::cube Ul,
                                          arma::colvec cl,
                                          arma::cube Dl,
                                          
                                          int const L,
                                          int const J,
                                          int const K){
  
  
  arma::field<arma::mat> XI_jil(J); // J different nj x L matrices
  // ElnOM_lk contains the matrix [g(\bar{a}_lh, \bar{b}_lh) + \sum_{t=1}^{l-1}g(\bar{b}_th, \bar{a}_th)]
  // For l = 1,...,L; h = 1,...,H
  // So ElnOM_lk is an LxH matrix
  // RHO_jk is a JxH matrix 
  arma::mat M_LJ  =  ElnOM_lk * RHO_jk.t(); // This gives \sum_{h=1}^{H} RHO_jh g(\bar{a}_lh, \bar{b}_lh) + \sum_{t=1}^{l-1}g(\bar{b}_th, \bar{a}_th)
  // M_LJ is an LxJ matrix
  arma::mat M_JL  = M_LJ.t(); // Transposing this gives M_JL is a JxL matrix
  
  arma::rowvec logunn(L);
  
  for(int j= 0; j<J; j++){
    int N_j = Y_grouped[j].n_rows;
    // arma::colvec TT(N_j);
    arma::mat TT(N_j,L);
    
    for(int l= 0; l<L; l++){
      // TT is an n_j x L matrix
      TT.col(l) = E_log_p_Y_X_Mtheta_cpp_mvt(Y_grouped[j],
             Ml.slice(l),
             X.row(j).t(),
             Ul.slice(l),
             cl(l),
             Dl.slice(l));
      
    }
    // This calculates 0.5[\sum_{d=1}^{p}\Psi(c_l^y - d + 1)/2 + plog2 + log|D_l^y| -
    //  p/t_l^y + c_l^y (y_ji - m_l^y)^T D_l^y (y_ji - m_l^y)]; for i = 1,...,n_j; l=1,..,L
    arma::mat temp1 = TT;  // Nj x L
    arma::mat tempres(N_j,L);
    
    for(int i=0; i<N_j; i++){
      // Adding jth row of M_JL and ith row of temp1
      logunn = M_JL.row(j) + temp1.row(i);
      // This gives an n_j x L matrix
      tempres.row(i) =  logunn - LogSumExp_cpp(logunn); // Normalizing by LogSumExp trick
    }
    // // For a fixed j, returning the exp of normalized XI_jil to gives the probabilities XI_jil
    XI_jil(j) = exp(tempres);
    
  }// Running over the loop of j
  return(XI_jil); // This is essentially a list of J probabilities with jth one being of dimension n_jxL
}

// --------------------------------------------------------------------------- 2
//THIS IS FOR SCALAR VALUES OF TAU0 AND GAMMA0
// // [[Rcpp::export]]
// Rcpp::List Update_THETAl_Y_X_cpp_mvt(arma::field<arma::mat> Y_grouped,
//                                      arma::field<arma::mat> XI_jil,
//                                      arma::mat X,
//                                      arma::mat B0,
//                                      arma::mat iG0,
//                                      double tau0,
//                                      double gamma0,
//                                      double nu0,
//                                      arma::mat iW0){
//   
//   int p = iW0.n_rows;
//   int q = X.n_cols;
//   int J = XI_jil.n_elem;
//   int L = XI_jil(0).n_cols;
//   
//   
//   arma::rowvec Nl(L, arma::fill::zeros);
//   arma::mat SumY_l(p, L, arma::fill::zeros);
//   arma::cube SumY_j_l(p, L, J, arma::fill::zeros);
//   arma::cube Ybar_jl(p, L, J, arma::fill::zeros);
//   arma::mat Ybarl(p, L, arma::fill::zeros);
//   
//   arma::cube Hl(q,q, L, arma::fill::zeros);
//   arma::cube Ul(q,q, L, arma::fill::zeros);
//   arma::cube Ql(q,p, L, arma::fill::zeros);
//   
//   arma::cube Sl(p,p, L, arma::fill::zeros);
//   arma::field<arma::cube> S_jl(J);
//   S_jl.fill(arma::cube(p, p, L, arma::fill::zeros));
//   arma::field<arma::cube> S_jl_over_N_jl(J);
//   S_jl_over_N_jl.fill(arma::cube(p, p, L, arma::fill::zeros));
//   arma::cube Z1l(p,p, L, arma::fill::zeros);
//   arma::cube Zl(p,p, L, arma::fill::zeros);
//   
//   arma::cube Wl(p,p, L, arma::fill::zeros);
//   
//   arma::cube Ml(q, p,L);
//   arma::rowvec lambdal(L);
//   arma::rowvec cl(L);
//   arma::cube Dl(p,p,L);
//   arma::mat Nl_j(J, L, arma::fill::zeros);
//   
//   
//   for(int j=0; j <J; j++){
//     // For a fixed j, XI_jil is an n_j x L matrix
//     arma::mat subY = Y_grouped[j];       // This gives Y_j, which is a n_j x p matrix
//     Nl      += arma::sum( XI_jil(j), 0); // This calculates [\sum_{i=1}^{n_j} XI_jil] for a fixed j for l = 1,...,L, and running the loop over j gives \sum_{j=1}^{J}\sum_{i=1}^{n_j} XI_jil. So Nl is 1xL vector.
//     Nl_j.row(j) = arma::sum( XI_jil(j), 0);
//     SumY_j_l.slice(j) = subY.t() * XI_jil(j);
//     SumY_l  += subY.t() * XI_jil(j);  // p x L // This gives Y_j^T  XI_jil, which is \sum_{j=1}^{J}\sum_{i=1}^{n_j} XI_jil y_{ji}, which is p x L
//   }
//   
//   cl             = Nl + nu0 ;    // This is c_l = N_l + \nu_0
//   
//   
//   for(int l=0; l<L; l++){
//     
//     if(Nl(l)>0){
//       Ybarl.col(l) = SumY_l.col(l)/Nl(l); // This calculates \bar{y}_l
//     }
//   }
//   for(int j=0; j <J; j++){
//     for(int l=0; l<L; l++){
//       if(Nl_j(j, l)>0){
//         Ybar_jl.slice(j).col(l) = SumY_j_l.slice(j).col(l)/Nl_j(j, l); // This calculates \bar{y}_l
//       }
//     }
//   }
//   
//   for(int l=0; l<L; l++){
//     for(int j=0; j <J; j++){
//       Hl.slice(l) +=Nl_j(j,l) * X.row(j).t() * X.row(j);
//       Ql.slice(l) +=Nl_j(j,l) * X.row(j).t() * Ybar_jl.slice(j).col(l).t();
//       Z1l.slice(l) +=Nl_j(j,l) * Ybar_jl.slice(j).col(l) * Ybar_jl.slice(j).col(l).t();
//     }
//     Hl.slice(l) = Hl.slice(l) + gamma0 * tau0 *iG0;
//     Ql.slice(l) = Ql.slice(l)  + gamma0 * tau0 * (iG0 * B0);
//     Ml.slice(l) = solve(Hl.slice(l), Ql.slice(l));
//     Ul.slice(l) = Hl.slice(l).i();
//     Zl.slice(l) = Z1l.slice(l) + (gamma0 * tau0 * B0.t() * iG0 * B0) - (Ql.slice(l).t() * solve(Hl.slice(l), Ql.slice(l)));
//   }
//   
//   for(int l = 0; l<L; l++){
//     for(int j=0; j <J; j++){
//       arma::mat subY = Y_grouped[j].t(); // p x Nj
//       arma::mat temp_j(p, p, arma::fill::zeros);
//       for(int ii=0; ii<subY.n_cols; ii++){
//         arma::mat temp = ( subY.col(ii) - Ybar_jl.slice(j).col(l)) *
//           (subY.col(ii) - Ybar_jl.slice(j).col(l)).t(); //This calculating (y_ji - \bar{y}_l)(y_ji - \bar{y}_l)^T for a fixed j, i, and l
//         double xi = (XI_jil(j))(ii,l); // This is XI_{jil} for a fixed j, i, and l
//         Sl.slice(l) += xi * temp; //This calculating XI_{jil}(y_ji - \bar{y}_l)(y_ji - \bar{y}_l)^T for a fixed j, i, and l
//         temp_j += xi * temp;
//       }// Finally this sums over i, for a fixed j, and fixed l
//       S_jl(j).slice(l) = temp_j;
//     }
//     // Wl.slice(l) = arma::inv_sympd(
//     //     (iW0)+
//     //       Sl.slice(l) + Zl.slice(l)
//     //   ); // This gives D_l^{-1}
//     Wl.slice(l) = arma::inv_sympd(
//       (iW0)+
//         arma::symmatu(Sl.slice(l)) + arma::symmatu(Zl.slice(l))
//     ); // This gives D_l^{-1}
//   }
//   //This sums over j, for a fixed l. Returning S_l as Sl.slice(l)
//   // //Rcpp::Rcout<< Sl_Nl.slice(l) << ".....\nNl = " << Nl(l) << "\n.....\n";
//   for(int j = 0; j<J; j++){
//     for(int l = 0; l<L;l++){
//       if(Nl_j(j,l) >0){
//         S_jl_over_N_jl(j).slice(l) = S_jl(j).slice(l)/Nl_j(j,l);
//       }
//     }
//   }
//   // if(Nl(l)>0){
//   //   Sl_over_Nl.slice(l) = Sl.slice(l)/Nl(l); // This calculates S_l/N_l, to be used in the calculation of the ELBO
//   // }
//   // 
//   // Wl.slice(l) = arma::inv_sympd(
//   //   (iW0)+
//   //     Sl.slice(l) +
//   //     lambda0*Nl(l)/(lambda0+Nl(l)) *
//   //     (Ybarl.col(l)-m0) * (Ybarl.col(l)-m0).t()
//   // ); // This gives D_l^{-1}
//   
//   // }// This completes the loop over l
//   //Rcpp::Rcout << "Sl = " << Sl;
//   
//   Rcpp::List results = Rcpp::List::create(
//     Rcpp::_["Ybar_jl"] = Ybar_jl,
//     Rcpp::_["Zl"] = Zl,
//     Rcpp::_["Sl"] = Sl,
//     Rcpp::_["Wl"] = Wl,
//     Rcpp::_["Ml"]  = Ml,
//     Rcpp::_["Ul"] = Ul,
//     Rcpp::_["cl"] = cl.t(),
//     Rcpp::_["Hl"] = Hl,
//     Rcpp::_["Ql"] = Ql,
//     Rcpp::_["Ybar"] = Ybarl,
//     Rcpp::_["Nl"] = Nl,
//     Rcpp::_["Nl_j"] = Nl_j,
//     Rcpp::_["S_jl_over_N_jl"] = S_jl_over_N_jl
//   );
//   return(results);
//   
// }
// --------------------------------------------------------------------------- 2
//THIS IS FOR GENERAL VECTOR OF TAU AND GAMMA
// [[Rcpp::export]]
Rcpp::List Update_THETAl_Y_X_cpp_mvt(arma::field<arma::mat> Y_grouped,
                                     arma::field<arma::mat> XI_jil,
                                     arma::mat X,
                                     arma::mat B0,
                                     arma::mat iG0,
                                     arma::vec taul,
                                     arma::vec gammal,
                                     double nu0,
                                     arma::mat iW0){
  
  int p = iW0.n_rows;
  int q = X.n_cols;
  int J = XI_jil.n_elem;
  int L = XI_jil(0).n_cols;
  
  
  arma::rowvec Nl(L, arma::fill::zeros);
  arma::mat SumY_l(p, L, arma::fill::zeros);
  arma::cube SumY_j_l(p, L, J, arma::fill::zeros);
  arma::cube Ybar_jl(p, L, J, arma::fill::zeros);
  arma::mat Ybarl(p, L, arma::fill::zeros);
  
  arma::cube Hl(q,q, L, arma::fill::zeros);
  arma::cube Ul(q,q, L, arma::fill::zeros);
  arma::cube Ql(q,p, L, arma::fill::zeros);
  
  arma::cube Sl(p,p, L, arma::fill::zeros);
  arma::field<arma::cube> S_jl(J);
  S_jl.fill(arma::cube(p, p, L, arma::fill::zeros));
  arma::field<arma::cube> S_jl_over_N_jl(J);
  S_jl_over_N_jl.fill(arma::cube(p, p, L, arma::fill::zeros));
  arma::cube Z1l(p,p, L, arma::fill::zeros);
  arma::cube Zl(p,p, L, arma::fill::zeros);
  
  arma::cube Wl(p,p, L, arma::fill::zeros);
  
  arma::cube Ml(q, p,L);
  arma::rowvec lambdal(L);
  arma::rowvec cl(L);
  arma::cube Dl(p,p,L);
  arma::mat Nl_j(J, L, arma::fill::zeros);
  
  
  for(int j=0; j <J; j++){
    // For a fixed j, XI_jil is an n_j x L matrix
    arma::mat subY = Y_grouped[j];       // This gives Y_j, which is a n_j x p matrix
    Nl      += arma::sum( XI_jil(j), 0); // This calculates [\sum_{i=1}^{n_j} XI_jil] for a fixed j for l = 1,...,L, and running the loop over j gives \sum_{j=1}^{J}\sum_{i=1}^{n_j} XI_jil. So Nl is 1xL vector.
    Nl_j.row(j) = arma::sum( XI_jil(j), 0);
    SumY_j_l.slice(j) = subY.t() * XI_jil(j);
    SumY_l  += subY.t() * XI_jil(j);  // p x L // This gives Y_j^T  XI_jil, which is \sum_{j=1}^{J}\sum_{i=1}^{n_j} XI_jil y_{ji}, which is p x L
  }
  
  cl             = Nl + nu0 ;    // This is c_l = N_l + \nu_0
  
  
  for(int l=0; l<L; l++){
    
    if(Nl(l)>0){
      Ybarl.col(l) = SumY_l.col(l)/Nl(l); // This calculates \bar{y}_l
    }
  }
  for(int j=0; j <J; j++){
    for(int l=0; l<L; l++){
      if(Nl_j(j, l)>0){
        Ybar_jl.slice(j).col(l) = SumY_j_l.slice(j).col(l)/Nl_j(j, l); // This calculates \bar{y}_l
      }
    }
  }
  
  for(int l=0; l<L; l++){
    for(int j=0; j <J; j++){
      Hl.slice(l) +=Nl_j(j,l) * X.row(j).t() * X.row(j);
      Ql.slice(l) +=Nl_j(j,l) * X.row(j).t() * Ybar_jl.slice(j).col(l).t();
      Z1l.slice(l) +=Nl_j(j,l) * Ybar_jl.slice(j).col(l) * Ybar_jl.slice(j).col(l).t();
    }
    Hl.slice(l) = Hl.slice(l) + gammal(l) * taul(l) *iG0;
    Ql.slice(l) = Ql.slice(l)  + gammal(l) * taul(l)* (iG0 * B0);
    Ml.slice(l) = solve(Hl.slice(l), Ql.slice(l));
    Ul.slice(l) = Hl.slice(l).i();
    Zl.slice(l) = Z1l.slice(l) + (gammal(l) * taul(l) * B0.t() * iG0 * B0) - (Ql.slice(l).t() * solve(Hl.slice(l), Ql.slice(l)));
  }
  
  for(int l = 0; l<L; l++){
    for(int j=0; j <J; j++){
      arma::mat subY = Y_grouped[j].t(); // p x Nj
      arma::mat temp_j(p, p, arma::fill::zeros);
      for(int ii=0; ii<subY.n_cols; ii++){
        arma::mat temp = ( subY.col(ii) - Ybar_jl.slice(j).col(l)) *
          (subY.col(ii) - Ybar_jl.slice(j).col(l)).t(); //This calculating (y_ji - \bar{y}_l)(y_ji - \bar{y}_l)^T for a fixed j, i, and l
        double xi = (XI_jil(j))(ii,l); // This is XI_{jil} for a fixed j, i, and l
        Sl.slice(l) += xi * temp; //This calculating XI_{jil}(y_ji - \bar{y}_l)(y_ji - \bar{y}_l)^T for a fixed j, i, and l
        temp_j += xi * temp;
      }// Finally this sums over i, for a fixed j, and fixed l
      S_jl(j).slice(l) = temp_j;
    }
    // Wl.slice(l) = arma::inv_sympd(
    //     (iW0)+
    //       Sl.slice(l) + Zl.slice(l)
    //   ); // This gives D_l^{-1}
    Wl.slice(l) = arma::inv_sympd(
      (iW0)+
        arma::symmatu(Sl.slice(l)) + arma::symmatu(Zl.slice(l))
    ); // This gives D_l^{-1}
  }
  //This sums over j, for a fixed l. Returning S_l as Sl.slice(l)
  // //Rcpp::Rcout<< Sl_Nl.slice(l) << ".....\nNl = " << Nl(l) << "\n.....\n";
  for(int j = 0; j<J; j++){
    for(int l = 0; l<L;l++){
      if(Nl_j(j,l) >0){
        S_jl_over_N_jl(j).slice(l) = S_jl(j).slice(l)/Nl_j(j,l);
      }
    }
  }// This completes the loop over l
  //Rcpp::Rcout << "Sl = " << Sl;
  
  Rcpp::List results = Rcpp::List::create(
    Rcpp::_["Ybar_jl"] = Ybar_jl,
    Rcpp::_["Zl"] = Zl,
    Rcpp::_["Sl"] = Sl,
    Rcpp::_["Wl"] = Wl,
    Rcpp::_["Ml"]  = Ml,
    Rcpp::_["Ul"] = Ul,
    Rcpp::_["cl"] = cl.t(),
    Rcpp::_["Hl"] = Hl,
    Rcpp::_["Ql"] = Ql,
    Rcpp::_["Ybar"] = Ybarl,
    Rcpp::_["Nl"] = Nl,
    Rcpp::_["Nl_j"] = Nl_j,
    Rcpp::_["S_jl_over_N_jl"] = S_jl_over_N_jl
  );
  return(results);
  
}
// --------------------------------------------------------------------------- 8
// [[Rcpp::export]]
arma::mat  Update_RHOjk_cpp(arma::field<arma::mat> XI_jil,
                            arma::mat X,
                            arma::mat ml,
                            arma::colvec tl,
                            arma::colvec cl,
                            arma::cube Dl,
                            
                            arma::colvec ElnPI_k,
                            arma::mat    ElnOM_lk,
                            int const L,
                            int const J,
                            int const H){


  int q = X.n_cols; //Dimension of covariate X
  arma::mat N_jl(J,L);

  // To find the matrix \sum_{i=1}^{n_j} XI_jil
  for(int j=0; j<J; j++){
    N_jl.row(j) = arma::sum(XI_jil(j),0);
  }// This gives a JxL Matrix


  arma::mat unn_log_RHO_jk(J,H);
  arma::mat log_RHO_jk(J,H);
  arma::mat RHO_jk(J,H);
  // ElnOM_lk contains the matrix [g(\bar{a}_lh, \bar{b}_lh) + \sum_{t=1}^{l-1}g(\bar{b}_th, \bar{a}_th)]
  // For l = 1,...,L; h = 1,...,H
  // So ElnOM_lk is an LxH matrix
  // Z is the matrix [\sum_{l=1}^{L}\sum_{i=1}^{n_j} XI_jil {g(\bar{a}_lh, \bar{b}_lh) + \sum_{t=1}^{l-1}g(\bar{b}_th, \bar{a}_th)}]
  // for j = 1,...,J; h=1,...,H
  arma::mat Z = N_jl * ElnOM_lk; // This gives a JxH matrix

  // ElnPI_k is a Hx1 matrix and ElnPI_k(k) contains g(\bar{a}_k, \bar{b}_k) + \sum_{r=1}^{k-1}g(\bar{b}_r, \bar{a}_r)
  for(int k = 0; k < H; k++){
    unn_log_RHO_jk.col(k) =  Z.col(k) +  ElnPI_k(k);
  }// This gives the matrix 
  // g(\bar{a}_k, \bar{b}_k) + \sum_{r=1}^{k-1}g(\bar{b}_r, \bar{a}_r) + 
  // \sum_{l=1}^{L}\sum_{i=1}^{n_j} XI_jil {g(\bar{a}_lh, \bar{b}_lh) + \sum_{t=1}^{l-1}g(\bar{b}_th, \bar{a}_th)}
  
  arma::mat TT(J,H);
  
  for(int k= 0; k<H; k++){
    
    TT.col(k) = E_log_p_X_cpp_mvt(X,
           ml.col(k), tl(k), cl(k), Dl.slice(k));
    
  }// This calculates 0.5[\sum_{i=1}^{q}\Psi(c_h^x - i + 1)/2 + qlog2 + log|D_h^x| -
  //  q/t_h^x + c_h^x (x_j - m_h^x)^T D_h^x (x_j - m_h^x)]; for j = 1,...,J; h=1,..,H
  
  arma::rowvec logunn(H);
  for(int j = 0; j < J; j++){
    // Adding the jth row of TT and unn_log_RHO_jk, which gives a JxH matrix
    logunn = unn_log_RHO_jk.row(j) + TT.row(j); ;
    log_RHO_jk.row(j) = logunn - LogSumExp_cpp(logunn); // Normalizing by LogSumExp trick
  }
  // Returning the exp of normalized log_RHO_jk to gives the probabilities RHO_jk
  return(exp(log_RHO_jk));
}

// --------------------------------------------------------------------------- 8
// [[Rcpp::export]]
arma::mat  Update_RHOjk_CAM_cpp(arma::field<arma::mat> XI_jil,
                                arma::colvec ElnPI_k,
                                arma::mat    ElnOM_lk,
                                int const L,
                                int const J,
                                int const H){
  
  
  arma::mat N_jl(J,L);
  
  // To find the matrix \sum_{i=1}^{n_j} XI_jil
  for(int j=0; j<J; j++){
    N_jl.row(j) = arma::sum(XI_jil(j),0);
  }// This gives a JxL Matrix
  
  
  arma::mat unn_log_RHO_jk(J,H);
  arma::mat log_RHO_jk(J,H);
  arma::mat RHO_jk(J,H);
  // ElnOM_lk contains the matrix [g(\bar{a}_lh, \bar{b}_lh) + \sum_{t=1}^{l-1}g(\bar{b}_th, \bar{a}_th)]
  // For l = 1,...,L; h = 1,...,H
  // So ElnOM_lk is an LxH matrix
  // Z is the matrix [\sum_{l=1}^{L}\sum_{i=1}^{n_j} XI_jil {g(\bar{a}_lh, \bar{b}_lh) + \sum_{t=1}^{l-1}g(\bar{b}_th, \bar{a}_th)}]
  // for j = 1,...,J; h=1,...,H
  arma::mat Z = N_jl * ElnOM_lk; // This gives a JxH matrix
  
  // ElnPI_k is a Hx1 matrix and ElnPI_k(k) contains g(\bar{a}_k, \bar{b}_k) + \sum_{r=1}^{k-1}g(\bar{b}_r, \bar{a}_r)
  for(int k = 0; k < H; k++){
    unn_log_RHO_jk.col(k) =  Z.col(k) +  ElnPI_k(k);
  }// This gives the matrix 
  // g(\bar{a}_k, \bar{b}_k) + \sum_{r=1}^{k-1}g(\bar{b}_r, \bar{a}_r) + 
  // \sum_{l=1}^{L}\sum_{i=1}^{n_j} XI_jil {g(\bar{a}_lh, \bar{b}_lh) + \sum_{t=1}^{l-1}g(\bar{b}_th, \bar{a}_th)}
  
  arma::rowvec logunn(H);
  for(int j = 0; j < J; j++){
    // Adding the jth row of TT and unn_log_RHO_jk, which gives a JxH matrix
    logunn = unn_log_RHO_jk.row(j);
    log_RHO_jk.row(j) = logunn - LogSumExp_cpp(logunn); // Normalizing by LogSumExp trick
  }
  // Returning the exp of normalized log_RHO_jk to gives the probabilities RHO_jk
  return(exp(log_RHO_jk));
}
