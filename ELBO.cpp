#include <RcppDist.h> 
#include <RcppArmadilloExtensions/sample.h> 
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
using namespace Rcpp;
using namespace RcppArmadillo;
//------------------------------------------------------------------------------
// HELPER FUNCTIONS
//------------------------------------------------------------------------------
// [[Rcpp::export]]
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
// [[Rcpp::export]]
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
// [[Rcpp::export]]
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
// [[Rcpp::export]]
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

double lbeta_normconst_cpp(double a, double b){
  double C_ab = lgamma(a) + lgamma(b) - lgamma(a+b);
  return( - (C_ab));
}

arma::colvec lbeta_normconst_vec_cpp(arma::colvec a, arma::colvec b){
  arma::colvec C_ab = lgamma(a) + lgamma(b) - lgamma(a+b);
  return( - (C_ab));
}

arma::mat lbeta_normconst_mat_cpp(arma::mat a, arma::mat b){
  arma::mat C_ab = lgamma(a) + lgamma(b) - lgamma(a+b);
  return( - C_ab);
}

//------------------------------------------------------------------------------
// MAIN FUNCTIONS
//------------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
double elbo_p_Y(arma::field<arma::mat> XI_ijl,
                arma::mat ml,
                arma::colvec tl,
                arma::colvec cl,
                arma::cube Dl,
                arma::cube Sl,    // D x D x L
                arma::mat Ybar,   // D x L
                arma::rowvec Nl){

  // Here Sl is infact Sl_over_Nl from the output of Update_THETAl_Y_cpp_mvt
  // Nl is the output Nl from Update_THETAl_Y_cpp_mvt
  // Nl basically gives \sum_{j=1}^{J}\sum_{i=1}^{n_j} XI_{jil}
  double Z = 0.0;
  int p = Dl.n_rows;
  int L = cl.n_rows;

  arma::vec CSD(L);
  arma::vec ELDL(L);
  arma::vec H(L);

  for(int l=0; l < L; l++){
    CSD(l)  = Const_sum_digamma(p, cl(l));
    ELDL(l) = ElogDetLambda(Dl.slice(l), CSD(l));

    Z += Nl(l) *
      (
          ELDL(l) - p * (1.0/tl(l)) -
            cl(l) * arma::trace( Sl.slice(l) * Dl.slice(l) ) -
            cl(l) * arma::dot( (Ybar.col(l)-ml.col(l)),
                Dl.slice(l)*(Ybar.col(l)-ml.col(l)) ) -
                  p * log( 2.0 * arma::datum::pi)
      );

  }

  return(.5*Z);

}

// [[Rcpp::export]]
double elbo_p_X(arma::mat RHO_jk,
                arma::mat mk,
                arma::colvec tk,
                arma::colvec ck,
                arma::cube Dk,
                arma::cube Sk,    // q x q x H
                arma::mat Xbar,   // q x H
                arma::rowvec Nk){
  
  // Here Sk is infact Sk_over_Nk from the output of Update_THETAl_X_cpp_mvt
  // Nk is the output Nk from Update_THETAl_X_cpp_mvt
  // Nk basically gives \sum_{j=1}^{J} RHO_{jk}
  double Z = 0.0;
  int q = Dk.n_rows;
  int H = ck.n_rows;
  
  arma::vec CSD(H);
  arma::vec ELDL(H);
  arma::vec S(H);
  
  for(int k=0; k < H; k++){
    CSD(k)  = Const_sum_digamma(q, ck(k));
    ELDL(k) = ElogDetLambda(Dk.slice(k), CSD(k));
    
    Z += Nk(k) *
      (
          ELDL(k) - q * (1.0/tk(k)) -
            tk(k) * arma::trace( Sk.slice(k) * Dk.slice(k) ) -
            tk(k) * arma::dot( (Xbar.col(k)-mk.col(k)),
               Dk.slice(k)*(Xbar.col(k)-mk.col(k)) ) -
                 q * log( 2.0 * arma::datum::pi)
      );
    
  }
  
  return(.5*Z);
  
}
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
double elbo_diff_v(arma::colvec a_tilde_k,
                      arma::colvec b_tilde_k,
                      arma::colvec S_concDP){
  
  int H = a_tilde_k.n_rows;
  
  a_tilde_k.shed_row(H-1);
  b_tilde_k.shed_row(H-1);
  
  arma::colvec Y =
    E_log_beta(b_tilde_k,a_tilde_k) * ( (S_concDP[0]/S_concDP[1]) - 1.0 );
  
  double p1 = (H-1.0)  * ( R::digamma(S_concDP[0]) - log(S_concDP[1]) ) + arma::accu(Y);  
  double p2 = arma::accu(lbeta_normconst_vec_cpp(a_tilde_k,b_tilde_k) +
                         E_log_beta(a_tilde_k,b_tilde_k) % (a_tilde_k - 1.0 )+
                         E_log_beta(b_tilde_k,a_tilde_k) % (b_tilde_k - 1.0 ));
  
  return(p1-p2);
}

// -----------------------------------------------------------------------------
// [[Rcpp::export]]
double elbo_diff_u(arma::mat a_bar_Ulk,
                   arma::mat b_bar_Ulk,
                   arma::colvec R_concDP){
  
  int L = a_bar_Ulk.n_rows;
  int H = a_bar_Ulk.n_cols;
  
  a_bar_Ulk.shed_row(L-1);
  b_bar_Ulk.shed_row(L-1);
  
  arma::colvec R(H);
  
  for(int k = 0; k < H; k ++){
    
    R[k] =  arma::accu( E_log_beta(b_bar_Ulk.col(k),a_bar_Ulk.col(k)) );
    
  }

  double p1 = (H * (L-1.0))  * ( R::digamma(R_concDP[0]) - log(R_concDP[1]) ) + ( (R_concDP[0]/R_concDP[1]) - 1.0 ) * arma::accu(R);  
  
  arma::colvec S(H);
  
  for(int k = 0; k < H; k ++){
    
    S[k] =  arma::accu(lbeta_normconst_vec_cpp(a_bar_Ulk.col(k),b_bar_Ulk.col(k)) +
      E_log_beta(a_bar_Ulk.col(k),b_bar_Ulk.col(k)) % (a_bar_Ulk.col(k) - 1.0 )+
      E_log_beta(b_bar_Ulk.col(k),a_bar_Ulk.col(k)) % (b_bar_Ulk.col(k) - 1.0 ));
    
  }
  double p2 = arma::accu(S);
  
  return(p1-p2);
}

// [[Rcpp::export]]
double elbo_diff_alpha(arma::colvec conc_hyper,
                           arma::colvec S_concDP){
  double pa =
    conc_hyper[0] * log(conc_hyper[1]) - lgamma(conc_hyper[0]) +
    (conc_hyper[0] - 1.0) * ( R::digamma(S_concDP[0]) - log(S_concDP[1]) ) -
    conc_hyper[1] * S_concDP[0]/S_concDP[1];

  double qa = S_concDP[0] * log(S_concDP[1]) - lgamma(S_concDP[0]) +
    (S_concDP[0] - 1.0) * (R::digamma(S_concDP[0])-log(S_concDP[1])) -
    S_concDP[0];

  return(pa - qa);
}

// [[Rcpp::export]]
double elbo_diff_beta(arma::colvec conc_hyper,
                       arma::colvec R_concDP){
  double pa =
    conc_hyper[0] * log(conc_hyper[1]) - lgamma(conc_hyper[0]) +
    (conc_hyper[0] - 1.0) * ( R::digamma(R_concDP[0]) - log(R_concDP[1]) ) -
    conc_hyper[1] * R_concDP[0]/R_concDP[1];
  
  double qa = R_concDP[0] * log(R_concDP[1]) - lgamma(R_concDP[0]) +
    (R_concDP[0] - 1.0) * (R::digamma(R_concDP[0])-log(R_concDP[1])) -
    R_concDP[0];
  
  return(pa - qa);
}

// -----------------------------------------------------------------------------
// [[Rcpp::export]]
double elbo_diff_S(arma::mat RHO_jk,
                   arma::colvec ElnPI){

  arma::colvec mdot_k = arma::sum(RHO_jk, 0).t();
  double Z1 = arma::accu(mdot_k % ElnPI);
  double Z2 = arma::accu(RHO_jk % log(RHO_jk + 1e-12));
  return(Z1-Z2);
}

// -----------------------------------------------------------------------------
// [[Rcpp::export]]
double elbo_diff_M(arma::field<arma::mat> XI_jil,
                         arma::mat RHO_jk,
                         arma::mat ElnOM_lk){

  int L = ElnOM_lk.n_rows;
  int J = RHO_jk.n_rows;
  int K = RHO_jk.n_cols;

  arma::mat N_jl(J,L);
  arma::mat E_ln_omega(L,K);
  arma::colvec Z2(J);

  // To find the matrix \sum_{i=1}^{n_j} XI_jil
  for(int j=0; j<J; j++){
    N_jl.row(j) = arma::sum(XI_jil(j),0);
  }// This gives a JxL Matrix
  
  // ElnOM_lk contains the matrix [g(\bar{a}_lh, \bar{b}_lh) + \sum_{t=1}^{l-1}g(\bar{b}_th, \bar{a}_th)]
  // For l = 1,...,L; h = 1,...,H
  // So ElnOM_lk is an LxH matrix
  // Z is the matrix [\sum_{l=1}^{L}\sum_{i=1}^{n_j} XI_jil {g(\bar{a}_lh, \bar{b}_lh) + \sum_{t=1}^{l-1}g(\bar{b}_th, \bar{a}_th)}]
  // for j = 1,...,J; h=1,...,H
  arma::mat Z = N_jl * ElnOM_lk; // This gives a JxH matrix
  
  for(int j=0; j<J; j++){
    Z2(j) = arma::accu(XI_jil[j] % log(XI_jil[j] + 1e-12) );
  }

  arma::mat X_jl  = RHO_jk % Z;
  double Z1 = accu(X_jl);
  return(Z1- arma::accu(Z2));
}

// ----------------------------------------------------------------------- 10.74
// [[Rcpp::export]]
double elbo_p_THETA_Y(arma::colvec m0, double t0, double c0,
                    arma::mat D0, arma::mat iD0,
                    double lCpl0, double LlogB0,
                    arma::mat ml,
                    arma::colvec tl,
                    arma::colvec cl, arma::cube Dl){

  int p = D0.n_rows;
  int L = tl.n_rows;

  arma::vec CSD(L);
  arma::vec LCPL(L);

  double p1   = 0.0;
  double ell1 = 0.0;
  double p4   = 0.0;


  for(int l = 0; l<L; l++){

    CSD(l)  = Const_sum_digamma(p, cl(l));
    LCPL(l) = log_Const_prod_gamma(p, cl(l));

    ell1   += ElogDetLambda(Dl.slice(l), CSD(l));
    p1     += cl(l) * arma::dot((ml.col(l)-m0) ,
                  Dl.slice(l) * ((ml.col(l)-m0)));
    p4     += cl(l) * arma::trace( iD0 * Dl.slice(l) );

  }

  double p0 = .5* (L * p * log( t0/(2.0*arma::datum::pi) ) +
                   ell1 -
                   p * t0 * arma::accu(1.0/(tl)) -
                   t0 * p1 );


  double  Z =  LlogB0 + p0 + (c0 - p - 1.0) * 0.5 * ell1 - 0.5 * p4;

  return(Z);
}

// ----------------------------------------------------------------------- 10.75
// [[Rcpp::export]]
double elbo_p_THETA_X(arma::colvec m0, double t0, double c0,
                      arma::mat D0, arma::mat iD0,
                      double lCpl0, double LlogB0,
                      arma::mat mk,
                      arma::colvec tk,
                      arma::colvec ck, arma::cube Dk){
  
  int q = D0.n_rows;
  int H = tk.n_rows;
  
  arma::vec CSD(H);
  arma::vec LCPL(H);
  
  double p1   = 0.0;
  double ell1 = 0.0;
  double p4   = 0.0;
  
  
  for(int k = 0; k<H; k++){
    
    CSD(k)  = Const_sum_digamma(q, ck(k));
    LCPL(k) = log_Const_prod_gamma(q, ck(k));
    
    ell1   += ElogDetLambda(Dk.slice(k), CSD(k));
    p1     += ck(k) * arma::dot((mk.col(k)-m0) ,
                 Dk.slice(k) * ((mk.col(k)-m0)));
    p4     += ck(k) * arma::trace( iD0 * Dk.slice(k) );
    
  }
  
  double p0 = .5* (H * q * log( t0/(2.0*arma::datum::pi) ) +
                   ell1 -
                   q * t0 * arma::accu(1.0/(tk)) -
                   t0 * p1 );
  
  
  double  Z =  LlogB0 + p0 + (c0 - q - 1.0) * 0.5 * ell1 - 0.5 * p4;
  
  return(Z);
}

// ----------------------------------------------------------------------- 10.77
// [[Rcpp::export]]
double elbo_q_THETA_Y(arma::mat ml, arma::colvec tl,
                    arma::colvec cl, arma::cube Dl){

  int p = Dl.slice(0).n_rows;
  int L = cl.n_rows;
  arma::vec CSD(L);
  arma::vec ELDL(L);
  arma::vec H(L);

  for(int l = 0; l<L; l++){
    CSD(l)  = Const_sum_digamma(p, cl(l));
    ELDL(l) = ElogDetLambda(Dl.slice(l), CSD(l));
    H(l)    = H_Lambda(Dl.slice(l),
      cl(l), CSD(l),
      log_Const_prod_gamma(p,cl(l)));
  }

  arma::vec Z = .5 * ELDL +
    p * 0.5 * ( log( tl / (2 * arma::datum::pi)) - 1.0) -
    H; 

  return(arma::accu(Z));
}

// ----------------------------------------------------------------------- 10.78
// [[Rcpp::export]]
double elbo_q_THETA_X(arma::mat mk, arma::colvec tk,
                      arma::colvec ck, arma::cube Dk){
  
  int q = Dk.slice(0).n_rows;
  int H = ck.n_rows;
  arma::vec CSD(H);
  arma::vec ELDL(H);
  arma::vec S(H);
  
  for(int k = 0; k<H; k++){
    CSD(k)  = Const_sum_digamma(q, ck(k));
    ELDL(k) = ElogDetLambda(Dk.slice(k), CSD(k));
    S(k)    = H_Lambda(Dk.slice(k),
      ck(k), CSD(k),
      log_Const_prod_gamma(q,ck(k)));
  }
  
  arma::vec Z = .5 * ELDL +
    q * 0.5 * ( log( tk / (2 * arma::datum::pi)) - 1.0) -
    S; 
  
  return(arma::accu(Z));
}









// -----------------------------------------------------------------------------
// FOR VI WITH REGRESSION OF Y ON X
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// [[Rcpp::export]]
double elbo_p_Y_X(arma::field<arma::mat> XI_ijl,
                  arma::mat X,
                  arma::cube Ml,
                  arma::cube Ul,
                  arma::colvec cl,
                  arma::cube Dl,
                  arma::field<arma::cube> S_jl_over_N_jl,    // J x p x p x L
                  arma::cube Ybar_jl,   // D x L
                  arma::mat N_jl){
  
  // Here Sl is infact Sl_over_Nl from the output of Update_THETAl_Y_cpp_mvt
  // Nl is the output Nl from Update_THETAl_Y_cpp_mvt
  // Nl basically gives \sum_{j=1}^{J}\sum_{i=1}^{n_j} XI_{jil}
  double Z = 0.0;
  int p = Dl.n_rows;
  int L = cl.n_rows;
  int J = N_jl.n_rows;
  
  for(int j = 0; j<J; j++){
    arma::vec CSD(L);
    arma::vec ELDL(L);
    arma::vec H(L);
    for(int l=0; l < L; l++){
      CSD(l)  = Const_sum_digamma(p, cl(l));
      ELDL(l) = ElogDetLambda(Dl.slice(l), CSD(l));
      
      Z += N_jl(j, l) *
        (
            ELDL(l) - p * arma::trace(X.row(j) * Ul.slice(l) * X.row(j).t()) -
              cl(l) * arma::trace(S_jl_over_N_jl(j).slice(l) * Dl.slice(l) ) -
              cl(l) * arma::dot( (Ybar_jl.slice(j).col(l)- Ml.slice(l).t() * X.row(j).t()),
                 Dl.slice(l)*(Ybar_jl.slice(j).col(l)- Ml.slice(l).t() * X.row(j).t()) ) -
                   p * log( 2.0 * arma::datum::pi)
        );
      
    }
  }
  
  return(.5*Z);
  
}

// --------------------------------------------------------------------------- 
//THIS IS FOR SCALAR VALUES OF TAU0 AND GAMMA0
// --------------------------------------------------------------------------- 
// // [[Rcpp::export]]
// double elbo_p_THETA_Y_X(arma::mat B0, 
//                         arma::mat iG0,
//                         double tau0, double gamma0,
//                         arma::mat iPsi0,
//                         double lG0,
//                         double LlogB0,
//                         arma::cube Ml,
//                         arma::cube Ul,
//                         double nu0,
//                         arma::colvec cl, arma::cube Dl){
//   
//   int p = B0.n_cols;
//   int q = B0.n_rows;
//   
//   int L = cl.n_rows;
//   
//   arma::vec CSD(L);
//   arma::vec LCPL(L);
//   
//   double p1   = 0.0;
//   double p2   = 0.0;
//   double ell1 = 0.0;
//   double p4   = 0.0;
//   
//   
//   for(int l = 0; l<L; l++){
//     
//     CSD(l)  = Const_sum_digamma(p, cl(l));
//     LCPL(l) = log_Const_prod_gamma(p, cl(l));
//     
//     ell1   += ElogDetLambda(Dl.slice(l), CSD(l));
//     p1     += cl(l) *  gamma0 * tau0* arma::trace(Dl.slice(l) * (Ml.slice(l)-B0).t() * iG0 * ((Ml.slice(l)-B0)));
//     p2     +=  gamma0 * tau0 * p * arma::trace( iG0 * Ul.slice(l) );
//     p4     += cl(l) * arma::trace( iPsi0 * Dl.slice(l) );
//     
//   }
//   
//   double p0 = .5* ((p * L * tau0)  + (q * L * gamma0) - (p * L * lG0) - (p * q * L * log(2.0*arma::datum::pi)) );
//   
//   
//   double  Z =  LlogB0 + p0 + ((nu0 - p - 1.0)*0.5*ell1) + (q*0.5*ell1) - (0.5*p4) - (0.5*p2) - (0.5*p1);
//   
//   return(Z);
// }
// --------------------------------------------------------------------------- 
//THIS IS FOR GENERAL VECTOR OF TAU AND GAMMA
// --------------------------------------------------------------------------- 
// [[Rcpp::export]]
double elbo_p_THETA_Y_X(arma::mat B0, 
                        arma::mat iG0,
                        arma::vec taul, arma::vec gammal,
                        arma::mat iPsi0,
                        double lG0,
                        double LlogB0,
                        arma::cube Ml,
                        arma::cube Ul,
                        double nu0,
                        arma::colvec cl, arma::cube Dl){
  
  int p = B0.n_cols;
  int q = B0.n_rows;
  
  int L = cl.n_rows;
  
  arma::vec CSD(L);
  arma::vec LCPL(L);
  
  double p1   = 0.0;
  double p2   = 0.0;
  double ell1 = 0.0;
  double p4   = 0.0;
  
  
  for(int l = 0; l<L; l++){
    
    CSD(l)  = Const_sum_digamma(p, cl(l));
    LCPL(l) = log_Const_prod_gamma(p, cl(l));
    
    ell1   += ElogDetLambda(Dl.slice(l), CSD(l));
    p1     += cl(l) *  gammal(l) * taul(l)* arma::trace(Dl.slice(l) * (Ml.slice(l)-B0).t() * iG0 * ((Ml.slice(l)-B0)));
    p2     +=  gammal(l) * taul(l) * arma::trace( iG0 * Ul.slice(l) );
    p4     += cl(l) * arma::trace( iPsi0 * Dl.slice(l) );
    
  }
  
  double p0 = .5* ((p * arma::accu(taul))  + (q * arma::accu(gammal)) - (p * L * lG0) - (p * q * L * log(2.0*arma::datum::pi)) );
  
  
  double  Z =  LlogB0 + p0 + ((nu0 - p - 1.0)*0.5*ell1) + (q*0.5*ell1) - (0.5*p4) - (0.5*p2) - (0.5*p1);
  
  return(Z);
}

// ----------------------------------------------------------------------- 10.77
// [[Rcpp::export]]
double elbo_q_THETA_Y_X(arma::cube Ul,
                        arma::colvec cl, arma::cube Dl){
  
  int p = Dl.slice(0).n_rows;
  int q = Ul.slice(0).n_rows;
  int L = cl.n_rows;
  arma::vec CSD(L);
  arma::vec ELDL(L);
  arma::vec H(L);
  arma::vec ldetU(L);
  
  for(int l = 0; l<L; l++){
    CSD(l)  = Const_sum_digamma(p, cl(l));
    ELDL(l) = ElogDetLambda(Dl.slice(l), CSD(l));
    ldetU(l) = log(arma::det(Ul.slice(l)));
    H(l)    = H_Lambda(Dl.slice(l),
      cl(l), CSD(l),
      log_Const_prod_gamma(p,cl(l)));
  }
  
  arma::vec Z = .5 * q * ELDL - .5* p * q * log((2 * arma::datum::pi) + 1) - 
                .5 * q * ldetU  -
  H; 
  
  return(arma::accu(Z));
}
