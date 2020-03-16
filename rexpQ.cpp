#include <iostream>

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>

#include <cmath>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <chrono>

using namespace std;

inline double h(double x) {
  return 1.0-x+x*log(x);
}

inline double hifunc(double rho, double B,double loge) {
  return rho +(B-loge)/3.0*(1.0+sqrt(1.0+18.0*rho/(B-loge)));
}
inline double lofunc(double rho, double A, double loge) {
  const double logroot2pi=0.5*log(2*3.14159265);
  return rho+sqrt(2*rho)*sqrt(-logroot2pi-loge-1.5*log(A)+log(A-1));
}
unsigned int get_mlo(unsigned int mhi, double rho) {
  unsigned int mlo;
  // Since using unsigned int, need to be careful of negative numbers
  double dmlo=double(2*(int)(rho-0.5))-(double)mhi;
  if (dmlo>0) {
    mlo=(unsigned int) dmlo;
  }
  else {
    mlo=0;
  }

  return mlo;
}
// Take the L1 norm of the matrix and return the mmax for the req precision
// mlo is lower bound on possible m value
unsigned int get_m(double rho, double prec, unsigned int mlo=0) {
  if (rho > 4.2e9) {
    return 0;
  }
  
  const double logprec=log(prec), pi=3.14159265;
  double dmhi, dmlo;
  unsigned int mhi;
  
  dmhi= hifunc(rho,0.0,logprec)-1;
  dmlo=lofunc(rho,2*rho*h(dmhi/rho),logprec);
  if ((unsigned int)dmlo > mlo) {
    mlo=(unsigned int)dmlo;
  }
  
  if (log(boost::math::gamma_p((double)(mlo+1),rho))<logprec) {
    return mlo; // lower bound is the value we want - no binary search needed
  }
  else {
    const double B=-0.5*log(4*pi*rho*h(dmlo/rho));
    if (B>logprec) {
      dmhi=hifunc(rho,B,logprec);
    }
    mhi=(unsigned int)(dmhi+1);
  
    //    cout<<mhi<<"-"<<mlo<<"=width: "<< mhi-mlo<<"\n";
  
    while (mhi-mlo>1) {
      unsigned int mmid=mlo+(mhi-mlo)/2; // =(mlo+mhi)/2, rounds down
      double dm=(double)mmid;
      double loginv;
      
      loginv=log(boost::math::gamma_p(dm+1,rho));
    //    cout <<mlo<<", "<<mmid <<", "<<mhi<<", "<<loginv<<"\n";

      if (loginv<logprec) {
	mhi=mmid;
      }
      else {
	mlo=mmid;
      }
    }
  }
  return mhi;
}



// Sanity checks on input parameters
bool commonChecks(double prec, int vn_rows, int vn_cols, int Qn_rows, int Qn_cols, string sinit, bool ignore_v=false) {
  bool bad=false;

  if ((prec>=1.0) || (prec < 1e-16)) {
    cout << sinit << "precision of "<<prec<<"\n";
    bad=true;
  }
  if (!ignore_v && (Qn_rows != vn_cols)) {
    cout << sinit<<"v and Q dimension mismatch "<<vn_rows<<"x"<<vn_cols<<" vs "<<Qn_rows<<"x"<<Qn_cols<<"\n";
    bad=true;
  }
  if (Qn_rows != Qn_cols) {
    cout << sinit<<"Q dimensions "<<Qn_rows<<"x"<<Qn_cols<<"\n";
    bad=true;
  }

  return bad;
}

bool sChecks(const arma::mat v, const arma::sp_mat Q, double prec, string caller, bool ignore_v=false) {
  string sinit="ERROR ("+caller+"): ";
  bool bad=commonChecks(prec, v.n_rows, v.n_cols, Q.n_rows, Q.n_cols, sinit,ignore_v);

  if (Q.n_rows == Q.n_cols) {
    arma::dvec dg(Q.diag());
    if (dg.max()>0) {
      cout << sinit <<"positive elements on Q diagonal, max= "<<dg.max()<<"\n";
      bad=true;
    }
    arma::sp_mat Qpd=Q-speye(arma::size(Q))*Q.min();
    if (Qpd.min()<0) {
      cout << sinit<<"negative elements off Q diagonal, min= "<<Qpd.min()<<"\n";
      bad=true;
    }
  }
  arma::dvec zeros=arma::zeros(Q.n_cols);
  arma::dvec ones=zeros+1;
  if (arma::norm(Q*ones)>1e-10) {
    cout << sinit<<"Q 1 != 0\n";
    bad=true;
  }
  if (!ignore_v && (v.min()<0)) {
      cout << sinit <<"negative elements in v, min= "<<v.min()<<"\n";
      bad=true;
  }

  return bad;
}
bool dChecks(const arma::mat v, const arma::dmat Q, double prec, string caller, bool ignore_v=false) {
  string sinit="ERROR ("+caller+"): ";
  bool bad=commonChecks(prec, v.n_rows, v.n_cols, Q.n_rows, Q.n_cols, sinit, ignore_v);

  if (Q.n_rows == Q.n_cols) {
    arma::vec dg=Q.diag();
    if (dg.max()>0) {
      cout << sinit <<"positive elements on Q diagonal, max= "<<dg.max()<<"\n";
      bad=true;
    }
    arma::dmat Qnd=Q-arma::diagmat(dg);
    if (Qnd.min()<0) {
      cout << sinit<< "negative elements off Q diagonal, min= "<<Qnd.min()<<"\n";
      bad=true;
    }
  }
  arma::dvec zeros=arma::zeros(Q.n_cols);
  arma::dvec ones=zeros+1;
  if (arma::norm(Q*ones)>1e-10) {
    cout << sinit<< "Q 1 != 0\n";
    bad=true;
  }
  if (!ignore_v && (v.min()<0)) {
      cout << sinit <<"negative elements in v, min= "<<v.min()<<"\n";
      bad=true;
  }
  
  return bad;
}
// return value when checking finds errors
static arma::mat drivel(void ) {
  return arma::zeros(1,1);
}

//////////////////////////////
/*** Scaling and squaring ***/
//////////////////////////////
// Choose number of squares, s according to
// Old code that was used in MESA
//inline int get_s_from_rho(double rho) {
//  return (rho>1) ? (int) (log(rho)/log(2)): 0;
//}
// Better algorithm
int get_s_from_rho_prec(double rho, double eps) {
  const double log2=log(2), logrho=log(rho);
  int smax=logrho/log2<-7?0:(int)(log(rho)/log2) +7;
  int smin=smax<7?0:smax-7, s=smin;
  double rhosmall=rho/exp((double)s*log2);
  int nops=get_m(rhosmall,eps)+s, bestnops=nops, bests=s;
  for (s=smin+1;s<=smax;s++) {
    rhosmall/=2.0;
    nops=get_m(rhosmall,eps)+s;
    //    cout <<"s="<<s<<", nops="<<nops<<"\n";
    if (nops<bestnops) {
      bests=s;
      bestnops=nops;
    }
  }
  // The above works for dense matrices.
  // However for sparse matrices nops \propto beta get_m + s, for beta <1.
  // The below typically adds one or two ops, so makes little diff for dense
  // but improves speed for sparse.
  return bests<2?0:bests-2; 
}

/////////////////
/* Full Matrix */
/////////////////

void renormalise_P(arma::dmat &P) {
  int j, d=P.n_rows;
  arma::dvec rowsums=arma::sum(P,1);
  for (j=0;j<d;j++) {
    P.row(j)/=rowsums(j);
  }  
}
void renormalise_v(arma::dmat &v, const arma::dvec &v0sums) {
  int j, d=v.n_rows;
  arma::dvec rowsums=arma::sum(v,1);
  for (j=0;j<d;j++) {
    v.row(j)*=(v0sums(j)/rowsums(j));
  }  
}

// Sparse Q
arma::mat sSS_exp_Q(const arma::sp_mat Q, double prec, bool renorm=true) {
  const int d=Q.n_cols;
  const double rho=-Q.min(), BIG=1e100;
  if (rho <= 0) {
    return arma::eye(size(Q));
  }
  const unsigned int s=get_s_from_rho_prec(rho,prec);
  const double rhosmall=rho/exp((double)s*log(2)); // rho for Q/2^s
  const arma::sp_mat Msmall=rhosmall*(Q/rho+speye(size(Q)));
  arma::mat Asum, Apro=arma::mat(Msmall);
  const unsigned int mhi=get_m(rhosmall,prec);
  double fac=2.0, szprod, rncumu=0;
  unsigned int j;
  
  szprod=rhosmall;
  if (szprod>BIG) {
    Apro/=szprod;
    szprod=1; rncumu=log(szprod);
  }
  Asum=Apro+arma::eye(d,d); // 0th and 1st term of exp(Msmall)

  j=2; // start at 2nd term
  while (j<=mhi) { // do vsum <- v^T e^P
    Apro=Apro*Msmall;
    Apro/=fac;
    Asum+=Apro; 
    szprod*=rhosmall/fac++; 
    if (szprod>BIG) {
      Apro/=szprod; 
      Asum/=szprod; 
      if (!renorm) {
	rncumu+=log(szprod);
      }
      szprod=1;
    }
    j++;
  }

  if (renorm) {
    renormalise_P(Asum);
  }
  else {
    Asum*=exp(rncumu-rhosmall);
  }

  for (j=0; j<s;j++) { // repeated squaring
    Asum *= Asum;
    if (renorm) {
      renormalise_P(Asum);
    }
  }

  return Asum;
}
// Dense Q
arma::mat dSS_exp_Q(const arma::dmat Q, double prec, bool renorm=true) {
  const int d=Q.n_cols;
  const double rho=-Q.min(), BIG=1e100;
  if (rho==0) {
    return arma::eye(size(Q));
  }
  const unsigned int s=get_s_from_rho_prec(rho,prec);
  const double rhosmall=rho/exp((double)s*log(2));
  const arma::dmat Msmall=rhosmall*(Q/rho+eye(size(Q)));
  arma::mat Asum;
  arma::mat Apro=arma::mat(Msmall);
  const unsigned int mhi=get_m(rhosmall,prec);
  double fac=2.0, szprod, rncumu=0;
  unsigned int j;
  
  szprod=rhosmall;
  if (szprod>BIG) {
    Apro/=szprod;
    szprod=1; rncumu=log(szprod);
  }
  Asum=Apro+arma::eye(d,d); // 0th and 1st term of exp(Msmall)

  j=2; // start at 2nd term
  while (j<=mhi) { // do vsum <- v^T e^P
    Apro=Apro*Msmall;
    Apro/=fac;
    Asum+=Apro; 
    szprod*=rhosmall/fac++; 
    if (szprod>BIG) {
      Apro/=szprod; 
      Asum/=szprod; 
      if (!renorm) {
	rncumu+=log(szprod);
      }
      szprod=1;
    }
    j++;
  }

  if (renorm) {
    renormalise_P(Asum);
  }
  else {
    Asum*=exp(rncumu-rhosmall);
  }

  for (j=0; j<s;j++) { // repeated squaring
    Asum *= Asum;
    if (renorm) {
      renormalise_P(Asum);
    }
  }

  return Asum;
}

//' Exponentiate a rate matrix
//'
//' Exponentiates a whole rate matrix by making all elements non-negative, then  scaling and squaring.
//'
//' @param Q Rate matrix (sparse or dense).
//' @param prec Required precision - missing mass in the Poisson distribution.
//' @param renorm Force elements of each row to sum to 1? Defaults to TRUE.
//' @return exp(Q) A dense matrix.
//'
//' @author Chris Sherlock
//'
//' @examples
//' Qd=matrix(nrow=2,ncol=2,data=c(-1,1,2,-2),byrow=TRUE); SS_exp_Q(Qd,1e-10)
//'
//' library("Matrix")
//' d=5; f=0.3; ones=rep(1,d)
//' Qs=abs(rsparsematrix(d,d,f))
//' diag(Qs)=0
//' Qsum=Qs%*%ones
//' diag(Qs)=-Qsum
//' SS_exp_Q(Qs,1e-15)
//'
//' \dontrun{
//' M=matrix(nrow=2,ncol=2,data=c(1,1,2,2),byrow=TRUE); SS_exp_Q(M,1e-10)
//' M=matrix(nrow=2,ncol=2,data=c(1,-1,-2,2),byrow=TRUE); SS_exp_Q(M,1e-10)
//' SS_exp_Q(Qs,1.5)
//' SS_exp_Q(Qs,-2.0)
//'}
// [[Rcpp::export]]
arma::mat SS_exp_Q(SEXP Q, double prec, bool renorm=true) {
  bool bad=false, checks=true; // always do checks as cost is negligible
  
  if (Rf_isS4(Q)) {
    if(Rf_inherits(Q, "dgCMatrix")) {
      if (checks) {
	arma::mat v=arma::zeros(1,1);
	bad=sChecks(v,Rcpp::as<arma::sp_mat>(Q),prec,"sSS_exp_Q",true);
      }
      if (!bad) {
	return sSS_exp_Q(Rcpp::as<arma::sp_mat>(Q),prec,renorm);
      }
    }
    else {
      cout << "ERROR (SS_exp_Q): Incompatible class of Q.\n";
    }
  }
  else {
    if (checks) {
      arma::mat v=arma::zeros(1,1);
      bad=dChecks(v,Rcpp::as<arma::dmat>(Q),prec,"dSS_exp_Q",true);
    }
    if (!bad) {
      return dSS_exp_Q(Rcpp::as<arma::dmat>(Q),prec,renorm);
    }
  }
 
  return drivel();
}


/////////////////////
/* vector x matrix */
/////////////////////

// Sparse version
arma::mat sSS_v_exp_Q(const arma::mat v, const arma::sp_mat Q, double prec,
		      bool renorm=true) {
  const int d=v.n_cols;
  const arma::vec v0sums=arma::sum(v,1);
  const double rho=-Q.min(), BIG=1e100;
  if (rho==0) {
    return v;
  }
  const unsigned int s=get_s_from_rho_prec(rho,prec);
  const unsigned int twops=(unsigned int)pow(2,s);
  const double rhosmall=rho/(double)twops;
  const arma::sp_mat Msmall=rhosmall*(Q/rho+speye(size(Q)));
  arma::mat Asum, Apro=arma::mat(Msmall);
  const unsigned int mhi=get_m(rhosmall,prec);
  double fac=2.0, szprod, rncumu=0;
  arma::mat vtmp=v;
  unsigned int j;
  
  szprod=rhosmall;
  if (szprod>BIG) {
    Apro/=szprod;
    szprod=1; rncumu=log(szprod);
  }
  Asum=Apro+arma::eye(d,d); // 0th and 1st term

  j=2;
  while (j<=mhi) { // do vsum <- v^T e^P
    Apro=Apro*Msmall; // dense * sparse
    Apro/=fac;
    Asum+=Apro; 
    szprod*=rhosmall/fac++; 
    if (szprod>BIG) {
      Apro/=szprod; 
      Asum/=szprod; 
      if (!renorm) {
	rncumu+=log(szprod);
      }
      szprod=1;
    }
    j++;
  }

  if (renorm) {
    renormalise_P(Asum);
  }
  else {
    Asum*=exp(rncumu-rhosmall);
  }

  // Minimise s1 d^3 + 2^s2 d^2 subject to s1+s2=s.
  const unsigned int s2=(unsigned int)min((log((double)d/log(2))/log(2)),(double)s); 
  const unsigned int s1=s-s2, np=(unsigned int)pow(2,s2);
  //    cout << "s, s1, s2, np: "<< s<< ", "<< s1 << ", "<< s2<<","<< np<<"\n";
  
  for (j=0; j<s1;j++) { // repeated squaring
    Asum *= Asum;
  }
  if (renorm && (s1>0)) {
    renormalise_P(Asum);
  }

  for (j=0;j<np;j++) { // vector-matrix product
    vtmp*=Asum;
  }
  if (renorm) {
    renormalise_v(vtmp,v0sums);
  }
  
  return vtmp;
}
// Dense version
arma::mat dSS_v_exp_Q(const arma::mat v, const arma::mat Q, double prec,
		      bool renorm=true) {
  const int d=v.n_cols;
  const arma::vec v0sums=arma::sum(v,1);
  const double rho=-Q.min(), BIG=1e100;
  if (rho==0) {
    return v;
  }
  const unsigned int s=get_s_from_rho_prec(rho,prec);
  const unsigned int twops=(unsigned int)pow(2,s);
  const double rhosmall=rho/(double)twops;
  const arma::mat Msmall=rhosmall*(Q/rho+eye(size(Q)));
  arma::mat Asum, Apro=arma::mat(Msmall);
  const unsigned int mhi=get_m(rhosmall,prec);
  double fac=2.0, szprod, rncumu=0;
  arma::mat vtmp=v;
  unsigned int j;

  szprod=rhosmall;
  if (szprod>BIG) {
    Apro/=szprod;
    szprod=1; rncumu=log(szprod);
  }
  Asum=Apro+arma::eye(d,d); // 0th and 1st term

  j=2;
  while (j<=mhi) { // do vsum <- v^T e^P
    Apro=Apro*Msmall; 
    Apro/=fac;
    Asum+=Apro; 
    szprod*=rhosmall/fac++; 
    if (szprod>BIG) {
      Apro/=szprod; 
      Asum/=szprod; 
      if (!renorm) {
	rncumu+=log(szprod);
      }
      szprod=1;
    }
    j++;
  }

  if (renorm) {
    renormalise_P(Asum);
  }
  else {
    Asum*=exp(rncumu-rhosmall);
  }

  // Minimise s1 d^3 + 2^s2 d^2 subject to s1+s2=s.
  const unsigned int s2=(unsigned int)min((log((double)d/log(2))/log(2)),(double)s); 
  const unsigned int s1=s-s2, np=(unsigned int)pow(2,s2);
  //    cout << "s, s1, s2, np: "<< s<< ", "<< s1 << ", "<< s2<<","<< np<<"\n";
  
  for (j=0; j<s1;j++) { // repeated squaring
    Asum *= Asum;
  }
  if (renorm && (s1>0)) {
    renormalise_P(Asum);
  }
  
  for (j=0;j<np;j++) { // vector-matrix product
    vtmp*=Asum;
  }
  if (renorm) {
    renormalise_v(vtmp,v0sums);
  }

  return vtmp;
}
//' Product of horizontal vector and exponential of a rate matrix
//'
//' Evaluates v Exp(Q) by making all elements of Q non-negative, then scaling down by a power of two, and finally using a mixture of squaring and vector-matrix products.
//'
//' @author Chris Sherlock
//'
//' @param v Non-negative horizontal vector (dense).
//' @param Q Rate matrix (sparse or dense).
//' @param prec Required precision - missing mass in the Poisson distribution.
//' @param renorm Force elements of each row to sum to 1? Defaults to TRUE.
//' @param checks Perform sanity checks on the arguments? Defaults to TRUE.
//' @return v exp(Q) A dense horizontal vector.
//'
//' @examples
//' v=runif(2); v=v/sum(d)
//' Qd=matrix(nrow=2,ncol=2,data=c(-1,1,2,-2),byrow=TRUE); SS_v_exp_Q(t(v),Qd,1e-10)
//'
//' library("Matrix")
//' d=5; f=0.3; ones=rep(1,d); v=runif(d)
//' Qs=abs(rsparsematrix(d,d,f))
//' diag(Qs)=0
//' Qsum=Qs%*%ones
//' diag(Qs)=-Qsum
//' SS_v_exp_Q(t(v),Qs,1e-15)
//'
//' \dontrun{
//' v=runif(2); 
//' M=matrix(nrow=2,ncol=2,data=c(1,1,2,2),byrow=TRUE); SS_v_exp_Q(t(v),M,1e-10)
//' M=matrix(nrow=2,ncol=2,data=c(1,-1,-2,2),byrow=TRUE); SS_v_exp_Q(t(v),M,1e-10)
//' d=5; f=0.3; ones=rep(1,d); v=runif(d)
//' Qs=abs(rsparsematrix(d,d,f))
//' diag(Qs)=0; Qsum=Qs%*%ones; diag(Qs)=-Qsum
//' SS_v_exp_Q(v,Qs,1e-15)
//' SS_v_exp_Q(t(v),Qs,1.5)
//' SS_v_exp_Q(t(v),Qs,-2.0)
//' v=-runif(d)
//' SS_v_exp_Q(t(v),Qs,1e-15)
//'}
// [[Rcpp::export]]
arma::mat SS_v_exp_Q(const arma::mat v, SEXP Q, double prec, bool renorm=true, bool checks=true) {
  bool bad=false;
  if (Rf_isS4(Q)) {
    if(Rf_inherits(Q, "dgCMatrix")) {
      if (checks) {
	bad=sChecks(v,Rcpp::as<arma::sp_mat>(Q),prec,"sSS_v_exp_Q");
      }
      if (!bad) {
	return sSS_v_exp_Q(v,Rcpp::as<arma::sp_mat>(Q),prec,renorm);
      }
    }
    else {
      cout << "ERROR (SS_v_exp_Q): Incompatible class of Q.\n";
    }
  }
  else {
    if (checks) {
      bad=dChecks(v,Rcpp::as<arma::dmat>(Q),prec,"dSS_v_exp_Q");
    }
    if (!bad) {
      return dSS_v_exp_Q(v,Rcpp::as<arma::dmat>(Q),prec,renorm);
    }
  }
 
  return drivel();
}

//////////////////////////
//*** Uniformisation ***//
//////////////////////////
// v is one or more row vectors

// Sparse Q
arma::mat sUnif_v_exp_Q(const arma::mat v, const arma::sp_mat Q, double prec, bool renorm=true, bool t2=true) {
  const double rho=-Q.min(), BIG=1e100;
  if (rho==0) {
    return v;
  }
  const int d=v.n_cols, nc=v.n_rows;
  const arma::dvec v0sums=arma::sum(v,1);
  const unsigned int mhi=get_m(rho,prec/(1.0+(double)t2)),mlo=t2? get_mlo(mhi,rho): 0;
  const double dnc=(double)nc;
  unsigned int j;
  const arma::sp_mat P=Q+rho*speye(size(Q));
  arma::mat vsum(nc,d),vpro=v;
  double fac=1.0, szvpro, rncumu=0;

  //  cout << "rho="<<rho<<", mlo="<< mlo<< " ,mhi="<<mhi << ", d="<<d<<"\n";

  szvpro=arma::accu(v)/dnc;
  if (szvpro>BIG) {
    vpro/=szvpro;
    szvpro=1; rncumu=log(szvpro);
  }
  if (mlo==0) {
    vsum=vpro;
  }
  else {
    vsum.zeros();
  }

  j=1;
  while (j<=mhi) {

    vpro=vpro*P;
    vpro/=fac;
    if (j>=mlo) {
      vsum+=vpro; 
    }
    szvpro*=rho/fac++; 
    if ((szvpro>BIG)){
      vpro/=szvpro; 
      if (j>=mlo) { 
	vsum/=szvpro; 
      }
      if (!renorm) {
	rncumu+=log(szvpro);
      }
      szvpro=1;
    }
    j++;
  }

  if (renorm) {
    renormalise_v(vsum,v0sums);
  }
  else {
    vsum*=exp(rncumu-rho);
  }

  return vsum;
}
// Dense Q
arma::mat dUnif_v_exp_Q(const arma::mat v, const arma::mat Q, double prec, bool renorm=true, bool t2=true) {
  const double rho=-Q.min(), BIG=1e100;
  if (rho==0) {
    return v;
  }
  const int d=v.n_cols, nc=v.n_rows;
  const arma::dvec v0sums=arma::sum(v,1);
  const unsigned int mhi=get_m(rho,prec/(1.0+(double)t2)),mlo=t2? get_mlo(mhi,rho): 0;
  const double dnc=(double)nc;
  unsigned int j;
  const arma::mat P=Q+rho*eye(size(Q));
  arma::mat vsum(nc,d),vpro=v;
  double fac=1.0, szvpro, rncumu=0;
  
  //cout << "rho="<<rho<<", mlo="<< mlo<< " ,mhi="<<mhi << ", d="<<d<<"\n";

  szvpro=arma::accu(v)/dnc;
  if (szvpro>BIG) {
    vpro/=szvpro;
    szvpro=1; rncumu=log(szvpro);
  }
  if (mlo==0) {
    vsum=vpro;
  }
  else {
    vsum.zeros();
  }

  j=1;
  while (j<=mhi) {//&& (!cvgd)) { // do vsum <- v^T e^P
    vpro=vpro*P;
    vpro/=fac;
    if (j>=mlo) {
      vsum+=vpro; 
    }
    szvpro*=rho/fac++; 
    if ((szvpro>BIG)){
      vpro/=szvpro; 
      if (j>=mlo) { 
	vsum/=szvpro; 
      }
      if (!renorm) {
	rncumu+=log(szvpro);
      }
      szvpro=1;
    }
    j++;
  }

  if (renorm) {
    renormalise_v(vsum,v0sums);
  }
  else {
    vsum*=exp(rncumu-rho);
  }

  return vsum;
}
//' Product of horizontal vector and exponential of rate matrix
//'
//' Evaluates v exp(Q) by making all elements of Q non-negative, then using the uniformisation method.
//'
//' @author Chris Sherlock
//'
//' @param v Non-negative horizontal vector (dense).
//' @param Q Rate matrix (sparse or dense).
//' @param prec Required precision - missing mass in the Poisson distribution.
//' @param renorm Force elements of each row to sum to 1? Defaults to TRUE.
//' @param t2 Perform two-tailed truncation? Defaults to TRUE.
//' @param checks Perform sanity checks on the arguments? Defaults to TRUE.
//' @return v exp(Q) Dense horizontal vector.
//'
//' @examples
//' v=runif(2); v=v/sum(d)
//' Qd=matrix(nrow=2,ncol=2,data=c(-1,1,2,-2),byrow=TRUE); Unif_v_exp_Q(t(v),Qd,1e-10)
//'
//' library("Matrix")
//' d=5; f=0.3; ones=rep(1,d); v=runif(d)
//' Qs=abs(rsparsematrix(d,d,f))
//' diag(Qs)=0
//' Qsum=Qs%*%ones
//' diag(Qs)=-Qsum
//' Unif_v_exp_Q(t(v),Qs,1e-15)
//'
//' \dontrun{
//' v=runif(2); 
//' M=matrix(nrow=2,ncol=2,data=c(1,1,2,2),byrow=TRUE); Unif_v_exp_Q(t(v),M,1e-10)
//' M=matrix(nrow=2,ncol=2,data=c(1,-1,-2,2),byrow=TRUE); Unif_v_exp_Q(t(v),M,1e-10)
//' d=5; f=0.3; ones=rep(1,d); v=runif(d)
//' Qs=abs(rsparsematrix(d,d,f))
//' diag(Qs)=0; Qsum=Qs%*%ones; diag(Qs)=-Qsum
//' Unif_v_exp_Q(v,Qs,1e-15)
//' Unif_v_exp_Q(t(v),Qs,1.5)
//' Unif_v_exp_Q(t(v),Qs,-2.0)
//' v=-runif(d)
//' Unif_v_exp_Q(t(v),Qs,1e-15)
//'}
// [[Rcpp::export]]
arma::mat Unif_v_exp_Q(const arma::mat v, SEXP Q, double prec, bool renorm=true, bool t2=true, bool checks=true) {
  bool bad=false;
  if (Rf_isS4(Q)) {
    if(Rf_inherits(Q, "dgCMatrix")) {
      if (checks) {
	bad=sChecks(v,Rcpp::as<arma::sp_mat>(Q),prec,"sUnif_v_exp_Q");
      }
      if (!bad) {
	return sUnif_v_exp_Q(v,Rcpp::as<arma::sp_mat>(Q),prec,renorm,t2);
      }
    }
    else {
      cout << "ERROR (Unif_v_exp_Q): Incompatible class of Q.\n";
    }
  }
  else {
    if (checks) {
      bad=dChecks(v,Rcpp::as<arma::dmat>(Q),prec,"dUnif_v_exp_Q");
    }
    if (!bad) {
      return dUnif_v_exp_Q(v,Rcpp::as<arma::dmat>(Q),prec,renorm,t2);
    }
  }

  return drivel();
}

/////////////////////////////////////////////////////////////////
// vector transpose * matrix exponential - choose which algorithm
// test for choice depends on input matrix type
/////////////////////////////////////////////////////////////////
// Sparse version

bool should_use_SS_sparse(int d, double rho, double f) {
  const double log10=log(10.0);
  double log10d=log((double)d)/log10;
  double log10f=log(f)/log10;
  double log10rho=log(rho)/log10;
  bool method=0; // default to Uniform unless clear advantage
  
  if (log10d<3.34) { // don't even consider unless dimension is small
    if (log10rho>4.3) {
      method=1;
    }
    else if ((log10rho>4.0) && (log10d<2.9)) {
      method=1;
    }
    else if ((log10rho>3.67) && (log10d<2.8) && (log10f>-3.0)) {
      method=1;
    }
    else if ((log10rho>3.33) && (log10d<2.6) && (log10f>-2.5)) {
      method=1;
    }
    else if ((log10rho>3.0) && (log10d<2.2) && (log10f>-2.0)) {
      method=1;
    }
    else if ((log10rho>2.67) && (log10d<1.8) && (log10f>-1.5)) {
      method=1;
    }
  }
  return method;
}
bool should_use_SS_dense(int d, double rho) {
  const double log10=log(10.0);
  double log10d=log((double)d)/log10;
  double log10rho=log(rho)/log10;
  bool method=0; // default to Uniform unless clear advantage
  
  if (log10d<3.34) { // don't even consider unless dimension is small
    if (log10rho>3.67) {
      method=1;
    }
    else if ((log10rho>3.33) && (log10d<2.5)) {
      method=1;
    }
    else if ((log10rho>3.0) && (log10d<2.2)) {
      method=1;
    }
    else if ((log10rho>2.67) && (log10d<1.8)) {
      method=1;
    }
    else if ((log10rho>2.33) && (log10d<1.4)) {
      method=1;
    }
  }
  return method;
}

arma::mat sv_exp_Q(const arma::mat v, const arma::sp_mat Q, double prec, bool renorm=true, bool t2=true) {
  const int d=v.n_cols;
  double rho=-Q.min();
  if (rho>4.2e9) {
    cout << "WARNING (sv_exp_Q): rho of "<<rho<<" set to 4.2e9\n";
    rho=4.2e9;
  }
  
  // Obvious formula for sparsity CPU cost, c = k_1 f d^2 doesnt hold well
  // so decision is empirical based on decision point in many runs

  double f=(double) Q.n_nonzero / (double) Q.n_elem;
  if (should_use_SS_sparse(d,rho,f)) {
    return sSS_v_exp_Q(v, Q, prec, renorm);
  }
  else {
    return sUnif_v_exp_Q(v, Q, prec, renorm, t2);
  }
}

// Dense version
arma::mat dv_exp_Q(const arma::mat v, const arma::mat Q, double prec, bool renorm=true, bool t2=true) {
  const int d=v.n_cols;
  double rho=-Q.min();
  if (rho>4.2e9) {
    cout << "WARNING (dv_exp_Q): rho of "<<rho<<" set to 4.2e9\n";
    rho=4.2e9;
  }

  if (should_use_SS_dense(d,rho)) {
    return dSS_v_exp_Q(v, Q, prec, renorm);
  }
  else {
    return dUnif_v_exp_Q(v, Q, prec, renorm, t2);
  }
}
//' Product of horizontal vector and exponential of rate matrix
//'
//' Evaluates v exp(Q); automatically chooses between scaling and squaring or uniformisation.
//'
//' @author Chris Sherlock
//'
//' @param v Non-negative horizontal vector (dense).
//' @param Q Rate matrix (sparse or dense).
//' @param prec Required precision - missing mass in the Poisson distribution.
//' @param renorm Force elements of each row to sum to 1? Defaults to TRUE.
//' @param t2 Perform two-tailed truncation? Defaults to TRUE.
//' @param checks Perform sanity checks on the arguments? Defaults to TRUE.
//' @return v exp(Q) Dense horizontal vector.
//'
//' @examples
//' v=runif(2); v=v/sum(d)
//' Qd=matrix(nrow=2,ncol=2,data=c(-1,1,2,-2),byrow=TRUE); v_exp_Q(t(v),Qd,1e-10)
//'
//' library("Matrix")
//' d=5; f=0.3; ones=rep(1,d); v=runif(d)
//' Qs=abs(rsparsematrix(d,d,f))
//' diag(Qs)=0
//' Qsum=Qs%*%ones
//' diag(Qs)=-Qsum
//' v_exp_Q(t(v),Qs,1e-15)
//'
//' \dontrun{
//' v=runif(2); 
//' M=matrix(nrow=2,ncol=2,data=c(1,1,2,2),byrow=TRUE); v_exp_Q(t(v),M,1e-10)
//' M=matrix(nrow=2,ncol=2,data=c(1,-1,-2,2),byrow=TRUE); v_exp_Q(t(v),M,1e-10)
//' d=5; f=0.3; ones=rep(1,d); v=runif(d)
//' Qs=abs(rsparsematrix(d,d,f))
//' diag(Qs)=0; Qsum=Qs%*%ones; diag(Qs)=-Qsum
//' v_exp_Q(v,Qs,1e-15)
//' v_exp_Q(t(v),Qs,1.5)
//' v_exp_Q(t(v),Qs,-2.0)
//' v=-runif(d)
//' v_exp_Q(t(v),Qs,1e-15)
//'}
// [[Rcpp::export]]
arma::mat v_exp_Q(const arma::mat v, SEXP Q, double prec, bool renorm=true, bool t2=true, bool checks=true) {
  bool bad=false;
  if (Rf_isS4(Q)) {
    if(Rf_inherits(Q, "dgCMatrix")) {
      if (checks) {
	bad=sChecks(v,Rcpp::as<arma::sp_mat>(Q),prec,"sv_exp_Q");
      }
      if (!bad) {
	return sv_exp_Q(v,Rcpp::as<arma::sp_mat>(Q),prec,renorm,t2);
      }
    }
    else {
      cout << "ERROR (v_exp_Q): Incompatible class of Q.\n";
    }
    return drivel();
  }
  else {
    if (checks) {
      bad=dChecks(v,Rcpp::as<arma::dmat>(Q),prec,"dv_exp_Q");
    }
    if (!bad) {
      return dv_exp_Q(v,Rcpp::as<arma::dmat>(Q),prec,renorm,t2);
    }
  }

  return drivel();
}
//' Vector transpose * exponential of a rate matrix, all transposed.
//'
//' Evaluates [v' exp(Q)]'; automatically chooses between scaling and squaring or uniformisation.
//'
//' @author Chris Sherlock
//'
//' @param v Non-negative vertical vector (dense).
//' @param Q Rate matrix (sparse or dense).
//' @param prec Required precision - missing mass in the Poisson distribution.
//' @param renorm Force elements of each row to sum to 1? Defaults to TRUE.
//' @param t2 Perform two-tailed truncation? Defaults to TRUE.
//' @param checks Perform sanity checks on the arguments? Defaults to TRUE.
//' @return [v' exp(Q)]' Dense vertical vector.
//'
//' @examples
//' v=runif(2); v=v/sum(d)
//' Qd=matrix(nrow=2,ncol=2,data=c(-1,1,2,-2),byrow=TRUE); vT_exp_Q(v,Qd,1e-10)
//'
//' library("Matrix")
//' d=5; f=0.3; ones=rep(1,d); v=runif(d)
//' Qs=abs(rsparsematrix(d,d,f))
//' diag(Qs)=0; Qsum=Qs%*%ones; diag(Qs)=-Qsum
//' vT_exp_Q(v,Qs,1e-15)
//'
//' \dontrun{
//' v=runif(2); 
//' M=matrix(nrow=2,ncol=2,data=c(1,1,2,2),byrow=TRUE); vT_exp_Q(v,M,1e-10)
//' M=matrix(nrow=2,ncol=2,data=c(1,-1,-2,2),byrow=TRUE); vT_exp_Q(v,M,1e-10)
//' d=5; f=0.3; ones=rep(1,d); v=runif(d)
//' Qs=abs(rsparsematrix(d,d,f))
//' diag(Qs)=0; Qsum=Qs%*%ones; diag(Qs)=-Qsum
//' vT_exp_Q(t(v),Qs,1e-15)
//' vT_exp_Q(v,Qs,1.5)
//' vT_exp_Q(v,Qs,-2.0)
//' v=-runif(d)
//' vT_exp_Q(v,Qs,1e-15)
//'}
// [[Rcpp::export]]
arma::dvec vT_exp_Q(const arma::dvec &v, SEXP Q, double prec, bool renorm=true, bool t2=true, bool checks=true) {
  arma::dmat vT=v.t();
  return v_exp_Q(vT, Q, prec, renorm, t2, checks).t();
}


