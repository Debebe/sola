
#include <TMB.hpp>
template<class Type>
Type bym2_conditional_lpdf(const vector<Type> b,
                           const vector<Type> u,
                           const Type sigma,
                           const Type phi,
                           const Eigen::SparseMatrix<Type> Q) {
  
  Type nll(0.0);
  nll += -0.5 * b.size() * (2 * log(sigma) + log(1 - phi));  // normalising constant
  nll += -0.5 * (b * b).sum() / (sigma * sigma * (1 - phi)) ;
  nll += (b * u).sum()*sqrt(phi) / (sigma * (1 - phi));
  nll += -0.5 * (u * (Q * u)).sum();
  nll += -0.5 * (u * u).sum()* phi / (1 - phi) ;
  
  return(nll);
}
// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  
  DATA_MATRIX(Y);
  DATA_MATRIX(X);

  PARAMETER_VECTOR(beta1);
  PARAMETER_VECTOR(beta2);
  PARAMETER_VECTOR(beta3);
  
  
  Type nll=0;
  
  //••••Space-BYM2- hiv model•••••//
  DATA_SPARSE_MATRIX(Q_space);
  DATA_SPARSE_MATRIX(Z_space);
  PARAMETER(log_sigma_space);                                                   // marginal standard deviation
  PARAMETER(logit_phi_space);
  PARAMETER_VECTOR(b);                                                          // combined spatial effect
  PARAMETER_VECTOR(u);

  Type sigma_space(exp(log_sigma_space));
  nll -= dnorm(sigma_space, Type(0.0), Type(1.0), true) + log_sigma_space;      // log_sigma: log-absolute Jacobian of exp(log_sigma)

  Type phi_space(invlogit(logit_phi_space));
  nll -= log(phi_space) +  log(1 - phi_space);                                  // change of variables: logit_phi -> phi
  nll -= dbeta(phi_space, Type(0.5), Type(0.5), true);
  nll -= dnorm(sum(b), Type(0.0), Type(0.001) * b.size(), true);                // soft sum-to-zero constraint
  nll -= bym2_conditional_lpdf(b, u, sigma_space, phi_space, Q_space);

// risk for district-age-sex combo
  Type risk1=0;
  Type risk2=0;
  Type risk3=0;
  Type level=0;
  
  vector<Type> y(Y.cols());
  vector<Type> prob(Y.cols());
  matrix<Type> proportion(Y.rows(), Y.cols());
  
  
  vector<Type> log_risk1(Y.rows());
  vector<Type> log_risk2(Y.rows());
  vector<Type> log_risk3(Y.rows());
  
  log_risk1=(X*beta1+Z_space * b);
  log_risk2=(X*beta2+Z_space * b);
  log_risk3=(X*beta3+Z_space * b);
  
  vector<Type> r1(exp(log_risk1));
  vector<Type> r2(exp(log_risk2));
  vector<Type> r3(exp(log_risk3));
  
  
  for (int i=0; i<Y.rows(); i++) {
    risk1=r1(i);
    risk2=r2(i);
    risk3=r3(i);
    
    level=1 + risk1 + risk2+ risk3;
    
    prob(0)=risk1/level;
    prob(1)=risk2/level;
    prob(2)=risk3/level;
    prob(3)=1 - prob(0) - prob(1)-prob(2);
    y=Y.row(i);
    if(!isNA(y(0))) nll -= dmultinom(y, prob, true);
    //nll -= dmultinom(y, prob, true);
    proportion.row(i) = prob;
    REPORT(prob);
    
  }
  
  REPORT(beta1);
  REPORT(beta2);
  REPORT(proportion);
  return nll;
}

// dyn.load(dynlib("multinom_test"))
//   Y= as.matrix(DT[, 7:10])
//   X = model.matrix(~1+ethnicity+survey+age_group, DT)
//   
//   tmbdata <- list(Y =Y,X=X)
//   
//   tmbpar <- list(
//       beta1 = numeric(ncol(X)),
//       beta2 = numeric(ncol(X)),
//       beta3 = numeric(ncol(X)))
//   
//   obj <- MakeADFun(tmbdata, tmbpar,  DLL="multinom_test2")
//   f <- nlminb(obj$par, obj$fn, obj$gr)
//   
//   rep2 <- obj$report()