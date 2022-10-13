

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

  //••••Space-BYM2- •••••//
  DATA_SPARSE_MATRIX(Q_space);
  DATA_SPARSE_MATRIX(Z_space1);
  DATA_SPARSE_MATRIX(Z_space2);
  DATA_SPARSE_MATRIX(Z_space3);

  //space-race interaction(ICAR x IID)
  DATA_SPARSE_MATRIX(Z_space_race1);               // space race interaction // n_obs by  n_space
  DATA_SPARSE_MATRIX(Z_space_race2);               // space race interaction // n_obs by  n_space
  DATA_SPARSE_MATRIX(Z_space_race3);               // space race interaction // n_obs by  n_space
  DATA_SPARSE_MATRIX(R_race);


  PARAMETER(log_sigma_space1);                                                   // marginal standard deviation
  PARAMETER(logit_phi_space1);
  PARAMETER(log_sigma_space2);                                                   // marginal standard deviation
  PARAMETER(logit_phi_space2);
  PARAMETER(log_sigma_space3);                                                   // marginal standard deviation
  PARAMETER(logit_phi_space3);

  PARAMETER_VECTOR(b1);                                                          // combined spatial effect
  PARAMETER_VECTOR(u1);
  PARAMETER_VECTOR(b2);                                                          // combined spatial effect
  PARAMETER_VECTOR(u2);
  PARAMETER_VECTOR(b3);                                                          // combined spatial effect
  PARAMETER_VECTOR(u3);

  Type sigma_space1(exp(log_sigma_space1));
  nll -= dnorm(sigma_space1, Type(0.0), Type(1.0), true) + log_sigma_space1;

  Type phi_space1(invlogit(logit_phi_space1));
  nll -= log(phi_space1) +  log(1 - phi_space1);
  nll -= dbeta(phi_space1, Type(0.5), Type(0.5), true);
  nll -= dnorm(sum(b1), Type(0.0), Type(0.001) * b1.size(), true);                // soft sum-to-zero constraint
  nll -= bym2_conditional_lpdf(b1, u1, sigma_space1, phi_space1, Q_space);


  Type sigma_space2(exp(log_sigma_space2));
  nll -= dnorm(sigma_space2, Type(0.0), Type(1.0), true) + log_sigma_space2;      // log_sigma: log-absolute Jacobian of exp(log_sigma)
  Type phi_space2(invlogit(logit_phi_space2));
  nll -= log(phi_space2) +  log(1 - phi_space2);                                  // change of variables: logit_phi -> phi
  nll -= dbeta(phi_space2, Type(0.5), Type(0.5), true);
  nll -= dnorm(sum(b2), Type(0.0), Type(0.001) * b2.size(), true);                // soft sum-to-zero constraint
  nll -= bym2_conditional_lpdf(b2, u2, sigma_space2, phi_space2, Q_space);


  Type sigma_space3(exp(log_sigma_space3));
  nll -= dnorm(sigma_space3, Type(0.0), Type(1.0), true) + log_sigma_space3;      // log_sigma: log-absolute Jacobian of exp(log_sigma)
  Type phi_space3(invlogit(logit_phi_space3));
  nll -= log(phi_space3) +  log(1 - phi_space3);                                  // change of variables: logit_phi -> phi
  nll -= dbeta(phi_space3, Type(0.5), Type(0.5), true);
  nll -= dnorm(sum(b3), Type(0.0), Type(0.001) * b3.size(), true);                // soft sum-to-zero constraint
  nll -= bym2_conditional_lpdf(b3, u3, sigma_space3, phi_space3, Q_space);



  ///////////////////////////////////////////////////////
  /////////   space(ICAR)-race -interaction-////////////
  ///////////////////////////////////////////////////////

  PARAMETER(log_sigma_space_race1);
  PARAMETER_ARRAY(u_raw_space_race1);

  Type sigma_space_race1(exp(log_sigma_space_race1));
  nll -= dnorm(sigma_space_race1, Type(0.0), Type(2.5), true) + log_sigma_space_race1;

  vector<Type> u_space_race1(u_raw_space_race1 * sigma_space_race1);
  if(u_raw_space_race1.size() > 0)
    nll += SEPARABLE(GMRF(R_race), GMRF(Q_space))(u_raw_space_race1);

  //sum-to-zero- constraint on interaction term
  for (int i = 0; i < u_raw_space_race1.cols(); i++) {
    nll -= dnorm(u_raw_space_race1.col(i).sum(), Type(0), Type(0.001) * u_raw_space_race1.col(i).size(), true);
  }


  PARAMETER(log_sigma_space_race2);
  PARAMETER_ARRAY(u_raw_space_race2);

  Type sigma_space_race2(exp(log_sigma_space_race2));
  nll -= dnorm(sigma_space_race2, Type(0.0), Type(2.5), true) + log_sigma_space_race2;

  vector<Type> u_space_race2(u_raw_space_race2 * sigma_space_race2);
  if(u_raw_space_race2.size() > 0)
    nll += SEPARABLE(GMRF(R_race), GMRF(Q_space))(u_raw_space_race2);

  //sum-to-zero- constraint on interaction term
  for (int i = 0; i < u_raw_space_race2.cols(); i++) {
    nll -= dnorm(u_raw_space_race2.col(i).sum(), Type(0), Type(0.001) * u_raw_space_race2.col(i).size(), true);
  }

  PARAMETER(log_sigma_space_race3);
  PARAMETER_ARRAY(u_raw_space_race3);

  Type sigma_space_race3(exp(log_sigma_space_race3));
  nll -= dnorm(sigma_space_race3, Type(0.0), Type(2.5), true) + log_sigma_space_race3;

  vector<Type> u_space_race3(u_raw_space_race3 * sigma_space_race3);
  if(u_raw_space_race3.size() > 0)
    nll += SEPARABLE(GMRF(R_race), GMRF(Q_space))(u_raw_space_race3);

  //sum-to-zero- constraint on interaction term
  for (int i = 0; i < u_raw_space_race3.cols(); i++) {
    nll -= dnorm(u_raw_space_race3.col(i).sum(), Type(0), Type(0.001) * u_raw_space_race3.col(i).size(), true);
  }


// risk for district-age-race-survey combo
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

  log_risk1=(X*beta1+Z_space1 * b1 + Z_space_race1*u_space_race1);
  log_risk2=(X*beta2+Z_space2 * b2 + Z_space_race2*u_space_race2);
  log_risk3=(X*beta3+Z_space3 * b3 + Z_space_race3*u_space_race3);

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
  REPORT(beta3);

  REPORT(b1);
  REPORT(b2);
  REPORT(b3);


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
