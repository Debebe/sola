#include <TMB.hpp>
template<class Type>
  Type bym2_conditional_lpdf(const vector<Type> b,
                             const vector<Type> u,
                             const Type sigma,
                             const Type phi,
                             const Eigen::SparseMatrix<Type> Q) {

    Type val(0.0);

    // constant terms omitted: -0.5 * (n + rank(Q)) * log(2*pi) + 0.5 * log|Q|
      val += -0.5 * b.size() * (2 * log(sigma) + log(1 - phi));  // normalising constant
      val += -0.5 * (b * b).sum() / (sigma * sigma * (1 - phi)) ;
      val += (b * u).sum()*sqrt(phi) / (sigma * (1 - phi));
      val += -0.5 * (u * (Q * u)).sum();
      val += -0.5 * (u * u).sum()* phi / (1 - phi) ;


      return(val);
  }

template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  DATA_VECTOR(Y);
  DATA_VECTOR(N);
  DATA_MATRIX(X);
  DATA_IVECTOR(obs_idx);                          // !Index to data with no NAs


  DATA_SPARSE_MATRIX(Q_space);
  DATA_SPARSE_MATRIX(Z_space);

  
  Type val(0);
  PARAMETER_VECTOR(beta);


  //••••Space-BYM2•••••//
  PARAMETER(log_sigma_space);                                                   // marginal standard deviation
  PARAMETER(logit_phi_space);
  PARAMETER_VECTOR(b);                                                          // combined spatial effect
  PARAMETER_VECTOR(u);                                                          // spatially correlated component

  Type sigma_space(exp(log_sigma_space));
  val -= dnorm(sigma_space, Type(0.0), Type(1.0), true) + log_sigma_space;      // log_sigma: log-absolute Jacobian of exp(log_sigma)

  Type phi_space(invlogit(logit_phi_space));
  val -= log(phi_space) +  log(1 - phi_space);                                  // change of variables: logit_phi -> phi
  val -= dbeta(phi_space, Type(0.5), Type(0.5), true);

  val -= dnorm(sum(b), Type(0.0), Type(0.001) * b.size(), true);                // soft sum-to-zero constraint
  val -= bym2_conditional_lpdf(b, u, sigma_space, phi_space, Q_space);

  
  // model
  vector<Type> mu(X*beta + 
                  Z_space *b);

  vector <Type> prevalence(invlogit(mu));

  for (int i = 0; i < obs_idx.size(); i++) {           // index to complete data
    val -= dbinom(Y[obs_idx[i]], N[obs_idx[i]], prevalence[obs_idx[i]], true);

  }          

  REPORT(beta);
  REPORT(prevalence);
  ADREPORT(prevalence);
  ADREPORT(mu);
  return val;
}
