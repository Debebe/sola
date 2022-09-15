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
  DATA_VECTOR(hiv_positive);
  DATA_VECTOR(hiv_tested);
  DATA_VECTOR(anc_private);
  DATA_VECTOR(anc_attended);
  DATA_MATRIX(X);
  DATA_IVECTOR(anc_obs_idx);                          // !Index to data with no NAs
  DATA_IVECTOR(hiv_obs_idx);                          // !Index to data with no NAs

  DATA_SPARSE_MATRIX(Q_space);
  DATA_SPARSE_MATRIX(Z_space_hiv);
  DATA_SPARSE_MATRIX(Z_space_anc);
  DATA_SPARSE_MATRIX(R_race);
  DATA_SPARSE_MATRIX(Z_space_race);


  Type val(0);
  PARAMETER_VECTOR(beta_hiv);
  PARAMETER_VECTOR(beta_anc);


  //••••Space-BYM2- hiv model•••••//
  PARAMETER(log_sigma_space_hiv);                                                   // marginal standard deviation
  PARAMETER(logit_phi_space_hiv);
  PARAMETER_VECTOR(b_hiv);                                                          // combined spatial effect
  PARAMETER_VECTOR(u_hiv);

  //••••Space-BYM2- anc model•••••//
    PARAMETER(log_sigma_space_anc);                                                   // marginal standard deviation
  PARAMETER(logit_phi_space_anc);
  PARAMETER_VECTOR(b_anc);                                                          // combined spatial effect
  PARAMETER_VECTOR(u_anc);


  //hiv-model- latent effect
  Type sigma_space_hiv(exp(log_sigma_space_hiv));
  val -= dnorm(sigma_space_hiv, Type(0.0), Type(1.0), true) + log_sigma_space_hiv;      // log_sigma: log-absolute Jacobian of exp(log_sigma)

  Type phi_space_hiv(invlogit(logit_phi_space_hiv));
  val -= log(phi_space_hiv) +  log(1 - phi_space_hiv);                                  // change of variables: logit_phi -> phi
  val -= dbeta(phi_space_hiv, Type(0.5), Type(0.5), true);

  val -= dnorm(sum(b_hiv), Type(0.0), Type(0.001) * b_hiv.size(), true);                // soft sum-to-zero constraint
  val -= bym2_conditional_lpdf(b_hiv, u_hiv, sigma_space_hiv, phi_space_hiv, Q_space);



  //anc-model-latent effect
  Type sigma_space_anc(exp(log_sigma_space_anc));
  val -= dnorm(sigma_space_anc, Type(0.0), Type(1.0), true) + log_sigma_space_anc;      // log_sigma: log-absolute Jacobian of exp(log_sigma)

  Type phi_space_anc(invlogit(logit_phi_space_anc));
  val -= log(phi_space_anc) +  log(1 - phi_space_anc);                                  // change of variables: logit_phi -> phi
  val -= dbeta(phi_space_anc, Type(0.5), Type(0.5), true);

  val -= dnorm(sum(b_anc), Type(0.0), Type(0.001) * b_anc.size(), true);                // soft sum-to-zero constraint
  val -= bym2_conditional_lpdf(b_anc, u_anc, sigma_space_anc, phi_space_anc, Q_space);

// space-race interaction

  PARAMETER(log_sigma_space_race);
  PARAMETER_ARRAY(u_raw_space_race);

  Type sigma_space_race(exp(log_sigma_space_race));
  val -= dnorm(sigma_space_race, Type(0.0), Type(2.5), true) + log_sigma_space_race;

  vector<Type> u_space_race(u_raw_space_race * sigma_space_race);
  if(u_raw_space_race.size() > 0)
    val += SEPARABLE(GMRF(R_race), GMRF(Q_space))(u_raw_space_race);

  //sum-to-zero- constraint on interaction term
  for (int i = 0; i < u_raw_space_race.cols(); i++) {
    val -= dnorm(u_raw_space_race.col(i).sum(), Type(0), Type(0.001) * u_raw_space_race.col(i).size(), true);
  }


  // hiv model -likelihood

  vector<Type> mu_hiv(X*beta_hiv +
                      Z_space_hiv *b_hiv +
                      Z_space_race*u_space_race);
  vector <Type> prevalence_hiv(invlogit(mu_hiv));

  for (int i = 0; i < hiv_obs_idx.size(); i++) {                                   // index to exclude 0 number tested
    val -= dbinom(hiv_positive[hiv_obs_idx[i]], hiv_tested[hiv_obs_idx[i]], prevalence_hiv[hiv_obs_idx[i]], true);
  }

  // anc model- likelihood
  vector<Type> mu_anc(X*beta_anc + Z_space_anc *b_anc);
  vector <Type> prevalence_anc(invlogit(mu_anc));

  //val -= dbinom(anc_private, anc_attended, prevalence_anc, true).sum();
  for (int i = 0; i < anc_obs_idx.size(); i++) {                                   // index to exclude 0 districts with 0 anc attendance (denom)
    val -= dbinom(anc_private[anc_obs_idx[i]], anc_attended[anc_obs_idx[i]], prevalence_anc[anc_obs_idx[i]], true);
  }

  REPORT(beta_hiv);
  REPORT(beta_anc);

  REPORT(prevalence_hiv);
  REPORT(prevalence_anc);

  ADREPORT(prevalence_hiv);
  ADREPORT(prevalence_anc);

  REPORT(b_hiv);
  REPORT(b_anc);
  REPORT(u_raw_space_race);

  ADREPORT(b_hiv);
  ADREPORT(b_anc);
  return val;
  }
