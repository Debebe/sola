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

  DATA_VECTOR(non_pregnant);
  DATA_VECTOR(all_women);

  DATA_MATRIX(X);      //hiv prev predictor
  DATA_MATRIX(Xanc);   //anc_place predictor
  DATA_MATRIX(Xpreg);   //pregnancy predictor


  DATA_IVECTOR(anc_obs_idx);                           // !Index to data with no NAs
  DATA_IVECTOR(hiv_obs_idx);                           // !Index to data with no NAs
  DATA_IVECTOR(preg_obs_idx);                          // !Index to data with no NAs



  DATA_SPARSE_MATRIX(Q_space);
  DATA_SPARSE_MATRIX(Z_space_hiv);
  DATA_SPARSE_MATRIX(Z_space_anc);
  DATA_SPARSE_MATRIX(Z_space_preg);

  //hiv model
  DATA_SPARSE_MATRIX(Z_space_race_hiv);               // space race interaction // n_obs by  n_space
  DATA_SPARSE_MATRIX(R_race_hiv);
  DATA_SPARSE_MATRIX(Z_space_ancplace_hiv);          // space race interaction // n_obs by  n_space
  DATA_SPARSE_MATRIX(R_ancplace_hiv);


  //anc model
  DATA_SPARSE_MATRIX(Z_space_race_anc);              // space race interaction // n_obs by  n_space
  DATA_SPARSE_MATRIX(R_race_anc);

  //preg model
  DATA_SPARSE_MATRIX(Z_space_race_preg);              // space race interaction // n_obs by  n_space
  DATA_SPARSE_MATRIX(R_race_preg);


  Type val(0);
  PARAMETER_VECTOR(beta_hiv);
  PARAMETER_VECTOR(beta_anc);
  PARAMETER_VECTOR(beta_preg);

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

  //••••Space-BYM2- preg model•••••//
  PARAMETER(log_sigma_space_preg);                                                   // marginal standard deviation
  PARAMETER(logit_phi_space_preg);
  PARAMETER_VECTOR(b_preg);                                                          // combined spatial effect
  PARAMETER_VECTOR(u_preg);


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




  //preg-model-latent effect
    Type sigma_space_preg(exp(log_sigma_space_preg));
    val -= dnorm(sigma_space_preg, Type(0.0), Type(1.0), true) + log_sigma_space_preg;      // log_sigma: log-absolute Jacobian of exp(log_sigma)

    Type phi_space_preg(invlogit(logit_phi_space_preg));
    val -= log(phi_space_preg) +  log(1 - phi_space_preg);                                  // change of variables: logit_phi -> phi
    val -= dbeta(phi_space_preg, Type(0.5), Type(0.5), true);

    val -= dnorm(sum(b_preg), Type(0.0), Type(0.001) * b_preg.size(), true);                // soft sum-to-zero constraint
    val -= bym2_conditional_lpdf(b_preg, u_preg, sigma_space_preg, phi_space_preg, Q_space);


  ///////////////////////////////////////////////////////
  /////////   ICAR-race(fixed) -interaction -hiv model  ////////////
  ///////////////////////////////////////////////////////

  PARAMETER(log_sigma_space_race_hiv);
  PARAMETER_ARRAY(u_raw_space_race_hiv);

  Type sigma_space_race_hiv(exp(log_sigma_space_race_hiv));
  val -= dnorm(sigma_space_race_hiv, Type(0.0), Type(2.5), true) + log_sigma_space_race_hiv;

  vector<Type> u_space_race_hiv(u_raw_space_race_hiv * sigma_space_race_hiv);
  if(u_raw_space_race_hiv.size() > 0)
    val += SEPARABLE(GMRF(R_race_hiv), GMRF(Q_space))(u_raw_space_race_hiv);

  //sum-to-zero- constraint on interaction term
  for (int i = 0; i < u_raw_space_race_hiv.cols(); i++) {
    val -= dnorm(u_raw_space_race_hiv.col(i).sum(), Type(0), Type(0.001) * u_raw_space_race_hiv.col(i).size(), true);
  }

   //Icar * anc-place (women group)- hiv model
  PARAMETER(log_sigma_space_ancplace_hiv);
  PARAMETER_ARRAY(u_raw_space_ancplace_hiv);

  Type sigma_space_ancplace_hiv(exp(log_sigma_space_ancplace_hiv));
  val -= dnorm(sigma_space_ancplace_hiv, Type(0.0), Type(2.5), true) + log_sigma_space_ancplace_hiv;

  vector<Type> u_space_ancplace_hiv(u_raw_space_ancplace_hiv * sigma_space_ancplace_hiv);
  if(u_raw_space_ancplace_hiv.size() > 0)
    val += SEPARABLE(GMRF(R_ancplace_hiv), GMRF(Q_space))(u_raw_space_ancplace_hiv);

  //sum-to-zero- constraint on interaction term
  for (int i = 0; i < u_raw_space_ancplace_hiv.cols(); i++) {
    val -= dnorm(u_raw_space_ancplace_hiv.col(i).sum(), Type(0), Type(0.001) * u_raw_space_ancplace_hiv.col(i).size(), true);
  }

  // anc_place(IID)*age(rw1)

  // DATA_SPARSE_MATRIX(Z_agegroup_ancplace);
  // DATA_SPARSE_MATRIX(R_agegroup);
  // DATA_SPARSE_MATRIX(R_ancplace);
  
  // PARAMETER(log_prec_agegroup_ancplace);
  // PARAMETER_ARRAY(u_raw_agegroup_ancplace);


  // val -= dlgamma(log_prec_agegroup_ancplace, Type(0.001), Type(1.0 / 0.001), true);
  // Type sigma_agegroup_ancplace(exp(-0.5 * log_prec_agegroup_ancplace));
  // vector<Type> u_agegroup_ancplace(u_raw_agegroup_ancplace * sigma_agegroup_ancplace);
  // val += SEPARABLE(GMRF(R_agegroup), GMRF(R_ancplace))(u_raw_agegroup_ancplace);



  ///////////////////////////////////////////////////////
  /////////   space(ICAR)-race -interaction-anc model////////////
  ///////////////////////////////////////////////////////

  PARAMETER(log_sigma_space_race_anc);
  PARAMETER_ARRAY(u_raw_space_race_anc);

  Type sigma_space_race_anc(exp(log_sigma_space_race_anc));
  val -= dnorm(sigma_space_race_anc, Type(0.0), Type(2.5), true) + log_sigma_space_race_anc;

  vector<Type> u_space_race_anc(u_raw_space_race_anc * sigma_space_race_anc);
  if(u_raw_space_race_anc.size() > 0)
    val += SEPARABLE(GMRF(R_race_anc), GMRF(Q_space))(u_raw_space_race_anc);

  //sum-to-zero- constraint on interaction term
  for (int i = 0; i < u_raw_space_race_anc.cols(); i++) {
    val -= dnorm(u_raw_space_race_anc.col(i).sum(), Type(0), Type(0.001) * u_raw_space_race_anc.col(i).size(), true);
  }



  ///////////////////////////////////////////////////////
  /////////   space(ICAR)-race -interaction-preg model////////////
  ///////////////////////////////////////////////////////

  PARAMETER(log_sigma_space_race_preg);
  PARAMETER_ARRAY(u_raw_space_race_preg);

  Type sigma_space_race_preg(exp(log_sigma_space_race_preg));
  val -= dnorm(sigma_space_race_preg, Type(0.0), Type(2.5), true) + log_sigma_space_race_preg;

  vector<Type> u_space_race_preg(u_raw_space_race_preg * sigma_space_race_preg);
  if(u_raw_space_race_preg.size() > 0)
    val += SEPARABLE(GMRF(R_race_preg), GMRF(Q_space))(u_raw_space_race_preg);

  //sum-to-zero- constraint on interaction term
  for (int i = 0; i < u_raw_space_race_preg.cols(); i++) {
    val -= dnorm(u_raw_space_race_preg.col(i).sum(), Type(0), Type(0.001) * u_raw_space_race_preg.col(i).size(), true);
  }


  // hiv model -likelihood

  vector<Type> mu_hiv(X*beta_hiv +
                      Z_space_hiv * b_hiv +
                      Z_space_race_hiv * u_space_race_hiv +
                      Z_space_ancplace_hiv * u_space_ancplace_hiv);
                      //Z_agegroup_ancplace * u_agegroup_ancplace);

  vector <Type> prevalence_hiv(invlogit(mu_hiv));

  for (int i = 0; i < hiv_obs_idx.size(); i++) {                                   // index to exclude 0 number tested
     val -= dbinom(hiv_positive[hiv_obs_idx[i]], hiv_tested[hiv_obs_idx[i]], prevalence_hiv[hiv_obs_idx[i]], true);
  }

  // anc model- likelihood
  vector<Type> mu_anc(Xanc*beta_anc +
                      Z_space_anc *b_anc +
                      Z_space_race_anc * u_space_race_anc);
  vector <Type> prevalence_anc(invlogit(mu_anc));

  for (int i = 0; i < anc_obs_idx.size(); i++) {                                   // index to exclude districts with 0|NA anc attendance (denom)
     val -= dbinom(anc_private[anc_obs_idx[i]], anc_attended[anc_obs_idx[i]], prevalence_anc[anc_obs_idx[i]], true);
  }

    // preg model- likelihood
    vector<Type> mu_preg(Xpreg*beta_preg +
                        Z_space_preg *b_preg +
                        Z_space_race_preg * u_space_race_preg);
    vector <Type> prevalence_preg(invlogit(mu_preg));

    for (int i = 0; i < anc_obs_idx.size(); i++) {                                   // index to exclude districts with 0|NA anc attendance (denom)
       val -= dbinom(non_pregnant[preg_obs_idx[i]], all_women[preg_obs_idx[i]], prevalence_preg[preg_obs_idx[i]], true);
    }


  REPORT(beta_hiv);
  REPORT(beta_anc);
  REPORT(beta_preg);

  REPORT(prevalence_hiv);
  REPORT(prevalence_anc);
  REPORT(prevalence_preg);

  REPORT(b_hiv);
  REPORT(b_anc);
  REPORT(u_raw_space_race_hiv);
  REPORT(u_raw_space_race_anc);

  // ADREPORT(b_hiv);
  // ADREPORT(b_anc);
  return val;
}
