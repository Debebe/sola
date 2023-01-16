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


  // Function for detecting NAs
  template<class Type>
  bool isNA(Type x){
    return R_IsNA(asDouble(x));
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
  DATA_SPARSE_MATRIX(Z_age);
  DATA_SPARSE_MATRIX(Z_age_sex);             // age sex interactrion // n_obs by  n_age
  DATA_SPARSE_MATRIX(Z_space_sex);          // space sex interactrion // n_obs by  n_space
  DATA_SPARSE_MATRIX(Z_time);
  DATA_SPARSE_MATRIX(Z_time_sex);         // time sex interactrion // n_obs by  n_age
  DATA_SCALAR(rankdef_Q_space);


  Type val(0);
  PARAMETER_VECTOR(beta);

   /////////////////////////////////////////////
   /////////   space-BYM2    ///////////////////
   /////////////////////////////////////////////

  PARAMETER(log_sigma_space);                                                   // marginal standard deviation
  PARAMETER(logit_phi_space);
  PARAMETER_VECTOR(b);                                                         // combined spatial effect
  PARAMETER_VECTOR(u);                                                         // spatially correlated component

  Type sigma_space(exp(log_sigma_space));
  val -= dnorm(sigma_space, Type(0.0), Type(1.0), true) + log_sigma_space;      // log_sigma: log-absolute Jacobian of exp(log_sigma)

  Type phi_space(invlogit(logit_phi_space));
  val -= log(phi_space) +  log(1 - phi_space);                                 // change of variables: logit_phi -> phi
  val -= dbeta(phi_space, Type(0.5), Type(0.5), true);

  val -= dnorm(sum(b), Type(0.0), Type(0.001) * b.size(), true);               // soft sum-to-zero constraint
  val -= bym2_conditional_lpdf(b, u, sigma_space, phi_space, Q_space);


   /////////////////////////////////////////////////////////
   /////////   AR1-age random effect    ///////////////////
   ///////////////////////////////////////////////////////

  PARAMETER(logit_phi_age);                                                       // correlation between ages
  PARAMETER(log_sigma_age);
  PARAMETER_VECTOR(u_age);

  val -= dnorm(logit_phi_age, Type(0.0), Type(2.582), true);                      // INLA default
  Type phi_age(2.0 * invlogit(logit_phi_age) - 1.0);
  Type sigma_age(exp(log_sigma_age));
  val -= dnorm(sigma_age, Type(0.0), Type(2.5), true) + log_sigma_age;
  if(u_age.size() > 0)
    val += SCALE(AR1(phi_age), sigma_age)(u_age);

  // Sum to zero constraint
  val -= dnorm(u_age.sum(), Type(0), Type(0.001) * u_age.size(), TRUE);

   ///////////////////////////////////////////////////////
   /////////   space-age -interaction   //////////////////
   ///////////////////////////////////////////////////////

   DATA_SPARSE_MATRIX(Z_space_age);
   PARAMETER(log_sigma_space_age);
   PARAMETER(logit_rho_space_age);
   PARAMETER_ARRAY(u_raw_space_age);

   // val -= dlgamma(log_sigma_space_age, Type(0.001), Type(1.0 / 0.001), true);
   // Type sigma_space_age(exp(-0.5 * log_sigma_space_age));
   // val -= dnorm(logit_rho_space_age, Type(0.0), Type(1.0 / sqrt(0.15)), true);
   // Type rho_age(2 * exp(logit_rho_space_age)/(1 + exp(logit_rho_space_age)) - 1);
   // vector<Type> u_space_age(u_raw_space_age * sigma_space_age);
   // val += SEPARABLE(AR1(rho_age), GMRF(Q_space))(u_raw_space_age);

    Type sigma_space_age(exp(log_sigma_space_age));
    val -= dnorm(sigma_space_age, Type(0.0), Type(2.5), true) + log_sigma_space_age;
    val -= dnorm(logit_rho_space_age, Type(0.0), Type(2.582), true);
    Type rho_space_age(2.0 * invlogit(logit_rho_space_age) - 1.0);
    vector<Type> u_space_age(u_raw_space_age * sigma_space_age);
    val += SEPARABLE(AR1(rho_space_age), GMRF(Q_space))(u_raw_space_age);

   // Adjust normalising constant for rank deficience of R_space
    Type log_det_Qar1_space_age((u_raw_space_age.cols() - 1) * log(1 - rho_space_age * rho_space_age));
    val -= rankdef_Q_space * 0.5 * (log_det_Qar1_space_age - log(2 * PI));

   // sum-to-zero- constraint on interaction term

   for (int i = 0; i < u_raw_space_age.cols(); i++) {
      val -= dnorm(u_raw_space_age.col(i).sum(), Type(0), Type(0.001) * u_raw_space_age.rows(), true); //gives NA/NAN warinng
   }

   ///////////////////////////////////////////////////////
   /////////   Age(AR1)-sex(fixed) -interaction   ////////
   ///////////////////////////////////////////////////////

   PARAMETER(logit_phi_age_sex);
   PARAMETER(log_sigma_age_sex);
   PARAMETER_VECTOR(u_age_sex);                                  //random effect for the interaction

   val -= dnorm(logit_phi_age_sex, Type(0.0), Type(2.582), true);  // INLA default
   Type phi_age_sex(2.0 * invlogit(logit_phi_age_sex) - 1.0);

   Type sigma_age_sex(exp(log_sigma_age_sex));
   val -= dnorm(sigma_age_sex, Type(0.0), Type(2.5), true) + log_sigma_age_sex;

   if(u_age_sex.size() > 0)
    val += SCALE(AR1(phi_age_sex), sigma_age_sex)(u_age_sex);

   // Sum to zero constraint
   val -= dnorm(u_age_sex.sum(), Type(0), Type(0.001) * u_age_sex.size(), TRUE);

     ///////////////////////////////////////////////////////
    /////////   ICAR-sex(fixed) -interaction   ////////////
   ///////////////////////////////////////////////////////

   PARAMETER(log_sigma_space_sex);
   PARAMETER_VECTOR(u_space_sex);

   Type sigma_space_sex(exp(log_sigma_space_sex));
   val -= dnorm(sigma_space_sex, Type(0.0), Type(2.5), true) + log_sigma_space_sex;

   if(u_space_sex.size() > 0)
    val += SCALE(GMRF(Q_space), sigma_space_sex)(u_space_sex);

    /*pay attention here*/
    // sum-to-zero- constraint on interaction term
   for (int i = 0; i < u_space_sex.cols(); i++) {
     val -= dnorm(u_space_sex.col(i).sum(), Type(0), Type(0.001) * u_space_sex.rows(), true);
  }

    ////////////////////////////////////////////////////////
    /////////   year - AR1       //////////////////////////
   ///////////////////////////////////////////////////////

   PARAMETER(logit_phi_time);                                                       // correlation between ages
   PARAMETER(log_sigma_time);
   PARAMETER_VECTOR(u_time);

   val -= dnorm(logit_phi_time, Type(0.0), Type(2.582), true);                      // INLA default
   Type phi_time(2.0 * invlogit(logit_phi_time) - 1.0);
   Type sigma_time(exp(log_sigma_time));
   val -= dnorm(sigma_time, Type(0.0), Type(2.5), true) + log_sigma_time;
   if(u_time.size() > 0)
    val += SCALE(AR1(phi_time), sigma_time)(u_time);

   val -= dnorm(u_time.sum(), Type(0), Type(0.001) * u_time.size(), TRUE);


   ////////////////////////////////////////////////////////
   /////////   year(AR1)-sex interaction       ///////////
   ///////////////////////////////////////////////////////

   PARAMETER(logit_phi_time_sex);
   PARAMETER(log_sigma_time_sex);
   PARAMETER_VECTOR(u_time_sex);                                                  //random effect for the interaction

   val -= dnorm(logit_phi_time_sex, Type(0.0), Type(2.582), true);                // INLA default
   Type phi_time_sex(2.0 * invlogit(logit_phi_time_sex) - 1.0);

   Type sigma_time_sex(exp(log_sigma_time_sex));
   val -= dnorm(sigma_time_sex, Type(0.0), Type(2.5), true) + log_sigma_time_sex;

   if(u_time_sex.size() > 0)
    val += SCALE(AR1(phi_time_sex), sigma_time_sex)(u_time_sex);

   val -= dnorm(u_time_sex.sum(), Type(0), Type(0.001) * u_time_sex.size(), TRUE);


    ////////////////////////////////////////////////////////
    /////////   year X age (AR1 X AR1)          ///////////
   ///////////////////////////////////////////////////////

   DATA_SPARSE_MATRIX(Z_time_age);       // time age interactrion //
   PARAMETER_ARRAY(u_raw_time_age);
   PARAMETER(log_sigma_time_age);
   Type sigma_time_age (exp(log_sigma_time_age));
   val -= dnorm(sigma_time_age,  Type(0.0), Type(2.5), true) + log_sigma_time_age;

   PARAMETER(logit_phi_timex);
   PARAMETER(logit_phi_agex);

   val -= dnorm(logit_phi_timex, Type(0.0), Type(2.582), TRUE);
   val -= dnorm(logit_phi_agex,  Type(0.0), Type(2.582), TRUE);

   Type phi_timex (2.0 * invlogit(logit_phi_timex) - 1.0);
   Type phi_agex (2.0 * invlogit(logit_phi_agex) - 1.0);


   vector<Type> u_time_age(u_raw_time_age * sigma_time_age);

   if(u_raw_time_age.size() > 0)
   val += SEPARABLE(AR1(phi_agex), AR1(phi_timex))(u_raw_time_age);
   
    // sum-to-0 constraint
   for (int i = 0; i < u_raw_time_age.cols(); i++) {
      val -= dnorm(u_raw_time_age.col(i).sum(), Type(0), Type(0.001) * u_raw_time_age.col(i).size(), true);
    }
     ///////////////////////////////////////////////////////
    /////////   space-time interaction          ///////////
   ///////////////////////////////////////////////////////

   DATA_SPARSE_MATRIX(Z_space_time);
   PARAMETER(log_sigma_space_time);
   PARAMETER(logit_rho_space_time);
   PARAMETER_ARRAY(u_raw_space_time);

   // val -= dlgamma(log_sigma_space_time, Type(0.001), Type(1.0/0.001), true);
   // Type sigma_space_time(exp(-0.5 * log_sigma_space_time));
   // val -= dnorm(logit_rho_space_time, Type(0.0), Type(1.0/sqrt(0.15)), true);
   // Type rho_time(2 * exp(logit_rho_space_time)/(1 + exp(logit_rho_space_time)) - 1);

   Type sigma_space_time(exp(log_sigma_space_time));
   val -= dnorm(sigma_space_time, Type(0.0), Type(2.5), true) + log_sigma_space_time;
   val -= dnorm(logit_rho_space_time, Type(0.0), Type(2.582), true);
   Type rho_space_time(2.0 * invlogit(logit_rho_space_time) - 1.0);
   vector<Type> u_space_time(u_raw_space_time * sigma_space_time);
   val += SEPARABLE(AR1(rho_space_time), GMRF(Q_space))(u_raw_space_time);

   // Adjust normalising constant for rank deficience of R_space
   Type log_det_Qar1_space_time((u_raw_space_time.cols() - 1) * log(1 - rho_space_time * rho_space_time));
   val -= rankdef_Q_space * 0.5 * (log_det_Qar1_space_time - log(2 * PI));

   //sum-to-zero- constraint on interaction term

   for (int i = 0; i < u_raw_space_time.cols(); i++) {
     val -= dnorm(u_raw_space_time.col(i).sum(), Type(0), Type(0.001) * u_raw_space_time.rows(), true);
   }

   /////////////////////////////////////////////
   /// Prior on the three way random effects ///
   /////////////////////////////////////////////

   DATA_SPARSE_MATRIX(R_space);                 // diagonal matrix for IID
   DATA_SPARSE_MATRIX(Z_age_space_time);        // Design matrix for the 3 way interaction
   PARAMETER_ARRAY(u_raw_age_space_time);
   PARAMETER(log_sigma_age_space_time);
   PARAMETER(logit_rho2_age);
   PARAMETER(logit_rho2_time);

   Type sigma_age_space_time(exp(log_sigma_age_space_time));
   val -= dnorm(sigma_age_space_time, Type(0.0), Type(2.5), true) + log_sigma_age_space_time;
   //val -= dgamma(sigma_age_space_time, Type(1), Type(2000), true);
   //Type rho2_age(exp(logit_rho2_age)/(1+exp(logit_rho2_age)));
   //val -= log(rho2_age) +  log(1 - rho2_age);
   //val -= dbeta(rho2_age, Type(0.5), Type(0.5), true);
   //Type rho2_time(exp(logit_rho2_time)/(1+exp(logit_rho2_time)));
   //val -= log(rho2_time) +  log(1 - rho2_time);
   //val -= dbeta(rho2_time, Type(0.5), Type(0.5), true);
   Type rho2_age(2.0 * invlogit(logit_rho2_age) - 1.0);
   val -= dnorm(logit_rho2_age, Type(0.0), Type(2.582), true);
   Type rho2_time(2.0 * invlogit(logit_rho2_time) - 1.0);
   val -= dnorm(logit_rho2_time, Type(0.0), Type(2.582), true);

   vector<Type> u_age_space_time(u_raw_age_space_time * sigma_age_space_time);
   val += SEPARABLE(AR1(Type(rho2_age)), SEPARABLE(AR1(Type(rho2_time)), GMRF(R_space)))(u_raw_age_space_time); //IID for space
   //val += SEPARABLE(AR1(Type(rho2_age)), SEPARABLE(AR1(Type(rho2_time)), GMRF(Q_space)))(u_raw_age_space_time); //ICAR for space

   // sum-to-zero constraint

   val -= dnorm(u_raw_age_space_time.sum(), Type(0), Type(0.001) * u_raw_age_space_time.size(), TRUE);


     ///////////////////////////////////////////////////////
    /////////   Space-age-sex                   ///////////
   ///////////////////////////////////////////////////////

   DATA_SPARSE_MATRIX(Z_space_sex_age);
   PARAMETER_ARRAY(u_raw_space_sex_age);
   PARAMETER(log_sigma_space_sex_age);
   PARAMETER(logit_rho3_age);

   Type sigma_space_sex_age(exp(log_sigma_space_sex_age));
   val -= dnorm(sigma_space_sex_age, Type(0.0), Type(2.5), true) + log_sigma_space_sex_age;
   val -= dnorm(logit_rho3_age, Type(0.0), Type(2.582), true);
   Type rho3_age(2 * exp(logit_rho3_age)/(1 + exp(logit_rho3_age)) - 1);

   vector<Type> u_space_sex_age(u_raw_space_sex_age * sigma_space_sex_age);

   if(u_raw_space_sex_age.size() > 0)
   val += SEPARABLE(AR1(rho3_age), GMRF(Q_space))(u_raw_space_sex_age);

   // Adjust normalising constant for rank deficience of R_space
   Type log_det_Qar1_space_age_sex((u_raw_space_sex_age.cols() - 1) * log(1 - rho3_age * rho3_age));
   val -= rankdef_Q_space * 0.5 * (log_det_Qar1_space_age_sex - log(2 * PI));

  //sum-to-zero- constraint on interaction term
   for (int i = 0; i < u_raw_space_sex_age.cols(); i++) {
     val -= dnorm(u_raw_space_sex_age.col(i).sum(), Type(0), Type(0.001) * u_raw_space_sex_age.rows(), true);
   }


     //////////////////////////////////////////////////////
    /////////   Space-time-sex                   /////////
   //////////////////////////////////////////////////////

   DATA_SPARSE_MATRIX(Z_space_sex_time);
   PARAMETER_ARRAY(u_raw_space_sex_time);
   PARAMETER(log_sigma_space_sex_time);
   PARAMETER(logit_rho3_time);

   Type sigma_space_sex_time(exp(log_sigma_space_sex_time));
   val -= dnorm(sigma_space_sex_time, Type(0.0), Type(2.5), true) + log_sigma_space_sex_time;
   val -= dnorm(logit_rho3_time, Type(0.0), Type(2.582), true);
   // val -= dnorm(logit_rho3_time, Type(0.0), Type(1.0 / sqrt(0.15)), true);
   Type rho3_time(2 * exp(logit_rho3_time)/(1 + exp(logit_rho3_time)) - 1);
   vector<Type> u_space_sex_time(u_raw_space_sex_time * sigma_space_sex_time);
   if(u_raw_space_sex_time.size() > 0)
   val += SEPARABLE(AR1(rho3_time), GMRF(Q_space))(u_raw_space_sex_time);

   // Adjust normalising constant for rank deficience of R_space
   Type log_det_Qar1_space_time_sex((u_raw_space_sex_time.cols() - 1) * log(1 - rho3_time * rho3_time));
   val -= rankdef_Q_space * 0.5 * (log_det_Qar1_space_time_sex - log(2 * PI));

  //sum-to-zero- constraint on interaction term
   for (int i = 0; i < u_raw_space_sex_time.cols(); i++) {
    val -= dnorm(u_raw_space_sex_time.col(i).sum(), Type(0), Type(0.001) * u_raw_space_sex_time.rows(), true);

   }

  // model
  vector<Type> mu(X*beta +
                  Z_space *b +
                  Z_age * u_age +

                  Z_age_sex * u_age_sex +
                  Z_space_sex * u_space_sex +
                  Z_space_age * u_space_age +

                  Z_time * u_time +

                  Z_time_sex * u_time_sex +
                  Z_time_age * u_time_age +

                  Z_space_time * u_space_time+
                  Z_age_space_time*u_age_space_time +

                  Z_space_sex_age  * u_space_sex_age +
                  Z_space_sex_time * u_space_sex_time);

  vector <Type> prevalence(invlogit(mu));

// for (int i = 0; i < obs_idx.size(); i++) {                                   // index to complete data
//    val -= dbinom(Y[obs_idx[i]], N[obs_idx[i]], prevalence[obs_idx[i]], true);
// }

   vector <Type> pointwise_ll(Y.size());

  for (int i = 0; i < Y.size(); i++) {
    if (Y[i] != Y[i]) { // f!= f returns true if and only if f is NaN.
      pointwise_ll(i) = Type(0.0);
    }
    else {
      pointwise_ll(i) = dbinom(Y[i], N[i], prevalence[i], true);
    }
  }
  val -= sum(pointwise_ll);

  // This one is also works- but if using the following comment out obs_idx in the data declaration section
  // for( int i=0; i<Y.size(); i++){
  //   if( !isNA(Y(i)) ) val -= dbinom(Y[i], N[i], prevalence[i], true);      //
  // }

  REPORT(beta);
  REPORT(prevalence);
  REPORT(pointwise_ll);

  REPORT(b);
  REPORT(u);
  REPORT(u_raw_space_age);
  REPORT(u_age);
  REPORT(u_age_sex);   
  REPORT(u_time_sex);
  REPORT(u_space_sex);
  REPORT(u_raw_time_age);
  REPORT(u_time);
  REPORT(u_raw_space_time);
  REPORT(u_raw_age_space_time);
  REPORT(u_raw_space_sex_age);
  REPORT(u_raw_space_sex_time);


  return val;
}
