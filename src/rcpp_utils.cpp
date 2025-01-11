#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <cmath>
#include <vector>

double hermite_he_function(double arg, int order, double gamma);
double exp_kernel_stat(arma::vec& resid, arma::mat& x_transformed);
double bierens_stat(arma::vec& resid, arma::mat& x_transformed, double tau);
arma::mat tuples(int location, int total_sum, arma::mat tmp_tuple);
arma::mat d_tuple_k_sum(int x_dim, int k);
double psi_he_k(arma::rowvec arg, arma::colvec order, double gamma);
arma::mat Q_dist_fft(arma::vec& eigen_values, double h, int N);

arma::mat sample_psi_he_k(arma::mat& arg, arma::mat& order, double gamma) {
  int n_obs = arg.n_rows;
  int n_order = order.n_cols;
  arma::mat res(n_obs, n_order);
  
  for(int i = 0; i < n_order; ++i) {
    for(int j = 0; j < n_obs; ++j) {
      res(j, i) = psi_he_k(arg.row(j), order.col(i), gamma);
    }
  }
  
  return res;
}

void rescale_cov(arma::mat& cov, arma::rowvec lambda){
  int cov_n_rows = cov.n_rows;
  for(int i = 0; i < cov_n_rows; ++i){
    for(int j = 0; j < cov_n_rows; ++j){
      cov(i,j) = cov(i,j)*sqrt(lambda(i)*lambda(j));
    }
  }
}


double psi_he_k(arma::rowvec arg, arma::colvec order, double gamma) {
  int nx = arg.n_elem;
  double psi = 1.0;

  for(int i = 0; i < nx; ++i) {
    psi *= hermite_he_function(arg(i), order(i), gamma);
  }

  return pow(1.0 - gamma*gamma, (double)nx/4.0)*psi;
}

double hermite_he_function(double arg, int order, double gamma) {
// A C++ translation of the code at https://www.numbercrunch.de/blog/2014/08/calculating-the-hermite-functions/
  if(order == 0){
    return 1.0*exp(-arg*arg*gamma/(2.0*(1.0 + gamma)));
  }

  if(order == 1){
    return arg*exp(-arg*arg*gamma/(2.0*(1.0 + gamma)));
  }

  if(arg == 0.0){
    arg = 0.000001; // the case arg = 0 to be handled separately, but it is not important
  }

  if(order > 1){
    double he_i_2 = 1.0;
    double he_i_1 = arg;
    double he_i, scale, log_scale, sum_log_scale = 0.0;
    for(int i = 2; i < order + 1; ++i){
      he_i = arg*he_i_1/sqrt(i) - sqrt((double)(i - 1)/(double)i)*he_i_2; // recursion for He_n/sqrt(n!)
      log_scale = round(log(fabs(he_i))); //rescale to avoid overflow, note the use of fabs instead of abs
      scale = exp(-log_scale);
      he_i_2 = he_i_1*scale;
      he_i_1 = he_i*scale;
      he_i = he_i*scale;
      sum_log_scale += log_scale;
    }
    return he_i*exp(-arg*arg*gamma/(2.0*(1.0 + gamma)) + sum_log_scale); //return the sum of all corrections to the exponent to prevent underflow
  }
  return 0.0; // if n meets none of the conditions
}


arma::mat Q_dist_fft(arma::vec& eigen_values, double h, int N){
  double max_ev = eigen_values(0);
  arma::vec eigen_values_scaled = eigen_values / max_ev;
  double s = 1.0 / (h*N);
  arma::vec t1 = arma::linspace<arma::vec>(1, N, N);
  arma::vec t2 = 2.0 * M_PI * (t1 - 1.0 - N / 2.0)*s;

  arma::cx_vec Q_chf(N);
  arma::vec denom_re(eigen_values_scaled.n_elem, arma::fill::ones);
  
  for(int i = 0; i < N; ++i){
    arma::vec denom_im = -2.0 * eigen_values_scaled * t2(i);
    arma::cx_vec Q_chf_denom_i(denom_re, denom_im);
    Q_chf(i) = arma::prod(1.0 / arma::sqrt(Q_chf_denom_i));
  }

  arma::vec ones_vec(N, arma::fill::ones);
  arma::cx_vec to_fft = arma::pow(-ones_vec, t1 - 1.0) % Q_chf;
  arma::vec scaler = s * arma::pow(-ones_vec, t1 - 1.0 - N / 2.0);
  arma::vec pdf = arma::real(arma::fft(to_fft)) % scaler / max_ev;

  arma::vec arg = (t1 - 1.0 - N / 2.0) * h * max_ev;
  pdf.shed_rows(arma::find(arg < 0.0));
  arg.shed_rows(arma::find(arg < 0.0));

  arma::vec cdf = arma::cumsum(arma::abs(pdf))/arma::sum(arma::abs(pdf));
  // organize the output
  arma::mat distrib(arg.n_elem, 3);
  distrib.col(0) = arg;
  distrib.col(1) = pdf;
  distrib.col(2) = cdf;
  return distrib;
}

// [[Rcpp::export]]
Rcpp::List lin_test_mc_cpp(arma::mat& x, arma::vec& y, int n_boot, int k_grid, Rcpp::String type, double bierens_tau){
  arma::vec y_centered = y - arma::mean(y);
  int x_dim = x.n_cols;
  int n_obs = x.n_rows;

  arma::mat x_centered(n_obs, x_dim);
  for(int i = 0; i < x_dim; ++i){
    x_centered.col(i) = x.col(i) - arma::mean(x.col(i));
  }

  arma::vec resid = y_centered - x_centered*arma::inv(x_centered.t()*x_centered)*x_centered.t()*y_centered;
  arma::mat C_half = arma::chol(arma::cov(x_centered));
  arma::mat x_std = x_centered*arma::inv(C_half);
  arma::mat u;
  u.set_size(x_dim, k_grid);
  double stat;
  
  if(type == "Bierens"){
    x_std = arma::atan(x_std);
    u.fill(arma::fill::randu);
    u = u * 2.0 * bierens_tau - bierens_tau;
    stat = bierens_stat(resid, x_std, bierens_tau);
  }else{
    u.fill(arma::fill::randn);
    stat = exp_kernel_stat(resid, x_std);
  }

  arma::vec vec_zeros(n_obs, arma::fill::zeros);
  arma::cx_mat rhs_u(n_obs, k_grid);
  arma::vec stat_sample(n_boot);

  arma::mat scale_x_inv = arma::inv(x_centered.t()*x_centered / n_obs);
  for(int i = 0; i < k_grid; ++i){
    arma::cx_vec dot_prod(vec_zeros, x_std * u.col(i));
    arma::cx_vec exp_ux = arma::exp(dot_prod);

    arma::cx_vec mean_x_exp_ux(x_dim);
    for(int j = 0; j < x_dim; ++j){
      mean_x_exp_ux(j) = arma::mean(x_centered.col(j) % exp_ux);
    }

    rhs_u.col(i) = resid % (exp_ux - x_centered * scale_x_inv * mean_x_exp_ux)/sqrt(n_obs);
  }

  for(int i = 0; i < n_boot; ++i){
    arma::rowvec delta(n_obs, arma::fill::randn);
    arma::cx_rowvec W_delta_u = delta * rhs_u;
    W_delta_u -= arma::mean(W_delta_u);
    stat_sample(i) = arma::mean(arma::pow(arma::abs(W_delta_u), 2));
  }

  return Rcpp::List::create(Rcpp::Named("stat") = stat, 
                            Rcpp::Named("stat_boot") = stat_sample, 
                            Rcpp::Named("rhs_u") = rhs_u, 
                            Rcpp::Named("resid") = resid);
}

double bierens_stat(arma::vec& resid, arma::mat& x_transformed, double tau){
  double int_sum = 0.0;
  int n_obs = x_transformed.n_rows;
  int x_dim = x_transformed.n_cols;

  for(int i = 0; i < n_obs - 1; ++i){
    for(int j = i + 1; j < n_obs; ++j){

      double prod = 1.0;
      arma::rowvec delta_x = x_transformed.row(i) - x_transformed.row(j);

      for(int k = 0; k < x_dim; ++k){
        if(delta_x(k) != 0.0){
          prod *= sin(tau*delta_x(k))/(tau*delta_x(k));
        }
      }

      int_sum += resid(i)*resid(j)*prod;
    }
  }
  int_sum *= 2.0/n_obs;
  int_sum += sum(resid % resid)/n_obs;
  return int_sum;
}

double exp_kernel_stat(arma::vec& resid, arma::mat& x_transformed){
  double stat = 0.0;
  int n_obs = x_transformed.n_rows;
  for(int i = 0; i < n_obs - 1; ++i){

    for(int j = i + 1; j < n_obs; ++j){
      arma::rowvec delta_x = x_transformed.row(i) - x_transformed.row(j);
      stat += resid(i)*resid(j)*exp(-sum(delta_x % delta_x)/2.0);
    }

  }
  stat *= 2.0/n_obs;
  stat += sum(resid % resid)/n_obs;
  return stat;
}

// [[Rcpp::export]]
Rcpp::List lin_test_cpp(arma::mat& x, arma::vec& y, int k, double gamma){
  arma::vec y_centered = y - arma::mean(y);
  int x_dim = x.n_cols;
  int n_obs = x.n_rows;

  arma::mat x_centered(n_obs, x_dim);
  for(int i = 0; i < x_dim; ++i){
    x_centered.col(i) = x.col(i) - arma::mean(x.col(i));
  }

  arma::vec resid = y_centered - x_centered*arma::inv(x_centered.t()*x_centered)*x_centered.t()*y_centered;
  arma::mat C_half = arma::chol(arma::cov(x_centered));
  arma::mat x_std = x_centered*arma::inv(C_half);

  double stat = exp_kernel_stat(resid, x_std);
  double sig_res = arma::stddev(resid);

  arma::mat poly_order = d_tuple_k_sum(x_dim, k);
  arma::rowvec lambdas(poly_order.n_cols);
  lambdas.fill(gamma);
  lambdas = arma::pow(lambdas, arma::sum(poly_order));

  arma::mat psi_x = sample_psi_he_k(x_std, poly_order, gamma);
  arma::mat scale_x_inv = arma::inv(x_centered.t()*x_centered / n_obs);
  arma::mat mean_xpsi_k = x_centered.t() * psi_x / n_obs;
  arma::mat f_val = psi_x - arma::trans(mean_xpsi_k.t() * scale_x_inv * x_centered.t());

  int C_rank_upper;
  if((int)lambdas.n_elem <= n_obs){
    C_rank_upper = (int)lambdas.n_elem;
  }else{
    C_rank_upper = n_obs;
  }

  // compute the covariance and rescale it
  arma::mat C = sig_res * sig_res * arma::cov(f_val.cols(0, C_rank_upper - 1));
  arma::rowvec lambdas_up_to_rank = lambdas.cols(0, C_rank_upper - 1);
  rescale_cov(C, lambdas_up_to_rank);

  // compute the eigenvalues and remove the unnecessary ones
  arma::vec eigen_values = arma::eig_sym(C);
  eigen_values = arma::sort(eigen_values, "descend");
  int C_rank = arma::rank(C);

  arma::vec chf_kappa;
  if((int)lambdas.n_elem > n_obs){
    eigen_values.shed_rows(C_rank, eigen_values.n_elem - 1);
    arma::rowvec eigen_tail = arma::var(f_val.cols(C_rank, f_val.n_cols - 1));
    arma::rowvec lambdas_tail = lambdas.cols(C_rank, lambdas.n_elem - 1);
    eigen_tail = sig_res * sig_res * (eigen_tail % lambdas_tail);
    chf_kappa = arma::join_vert(eigen_values, eigen_tail.t());
  }else{
    chf_kappa = eigen_values;
  }

  double h = 0.005;
  int K = 2 + log(2.0*sum(chf_kappa)/(h*chf_kappa(0)))/log(2.0);
  if(K < 13){
    K = 13;
  }

  arma::mat distrib = Q_dist_fft(chf_kappa, h, pow(2, K));
//  arma::uvec indices = arma::find(distrib.col(0) > stat, 1);
//  double prob = distrib(indices[0] - 1, 2);
  double prob;
  if(distrib.col(0).max() <= stat){
    prob = 1.0;
  }else{
    arma::uvec indices = arma::find(distrib.col(0) > stat, 1);
    prob = distrib(indices(0) - 1, 2);
  }
  
  return Rcpp::List::create(Rcpp::Named("stat") = stat, 
                            Rcpp::Named("dist") = distrib, 
                            Rcpp::Named("prob") = prob, 
                            Rcpp::Named("p_value") = 1.0 - prob);
}


arma::mat tuples(int location, int total_sum, arma::mat tmp_tuple){

  if(location == 0){
    tmp_tuple(0, 0) = total_sum;
    return tmp_tuple;
  }

  arma::mat res;
  arma::mat out;
  for(int i = 0; i <= total_sum; ++i){
    tmp_tuple(location, 0) = i;
    out = tuples(location - 1, total_sum - i, tmp_tuple);
    res = arma::join_horiz(res, out);
  }
  return res;
}

arma::mat d_tuple_k_sum(int x_dim, int k){
  arma::mat res(x_dim, 1, arma::fill::zeros);
  arma::mat initialize = res;
  for(int i = 1; i <= k; ++i){
    arma::mat v_i = tuples(x_dim - 1, i, initialize);
    res = arma::join_horiz(res, v_i);
  }
  return res;
}