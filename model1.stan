data{
  int<lower=1> N; // observations
  int<lower=1> A; // ages
  int<lower=1> L; // countries
  int<lower=1> T; // years
  int<lower=1> K; // spline basis functions
  int<lower=1> Q; // K-2
  int<lower=1> P; // covariates
  
  matrix[T,K] B; // B-spline basis matrix
  matrix[K,Q] Pmat; // projection matrix: D'(DD)'^-1
  
  vector[N] logit_qx; // mortality
  matrix[N,P] X; // covariates
  
  array[N] int<lower=1,upper=A> xid; // age index
  array[N] int<lower=1,upper=T> tid; // year index
  array[N] int<lower=1,upper=L> lid; // country index
  
  vector[P] mu_beta; // infromative prior on beta
  vector<lower=0>[P] sigma_beta;
  
  real<lower=0> scale_global; // tau_0 scale
  real<lower=1> nu_global; // df for half-t prior on tau
  real<lower=1> nu_local; // df for half-t prior on lambda
  real<lower=0> slab_scale; // slab scale for RHS
  real<lower=0> slab_df; // df for RHS
}
parameters{
  matrix[A,L] alpha0; // spline level
  matrix[A,L] alpha1; // spline slope
  array[A,L] vector[Q] eps_raw; // second-order difference
  
  real chi; // mean sigma
  real<lower=0> phi_sigma; // sd sigma
  matrix<lower=0>[A,L] sigma_xl;
  
  vector[P] beta;
  real<lower=0> sigma_q;
  
  matrix[T,L] z;
  matrix<lower=0>[T,L] lambda_d;
  real<lower=0> tau;
  real<lower=0> caux;
}
transformed parameters{
  real<lower=0> c;
  matrix<lower=0>[T,L] lambda_tilde;
  matrix[T,L] delta;
  array[A,L] vector[Q] eps;
  array[A,L] vector[K] lambda;
  vector[K] linear_trend;
  vector[N] mu_logit;
  
  c = slab_scale * sqrt(caux);
  for(t in 1:T) for(l in 1:L){
    lambda_tilde[t,l] = sqrt(c^2 * square(lambda_d[t,l])/(c^2 + tau^2 * square(lambda_d[t,l])));
    delta[t,l] = z[t,l] * lambda_tilde[t,l] * tau;
  }
  
  for(a in 1:A) for(l in 1:L){
    eps[a,l] = eps_raw[a,l] * sigma_xl[a,l];
    for(k in 1:K){
      linear_trend[k] = alpha0[a,l] + alpha1[a,l]*(k-K/2.0);
    }
    lambda[a,l] = linear_trend + Pmat * eps[a,l];
  }
  for(n in 1:N){
    mu_logit[n] = dot_product(B[tid[n]],lambda[xid[n], lid[n]]) + dot_product(X[n], beta) + delta[tid[n],lid[n]];
  }
}
model{
  to_vector(alpha0) ~ normal(0,5);
  to_vector(alpha1) ~ normal(0,1);
  chi ~ normal(0,1);
  phi_sigma ~ normal(0,1);
  for(a in 1:A)for(l in 1:L){
    sigma_xl[a,l] ~ lognormal(chi, phi_sigma);
    eps_raw[a,l] ~ normal(0, 1);
  }
  beta ~ normal(mu_beta, sigma_beta);
  sigma_q ~ normal(0,0.5);
  
  to_vector(z) ~ normal(0,1);
  to_vector(lambda_d) ~ student_t(nu_local, 0, 1);
  tau ~ student_t(nu_global, 0, scale_global*sigma_q);
  caux ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df);
  
  logit_qx ~ normal(mu_logit, sigma_q);
}
generated quantities{
  vector[N] qx_rep;
  vector[N] log_lik;
  for(n in 1:N){
    log_lik[n] = normal_lpdf(logit_qx[n] | mu_logit[n], sigma_q);
    qx_rep[n] = inv_logit(normal_rng(mu_logit[n], sigma_q));
  }
}

