

### Stan models
### this file has the Stan code for the 2Pl, 3PL and GRM models

## add
## Partial credit and nominal later


models_stan <- list(
  ##### 1 dimension, 2PL
  irt_2pl_1d = "
data{
int N; //number of subjects
int n; // number of items
int<lower=0, upper=1> x[N,n]; // binary data matrix
}

parameters{
real<lower=0, upper=5> a[n];
real<lower=-4, upper=4> b[n]; // item difficulties
real theta[N]; // latent factor scores
}

model{

for(i in 1:N){// MIRT formulation
for(j in 1:n){
x[i,j] ~ bernoulli_logit(a[j] * (theta[i] - b[j]) );
}
}

b ~ normal(0,1);  /// difficulty prior
a ~ lognormal(0,1); /// or lognormal(0,2)
theta ~ normal(0,1); // latent factor

}

// calculate the correlation from the cholesky
generated quantities {
vector[n] log_lik[N]; ///// row log likelihood
real dev; /////// deviance
real log_lik0; ///// global log likelihood
vector[N] log_lik_row;

for(i in 1:N){// MIRT formulation
for(j in 1:n){
log_lik[i,j] = bernoulli_logit_lpmf(x[i,j] | a[j] * (theta[i] - b[j]) );
}
}

for(i in 1:N){
log_lik_row[i] = sum(log_lik[i]);
}
log_lik0 = sum(log_lik_row); // global log-likelihood
dev = -2*log_lik0; // model deviance

}
"
,
### multiple Dimensions, 2PL
irt_2pl_mD = "
data{
int N; //number of subjects
int n; // number of items
int<lower=0, upper=1> x[N,n]; // binary data matrix
int D; // number of dimensions/factors
int it_d[n]; // which factor each item goes
int a_sign[n]; // sign for items, first item for each factor is constraint positive
}

parameters{
real<lower=0, upper=1.7> a_norm_pos[n]; // positively constraint discrimitaions
real<upper=1.7> a_norm_rest[n]; // not positivily constraint discrimitaions
real<lower=-4, upper=4> b[n]; // item difficulties
vector[D] theta[N]; // latent factor scores
cholesky_factor_corr[D] L_corr_d; // cholesky correlation between factors
}

transformed parameters{
vector[D] mu_d; // mean of factors
real a_norm[n]; // normalized discriminations, log-discriminations
real a[n]; // discrimination parameter in IRT parameterization

for(j in 1:D){ // fixed factor means
mu_d[j] = 0;}

for(k in 1:n){ // assign the positiveley constrant discrimination to the respective item
if(a_sign[k] > 0){a_norm[k] = a_norm_pos[k];}
if(a_sign[k] == 0){a_norm[k] = a_norm_rest[k];}
}

a = exp(a_norm); // transform the normalized discrimination into IRT

}

model{

for(i in 1:N){// MIRT formulation
for(j in 1:n){
x[i,j] ~ bernoulli_logit(a[j] * (theta[i,it_d[j]] - b[j]) );
}
}

L_corr_d ~ lkj_corr_cholesky(1);
b ~ normal(0,1);  /// difficulty prior
a_norm_pos ~ normal(0,1); // discrimination prior
a_norm_rest ~ normal(0,1);// discrimination prior
theta ~ multi_normal_cholesky(mu_d, L_corr_d); // latent factor multivariate normal distribution

}

// calculate the correlation from the cholesky
generated quantities {
corr_matrix[D] corr_d;
vector[n] log_lik[N]; ///// row log likelihood
real dev; /////// deviance
real log_lik0; ///// global log likelihood
vector[N] log_lik_row;

corr_d = multiply_lower_tri_self_transpose(L_corr_d); // tranform cholesky into correlation matrix

for(i in 1:N){// MIRT formulation
for(j in 1:n){
log_lik[i,j] = bernoulli_logit_lpmf(x[i,j] | a[j] * (theta[i,it_d[j]] - b[j]) );
}
}

for(i in 1:N){
log_lik_row[i] = sum(log_lik[i]);
}
log_lik0 = sum(log_lik_row); // global log-likelihood
dev = -2*log_lik0; // model deviance

}
"
,
##### 1 dimension, 3PL
irt_3pl_1d = "
data{
int N; //number of subjects
int n; // number of items
int<lower=0, upper=1> x[N,n]; // binary data matrix
}

parameters{
real<lower=0, upper=5> a[n];
real<lower=-4, upper=4> b[n]; // item difficulties
real theta[N]; // latent factor scores
real<lower=0, upper=1> c[n];
}

model{

for(i in 1:N){//
for(j in 1:n){
x[i,j] ~ bernoulli( c[j] + (1-c[j]) * inv_logit( a[j] * (theta[i] - b[j]) ) );
}
}

b ~ normal(0,1);  /// difficulty prior
a ~ lognormal(0,1); /// or lognormal(0,2)
c ~ beta(5,23);  /// beta(5,17)
theta ~ normal(0,1); // latent factor

}

// calculate the correlation from the cholesky
generated quantities {
vector[n] log_lik[N]; ///// row log likelihood
real dev; /////// deviance
real log_lik0; ///// global log likelihood
vector[N] log_lik_row;

for(i in 1:N){// MIRT formulation
for(j in 1:n){
log_lik[i,j] = bernoulli_lpmf(x[i,j] | c[j] + (1-c[j]) * inv_logit( a[j] * (theta[i] - b[j]) ) );
}
}

for(i in 1:N){
log_lik_row[i] = sum(log_lik[i]);
}
log_lik0 = sum(log_lik_row); // global log-likelihood
dev = -2*log_lik0; // model deviance

}
"
,
#### multiple D, 3PL
irt_3pl_mD = "
data{
int N; //number of subjects
int n; // number of items
int<lower=0, upper=1> x[N,n]; // binary data matrix
int D; // number of dimensions/factors
int it_d[n]; // which factor each item goes
int a_sign[n]; // sign for items, first item for each factor is constraint positive
}

parameters{
real<lower=0, upper=1.7> a_norm_pos[n]; // positively constraint discrimitaions
real<upper=1.7> a_norm_rest[n]; // not positivily constraint discrimitaions
real<lower=-4, upper=4> b[n]; // item difficulties
vector[D] theta[N]; // latent factor scores
cholesky_factor_corr[D] L_corr_d; // cholesky correlation between factors
real<lower=0, upper=1> c[n];
}

transformed parameters{
vector[D] mu_d; // mean of factors
real a_norm[n]; // normalized discriminations, log-discriminations
real a[n]; // discrimination parameter in IRT parameterization

for(j in 1:D){ // fixed factor means
mu_d[j] = 0;}

for(k in 1:n){ // assign the positiveley constrant discrimination to the respective item
if(a_sign[k] > 0){a_norm[k] = a_norm_pos[k];}
if(a_sign[k] == 0){a_norm[k] = a_norm_rest[k];}
}

a = exp(a_norm); // transform the normalized discrimination into IRT

}

model{

for(i in 1:N){// MIRT formulation
for(j in 1:n){
x[i,j] ~ bernoulli( c[j] + (1-c[j]) * inv_logit( a[j] * (theta[i,it_d[j]] - b[j]) ) );
}
}

L_corr_d ~ lkj_corr_cholesky(1);
b ~ normal(0,1);  /// difficulty prior
c ~ beta(5,23);  /// beta(5,17)
a_norm_pos ~ normal(0,1); // discrimination prior
a_norm_rest ~ normal(0,1);// discrimination prior
theta ~ multi_normal_cholesky(mu_d, L_corr_d); // latent factor multivariate normal distribution

}

// calculate the correlation from the cholesky
generated quantities {
corr_matrix[D] corr_d;
vector[n] log_lik[N]; ///// row log likelihood
real dev; /////// deviance
real log_lik0; ///// global log likelihood
vector[N] log_lik_row;

corr_d = multiply_lower_tri_self_transpose(L_corr_d); // tranform cholesky into correlation matrix

for(i in 1:N){// MIRT formulation
for(j in 1:n){
log_lik[i,j] = bernoulli_lpmf(x[i,j] | c[j] + (1-c[j]) * inv_logit( a[j] * (theta[i,it_d[j]] - b[j]) ) );
}
}

for(i in 1:N){
log_lik_row[i] = sum(log_lik[i]);
}
log_lik0 = sum(log_lik_row); // global log-likelihood
dev = -2*log_lik0; // model deviance

}
"
,
## GRM 1 dimension
grm_1d = "
data{
int n; //number of subjects
int p; // number of items
int K; // number of categories
int x[n,p]; // categorical data matrix
}

parameters{
real<lower=0, upper=5> alpha[p]; // slopes for each variable
ordered[K-1] kappa[p]; // thresholds for each variable
real theta[n]; // factor scores
real mu_kappa; //mean of the prior distribution of category difficulty
real<lower=0> sigma_kappa; //sd of the prior distribution of category difficulty
real<lower=0> sigma_alpha;
}

transformed parameters{
vector[K-1] thresholds[p];

for(i in 1:p){
for(k in 1:(K-1)){
thresholds[i,k] = kappa[i,k]/alpha[i];
}
}

}

model{

for(i in 1:n){
for(j in 1:p){
x[i,j] ~ ordered_logistic(alpha[j] * theta[i], kappa[j]);
}
}

for(i in 1:p){
for(k in 1:(K-1)){
kappa[i,k] ~ normal(mu_kappa, sigma_kappa);
}
}

mu_kappa ~ normal(0,5);
sigma_kappa ~ cauchy(0,5);

alpha ~ lognormal(0,sigma_alpha);
sigma_alpha ~ cauchy(0,5);

theta ~ normal(0, 1);

}

//
generated quantities {
vector[p] log_lik[n]; ///// row log likelihood
real dev; /////// deviance
real log_lik0; ///// global log likelihood
vector[n] log_lik_row;

for(i in 1:n){//
for(j in 1:p){
log_lik[i,j] = ordered_logistic_lpmf(x[i,j] | alpha[j] * theta[i], kappa[j] );
}
}

for(i in 1:n){
log_lik_row[i] = sum(log_lik[i]);
}
log_lik0 = sum(log_lik_row); // global log-likelihood
dev = -2*log_lik0; // model deviance

}
"
,
## GMR multiple Dim
grm_stan_md = "
data{
int n; //number of subjects
int p; // number of items
int K; // number of categories
int x[n,p]; // categorical data matrix
int D; // number of dimensions/factors
int it_d[p]; // which factor each item goese
}

parameters{
real<lower=0, upper=5> alpha[p]; // slopes for each variable
ordered[K-1] kappa[p]; // thresholds for each variable
vector[D] theta[n]; // factor scores
real mu_kappa; //mean of the prior distribution of category difficulty
real<lower=0> sigma_kappa; //sd of the prior distribution of category difficulty
real<lower=0> sigma_alpha;
cholesky_factor_corr[D] L_corr_d; // cholesky correlation between factors
}

transformed parameters{
vector[D] mu_d; // mean of factors
vector[K-1] thresholds[p];


for(j in 1:D){ // fixed factor means
mu_d[j] = 0;}

for(i in 1:p){
for(k in 1:(K-1)){
thresholds[i,k] = kappa[i,k]/alpha[i];
}
}

}

model{

for(i in 1:n){
for(j in 1:p){
x[i,j] ~ ordered_logistic(alpha[j] * theta[i,it_d[j]] , kappa[j]);
}
}

for(i in 1:p){
for(k in 1:(K-1)){
kappa[i,k] ~ normal(mu_kappa, sigma_kappa);
}
}

mu_kappa ~ normal(0,5);
sigma_kappa ~ cauchy(0,5);

alpha ~ lognormal(0,sigma_alpha);
sigma_alpha ~ cauchy(0,5);


L_corr_d ~ lkj_corr_cholesky(1);
theta ~ multi_normal_cholesky(mu_d, L_corr_d); // latent factor multivariate normal distribution

}

//
generated quantities {
corr_matrix[D] corr_d;
vector[p] log_lik[n]; ///// row log likelihood
real dev; /////// deviance
real log_lik0; ///// global log likelihood
vector[n] log_lik_row;

corr_d = multiply_lower_tri_self_transpose(L_corr_d); // tranform cholesky into correlation matrix

for(i in 1:n){//
for(j in 1:p){
log_lik[i,j] = ordered_logistic_lpmf(x[i,j] | alpha[j] * theta[i,it_d[j]] , kappa[j] );
}
}

for(i in 1:n){
log_lik_row[i] = sum(log_lik[i]);
}
log_lik0 = sum(log_lik_row); // global log-likelihood
dev = -2*log_lik0; // model deviance

}
"
)

#models_stan
#names(models_stan)



params <- list(
  irt_2pl_1d = c("a","b","log_lik0","dev"),
  irt_2pl_mD = c("a","b","corr_d","log_lik0","dev"),
  irt_3pl_1d = c("a","b","c","log_lik0","dev"),
  irt_3pl_mD = c("a","b","c","corr_d","log_lik0","dev"),
  grm_1d = c("alpha","thresholds","log_lik0","dev"),
  grm_stan_md = c("alpha","thresholds","corr_d","log_lik0","dev")
)


