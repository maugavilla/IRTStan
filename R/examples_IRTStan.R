

library(IRTStan)
library(rstan)
library(loo)

lsat <- ltm::LSAT
head(lsat)
dim(lsat)


#### 1 Dimension 2PL
ex_2pl_1d <- run_IRT_Stan(data=lsat,
                         model = c("2PL"),
                         D = 1,
                         burnin = 3000,
                         sample = 1000,
                         chains=3,
                         cores = parallel::detectCores(),
                         log_lik = T,
                         theta = T,
                         convergence_loop = T,
                         Rhat_criteria = 1.05)
ex_2pl_1d

print(ex_2pl_1d, pars=c("a","b"))
print(ex_2pl_1d, pars=c("theta"))

ex_ll_2pl <- extract_log_lik(ex_2pl_1d, merge_chains = FALSE)
loo_2pl <- loo(ex_ll_2pl, r_eff = relative_eff(exp(ex_ll_2pl)) )
loo_2pl



#### 1 Dimension 3PL
ex_3pl_1d <- run_IRT_Stan(data=lsat,
                          model = c("3PL"),
                          D = 1,
                          burnin = 3000,
                          sample = 1000,
                          chains=3,
                          cores = parallel::detectCores(),
                          log_lik = T,
                          theta = T,
                          convergence_loop = T,
                          Rhat_criteria = 1.05)
ex_3pl_1d

print(ex_3pl_1d, pars=c("a","b","c"))
print(ex_3pl_1d, pars=c("theta"))

ex_ll_3pl <- extract_log_lik(ex_3pl_1d, merge_chains = FALSE)
loo_3pl <- loo(ex_ll_3pl, r_eff = relative_eff(exp(ex_ll_3pl)) )
loo_3pl


### compare LOO: 2PL vs 3PL
loo_compare(loo_2pl, loo_3pl)



#### 3 Dimesional 2PL
hzs <- lavaan::HolzingerSwineford1939

dat_bin <- rockchalk::mvrnorm(n = 500,
                               mu = rep(0, 9),
                               Sigma = cov(hzs[,7:15]) )
psych::describe(dat_bin)
dat_bin <- apply(dat_bin, 2, cut, breaks = 2, labels = FALSE) - 1
apply(dat_bin, 2, table)


ex_2pl_3d <- run_IRT_Stan(data=dat_bin,
                          model = c("2PL"),
                          D = 3,
                          it_d=c(1,1,1,2,2,2,3,3,3),
                          burnin = 3000,
                          sample = 2000,
                          chains=3,
                          cores = parallel::detectCores(),
                          log_lik = T,
                          theta = T,
                          convergence_loop = T,
                          Rhat_criteria = 1.05)
ex_2pl_3d

print(ex_2pl_3d, pars=c("a","b","corr_d"))
print(ex_2pl_3d, pars=c("theta"))

ex_ll_2pl_3d <- extract_log_lik(ex_2pl_3d, merge_chains = FALSE)
loo_2pl_3d <- loo(ex_ll_2pl_3d, r_eff = relative_eff(exp(ex_ll_2pl_3d)) )
loo_2pl_3d



#### 3 Dimesional 3PL

ex_3pl_3d <- run_IRT_Stan(data=dat_bin,
                          model = c("3PL"),
                          D = 3,
                          it_d=c(1,1,1,2,2,2,3,3,3),
                          burnin = 3000,
                          sample = 2000,
                          chains=3,
                          cores = parallel::detectCores(),
                          log_lik = T,
                          theta = T,
                          convergence_loop = T,
                          Rhat_criteria = 1.05)
ex_3pl_3d

print(ex_3pl_3d, pars=c("a","b","c","corr_d"))
print(ex_3pl_3d, pars=c("theta"))

ex_ll_3pl_3d <- extract_log_lik(ex_3pl_3d, merge_chains = FALSE)
loo_3pl_3d <- loo(ex_ll_3pl_3d, r_eff = relative_eff(exp(ex_ll_3pl_3d)) )
loo_3pl_3d





###### 1 Dimension GRM

sci <- ltm::Science
head(sci)
dim(sci)

for(j in 1:ncol(sci)){
  sci[,j] <- car::recode(sci[,j],
                         " 'strongly disagree'=1;'disagree'=2;
                         'agree'=3;'strongly agree'=4 ")
}

apply(sci, 2, table)
apply(sci[,c(1,3,4,7)], 2, table) ## positively worded
apply(sci[,c(2,5,6)], 2, table) ## negatively worded


ex_grm_1d <- run_IRT_Stan(data=sci[,c(1,3,4,7)],
                          model = c("GRM"),
                          D = 1,
                          burnin = 3000,
                          sample = 1000,
                          chains=3,
                          cores = parallel::detectCores(),
                          log_lik = T,
                          theta = T,
                          convergence_loop = T,
                          Rhat_criteria = 1.05)
ex_grm_1d

print(ex_grm_1d, pars=c("alpha","thresholds"))
print(ex_grm_1d, pars=c("theta"))

ex_ll_grm <- extract_log_lik(ex_grm_1d, merge_chains = FALSE)
loo_grm <- loo(ex_ll_grm, r_eff = relative_eff(exp(ex_ll_grm)) )
loo_grm



###### 3 Dimension GRM
hzs <- lavaan::HolzingerSwineford1939

dat_grmd <- rockchalk::mvrnorm(n = 500,
                          mu = rep(0, 9),
                          Sigma = cov(hzs[,7:15]) )
psych::describe(dat_grmd)
dat_grmd <- apply(dat_grmd, 2, cut, breaks = 4, labels = FALSE)
apply(dat_grmd, 2, table)


ex_grm_3d <- run_IRT_Stan(data=dat_grmd,
                          model = c("GRM"),
                          D = 3,
                          it_d=c(1,1,1,2,2,2,3,3,3),
                          burnin = 3000,
                          sample = 2000,
                          chains=3,
                          cores = parallel::detectCores(),
                          log_lik = T,
                          theta = T,
                          convergence_loop = T,
                          Rhat_criteria = 1.05)
ex_grm_3d

print(ex_grm_3d, pars=c("alpha","thresholds","corr_d"))
print(ex_grm_3d, pars=c("theta"))

# monitor(ex_grm_3d)

ex_ll_grm_3d <- extract_log_lik(ex_grm_3d, merge_chains = FALSE)
loo_grm_3d <- loo(ex_ll_grm_3d, r_eff = relative_eff(exp(ex_ll_grm_3d)) )
loo_grm_3d



