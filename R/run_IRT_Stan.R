
### run IRT Stan models


run_IRT_Stan <- function(data,
                         model = c("2PL","3PL","GRM"),
                         D = 1,
                         it_d = NULL,
                         burnin = 3000,
                         sample = 1000,
                         chains=3,
                         cores = parallel::detectCores(),
                         log_lik = F,
                         theta = F,
                         convergence_loop = T,
                         Rhat_criteria = 1.05){


  model <- tolower(model)

  if(model == "grm"){
    if(sum(data == 0) > 0){stop("For GRM the ordinal data most start at 1 instead of 0. \n
                                Recode the data to have no zeros")}
  }

  ### select desired model
  if(D == 1){ mod <- paste0("irt_", model, "_1d" ) }
  if( D > 1){ mod <- paste0("irt_", model, "_md" ) }

  ### get parameters to extract
  pars <- params[[mod]]

  if(log_lik){pars <- c(pars, "log_lik")}
  if(theta){pars <- c(pars, "theta")}

  dat_stan <- data_prep(data,
                        model = model,
                        D = D,
                        it_d = it_d)

  if(convergence_loop){

    iters <- burnin
    keep <- sample
    rhat <- 20
    while(rhat > Rhat_criteria){
      iters <- iters + sample
      burn <- iters - keep

      out <- stan(data=dat_stan,
                  pars=pars,
                  model_code=mod,
                  chains=chains,
                  iter=iters,
                  warmup=burn,
                  cores=cores)

      ss <- summary(OUT_grm)$summary
      rhat <- max(ss[,"Rhat"], na.rm=T)
      print(paste0("Rhat=",rhat))
    }

  }else{
    iter <- burnin + sample

    out <- stan(data=dat_stan,
                pars=pars,
                model_code=mod,
                chains=chains,
                iter=iter,
                warmup=burnin,
                cores=cores)
  }

  return(out)

}




