

##### function for data preparation


data_prep <- function(data,
                      model = c("2PL","3PL","GRM"),
                      D = 1,
                      it_d = NULL){

  model <- tolower(model)

  if(model %in% c("2pl","3pl") & D == 1){
    N <- nrow(data) ## number of subjects
    n <- ncol(data) ## number of items

    dat <- list(N=N,n=n, x=data)
  }

  if(model %in% c("2pl","3pl") & D > 1){
    N <- nrow(data) ## number of subjects
    n <- ncol(data) ## number of items

    its <- 1:length(it_d)
    loc_first <- NULL
    uqs <- unique(it_d)
    for(i in 1:length(uqs)){
      tt <- its[it_d == uqs[i]]
      tt_f <- tt[1]
      loc_first[[i]] <- tt_f
    }

    a_sign <- as.numeric(its %in% loc_first)

    dat <- list(N=N,n=n, x=data, D=D, it_d=it_d, a_sign=a_sign)
  }

  if(model %in% c("grm") & D == 1){
    n <- nrow(data) ## number of subjects
    p <- ncol(data) ## number of variables
    cats <- length(table(data[,1]))

    dat <- list(n=n, p=p, x=data, K = cats)
  }

  if(model %in% c("grm") & D > 1){
    n <- nrow(data) ## number of subjects
    p <- ncol(data) ## number of variables
    cats <- length(table(data[,1]))

    dat <- list(n=n, p=p, x=data, K = cats, D=D, it_d=it_d)
  }


  return(dat)

}

