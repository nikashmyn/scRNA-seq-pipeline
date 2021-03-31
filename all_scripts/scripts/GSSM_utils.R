

getMyGSSM <- function(mv1) {
  require(KFAS)
  #mv1 <- adt[,c("seqnames","160714_A3")] #seqnames == "chr6"
  mv1 <- ts(mv1)
  Zt <- matrix(c(1, 0), 1, 2)
  Ht <- matrix(NA)
  Tt <- matrix(c(1, 0, 1, 1), 2, 2)
  Rt <- matrix(c(1, 0), 2, 1)
  Qt <- matrix(NA)
  a1 <- matrix(c(1, 0), 2, 1)
  P1 <- matrix(0, 2, 2)
  P1inf <- diag(2)
  
  if(sum(mv1, na.rm = T) == 0) { mv1[sample(length(mv1), 2)] <- 0.00001 }
  
  model_gaussian <- SSModel(mv1 ~ -1 +
                              SSMcustom(Z = Zt, T = Tt, R = Rt, Q = Qt, a1 = a1, P1 = P1, P1inf = P1inf), 
                            H = Ht)
  fit_gaussian <- fitSSM(model_gaussian, inits = c(0, 0), method = "BFGS") #inits = c(0, 0)
  out_gaussian <- KFS(fit_gaussian$model)
  return(as.numeric(coef(out_gaussian)[,1]))
}

getMyGSSM2 <- function(mv1) {
  require(KFAS)
  #mv1 <- adt[,c("seqnames","160714_A3")] #seqnames == "chr6"
  mv1 <- ts(mv1)
  Zt <- matrix(c(1, 0), 1, 2)
  Ht <- matrix(NA)
  Tt <- matrix(c(1, 0, 1, 1), 2, 2)
  Rt <- matrix(c(1, 0), 2, 1)
  Qt <- matrix(NA)
  a1 <- matrix(c(1, 0), 2, 1)
  P1 <- matrix(0, 2, 2)
  P1inf <- diag(2)
  
  if(sum(mv1, na.rm = T) == 0) { mv1[sample(length(mv1), 2)] <- 0.00001 }
  
  model_gaussian <- SSModel(mv1 ~ -1 +
                              SSMcustom(Z = Zt, T = Tt, R = Rt, Q = Qt, a1 = a1, P1 = P1, P1inf = P1inf), 
                            H = Ht)
  fit_gaussian <- fitSSM(model_gaussian, inits = c(0, 0), method = "L-BFGS-B") #inits = c(0, 0)
  out_gaussian <- KFS(fit_gaussian$model)
  return(as.numeric(coef(out_gaussian)[,1]))
}

getMy3VariateGSSM <- function(mv3, polyDegree = 1) {
  N <- 3
  #In cases where all signal is zero we seed some random values to avoid crash
  if(sum(mv3[,1], na.rm = T) == 0)
    mv3[sample(nrow(mv3), 2),1] <- 0.00001
  if(sum(mv3[,2], na.rm = T) == 0)
    mv3[sample(nrow(mv3), 2),2] <- 0.00001
  if(sum(mv3[,3], na.rm = T) == 0)
    mv3[sample(nrow(mv3), 2),3] <- 0.00001
  
  if(polyDegree == 2) {
    model <- SSModel(mv3 ~ -1 + SSMtrend(degree = 2, Q = list(matrix(NA, N, N), matrix(0, N, N))) +
                       SSMcustom(Z = diag(1, N), T = diag(0, N), Q = matrix(NA, N, N), P1 = matrix(NA, N, N)))
  } else if(polyDegree == 1) {
    model <- SSModel(mv3 ~ -1 + SSMtrend(degree = 1, Q = list(matrix(NA, N, N))) +
                       SSMcustom(Z = diag(1, N), T = diag(0, N), Q = matrix(NA, N, N), P1 = matrix(NA, N, N)))
  }
  
  init <- chol(cov(mv3), pivot = T)
  inits <- rep(c(log(diag(init)), init[upper.tri(init)]), 2)
  inits[is.nan(inits)] <- mean(inits, na.rm=T)
  inits[is.infinite(inits)] <- mean(inits[!is.infinite(inits)], na.rm=T)
  
  fitinit1 <- fitSSM(model, updatefn = updatefn3,
                     inits = inits,
                     method = "BFGS")
  fit1 <- fitSSM(model, updatefn = updatefn3, inits = fitinit1$optim.out$par,
                 method = "BFGS")
  
  out <- KFS(fit1$model) 
  
  return(coef(out))
  
}

updatefn3 <- function(pars, model, ...) {
  Q <- diag(exp(pars[1:3]))
  Q[upper.tri(Q)] <- pars[4:6]  
  model["Q", etas = "level"] <- crossprod(Q)
  Q <- diag(exp(pars[7:9]))
  Q[upper.tri(Q)] <- pars[10:12]
  model["Q", etas = "custom"] <- model["P1", states = "custom"] <- crossprod(Q)
  model
}

getMy2VariateGSSM <- function(mv2, polyDegree = 1) {
  require(KFAS)
  
  N <- 2
  #In cases where all signal is zero we seed some random values to avoid crash
  if(sum(mv2[,1], na.rm = T) == 0)
    mv2[sample(nrow(mv2), 2),1] <- 0.00001
  if(sum(mv2[,2], na.rm = T) == 0)
    mv2[sample(nrow(mv2), 2),2] <- 0.00001
  
  if(polyDegree == 2) { #poly = 2
    model <- SSModel(mv2 ~ -1 + SSMtrend(degree = 2, Q = list(matrix(NA, N, N), matrix(0, N, N))) +
                       SSMcustom(Z = diag(1, N), T = diag(0, N), Q = matrix(NA, N, N), P1 = matrix(NA, N, N)))
  } else {
    model <- SSModel(mv2 ~ -1 + SSMtrend(degree = 1, Q = list(matrix(NA, N, N))) +
                       SSMcustom(Z = diag(1, N), T = diag(0, N), Q = matrix(NA, N, N), P1 = matrix(NA, N, N)))
  }
  
  init <- chol(cov(mv2), pivot = T)
  inits <- rep(c(log(diag(init)), init[upper.tri(init)]), 2)
  inits[is.nan(inits)] <- mean(inits, na.rm=T)
  inits[is.na(inits)] <- mean(inits, na.rm=T)
  inits[is.infinite(inits)] <- mean(inits[!is.infinite(inits)], na.rm=T)
  
  
  fitinit <- fitSSM(model, updatefn = updatefn2,
                    inits = inits,
                    method = "BFGS")
  fit <- fitSSM(model, updatefn = updatefn2, inits = fitinit$optim.out$par,
                method = "BFGS")
  out <- KFS(fit$model) 
  return(coef(out)[,1])
}

getMy2VariateGSSM.raw <- function(mv2, polyDegree = 1) {
  require(KFAS)
  
  N <- 2
  #In cases where all signal is zero we seed some random values to avoid crash
  if(sum(mv2[,1], na.rm = T) == 0)
    mv2[sample(nrow(mv2), 2),1] <- 0.00001
  if(sum(mv2[,2], na.rm = T) == 0)
    mv2[sample(nrow(mv2), 2),2] <- 0.00001
  
  if(polyDegree == 2) { #poly = 2
    model <- SSModel(mv2 ~ -1 + SSMtrend(degree = 2, Q = list(matrix(NA, N, N), matrix(0, N, N))) +
                       SSMcustom(Z = diag(1, N), T = diag(0, N), Q = matrix(NA, N, N), P1 = matrix(NA, N, N)))
  } else {
    model <- SSModel(mv2 ~ -1 + SSMtrend(degree = 1, Q = list(matrix(NA, N, N))) +
                       SSMcustom(Z = diag(1, N), T = diag(0, N), Q = matrix(NA, N, N), P1 = matrix(NA, N, N)))
  }
  
  init <- chol(cov(mv2), pivot = T)
  inits <- rep(c(log(diag(init)), init[upper.tri(init)]), 2)
  inits[is.nan(inits)] <- mean(inits, na.rm=T)
  inits[is.na(inits)] <- mean(inits, na.rm=T)
  inits[is.infinite(inits)] <- mean(inits[!is.infinite(inits)], na.rm=T)
  
  
  fitinit <- fitSSM(model, updatefn = updatefn2,
                    inits = inits,
                    method = "BFGS")
  fit <- fitSSM(model, updatefn = updatefn2, inits = fitinit$optim.out$par,
                method = "BFGS")
  out <- KFS(fit$model) 
  return(coef(out))
}

updatefn2 <- function(pars, model, ...) {
  Q <- diag(exp(pars[1:2]))
  Q[upper.tri(Q)] <- pars[3]  
  model["Q", etas = "level"] <- crossprod(Q)
  Q <- diag(exp(pars[4:5]))
  Q[upper.tri(Q)] <- pars[6]
  model["Q", etas = "custom"] <- model["P1", states = "custom"] <- crossprod(Q)
  model
}


