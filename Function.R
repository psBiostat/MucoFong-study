Boots_Phyglm_LOO_CV_Cont <- function(Data, B){
  
  indseulOTU <- which(apply(Data$X, 2, function(x)sum(x!=0))==1) #Chercher les OTUs qui sont que chez un seul patient
  
  #Imposer que ces individus soit dans l'échantillon boot
  k=0
  indseul <- c()
  for(j in indseulOTU){
    k=k+1
    indseul[k] <- which(Data$X[,j]!=0) #On repère les individus
  }
  
  Res_Beta_Boot <- list()
  #Bootstrap procedure
  for(b in 1:B){
    
    bootid <- sample(seq(1:dim(Data$X)[1])[-indseul], dim(Data$X)[1] - length(indseul), replace = TRUE)
    fit_PhyLasso_Boot <- phyglm_LOO_CV_Cont(x = Data$X[bootid,], y = Data$y[bootid,], taxonomy = Data$Taxonomy,
                                        family = "gaussian")
    Res_Beta_Boot[[b]] <- c(fit_PhyLasso_Boot$a0, fit_PhyLasso_Boot$beta)
    
  }
  return(Res_Beta_Boot)
}
                            
phyglm_LOO_CV_Log <- function (x, y, taxonomy, family = "binomial", weights,
                               offset = NULL, alpha = 1, nlambda = 100,
                               lambda.min.ratio = ifelse(nobs < nvars, 0.01, 1e-08), lambda = NULL,
                               standardize = FALSE, intercept = TRUE, thresh = 1e-06,
                               dfmax = nvars + 1, pmax = min(dfmax * 2 + 20, nvars), exclude,
                               penalty.factor = rep(1, nvars), lower.limits = -Inf,
                               upper.limits = Inf, maxit = 1e+05,
                               type.gaussian = ifelse(nvars < 500, "covariance", "naive"),
                               type.logistic = c("Newton", "modified.Newton"),
                               standardize.response = FALSE,
                               type.multinomial = c("ungrouped", "grouped"),
                               Covariable = NULL)
{
  # browser()
  family = match.arg(family)
  alpha = as.double(alpha)
  this.call = match.call()
  nlam = as.integer(nlambda)
  y = drop(y)
  np = dim(x)
  nobs = as.integer(np[1])
  nvars = as.integer(np[2])
  taxonomy = as.data.frame(taxonomy)
  dimtax <- dim(taxonomy)
  if (dimtax[1] != nvars) stop(paste("number of variables does not match
                                     size of taxonomy"))
  if (!(dimtax[2] > 1)) stop(paste("taxonomy must consist of the otu
                                   level and at least one taxon level
                                   above"))
  if (family == "binomial" & is.null(dim(y))) {y = cbind(y, 1-y)}
  # Normalize covariates.
  if (standardize == TRUE) {x = scale(x)}
  else {x = scale(x, scale = rep(1,nvars))}
  if(is.null(Covariable)){
    modb = glmnet(x, y, family = family, offset = offset,
                  alpha = alpha, lambda = lambda, standardize = FALSE,
                  nlambda = nlam, thresh = thresh * 1e-01)
  }else {
    modb = glmnet(as.matrix(cbind(Covariable, x)), y, family = family, offset = offset,
                  alpha = alpha, lambda = lambda, standardize = FALSE,
                  nlambda = nlam, thresh = thresh * 1e-01,
                  penalty.factor = c(rep(0, dim(Covariable)[2]), rep(1, dim(x)[2])))
  }
  ulam = modb$lambda
  nlam = length(ulam)
  # Initialize intermediates
  ntl = ncol(taxonomy)
  vnames = taxonomy[, ntl]
  taxa = lapply(taxonomy[, -ntl, drop = FALSE], levels)
  b_old = Matrix(rep(0,nvars), sparse = TRUE,
                 dimnames = list(vnames, NULL))
  beta = list()
  b0 = c()
  EC_lam <- c()
  Pred_lamb <- c()
  
  ###############################################################################
  #Step 1: Lambda choice by LOO-CV based on classification error
  ###############################################################################
  for (i in 1:nlam) {
    cat("lambda: ", i, " sur ", nlam, "\n")
    #LOO_CV
    # prediction <- c()
    EC_CV <- c()
    for(k in 1:nrow(x)){
      cat("patient: ", k, " sur ", nrow(x), "\n")
      
      b_old = modb$beta[,i]
      a_old = Matrix(rep(1,nvars), sparse = TRUE,
                     dimnames = list(vnames, NULL))
      d_old = lapply(taxa, function(x) data.frame(d=rep(1,length(x)),
                                                  row.names=x))
      if(is.null(Covariable)){
        red = redistribute6(taxonomy, d_old, a_old, b_old)
      }else{
        red = redistribute6(taxonomy, d_old, a_old, b_old[-c(1:dim(Covariable)[2])])
      }
      a_old = red$ao
      d_old = red$do
      w = phi(taxonomy, d_old) * abs(a_old)
      w[w < 1e-2] = 0.5 * w[w < 1e-2] + 0.5 * 1e-2
      for (j in 1:1000) {
        cat("Iter:", j, "\n")
        # browser()
        if(is.null(Covariable)){
          mod = glmnet(x[-k,], y[-k], alpha = alpha, family = family,
                       standardize = FALSE, lambda = ulam[i], penalty.factor = w^-1,
                       thresh = thresh * 1e-01, maxit = maxit)
        }else{
          
          mod = glmnet(as.matrix(cbind(Covariable,x)[-k,]), y[-k,], alpha = alpha, family = family,
                       standardize = FALSE, lambda = ulam[i], 
                       penalty.factor = as.vector(c(rep(0, dim(Covariable)[2]),as.vector(w)^-1)),
                       thresh = thresh * 1e-01, maxit = maxit)
        }
        b_new = mod$beta
        if(is.null(Covariable)){
          a_new = b_new * w^-1
        }else{
          a_new = b_new[-c(1:dim(Covariable)[2])] * w^-1
        }
        red = redistribute6(taxonomy, d_old, a_old, a_new)
        a_new = red$an
        a_old = red$ao
        d_new = red$do
        w = phi(taxonomy, d_new) * abs(a_new)
        expt = min(50, j + 2)
        w[w < 10^-expt] = 0.5 * w[w < 10^-expt] + 0.5 * 10^-expt
        diff1 = sum(abs(b_new - b_old))
        diff2 = sum(abs(b_new - b_old) / max(abs(b_old), thresh))
        b_old = b_new; d_old = d_new;
        # print(paste("Lambda ", i, " = ", ulam[i], ", Iter ", j, ": ",
        # diff1, diff2, sep=" "))
        if (diff2 < thresh) break
      }
      # beta = c(beta, b_old)
      # b0 = c(b0, mod$a0)
      
      #prediction sur l'obs mise de côté
      beta <- b_old
      if(is.null(Covariable)){
        pred <- mod$a0 + as.numeric(attr(x, "scaled:center") %*% as.vector(beta))
      } else {
        pred <- mod$a0 + as.numeric(c( attr(scale(Covariable),"scaled:center"),attr(x, "scaled:center")) %*% as.vector(beta))
      }
      pi <- 1/(1+exp(-pred))
      if (pi > 0.5) prediction = 1
      else prediction=0
      EC_CV[k] <- abs(y[k] - prediction)
    }
    EC_lam <- cbind(EC_lam, EC_CV)
    # Pred_lam <- cbind(Pred_lam, prediction)
  }
  
  EC <- apply(EC_lam, 2, sum)/length(y)
  lambda_opt <- min(ulam[which.min(EC)])
  ###############################################################################
  
  ###############################################################################
  #Step2 2: Parameter estimation on overall data
  ###############################################################################
  if(is.null(Covariable)){
    modb = glmnet(x, y, family = family, offset = offset,
                  alpha = alpha, lambda = lambda_opt, standardize = FALSE,
                  nlambda = nlam, thresh = thresh * 1e-01)
    
  }else {
    modb = glmnet(as.matrix(cbind(Covariable, x)), y, family = family, offset = offset,
                  alpha = alpha, lambda = lambda_opt, standardize = FALSE,
                  nlambda = nlam, thresh = thresh * 1e-01,
                  penalty.factor = c(rep(0, dim(Covariable)[2]), rep(1, dim(x)[2])))
  }
  
  b_old = modb$beta
  a_old = Matrix(rep(1,nvars), sparse = TRUE,
                 dimnames = list(vnames, NULL))
  d_old = lapply(taxa, function(x) data.frame(d=rep(1,length(x)),
                                              row.names=x))
  if(is.null(Covariable)){
    red = redistribute6(taxonomy, d_old, a_old, b_old)
  }else{
    red = redistribute6(taxonomy, d_old, a_old, b_old[-c(1:dim(Covariable)[2])])
  }
  a_old = red$ao
  d_old = red$do
  w = phi(taxonomy, d_old) * abs(a_old)
  w[w < 1e-2] = 0.5 * w[w < 1e-2] + 0.5 * 1e-2
  for (j in 1:1000) {
    if(is.null(Covariable)){
      mod = glmnet(x, y, alpha = alpha, family = family,
                   standardize = FALSE, lambda = lambda_opt, penalty.factor = w^-1,
                   thresh = thresh * 1e-01, maxit = maxit)
    }else{
      mod = glmnet(as.matrix(cbind(Covariable,x)), y, alpha = alpha, family = family,
                   standardize = FALSE, lambda = lambda_opt, 
                   penalty.factor = c(rep(0, dim(Covariable)[2]),as.vector(w)^-1),
                   thresh = thresh * 1e-01, maxit = maxit)
      
    }
    b_new = mod$beta
    if(is.null(Covariable)){
      a_new = b_new * w^-1
    }else{
      a_new = b_new[-c(1:dim(Covariable)[2])] * w^-1
    }
    red = redistribute6(taxonomy, d_old, a_old, a_new)
    a_new = red$an
    a_old = red$ao
    d_new = red$do
    w = phi(taxonomy, d_new) * abs(a_new)
    expt = min(50, j + 2)
    w[w < 10^-expt] = 0.5 * w[w < 10^-expt] + 0.5 * 10^-expt
    diff1 = sum(abs(b_new - b_old))
    diff2 = sum(abs(b_new - b_old) / max(abs(b_old), thresh))
    b_old = b_new; d_old = d_new;
    # print(paste("Lambda ", i, " = ", ulam[i], ", Iter ", j, ": ",
    # diff1, diff2, sep=" "))
    if (diff2 < thresh) break
  }
  beta = b_old
  b0 = mod$a0
  
  ###############################################################################
  
  
  
  # beta = do.call(cBind, beta)
  beta = beta / attr(x, "scaled:scale")
  if(is.null(Covariable)){
    b0 = b0 - as.numeric(attr(x, "scaled:center") %*% as.vector(beta))
  }else{
    b0 = b0 - as.numeric(c( attr(scale(Covariable),"scaled:center"),attr(x, "scaled:center")) %*% as.vector(beta))
    # b0 = b0 - as.numeric(c( attr(scale(Covariable),"scaled:center"),attr(x, "scaled:center")) %*% as.vector(beta))
  }
  return(list(beta=beta, a0 = b0, lambda = ulam, lambda = lambda_opt, EC_lam = EC_lam))
}





phyglm_LOO_CV_Cont <- function (x, y, taxonomy, family = "gaussian", weights,
                               offset = NULL, alpha = 1, nlambda = 100,
                               lambda.min.ratio = ifelse(nobs < nvars, 0.01, 1e-08), lambda = NULL,
                               standardize = FALSE, intercept = TRUE, thresh = 1e-06,
                               dfmax = nvars + 1, pmax = min(dfmax * 2 + 20, nvars), exclude,
                               penalty.factor = rep(1, nvars), lower.limits = -Inf,
                               upper.limits = Inf, maxit = 1e+05,
                               type.gaussian = ifelse(nvars < 500, "covariance", "naive"),
                               type.logistic = c("Newton", "modified.Newton"),
                               standardize.response = FALSE,
                               type.multinomial = c("ungrouped", "grouped"),
                               Covariable = NULL)
{
  # browser()
  family = match.arg(family)
  alpha = as.double(alpha)
  this.call = match.call()
  nlam = as.integer(nlambda)
  y = drop(y)
  np = dim(x)
  # browser()
  nobs = as.integer(np[1])
  nvars = as.integer(np[2])
  taxonomy = as.data.frame(taxonomy)
  dimtax <- dim(taxonomy)
  # browser()
  if (dimtax[1] != nvars) stop(paste("number of variables does not match
                                     size of taxonomy"))
  if (!(dimtax[2] > 1)) stop(paste("taxonomy must consist of the otu
                                   level and at least one taxon level
                                   above"))
  if (family == "binomial" & is.null(dim(y))) {y = cbind(y, 1-y)}
  # browser()
  # Normalize covariates.
  if (standardize == TRUE) {x = scale(x)}
  else {x = scale(x, scale = rep(1,nvars))}
  # browser()
  if(is.null(Covariable)){
    modb = glmnet(as.matrix(x), y, family = family, offset = offset,
                  alpha = alpha, lambda = lambda, standardize = FALSE,
                  nlambda = nlam, thresh = thresh * 1e-01)
  }else {
    modb = glmnet(as.matrix(cbind(Covariable, x)), y, family = family, offset = offset,
                  alpha = alpha, lambda = lambda, standardize = FALSE,
                  nlambda = nlam, thresh = thresh * 1e-01,
                  penalty.factor = c(rep(0, dim(Covariable)[2]), rep(1, dim(x)[2])))
  }
  ulam = modb$lambda
  nlam = length(ulam)
  # Initialize intermediates
  ntl = ncol(taxonomy)
  vnames = taxonomy[, ntl]
  taxa = lapply(taxonomy[, -ntl, drop = FALSE], levels)
  b_old = Matrix(rep(0,nvars), sparse = TRUE,
                 dimnames = list(vnames, NULL))
  beta = list()
  b0 = c()
  MSE_lam <- c()
  Pred_lamb <- c()

  
  ###############################################################################
  #Etape 1: Choix du lambda LOO-CV basé sur l'EC
  ###############################################################################
  for (i in 1:nlam) {
    cat("lambda: ", i, " sur ", nlam, "\n")
    #LOO_CV
    # prediction <- c()
    EC_MSE <- c()
    for(k in 1:nrow(x)){
      # browser()
      cat("patient: ", k, " sur ", nrow(x), "\n")
      
      b_old = modb$beta[,i]
      a_old = Matrix(rep(1,nvars), sparse = TRUE,
                     dimnames = list(vnames, NULL))
      d_old = lapply(taxa, function(x) data.frame(d=rep(1,length(x)),
                                                  row.names=x))
      if(is.null(Covariable)){
        red = redistribute6(taxonomy, d_old, a_old, b_old)
      }else{
        red = redistribute6(taxonomy, d_old, a_old, b_old[-c(1:dim(Covariable)[2])])
      }
      a_old = red$ao
      d_old = red$do
      w = phi(taxonomy, d_old) * abs(a_old)
      w[w < 1e-2] = 0.5 * w[w < 1e-2] + 0.5 * 1e-2
      for (j in 1:1000) {
        cat("Iter:", j, "\n")
        # browser()
        if(is.null(Covariable)){
          mod = glmnet(as.matrix(x)[-k,], y[-k], alpha = alpha, family = family,
                       standardize = FALSE, lambda = ulam[i], penalty.factor = w^-1,
                       thresh = thresh * 1e-01, maxit = maxit)
        }else{
          
          mod = glmnet(as.matrix(cbind(Covariable,x)[-k,]), y[-k], alpha = alpha, family = family,
                       standardize = FALSE, lambda = ulam[i], 
                       penalty.factor = as.vector(c(rep(0, dim(Covariable)[2]),as.vector(w)^-1)),
                       thresh = thresh * 1e-01, maxit = maxit)
        }
      
        b_new = mod$beta
        
        
        if(is.null(Covariable)){
          a_new = b_new * w^-1
        }else{
          a_new = b_new[-c(1:dim(Covariable)[2])] * w^-1
        }
        red = redistribute6(taxonomy, d_old, a_old, a_new)
        a_new = red$an
        a_old = red$ao
        d_new = red$do
        w = phi(taxonomy, d_new) * abs(a_new)
        expt = min(50, j + 2)
        w[w < 10^-expt] = 0.5 * w[w < 10^-expt] + 0.5 * 10^-expt
        diff1 = sum(abs(b_new - b_old))
        diff2 = sum(abs(b_new - b_old) / max(abs(b_old), thresh))
        b_old = b_new; d_old = d_new;
        # print(paste("Lambda ", i, " = ", ulam[i], ", Iter ", j, ": ",
        # diff1, diff2, sep=" "))
        if (diff2 < thresh) break
      }
      # browser()
      beta <- b_old
      # if(is.null(Covariable)){
      #   pred <- mod$a0 + as.numeric(attr(x, "scaled:center") * as.vector(beta))
      # } else {
      #   pred <- mod$a0 + as.numeric(c( attr(scale(Covariable),"scaled:center"),attr(x, "scaled:center")) %*% as.vector(beta))
      # }
      # 
      if(is.null(Covariable)){
        pred <- mod$a0 + as.matrix(x)[k,] %*% as.vector(beta)
      } else {
        pred <- mod$a0 + as.numeric(c( attr(scale(Covariable),"scaled:center"),attr(x, "scaled:center")) %*% as.vector(beta))
      }
      

      EC_MSE[k] <-abs(y[k] - pred)
    }
    
    MSE_lam <- cbind(MSE_lam, EC_MSE)
    # Pred_lam <- cbind(Pred_lam, prediction)
  }
  browser()
  MSE <- apply(MSE_lam, 2, sum)/length(y)
  lambda_opt <- min(ulam[which.min(MSE)])
  # browser()
  ###############################################################################
  
  ###############################################################################
  #Etape 2: Estimation des paramètres sur l'ensemble des données
  ###############################################################################
  if(is.null(Covariable)){
    modb = glmnet(as.matrix(x), y, family = family, offset = offset,
                  alpha = alpha, lambda = lambda_opt, standardize = FALSE,
                  nlambda = nlam, thresh = thresh * 1e-01)
    
  }else {
    modb = glmnet(as.matrix(cbind(Covariable, x)), y, family = family, offset = offset,
                  alpha = alpha, lambda = lambda_opt, standardize = FALSE,
                  nlambda = nlam, thresh = thresh * 1e-01,
                  penalty.factor = c(rep(0, dim(Covariable)[2]), rep(1, dim(x)[2])))
  }
  
  b_old = modb$beta
  a_old = Matrix(rep(1,nvars), sparse = TRUE,
                 dimnames = list(vnames, NULL))
  d_old = lapply(taxa, function(x) data.frame(d=rep(1,length(x)),
                                              row.names=x))
  if(is.null(Covariable)){
    red = redistribute6(taxonomy, d_old, a_old, b_old)
  }else{
    red = redistribute6(taxonomy, d_old, a_old, b_old[-c(1:dim(Covariable)[2])])
  }
  a_old = red$ao
  d_old = red$do
  w = phi(taxonomy, d_old) * abs(a_old)
  w[w < 1e-2] = 0.5 * w[w < 1e-2] + 0.5 * 1e-2
  for (j in 1:1000) {
    if(is.null(Covariable)){
      mod = glmnet(as.matrix(x), y, alpha = alpha, family = family,
                   standardize = FALSE, lambda = lambda_opt, penalty.factor = w^-1,
                   thresh = thresh * 1e-01, maxit = maxit)
    }else{
      mod = glmnet(as.matrix(cbind(Covariable,x)), y, alpha = alpha, family = family,
                   standardize = FALSE, lambda = lambda_opt, 
                   penalty.factor = c(rep(0, dim(Covariable)[2]),as.vector(w)^-1),
                   thresh = thresh * 1e-01, maxit = maxit)
     
    }
    b_new = mod$beta
    if(is.null(Covariable)){
      a_new = b_new * w^-1
    }else{
      a_new = b_new[-c(1:dim(Covariable)[2])] * w^-1
    }
    red = redistribute6(taxonomy, d_old, a_old, a_new)
    a_new = red$an
    a_old = red$ao
    d_new = red$do
    w = phi(taxonomy, d_new) * abs(a_new)
    expt = min(50, j + 2)
    w[w < 10^-expt] = 0.5 * w[w < 10^-expt] + 0.5 * 10^-expt
    diff1 = sum(abs(b_new - b_old))
    diff2 = sum(abs(b_new - b_old) / max(abs(b_old), thresh))
    b_old = b_new; d_old = d_new;
    # print(paste("Lambda ", i, " = ", ulam[i], ", Iter ", j, ": ",
    # diff1, diff2, sep=" "))
    if (diff2 < thresh) break
  }
  beta = b_old
  b0 = mod$a0
  
  ###############################################################################
  
  # browser()
  
  # beta = do.call(cBind, beta)
  # beta = beta / attr(x, "scaled:scale")
  # if(is.null(Covariable)){
  #   b0 = b0 - as.numeric(attr(x, "scaled:center") %*% as.vector(beta))
  # }else{
  #   b0 = b0 - as.numeric(c( attr(scale(Covariable),"scaled:center"),attr(x, "scaled:center")) %*% as.vector(beta))
  #   # b0 = b0 - as.numeric(c( attr(scale(Covariable),"scaled:center"),attr(x, "scaled:center")) %*% as.vector(beta))
  # }
  # browser()
  # #beta = beta / attr(x, "scaled:scale")
  # if(is.null(Covariable)){
  #   b0 = b0 - as.matrix(x) %*% as.vector(beta)
  # }else{
  #   b0 = b0 - as.numeric(c( attr(scale(Covariable),"scaled:center"),attr(x, "scaled:center")) %*% as.vector(beta))
  #   # b0 = b0 - as.numeric(c( attr(scale(Covariable),"scaled:center"),attr(x, "scaled:center")) %*% as.vector(beta))
  # }
  return(list(beta=beta, a0 = b0, lambda = ulam, lambda = lambda_opt, MSE_lam = MSE_lam))
}
