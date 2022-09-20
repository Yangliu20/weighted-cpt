# Define functions to simulate data

library(mvtnorm)
library(glmnet)

#' Generate auto-regressive covariance matrix
#' @param ar autocorrelation
#' 
AR_cov <- function(p, ar){
  ar_series <- ar^(c(1:p) - 1)
  cov_mat <- ar_series
  for (i in 1:(p - 1)){
    cov_mat <- cbind(cov_mat, ar_series[c((p - i + 1):p ,1:(p - i))])
  }
  for (i in 1:(p - 1)){
    for (j in (i + 1):p){
      cov_mat[i, j] <- cov_mat[j, i]
    }
  }
  rownames(cov_mat) <- colnames(cov_mat)
  return(cov_mat)
}


#' Compute conditional guassian distribution
Creat_condition_gaussian <- function(XY, indx, X.type, model.X.known, prob.mixture, X_unlabel=NULL){
  X = XY$X
  num.mixture = length(prob.mixture)
  if(model.X.known){
    if(X.type == 'Zmixture'){
      cond = list(mean_x = XY$cond.mean_x, sigma2_x = XY$cond.sigma2_x, prob.mixture = XY$cond.prob.mixture)
      return(cond)
    }
    if(X.type == 'linear'){
      cond = list(mean_x = XY$model.X.intercept+X[,-j]%*%XY$model.X.coeff, sigma2_x = XY$model.X.sd^2)
      return(cond)
    }
    if((X.type == 'AR') | (X.type == 'Equicorr')){
      mu = XY$mu_true; Sigma = XY$Sigma_true
      beta_x <- solve(Sigma[-indx, -indx], Sigma[-indx, indx])
      sigma2_x <- Sigma[indx, indx] - Sigma[indx, -indx] %*% beta_x
      if(is.null(prob.mixture)){ # single normal
        X_bar <- X[ ,-indx] %*% beta_x - as.numeric(mu[-indx] %*% beta_x) + mu[indx]
      }else{ # mixture of normals
        X_bar = c()
        for(k in 1:num.mixture){
          X_bar = cbind(X_bar, X[ ,-indx] %*% beta_x - as.numeric(mu[-indx, k] %*% beta_x) + mu[indx, k])
        }
      }
      return(list(mean_x = X_bar, sigma2_x = sigma2_x, gamma = beta_x, prob.mixture = prob.mixture))
    }
  }else if(!model.X.known){
    if((X.type == 'linear') | (is.null(prob.mixture))){
      # linear or single normal
      
      # cv_lasso_x <- cv.glmnet(X_unlabel$X[,-indx], X_unlabel$X[,indx], alpha=1, intercept=T, dfmax=as.integer(p/2))
      # lambda_x <- cv_lasso_x$lambda.min
      # opt_model_x <- glmnet(X_unlabel$X[,-indx], X_unlabel$X[,indx], alpha=1, lambda=lambda_x, intercept=T, dfmax=as.integer(p/2))
      # beta_x <- opt_model_x$beta
      # X_bar <- predict(opt_model_x, X[,-indx])
      # x_res <- X_unlabel$X[,indx] - predict(opt_model_x, X_unlabel$X[,-indx]) 
      # sigma2_x <- mean(x_res^2)
      
      model = lm(X_unlabel$X[,indx] ~ X_unlabel$X[,-indx])
      beta_x = model$coefficients[-1]
      intercept = model$coefficients[1]
      sigma2_x = sum((model$residuals)^2)/model$df.residual
      X_bar = X[,-indx] %*% beta_x + intercept
      
      return(list(mean_x = X_bar, sigma2_x = sigma2_x, gamma = beta_x, prob.mixture = prob.mixture))
    }
    if((X.type == 'Zmixture') | (X.type == 'AR') | (X.type == 'Equicorr')){
      # mixture model
      dummy = model.matrix(~as.factor(X_unlabel$cluster)-1)
      x = cbind(dummy, X_unlabel$X[,-indx])
      y = X_unlabel$X[,indx]
      cv_lasso_x = cv.glmnet(x, y, alpha=1, intercept=F)
      lambda_x = cv_lasso_x$lambda.min
      opt_model_x = glmnet(x, y, alpha=1, intercept=F, lambda = lambda_x)
      
      # intercept.mixture = opt_model_x$beta[1:num.mixture]
      # beta.x.z = opt_model_x$beta[-(1:num.mixture)]
      prob.mixture.hat = as.numeric(table(X_unlabel$cluster)/length(X_unlabel$cluster))
      # mean_x = t(apply(X[,-indx] %*% beta.x.z, 1, function(x){x + intercept.mixture}))
      sigma2_x = mean((y-predict(opt_model_x, x))^2)
      mean_x = c()
      for(i in 1:num.mixture){
        x.mixture = cbind(0,0,0,0,X[,-indx])
        x.mixture[,i] = 1
        mean_x = cbind(mean_x, predict(opt_model_x, x.mixture))
      }
      
      return(list(mean_x=mean_x, sigma2_x=sigma2_x, prob.mixture=prob.mixture.hat, 
                  param=opt_model_x$beta))
    }
  }
}


generate_covariates = function(n, p, X.type, ar_coeff = NULL, prob.mixture = NULL, sd.a = NULL, mix.mean.sd = NULL, mis.specify = NULL){
  
  if(X.type=='linear'){
    Z = matrix(rnorm(n*(p-1)), n, p-1)
    a = rnorm(p-1) * sd.a
    X = Z%*%a + rnorm(n)
    return(list(X=cbind(X,Z), model.X.coeff=a, model.X.intercept=0, model.X.sd=1)) 
  }
  
  if(X.type=='Zmixture'){
    num.mixture = length(prob.mixture)
    mu.Z = matrix(rnorm(num.mixture*(p-1), sd=mix.mean.sd), nrow=p-1, ncol=num.mixture)
    g = sample(1:num.mixture, n, replace=TRUE, prob=prob.mixture)
    Z = t(mu.Z[,g]) + matrix(rnorm(n*(p-1)), n, p-1) # dim=n*(p-1)
    beta.x.z = rnorm(p-1)*sd.a
    intercept.mixture = rnorm(num.mixture, sd=mix.mean.sd)
    intercept.x = intercept.mixture[g]
    X = intercept.x + Z %*% beta.x.z + rnorm(n)
    cond.mean_x = t(apply(Z %*% beta.x.z, 1, function(x){x + intercept.mixture}))
    return(list(X=cbind(X,Z), cluster=g, 
                cond.mean_x=cond.mean_x, 
                cond.sigma2_x=1, cond.prob.mixture=prob.mixture, 
                intercept.true = intercept.mixture, 
                beta.x.z.true = beta.x.z))
  }
  
  if(X.type=='AR'){ # AR covariance
    Sigma_true = AR_cov(p, ar_coeff)
  }else if(X.type=='Equicorr'){ # Equicorrelation
    Sigma_true = matrix(ar_coeff, p, p)
    diag(Sigma_true) = 1
  }
  
  if(is.null(prob.mixture)){ # single normal
    X = rmvnorm(n, rep(0,p), Sigma_true)
    mu = rep(0,p)
    return(list(X = X, Sigma_true = Sigma_true, mu_true = mu))
  }else{ # mixture of normals
    num.mixture = length(prob.mixture)
    g = sample(1:num.mixture, n, replace=TRUE, prob=prob.mixture)
    mu = matrix(rnorm(num.mixture*p, sd=mix.mean.sd), nrow=p, ncol=num.mixture) # p*num.mixture
    X = t(mu[,g]) + rmvnorm(n, sigma=Sigma_true)
    return(list(X = X, Sigma_true = Sigma_true, mu_true = mu, cluster=g))
  }
}

generate_coeff = function(p, eff_size, signal){
  beta_true = rnorm(p, 0, 1) * eff_size
  beta_true[1] = signal
  return(beta_true)
}

generate_response = function(X, beta, intercept, sd.err){
  intercept + X %*% beta + rnorm(n, 0, sd.err)
}

generate_XY = function(n, p, X.type, eff_size, signal, 
                       ar_coeff = NULL, prob.mixture = NULL, sd.a = NULL, mix.mean.sd = NULL, 
                       mis.specify = NULL, n_unlabel = 0){
  covariate = generate_covariates(n+n_unlabel, p, X.type, ar_coeff, prob.mixture, sd.a, mix.mean.sd, mis.specify)
  X = covariate$X[1:n,]
  beta_true = generate_coeff(p, eff_size, signal) 
  Y = generate_response(X, beta_true, 0, 1)
  if(X.type == 'linear'){
    return(list(X=X, Y=Y, X_unlabel = covariate$X[-(1:n),],
                model.X.coeff = covariate$model.X.coeff, 
                model.X.intercept = covariate$model.X.intercept, 
                model.X.sd = covariate$model.X.sd))}
  if(X.type == 'Zmixture'){
    return(list(X=X, Y=Y, cluster_true=covariate$cluster[1:n],
                cond.mean_x=covariate$cond.mean_x[1:n,], 
                cond.sigma2_x=covariate$cond.sigma2_x, 
                cond.prob.mixture=covariate$cond.prob.mixture, 
                X_unlabel=covariate$X[-(1:n),],
                cluster_true_unlabel=covariate$cluster[-(1:n)], 
                cond.mean_x_unlabel=covariate$cond.mean_x[-(1:n),], 
                intercept.mixture.true=covariate$intercept.true, 
                beta.x.z.true=covariate$beta.x.z.true))
  }
  if((X.type == 'AR') | (X.type == 'Equicorr')){
    Sigma = covariate$Sigma_true; mu = covariate$mu_true
    if(is.null(prob.mixture)){
      # single normal
      return(list(X=X, Y=Y, mu_true=mu, beta_true=beta_true, Sigma_true=Sigma,
                  X_unlabel=covariate$X[-(1:n),]))
    }
    return(list(X=X, Y=Y, mu_true=mu, beta_true=beta_true, Sigma_true=Sigma,
                cluster_true=covariate$cluster[1:n], 
                X_unlabel=covariate$X[-(1:n),],
                cluster_true_unlabel=covariate$cluster[-(1:n)]))
  }
}



# Define functions of efficient permutation

#' Propose a permutation efficiently
#' @param prob probability of permuting within group
#' @param g a vector of group assignment
#' 
propose.perm = function(prob, g){
  k.arr = unique(g)
  n = length(g)
  if(runif(1) < prob){ # permute within groups
    perm = numeric(n)
    for(c in k.arr){
      if(sum((g == c)) == 1){
        perm[which(g == c)] = which(g == c)
      }else{
        perm[which(g == c)] = sample(which(g == c))
      }
    }
    return(list(within.group=1, perm=perm))
  }else{ # permute completely at random
    return(list(within.group=0, perm=sample(1:n)))
  }
}

#' Propose M permutations efficiently
#' @return matrix of dim n*(M+1), with first column being 1:n
#' 
propose.perm.M = function(prob, g, M, type='naive'){
  if(type=='naive'){
    n = length(g)
    Xcopy = matrix(rep(1:n, M+1), nrow=n, ncol=M+1)
    n_within = rbinom(1, M, prob)
    temp = c(rep(1, n_within+1), rep(0, M-n_within))
    k.arr = unique(g)
    for(c in k.arr){
      if(sum((g == c)) > 1){
        Xcopy[which(g == c), 2:(n_within+1)] = apply(Xcopy[which(g == c), 2:(n_within+1)], 2, function(x){sample(x)})
      }
    }
    if(n_within < M){
      Xcopy[,-(2:(n_within+1))] = apply(Xcopy[,-(2:(n_within+1))], 2, function(x){sample(x)})
    }
    # for(m in 2:(M+1)){
    #   proposal = propose.perm(prob, g)
    #   Xcopy[,m] = proposal$perm
    #   if(proposal$within.group == 1){
    #     temp[m] = 1
    #   }else{
    #     temp[m] = all(g[proposal$perm] == g)*1
    #   }
    # }
    const = multicool::multinom(as.numeric(table(g)))
    log.prop.den = log(prob * const * temp + (1-prob))
  }
  list(Xcopy=Xcopy, log.prop.den=log.prop.den, within.group=temp)
}

#' Compute the log target density
#' @param x a vector of observations
#' @param permutations matrix of dim n*M
#' @param mu conditional mean (matrix if mixture model)
#' @param sigma2 conditional variance
#' @param prob.mixture probability of mixture components if applicable
#' @return log target density for each permutation, up to constant
#' 
log.target.den.perm = function(x, permutations, mu, sigma2, prob.mixture){
  if(is.null(prob.mixture)){
    return(apply(permutations, 2, 
                 function(z){sum(dnorm(x[z], mean=mu, sd=sqrt(sigma2), log=TRUE))}))
  }else{
    M = length(permutations[1,])
    num.mixture = length(mu[1,])
    
    lik_mat = matrix(0, length(x), length(x))
    sig2 = rep(sigma2, length(x))
    for(k in 1:num.mixture){
      lik_mat = lik_mat + prob.mixture[k] * exp(-(x^2)%*%t(1/2/sig2) + x%*%t(mu[,k]/sig2) - rep(1,length(x))%*% t(mu[,k]^2/sig2/2))
    }
    log_lik_mat = log(lik_mat) # log_lik_mat[i,j] is the log lik of X[i]|Z[j]
    
    # log.den = numeric(M)
    # for(m in 1:M){ # for each copy
    #   log.den[m] = sum(log_lik_mat[cbind(permutations[,m], 1:length(x))])
    # }
    log.den = apply(permutations, 2, function(z){sum(log_lik_mat[cbind(z, 1:length(x))])})
    return(log.den)
  }
}

#' Create groups
#' @param X data matrix of dimension n*p
#' @param j index of interest
#' @param k number of groups
#' @param type kmeans or gmm
#' @param prob probability of components, length k, if applicable
#' @param mu mean matrix of dimension p*k
#' @param Sigma covariance matrix (common) of dim p*p
#' 
create_groups = function(x, k, type='kmeans', prob=NULL, mu=NULL, Sigma=NULL){
  if(is.null(dim(x))){
    return(kmeans(x, k))
  }else{
    return(kmeans(x[,1], k))
  }
}


weighted_CPT = function(Y, X, indx, prob, g, M, cond, test.statistic = 'cor.X.res', res=NULL){
  starttime = Sys.time()
  proposals = propose.perm.M(prob, g$cluster, M) # generate copies
  Xcopy = proposals$Xcopy; log.prop.den = proposals$log.prop.den 
  log.target.den = log.target.den.perm(X[,indx], Xcopy, cond$mean_x, cond$sigma2_x, cond$prob.mixture)
  log.weight = log.target.den - log.prop.den
  weight = exp(log.weight - max(log.weight))/sum(exp(log.weight - max(log.weight)))
  
  # test statistic
  if(test.statistic == 'cor.X.res'){
    t = apply(Xcopy, 2, function(z){abs(cor(X[z, indx], res))})
  }else if(test.statistic == 'cor.X.Y'){
    t = apply(Xcopy, 2, function(z){abs(cor(X[z, indx], Y))})
  }
  
  endtime = Sys.time()
  list(pval = sum(weight * (t >= t[1])), ess = 1/sum(weight^2), comp_time = as.numeric(endtime-starttime))
}




# Define functions of original CPT

generate_X_CPT_gaussian = function(nstep,M,X0,mu,sig2,prob.mixture){
  # Runs the conditional permutation test using the distribution X | Z=Z[i] ~ N(mu[i],sig2[i])
  if(is.null(prob.mixture)){
    log_lik_mat = -(X0^2)%*%t(1/2/sig2) + X0%*%t(mu/sig2) - rep(1,length(X0))%*% t(mu^2/sig2/2)
  }else{
    num.mixture = length(mu[1,])
    lik_mat = matrix(0, length(X0), length(X0))
    for(k in 1:num.mixture){
      lik_mat = lik_mat + prob.mixture[k] * exp(-(X0^2)%*%t(1/2/sig2) + X0%*%t(mu[,k]/sig2) - rep(1,length(X0))%*% t(mu[,k]^2/sig2/2))
    }
    log_lik_mat = log(lik_mat)
  }
  # log_lik_mat[i,j] = density at X=X0[i] when Z=Z[j]
  Pi_mat = generate_X_CPT(nstep,M,log_lik_mat)
  # X_mat = X0[Pi_mat]
  # dim(X_mat) = c(M,length(X0))
  return(cbind(1:n, t(Pi_mat)))
}

generate_X_CPT = function(nstep,M,log_lik_mat,Pi_init=NULL){
  # log_lik_mat is the n-by-n matrix with entries log(q(X_i|Z_j))
  # this function produces M exchangeable permutations, initialized with permutation Pi_init
  n = dim(log_lik_mat)[1]
  if(length(Pi_init)==0){
    Pi_init = 1:n
  }
  Pi_ = generate_X_CPT_MC(nstep,log_lik_mat,Pi_init)
  Pi_mat = matrix(0,M,n)
  for(m in 1:M){
    Pi_mat[m,] = generate_X_CPT_MC(nstep,log_lik_mat,Pi_)
  }
  return(Pi_mat)
}

generate_X_CPT_MC = function(nstep,log_lik_mat,Pi_){
  # log_lik_mat is the n-by-n matrix with entries log(q(X_i|Z_j))
  # this function runs the MC sampler, initialized with permutation Pi_
  n = length(Pi_)
  npair = floor(n/2)
  for(istep in 1:nstep){
    perm = sample(n)
    inds_i = perm[1:npair]
    inds_j = perm[(npair+1):(2*npair)]
    # for each k=1,...,npair, deciding whether to swap Pi_[inds_i[k]] with Pi_[inds_j[k]]
    log_odds = (log_lik_mat[cbind(Pi_[inds_i],inds_j)] + log_lik_mat[cbind(Pi_[inds_j],inds_i)]
                - log_lik_mat[cbind(Pi_[inds_i],inds_i)] - log_lik_mat[cbind(Pi_[inds_j],inds_j)])
    swaps = rbinom(npair,1,1/(1+exp(-pmax(-500,log_odds))))
    Pi_[c(inds_i,inds_j)] = Pi_[c(inds_i,inds_j)] + swaps*(Pi_[c(inds_j,inds_i)]-Pi_[c(inds_i,inds_j)])
  }
  return(Pi_)
}

CPT = function(Y, X, indx, nstep, M, cond, test.statistic='cor.X.res', res=NULL){
  starttime = Sys.time()
  Xcopy_CPT = generate_X_CPT_gaussian(nstep,M,X[,indx],cond$mean_x,rep(cond$sigma2_x, n),cond$prob.mixture)
  if(test.statistic=='cor.X.res'){
    # res = lm(Y~X[,-indx])$residuals
    t_CPT = apply(Xcopy_CPT, 2, function(z){abs(cor(X[z, indx], res))})
  }else if(test.statistic=='cor.X.Y'){
    t_CPT = apply(Xcopy_CPT, 2, function(z){abs(cor(X[z, indx], Y))})
  }
  endtime = Sys.time()
  list(pval = mean(t_CPT >= t_CPT[1]), comp_time = as.numeric(endtime-starttime))
}




# Define plotting function

plot_rejection_rate = function(pval_weighted, pval_CPT, x.list, x.label, cpt.run, alpha=0.05, legend.bool=TRUE){
  rej_rate = apply(pval_weighted, c(1,2,3), function(x){mean(x<=alpha)})
  
  if(X.type=='AR'){
    title = paste0('n=',n,'; p=',p,'; ',group.method,'; ar_coeff=',ar_coeff,'; mixtures=',length(prob.mixture))
    fig.name = paste0('n',n,'_p',p,'_',X.type,'_ar',as.integer(ar_coeff*10),'_mix',length(prob.mixture))
    if(!is.null(mix.mean.sd)){
      fig.name = paste0(fig.name, '_mixMeanSd',mix.mean.sd)
    }
    fig.name = paste0(fig.name, '.png')
  }
  if(X.type=='linear'){
    title = paste0('n=',n,'; p=',p,'; ',group.method,'; sd.a=',sd.a)
    fig.name = paste0('n',n,'_p',p,'_',X.type,'_sda', sd.a,'.png')
  }
  if(X.type=='Zmixture'){
    title = paste0('n=',n,'; p=',p,'; ',group.method,'; sd_a=',sd.a,'; mix_mean_sd=',mix.mean.sd,'; mixtures=',length(prob.mixture))
    # fig.name = paste0('n',n,'_p',p,'_',X.type,'_ar',as.integer(ar_coeff*10),'_mix',length(prob.mixture),'.png')
  }
  if(X.type=='Equicorr'){
    title = paste0('n=',n,'; p=',p,'; ',group.method,'; equicorr=',ar_coeff,'; mixtures=',length(prob.mixture))
    fig.name = paste0('n',n,'_p',p,'_',X.type,'_',ar_coeff,'_mix',length(prob.mixture))
    if(!is.null(mix.mean.sd)){
      fig.name = paste0(fig.name, '_mixMeanSd',mix.mean.sd)
    }
    fig.name = paste0(fig.name, '.png')
  }
  
  pch.list = c(15, 16, 17, 18, 19) # k
  col.list = c(2, 3, 4, 5, 6) # prob
  label = c(); pchs = c(); cols = c()
  y.max = min(1, max(max(rej_rate), max(rowMeans(pval_CPT <= alpha)))*1.2)
  # y.max = 1
  
  png(fig.name, width = 3.25, height = 3.25, units = "in", res = 1200)
  
  par(cex.main = 0.8, mai=c(.6,.6,.25,.05),mgp=c(2,.7,0))
  if(cpt.run){
    plot(x.list, rowMeans(pval_CPT <= alpha), type = 'b', pch=22, col=1, lwd=1., lty=2,
         xlab = x.label, ylab='Rejection rate', ylim=c(0,y.max))#, main=title)
  }else{
    plot(x.list, rep(0, length(x.list)), col=0, 
         xlab = x.label, ylab='Rejection rate', ylim=c(0,y.max))#, main=title)
  }
  
  for(k.i in 1:length(k.list)){
    for(prob.i in 1:length(prob.list)){
      points(x.list, rej_rate[, k.i, prob.i], type = 'b',
             pch=pch.list[k.i], col=col.list[k.i])
      label = c(label, paste0('K=',k.list[k.i]))#,', prob=',prob.list[prob.i]))
      pchs = c(pchs, pch.list[k.i]); cols = c(cols, col.list[k.i])
    }
  }
  abline(h=alpha, lty=3, lwd=1.)
  if(legend.bool){
    legend('topleft', legend=c(label,paste0('CPT')),
           lty = c(rep(1, length(label)), 2), lwd=rep(1., length(k.list)+1),
           pch = c(pchs, 22), col = c(cols, 1), cex=0.8)
  }
  
  dev.off()
  
}



# Small size data simulation

## Assume model-X is known

compare_efficient_CPT = function(n, p, M_CPT, M_efficient, 
                                 X.type, ar_coeff, prob.mixture, sd.a,mix.mean.sd,
                                 signal.list, k.list, prob.list, nstep, num.sim, test.statistic='cor.X.res',
                                 cpt.run=TRUE, alpha=0.05, legend.bool=TRUE){
  
  pval = array(0,dim=c(length(signal.list),length(k.list),length(prob.list),num.sim))
  ess = array(0,dim=c(length(signal.list),length(k.list),length(prob.list),num.sim))
  comp_time = array(0,dim=c(length(signal.list),length(k.list),length(prob.list),num.sim))
  pval_CPT = array(0,dim=c(length(signal.list),num.sim))
  comp_time_CPT = array(0,dim=c(length(signal.list),num.sim))
  
  
  for(signal.i in 1:length(signal.list)){
    signal = signal.list[signal.i]
    message('\nSignal = ', signal)
    for(num in 1:num.sim){
      if(num%%10==0){
        message(num,' of ', num.sim)
      }
      
      XY = generate_XY(n, p, X.type, eff_size, signal, ar_coeff, prob.mixture, sd.a, mix.mean.sd)
      X = XY$X; Y = XY$Y
      cond = Creat_condition_gaussian(XY, j, X.type, TRUE, prob.mixture)
      res = lm(Y~X[,-j])$residuals
      
      ### efficient sampler + weighting
      for(k.i in 1:length(k.list)){
        for(prob.i in 1:length(prob.list)){
          k = k.list[k.i]
          prob = prob.list[prob.i]
          # g = create_groups(X, j, k, 'kmeans')
          g = create_groups(cond$mean_x, k, 'kmeans')
          
          result = weighted_CPT(Y, X, j, prob, g, M_efficient, cond, test.statistic, res)
          ess[signal.i, k.i, prob.i, num] = result$ess
          pval[signal.i, k.i, prob.i, num] = result$pval
          comp_time[signal.i, k.i, prob.i, num] = result$comp_time
        }
      }
      ### original CPT
      if(cpt.run){
        result_CPT = CPT(Y, X, j, nstep, M_CPT, cond, test.statistic, res)
        pval_CPT[signal.i, num] = result_CPT$pval
        comp_time_CPT[signal.i, num] = result_CPT$comp_time
      }
    }
  }
  
  plot_rejection_rate(pval, pval_CPT, signal.list, 'Signal strength', cpt.run, alpha, legend.bool)
  
  ess_k_prob = apply(ess, c(2,3), mean)
  dimnames(ess_k_prob) = list(paste0('k=',k.list), paste0('prob=',prob.list))
  rej = apply(pval, c(1,2,3), function(x){mean(x<=alpha)})
  dimnames(rej) = list(paste0('signal=',signal.list), paste0('k=',k.list), paste0('prob=',prob.list))
  
  cat(paste('\n\nComputational time of efficient method:',round(mean(comp_time),2)))
  if(cpt.run){
    cat(paste('\nComputational time of original CPT:',round(mean(comp_time_CPT),2)))
  }
  
  cat(paste('\n\nEffective sample size:\n'))
  print(ess_k_prob)
  cat('\n\nRejection rate of efficient method:\n')
  print(rej)
  if(cpt.run){
    cat('\n\nRejection rate of CPT:\n')
    print(rowMeans(pval_CPT <= alpha))
  }
  
  list(pval=pval, pval_CPT=pval_CPT, ess=ess, comp_time=comp_time, comp_time_CPT=comp_time_CPT)
}



n = 50
p = 20
j = 1 # our target is the first covariate
num.sim = 500 # number of repeated experiments

### parameters for the response model
eff_size = 1

### parameters for the efficient method
model.X.known = T; n_unlabel = n*10
group.method = 'kmeans'
k.list = c(5,10,15,20,25)
prob.list = c(1)

### parameters for CPT
nstep = 50
M_CPT = 500
M_efficient = M_CPT 

signal.list = seq(0, 1, 0.2)


### Single normal
### parameters for the covariate model
X.type = 'AR'; ar_coeff = 0.3; sd.a = NULL; mix.mean.sd = NULL
prob.mixture = NULL
sim1 = compare_efficient_CPT(n, p, M_CPT, M_efficient,
                             X.type, ar_coeff, prob.mixture, sd.a, mix.mean.sd,
                             signal.list, k.list, prob.list, nstep, num.sim)
saveRDS(sim1, file='sim1.rds')

X.type = 'AR'; ar_coeff = 0.5; sd.a = NULL; mix.mean.sd = NULL
prob.mixture = NULL
sim2 = compare_efficient_CPT(n, p, M_CPT, M_efficient,
                             X.type, ar_coeff, prob.mixture, sd.a, mix.mean.sd,
                             signal.list, k.list, prob.list, nstep, num.sim)
saveRDS(sim2, file='sim2.rds')

X.type = 'AR'; ar_coeff = 0.7; sd.a = NULL; mix.mean.sd = NULL
prob.mixture = NULL
sim3 = compare_efficient_CPT(n, p, M_CPT, M_efficient,
                             X.type, ar_coeff, prob.mixture, sd.a, mix.mean.sd,
                             signal.list, k.list, prob.list, nstep, num.sim)
saveRDS(sim3, file='sim3.rds')


### mixture of normals
X.type = 'AR'; ar_coeff = 0.3; sd.a = NULL; mix.mean.sd = 3
prob.mixture = rep(1,4)/4
sim4 = compare_efficient_CPT(n, p, M_CPT, M_efficient, 
                             X.type, ar_coeff, prob.mixture, sd.a, mix.mean.sd,
                             signal.list, k.list, prob.list, nstep, num.sim)
saveRDS(sim4, file='sim4.rds')

X.type = 'AR'; ar_coeff = 0.5; sd.a = NULL; mix.mean.sd = 3
prob.mixture = rep(1,4)/4
sim5 = compare_efficient_CPT(n, p, M_CPT, M_efficient, 
                             X.type, ar_coeff, prob.mixture, sd.a, mix.mean.sd,
                             signal.list, k.list, prob.list, nstep, num.sim)
saveRDS(sim5, file='sim5.rds')

X.type = 'AR'; ar_coeff = 0.7; sd.a = NULL; mix.mean.sd = 3
prob.mixture = rep(1,4)/4
sim6 = compare_efficient_CPT(n, p, M_CPT, M_efficient, 
                             X.type, ar_coeff, prob.mixture, sd.a, mix.mean.sd,
                             signal.list, k.list, prob.list, nstep, num.sim)
saveRDS(sim6, file='sim6.rds')


### linear
X.type = 'linear'; ar_coeff = NULL; sd.a = 0.5; mix.mean.sd = NULL
prob.mixture = NULL
sim7 = compare_efficient_CPT(n, p, M_CPT, M_efficient, 
                             X.type, ar_coeff, prob.mixture, sd.a, mix.mean.sd,
                             signal.list, k.list, prob.list, nstep, num.sim)
saveRDS(sim7, file='sim7.rds')

X.type = 'linear'; ar_coeff = NULL; sd.a = 1; mix.mean.sd = NULL
prob.mixture = NULL
sim8 = compare_efficient_CPT(n, p, M_CPT, M_efficient, 
                             X.type, ar_coeff, prob.mixture, sd.a, mix.mean.sd,
                             signal.list, k.list, prob.list, nstep, num.sim)
saveRDS(sim8, file='sim8.rds')


### mixture of normals
X.type = 'AR'; ar_coeff = 0.3; sd.a = NULL; mix.mean.sd = 1
prob.mixture = rep(1,4)/4
sim9 = compare_efficient_CPT(n, p, M_CPT, M_efficient, 
                             X.type, ar_coeff, prob.mixture, sd.a, mix.mean.sd,
                             signal.list, k.list, prob.list, nstep, num.sim)
saveRDS(sim9, file='sim9.rds')

X.type = 'AR'; ar_coeff = 0.5; sd.a = NULL; mix.mean.sd = 1
prob.mixture = rep(1,4)/4
sim10 = compare_efficient_CPT(n, p, M_CPT, M_efficient, 
                              X.type, ar_coeff, prob.mixture, sd.a, mix.mean.sd,
                              signal.list, k.list, prob.list, nstep, num.sim)
saveRDS(sim10, file='sim10.rds')

X.type = 'AR'; ar_coeff = 0.7; sd.a = NULL; mix.mean.sd = 1
prob.mixture = rep(1,4)/4
sim11 = compare_efficient_CPT(n, p, M_CPT, M_efficient, 
                              X.type, ar_coeff, prob.mixture, sd.a, mix.mean.sd,
                              signal.list, k.list, prob.list, nstep, num.sim)
saveRDS(sim11, file='sim11.rds')


### Equicorrelations
X.type = 'Equicorr'; ar_coeff = 0.15; sd.a = NULL; mix.mean.sd = NULL
prob.mixture = NULL
sim12 = compare_efficient_CPT(n, p, M_CPT, M_efficient, 
                              X.type, ar_coeff, prob.mixture, sd.a, mix.mean.sd,
                              signal.list, k.list, prob.list, nstep, num.sim)
saveRDS(sim12, file='sim12.rds')

X.type = 'Equicorr'; ar_coeff = 0.3; sd.a = NULL; mix.mean.sd = NULL
prob.mixture = NULL
sim13 = compare_efficient_CPT(n, p, M_CPT, M_efficient, 
                              X.type, ar_coeff, prob.mixture, sd.a, mix.mean.sd,
                              signal.list, k.list, prob.list, nstep, num.sim)
saveRDS(sim13, file='sim13.rds')

X.type = 'Equicorr'; ar_coeff = 0.45; sd.a = NULL; mix.mean.sd = NULL
prob.mixture = NULL
sim14 = compare_efficient_CPT(n, p, M_CPT, M_efficient, 
                              X.type, ar_coeff, prob.mixture, sd.a, mix.mean.sd,
                              signal.list, k.list, prob.list, nstep, num.sim)
saveRDS(sim14, file='sim14.rds')





## Robustness to model misspecification

compare_efficient_CPT_robust = function(n, p, M_CPT, M_efficient, 
                                        X.type, prob.mixture, ar_coeff, sd.a, mix.mean.sd,
                                        mis.specify.list, k.list, prob.list, nstep, num.sim, test.statistic='cor.X.res',
                                        cpt.run=TRUE, alpha=0.05, legend.bool=TRUE){
  
  pval = array(0,dim=c(length(mis.specify.list),length(k.list),length(prob.list),num.sim))
  ess = array(0,dim=c(length(mis.specify.list),length(k.list),length(prob.list),num.sim))
  comp_time = array(0,dim=c(length(mis.specify.list),length(k.list),length(prob.list),num.sim))
  pval_CPT = array(0,dim=c(length(mis.specify.list),num.sim))
  comp_time_CPT = array(0,dim=c(length(mis.specify.list),num.sim))
  
  
  for(mis.specify.i in 1:length(mis.specify.list)){
    mis.specify = mis.specify.list[mis.specify.i]
    message('\nMisspecification = ', mis.specify)
    for(num in 1:num.sim){
      if(num%%10==0){
        message(num,' of ', num.sim)
      }
      
      XY = generate_XY(n, p, X.type, eff_size, 0, ar_coeff, prob.mixture, sd.a, mix.mean.sd, mis.specify)
      X = XY$X; Y = XY$Y
      cond = Creat_condition_gaussian(XY, j, X.type, TRUE, prob.mixture)
      cond$mean_x = cond$mean_x * (1+mis.specify)
      # cond$sigma2_x = cond$sigma2_x * ((1-mis.specify)^2+1)/2
      res = lm(Y~X[,-j])$residuals
      
      ### efficient sampler + weighting
      for(k.i in 1:length(k.list)){
        for(prob.i in 1:length(prob.list)){
          k = k.list[k.i]
          prob = prob.list[prob.i]
          g = create_groups(cond$mean_x, k, 'kmeans')
          
          result = weighted_CPT(Y, X, j, prob, g, M_efficient, cond, test.statistic, res)
          ess[mis.specify.i, k.i, prob.i, num] = result$ess
          pval[mis.specify.i, k.i, prob.i, num] = result$pval
          comp_time[mis.specify.i, k.i, prob.i, num] = result$comp_time
        }
      }
      ### original CPT
      if(cpt.run){
        result_CPT = CPT(Y, X, j, nstep, M_CPT, cond, test.statistic, res)
        pval_CPT[mis.specify.i, num] = result_CPT$pval
        comp_time_CPT[mis.specify.i, num] = result_CPT$comp_time
      }
    }
  }
  
  plot_rejection_rate(pval, pval_CPT, mis.specify.list, 'Misspecification', cpt.run, alpha, legend.bool)
  
  ess_k_prob = apply(ess, c(2,3), mean)
  dimnames(ess_k_prob) = list(paste0('k=',k.list), paste0('prob=',prob.list))
  rej = apply(pval, c(1,2,3), function(x){mean(x<=alpha)})
  dimnames(rej) = list(paste0('Misspec=',mis.specify.list), paste0('k=',k.list), paste0('prob=',prob.list))
  
  cat(paste('\n\nComputational time of efficient method:',round(mean(comp_time),2)))
  if(cpt.run){
    cat(paste('\nComputational time of original CPT:',round(mean(comp_time_CPT),2)))
  }
  
  cat(paste('\n\nEffective sample size:\n'))
  print(round(ess_k_prob,2))
  cat('\n\nRejection rate of efficient method:\n')
  print(rej)
  cat('\n\nRejection rate of CPT:\n')
  print(rowMeans(pval_CPT <= alpha))
  
  list(pval=pval, pval_CPT=pval_CPT, ess=ess, comp_time=comp_time, comp_time_CPT=comp_time_CPT)
}



n = 50
p = 20
j = 1 # our target is the first covariate
num.sim = 500 # number of repeated experiments

### parameters for the response model
eff_size = 1

### parameters for the efficient method
model.X.known = T; n_unlabel = n*10
group.method = 'kmeans'
k.list = c(5,10,15,20,25)
prob.list = c(1)

### parameters for CPT
nstep = 50
M_CPT = 500
M_efficient = M_CPT #*30

mis.specify.list = seq(0,1,0.2)

### Single normal
X.type = 'AR'; ar_coeff = 0.3; sd.a = NULL; mix.mean.sd = NULL
prob.mixture = NULL
sim.robust.1 = compare_efficient_CPT_robust(n, p, M_CPT, M_efficient, 
                                            X.type, prob.mixture, ar_coeff, sd.a, mix.mean.sd,
                                            mis.specify.list, k.list, prob.list, nstep, num.sim, 
                                            cpt.run=TRUE, legend.bool=TRUE)
saveRDS(sim.robust.1, file='sim_robust_1.rds')

X.type = 'AR'; ar_coeff = 0.5; sd.a = NULL; mix.mean.sd = NULL
prob.mixture = NULL
sim.robust.2 = compare_efficient_CPT_robust(n, p, M_CPT, M_efficient, 
                                            X.type, prob.mixture, ar_coeff, sd.a, mix.mean.sd,
                                            mis.specify.list, k.list, prob.list, nstep, num.sim, 
                                            cpt.run=TRUE, legend.bool=TRUE)
saveRDS(sim.robust.2, file='sim_robust_2.rds')

X.type = 'AR'; ar_coeff = 0.7; sd.a = NULL; mix.mean.sd = NULL
prob.mixture = NULL
sim.robust.3 = compare_efficient_CPT_robust(n, p, M_CPT, M_efficient, 
                                            X.type, prob.mixture, ar_coeff, sd.a, mix.mean.sd,
                                            mis.specify.list, k.list, prob.list, nstep, num.sim, 
                                            cpt.run=TRUE, legend.bool=TRUE)
saveRDS(sim.robust.3, file='sim_robust_3.rds')


### mixture of normal
X.type = 'AR'; ar_coeff = 0.3; sd.a = NULL; mix.mean.sd = 1
prob.mixture = rep(1,4)/4
sim.robust.4 = compare_efficient_CPT_robust(n, p, M_CPT, M_efficient, 
                                            X.type, prob.mixture, ar_coeff, sd.a, mix.mean.sd,
                                            mis.specify.list, k.list, prob.list, nstep, num.sim, 
                                            cpt.run=TRUE, legend.bool=TRUE)
saveRDS(sim.robust.4, file='sim_robust_4.rds')

X.type = 'AR'; ar_coeff = 0.5; sd.a = NULL; mix.mean.sd = 1
prob.mixture = rep(1,4)/4
sim.robust.5 = compare_efficient_CPT_robust(n, p, M_CPT, M_efficient, 
                                            X.type, prob.mixture, ar_coeff, sd.a, mix.mean.sd,
                                            mis.specify.list, k.list, prob.list, nstep, num.sim, 
                                            cpt.run=TRUE, legend.bool=TRUE)
saveRDS(sim.robust.5, file='sim_robust_5.rds')

X.type = 'AR'; ar_coeff = 0.7; sd.a = NULL; mix.mean.sd = 1
prob.mixture = rep(1,4)/4
sim.robust.6 = compare_efficient_CPT_robust(n, p, M_CPT, M_efficient, 
                                            X.type, prob.mixture, ar_coeff, sd.a, mix.mean.sd,
                                            mis.specify.list, k.list, prob.list, nstep, num.sim, 
                                            cpt.run=TRUE, legend.bool=TRUE)
saveRDS(sim.robust.6, file='sim_robust_6.rds')


### linear
X.type = 'linear'; ar_coeff = NULL; sd.a = 0.5; mix.mean.sd = NULL
prob.mixture = NULL
sim.robust.linear1 = compare_efficient_CPT_robust(n, p, M_CPT, M_efficient, 
                                                  X.type, prob.mixture, ar_coeff, sd.a, mix.mean.sd,
                                                  mis.specify.list, k.list, prob.list, nstep, num.sim, 
                                                  cpt.run=TRUE, legend.bool=TRUE)
saveRDS(sim.robust.linear1, file='sim_robust_linear1_more.rds')

X.type = 'linear'; ar_coeff = NULL; sd.a = 1; mix.mean.sd = NULL
prob.mixture = NULL
sim.robust.linear2 = compare_efficient_CPT_robust(n, p, M_CPT, M_efficient,
                                                  X.type, prob.mixture, ar_coeff, sd.a, mix.mean.sd,
                                                  mis.specify.list, k.list, prob.list, nstep, num.sim,
                                                  cpt.run=TRUE, legend.bool=TRUE)
saveRDS(sim.robust.linear2, file='sim_robust_linear2.rds')




## Robustness to covariate model estimated on unlabeled data

compare_efficient_CPT_estimate = function(n, p, M_CPT, M_efficient, 
                                          X.type, prob.mixture, ar_coeff, sd.a, mix.mean.sd,
                                          n.unlabel.list, k.list, prob.list, nstep, num.sim, test.statistic='cor.X.res',
                                          cpt.run=TRUE, alpha=0.05, legend.bool=TRUE){
  
  pval = array(0,dim=c(length(n.unlabel.list),length(k.list),length(prob.list),num.sim))
  ess = array(0,dim=c(length(n.unlabel.list),length(k.list),length(prob.list),num.sim))
  comp_time = array(0,dim=c(length(n.unlabel.list),length(k.list),length(prob.list),num.sim))
  pval_CPT = array(0,dim=c(length(n.unlabel.list),num.sim))
  comp_time_CPT = array(0,dim=c(length(n.unlabel.list),num.sim))
  
  
  for(n.unlabel.i in 1:length(n.unlabel.list)){
    n.unlabel = n.unlabel.list[n.unlabel.i]
    message('\nN_unlabel = ', n.unlabel)
    for(num in 1:num.sim){
      if(num%%10==0){
        message(num,' of ', num.sim)
      }
      
      XY = generate_XY(n, p, X.type, eff_size, 0, ar_coeff, prob.mixture, sd.a, mix.mean.sd, n_unlabel = n.unlabel)
      X = XY$X; Y = XY$Y
      X_unlabel = list(X=XY$X_unlabel, cluster=XY$cluster_true_unlabel)
      cond = Creat_condition_gaussian(XY, j, X.type, FALSE, prob.mixture, X_unlabel)
      res = lm(Y~X[,-j])$residuals
      
      ### efficient sampler + weighting
      for(k.i in 1:length(k.list)){
        for(prob.i in 1:length(prob.list)){
          k = k.list[k.i]
          prob = prob.list[prob.i]
          g = create_groups(cond$mean_x, k, 'kmeans')
          
          result = weighted_CPT(Y, X, j, prob, g, M_efficient, cond, test.statistic, res)
          ess[n.unlabel.i, k.i, prob.i, num] = result$ess
          pval[n.unlabel.i, k.i, prob.i, num] = result$pval
          comp_time[n.unlabel.i, k.i, prob.i, num] = result$comp_time
        }
      }
      ### original CPT
      if(cpt.run){
        result_CPT = CPT(Y, X, j, nstep, M_CPT, cond, test.statistic, res)
        pval_CPT[n.unlabel.i, num] = result_CPT$pval
        comp_time_CPT[n.unlabel.i, num] = result_CPT$comp_time
      }
    }
  }
  
  plot_rejection_rate(pval, pval_CPT, n.unlabel.list, 'N_unlabel', cpt.run, alpha, legend.bool)
  
  ess_k_prob = apply(ess, c(2,3), mean)
  dimnames(ess_k_prob) = list(paste0('k=',k.list), paste0('prob=',prob.list))
  rej = apply(pval, c(1,2,3), function(x){mean(x<=alpha)})
  dimnames(rej) = list(paste0('N_unlabel=',n.unlabel.list), paste0('k=',k.list), paste0('prob=',prob.list))
  
  cat(paste('\n\nComputational time of efficient method:',round(mean(comp_time),2)))
  if(cpt.run){
    cat(paste('\nComputational time of original CPT:',round(mean(comp_time_CPT),2)))
  }
  
  cat(paste('\n\nEffective sample size:\n'))
  print(round(ess_k_prob,2))
  cat('\n\nRejection rate of efficient method:\n')
  print(rej)
  if(cpt.run){
    cat('\n\nRejection rate of CPT:\n')
    print(rowMeans(pval_CPT <= alpha))
  }
  
  list(pval=pval, pval_CPT=pval_CPT, ess=ess, comp_time=comp_time, comp_time_CPT=comp_time_CPT)
}

n.unlabel.list = c(50, 100, 150, 200)

X.type = 'linear'; ar_coeff = NULL; sd.a = 0.5; mix.mean.sd = NULL
prob.mixture = NULL
sim.estimate.linear1.more = compare_efficient_CPT_estimate(n, p, M_CPT, M_efficient, 
                                                           X.type, prob.mixture, ar_coeff, sd.a, mix.mean.sd,
                                                           n.unlabel.list, k.list, prob.list, nstep, num.sim, test.statistic='cor.X.res',
                                                           cpt.run=TRUE, alpha=0.05, legend.bool=TRUE)
saveRDS(sim.estimate.linear1.more, file='sim_estimate_linear1_more.rds')

