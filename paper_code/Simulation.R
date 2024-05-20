#Run PrepData.R, but filter on Allston Neighborhood

library(tidyverse)
library(sp)
library(rgeos)
library(igraph)
library(Matrix)
library(banner)
library(MASS)
library(ggplot2)
library(tidyr)
library(gridExtra)

logit <- qlogis;  expit <- plogis

Incidence <- function(g)  {
  #This function computes the weighted incidence matrix of a graph.
  #Args: 
  #   g: an igraph with edge attributes (weight)
  #Returns:
  #   M: the corresponding weighted incidence matrix
  #######################################################################
  
  nrows <- length(E(g))
  ncols <- length(V(g))
  edgelist <- get.edges(g,seq(1, nrows, by=1))
  M <- matrix(rep(0,nrows*ncols),nrows,ncols)
  
  for(i in 1:nrows)   {
    a <- edgelist[i, 1]
    b <- edgelist[i, 2]
    M[i, a] <- 1*sqrt(E(g)[i]$weight)
    M[i, b] <- -1*sqrt(E(g)[i]$weight)
  }
  M
}


M <- Incidence(g)

size <- vcount(g)
lambda <- 2.477076 #flexible
n.e <- 20 #basis expansion- change this for exploratory 

rank <- c(rep(n.e, ncol(preds) - 1), 1) #based on int, wealth, pop, pol, college-b and n.e above.  This is flexible
X <- expand_design(preds, e.vectors, rank = rank)

qr.m <- qr(M %*% X)
R1 <- size/10*qr.R(qr.m) #cholesky factor of precision, R1-1. 
mu <- c(rep(0, dim(X)[2]))

r.rank <- qr.m$rank

#Choosing parameters that control basis rank expansion.  
inveigs <- 1/e.values
Kmax <- 100 
ratio <- cumsum(inveigs[2:Kmax])/sum(inveigs[2:Kmax])
diff <- abs(ratio-.5) #50%
tau_bar <- which.min(diff); 

tone <- tau_prior(K = Kmax, tau_bar = tau_bar, alpha_0 = 0.99, alpha_1 = 0.90) #flexible
log_tau_p <- log(tone)

ppl_p<- function (y, yhat) { # posterior predictive loss for Poisson
  sum((y - yhat) ^ 2 + yhat)
}

###, diagnostic two
ppl_em <- function (y, zeta, eta, pg) { # posterior predictive loss for EM
  muz <- exp(zeta); mue <- exp(eta)
  mu <- pg * muz + (1 - pg) * mue
  v <- mu + pg * (1 - pg) * (muz - mue) ^ 2
  sum((y - mu) ^ 2 + v)
}

preds2 <- preds[,c("intercept", "log_population")]

sim_mods <- function(type, mu, lambda, R1, g, e.vectors, diagnostic = ppl_p, ppl_em, X = X, preds2 = preds2){
  
  if(type == "simple"){
    Y <- rpois(size, 5)}
  
  if(type == "complex"){
  z <- rnorm(length(mu), mu, sd = 1/sqrt(lambda) )
  betas <- backsolve(R1, z)  + 3
  betas[1] <- .5 #FIXME? Not sure why this happens if fullrank

  #3 hotzones
  hotzone <- function(hzsize, g)  {
    initiate <- sample(seq(1, length(V(g)), by =1), size = 1)
    hzsg <- bfs(g, root = initiate, "all", order = TRUE)
    z <- hzsg$order[1:hzsize]
  }

  marker1 <- hotzone(0.1*vcount(g), g) #10% of nodes
  marker2 <- hotzone(0.1*vcount(g), g)
  marker3 <- hotzone(0.1*vcount(g), g) #identifies hot zones.  In the hot zone, z = 0.
  zind <- rep(1, length(V(g)))
  zind[c(marker1, marker2, marker3)] <- 0

  ####sampling from neg bin#############
  eta <- 2 #use to match distribution of counts on actual network
  mu <- exp(zind*(eta)+(1-zind)*(X%*%betas))
  Y <- rnegbin(length(mu),mu, theta = 1)
  }

  Y[which(Y > 30)] <- 15 #may need to control to match crime distribution more closely
  
  Kmax <- 100
  pl <- c()
  
  ###Find max basis expansion
  lambda_values <- 10
  for (j in 1:Kmax){  
    X <- e.vectors[,1:j]
    prec <- lambda_values * e.values[1:j]
    pc <- list(mean = rep(0, size = j), precision = prec)
    fit <- try(glm(Y ~ X - 1, family = poisson,
                   method = bsglm::fitter(prior_coef = pc)))
    if (class(fit)[1] != "try-error") pl[j] <- sum(press_loop(fit))
  }
  
  Kmax <- which.min(pl)
  
  #Model 1, kernel regression (intercept only)
  dim_mod1 <- Kmax
  X1 <- e.vectors[,1:Kmax]
  
  #find lambda
  l <- 10 #how many lambda values to use in search for best lambda 
  lambda_values <- 10^(seq(-4,4,length=l))
  pl <- c()
  
  for (j in seq_along(lambda_values)) { 
    prec <- lambda_values[j] * e.values[1:Kmax]
    pc <- list(mean = rep(0, Kmax), precision = prec)
    fit <- try(glm(Y ~ X1 - 1, family = poisson,
                   method = bsglm::fitter(prior_coef = pc)))
    if (class(fit)[1] != "try-error") pl[j] <- sum(press_loop(fit))
  }
  
  index <- which.min(pl)
  prec <- lambda_values[index] * e.values[1:Kmax]
  pc1 <- list(mean = rep(0, Kmax), precision = prec)
  model1 <- glm(Y ~ X1 - 1, family = poisson, method = bsglm::fitter(prior_coef = pc1))
  yhat1 <- exp(as.matrix(X1)%*%model1$coef)
  
  #Model 2- kernel regression with diffusion kernel
  #find lambda
  for (j in seq_along(lambda_values)) { 
    prec <- exp(-lambda_values[j] * e.values[1:Kmax])
    pc <- list(mean = rep(0, Kmax), precision = prec)
    fit <- try(glm(Y ~ X1 - 1, family = poisson,
                   method = bsglm::fitter(prior_coef = pc)))
    if (class(fit)[1] != "try-error") pl[j] <- sum(press_loop(fit))
  }
  
  index <- which.min(pl)
  prec <- exp(-lambda_values[index] * e.values[1:Kmax])
  pc2 <- list(mean = rep(0, Kmax), precision = prec)
  model2 <- glm(Y ~ X1 - 1, family = poisson, method = bsglm::fitter(prior_coef = pc2))
  yhat2 <- exp(as.matrix(X1)%*%model2$coef)
  
  #Model 3- regular Poisson
  model3 <- glm(Y~preds, family=poisson)
  yhat3 <- model3$fitted.values
  

  #Model 4- intercept with hotzones
  z <- c()
  z[which(yhat1<=quantile(yhat1)[4])] <- 1 
  z[which(yhat1>quantile(yhat1)[4])] <- 0
  pm.z1 <- glm(Y ~ 1, family = poisson, subset = (z==1)) #z = 1
  pm.z0 <- glm(Y ~ X1 - 1, family = poisson, subset = (z==0), method = bsglm::fitter(prior_coef = pc1))
  bern <- glm(z ~ X1 - 1, family = binomial, method = bsglm::fitter(prior_coef = pc1))
  
  etaknot <- pm.z1$coef ##poisson, z=1 (zero crime)
  beta <- pm.z0$coef ##poisson, z=0  (a NON zero crime event)
  gamma <- bern$coef ##bernoulli to calculate pzero
  
  maxit <- 100
  D <- 3000 
  tol <- 1e-4
  
  for (iter in 1:maxit){
    
    nlpi <- log(1+exp(-Y*(etaknot)+exp(etaknot)+as.matrix(X1)%*%beta*Y-exp(as.matrix(X1)%*%beta)-as.matrix(X1)%*%gamma))
    Pi <- exp(-nlpi)
    
    nzeros <- which(Pi!=1)
    
    mod1 <- glm(Pi ~ X1 - 1, subset = nzeros, family = quasibinomial, method = bsglm::fitter(prior_coef = pc1))
    gamma <- mod1$coef
    
    mod2 <- glm(Pi*Y ~ 1, family = quasipoisson, offset = log(Pi+.01))
    etaknot <- mod2$coef
    
    mod3 <- glm((1-Pi)*Y ~ X1 - 1, subset = nzeros, family = quasipoisson, method = bsglm::fitter(prior_coef = pc1), offset = log(1-Pi))
    beta <- mod3$coef
    
    X1 <- as.matrix(X1)
    phat <- expit(as.matrix(X1[nzeros,])%*%gamma)
    muhat <- exp(as.matrix(X1[nzeros,])%*%beta+log(1-Pi[nzeros]))
    
    D.new <- mod1$deviance + mod2$deviance + mod3$deviance
    
    if (abs(D.new-D)/abs(D) < tol)
      break
    
    D <- D.new
    
  }
  
  yhat4 <- Pi*(exp(etaknot))+(1-Pi)*exp(as.matrix(X1)%*%beta)
  eta <- as.matrix(X1)%*%beta
 

#diagnostics for Model 4 only (requires different diagnostic)  
  
  if(type == "simple"){
    res_mod4 <- ppl_em(y = Y, zeta = etaknot, eta = eta, pg = Pi)
  }
    
  if(type == "complex"){
    bg <- seq(1, length(Y))
    hz <- c(marker1, marker2, marker3)
    bg <- bg[-hz]
  
    hot_mod4 <- ppl_em(y = Y[hz], zeta = etaknot, eta = eta[hz], pg = Pi[hz])
    bg_mod4 <-  ppl_em(y = Y[bg], zeta = etaknot, eta = eta[bg], pg = Pi[bg])
  }
  
  #Model 5- Random intercept, set rank to be one for all preds
  rank <- c(Kmax, rep(1,ncol(preds) - 1))
  X5 <- expand_design(X = preds, Phi = e.vectors, rank = rank)
  
  pl <- rep(NA, length(lambda_values))
  
  for (j in seq_along(lambda_values)) { 
      prec <- lambda_values[j] * crossprod(X5, as.matrix(mL %*% X5))
      pc5 <- list(mean = rep(0, dim(X5)[2]), precision = prec)
      fit <- try(glm(Y ~ X5 - 1, family = poisson,
                     method = bsglm::fitter(prior_coef = pc5)))
      if (class(fit)[1] != "try-error") pl[j] <- sum(press_loop(fit))
    }

  index <- which.min(pl)
  prec <- lambda_values[index] * crossprod(X5, as.matrix(mL %*% X5))
  pc5 <- list(mean = rep(0, dim(X5)[2]), precision = prec)
  model5 <- glm(Y ~ X5 - 1, family = poisson, method = bsglm::fitter(prior_coef = pc5))
  yhat5 <- exp(X5%*%model5$coef)
  
  
  #Model 6- Poisson with basis expansion. Find basis expansion for each predictor 
  rank <- c(rep(Kmax, ncol(preds) - 1), 1)
  X6 <- expand_design(X = preds, Phi = e.vectors, rank = rank)
  
  pl <- rep(NA, length(lambda_values))
  
  for (j in seq_along(lambda_values)) { 
    res <- try(em_network_rank(y = Y, X = X6, L = mL, lambda = lambda_values[j], 
                               log_tau_p = log_tau_p,control = list(trace = FALSE), v0 = 0.2))
    if (class(res)[1] != "try-error") {
      tau <- centroid_network_rank(exp(res$logptau), kappa = 4)
      tau[5] <- 1
      Xnew <- expand_design(preds, e.vectors, rank = tau) 
      prec <- lambda_values[j] * crossprod(Xnew, as.matrix(mL %*% Xnew))
      pc <- list(mean = rep(0, sum(tau)), precision = prec)
      fit <- try(glm(Y ~ Xnew - 1, family = res$family,
                     method = bsglm::fitter(prior_coef = pc)))
      if (class(fit)[1] != "try-error") pl[j] <- sum(press_loop(fit))
    }
  }
  index <- which.min(pl)
  
  res <- em_network_rank(y = Y, X = X6, L = mL, lambda = lambda_values[index], 
                         log_tau_p = log_tau_p,control = list(trace = FALSE), v0 = 0.2) 
  tau <- centroid_network_rank(exp(res$logptau), kappa = 4) 
  Xnew <- expand_design(preds, e.vectors, rank = tau) 
  prec <- lambda_values[index] * crossprod(Xnew, as.matrix(mL %*% Xnew))
  pc6 <- list(mean = rep(0, sum(tau)), precision = prec)
  model6 <- glm(Y ~ Xnew - 1, family = res$family,method = bsglm::fitter(prior_coef = pc6))
  yhat6 <- model6$fitted.values
  
  #Model 7- ZIP model.  Use results from model 6 above. 
  z <- c()
  z[which(yhat6<=quantile(yhat6)[4])] <- 1 
  z[which(yhat6>quantile(yhat6)[4])] <- 0
  pm.z1 <- glm(Y ~ preds2 - 1, family = poisson, subset = (z==1)) #z = 1
  pm.z0 <- glm(Y ~ Xnew - 1, family = poisson, subset = (z==0), method = bsglm::fitter(prior_coef = pc6))
  bern <- glm(z ~ Xnew - 1, family = binomial, method = bsglm::fitter(prior_coef = pc6))
  
  etaknot <- pm.z1$coef ##poisson, z=1 (zero crime)
  beta <- pm.z0$coef ##poisson, z=0  (a NON zero crime event)
  gamma <- bern$coef ##bernoulli to calculate pzero
  
  for (iter in 1:maxit){
    
    nlpi <- log(1+exp(-Y*(preds2%*%etaknot)+exp(preds2%*%etaknot)+Xnew%*%beta*Y-exp(Xnew%*%beta)-Xnew%*%gamma))
    Pi <- exp(-nlpi)
    
    nzeros <- which(Pi!=1)
    
    mod1 <- glm(Pi ~ Xnew - 1, subset = nzeros, family = quasibinomial, method = bsglm::fitter(prior_coef = pc6))
    gamma <- mod1$coef
    
    mod2 <- glm(Pi*Y ~ preds2-1, family = quasipoisson, offset = log(Pi))
    etaknot <- mod2$coef
    
    mod3 <- glm((1-Pi)*Y ~ Xnew - 1, subset = nzeros, family = quasipoisson, method = bsglm::fitter(prior_coef = pc6), offset = log(1-Pi))
    beta <- mod3$coef
    
    phat <- expit(Xnew[nzeros,]%*%gamma)
    muhat <- exp(Xnew[nzeros,]%*%beta+log(1-Pi[nzeros]))
    
    D.new <- mod1$deviance + mod2$deviance + mod3$deviance
    
    if (abs(D.new-D)/abs(D) < tol)
      break
    
    D <- D.new
  }
  
  yhat7 <- Pi*(exp(preds2%*%etaknot))+(1-Pi)*exp(Xnew%*%beta)
  eta <- Xnew%*%beta
  zeta <- preds2%*%etaknot
  
  yhats <- cbind(yhat1, yhat2, yhat3, yhat5, yhat6)
  
  if(type == c("simple")){
    res <-  apply(yhats, 2, diagnostic, y = Y)
    res_mod7 <-  ppl_em(y = Y, zeta = zeta, eta = eta, pg = Pi)
    all <- c(res, res_mod4, res_mod7)
    }
 
  if(type == c("complex")){
    hot_mod <- apply(yhats[hz,], 2, diagnostic, y = Y[hz])
    bg_mod <-  apply(yhats[bg,], 2, diagnostic, y = Y[bg])
  
    hot_mod7 <-  ppl_em(y = Y[hz], zeta = zeta[hz], eta = eta[hz], pg = Pi[hz])
    bg_mod7 <-  ppl_em(y = Y[bg], zeta = zeta[bg], eta = eta[bg], pg = Pi[bg])
    all <- c(hot_mod, hot_mod4, hot_mod7, bg_mod, bg_mod4, bg_mod7)
   }
  
  dim <- c(dim_mod1, tau)
  report <- list(PPL = all,SIZE = dim)
}

#if type is simple, then PPL will have 7 entries per simulation (for the 7 models)
#if type is complex, then PPL will have 14 entries per simulation (two for each of the 7 models)

#simple
set.seed(10)
maxiter <- 200
results <- matrix(nrow = 7, ncol = 200)
dims <- matrix(nrow = 200, ncol = 6)
for (iter in 1:maxiter){
  print(iter)
  res <- try(sim_mods(type = c("simple"), mu = mu, lambda, R1, g, e.vectors, diagnostic = ppl_p, ppl_em, X, preds2 = preds2))
  if (class(res)!="try-error")  {
    results[,iter] <- res$PPL
    dims[iter,] <- res$SIZE
    
  }
}
#save.image(file = "SimpleSims.Rdata") 

#complex
set.seed(10)
maxiter <- 200
results2 <- matrix(nrow = 14, ncol = 200)
dims <- matrix(nrow = 200, ncol = 6)
for (iter in 1:maxiter){
  print(iter)
  res <- try(sim_mods(type = c("complex"), mu = mu, lambda, R1, g, e.vectors, diagnostic = ppl_p, ppl_em, X, preds2 = preds2))

  if (class(res)!="try-error")  {
    results2[,iter] <- res$PPL
    dims[iter,] <- res$SIZE
    
  }
}
#save.image(file = "ComplexSims.Rdata") 


#SimpleSims visualization [ Figure 5 ]
test <- data.frame(
  value = c(results),
  group = 1:7
)

p1<- ggplot(test, aes(factor(group), value)) + 
  geom_boxplot(color="#aa0022") +
  ggtitle("A") +
  scale_x_discrete() +
  labs(x=" ",y="") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14), 
        axis.text.x = element_blank())

#Complex sims visualization  [ Figure 5 ]
test2 <- data.frame(
  value = c(results2),
  group = 1:14
)

test3 <- test2 %>% filter(group < 8)

p2 <- ggplot(test3, aes(factor(group), value)) + 
  geom_boxplot(color="#aa0022") +
  ggtitle("B, Hot Zone") +
  scale_x_discrete() +
  labs(x= "", y ="PPL")+
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.text.x = element_blank())

test4 <- test2 %>% filter(group > 7)

p3 <- ggplot(test4, aes(factor(group), value)) + 
  geom_boxplot(color="#aa0022") +
  ggtitle("B, Background") +
  scale_x_discrete(labels=rep(c("Mod1", "Mod2","Mod3", "Mod4", "Mod5" , "Mod6", "Mod7"),2)) +
  labs(x= "", y ="")+
  theme(axis.title = element_text(size = 16),axis.text = element_text(size = 14))

grid.arrange(p1,p2,p3, heights = c(3, 3, 3))






