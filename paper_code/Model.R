# Exploratory analysis to determine max rank (K)

# [ find max rank: use single intercept model ]
# first focus on surrogate "hot zones" where crimes > 1

covars1 <- covars[covars$crime > 1, ]

rem <- setdiff(V(gl)$id, covars1$seg_id) # to be removed
rem2 <- which(V(gl)$id %in% rem == TRUE)

mL1 <- laplacian_marginal(gl, rem2, weights = sim)
e.mL1 <- eigen(mL1, symmetric = TRUE)

m <- nrow(e.mL1$vectors)
e.values1 <- e.mL1$values[m:1] # reversing eigenvalues
e.vectors1 <- e.mL1$vectors[, m:1] # reversing eigenvectors


preds1 <- preds[covars$crime > 1, 1, drop = FALSE] # intercept-only
y1 <- y[covars$crime > 1]

K_values <- seq(100, 500, by = 100)
lambda_values <- c(5, 10, 20, 30, 40)
pl <- matrix(nrow = length(K_values), ncol = length(lambda_values))
rownames(pl) <- paste0("K", K_values)
colnames(pl) <- paste0("l", lambda_values)

for (i in seq_along(K_values)) {
  K <- K_values[i]
  X <- expand_design(preds1, e.vectors1, rank = K)
  log_tau_p <- c(-Inf, rep(0, K)) # non-informative
  for (j in seq_along(lambda_values)) {
    lambda <- lambda_values[j]
    res <- em_network_rank(y1, X, mL1, lambda, log_tau_p,
                           control = list(maxit = 50, epsilon = 1e-4,
                                          trace = TRUE))
    tau <- centroid_network_rank(exp(res$logptau))
    X1 <- expand_design(preds1, e.vectors1, rank = tau)
    fit <- try(glm(y1 ~ X1 - 1, family = res$family,
                   method = make_fitter(lambda, X1, mL1)))
    if (class(fit)[1] != "try-error") pl[i, j] <- sum(press_loop(fit))
    message("[", i, ", ", j, "] K = ", K, ", l = ", lambda,
            " | press loop = ", round(pl[i, j], 2))
  }
}

# NOTES
# - lowest press-loop comes from larger K and larger lambda (diagonal effect)
# - press-loop seems to saturate for fixed K and as lambda increases

# given that PL is not perfect and adding parsimony, here's a reasonable value:

K <- 200

# [ base model without hot zones ]
rank <- c(rep(K, ncol(preds) - 2), 1, 10) #can change this to reflect max rank of each pred
X <- expand_design(preds, e.vectors, rank = rank)

#log_tau_p <- c(-Inf, rep(0, K)) # non-informative

# [hyperpriors, following section 3 of paper]
perc_var_target <- .9
inveigs <- 1 / e.values
perc_var <- cumsum(inveigs[2:K]) / sum(inveigs[2:K])
(tau_bar <- which.min(abs(perc_var - perc_var_target)))
plot(perc_var[1:K], xlab = "Number of Eigenvalues",
     ylab = "Partial Sum / Total Sum",
     main = "Variation Described by Eigenvalues (1/e)")
abline(h = perc_var_target); abline(v = tau_bar)

# at least one expansion (include pred), good chance of selecting higher ranks:
alpha0 <- 1; alpha1 <- .99
log_tau_p <- log(tau_prior(K, tau_bar, alpha0, alpha1))


lambda_values <- c(5, 10, 20, 30, 40) #can change based on prior exploratory analysis
pl <- rep(NA, length(lambda_values)); names(pl) <- paste0("l", lambda_values)
for (j in seq_along(lambda_values)) {
  lambda <- lambda_values[j]
  res <- em_network_rank(y, X, mL, lambda, log_tau_p,
                         control = list(epsilon = 1e-4, trace = TRUE))
  tau <- centroid_network_rank(exp(res$logptau))
  X1 <- expand_design(preds, e.vectors, rank = tau)
  fit <- try(glm(y ~ X1 - 1, family = res$family,
                 method = make_fitter(lambda, X1, mL)))
  if (class(fit)[1] != "try-error") pl[j] <- sum(press_loop(fit))
  message("[", j, "] l = ", lambda, " | press loop = ", round(pl[j], 2))
}


#lambda_nohz <- lambda_values[which.min(pl)] 
lambda_nohz <- 5

res_nohz <- em_network_rank(y, X, mL, lambda_nohz, log_tau_p,
                            control = list(epsilon = 1e-4, trace = TRUE))
tau <- centroid_network_rank(exp(res_nohz$logptau))
message_vec("tau", tau)
X1 <- expand_design(preds, e.vectors, rank = tau)
fit_nohz <- try(glm(y ~ X1 - 1, family = res_nohz$family,
                    method = make_fitter(lambda_nohz, X1, mL),
                    control = list(trace = TRUE)))
op <- par(mfrow = c(1, 2))
plot(y, fitted(fit_nohz)); abline(0, 1, col = "red")
mu <- fitted(fit_nohz); r <- (y - mu) / sqrt(mu)
plot(mu, r); lines(lowess(mu, r), col = "red", lwd = 2)
par(op)


# [ full model with hot zones ]

rank <- c(rep(K, ncol(preds) - 2), 1, 10)
X <- X_hz <- expand_design(preds, e.vectors, rank = rank)
X_bg <- X[, c(1, K * (ncol(preds) - 2) + 1)] # simple, no network expansion. can explore/change
p_bg_init <- as.integer(y <= 1) # simple threshold

l_values <- c(5, 10, 20, 30, 40)
pl <- matrix(nrow = length(l_values), ncol = length(l_values))
for (i in seq_along(l_values)) {
  for (j in seq_along(l_values)) {
    lambda <- l_values[i]; lambda_hz <- l_values[j]
    res <- em_network_rank_bg(y, X, X_bg, X, mL, lambda, lambda_hz,
                              log_tau_p, p_bg_init,
                              control = list(maxit = 20, epsilon = 1e-4, trace = TRUE))
    
    tau <- centroid_network_rank(exp(res$logptau))
    X1 <- X[, get_index_rank(attr(X, "rank"), tau), drop = FALSE]
    tau_hz <- centroid_network_rank(exp(res$logptau_hz))
    X1_hz <- X_hz[, get_index_rank(attr(X_hz, "rank"), tau_hz), drop = FALSE]

    fit <- em_network_bg(y, X1, X_bg, X1_hz, mL, lambda, lambda_hz, p_bg_init,
                         control = list(maxit = 50, epsilon = 1e-4, trace = TRUE))
    plot_em_fit(y, X1, X_bg, fit)
    pl[i, j] <- fit$pl
    message(">>> [", i, ", ", j, "] PL = ", round(pl[i, j], 2))
  }
}

lambda <- 10; lambda_hz <- 40 

# Notes:
# - As regularization penalty increases, the rank of predictors increases (as opposed to the intercept)

#And, putting it all together:

res_em <- em_network_rank_bg(y, X, X_bg, X, mL, lambda, lambda_hz,
                             log_tau_p, p_bg_init,
                             control = list(maxit = 50, epsilon = 1e-4, trace = TRUE))

tau <- centroid_network_rank(exp(res_em$logptau))
X1 <- X[, get_index_rank(attr(X, "rank"), tau), drop = FALSE]
tau_hz <- centroid_network_rank(exp(res_em$logptau_hz))
X1_hz <- X_hz[, get_index_rank(attr(X_hz, "rank"), tau_hz), drop = FALSE]

fit_em <- em_network_bg(y, X1, X_bg, X1_hz, mL, lambda, lambda_hz, p_bg_init,
                        control = list(maxit = 50, epsilon = 1e-4, trace = TRUE))
plot_em_fit(y, X1, X_bg, fit_em)
message("PL = ", fit_em$pl)


#save(K, lambda, lambda_hz,
#     res_nohz, fit_nohz, res_em, fit_em, file="boston-res.RData")


# [ HMC Gibbs sampler ]
K <- 200; lambda <- 40; lambda_hz <- 10

rank <- c(rep(K, ncol(preds) - 2), 1, 1)
X <- X_hz <- expand_design(preds, e.vectors, rank = rank)
X_bg <- X[, c(1, K * (ncol(preds) - 2) + 1)] # non-expanded

tau <- centroid_network_rank(exp(res_em$logptau))
X <- X[, get_index_rank(attr(X, "rank"), tau), drop = FALSE]
tau_hz <- centroid_network_rank(exp(res_em$logptau_hz))
X_hz <- X_hz[, get_index_rank(attr(X_hz, "rank"), tau_hz), drop = FALSE]


ns <- 1000 # final number of samples
burn_in <- .2; thin <- 10
control <- list(nchains = 4, nsamples = thin * ns / (1 - burn_in), trace = TRUE, epsilon = 1e-4)
simplified <- TRUE
L <- mL 

n <- nrow(X) 
nt <- ncol(X); nz <- ncol(X_bg); nw <- ncol(X_hz)
step_size <- max(nt, nz, nw) ^ (- 1. / 3)
sims <- bsglm::mcmc_init_array(control$nsamples, control$nchains,
                               c(colnames(X), paste0(colnames(X_bg), "_bg"),
                                 paste0(colnames(X_hz), "_hz"),
                                 paste0("Z", 1:nrow(X)), "__lq"))

theta_pc <- list(mean = rep(0, ncol(X)),
                 precision = lambda * crossprod(X, as.matrix(L %*% X)))
sampler_theta <- bsglm::sampler(prior_coef = theta_pc,
                                step_size = step_size,
                                simplified = simplified)

omega_pc <- list(mean = rep(0, ncol(X_hz)),
                 precision = lambda_hz * crossprod(X_hz, as.matrix(L %*% X_hz)))
sampler_omega <- bsglm::sampler(prior_coef = omega_pc,
                                step_size = step_size,
                                simplified = simplified)

sampler_zeta <- bsglm::sampler(step_size = step_size,
                               simplified = simplified)
one_sample <- list(nchains = 1, nsamples = 1)

for (chain in 1:control$nchains) {
  # init chain at EM MAP:
  theta <- fit_em$coef; zeta <- fit_em$coef_bg; omega <- fit_em$coef_hz #can change init values to explore
 
  for (iter in 1:control$nsamples) {
    lp <- 0
    # [ Z | theta, omega, zeta, y ]
    eta <- drop(X %*% theta); eta_hz <- drop(X_hz %*% omega)
    eta_bg <- drop(X_bg %*% zeta)
    p_bg <- plogis(eta_hz + y * eta_bg - exp(eta_bg) - (y * eta - exp(eta)))
    # p_bg in (eps, 1 - eps), avoid infs:
    p_bg <- pmin(pmax(p_bg, control$epsilon), 1 - control$epsilon)
    Z <- rbinom(n, 1, p_bg)
    bg_vertices <- Z == 1
    # [ theta | Z, zeta, omega, y ]
    ft <- glm(y[!bg_vertices] ~ X[!bg_vertices, ] - 1, family = poisson,
              start = theta, method = sampler_theta, control = one_sample)
    theta <- ft$samples[1, 1, 1:nt]
    lp <- lp + ft$samples[1, 1, nt + 2]
    # [ zeta | Z, theta, omega, y ]
    fz <- glm(y[bg_vertices] ~ X_bg[bg_vertices, ] - 1, family = poisson,
              start = zeta, method = sampler_zeta, control = one_sample)
    zeta <- fz$samples[1, 1, 1:nz]
    lp <- lp + fz$samples[1, 1, nz + 2]
    # [ omega | Z, theta, zeta, y ]
    fw <- glm(Z ~ X_hz - 1, family = binomial, start = omega,
              method = sampler_omega, control = one_sample)
    omega <- fw$samples[1, 1, 1:nw]
    lp <- lp + fw$samples[1, 1, nw + 2]
    
    sims[iter, chain, 1:nt] <- theta
    sims[iter, chain, nt + 1:nz] <- zeta
    sims[iter, chain, nt + nz + 1:nw] <- omega
    sims[iter, chain, nt + nz + nw + 1:n] <- Z
    sims[iter, chain, nt + nz + nw + n + 1] <- lp
    if (control$trace) message("[", chain, ", ", iter, "] lp = ", lp)
    plot(sims[1:iter, chain, nt + nz + nw + n + 1],
         xlim = c(1, control$nsamples), xlab = "iter",
         ylab = "lp", type = "l", main = paste0("Chain ", chain))
  }
  #save(sims, control, K, lambda, lambda_hz, file = "boston-mcmc.RData")
}

#save(sims, control, K, lambda, lambda_hz, file = "boston-mcmc.RData")



library(bayesplot)
color_scheme_set("mix-blue-red")
mcmc_trace(sims, pars = "__lq")
mcmc_trace(sims[valid, , ], pars = "__lq")

pars <- c("intercept_1")
pars <- c("__lq")

sims1 <- sims[, 1, , drop = TRUE]
#mcmc_intervals(sims1, pars = pars)
mcmc_areas(sims1, pars = pars,
           prob = .8, # .5
           prob_outer = .99, # .9
           point_est = "mean")

#mcmc_hist(sims1, pars = pars)
#mcmc_hist_by_chain(sims, pars = pars)
#mcmc_dens(sims1, pars = pars)
mcmc_dens_overlay(sims, pars = pars)
mcmc_violin(sims, pars = pars)

library(rstan) # for `monitor`
rhat <- function (sims, ...) monitor(sims, print = FALSE, ...)[, "Rhat"]
neff_ratio <- function (sims, ...) monitor(sims, print = FALSE, ...)[, "n_eff"]
mcmc_rhat(rhat(sims))
mcmc_rhat(rhat(sims[valid, , ]))
#mcmc_neff(neff_ratio(sims), size = 2)
mcmc_acf(sims, pars = pars, lags = 20)
mcmc_acf(sims[valid, , ], pars = pars, lags = 20)


# [ Visualizations (summaries) ]
library(banner)
library(Matrix)
setwd("/Users/emu1/Dropbox/BANNER_git/banner_nov21")
#setwd("~/git/banner/work/data")
load("boston.RData")
load("boston-res-noexpbg.RData")
load("boston-mcmc.RData")

# [1] choosing tau based on a threshold (visualization of the cumulative posterior of tau for the different predictors)
rank <- c(rep(K, ncol(preds) - 2), 1, 1)
plot(1:max(rank), type = "n", xlab = "rank", ylab = "cumulative prob", ylim = c(0, 1),
     main = "Main effect (counts)")
jpos <- 0
for (j in seq_along(rank)) {
  lines(1:rank[j], cumsum(exp(res_em$logptau[jpos + 1:rank[j]])), col = j)
  jpos <- jpos + rank[j]
}
abline(h = 0.5, lty = 2)

plot(1:max(rank), type = "n", xlab = "rank", ylab = "cumulative prob", ylim = c(0, 1),
     main = "Latent effect (BG indicators)")
jpos <- 0
for (j in seq_along(rank)) {
  lines(1:rank[j], cumsum(exp(res_em$logptau_hz[jpos + 1:rank[j]])), col = j)
  jpos <- jpos + rank[j]
}
abline(h = 0.5, lty = 2)


# [2] posterior distribution of coefficients (selection of a few) from the Gibbs sampler with the mean and EM-MAP estimates overlaid
X <- X_hz <- expand_design(preds, e.vectors, rank = rank)
#X_bg <- X[, c(1:K, K * (ncol(preds) - 2) + 1)]
X_bg <- X[, c(1, K * (ncol(preds) - 2) + 1)] # non-expanded

tau <- centroid_network_rank(exp(res_em$logptau))
X <- X[, get_index_rank(attr(X, "rank"), tau), drop = FALSE]
tau_hz <- centroid_network_rank(exp(res_em$logptau_hz))
X_hz <- X_hz[, get_index_rank(attr(X_hz, "rank"), tau_hz), drop = FALSE]

burn_in <- .2; thin <- 10
n <- nrow(X) # == nrow(X_bg) == nrow(X_hz)
nt <- ncol(X); nz <- ncol(X_bg); nw <- ncol(X_hz)

# quick check:
p_bg <- colMeans(sims[, 1, nt + nz + nw + 1:n], na.rm = TRUE)
plot(fit_em$p_bg, p_bg) # good agreement

valid <- (control$nsamples * burn_in):control$nsamples # discard burn-in
valid <- valid[valid %% thin == 0] # thin

# linear predictor plots
THOLD <- 5
op <- par(mfrow = c(control$nchains, 3))
for (chain in 1:control$nchains) {
  theta_s <- tcrossprod(X, sims[valid, chain, 1:nt])
  theta_f <- X %*% fit_em$coef[1:nt]
  hist(theta_s, main = paste("chain", chain)); abline(v = mean(theta_f), col = "red", lwd = 2)
  
  zeta_s <- tcrossprod(X_bg, sims[valid, chain, nt + 1:nz])
  zeta_f <- X_bg %*% fit_em$coef_bg
  hist(zeta_s, main = paste("chain", chain)); abline(v = mean(zeta_f), col = "red", lwd = 2)
  
  omega_s <- tcrossprod(X_hz, sims[valid, chain, nt + nz + 1:nw])
  omega_f <- X_hz %*% fit_em$coef_hz
  hist(omega_s, main = paste("chain", chain)); abline(v = mean(omega_f), col = "red", lwd = 2)
}
par(op)

# averaged main effect (original scale)
op <- par(mfrow = c(control$nchains, length(tau)))
for (chain in 1:control$nchains) {
  ipos <- 0
  for (j in seq_along(tau)) {
    ij <- ipos + 1:tau[j]
    betaj_s <- colMeans(tcrossprod(e.vectors[, 1:tau[j]], sims[valid, chain, ij]))
    betaj_f <- mean(e.vectors[, 1:tau[j], drop = FALSE] %*% fit_em$coef[ij])
    hist(betaj_s, main = paste("chain", chain), xlab = paste0("avg beta_", j))
    abline(v = betaj_f, col = "red", lwd = 2)
    ipos <- ipos + tau[j]
  }
}
par(op)

# theta (linear predictor)
THOLD <- 5
op <- par(mfrow = c(control$nchains, length(tau)))
for (chain in 1:control$nchains) {
  ipos <- 0
  for (j in seq_along(tau)) {
    ij <- ipos + 1:tau[j]
    thetaj_s <- tcrossprod(X[, ij], sims[valid, chain, ij])
    thetaj_f <- X[, ij, drop = FALSE] %*% fit_em$coef[ij]
    hist(thetaj_s[abs(thetaj_s) < THOLD],
         main = paste("chain", chain), xlab = paste0("theta_", j))
    abline(v = mean(thetaj_f), col = "red", lwd = 2)
    ipos <- ipos + tau[j]
  }
}
par(op)


# zeta (linear predictor)
THOLD <- 10
op <- par(mfrow = c(control$nchains, ncol(X_bg)))
for (chain in 1:control$nchains) {
  for (j in seq_len(ncol(X_bg))) {
    zetaj_s <- tcrossprod(X_bg[, j], sims[valid, chain, nt + j])
    zetaj_f <- X_bg[, j, drop = FALSE] %*% fit_em$coef_bg[j]
    hist(zetaj_s[abs(zetaj_s) < THOLD],
         main = paste("chain", chain), xlab = paste0("zeta_", j))
    abline(v = mean(zetaj_f), col = "red", lwd = 2)
  }
}
par(op)


# averaged BG indicator effect (original scale)
op <- par(mfrow = c(control$nchains, length(tau_hz)))
for (chain in 1:control$nchains) {
  ipos <- 0
  for (j in seq_along(tau_hz)) {
    ij <- ipos + 1:tau_hz[j]
    gammaj_s <- colMeans(tcrossprod(e.vectors[, 1:tau_hz[j]], sims[valid, chain, nt + nz + ij]))
    gammaj_f <- mean(e.vectors[, 1:tau_hz[j], drop = FALSE] %*% fit_em$coef_hz[ij])
    hist(gammaj_s, main = paste("chain", chain), xlab = paste0("gamma_", j))
    abline(v = gammaj_f, col = "red", lwd = 2)
    ipos <- ipos + tau_hz[j]
  }
}
par(op)

# omega (linear predictor)
THOLD <- 5
op <- par(mfrow = c(control$nchains, length(tau_hz)))
for (chain in 1:control$nchains) {
  ipos <- 0
  for (j in seq_along(tau_hz)) {
    ij <- ipos + 1:tau_hz[j]
    omegaj_s <- tcrossprod(X_hz[, ij], sims[valid, chain, nt + nz + ij])
    omegaj_f <- X_hz[, ij, drop = FALSE] %*% fit_em$coef_hz[ij]
    hist(omegaj_s[abs(omegaj_s) < THOLD],
         main = paste("chain", chain), xlab = paste0("omega_", j))
    abline(v = mean(omegaj_f), col = "red", lwd = 2)
    ipos <- ipos + tau_hz[j]
  }
}
par(op)


# 3 - binned deviance residual plot
p_bg <- plogis(X_hz %*% fit_em$coef_hz)
mu <- drop(p_bg * exp(X_bg %*% fit_em$coef_bg) + (1 - p_bg) * exp(X %*% fit_em$coef))
nc <- ceiling(.05 * length(y) / 50) * 50 # 5% of obs, rounded up by 50

ii <- factor(cutree(hclust(dist(preds[, -1])), nc)) # cluster by predictors
wf <- tapply(y, ii, length)
yf <- tapply(y, ii, sum) / wf # average response

xf <- (1 / as.vector(wf)) * rowsum(X_hz, ii) # average preds (hz)
pf <- drop(plogis(xf %*% fit_em$coef_hz))
xf <- (1 / as.vector(wf)) * rowsum(X_bg, ii) # average preds (bg)
mf0 <- drop(exp(xf %*% fit_em$coef_bg))
xf <- (1 / as.vector(wf)) * rowsum(X, ii) # average preds (main)
mf1 <- drop(exp(xf %*% fit_em$coef))
mf <- pf * mf0 + (1 - pf) * mf1 # E[y]
sf <- sqrt(mf + pf * (1 - pf) * (mf0 - mf1) ^ 2) # SD[y]
rf <- sqrt(wf) * (yf - mf) / sf # Pearson residuals

op <- par(mfrow = c(1, 2))
plot(mf, rf, xlab = "pooled fitted", ylab = "pooled Pearson residuals")
abline(h = 0, lty = 2)
lines(lowess(mf, rf), col = "red")
qqnorm(rf); qqline(rf)
par(op)



# 4 - wealth effect (beta) varying over the city (or different coefficient varying over city), and probability of hot zone for each intersection
library(tmap)

ii <- c(1:tau[1], sum(tau[1:3]) + 1) # intercept + log_population
covars$intercept <- drop(e.vectors[, ii] %*% fit_em$coef[ii])
covars$logtax_effect <- drop(e.vectors[, 1:tau[2], drop = FALSE] %*%
                               fit_em$coef[tau[1] + 1:tau[2]])
covars$hz_prob <- drop(plogis(-X_hz %*% fit_em$coef_hz))
covars$hz_em_prob <- 1 - res_em$p_bg

tmap_arrange(
  tm_shape(covars) +
    tm_bubbles(col = "intercept", palette = "-RdYlBu", size = 0.2, style = "quantile"),
  tm_shape(covars) +
    tm_bubbles(col = "logtax_effect", palette = "Greens", size = 0.2, style = "quantile"),
  tm_shape(covars) +
    tm_bubbles(col = "hz_prob", palette = "Reds", size = 0.2, style = "quantile"),
  tm_shape(covars) +
    tm_bubbles(col = "hz_em_prob", palette = "Reds", size = 0.2, style = "quantile"),
  nrow = 2
)