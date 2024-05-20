# banner: BAyesiaN NEtwork Regularization

# [ Helper functions ]

LOGEPS <- log(.Machine$double.eps)
lse <- function (x, y) {
  if (is.infinite(x)) return(y)
  if (is.infinite(y)) return(x)
  w <- max(x, y); d <- -abs(x - y)
  if (d < LOGEPS) w else w + log1p(exp(d))
}
log_sum_exp <- function (x) Reduce(lse, x, init = -Inf)

mat_mult <- function (A, b) if (length(b) == 1) A * b else A %*% b
mat_quad <- function (v, L) sum(v * as.matrix(L %*% v))
mat_quad_inv <- function (v, C) drop(crossprod(backsolve(C, v, transpose = TRUE)))
chol_solve <- function (C, v) backsolve(C, backsolve(C, v, transpose = TRUE))

check_canonical <- function (fam) {
  ((!is.null(fam$linkcan)) && fam$linkcan) ||
    ((fam$family == "poisson" || fam$family == "quasipoisson") &&
        fam$link == "log") ||
    ((fam$family == "binomial" || fam$family == "quasibinomial") &&
        fam$link == "logit") ||
    (fam$family == "gaussian" && fam$link == "identity")
}

# used in `em_network_rank`:
Etau <- function (v0, nu) function (l, m) {
  nu_min <- nu[pmin(l, m)]; nu_max <- nu[pmax(l, m)]
  nu_min * v0 + (nu_max - nu_min) * sqrt(v0) + 1 - nu_max
}


#' Create network regularized fitter.
#'
#' Create a regularized fitter based on a normal prior with zero mean and
#' precision \code{lambda * t(X) * L * X}.
#'
#' @param lambda precision scale or regularization penalty
#' @param X design matrix
#' @param L Laplacian matrix
#' @return fitter object from package \code{bsglm}
#' @export
make_fitter <- function (lambda, X, L) {
  prec <- lambda * crossprod(X, as.matrix(L %*% X))
  pc <- list(mean = rep(0, ncol(X)), precision = prec)
  bsglm::fitter(prior_coef = pc)
}


# [ Hyper-prior functions ]

#' Find hyper-prior network range parameter based on network distances.
#'
#' For the similarity score defined as \code{S = exp(-dist / phi)}, find
#' \code{phi} such that the median similarity score computed from distances
#' \code{dist} is \code{median_similarity}.
#'
#' @param dist distance matrix
#' @param median_similarity median similarity distance
#' @param interval root search interval, default (0, 1)
#' @return The hyper-prior network range parameter \code{phi}
#' @export
similarity_phi <- function (dist, median_similarity, interval = c(0, 1), ...) {
  min_med_range <- (median(dist) - min(dist)) / max(dist) # range in 0-1 scale
  uniroot(function (phi) exp(-min_med_range / phi) - median_similarity,
          interval, ...)$root
}

#' Laplacian marginalization
#'
#' Given a set \code{v} of vertices to be marginalized in a weighted
#' \code{graph} with edge \code{weights}, compute the Laplacian of
#' \code{G[-v]}.
#'
#' @param graph igraph object
#' @param v vertex set to be marginalized
#' @param weights edge weights
#' @return Laplacian matrix.
#' @export
laplacian_marginal <- function (graph, v, weights) {
  L <- graph.laplacian(graph, sparse = TRUE, weights = weights)
  C <- Matrix::Cholesky(L[v, v])
  L[-v, -v] - Matrix::crossprod(L[v, -v], Matrix::solve(C, L[v, -v]))
}


#' Prior for basis indicator \code{tau}, using hyper-parameter \code{rho}.
#'
#' Computes the prior rank probabilities \code{P(tau | K)} where \code{K} is
#' the maximum basis expansion based on parameters \code{alpha_0 = P(tau > 0)},
#' \code{alpha_1 = P(tau > 1 | tau > 0)}, and \code{rho}, where
#' \code{P(tau | tau > 1)} is proportional to \code{rho^{tau - 1}}.
#'
#' @param K maximum basis expansion
#' @param rho probability ratio \code{P(tau = k + 1 | tau > 1) / P(tau = k | tau > 1)}
#' @param alpha_0 probability of rank > 0
#' @return alpha_1 probability of rank > 1 given rank > 0
#' @return vector with probabilities \code{P(tau | K)}
#' @export
tau_prior_rho <- function (K, rho, alpha_0, alpha_1) {
  t <- rho ^ (1:(K - 1))
  c(1 - alpha_0, alpha_0 * (1 - alpha_1), alpha_0 * alpha_1 * t / sum(t))
}

#' Prior for basis indicator \code{tau}, using hyper-parameter \code{tau_bar}.
#'
#' Computes the prior rank probabilities \code{P(tau | K)} where \code{K} is
#' the maximum basis expansion based on parameters \code{alpha_0 = P(tau > 0)},
#' \code{alpha_1 = P(tau > 1 | tau > 0)}, and \code{rho}, where
#' \code{P(tau | tau > 1)} is proportional to \code{rho^{tau - 1}}.
#' Parameter \code{rho} is computed from \code{tau_bar}.
#'
#' @param K maximum basis expansion
#' @param tau_bar prior expected rank \code{E(tau)}
#' @param alpha_0 probability of rank > 0
#' @return alpha_1 probability of rank > 1 given rank > 0
#' @return vector with probabilities \code{P(tau | K)}
#' @export
tau_prior <- function (K, tau_bar, alpha_0, alpha_1, ...) {
  # find rho first:
  L <- (tau_bar - alpha_0) / (alpha_0 * alpha_1)
  f <- function (p) 1 / (1 - p) - (K - 1) * p ^ (K - 1) / (1 - p ^ (K - 1)) - L
  rho <- uniroot(f, c(0, 1 - 1e-6), ...)$root
  tau_prior_rho(K, rho, alpha_0, alpha_1)
}


# [ Cross-validation ]

#' Leave-One-Out-Proxy (LOOP) for PRESS score
#'
#' Given a model \code{fit}, compute PRESS LOOP score.
#'
#' @param fit model fit object
#' @param X design matrix to use in PRESS evaluation
#' @return PRESS score
#' @export
press_loop <- function (fit, X = model.matrix(fit)) {
  eta <- drop(X %*% fit$coef); W <- fit$family$mu.eta(eta)
  mu <- fit$family$linkinv(eta); V <- fit$family$variance(mu)
  if (!check_canonical(fit$family)) W <- (W / V) * W # adjust?
  h <- W * apply(X, 1, mat_quad_inv, fit$qr$qr) # leverages
  (fit$y - mu) ^ 2 / V / (1 - h) ^ 2
}


# [ EM for network basis rank ]

#' Network-expanded design.
#'
#' Given a base design matrix \code{X} and the Laplacian eigenvectors
#' \code{Phi}, computes the expanded design
#' \code{[Diag(X_i) * Phi_[1:rank[i]]]_i} for the i-th predictor in \code{X}.
#'
#' @param X base design matrix, including intercept
#' @param Phi Laplacian eigenvectors
#' @param rank expansion ranks, a vector with the same length as \code{ncol(X)}; defaults to a full expansion with \code{rep(ncol(Phi), ncol(X))}
#' @return the expanded design in a matrix with \code{rank} as attribute
#' @export
expand_design <- function (X, Phi, rank = rep(ncol(Phi), ncol(X))) {
  xnames <- colnames(X)
  if (is.null(xnames)) xnames <- paste0("X", 1:ncol(X))
  D <- matrix(nrow = nrow(X), ncol = sum(rank))
  xpos <- c(0, cumsum(rank[-length(rank)])) # column offsets
  dnames <- character(sum(rank))
  for (k in seq_along(rank)) {
    Phi_k <- Phi[, 1:rank[k], drop = FALSE]
    D[, xpos[k] + 1:rank[k]] <- sweep(Phi_k, 1, X[, k], `*`)
    dnames[xpos[k] + 1:rank[k]] <- paste0(xnames[k], "_", 1:rank[k])
  }
  attr(D, "rank") <- rank
  colnames(D) <- dnames
  D
}


#' Forward Stagewise Additive Modeling (FSAM).
#'
#' Performs the FSAM regression \code{y ~ X} under family \code{fam}, with
#' coefficients having prior distribution \code{N(0, (lambda t(X) * L * X)^-)},
#' with \code{L} taken as a differential operator such as a Laplacian.
#'
#' @param y response vector
#' @param X expanded design matrix, including intercept (see \code{expand_design})
#' @param L Laplacian matrix
#' @param lambda prior precision scale (or regularization penalty)
#' @param fam exponential family distribution, such as \code{family} in \code{glm}; defaults to \code{poisson}
#' @param maxit maximum number of iterations, defaults to maximum rank in \code{X}
#' @param trace print iteration progress messages?
#' @return list containing coefficients \code{theta} and \code{deviance}
#' @export
fsam_network <- function (y, X, L, lambda, fam = poisson(),
                          maxit = max(attr(X, "rank")), trace = FALSE) {
  rank <- attr(X, "rank")
  p <- length(rank)
  xpos <- c(0, cumsum(rank[-p])) # column offsets
  theta <- rep(NA, ncol(X))
  dev <- numeric(maxit)
  offset <- rep(0, nrow(X))

  P <- vector("list", length = p) # precision matrices
  for (j in 1:p) {
    Xj <- X[, xpos[j] + 1:rank[j]]
    P[[j]] <- crossprod(Xj, as.matrix(L %*% Xj))
  }

  # for each order k, regress theta_[i, k] given theta_[i, 1:(k - 1)]:
  for (k in 1:maxit) {
    valid <- which(k <= rank)
    pk <- length(valid)
    mu <- numeric(pk); pr <- numeric(pk) # conditional mean and precision
    for (j in 1:pk) {
      Pj <- P[[valid[j]]]; pr[j] <- Pj[k, k]
      if (k > 1) {
        theta_k <- theta[xpos[valid[j]] + 1:(k - 1)]
        mu[j] <- -sum(Pj[k, 1:(k - 1)] / Pj[k, k] * theta_k)
      }
    }
    pt <- list(mean = mu, precision = lambda * pr)

    Xk <- X[, k + xpos[valid], drop = FALSE]
    tk <- glm(y ~ Xk - 1, family = fam, offset = offset,
              method = bsglm::fitter(prior_coef = pt))
    theta[k + xpos[valid]] <- tk$coef
    dev[k] <- deviance(tk)
    offset <- offset + Xk %*% tk$coef
    if (trace) message("[", k, "] FSAM deviance = ", dev[k])
  }
  list(theta = theta, deviance = dev)
}


#' Expectation-Maximization (EM) for network basis rank determination.
#'
#' Given the rank-selection prior \code{theta | tau ~ N(0, (lambda Omega)^-)},
#' with \code{Omega = direct-sum_i M_[tau[i]] * t(X_i) * L * X_i M_[tau[i]]},
#' regresses \code{y ~ X}, \code{X} an expanded design, under family \code{fam}
#' via conditional EM with ranks \code{tau} as a latent. \code{tau} has log
#' prior distribution specified in \code{log_tau_p}, defined with respect to
#' the maximum rank of \code{X}. FSAM is used to initialize the procedure.
#'
#' @param y response vector
#' @param X expanded design matrix, including intercept (see \code{expand_design})
#' @param L Laplacian matrix
#' @param lambda prior precision scale (or regularization penalty)
#' @param log_tau_p log prior for \code{tau} (see \code{tau_prior})
#' @param v0 prior variance for non-selected basis, defaults to zero
#' @param fam exponential family distribution, such as \code{family} in \code{glm}; defaults to \code{poisson}
#' @param offset linear predictor offset, defaults to zero
#' @param start initial coefficient estimates
#' @param control list specifying running parameters (see \code{glm.control})
#' @return list containing coefficients \code{coef}, estimated log posterior \code{log P(tau | theta, y)} in \code{logptau} as (log) E-weights, and \code{deviance}
#' @export
em_network_rank <- function (y, X, L, lambda, log_tau_p, v0 = 0,
                             fam = poisson(), offset = 0,
                             start = NULL, control = list()) {
  control <- do.call("glm.control", control)
  link_is_canonical <- check_canonical(fam)

  rank <- attr(X, "rank")
  p <- length(rank)
  xpos <- c(0, cumsum(rank[-p])) # column offsets
  if (!is.null(start)) {
    theta <- start
    offset <- X %*% theta + offset
    dev <- sum(fam$dev.resids(y, fam$linkinv(offset), 1))
  } else {
    fn <- fsam_network(y, X, L, lambda, fam) # FIXME: use offset
    theta <- fn$theta; dev <- tail(fn$dev, 1)
    offset <- X %*% theta + offset
  }
  dev <- dev + lambda * mat_quad(theta, crossprod(X, as.matrix(L %*% X)))
  logptau <- numeric(sum(rank)) # log E-weights
  attr(logptau, "rank") <- rank
  if (control$trace) message("[0] deviance = ", dev)

  # cache log determinants
  log_det_prec <- vector(mode = "list", length = p)
  for (j in 1:p) {
    Xj <- X[, xpos[j] + 1:rank[j]]
    cj <- if (j == 1) rep(0, rank[j]) else # intercept?
      log(diag(chol(lambda * crossprod(Xj, as.matrix(L %*% Xj)))))
    log_det_prec[[j]] <- c(0, 2 * cumsum(cj))
  }

  for (it in 1:control$maxit) {
    for (j in 1:p) { # cycle
      rangej <- xpos[j] + 1:rank[j]
      Xj <- X[, rangej, drop = FALSE]; thetaj <- theta[rangej]
      offset <- offset - mat_mult(Xj, thetaj)

      # [ E-step ]
      tj <- numeric(rank[j] + 1)
      ldpj <- log_det_prec[[j]]
      for (tau in 0:rank[j]) {
        Mtau <- ifelse(1:rank[j] > tau, v0, 1)
        eta <- mat_mult(Xj, thetaj * Mtau) # select
        devj <- sum(fam$dev.resids(y, fam$linkinv(offset + eta), 1))
        u <- mat_mult(Xj, thetaj * sqrt(Mtau))
        tj[tau + 1] <- log_tau_p[tau + 1] +
          .5 * (ldpj[tau + 1] - lambda * mat_quad(u , L) - devj)
      }
      tj <- tj - log_sum_exp(tj) # normalize
      logptau[rangej] <- tj[-1]
      pij <- exp(tj)
      nu <- cumsum(pij[-(rank[j] + 1)])
      Tj <- outer(1:rank[j], 1:rank[j], Etau(v0, nu))
      Omega <- lambda * Tj * crossprod(Xj, as.matrix(L %*% Xj)) # prior precision

      # [ M-step ]
      pij <- pij[-1]
      ut <- -Omega %*% theta[rangej]; Ht <- Omega
      for (tau in 1:rank[j]) {
        tauj <- 1:tau
        Xjt <- Xj[, tauj, drop = FALSE]
        etaj <- drop(mat_mult(Xjt, thetaj[tauj]) + offset)
        wj <- fam$mu.eta(etaj)
        muj <- fam$linkinv(etaj); rj <- y - muj
        if (!link_is_canonical) { # adjust?
          aj <- wj / fam$variance(muj)
          wj <- aj * wj; rj <- aj * rj
        }
        ut[tauj] <- ut[tauj] + crossprod(Xjt, rj) * pij[tau]
        Ht[tauj, tauj] <- Ht[tauj, tauj] + crossprod(Xjt, Xjt * wj) * pij[tau]
      }
      #theta[rangej] <- drop(chol_solve(chol(Ht), Ht %*% thetaj + ut))
      eht <- eigen(Ht, symm = TRUE)
      eh_values <- eht$values
      ind <- which(eh_values > 1e-4) # FIXME
      eh_values[-ind] <- 0; eh_values[ind] <- 1 / eh_values[ind] # pseudo-inverse
      ht <- crossprod(eht$vectors, Ht %*% thetaj + ut) * eh_values
      theta[rangej] <- drop(eht$vectors %*% ht)

      offset <- offset + mat_mult(Xj, theta[rangej])
    }
    dev_new <- sum(fam$dev.resids(y, fam$linkinv(offset), 1)) +
      lambda * mat_quad(theta, crossprod(X, as.matrix(L %*% X)))
    if (control$trace) message("[", it, "] deviance = ", dev_new)
    if ((dev - dev_new) / (dev + .1) < control$epsilon) break
    dev <- dev_new
  }
  list(coef = theta, logptau = logptau, deviance = dev, family = fam)
}


#' Sequential centroid estimation for network basis rank.
#'
#' Computes the sequential centroid estimator for the basis rank \code{tau}
#' based on the probabilities in \code{ptau}. The threshold is based on the
#' true positive gain \code{kappa}, given by \code{kappa / (1 + kappa)}.
#'
#' @param ptau probabilities of tau = 1, ..., rank[i] for each predictor i; the vector needs an attribute \code{rank}, as in the expanded design
#' @param kappa true positive gain, defaults to unit
#' @return the sequential centroid estimator
#' @export
centroid_network_rank <- function (ptau, kappa = 1) {
  kt <- kappa / (1 + kappa)
  rank <- attr(ptau, "rank")
  xpos <- c(0, cumsum(rank)[-length(rank)])
  ce <- numeric(length(rank))
  for (j in seq_along(rank)) {
    ptj <- ptau[xpos[j] + 1:rank[j]]
    ptj <- c(1 - sum(ptj), ptj) # add P(tau = 0)
    ce[j] <- min(which(cumsum(ptj) - kt > 0)) - 1
  }
  ce
}


#' Sub-index rank selection.
#'
#' Gets the index range from \code{rank} based on the rank selection specified
#' in \code{tau}.
#'
#' @param rank full rank indices, as returned from an expanded design \code{X}
#' with \code{attr(X, "rank")} (see also \code{expand_design})
#' @param tau rank subsets, with same length as \code{rank}
#' @return index range
#' @export
get_index_rank <- function (rank, tau) {
  p <- length(rank) # == length(tau)
  rpos <- c(0, cumsum(rank[-p]))
  ipos <- c(0, cumsum(tau[-p]))
  ind_range <- integer(sum(tau))
  for (j in seq_len(p))
    ind_range[ipos[j] + 1:tau[j]] <- rpos[j] + 1:tau[j]
  ind_range
}

#' Expectation-Maximization (EM) for network basis rank determination and hot
#' zones.
#'
#' Given the rank-selection priors \code{theta | tau ~ N(0, (lambda Omega)^-)} and
#' \code{omega | tau_hz ~ N(0, (lambda_hz Omega_hz)^-)}, with
#' \code{Omega = direct-sum_i M_[tau[i]] * t(X_i) * L * X_i M_[tau[i]]} and
#' \code{Omega_hz = direct-sum_i M_[tau_hz[i]] * t(X_hz_i) * L * X_hz_i M_[tau_hz[i]]},
#' Poisson regresses \code{y | Z = 1 ~ X} and \code{y | Z = 0 ~ X_bg} with
#' \code{Z ~ X_hz}, the "hot zone" indicators, also network regularized via
#' \code{omega}. Here \code{X} and \code{X_hz} are network-expanded designs
#' (see \code{expand_design}). The conditional regressions are performed via
#' conditional EM with ranks \code{tau} and \code{tau_hz} and hot zone
#' indicators as latent. Both \code{tau} and \code{tau_hz} have log prior
#' distribution specified in \code{log_tau_p}, defined with respect to the
#' maximum rank of \code{X}.
#'
#' @param y response vector
#' @param X expanded design matrix for hot zone counts, including intercept (see \code{expand_design})
#' @param X_bg design matrix for background counts, including intercept
#' @param X_hz expanded design matrix for hot zone indicators, including intercept (see \code{expand_design}); defaults to \code{X}
#' @param L Laplacian matrix
#' @param lambda prior precision scale for hot zone counts (or regularization penalty)
#' @param lambda_hz prior precision scale for hot zone indicators (or regularization penalty)
#' @param log_tau_p log prior for \code{tau} (see \code{tau_prior})
#' @param p_bg_init initial background probabilities, assumed binary
#' @param control list specifying running parameters (see \code{glm.control})
#' @return list containing coefficients \code{coef}, \code{coef_bg}, and \code{coef_hz} for the count regressions on hot zones and background and logistic regression on hot zone indicators, respectively; estimated log posteriors \code{log P(tau | theta, y)} and \code{log P(tau_hz | omega, y)} in \code{logptau} and \code{logptau_hz} as (log) E-weights; and \code{deviance}
#' @export
em_network_rank_bg <- function (y, X, X_bg, X_hz = X, L, lambda, lambda_hz,
                                log_tau_p, p_bg_init, control = list()) {
  control <- do.call("glm.control", control)
  one_iteration <- list(maxit = 1)
  prec <- lambda * crossprod(X, as.matrix(L %*% X))
  prec_hz <- lambda_hz * crossprod(X_hz, as.matrix(L %*% X_hz))
  p_bg <- p_bg_init

  o <- p_bg; o[p_bg == 1] <- log(1 - control$epsilon)
  ft <- em_network_rank((1 - p_bg) * y, X, L, lambda, log_tau_p,
                        fam = quasipoisson(), offset = o,
                        control = one_iteration)
  tau <- centroid_network_rank(exp(ft$logptau))
  fz <- glm(y ~ X_bg - 1, family = quasipoisson, subset = p_bg == 1,
            control = one_iteration)
  fw <- em_network_rank(p_bg, X_hz, L, lambda_hz, log_tau_p, fam = quasibinomial(),
                        control = one_iteration)
  tau_hz <- centroid_network_rank(exp(fw$logptau))

  dev <- ft$deviance + fz$deviance + fw$deviance +
    mat_quad(ft$coef, prec) + mat_quad(fw$coef, prec_hz)
  if (control$trace) message("[0] deviance = ", round(dev, 2))

  for (it in 1:control$maxit) {
    # [ E-step ]
    theta <- ft$coef; zeta <- fz$coef; omega <- fw$coef
    logptau <- ft$logptau; logptau_hz <- fw$logptau
    ind <- get_index_rank(attr(X, "rank"), tau)
    eta <- drop(X[, ind, drop = FALSE] %*% theta[ind])
    ind <- get_index_rank(attr(X_hz, "rank"), tau_hz)
    eta_hz <- drop(X_hz[, ind, drop = FALSE] %*% omega[ind])
    eta_bg <- drop(X_bg %*% zeta)
    p_bg <- plogis(eta_hz + y * eta_bg - exp(eta_bg) - (y * eta - exp(eta)))
    # p_bg in (eps, 1 - eps), avoid infs:
    p_bg <- pmin(pmax(p_bg, control$epsilon), 1 - control$epsilon)

    # [ M-step theta ]
    ft <- em_network_rank((1 - p_bg) * y, X, L, lambda, log_tau_p,
                          fam = quasipoisson(), offset = log(1 - p_bg),
                          start = theta,
                          control = one_iteration)
    tau <- centroid_network_rank(exp(ft$logptau))
    # [ M-step zeta ]
    fz <- suppressWarnings(glm(p_bg * y ~ X_bg - 1, offset = log(p_bg),
                               family = quasipoisson, start = zeta,
                               control = one_iteration))
    # [ M-step omega ]
    fw <- em_network_rank(p_bg, X_hz, L, lambda_hz, log_tau_p,
                          fam = quasibinomial(), start = omega,
                          control = one_iteration)
    tau_hz <- centroid_network_rank(exp(fw$logptau))

    dev_new <- ft$deviance + fz$deviance + fw$deviance +
      mat_quad(ft$coef, prec) + mat_quad(fw$coef, prec_hz)
    if (control$trace) message("[", it, "] deviance = ", round(dev_new, 2))
    if ((dev - dev_new) / dev < control$epsilon) break
    dev <- dev_new
  }
  list(lambda = lambda, lambda_hz = lambda_hz,
       coef = theta, coef_bg = zeta, coef_hz = omega,
       logptau = logptau, logptau_hz = logptau_hz,
       p_bg = p_bg, deviance = dev)
}


#' Expectation-Maximization (EM) for hot zones.
#'
#' Given the priors \code{theta ~ N(0, (lambda Omega)^-)} and
#' \code{omega ~ N(0, (lambda_hz Omega_hz)^-)}, with
#' \code{Omega = direct-sum_i t(X_i) * L * X_i} and
#' \code{Omega_hz = direct-sum_i t(X_hz_i) * L * X_hz_i},
#' Poisson regresses \code{y | Z = 1 ~ X} and \code{y | Z = 0 ~ X_bg} with
#' \code{Z ~ X_hz}, the "hot zone" indicators, also network regularized via
#' \code{omega}. Here \code{X} and \code{X_hz} are network-expanded designs
#' (see \code{expand_design}). The conditional regressions are performed via
#' conditional EM with hot zone indicators as latent.
#'
#' @param y response vector
#' @param X expanded design matrix for hot zone counts, including intercept (see \code{expand_design})
#' @param X_bg design matrix for background counts, including intercept
#' @param X_hz expanded design matrix for hot zone indicators, including intercept (see \code{expand_design}); defaults to \code{X}
#' @param L Laplacian matrix
#' @param lambda prior precision scale for hot zone counts (or regularization penalty)
#' @param lambda_hz prior precision scale for hot zone indicators (or regularization penalty)
#' @param p_bg_init initial background probabilities, assumed binary
#' @param control list specifying running parameters (see \code{glm.control})
#' @return list containing coefficients \code{coef}, \code{coef_bg}, and \code{coef_hz} for the count regressions on hot zones and background and logistic regression on hot zone indicators, respectively; estimated \code{P(Z_i = 1 | y)} as E-weights \code{p_bg}; and \code{deviance} and press-loop cross-validation loss \code{pl}.
#' @export
em_network_bg <- function (y, X, X_bg, X_hz = X, L, lambda, lambda_hz,
                           p_bg_init, control = list()) {
  control <- do.call("glm.control", control)
  one_iteration <- list(maxit = 1)
  prec <- lambda * crossprod(X, as.matrix(L %*% X))
  prec_hz <- lambda_hz * crossprod(X_hz, as.matrix(L %*% X_hz))
  net_fitter <- make_fitter(lambda, X, L)
  net_fitter_hz <- make_fitter(lambda_hz, X_hz, L)
  p_bg <- p_bg_init

  ft <- glm(y ~ X - 1, family = quasipoisson, subset = p_bg == 0,
            method = net_fitter, control = one_iteration)
  fz <- suppressWarnings(glm(y ~ X_bg - 1, family = quasipoisson,
                             subset = p_bg == 1, control = one_iteration))
  fw <- glm(p_bg ~ X_hz - 1, family = quasibinomial,
            method = net_fitter_hz)

  dev <- ft$deviance + fz$deviance + fw$deviance +
    mat_quad(ft$coef, prec) + mat_quad(fw$coef, prec_hz)
  if (control$trace) message("[0] deviance = ", round(dev, 2))

  for (it in 1:control$maxit) {
    # [ E-step ]
    theta <- ft$coef; zeta <- fz$coef; omega <- fw$coef
    eta <- drop(X %*% theta); eta_hz <- drop(X_hz %*% omega)
    eta_bg <- drop(X_bg %*% zeta)
    p_bg <- plogis(eta_hz + y * eta_bg - exp(eta_bg) - (y * eta - exp(eta)))
    # p_bg in (eps, 1 - eps), avoid infs:
    p_bg <- pmin(pmax(p_bg, control$epsilon), 1 - control$epsilon)
    # press loop, for convenience:
    if (it > 1) {
      pl <- 0
      if (class(ft)[1] != "try-error")
        pl <- pl + sum((1 - p_bg) * press_loop(ft))
      if (class(fz)[1] != "try-error")
        pl <- pl + sum(p_bg * press_loop(fz))
    }

    # [ M-step theta ]
    ft <- glm((1 - p_bg) * y ~ X - 1, family = quasipoisson,
              offset = log(1 - p_bg), start = theta,
              method = net_fitter, control = one_iteration)
    # [ M-step zeta ]
    fz <- suppressWarnings(glm(p_bg * y ~ X_bg - 1, offset = log(p_bg),
                               family = quasipoisson, start = zeta,
                               control = one_iteration))
    # [ M-step omega ]
    fw <- glm(p_bg ~ X_hz - 1, family = quasibinomial, start = omega,
              method = net_fitter_hz, control = one_iteration)

    dev_new <- ft$deviance + fz$deviance + fw$deviance +
      mat_quad(ft$coef, prec) + mat_quad(fw$coef, prec_hz)
    if (control$trace) message("[", it, "] deviance = ", round(dev_new, 2))
    if ((dev - dev_new) / dev < control$epsilon) break
    dev <- dev_new
  }
  list(lambda = lambda, lambda_hz = lambda_hz,
       coef = theta, coef_bg = zeta, coef_hz = omega,
       p_bg = p_bg, dev = dev, pl = pl)
}


#' Plot EM hot-zone network fit.
#'
#' Given the output of a EM hot-zone network fit, plots observed vs espected counts, estimated probabilities of background vs observed counts, and residuals vs expected counts.
#'
#' @param y response vector
#' @param X expanded design matrix for hot zone counts, including intercept (see \code{expand_design})
#' @param X_bg design matrix for background counts, including intercept
#' @param res list as returned from \code{em_network_bg}
#' @return list with \code{fitted} values and \code{residuals}, as side effect (invisible)
#' @export
plot_em_fit <- function (y, X, X_bg, res) {
  eta_bg <- X_bg %*% res$coef_bg
  eta <- X %*% res$coef
  op <- par(mfrow = c(1, 3))
  muz <- exp(eta_bg); mue <- exp(eta) # Poisson family
  mu <- res$p_bg * muz + (1 - res$p_bg) * mue
  plot(y, mu); abline(0, 1)
  plot(y, res$p_bg)
  v <- mu + res$p_bg * (1 - res$p_bg) * (muz - mue) ^ 2
  r <- (y - mu) / sqrt(v)
  plot(mu, r); lines(lowess(mu, r), col = "red", lwd = 2)
  par(op)
  invisible(list(fitted = mu, residuals = r))
}

