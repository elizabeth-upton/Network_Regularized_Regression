% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/banner.R
\name{tau_prior}
\alias{tau_prior}
\title{Prior for basis indicator \code{tau}, using hyper-parameter \code{tau_bar}.}
\usage{
tau_prior(K, tau_bar, alpha_0, alpha_1, ...)
}
\arguments{
\item{K}{maximum basis expansion}

\item{tau_bar}{prior expected rank \code{E(tau)}}

\item{alpha_0}{probability of rank > 0}
}
\value{
alpha_1 probability of rank > 1 given rank > 0

vector with probabilities \code{P(tau | K)}
}
\description{
Computes the prior rank probabilities \code{P(tau | K)} where \code{K} is
the maximum basis expansion based on parameters \code{alpha_0 = P(tau > 0)},
\code{alpha_1 = P(tau > 1 | tau > 0)}, and \code{rho}, where
\code{P(tau | tau > 1)} is proportional to \code{rho^{tau - 1}}.
Parameter \code{rho} is computed from \code{tau_bar}.
}
