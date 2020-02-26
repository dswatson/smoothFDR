#' Smooth FDR workhorse
#' 
#' This is the workhorse of the smoothFDR package. The eponymous function
#' is just a wrapper.
#'
#' @import FDRreg
#' @import glmgen
#' @import dplyr
#' 

smooth_fdr <- function(probe, z, pos, chr, nulltype, nlambda, tol, maxit) {
  ### Part I: Estimate null and alternative distros
  n <- length(z)
  x <- pos - min(pos) + 1
  fdr1 <- FDRreg2(z, x, nulltype = nulltype)
  f0 <- fdr1$M0
  f1 <- fdr1$M1
  p0 <- fdr1$p0
  ### Part II: Expectation-Maximization
  # Initial parameters for EM algorithm
  drift <- 1
  iter <- 1
  prior_prob <- rep(1 - p0, n)
  beta_hat <- qlogis(prior_prob)
  m1 <- prior_prob * f1
  m0 <- (1 - prior_prob) * f0
  post_prob <- m1 / (m1 + m0)
  wts <- prior_prob * (1 - prior_prob)
  y <- beta_hat - (prior_prob - post_prob) / wts
  fl0 <- trendfilter(y, weights = wts, k = 0, nlambda = nlambda)
  obj_old <- sum(log(prior_prob * f1 + (1 - prior_prob) * f0))
  # Data frame to keep track of BIC, matrix for pi vectors
  res <- data.frame('lambda' = fl0$lambda, 'BIC' = NA_real_)
  pi <- matrix(nrow = n, ncol = nlambda + 1)
  pi[, 1] <- prior_prob
  # EM algorithm
  em <- function(lambda, prior_prob) {
    idx <- which(res$lambda == lambda)
    beta_hat <- qlogis(prior_prob)
    while(drift >= tol && iter <= maxit) {
      if (iter == 1) {
        fl <- trendfilter(y, weights = wts, k = 0, lambda = lambda)
      } else {
        # E step
        m1 <- prior_prob * f1
        m0 <- (1 - prior_prob) * f0
        post_prob <- m1 / (m1 + m0)
        # M step
        wts <- prior_prob * (1 - prior_prob)
        y <- beta_hat - (prior_prob - post_prob) / wts
        fl <- trendfilter(y, weights = wts, k = 0, lambda = lambda)
      }
      # Update  
      beta_hat <- as.numeric(fl$beta)
      prior_prob <- plogis(beta_hat)
      obj_new <- sum(log(prior_prob * f1 + (1 - prior_prob) * f0))
      drift <- abs(obj_old - obj_new) / (abs(obj_old) + tol)
      iter <- iter + 1
      obj_old <- obj_new
    }
    out <- list('BIC' = -2 * obj_new + log(n) * fl$df, 'pi1' = prior_prob)
    return(out)
  }
  # Loop over lambdas using previous outputs as warm starts
  for (i in seq_len(nlambda)) {
    tmp <- em(res$lambda[i], pi[, i])
    res$BIC[i] <- tmp$BIC
    pi[, (i + 1)] <- tmp$pi1
  }
  pi <- pi[, -1]
  ### Part III: Use posteriors to compute FDR
  j <- which.min(res$BIC)
  m1 <- pi[, j] * f1
  m0 <- (1 - pi[, j]) * f0
  post_prob <- m1 / (m1 + m0)
  # Export data frame
  data.frame(Idx = seq_len(n), probe = probe, z = z, p.value = pnorm(z), 
             pos = pos, chr = chr) %>%
    mutate(lfdr = 1 - post_prob) %>%
    arrange(lfdr) %>%
    mutate(q.value = cummean(lfdr)) %>%
    arrange(Idx) %>%
    select(-Idx) %>%
    return(.)
}






