#' Smooth FDR workhorse
#' 
#' This is the workhorse of the smoothFDR package. The eponymous function
#' is just a wrapper.
#'
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

#' The following functions (and RcppExports.R) are copy/pasted from 
#' jgscott's FDRreg package:
#' https://github.com/jgscott/FDRreg/tree/master/R_pkg/R
#' I've embedded the code and changed function names to distinguish between
#' this and the CRAN package of the same name.
#' 

FDRreg2 = function(z, features, nulltype='theoretical', method='pr', stderr = NULL, control=list()) {
  # False discovery rate regression
  # z = vector of z scores
  # features = design matrix of covariates, assumed NOT to have an intercept just as in vanilla lm()
  # nulltype = flag for what kind of null hypothesis to assume, theoretical/empirical/heteroscedastic
  
  stopifnot(any(method=='pr', method=='efron'))
  
  # Set up control parameters
  mycontrol = list(center=TRUE, scale=TRUE)
  if(method=='pr') {
    mycontrol$gridsize = 300
    mycontrol$decay = -0.67
    mycontrol$npasses = 10
    mycontrol$lambda = 0.01
    mycontrol$densknots=10
    mycontrol$nmids=150
  } else if(method=='efron') {
    mycontrol$gridsize = 150
    mycontrol$nmids=150
    mycontrol$densknots=10
  }
  # Overwrite with user choices
  mycontrol[(namc <- names(control))] <- control
  
  # Matrix of regressors, centered and scaled as requested
  N = length(z)
  X = scale(features, center=mycontrol$center, scale=mycontrol$scale)
  P = ncol(X)
  
  # Estimate the marginal density
  if(method=='pr') {
    
    # Compute M0 and M1, the marginals under null and alternative for each observation
    if(nulltype=='empirical') {
      
      # Currently using my implementation of Efron's central matching estimator
      l1 = efron(z, nmids=mycontrol$nmids, df=mycontrol$densknots, nulltype=nulltype)
      mu0 = l1$mu0
      sig0 = l1$sig0
      
      prfit = prfdr(z, mu0, sig0, control=mycontrol)
      fmix_grid = prfit$fmix_grid
      f0_grid = dnorm(prfit$x_grid, mu0, sig0)
      f1_grid = prfit$f1_grid
    } else if(nulltype=='heteroscedastic') {
      if(missing(stderr)) {
        stop("Must specify standard error (stderr) if assuming heteroscedastic null.")
      }
      mu0 = 0.0
      sig0 = stderr
      prfit = prfdr_het(z, mu0, sig0, control=mycontrol)
      fmix_grid = NULL
      f0_grid = NULL
      f1_grid = NULL
    } else {
      mu0 = 0.0
      sig0 = 1.0
      prfit = prfdr(z, mu0, sig0, control=mycontrol)
      fmix_grid = prfit$fmix_grid
      f0_grid = dnorm(prfit$x_grid, mu0, sig0)
      f1_grid = prfit$f1_grid
    }
    
    # Extract marginal densities and fit regression
    p0 = prfit$pi0
    M0 = dnorm(z, mu0, sig0)
    M1 = prfit$f1_z
    m1zeros = which(M1 < .Machine$double.eps)
    if(length(m1zeros > 0)) {
      M1[m1zeros] = min(M1[-m1zeros]) # substitute in the smallest nonzero value
      M1 = pmax(M1, .Machine$double.eps) # shouldn't happen but just in case!
    }
    m0zeros = which(M0 < .Machine$double.eps)
    if(length(m0zeros > 0)) {
      M0[m0zeros] = min(M0[-m0zeros]) # substitute in the smallest nonzero value
      M0 = pmax(M0, .Machine$double.eps) # shouldn't happen but just in case!
    }
    x_grid = prfit$x_grid
    regressfit = fdrr_regress_pr(M0, M1, X, 1-p0, lambda=mycontrol$lambda)
    
  } else if(method=='efron') {
    if(nulltype=='heteroscedastic') {
      stop("Cannot use Efron's method under a heteroscedastic null.")
    }
    l1 = efron(z, nmids=mycontrol$nmids, df=mycontrol$densknots, nulltype=nulltype)
    mu0 = l1$mu0
    sig0 = l1$sig0
    p0 = l1$p0
    M0 = dnorm(z, mu0, sig0)
    M1 = NULL
    MTot = l1$fz	
    x_grid = l1$mids
    fmix_grid = l1$zdens
    f0_grid = dnorm(x_grid, mu0, sig0)
    f1_grid = NULL
    regressfit = fdrr_regress_efron(M0, MTot, X, 1-p0, N)
  }
  
  out2 = getFDR(regressfit$PostProb)
  list(	z=z, X=X, localfdr=out2$localfdr, FDR=out2$FDR, x_grid = x_grid,
        M0 = M0, M1 = M1,
        fmix_grid=fmix_grid, f0_grid = f0_grid, f1_grid = f1_grid, 
        mu0=mu0, sig0=sig0, p0=p0, priorprob = regressfit$W,
        postprob = regressfit$PostProb, model=regressfit$model
  )
  
}



BayesFDRreg = function(z, features, mu0=NULL, sig0 = NULL, empiricalnull=FALSE, nmc=5000, nburn=1000,
                       control=list(), ncomps=NULL, priorpars = NULL) {
  # Fully Bayesian version of false discovery rate regression
  # z = vector of z scores
  # features = design matrix of covariates, assumed NOT to have an intercept just as in vanilla lm()
  # nulltype = flag for what kind of null hypothesis to assume, theoretical or empirical
  # ncomps = how many mixture components for the alternative hypothesis
  
  mycontrol = list(center=FALSE, scale=FALSE, verbose=nmc+nburn+1)
  mycontrol[(namc <- names(control))] <- control
  
  N = length(z)
  X = cbind(1,scale(features, center= mycontrol$center, scale= mycontrol$scale))
  P = ncol(X)
  
  if(empiricalnull) {
    l1 = efron(z, nmids=150, df=15, nulltype='empirical')
    mu0 = l1$mu0
    sig0 = l1$sig0
    p0 = l1$p0
  } else {
    if(missing(sig0)) sig0 = rep(1,N)
    if(missing(mu0)) mu0 = 0
    p0 = NULL
  }
  sig0squared = sig0^2
  M0 = dnorm(z, mu0, sig0)
  
  # Initialize MCMC
  if(missing(priorpars)) {
    PriorPrec = diag(rep(1/25, P))
    PriorMean = rep(0,P)
  } else{
    PriorPrec = priorpars$PriorPrec
    PriorMean = priorpars$PriorMean
  }
  
  if(missing(ncomps)) {
    foundfit = FALSE
    ncomps = 1
    emfit = deconvolveEM(z, ncomps)
    while(!foundfit) {
      newfit = deconvolveEM(z, ncomps+1)
      if(newfit$AIC > emfit$AIC) {
        foundfit = TRUE
      } else {
        emfit = newfit
        ncomps = ncomps+1
      }
    }
    M1 = dnormix(z, emfit$weights[1:ncomps]/sum(emfit$weights[1:ncomps]),
                 emfit$means[1:ncomps], emfit$vars[1:ncomps])
  } else 	M1 = dnorm(z, 0, 4)
  
  PriorPrecXMean = PriorPrec %*% PriorMean
  Beta = rep(0,P)
  Beta[1] = -3
  BetaSave = matrix(0, nrow=nmc, ncol=P)
  MuSave = matrix(0, nrow=nmc, ncol=ncomps)
  VarSave = matrix(0, nrow=nmc, ncol=ncomps)
  WeightsSave = matrix(0, nrow=nmc, ncol=ncomps)
  PostProbSave = 0
  PriorProbSave = 0
  M1Save = 0
  
  # Alternative hypothesis
  comp_weights = rep(1/ncomps, ncomps)
  comp_means = quantile(z[abs(z/sig0)>2], probs=seq(0.025,0.975,length=ncomps))
  comp_variance = rep(1, ncomps)
  myvar = comp_variance + 1
  
  # Main MCMC
  for(t in 1:(nmc+nburn)) {
    
    if(t %% mycontrol$verbose == 0) cat(t, "\n")
    
    ### Update indicators
    Psi = drop(X %*% Beta)
    W = ilogit(Psi)
    PostProb = W*M1/{(1-W)*M0 + W*M1}		
    Gamma = rbinom(N,1,PostProb)
    
    
    ### Update mixture of normals model
    cases = which(Gamma==1)
    signals = z[cases]
    components = draw_mixture_component(signals, sig0[cases], weights=comp_weights, mu = comp_means, tau2 = comp_variance) + 1
    
    # Draw latent means
    if(length(cases) > 0) {
      latentmeans.var = 1.0/(1.0/sig0[cases]^2 + 1.0/comp_variance[components])
      latentmeans.mu = latentmeans.var*(signals/(sig0[cases]^2) + (comp_means/comp_variance)[components])
      latentmeans = rnorm(length(cases), latentmeans.mu, sqrt(latentmeans.var))
      nsig = mosaic::maggregate(signals ~ factor(components, levels=1:ncomps), FUN='length')
      tss_thetai = mosaic::maggregate((latentmeans-comp_means[components])^2 ~ factor(components, levels=1:ncomps), FUN='sum')
      sum_thetai = mosaic::maggregate(latentmeans ~ factor(components, levels=1:ncomps), FUN='sum')
    } else {
      nsig = rep(0,ncomps)
      tss_thetai = rep(0,ncomps)
      mean_thetai = rep(0,ncomps)
    }
    
    # Actual updates
    for(k in 1:ncomps) comp_variance[k] = 1.0/rgamma(1, {nsig[k]+2}/2, rate={tss_thetai[k]+2}/2)
    muvar = comp_variance/{nsig + comp_variance*0.1}
    muhat = sum_thetai/{nsig + comp_variance*0.1}
    comp_means = rnorm(ncomps, muhat, sqrt(muvar))	
    comp_weights = rdirichlet_once(rep(5, ncomps) + nsig)
    M1 = marnormix(z, sig0squared, comp_weights, comp_means, comp_variance)
    
    ### Update latent variables in logit likelihood
    Om = as.numeric(BayesLogit::rpg(N,rep(1,N),Psi))
    
    ### Update regression parameters
    Kap = PostProb - 0.5
    PrecMat = t(X) %*% {Om * X} + PriorPrec
    Beta.V = solve(PrecMat)
    Beta.mu = Beta.V %*% {t(X) %*% Kap + PriorPrecXMean}
    Beta = t(mvtnorm::rmvnorm(1,mean=Beta.mu,sigma=Beta.V))	
    if(t > nburn) {
      BetaSave[t-nburn,] = Beta
      MuSave[t-nburn,] = comp_means
      VarSave[t-nburn,] = comp_variance
      WeightsSave[t-nburn,] = comp_weights
      PostProbSave = PostProbSave + (1.0/nmc)*PostProb
      PriorProbSave = PriorProbSave + (1.0/nmc)*W
      M1Save = M1Save + (1.0/nmc)*M1
    }
  }
  out2 = getFDR(PostProbSave)
  
  mylist = list(z=z, localfdr=out2$localfdr, FDR=out2$FDR, X=X,
                M0 = M0, M1 = M1Save, mu0=mu0, sig0=sig0, p0=p0, ncomps=ncomps,
                priorprob = PriorProbSave, postprob = PostProbSave, 
                coefficients = BetaSave, weights = WeightsSave, means=MuSave, vars = VarSave
  )
  return(mylist);
}

ilogit = function(x) 1/{1+exp(-x)}
flogit = function(x) log(x/{1-x})

# Truncated gamma random draws
rtgamma = function(n, a, b, lb, ub) {
  lub = pgamma(lb, a, rate=b)
  uub = pgamma(ub, a, rate=b)
  u = runif(n, lub, uub)
  qgamma(u, a, rate=b)
}


SoftLogitFit = function(y, X, lambda=1e-6, start=NULL) {
  # y: a vector of "fuzzy" ones and zeros, i.e. expected binary outcomes
  # X: design matrix assumed to include a column of 1's for an intercept
  # lambda is a ridge penalty parameter, defaulting to 1e-6
  #	(basically no regularization, only there to ensure numerical stability)
  if(missing(start)) {
    start = rep(0, ncol(X))
  }
  mymax = optim(start, fn = SoftLogitLoss, gr = SoftLogitGradient,
                method='BFGS', y=y, X = X, lambda=lambda, hessian=TRUE)
  list(coef = mymax$par, value = mymax$value, hessian = mymax$hessian)
}


efron = function(z, nmids=150, pct=-0.01, pct0=0.25, df=10, nulltype='theoretical') {
  # estimate f(z) and f_0(z) using Efron (2004)'s method
  stopifnot(any(nulltype == 'theoretical', nulltype=='empirical'))
  
  result = tryCatch({
    
    N = length(z)
    med = median(z)
    myrange = med + (1 - pct) * (range(z) - med)
    lb = myrange[1]
    ub = myrange[2]
    
    breaks = seq(lb, ub, length= nmids +1)
    h1 = hist(z, breaks = breaks, plot = FALSE)
    mids = (breaks[-1] + breaks[-(nmids+1)])/2
    zcounts = h1$counts
    glm1 = glm(zcounts ~ splines::ns(mids, df = df), family=poisson)
    zrate = glm1$fit
    D = (zcounts - zrate)/sqrt(zrate+1)
    D = D[-c(1,nmids)]
    if (sum(D^2) > qchisq(0.9, nmids-2-df)) {
      warning(paste0("f(z) misfit = ", round(D, 1), ".  Rerun with increased df."))
    }	
    
    zdens = {zrate/sum(zrate)}/diff(breaks)
    
    # Now do spline interpolation for the density at the observed points
    ispl2 = splines::interpSpline( zdens ~ mids )
    fz = predict(ispl2, z)$y
    
    # Pick out the middle of the data points
    ql = quantile(z, pct0)
    qu = quantile(z, 1-pct0)
    ind0 = intersect(which(z > ql), which(z<qu))
    if(nulltype=='empirical') {
      # empirical null by central moment matching
      z0 = z[ind0]
      l0 = log(fz[ind0])
      zmax = z[which.max(l0)]
      lm0 = lm(l0~I(z0-zmax) + I((z0-zmax)^2))
      b0 = coef(lm0)
      sig = as.numeric(sqrt(-1.0/{2.0*b0[3]}))
      mu = as.numeric(-b0[2]/(2.0*b0[3]) + zmax)
      # lm0 = lm(l0 ~ z0 + I(z0^2))
      # b0 = coef(lm0)
      # sig = as.numeric(sqrt(-1/{2*b0[3]}))
      # mu = as.numeric(b0[2] * sig^2)
    } else {
      # theoretical null
      sig = 1
      mu = 0
    }
    p0 = sum(fz[ind0])/sum(dnorm(z[ind0], mu, sig))
    localfdr = pmin(1, p0*dnorm(z, mu, sig)/fz)
    list(mids=mids, breaks=breaks, zcounts=zcounts, zdens=zdens,
         z=z, fz=fz, mu0=mu, sig0=sig, p0=p0, fdr=localfdr)
  }, error = function(err) {
    print(err)
    list(mids=NULL, breaks=NULL, zcounts=NULL, zdens=NULL,
         z=NULL, fz=NULL, mu0=NULL, sig0=NULL, p0=NULL, fdr=NULL)
  }, finally = {
    # Nothing to clean up
  })
  return(result)
}

getFDR = function(postprob) {
  # postprob is a vector of posterior probabilities
  # from which local fdr and (Bayesian) FDR are extracted
  indices = 1:length(postprob)
  iorder = order(postprob, decreasing=TRUE)
  porder = postprob[iorder]
  localfdr.order = 1-porder
  FDR.order = cumsum(localfdr.order)/indices
  localfdr = indices  # placeholder
  localfdr[iorder] = localfdr.order
  FDR = indices  # placeholder
  FDR[iorder] = FDR.order
  
  # Where local fdr is 1, report the most conservative FDR 
  fdrmax = which(localfdr == 1)
  FDR[fdrmax] = max(FDR)
  list(localfdr=localfdr, FDR=FDR)
}

plotFDR = function(fdrr, Q=0.1, showrug=TRUE, showsub=TRUE, breaks=150, ...) {
  N = length(fdrr$z)
  mytitle = paste0('')
  par(mar=c(5,4,1,1))
  par(...)
  h1 = hist(fdrr$z, breaks, plot=FALSE)
  plot(h1, freq=FALSE, col='lightgrey', border='grey',
       main=mytitle, xlab='', ylab='', axes=FALSE, ylim=c(0, 1.05*max(h1$density)))
  mysub = paste0('Grey bars: original z scores\nBlue bars: fraction signals in each bin')
  axis(1, pos=0, tick=FALSE, cex.axis=0.9)
  axis(2, tick=FALSE, las=1, cex.axis=0.9)
  zcut = data.frame(prob=fdrr$postprob, bucket=cut(fdrr$z, h1$breaks))
  pmean = mosaic::maggregate(prob~bucket, data=zcut, FUN=mean)
  pmean[is.na(pmean)] = 0
  par(new=TRUE)
  h2 = h1
  h2$density = h2$density * pmean
  plot(h2, freq=FALSE, col='blue', border='grey', axes=FALSE, xlab='', ylab='', main='', ylim=c(0, 1.05*max(h1$density)))
  lines(fdrr$x_grid, fdrr$fmix_grid, col='black')
  lines(fdrr$x_grid, fdrr$p0 * fdrr$f0_grid, col='red', lty='dotted', lwd=1)
  legend('topright', c(expression(f(z)), expression(pi[0] %.% f[0](z))),
         lty=c('solid', 'dotted'), col=c('black', 'red'), bty='n')
  if(showrug) {
    rug( fdrr$z[fdrr$FDR < Q], ticksize=0.03, col='black')
    mysub = paste0(mysub, '\nBlack rug: discoveries at FDR = ', Q)
  }
  if(showsub) title(sub=mysub)
  par(new=FALSE)
}


# Benjamini-Hochberg with two-sided p-values
BenjaminiHochberg = function(zscores, fdr_level) {
  # zscores is a vector of z scores
  # fdr_level is the desired level (e.g. 0.1) of control over FDR
  # returns a binary vector where 0=nofinding, 1=finding at given FDR level
  N = length(zscores)
  pval2 = 2*pmin(pnorm(zscores), 1- pnorm(zscores))
  cuts = (1:N)*fdr_level/N
  bhdiff = sort(pval2)-cuts
  bhcutind2 = max(which(bhdiff < 0))
  bhcut2 = sort(pval2)[bhcutind2]
  0+{pval2 <= bhcut2}
}

# Utility function for extracting error rates and a confusion matrix
GetErrorRates = function(truth, guess) {
  # truth is a binary vector saying which cases are signals
  # guess is a binary vector saying which cases are "findings" from a given procedure
  confusion_matrix = table(factor(truth, levels=c(0,1)), factor(guess, levels=c(0,1)))
  
  # Need to catch the corner case: no signals
  if( sum(truth) == 0 ) {
    true_positive_rate = 0
  } else {
    true_positive_rate = confusion_matrix[2,2]/sum(truth)
  }
  
  # Need to catch the corner case: no discoveries
  if( sum(guess) == 0 )  {
    false_discovery_rate = 0
  } else {
    false_discovery_rate = confusion_matrix[1,2]/sum(guess)
  }
  
  list(tpr = true_positive_rate, fdr = false_discovery_rate, confusion = confusion_matrix)
}


# Iteratively fit the regression piece of the PR-based FDR regression
fdrr_regress_pr = function(M0, M1, X, W_initial, maxit=2500, abstol = 1e-6, lambda=0.01) {
  stopifnot( M0 > 0, M1 > 0 )
  P = ncol(X)
  result = tryCatch({
    
    # Initialize the regression fit
    travel=1
    Xs = cbind(1,X)
    PostProb = W_initial*M1/(W_initial*M1 + (1-W_initial)*M0)
    suppressWarnings(lm1 <- SoftLogitFit(PostProb, Xs, lambda=lambda))
    betaguess = lm1$coef
    W = ilogit(drop(Xs %*% betaguess))
    PostProb = W*M1/(W*M1 + (1-W)*M0)
    oldval = lm1$value
    
    # Iterate until convergence, each time with a warm start
    passcounter = 0
    while(abs(travel/oldval) > abstol && passcounter <= maxit) {
      suppressWarnings(lm1 <- SoftLogitFit(PostProb, Xs, lambda=lambda, start=betaguess))
      newval = lm1$value
      travel = abs((oldval - newval)/(oldval + abstol))
      oldval = newval
      betaguess = lm1$coef
      W = ilogit(drop(Xs %*% betaguess))
      PostProb = W*M1/(W*M1 + (1-W)*M0)
      passcounter = passcounter + 1
    }
    if(abs(travel/oldval) > abstol) {
      mywarning = paste0('\nMaximum FDRR iteration (maxit) reached for PR method. ',
                         'Try re-running with a weaker tolerance or larger maxit.')
      warning(mywarning, immediate.=FALSE)
    }
    list(PostProb = PostProb, W = W, model = lm1)
  }, error = function(err) {
    print(err)
    print("An error was encountered in fitting the regression.  Reverting to the PR no-covariates model.")
    list(PostProb = W_initial*M1/(W_initial*M1 + (1-W_initial)*M0), W = W_initial, model = NULL)
  }, finally = {
    # No cleanup necessary
  })
  return(result)
}

# Iteratively fit the regression piece of the FDR regression using Efron's estimator
fdrr_regress_efron = function(M0, MTot, X, W_initial, maxit=500, abstol = 1e-6) {
  stopifnot( M0 > 0, MTot > 0 )
  P = ncol(X)
  result = tryCatch({
    
    # Initialize the regression fit
    travel=1
    nullcases = which((1-W_initial)*M0>MTot)	
    PostProb = pmax(0, pmin(1, 1 - (1-W_initial)*M0/MTot))
    PostProb[nullcases] = 0
    suppressWarnings(lm1 <- glm(PostProb ~ X, family=binomial))
    W = fitted(lm1)
    PostProb = pmax(0, pmin(1, 1 - (1-W)*M0/MTot))
    PostProb[nullcases] = 0
    betaguess = coef(lm1)
    
    
    # Iteratively refine estimate for beta
    travel=1
    passcounter = 0
    while(travel > abstol && passcounter <= maxit) {
      suppressWarnings(lm1 <- glm(PostProb ~ X, family=binomial, start=betaguess))
      newbeta = coef(lm1)
      W = fitted(lm1)
      PostProb = pmax(0, pmin(1, 1 - (1-W)*M0/MTot))
      PostProb[nullcases] = 0
      travel = sum(abs(betaguess-newbeta))
      betaguess = newbeta
      passcounter = passcounter + 1
    }
    
    if(travel > abstol) {
      myerror = paste0('Maximum FDRR iteration (maxit) reached for Efron\'s method. ',
                       'Try re-running with a weaker tolerance or larger maxit.')
      warning(mywarning, immediate. = TRUE)
    }
    
    list(PostProb = PostProb, W = W, model = lm1)
    
  }, error = function(err) {
    #print(err)
    print("An error was encountered in fitting the regression.  Reverting to the Efron no-covariates model.")
    list(PostProb = pmax(0, pmin(1, 1 - (1-W_initial)*M0/MTot)), W = W_initial, model = NULL)
  }, finally = {
    # No cleanup necessary
  })
  return(result)
}




# EM deconvolution for a Gaussian mixture model for mu_i with N(mu_i,1) observations
deconvolveEM = function(z, ncomps, rel.tol=1e-7, plotit=FALSE) {
  N = length(z)
  
  comp_weights = c(0.5*rep(1/ncomps, ncomps), 0.5)
  comp_means = c(quantile(z[abs(z)>2], probs=seq(.025, 0.975, length=ncomps)), 0)
  comp_variance = c(rep(2, ncomps), 1)
  
  loglike = sum(log(dnormix(z, comp_weights, comp_means, comp_variance)))
  converged=FALSE
  
  marg_like = matrix(0, nrow=N, ncol=ncomps+1)
  marg_like[,ncomps+1] = dnorm(z)
  while(!converged) {
    
    # E step
    for(j in 1:ncomps) {
      marg_like[,j] = dnorm(z, comp_means[j], sqrt(comp_variance[j]))
    }
    post_probs = scale(marg_like, center=FALSE, scale=1/comp_weights)
    post_probs = post_probs/rowSums(post_probs)
    
    # M step
    comp_weights = colSums(post_probs)/N
    for(j in 1:ncomps) {
      comp_means[j] = sum(z*post_probs[,j])/sum(post_probs[,j])
      # The component-level variances are constrained to be >=1, consistent with a deconvolution
      comp_variance[j] = max(1,  sum( {(z-comp_means[j])^2} * post_probs[,j])/sum(post_probs[,j]) )
    }
    
    loglikenew = sum(log(dnormix(z, comp_weights, comp_means, comp_variance)))
    relative_change = abs(loglikenew - loglike)/abs(loglike + rel.tol)
    converged = {relative_change < rel.tol}
    loglike = loglikenew
  }
  if(plotit) {
    par(mfrow=c(1,2))
    hist(z, 100, prob=TRUE, col='grey', border=NA)
    curve(marnormix(x, rep(1, N), comp_weights, comp_means, comp_variance-1), add=TRUE, n=1000)
    hist(z, 100, prob=TRUE, col='grey', border=NA)
    curve(dnormix(x, comp_weights[1:ncomps]/sum(comp_weights[1:ncomps]), comp_means[1:ncomps], comp_variance[1:ncomps]), add=TRUE, n=1000)
    
  }
  list(weights = comp_weights, means = comp_means, vars = comp_variance, loglike = loglike, AIC = -2*loglike + 2*{3*ncomps -1 })
}

FDRsmooth1D = function(z, x=NULL, lambda=NULL,
                       nulltype='theoretical', method='pr', stderr = NULL, control=list()) {
  # False discovery rate smoothing
  # z = vector of z scores
  # x = integer locations of entries in z along 1D chain graph
  # 		if x is missing, z is assumed to be in order, i.e. x=(1,...,N)
  # nulltype = flag for what kind of null hypothesis to assume, theoretical/empirical/heteroscedastic
  
  stopifnot(method=='pr')
  
  # Set up control parameters
  mycontrol = list()
  if(method=='pr') {
    mycontrol$reltol = 1e-6
    mycontrol$maxit = 1000
    mycontrol$gridsize = 300
    mycontrol$decay = -0.67
    mycontrol$npasses = 10
    mycontrol$densknots=10
  }
  # Overwrite with user choices
  mycontrol[(namc <- names(control))] <- control
  
  # Sample size
  n = length(z)
  
  
  # Deal with ordering of z scores
  if(missing(x)) {
    z_sort = z
    x_order = 1:n
  } else {
    x_order = order(x)
    z_sort = z[x_order]
  }
  
  
  # Estimate the marginal density
  if(method=='pr') {
    
    # Compute f0 and f1, the marginals under null and alternative for each observation
    if(nulltype=='empirical') {
      
      # Currently using my implementation of Efron's central matching estimator
      l1 = efron(z, nmids=150, df=mycontrol$densknots, nulltype=nulltype)
      mu0 = l1$mu0
      sig0 = l1$sig0
      
      prfit = prfdr(z_sort, mu0, sig0, control=mycontrol)
      fmix_grid = prfit$fmix_grid
      f0_grid = dnorm(prfit$x_grid, mu0, sig0)
      f1_grid = prfit$f1_grid
    } else if(nulltype=='heteroscedastic') {
      if(missing(stderr)) {
        stop("Must specify standard error (stderr) if assuming heteroscedastic null.")
      }
      mu0 = 0.0
      sig0 = stderr
      prfit = prfdr_het(z_sort, mu0, sig0[x_order], control=mycontrol)
      fmix_grid = NULL
      f0_grid = NULL
      f1_grid = NULL
    } else {
      mu0 = 0.0
      sig0 = 1.0
      prfit = prfdr(z_sort, mu0, sig0, control=mycontrol)
      fmix_grid = prfit$fmix_grid
      f0_grid = dnorm(prfit$x_grid, mu0, sig0)
      f1_grid = prfit$f1_grid
    }
    # Extract marginal densities and fit regression
    p0 = prfit$pi0
    f0 = dnorm(z_sort, mu0, sig0)
    f1 = prfit$f1_z
    f1zeros = which(f1 < .Machine$double.eps)
    if(length(f1zeros > 0)) {
      f1[f1zeros] = min(f1[-m1zeros]) # substitute in the smallest nonzero value
      f1 = pmax(f1, .Machine$double.eps) # shouldn't happen but just in case!
    }
    x_grid = prfit$x_grid
  }
  
  # Now run FDR smoothing
  if(missing(lambda)) {
    lambda_grid = 10^seq(2,-0.5,length=50)
  } else {
    lambda_grid = lambda
  }
  n_lambda = length(lambda_grid)
  dof_grid = rep(0, n_lambda)
  dev_grid = rep(0, n_lambda)
  theta = matrix(0, nrow=n, ncol=n_lambda)
  theta_start = rep(0, n)
  fit_list = list()
  for(i in seq_along(lambda_grid)) {
    lambda = lambda_grid[i]
    this_fit = fdrs1D(f0, f1, lambda=lambda, mycontrol$reltol, mycontrol$maxit)
    dof_grid[i] = this_fit$n_jumps
    dev_grid[i] = 2*this_fit$val
    # update warm-start value and store result
    theta_start = this_fit$beta_hat
    theta[,i] = theta_start
  }
  bic_grid = dev_grid + sqrt(n) * dof_grid
  aic_grid = dev_grid + 2 * dof_grid
  list(theta = theta, prfit = prfit,
       lambda = lambda_grid, dof = dof_grid, dev = dev_grid, bic = bic_grid, aic = aic_grid)
}


fdrs1D = function(f0, f1, lambda, reltol, maxit, beta_start = NULL) {
  # Used internally by the main wrapper function FDRsmooth1D
  # f0 = vector of densities of z[i] under null
  # f1 = vector of densities of z[i] under alternative
  
  n = length(f0)
  if(missing(beta_start)) {
    beta_hat = rep(0, n)
  } else {
    beta_hat = beta_start
  }
  
  prior_prob = ilogit(beta_hat)
  objective_old = -sum(log(prior_prob*f1 + (1-prior_prob)*f0))
  
  # Start iterates
  drift = max(1, 10*reltol)
  step_counter = 0
  while( {drift > reltol} & {step_counter <= maxit} ) {
    
    # E step
    prior_prob = ilogit(beta_hat)
    m1 = prior_prob*f1
    m0 = (1-prior_prob)*f0
    post_prob = m1/(m1+m0)
    
    # M step
    ebeta = exp(beta_hat)
    weights = ebeta/{(1+ebeta)^2}
    y = {(1+ebeta)^2}*post_prob/ebeta + beta_hat - (1+ebeta)
    weights = prior_prob*(1-prior_prob)
    y = beta_hat - (prior_prob - post_prob)/weights
    
    # Fit the 1-D fused lasso using the GLM pseudo-responses and weights
    beta_hat = fl_dp_weight(drop(y), w=drop(weights), lam=lambda)
    prior_prob = ilogit(beta_hat)
    
    objective_new = -sum(log(prior_prob*f1 + (1.0-prior_prob)*f0))
    drift = abs(objective_old - objective_new)/(abs(objective_new) + reltol)
    objective_old = objective_new
    step_counter = step_counter + 1
  }
  
  # Finish up and return fitted prior/posterior probabilities
  n_jumps = sum(abs(diff(beta_hat)) > 1e-6)
  m1 = prior_prob*f1
  m0 = (1-prior_prob)*f0
  post_prob = m1/(m1+m0)
  list(prior_prob = prior_prob, post_prob = post_prob,
       beta_hat = beta_hat, n_jumps = n_jumps, val = objective_old)
}

# April 11, 2014: migrated this functionality from a separate package over here

prfdr = function(z, mu0=0.0, sig0=1.0, control=list()) {
  
  mycontrol=list(gridsize = 500, npasses=10, decay=-0.67)
  mycontrol[(namc <- names(control))] <- control
  stopifnot(mycontrol$decay < -2/3, mycontrol$decay > -1, mycontrol$npasses >= 1)
  if(mycontrol$npasses < 5) warning("It is not recommended to use predictive recursion with < 5 passes through the data.\n")
  
  # Initial guess for alternative density and pi0
  x_grid = seq(min(z), max(z), length=mycontrol$gridsize)
  nullprob = 0.95
  # theta_guess = pmin(1, 0.1+seq(-2,2, length=mycontrol$gridsize)^2)
  theta_guess = rep(1, length(x_grid))
  theta_guess = (1.0-nullprob)* theta_guess/trapezoid(x_grid, theta_guess)
  
  # We sweep through the data npasses times in random order
  # Set up the vector of indices
  N = length(z)
  sweeporder = rep(0L, mycontrol$npasses*N)
  for(k in 1:mycontrol$npasses) {
    sweeporder[1:N + (k-1)*N] = sample(0:(N-1))
  }
  
  out1 = PredictiveRecursionFDR(z, sweeporder, x_grid, theta_guess,
                                mu0, sig0, nullprob=nullprob, decay=mycontrol$decay)
  out2 = eval_pr_dens(z, mu0, rep(sig0, N), x_grid, out1$theta_subdens)
  
  # pi(theta), the prior (normalized)
  pitheta_grid = out1$theta_subdens / (1-out1$pi0)
  
  # The mixture density
  f1_z = out2$fsignal_z
  f0_z = dnorm(z, mu0, sig0)
  fmix_z = (1.0-out1$pi0) * f1_z + out1$pi0*f0_z
  postprob = (1.0-out1$pi0) * f1_z / fmix_z
  
  out3 <- list(x_grid = x_grid, pi0 = out1$pi0, pitheta_grid = pitheta_grid, 
               f1_grid = out1$y_signal, fmix_grid = out1$y_mix,
               f0_z = f0_z, f1_z = out2$fsignal_z, fmix_z = fmix_z, postprob = postprob)
  out3
}


prfdr_het = function(z, mu0=0.0, sig0, control=list()) {
  # Predictive recursion estimate of a mixing density under heteroscedastic Gaussian error
  
  # Set up control parameters
  mycontrol=list(gridsize = 500, npasses=20, decay=-0.67)
  mycontrol[(namc <- names(control))] <- control
  stopifnot(mycontrol$decay < -.5, mycontrol$decay > -1, mycontrol$npasses >= 1)
  if(mycontrol$npasses < 5) warning("It is not recommended to use predictive recursion with < 5 passes through the data.\n")
  
  # Initial guess for alternative density and pi0
  x_grid = seq(min(z), max(z), length=mycontrol$gridsize)
  nullprob = 0.95
  # theta_guess = pmin(1, 0.1+seq(-2,2, length=mycontrol$gridsize)^2)
  theta_guess = rep(1, length(x_grid))
  theta_guess = (1.0-nullprob)* theta_guess/trapezoid(x_grid, theta_guess)
  
  # We sweep through the data npasses times in random order
  # Set up the vector of indices
  N = length(z)
  sweeporder = rep(0L, mycontrol$npasses*N)
  for(k in 1:mycontrol$npasses) {
    sweeporder[1:N + (k-1)*N] = sample(0:(N-1))
  }
  
  # Call cpp routine
  out1 = PredictiveRecursion_DifferentSigma(z, mu0, sig0, sweeporder, x_grid, theta_guess,
                                            nullprob=nullprob, decay=mycontrol$decay)
  out2 = eval_pr_dens(z, mu0, sig0, x_grid, out1$theta_subdens)
  
  out3 <- list(x_grid = x_grid, pi0 = out1$pi0, ftheta_grid = out1$theta_subdens, 
               fsignal_z=out2$fsignal_z)
  out3
}




























