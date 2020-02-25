#' FDR Smoothing
#'
#' This function implements Tansey et al.'s (2018) FDR smoothing algorithm.
#'
#' @param dat Data frame with columns for probe name, z-statistics, chromosomal 
#'   position, and chromosome.
#' @param probe String denoting column name for probes (CpGs, SNPs, etc.).
#' @param z String denoting column name for z-scores.
#' @param pos String denoting column name for chromosomal positions.
#' @param chr String denoting column name for chromosomes.
#' @param nulltype How should the null distribution be estimated? Choose 
#'   \code{"empirical"} for Efron's central-matching method (default). Choose
#'   \code{"theoretical"} for a standard normal null. 
#' @param nlambda Length of the lambda sequence for the fused lasso subroutine.
#' @param tol Convergence tolerance for the expectation maximization (EM) 
#'   algorithm.
#' @param maxit Maximum number of iterations for the EM algorithm.
#' @param parallel Process in parallel? Only relevant if data spans multiple 
#'   chromosomes. If \code{TRUE}, backend must be registered beforehand.
#'
#' @details
#' FDR smoothing is an empirical Bayes method for exploiting spatial structure 
#' in large multiple-testing problems. The method automatically finds 
#' spatially localized regions of significant test statistics. It then relaxes 
#' the threshold of statistical significance within these regions and tightens 
#' it elsewhere, in a manner that controls the overall false discovery rate at a 
#' given level. This results in increased power and cleaner spatial separation 
#' of signals from noise. The approach requires solving a nonstandard 
#' high-dimensional optimization problem, for which an efficient 
#' augmented-Lagrangian algorithm is implemented. See (Tansey et al., 2018) 
#' for details.
#'
#' @references
#' Efron, B. (2004). \href{https://bit.ly/2kYx5AR}{Large-Scale Simultaneous 
#' Hypothesis Testing: The Choice of a Null Hypothesis}. \emph{JASA, 99}(465),
#' 96-104.
#' 
#' Newton, M.A. (2002). \href{https://bit.ly/2kYjCce}{On a Nonparametric 
#' Recursive Estimator of the Mixing Distribution}. \emph{Sankhyā: The Indian
#' Journal of Statistics, Series A, 64}(2), 306-322. 
#' 
#' Tansey, W., Koyejo, O., Poldrack, R.A., & Scott, J.G. (2018). 
#' \href{https://bit.ly/2msEdWH}{False Discovery Rate Smoothing}. \emph{JASA, 
#' 113}(523), 1156-1171.
#' 
#' Tibshirani, R., Saunders, M., Rosset, S., Zhu, J., & Knight, K. (2005). 
#' \href{https://stanford.io/2lYzEmJ}{Sparsity and Smoothness via the Fused
#' Lasso}. \emph{J. R. Statist. Soc. B, 67}(1), 91-108.
#' 
#' @examples
#' # Import data
#' data('DNAm')
#' 
#' # Set seed
#' set.seed(123)
#' 
#' # Run FDR smoothing
#' res <- smoothFDR(DNAm, probe = 'cpg', parallel = FALSE)
#' 
#' # Compare q-values to Benjamini-Hochberg estimates
#' sum(res$BH_q.value <= 0.05)
#' sum(res$q.value <= 0.05)
#' 
#' @export
#' @importFrom data.table data.table
#' @import foreach
#' @import dplyr
#' 


smoothFDR <- function(dat, 
                      probe = 'probe', 
                      z = 'z', 
                      pos = 'pos', 
                      chr = 'chr', 
                      nulltype = 'empirical', 
                      nlambda = 30, 
                      tol = 1e-6, 
                      maxit = 100, 
                      parallel = TRUE) {
  # Check for errors, etc.
  if (!is.data.frame(dat)) {
    stop('dat must be a data.frame.')
  }
  if (!probe %in% colnames(dat)) {
    stop('Column \"', probe, '\" not found.')
  }
  if (!z %in% colnames(dat)) {
    stop('Column \"', z, '\" not found.')
  }
  if (!pos %in% colnames(dat)) {
    stop('Column \"', pos, '\" not found.')
  }
  if (!chr %in% colnames(dat)) {
    stop('Column \"', chr, '\" not found.')
  }
  # Prepare data
  dat <- data.table(
      'Idx' = seq_len(nrow(dat)), 
    'probe' = dat[[probe]], 
        'z' = dat[[z]],
      'pos' = dat[[pos]], 
      'chr' = dat[[chr]]
  )
  # Execute in parallel or serial
  if (parallel) {
    out <- foreach(ch = dat[, unique(chr)], .combine = rbind) %dopar%
      smooth_fdr(probe = dat[chr == ch, probe], z = dat[chr == ch, z], 
                 pos = dat[chr == ch, pos], chr = ch, nulltype = nulltype,
                 nlambda = nlambda, tol = tol, maxit = maxit)
  } else {
    out <- foreach(ch = dat[, unique(chr)], .combine = rbind) %do%
      smooth_fdr(probe = dat[chr == ch, probe], z = dat[chr == ch, z], 
                 pos = dat[chr == ch, pos], chr = ch, nulltype = nulltype,
                 nlambda = nlambda, tol = tol, maxit = maxit)
  }
  # Rearrange and export
  dat %>% 
    select(Idx, probe) %>%
    inner_join(out, by = 'probe') %>% 
    arrange(Idx) %>%
    select(-Idx) %>%
    return(.)
}


