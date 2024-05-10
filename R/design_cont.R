#' @title  Design the clinical trial for continuous outcome
#' @description This function generates the design parameters of a clinical trial for continuous outcome.
#' @param delta0 numeric Standardized effect size in ineffective arm.
#' @param delta1 numeric Standardized effect size in effective arm.
#' @param alpha numeric Type I error.
#' @param beta numeric Type II error.
#' @param k numeric Number of treatment arms.
#' @param frac numeric Vector of fractions for information time at each look.
#' @return List of cumulative sample size for each stage of treatment and control groups along with maximum total sample size of the trial. It also provides efficacy and futility boundaries of the trial.
#' @examples
#' design_cont(delta0 = 0.178, delta1 = 0.545, alpha = 0.05, beta = 0.1, k = 4, frac = c(1 / 2, 1))
#' @export


design_cont <- function(delta0, delta1, alpha, beta, k, frac) {
  n <- size_cont(delta0 = delta0, delta1 = delta1, alpha = alpha, beta = beta, k = k)
  mat <- matrix(NA, nrow = 2, ncol = length(frac))
  rownames(mat) <- c("Cumulative sample size for treatment group", "Cumulative sample size for control group")
  colnames(mat) <- paste("Stage", seq_len(length(frac)))
  mat[1, ] <- ceiling(n * frac)
  mat[2, ] <- ceiling(n * frac)

  mat1 <- matrix(NA, nrow = 2, ncol = length(frac))
  rownames(mat1) <- c("Lower bound", "Upper bound")
  colnames(mat1) <- paste("Stage", seq_len(length(frac)))
  bv <- scprt(alpha = alpha, k = k, frac = frac)
  mat1[1, ] <- bv$lshape
  mat1[2, ] <- bv$ushape
  p <- list("Sample size" = mat, "Maximum total sample size for the trial" = (k + 1) * n, "Boundary values" = mat1)
  return(p)
}
