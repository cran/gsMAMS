#' @title  Design the clinical trial for ordinal outcome
#' @description This function generates the design parameters of a clinical trial for ordinal outcome.
#' @param alpha numeric Type I error.
#' @param beta numeric Type II error.
#' @param k numeric Number of treatment arms.
#' @param prob numeric Probability of ordinal outcomes in control group.
#' @param or0 numeric Odds ratio of ineffective treatment group vs control.
#' @param or numeric Odds ratio of effective treatment group vs control.
#' @param frac numeric Vector of fractions for information time at each look.
#' @return List of cumulative sample size for each stage of treatment and control groups along with maximum total sample size of the trial. It also provides efficacy and futility boundaries of the trial.
#' @examples
#' design_ord(alpha = 0.05,
#'            beta = 0.1,
#'            k = 4,
#'            prob = c(0.075, 0.182, 0.319, 0.243, 0.015, 0.166),
#'            or = 3.06,
#'            or0 = 1.32,
#'            frac = c(1 / 2, 1))
#' @export

design_ord <- function(alpha, beta, k, prob, or0, or, frac) {
  n <- size_ord(alpha, beta, k=k, prob, or0, or)
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
