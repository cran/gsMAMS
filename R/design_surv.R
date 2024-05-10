#' @title  Design the clinical trial for survival outcome
#' @description This function generates the design parameters of a clinical trial for survival outcome.
#' @param m0 numeric Median survival time of control group.
#' @param alpha numeric Type I error.
#' @param beta numeric Type II error.
#' @param k numeric Number of treatment arms.
#' @param hr0 numeric Hazard ratio of ineffective treatment group vs control.
#' @param hr1 numeric Hazard ratio of effective treatment group vs control.
#' @param ta numeric Accrual time.
#' @param tf numeric Follow-up time.
#' @param kappa numeric Shape parameter (kappa=1 for exponential distribution).
#' @param eta numeric Rate of loss to follow-up.
#' @param frac numeric Vector of fractions for information time at each look.
#' @return List of cumulative number of events for each stage of combined treatment and control groups along with total number of subjects and maximum total number of events for the trial. It also provides efficacy and futility boundaries of the trial.
#' @examples
#' design_surv(m0 = 20,
#'             hr0 = 1,
#'             hr1 = 0.65,
#'             ta = 20,
#'             tf = 40,
#'             alpha = 0.05,
#'             beta = 0.1,
#'             k = 3,
#'             kappa = 1,
#'             eta = 0,
#'             frac = c(1 / 2, 1))

#' @export

design_surv <- function(m0, alpha, beta, k, hr0, hr1, ta, tf, kappa, eta, frac) {
  n <- size_surv(m0, alpha, beta, k=k, hr0=hr0, hr1=hr1, ta, tf, kappa, eta, frac)
  mat <- matrix(NA, nrow = 1, ncol = length(frac))
  rownames(mat) <- c("Cumulative number of events for combined treatment & control")
  colnames(mat) <- paste("Stage", seq_len(length(frac)))
  mat[1, ] <- ceiling(2 * n[1] * frac)

  mat1 <- matrix(NA, nrow = 2, ncol = length(frac))
  rownames(mat1) <- c("Lower bound", "Upper bound")
  colnames(mat1) <- paste("Stage", seq_len(length(frac)))
  bv <- scprt(alpha = alpha, k = k, frac = frac)
  mat1[1, ] <- bv$lshape
  mat1[2, ] <- bv$ushape
  p <- list("Sample size" = mat, "Maximum total number of events for the trial" = (k + 1) * n[1], "Total number of subjects required for the trial" = n[2] * (k + 1), "Boundary values" = mat1)
  return(p)
}
