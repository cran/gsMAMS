#' @title  Calculates the Sample Size for a Clinical Trial
#' @description This function calculates the sample size per arm of a clinical trial for continuous outcome.
#' @param delta0 numeric Standardized effect size in ineffective arm.
#' @param delta1 numeric Standardized effect size in effective arm.
#' @param alpha numeric Type I error.
#' @param beta numeric Type II error.
#' @param k numeric Number of treatment arms.
#' @return A numeric value indicating the sample size per arm.
#' @examples
#' size_cont(delta0 = 0.178, delta1 = 0.545, alpha = 0.05, beta = 0.1, k = 4)
#' @keywords internal
#' @noRd


size_cont <- function(delta0, delta1, alpha, beta, k) {
  r <- 1
  delta <- c(delta1, rep(delta0, k - 1))
  if (k == 1) {
    z0 <- stats::qnorm(1 - alpha)
    z1 <- stats::qnorm(1 - beta)
    n <- ceiling((z0 + z1)^2 * (1 + r) / (r * delta^2))
  }
  if (k >= 2) {
    Sigma <- matrix((1 / (1 + r)), k, k)
    diag(Sigma) <- 1
    root <- function(c) {
      alpha - (1 - mvtnorm::pmvnorm(
        lower = rep(-Inf, k), upper = rep(c, k),
        mean = rep(0, k), sigma = Sigma
      )[1])
    }
    c <- stats::uniroot(root, lower = 0, upper = 999)$root
    2
    Sigma11 <- Sigma[1:(k - 1), 1:(k - 1)]
    Sigma12 <- Sigma[1:(k - 1), k]
    Sigma21 <- Sigma[k, 1:(k - 1)]
    A <- (-1) * diag(k)
    A[, 1] <- 1
    A <- rbind(A[-1, ], A[1, ])
    B <- A %*% Sigma %*% t(A)
    root1 <- function(n) {
      mu <- sqrt(r * n / (1 + r)) * delta
      b <- as.numeric(A %*% mu)
      int <- mvtnorm::pmvnorm(
        lower = c(rep(0, k - 1), c), upper = rep(Inf, k),
        mean = b, sigma = B
      )[1]
      1 - beta - as.double(int)
    }
    n <- ceiling(stats::uniroot(root1, lower = 1, upper = 999)$root)
  }
  return(n)
}
