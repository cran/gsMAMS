#' @title  Calculates the Sample Size for a Clinical Trial
#' @description This function calculates the sample size per arm of a clinical trial of ordinal outcome.
#' @param alpha numeric Type I error.
#' @param beta numeric Type II error.
#' @param k numeric Number of treatment arms.
#' @param prob numeric Probability of ordinal outcomes in control group.
#' @param or0 numeric Odds ratio of ineffective treatment group vs control.
#' @param or numeric Odds ratio of effective treatment group vs control.
#' @return A numeric value indicating the sample size per arm.
#' @examples
#' size_ord(prob = c(0.075, 0.182, 0.319, 0.243, 0.015, 0.166), or = 3.06, or0 = 1.32, alpha = 0.05, beta = 0.1, k = 4)
#' @keywords internal
#' @noRd

size_ord <- function(alpha, beta, k, prob, or0, or) {
   if (length(prob) > 2) {
    p0 <- prob
    Q0 <- cumsum(p0)
    Qk <- 1 / (1 + ((1 - Q0) / Q0) / or)
    pk <- c(Qk[1], diff(Qk, lag = 1))
    pbar <- (p0 + k * pk) / (k + 1)
    q <- (1 - sum(pbar^3)) / 3
    sigma <- 1 / sqrt(q)
  }
  if (length(prob) == 2) {
    p0 <- prob[1]
    p1 <- (1 + (1 - p0) / (p0 * or))^(-1)
    pbar <- c((p0 + p1) / 2, 1 - (p0 + p1) / 2)
    2
    q <- (1 - sum(pbar^3)) / 3
    sigma <- 1 / sqrt(q)
  }
  delta <- c(log(or), rep(log(or0), k - 1))
  if (k == 1) {
    V <- (stats::qnorm(1 - alpha) + stats::qnorm(1 - beta))^2 / log(or)^2
    n <- ceiling(2 * V / q) ## sample size per arm ###
    return(n)
  }
  if (k >= 2) {
    Sigma <- matrix(0.5, k, k)
    diag(Sigma) <- 1
    root <- function(c) {
      alpha - (1 - mvtnorm::pmvnorm(lower = rep(-Inf, k), upper = rep(c, k), mean = rep(0, k), sigma = Sigma)[1])
    }
    c <- stats::uniroot(root, lower = 0, upper = 999)$root
    Sigma11 <- Sigma[1:(k - 1), 1:(k - 1)]
    Sigma12 <- Sigma[1:(k - 1), k]
    Sigma21 <- Sigma[k, 1:(k - 1)]
    A <- (-1) * diag(k)
    A[, 1] <- 1
    A <- rbind(A[-1, ], A[1, ])
    B <- A %*% Sigma %*% t(A)
    root1 <- function(n) {
      mu <- sqrt(n / (2 * sigma^2)) * delta
      b <- as.numeric(A %*% mu)
      int <- mvtnorm::pmvnorm(lower = c(rep(0, k - 1), c), upper = rep(Inf, k), mean = b, sigma = B)[1]
      3
      1 - beta - as.double(int)
    }
    n <- ceiling(stats::uniroot(root1, lower = 1, upper = 999)$root)
  }
  return(n)
}
