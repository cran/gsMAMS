#' @title Provides operating characteristics of group sequential MAMS trial for survival outcome under null hypothesis
#' @description Computes FWER and other characteristics for group-sequential MAMS trial for survival outcome.
#' @param m0 numeric Median survival time in control group.
#' @param alpha numeric Type I error.
#' @param beta numeric Type II error.
#' @param p numeric Number of treatment arms.
#' @param frac numeric Vector of fractions for information time at each look.
#' @param hr0 numeric Hazard ratio of ineffective treatment group vs control.
#' @param hr1 numeric Hazard ratio of effective treatment group vs control.
#' @param nsim numeric Number of simulations.
#' @param ta numeric Accrual time.
#' @param tf numeric Follow-up time.
#' @param kappa numeric Shape parameter (Kappa=1 for exponential distribution).
#' @param eta numeric  Rate of loss to follow-up.
#' @param seed numeric Random seed number.
#' @return A list of FWER, stage-wise type I error, stopping probability, probability of futility, average number of events happened per arm, average duration of trial.
#' @examples
#' op_fwer_surv(m0 = 20,
#'              alpha = 0.05,
#'              beta = 0.1,
#'              p = 4,
#'              frac = c(1 / 2, 1),
#'              hr0 = 1,
#'              hr1 = 0.75,
#'              nsim = 12,
#'              ta = 40,
#'              tf = 20,
#'              kappa = 1,
#'              eta = 0,
#'              seed = 12)
#' @export



op_fwer_surv <- function(m0, alpha, beta, p, frac, hr0, hr1, nsim, ta, tf, kappa, eta, seed) {
  HR0 <- hr0
  HR1 <- hr1
  K <- p
  if (K <= 1) {
    stop("K should be greater than 1.")
  }
  if (length(frac) == 1) {
    stop("The length of frac should be greater than 1.")
  }
  j <- length(frac)

  bound <- scprt(alpha = alpha, k = K, frac = frac)



  hr0 <- HR0 # exp(-delta0)
  hr1 <- HR1 # exp(-delta1)

  lambda0 <- log(2) / m0^kappa

  lambda <- numeric()



  for (i in 1:K) {
    lambda[i] <- lambda0^(1 / kappa)
  }


  tau <- ta + tf

  scale <- vector()
  scale[1] <- 1 / lambda0^(1 / kappa)


  for (i in 2:(K + 1)) {
    scale[i] <- 1 / lambda[i - 1]
  }





  n <- size_surv(m0 = m0, hr0 = HR0, hr1 = HR1, ta = ta, tf = tf, k = K, beta = beta, alpha = alpha, kappa = kappa, eta = eta, frac = frac)[2]
  d <- size_surv(m0 = m0, hr0 = HR0, hr1 = HR1, ta = ta, tf = tf, k = K, beta = beta, alpha = alpha, kappa = kappa, eta = eta, frac = frac)[1]
  tstar <- 1 / j


  tau <- ta + tf



  s1 <- 0
  s2 <- numeric(length = j)

  frac <- frac
  d1 <- numeric(length = j)
  for (i in 1:j) {
    d1[i] <- ceiling(d * frac[i])
  }


   
  a <- data.frame()
  for (i in 1:K) {
    a <- rbind.data.frame(a, bound$lshape)
  }

  b <- data.frame()
  for (i in 1:K) {
    b <- rbind.data.frame(b, bound$ushape)
  }


  sp <- numeric()
  pf <- numeric()
  smn <- numeric()
  smd <- numeric()
  dur <- numeric()
   
  set.seed(seed)
  for (e in 1:nsim) {
    Q <- rep(0, j)
    # ASN=0
    datagen <- NULL
    for (i in 1:(K + 1)) {
      w <- stats::rweibull(n, kappa, scale[i])
      u <- stats::runif(n, 0, ta) ## generate accrual time
      if (eta != 0) {
        g <- stats::rexp(n, rate = eta)
      }
      if (eta == 0) {
        g <- tau - u
      }
      x <- pmax(0, pmin(w, tau - u, g)) ## observed survival time
      cens <- as.numeric(w < pmin((tau - u), g)) ## censoring indicator
      group <- rep((i - 1), n)
      datagen[[i]] <- cbind(x, cens, group, u)
    }


    data <- NULL
    for (i in 1:K) {
      dat <- rbind.data.frame(datagen[[1]], datagen[[i + 1]])
      colnames(dat) <- c("time", "event", "group", "accrualtime")
  
      dat$calendar_T <- dat$accrualtime + dat$time
      dat <- dat[order(dat$calendar_T), ]
      dat$cum_event <- cumsum(dat$event)
      data[[i]] <- dat
    }


    loc <- matrix(NA, K, j)
    z <- matrix(NA, K, j)
    var_result <- matrix(NA, K, j)
    N_result <- matrix(NA, K + 1, j)
    d_result <- matrix(NA, K + 1, j)
    calendarT_look <- matrix(NA, K, j)
    stagensizeT <- matrix(NA, K, j)
    stagensizeC <- matrix(NA, K, j)
    stagedsizeT <- matrix(NA, K, j)
    stagedsizeC <- matrix(NA, K, j)


    for (h in (1:j)) {
      for (k in (1:K)) {
        if (h < j) {
      
          dlook <- ceiling(tstar * h * 2 * d)
          loc <- min(which(data[[k]]$cum_event == dlook))
          calendarT_look[k, h] <- data[[k]]$calendar_T[loc]
          data_compare_part1 <- data[[k]][1:loc, ]
          data_compare_part2 <- data[[k]][(loc + 1):(2 * n), ]
          data_compare_part2 <- data_compare_part2[
            data_compare_part2$accrualtime < data_compare_part1$calendar_T[loc],
          ]
          data_compare_part2$event <- 0
          data_compare_part2$time <- data_compare_part1$calendar_T[loc] - data_compare_part2$accrualtime
    
          data_compare <- rbind(data_compare_part1, data_compare_part2)
          stagensizeT[k, h] <- nrow(data_compare[data_compare$group != 0, ])
          stagensizeC[k, h] <- nrow(data_compare[data_compare$group == 0, ])

          stagedsizeT[k, h] <- nrow(data_compare[data_compare$group != 0 & data_compare$event == 1, ])
          stagedsizeC[k, h] <- nrow(data_compare[data_compare$group == 0 & data_compare$event == 1, ])
        }

      
        if (h == j) {
          calendarT_look[k, h] <- tau
          data_compare <- data[[k]]
          stagensizeT[k, h] <- nrow(data_compare[data_compare$group != 0, ])
          stagensizeC[k, h] <- nrow(data_compare[data_compare$group == 0, ])
          stagedsizeT[k, h] <- nrow(data_compare[data_compare$group != 0 & data_compare$event == 1, ])
          stagedsizeC[k, h] <- nrow(data_compare[data_compare$group == 0 & data_compare$event == 1, ])
        }
         
        temp <- survival::survdiff(survival::Surv(time, event) ~ group, data = data_compare)
        z[k, h] <- sign(temp$obs[1] - temp$exp[1]) * sqrt(temp$chisq)
        
      }
    }
    

    w <- K
    g <- data.frame(matrix(ncol = w, nrow = 0))
    mp <- data.frame(matrix(ncol = w, nrow = 0))
    for (q in 1:(length(g))) {
      # j<-3
      p <- numeric(length = (j - 1) * 3)
      v <- numeric(length = (j - 1) * 3)
      k <- seq(1, 25, by = 3)[1:(j - 1)]
      sk <- list(a = 1, b = c(2, 3))
      

      for (i in 2:(j)) {
        if (i == 2) {
      
          p[k[i - 1]] <- z[q, (i - 1)] < a[q, (i - 1)]
          p[k[i - 1] + 1] <- z[q, (i - 1)] > a[q, (i - 1)] & z[q, (i - 1)] < b[q, (i - 1)]
          p[k[i - 1] + 2] <- z[q, (i)] < b[q, (i)]

          v[k[i - 1]] <- z[q, (i - 1)] > b[q, (i - 1)]
          v[k[i - 1] + 1] <- z[q, (i - 1)] > a[q, (i - 1)] & z[q, (i - 1)] < b[q, (i - 1)]
          v[k[i - 1] + 2] <- z[q, (i)] > b[q, (i)]

          if (i == j) {
            g[1, q] <- 0
            for (o in seq_len(length(sk))) {
              g[1, q] <- g[1, q] + prod(p[sk[[o]]])
              mp[o, q] <- prod(v[sk[[o]]])
            }
          }
          next
        }
      

        tk <- sk[[length(sk)]]

        sk[[length(sk)]][length(sk[[length(sk)]])] <- sk[[length(sk)]][length(sk[[length(sk)]])] + 1

        sk[[length(sk) + 1]] <- tk
        sk[[length(sk)]][length(sk[[length(sk)]])] <- sk[[length(sk)]][length(sk[[length(sk)]])] + 2
        sk[[length(sk)]][length(sk[[length(sk)]]) + 1] <- sk[[length(sk)]][length(sk[[length(sk)]])] + 1

        p[k[i - 1]] <- z[q, (i - 1)] < a[q, (i - 1)]
        p[k[i - 1] + 1] <- z[q, (i - 1)] > a[q, (i - 1)] & z[q, (i - 1)] < b[q, (i - 1)]
        p[k[i - 1] + 2] <- z[q, (i)] < b[q, (i)]

        v[k[i - 1]] <- z[q, (i - 1)] > b[q, (i - 1)]
        v[k[i - 1] + 1] <- z[q, (i - 1)] > a[q, (i - 1)] & z[q, (i - 1)] < b[q, (i - 1)]
        v[k[i - 1] + 2] <- z[q, (i)] > b[q, (i)]

        if (i == j) {
          g[1, q] <- 0
          for (o in seq_len(length(sk))) {
            g[1, q] <- g[1, q] + prod(p[sk[[o]]])
            mp[o, q] <- prod(v[sk[[o]]])
          }
        }
      }
    }


    if (prod(g[1, ]) >= 1) {
      s1 <- s1 + 1
    }

    for (y in 1:j) {
      s2[y] <- s2[y] + ifelse(sum(mp[y, ]) >= 1, 1, 0)
    }

    


    

    stopprob <- rep(0, j)
    probfut <- rep(0, j)
    stop <- rep(0, j)
    stopF <- matrix(NA, K, j)
    stopE <- matrix(NA, K, j)
    Fflag <- rep(0, j)
    Eflag <- rep(0, j)

    for (k in 1:K) {
      for (h in 1:j) {
        stopF[k, h] <- z[k, h] <= a[k, h]
        stopE[k, h] <- z[k, h] >= b[k, h]
      }
    }


    for (k in 1:K) {
      if (any(stopF[k, ] == TRUE)) {
        s <- min(which(stopF[k, ] == TRUE))
        stopF[k, s:j] <- TRUE
      }
    }


    for (h in 1:j) {
      if (sum(stopF[, h] == TRUE) == K) {
        Fflag[h] <- 1
      }
      if (sum(stopE[, h] == TRUE) >= 1) {
        Eflag[h] <- 1
      }
      stop[h] <- Eflag[h] + Fflag[h] # stop is rep (0,J) originally, if we stop at stage 3, we will have (0,0,1,0), thus stop records at which stage we stop for current simulation
    }



    stopstage <- min(which(stop >= 1)) # the stage we should stop, no matter it is due to efficacy or futility
    stopFstage <- ifelse(all(Fflag == 0), 0, min(which(Fflag == 1))) # the stage we  stop due to futility

    if (stopFstage == stopstage) {
      probfut[stopFstage] <- 1
    }

    stopprob[stopstage] <- 1

    samplesize <- rep(0, K + 1)
    samplesized <- rep(0, K + 1)

    for (k in 1:K) {
      stopFstage1 <- ifelse(any(stopF[k, ] == TRUE), min(which(stopF[k, ] == TRUE)), Inf)

      stopstagefork <- min(stopstage, stopFstage1)
      samplesize[k] <- stagensizeT[k, stopstagefork]
      samplesized[k] <- stagedsizeT[k, stopstagefork]
    }
    samplesize[K + 1] <- max(stagensizeC[, stopstage])
    samplesized[K + 1] <- max(stagedsizeC[, stopstage])



    duration <- max(calendarT_look[, stopstage])




    sp <- rbind.data.frame(stopprob, sp)
    pf <- rbind.data.frame(probfut, pf)
    smn[e] <- sum(samplesize) / (K + 1)
    smd[e] <- sum(samplesized) / (K + 1)
    dur[e] <- duration
  }



  names(sp) <- paste0("look", 1:j)
  names(pf) <- paste0("look", 1:j)

  fwer <- round(1 - (s1 / nsim), 3)
  fwg <- rbind.data.frame(s2 / nsim)
  names(fwg) <- paste0("look", 1:j)


  p <- list(FWER = fwer, "Stagewise FWER" = colMeans(fwg), "Stopping probability under null" = colMeans(sp), "Probability of futility under null" = colMeans(pf), "Average number of events happened per arm under null" = mean(smd), "Average duration of trial(months)" = mean(dur))
  return(p)
}
