#' Log Rank Test
#'
#' A function to carry out the log rank test.
#'
#' @param time Numeric vector indicating the time of event or censoring.
#' @param event Numeric vector indicating whether the event has occurred.
#' @param group Factor vector indicating the group membership.
#' @return Z test statistic.
#' @keywords internal
#' @noRd
#' @examples
#' logranktest(time = time, event = event, group = group)
logranktest <- function(time, event, group) {
  n <- length(time)
  group <- factor(group)
  Ag <- stats::aggregate(event,
    by = list(time = time, group = group),
    FUN = sum, drop = FALSE
  )
  Ag$x <- ifelse(is.na(Ag$x), 0, Ag$x)
  tab <- data.frame(
    time = Ag$time[Ag$group == levels(group)[1]],
    event1 = Ag$x[Ag$group == levels(group)[1]],
    event2 = Ag$x[Ag$group == levels(group)[2]]
  )
  tab$atrisk1 <- NA
  tab$atrisk2 <- NA

  tab$atrisk1 <- unlist(lapply(tab$time, function(x) sum(time[group == levels(group)[1]] >= x)))
  tab$atrisk2 <- unlist(lapply(tab$time, function(x) sum(time[group == levels(group)[2]] >= x)))
  tab$atrisk <- tab$atrisk1 + tab$atrisk2
  tab$event <- tab$event1 + tab$event2
  D <- tab[tab$event > 0, ]
  D$expected1 <- D$event * D$atrisk1 / D$atrisk 
  D$diff1 <- D$event1 - D$expected1
  D$var <- D$atrisk1 * D$event * D$atrisk2 * (D$atrisk - D$event) / (D$atrisk^2 *
    (D$atrisk - 1))
  D <- D[D$atrisk1 > 0 & D$atrisk2 > 0, ]
  z <- sum(D$diff1) / sqrt(sum(D$var)) 
  return(z)
}
