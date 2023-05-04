#' Log Rank Test
#' 
#' A function to carry out the log rank test.
#' 
#' @param time Numeric vector indicating the time of event or censoring.
#' @param event Numeric vector indicating whether the event has occurred.
#' @param group Factor vector indicating the group membership. 
#' @return Z test statistic.
#' @keywords internal
#' @import stats
#' @noRd
#' @examples 
#' logranktest(time=time, event=event, group=group)

logranktest<-function (time, event, group)
{
  
  n <- length(time)
  #ng <- table(group)
  group <- factor(group)
  Ag <- aggregate(event, by = list(time = time, group = group),
                  FUN = sum, drop = FALSE)
  Ag$x <- ifelse(is.na(Ag$x), 0, Ag$x)
  tab <- data.frame(time = Ag$time[Ag$group == levels(group)[1]],
                    event1 = Ag$x[Ag$group == levels(group)[1]],
                    event2 = Ag$x[Ag$group == levels(group)[2]])
  tab$atrisk1 <- NA
  tab$atrisk2 <- NA
  
  #browser()
  #for (i in 1:dim(tab)[1]) {
  tab$atrisk1 <-unlist(lapply(tab$time,function(x) sum(time[group==levels(group)[1]]>=x))) 
  #tab$atrisk1[i]<-sum(time[group == unique(group)[1]] >=tab$time[i])
  tab$atrisk2<-unlist(lapply(tab$time,function(x) sum(time[group==levels(group)[2]]>=x)))
  # tab$atrisk2[i]<-sum(time[group == unique(group)[2]] >=tab$time[i])
  #}
  #nz <- dim(tab)[1]
  tab$atrisk <- tab$atrisk1 + tab$atrisk2
  tab$event <- tab$event1 + tab$event2
  D <- tab[tab$event > 0, ]
  D$expected1 <- D$event * D$atrisk1/D$atrisk
  #D$expected2 <- D$event * D$atrisk2/D$atrisk
  D$diff1 <- D$event1 - D$expected1
  D$var <- D$atrisk1 * D$event * D$atrisk2 * (D$atrisk - D$event)/(D$atrisk^2 *
                                                                     (D$atrisk - 1))
  #D$S <- cumprod((D$atrisk - D$event)/D$atrisk)
  #Dvoll <- D
  D <- D[D$atrisk1 > 0 & D$atrisk2 > 0, ]
  z <- sum( D$diff1)/sqrt(sum(D$var))
  #Chisq <- (sum( D$diff1))^2/sum( D$var)
  #27
  #df <- nlevels(group) - 1
  #p <- 1 - pchisq(Chisq, df = df)
  #out = list(D = D, test = data.frame( z, Chisq, df, p),
  #          var = sum(D$var), obs = c(sum(D$event1),
  #                                   sum(D$event2)), exp = c(sum(D$expected1), sum(D$expected2)),
  #        n = c(sum(D$atrisk1[1]), sum(D$atrisk2[1])))
  #class(out) = "logrank"
  return(z)
}
