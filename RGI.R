library(survival)
library(dplyr)

resample <- function(x, ...)
  x[sample.int(length(x), ...)]

dc = u[u$ftype != 2, ] ## complete cases
dd = u[u$ftype == 2, ] ##cases that need imputation

if (length(which(u$ftype == 2)) > 0) {
  xt = dd$ftime ################################################
  
  n = nrow(u) # total number cases
  b = v = NULL # objects to store results
  log_b <- log_v <- NULL
  
  g = summary(survfit(Surv(ftime, ftype == 0) ~ 1, data = u))
  w = g$time
  gm = g$surv[length(g$surv)]
  wp = -diff(c(1, g$surv))
  if (length(gm) == 0) {
    wp <- 1
    w <- max(u$ftime) + 1
  } else if (gm > 0) {
    # add mass point to censoring distribution, if needed
    wp = c(wp, gm)
    w = c(w, max(u$ftime) + 1)
  }
  
  for (j in 1:m) {
    # b). obtain the jth imputed dataset, randomly sampling
    # censoring times from KM estimate of censoring dist
    cts = NULL
    for (jj in 1:length(xt)) {
      sub = w > xt[jj]
      if (length(w[sub]) == 0) {
        cts[jj] = xt[jj]
      }
      else if (length(w[sub]) == 1) {
        cts[jj] = w[sub]
      }
      else
        cts[jj] = sample(w[sub], 1, replace = TRUE, prob = wp[sub])
    }
    ## replace competing risks failure (’failure type=2’) ## by imputed censoring times
    dd$ftime = cts
    ipd = rbind(dd, dc) ## imputed full dataset (’censoring complete data’)
    outs = summary(survfit(Surv(ftime, ftype == 1) ~ 1, data = ipd), times = new_times)
    b = rbind(b, 1 - outs$surv)
    v = rbind(v, outs$std.err ^ 2)
    log_b <- rbind(log_b, log(-log(1 - outs$surv)))
    log_v <- rbind(log_v, (outs$std.err / ((1 - outs$surv) * log(1 - outs$surv))) ^ 2)
  }
  
  if (m > 1){
    est.cuminc2 = apply(b, 2, mean)
    va = apply(v, 2, mean)
    var.coef = va + (1 + 1 / m) * apply(b, 2, var) # estimated variance-covariance matrix
    
    va_log <- apply(log_v, 2, mean)
    var.log = va_log + (1 + 1 / m) * apply(log_b, 2, var)
  } else{
    est.cuminc2 <- as.numeric(b)
    var.coef <- as.numeric(v)
    var.log <- as.numeric(log_v)
  }
} else{
  outs = summary(survfit(Surv(ftime, ftype == 1) ~ 1, data = u), times = new_times)
  est.cuminc2 <- 1-outs$surv
  var.coef <- outs$std.err ^ 2
  var.log <- (outs$std.err / ((1 - outs$surv) * log(1 - outs$surv))) ^ 2
}
#

mycurvefits <- data.frame("Incidence" = est.cuminc2, "Time" = times, "SE" = sqrt(var.coef),
                          "Lower" = est.cuminc2^(exp(1.96*sqrt(var.log))), 
                          "Upper" = est.cuminc2^(exp(-1.96*sqrt(var.log))))
