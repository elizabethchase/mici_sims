library(dplyr)

true_dat <- data.frame("Death" = rweibull(n, shape = 1, scale = 10), 
  "Relapse" = rweibull(n, shape = kappa_r, scale = 1/rho),
  "aGVHD" = rweibull(n, shape = kappa_a, scale = 1/rho),
  "Censor" = runif(n, 0, censor_bound))

true_dat$X <- pmin(true_dat$Death, true_dat$Relapse, true_dat$aGVHD, true_dat$Censor)
  true_dat$delta <- case_when(
  true_dat$X==true_dat$Relapse ~ 1,
  true_dat$X==true_dat$Censor ~ 0,
  (true_dat$X==true_dat$aGVHD | true_dat$X==true_dat$Death) ~ 2
)

mydat <- data.frame("ftime" = true_dat$X, "ftype" = true_dat$delta)

integrand1 <- function(t){kappa_r*rho*(rho*t)^(kappa_r - 1)*
   exp(-((rho*t)^kappa_r + (rho*t)^kappa_a + 0.1*t))}
integrand2 <- function(t){kappa_a*rho*(rho*t)^(kappa_a - 1)*
  exp(-((rho*t)^kappa_r + (rho*t)^kappa_a + 0.1*t))}
integrand3 <- function(t){0.1*exp(-((rho*t)^kappa_r + (rho*t)^kappa_a + 0.1*t))}

cuminc1 <- rep(NA, length(times))
cuminc2 <- rep(NA, length(times))
cuminc3 <- rep(NA, length(times))

for (i in 1:length(times)){
  cuminc1[i] <- integrate(integrand1, lower = 0, upper = times[i])$value
  cuminc2[i] <- integrate(integrand2, lower = 0, upper = times[i])$value
  cuminc3[i] <- integrate(integrand3, lower = 0, upper = times[i])$value
}

true_dat <- data.frame("Time" = times, "cuminc1" = cuminc1, "cuminc2" = cuminc2,
"cuminc3" = cuminc3)
