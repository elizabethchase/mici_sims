library(tidyverse);
library(survival);
library(cmprsk);
library(etm);
library(mici);

rm(list=ls())

write_to_folder = "/Users/echase/Library/CloudStorage/OneDrive-RANDCorporation/UMich Research/mici/mici_sims/Results/";
setwd("/Users/echase/Library/CloudStorage/OneDrive-RANDCorporation/UMich Research/mici/mici_sims")

rho <- 0.2
kappa_a <- 0.5
kappa_r <- 0.5

n_list <- c(25, 100, 500)
c_list <- c("low", "high")

nrep <- 6000

finaldat <- data.frame("Incidence" = NA,
                       "Time" = NA,
                       "SE" = NA,
                       "Lower" = NA,
                       "Upper" = NA,
                       "time" = NA,
                       "truth" = NA,
                       "Method" = NA,
                       "n" = NA,
                       "m" = NA,
                       "actual_m" = NA,
                       "censoring" = NA)

for (array_id in 1:nrep){
  
n <- n_list[1 + array_id%%length(n_list)]
c <- c_list[1 + array_id%%length(c_list)]

if (c=="low"){
  censor_bound <- 10
  times <- c(1, 2)
} else {
  censor_bound <- 0.2
  times <- c(0.01, 0.02)
}

source('data_gen.R')

u <- mydat

prob_inds <- which(times > max(u$ftime[u$ftype==1]))
if (length(prob_inds)>0){
  if (length(which(u$ftype==1))==0){
    new_times <- rep(0, length(times))
  } else{
    new_times <- c(times[-prob_inds], rep(max(u$ftime[u$ftype==1]), length(prob_inds)))
  }
} else{
  new_times <- times
}

########## m = 150

m <- 150

start.time <- Sys.time()
myimps <- mici.impute(ftime = ftime, ftype = ftype, data = u, M = m, scheme = "RSI")
myres <- mici.cuminc(myimps, times = times, int.type = "Wald")
mycurvefits <- data.frame("Incidence" = myres$cuminc, "Time" = myres$time, "SE" = NA,
                          "Lower" = myres$lower, "Upper" = myres$upper)
end.time <- Sys.time()
time.taken <- difftime(end.time, start.time, units="secs")

mycurvefits$time <- time.taken
mycurvefits$truth <- true_dat$cuminc1
mycurvefits$Method <- "RSI"
mycurvefits$n <- n
mycurvefits$m <- "150"
mycurvefits$actual_m <- 150
mycurvefits$censoring <- c
finaldat <- rbind(finaldat, mycurvefits)

start.time <- Sys.time()
myimps <- mici.impute(ftime = ftime, ftype = ftype, data = u, M = m, scheme = "RSI")
myres <- mici.cuminc(myimps, times = times)
mycurvefits <- data.frame("Incidence" = myres$cuminc, "Time" = myres$time, "SE" = NA,
                          "Lower" = myres$lower, "Upper" = myres$upper)
end.time <- Sys.time()
time.taken <- difftime(end.time, start.time, units="secs")

mycurvefits$time <- time.taken
mycurvefits$truth <- true_dat$cuminc1
mycurvefits$Method <- "RSI_lott"
mycurvefits$n <- n
mycurvefits$m <- "150"
mycurvefits$actual_m <- 150
mycurvefits$censoring <- c
finaldat <- rbind(finaldat, mycurvefits)

n_samps <- 1000
shape1 <- 0.8
shape2 <- 1.2

start.time <- Sys.time()
myimps <- mici.impute(ftime = ftime, ftype = ftype, data = u, M = m, scheme = "RSI")
myres <- mici.cuminc(myimps, times = times, int.type = "Bayes", bayes_alpha = shape1,
                     bayes_beta = shape2, bayes_samps = n_samps)
mycurvefits <- data.frame("Incidence" = myres$cuminc, "Time" = myres$time, "SE" = NA,
                          "Lower" = myres$lower, "Upper" = myres$upper)
end.time <- Sys.time()
time.taken <- difftime(end.time, start.time, units="secs")

mycurvefits$time <- time.taken
mycurvefits$truth <- true_dat$cuminc1
mycurvefits$Method <- "RSI_bayes"
mycurvefits$n <- n
mycurvefits$m <- "150"
mycurvefits$actual_m <- 150
mycurvefits$censoring <- c
finaldat <- rbind(finaldat, mycurvefits)

start.time <- Sys.time()
myimps <- mici.impute(ftime = ftime, ftype = ftype, data = u, M = m, scheme = "KMI")
myres <- mici.cuminc(myimps, times = times, int.type = "Wald")
mycurvefits <- data.frame("Incidence" = myres$cuminc, "Time" = myres$time, "SE" = NA,
                          "Lower" = myres$lower, "Upper" = myres$upper)
end.time <- Sys.time()
time.taken <- difftime(end.time, start.time, units="secs")

mycurvefits$time <- time.taken
mycurvefits$truth <- true_dat$cuminc1
mycurvefits$Method <- "KMI"
mycurvefits$n <- n
mycurvefits$m <- "150"
mycurvefits$actual_m <- 150
mycurvefits$censoring <- c
finaldat <- rbind(finaldat, mycurvefits)

start.time <- Sys.time()
myimps <- mici.impute(ftime = ftime, ftype = ftype, data = u, M = m, scheme = "RSI")
myres <- mici.cuminc(myimps, times = times)
mycurvefits <- data.frame("Incidence" = myres$cuminc, "Time" = myres$time, "SE" = NA,
                          "Lower" = myres$lower, "Upper" = myres$upper)
end.time <- Sys.time()
time.taken <- difftime(end.time, start.time, units="secs")

mycurvefits$time <- time.taken
mycurvefits$truth <- true_dat$cuminc1
mycurvefits$Method <- "KMI_lott"
mycurvefits$n <- n
mycurvefits$m <- "150"
mycurvefits$actual_m <- 150
mycurvefits$censoring <- c
finaldat <- rbind(finaldat, mycurvefits)

start.time <- Sys.time()
myimps <- mici.impute(ftime = ftime, ftype = ftype, data = u, M = m, scheme = "KMI")
myres <- mici.cuminc(myimps, times = times, int.type = "Bayes", bayes_alpha = shape1, 
                     bayes_beta = shape2, bayes_samps = n_samps)
mycurvefits <- data.frame("Incidence" = myres$cuminc, "Time" = myres$time, "SE" = NA,
                          "Lower" = myres$lower, "Upper" = myres$upper)
end.time <- Sys.time()
time.taken <- difftime(end.time, start.time, units="secs")

mycurvefits$time <- time.taken
mycurvefits$truth <- true_dat$cuminc1
mycurvefits$Method <- "KMI_bayes"
mycurvefits$n <- n
mycurvefits$m <- "150"
mycurvefits$actual_m <- 150
mycurvefits$censoring <- c
finaldat <- rbind(finaldat, mycurvefits)

start.time <- Sys.time()
source('RGI.R')
end.time <- Sys.time()
time.taken <- difftime(end.time, start.time, units="secs")

mycurvefits$time <- time.taken
mycurvefits$truth <- true_dat$cuminc1
mycurvefits$Method <- "RGI"
mycurvefits$n <- n
mycurvefits$m <- "150"
mycurvefits$actual_m <- 150
mycurvefits$censoring <- c
finaldat <- rbind(finaldat, mycurvefits)

########## m = 50 recommendations

m <- 50

start.time <- Sys.time()
myimps <- mici.impute(ftime = ftime, ftype = ftype, data = u, M = m, scheme = "RSI")
myres <- mici.cuminc(myimps, times = times, int.type = "Wald")
mycurvefits <- data.frame("Incidence" = myres$cuminc, "Time" = myres$time, "SE" = NA,
                          "Lower" = myres$lower, "Upper" = myres$upper)
end.time <- Sys.time()
time.taken <- difftime(end.time, start.time, units="secs")

mycurvefits$time <- time.taken
mycurvefits$truth <- true_dat$cuminc1
mycurvefits$Method <- "RSI"
mycurvefits$n <- n
mycurvefits$m <- "50"
mycurvefits$actual_m <- 50
mycurvefits$censoring <- c
finaldat <- rbind(finaldat, mycurvefits)

start.time <- Sys.time()
myimps <- mici.impute(ftime = ftime, ftype = ftype, data = u, M = m, scheme = "RSI")
myres <- mici.cuminc(myimps, times = times)
mycurvefits <- data.frame("Incidence" = myres$cuminc, "Time" = myres$time, "SE" = NA,
                          "Lower" = myres$lower, "Upper" = myres$upper)
end.time <- Sys.time()
time.taken <- difftime(end.time, start.time, units="secs")

mycurvefits$time <- time.taken
mycurvefits$truth <- true_dat$cuminc1
mycurvefits$Method <- "RSI_lott"
mycurvefits$n <- n
mycurvefits$m <- "50"
mycurvefits$actual_m <- 50
mycurvefits$censoring <- c
finaldat <- rbind(finaldat, mycurvefits)

start.time <- Sys.time()
myimps <- mici.impute(ftime = ftime, ftype = ftype, data = u, M = m, scheme = "RSI")
myres <- mici.cuminc(myimps, times = times, int.type = "Bayes", bayes_alpha = shape1,
                     bayes_beta = shape2, bayes_samps = n_samps)
mycurvefits <- data.frame("Incidence" = myres$cuminc, "Time" = myres$time, "SE" = NA,
                          "Lower" = myres$lower, "Upper" = myres$upper)
end.time <- Sys.time()
time.taken <- difftime(end.time, start.time, units="secs")

mycurvefits$time <- time.taken
mycurvefits$truth <- true_dat$cuminc1
mycurvefits$Method <- "RSI_bayes"
mycurvefits$n <- n
mycurvefits$m <- "50"
mycurvefits$actual_m <- 50
mycurvefits$censoring <- c
finaldat <- rbind(finaldat, mycurvefits)

start.time <- Sys.time()
myimps <- mici.impute(ftime = ftime, ftype = ftype, data = u, M = m, scheme = "KMI")
myres <- mici.cuminc(myimps, times = times, int.type = "Wald")
mycurvefits <- data.frame("Incidence" = myres$cuminc, "Time" = myres$time, "SE" = NA,
                          "Lower" = myres$lower, "Upper" = myres$upper)
end.time <- Sys.time()
time.taken <- difftime(end.time, start.time, units="secs")

mycurvefits$time <- time.taken
mycurvefits$truth <- true_dat$cuminc1
mycurvefits$Method <- "KMI"
mycurvefits$n <- n
mycurvefits$m <- "50"
mycurvefits$actual_m <- 50
mycurvefits$censoring <- c
finaldat <- rbind(finaldat, mycurvefits)

start.time <- Sys.time()
myimps <- mici.impute(ftime = ftime, ftype = ftype, data = u, M = m, scheme = "RSI")
myres <- mici.cuminc(myimps, times = times)
mycurvefits <- data.frame("Incidence" = myres$cuminc, "Time" = myres$time, "SE" = NA,
                          "Lower" = myres$lower, "Upper" = myres$upper)
end.time <- Sys.time()
time.taken <- difftime(end.time, start.time, units="secs")

mycurvefits$time <- time.taken
mycurvefits$truth <- true_dat$cuminc1
mycurvefits$Method <- "KMI_lott"
mycurvefits$n <- n
mycurvefits$m <- "50"
mycurvefits$actual_m <- 50
mycurvefits$censoring <- c
finaldat <- rbind(finaldat, mycurvefits)

start.time <- Sys.time()
myimps <- mici.impute(ftime = ftime, ftype = ftype, data = u, M = m, scheme = "KMI")
myres <- mici.cuminc(myimps, times = times, int.type = "Bayes", bayes_alpha = shape1, 
                     bayes_beta = shape2, bayes_samps = n_samps)
mycurvefits <- data.frame("Incidence" = myres$cuminc, "Time" = myres$time, "SE" = NA,
                          "Lower" = myres$lower, "Upper" = myres$upper)
end.time <- Sys.time()
time.taken <- difftime(end.time, start.time, units="secs")

mycurvefits$time <- time.taken
mycurvefits$truth <- true_dat$cuminc1
mycurvefits$Method <- "KMI_bayes"
mycurvefits$n <- n
mycurvefits$m <- "50"
mycurvefits$actual_m <- 50
mycurvefits$censoring <- c
finaldat <- rbind(finaldat, mycurvefits)

start.time <- Sys.time()
source('RGI.R')
end.time <- Sys.time()
time.taken <- difftime(end.time, start.time, units="secs")

mycurvefits$time <- time.taken
mycurvefits$truth <- true_dat$cuminc1
mycurvefits$Method <- "RGI"
mycurvefits$n <- n
mycurvefits$m <- "50"
mycurvefits$actual_m <- 50
mycurvefits$censoring <- c
finaldat <- rbind(finaldat, mycurvefits)

########## m = 10 recommendations

m <- 10

start.time <- Sys.time()
myimps <- mici.impute(ftime = ftime, ftype = ftype, data = u, M = m, scheme = "RSI")
myres <- mici.cuminc(myimps, times = times, int.type = "Wald")
mycurvefits <- data.frame("Incidence" = myres$cuminc, "Time" = myres$time, "SE" = NA,
                          "Lower" = myres$lower, "Upper" = myres$upper)
end.time <- Sys.time()
time.taken <- difftime(end.time, start.time, units="secs")

mycurvefits$time <- time.taken
mycurvefits$truth <- true_dat$cuminc1
mycurvefits$Method <- "RSI"
mycurvefits$n <- n
mycurvefits$m <- "10"
mycurvefits$actual_m <- 10
mycurvefits$censoring <- c
finaldat <- rbind(finaldat, mycurvefits)

start.time <- Sys.time()
myimps <- mici.impute(ftime = ftime, ftype = ftype, data = u, M = m, scheme = "RSI")
myres <- mici.cuminc(myimps, times = times)
mycurvefits <- data.frame("Incidence" = myres$cuminc, "Time" = myres$time, "SE" = NA,
                          "Lower" = myres$lower, "Upper" = myres$upper)
end.time <- Sys.time()
time.taken <- difftime(end.time, start.time, units="secs")

mycurvefits$time <- time.taken
mycurvefits$truth <- true_dat$cuminc1
mycurvefits$Method <- "RSI_lott"
mycurvefits$n <- n
mycurvefits$m <- "10"
mycurvefits$actual_m <- 10
mycurvefits$censoring <- c
finaldat <- rbind(finaldat, mycurvefits)

start.time <- Sys.time()
myimps <- mici.impute(ftime = ftime, ftype = ftype, data = u, M = m, scheme = "RSI")
myres <- mici.cuminc(myimps, times = times, int.type = "Bayes", bayes_alpha = shape1,
                     bayes_beta = shape2, bayes_samps = n_samps)
mycurvefits <- data.frame("Incidence" = myres$cuminc, "Time" = myres$time, "SE" = NA,
                          "Lower" = myres$lower, "Upper" = myres$upper)
end.time <- Sys.time()
time.taken <- difftime(end.time, start.time, units="secs")

mycurvefits$time <- time.taken
mycurvefits$truth <- true_dat$cuminc1
mycurvefits$Method <- "RSI_bayes"
mycurvefits$n <- n
mycurvefits$m <- "10"
mycurvefits$actual_m <- 10
mycurvefits$censoring <- c
finaldat <- rbind(finaldat, mycurvefits)

start.time <- Sys.time()
myimps <- mici.impute(ftime = ftime, ftype = ftype, data = u, M = m, scheme = "KMI")
myres <- mici.cuminc(myimps, times = times, int.type = "Wald")
mycurvefits <- data.frame("Incidence" = myres$cuminc, "Time" = myres$time, "SE" = NA,
                          "Lower" = myres$lower, "Upper" = myres$upper)
end.time <- Sys.time()
time.taken <- difftime(end.time, start.time, units="secs")

mycurvefits$time <- time.taken
mycurvefits$truth <- true_dat$cuminc1
mycurvefits$Method <- "KMI"
mycurvefits$n <- n
mycurvefits$m <- "10"
mycurvefits$actual_m <- 10
mycurvefits$censoring <- c
finaldat <- rbind(finaldat, mycurvefits)

start.time <- Sys.time()
myimps <- mici.impute(ftime = ftime, ftype = ftype, data = u, M = m, scheme = "RSI")
myres <- mici.cuminc(myimps, times = times)
mycurvefits <- data.frame("Incidence" = myres$cuminc, "Time" = myres$time, "SE" = NA,
                          "Lower" = myres$lower, "Upper" = myres$upper)
end.time <- Sys.time()
time.taken <- difftime(end.time, start.time, units="secs")

mycurvefits$time <- time.taken
mycurvefits$truth <- true_dat$cuminc1
mycurvefits$Method <- "KMI_lott"
mycurvefits$n <- n
mycurvefits$m <- "10"
mycurvefits$actual_m <- 10
mycurvefits$censoring <- c
finaldat <- rbind(finaldat, mycurvefits)

start.time <- Sys.time()
myimps <- mici.impute(ftime = ftime, ftype = ftype, data = u, M = m, scheme = "KMI")
myres <- mici.cuminc(myimps, times = times, int.type = "Bayes", bayes_alpha = shape1, 
                     bayes_beta = shape2, bayes_samps = n_samps)
mycurvefits <- data.frame("Incidence" = myres$cuminc, "Time" = myres$time, "SE" = NA,
                          "Lower" = myres$lower, "Upper" = myres$upper)
end.time <- Sys.time()
time.taken <- difftime(end.time, start.time, units="secs")

mycurvefits$time <- time.taken
mycurvefits$truth <- true_dat$cuminc1
mycurvefits$Method <- "KMI_bayes"
mycurvefits$n <- n
mycurvefits$m <- "10"
mycurvefits$actual_m <- 10
mycurvefits$censoring <- c
finaldat <- rbind(finaldat, mycurvefits)

start.time <- Sys.time()
source('RGI.R')
end.time <- Sys.time()
time.taken <- difftime(end.time, start.time, units="secs")

mycurvefits$time <- time.taken
mycurvefits$truth <- true_dat$cuminc1
mycurvefits$Method <- "RGI"
mycurvefits$n <- n
mycurvefits$m <- "10"
mycurvefits$actual_m <- 10
mycurvefits$censoring <- c
finaldat <- rbind(finaldat, mycurvefits)

########## m = 1 recommendations

m <- 1

start.time <- Sys.time()
myimps <- mici.impute(ftime = ftime, ftype = ftype, data = u, M = m, scheme = "RSI")
myres <- mici.cuminc(myimps, times = times, int.type = "Wald")
mycurvefits <- data.frame("Incidence" = myres$cuminc, "Time" = myres$time, "SE" = NA,
                          "Lower" = myres$lower, "Upper" = myres$upper)
end.time <- Sys.time()
time.taken <- difftime(end.time, start.time, units="secs")

mycurvefits$time <- time.taken
mycurvefits$truth <- true_dat$cuminc1
mycurvefits$Method <- "RSI"
mycurvefits$n <- n
mycurvefits$m <- "1"
mycurvefits$actual_m <- 1
mycurvefits$censoring <- c
finaldat <- rbind(finaldat, mycurvefits)

start.time <- Sys.time()
myimps <- mici.impute(ftime = ftime, ftype = ftype, data = u, M = m, scheme = "RSI")
myres <- mici.cuminc(myimps, times = times)
mycurvefits <- data.frame("Incidence" = myres$cuminc, "Time" = myres$time, "SE" = NA,
                          "Lower" = myres$lower, "Upper" = myres$upper)
end.time <- Sys.time()
time.taken <- difftime(end.time, start.time, units="secs")

mycurvefits$time <- time.taken
mycurvefits$truth <- true_dat$cuminc1
mycurvefits$Method <- "RSI_lott"
mycurvefits$n <- n
mycurvefits$m <- "1"
mycurvefits$actual_m <- 1
mycurvefits$censoring <- c
finaldat <- rbind(finaldat, mycurvefits)

start.time <- Sys.time()
myimps <- mici.impute(ftime = ftime, ftype = ftype, data = u, M = m, scheme = "RSI")
myres <- mici.cuminc(myimps, times = times, int.type = "Bayes", bayes_alpha = shape1,
                     bayes_beta = shape2, bayes_samps = n_samps)
mycurvefits <- data.frame("Incidence" = myres$cuminc, "Time" = myres$time, "SE" = NA,
                          "Lower" = myres$lower, "Upper" = myres$upper)
end.time <- Sys.time()
time.taken <- difftime(end.time, start.time, units="secs")

mycurvefits$time <- time.taken
mycurvefits$truth <- true_dat$cuminc1
mycurvefits$Method <- "RSI_bayes"
mycurvefits$n <- n
mycurvefits$m <- "1"
mycurvefits$actual_m <- 1
mycurvefits$censoring <- c
finaldat <- rbind(finaldat, mycurvefits)

start.time <- Sys.time()
myimps <- mici.impute(ftime = ftime, ftype = ftype, data = u, M = m, scheme = "KMI")
myres <- mici.cuminc(myimps, times = times, int.type = "Wald")
mycurvefits <- data.frame("Incidence" = myres$cuminc, "Time" = myres$time, "SE" = NA,
                          "Lower" = myres$lower, "Upper" = myres$upper)
end.time <- Sys.time()
time.taken <- difftime(end.time, start.time, units="secs")

mycurvefits$time <- time.taken
mycurvefits$truth <- true_dat$cuminc1
mycurvefits$Method <- "KMI"
mycurvefits$n <- n
mycurvefits$m <- "1"
mycurvefits$actual_m <- 1
mycurvefits$censoring <- c
finaldat <- rbind(finaldat, mycurvefits)

start.time <- Sys.time()
myimps <- mici.impute(ftime = ftime, ftype = ftype, data = u, M = m, scheme = "RSI")
myres <- mici.cuminc(myimps, times = times)
mycurvefits <- data.frame("Incidence" = myres$cuminc, "Time" = myres$time, "SE" = NA,
                          "Lower" = myres$lower, "Upper" = myres$upper)
end.time <- Sys.time()
time.taken <- difftime(end.time, start.time, units="secs")

mycurvefits$time <- time.taken
mycurvefits$truth <- true_dat$cuminc1
mycurvefits$Method <- "KMI_lott"
mycurvefits$n <- n
mycurvefits$m <- "1"
mycurvefits$actual_m <- 1
mycurvefits$censoring <- c
finaldat <- rbind(finaldat, mycurvefits)

start.time <- Sys.time()
myimps <- mici.impute(ftime = ftime, ftype = ftype, data = u, M = m, scheme = "KMI")
myres <- mici.cuminc(myimps, times = times, int.type = "Bayes", bayes_alpha = shape1, 
                     bayes_beta = shape2, bayes_samps = n_samps)
mycurvefits <- data.frame("Incidence" = myres$cuminc, "Time" = myres$time, "SE" = NA,
                          "Lower" = myres$lower, "Upper" = myres$upper)
end.time <- Sys.time()
time.taken <- difftime(end.time, start.time, units="secs")

mycurvefits$time <- time.taken
mycurvefits$truth <- true_dat$cuminc1
mycurvefits$Method <- "KMI_bayes"
mycurvefits$n <- n
mycurvefits$m <- "1"
mycurvefits$actual_m <- 1
mycurvefits$censoring <- c
finaldat <- rbind(finaldat, mycurvefits)

start.time <- Sys.time()
source('RGI.R')
end.time <- Sys.time()
time.taken <- difftime(end.time, start.time, units="secs")

mycurvefits$time <- time.taken
mycurvefits$truth <- true_dat$cuminc1
mycurvefits$Method <- "RGI"
mycurvefits$n <- n
mycurvefits$m <- "1"
mycurvefits$actual_m <- 1
mycurvefits$censoring <- c
finaldat <- rbind(finaldat, mycurvefits)

########## FCI recommendations

u2 <- filter(u, ftime <= max(u$ftime[u$ftype != 0]))

if (nrow(u2) > 0){
  censor_model <- survfit(Surv(ftime, as.numeric(ftype == 0)) ~ 1, data = u2)
  allcause_model <- survfit(Surv(ftime, as.numeric(ftype != 0)) ~ 1, data = u2)
  
  censor_time <- summary(censor_model)$table["median"]
  ac_time <- summary(allcause_model)$table["median"]
  
  if (is.na(ac_time)){
    ac_time <- max(u2$ftime[u2$ftype!= 0])
  }
  
  if (is.na(censor_time)){
    censor_time <- max(u2$ftime[u2$ftype==0])
  }
  
  fmi_imps <- ceiling((ac_time/censor_time)*length(which(u$ftype==0)))
  
  fmi <- max(c(fmi_imps, 10))
  m <- min(c(fmi, 200))
} else{
  m <- 1
}

if (is.na(m)){m <- 1}

start.time <- Sys.time()
myimps <- mici.impute(ftime = ftime, ftype = ftype, data = u, M = m, scheme = "RSI")
myres <- mici.cuminc(myimps, times = times, int.type = "Wald")
mycurvefits <- data.frame("Incidence" = myres$cuminc, "Time" = myres$time, "SE" = NA,
                          "Lower" = myres$lower, "Upper" = myres$upper)
end.time <- Sys.time()
time.taken <- difftime(end.time, start.time, units="secs")

mycurvefits$time <- time.taken
mycurvefits$truth <- true_dat$cuminc1
mycurvefits$Method <- "RSI"
mycurvefits$n <- n
mycurvefits$m <- "FCI"
mycurvefits$actual_m <- m
mycurvefits$censoring <- c
finaldat <- rbind(finaldat, mycurvefits)

start.time <- Sys.time()
myimps <- mici.impute(ftime = ftime, ftype = ftype, data = u, M = m, scheme = "RSI")
myres <- mici.cuminc(myimps, times = times)
mycurvefits <- data.frame("Incidence" = myres$cuminc, "Time" = myres$time, "SE" = NA,
                          "Lower" = myres$lower, "Upper" = myres$upper)
end.time <- Sys.time()
time.taken <- difftime(end.time, start.time, units="secs")

mycurvefits$time <- time.taken
mycurvefits$truth <- true_dat$cuminc1
mycurvefits$Method <- "RSI_lott"
mycurvefits$n <- n
mycurvefits$m <- "FCI"
mycurvefits$actual_m <- m
mycurvefits$censoring <- c
finaldat <- rbind(finaldat, mycurvefits)

start.time <- Sys.time()
myimps <- mici.impute(ftime = ftime, ftype = ftype, data = u, M = m, scheme = "RSI")
myres <- mici.cuminc(myimps, times = times, int.type = "Bayes", bayes_alpha = shape1,
                     bayes_beta = shape2, bayes_samps = n_samps)
mycurvefits <- data.frame("Incidence" = myres$cuminc, "Time" = myres$time, "SE" = NA,
                          "Lower" = myres$lower, "Upper" = myres$upper)
end.time <- Sys.time()
time.taken <- difftime(end.time, start.time, units="secs")

mycurvefits$time <- time.taken
mycurvefits$truth <- true_dat$cuminc1
mycurvefits$Method <- "RSI_bayes"
mycurvefits$n <- n
mycurvefits$m <- "FCI"
mycurvefits$actual_m <- m
mycurvefits$censoring <- c
finaldat <- rbind(finaldat, mycurvefits)

start.time <- Sys.time()
myimps <- mici.impute(ftime = ftime, ftype = ftype, data = u, M = m, scheme = "KMI")
myres <- mici.cuminc(myimps, times = times, int.type = "Wald")
mycurvefits <- data.frame("Incidence" = myres$cuminc, "Time" = myres$time, "SE" = NA,
                          "Lower" = myres$lower, "Upper" = myres$upper)
end.time <- Sys.time()
time.taken <- difftime(end.time, start.time, units="secs")

mycurvefits$time <- time.taken
mycurvefits$truth <- true_dat$cuminc1
mycurvefits$Method <- "KMI"
mycurvefits$n <- n
mycurvefits$m <- "FCI"
mycurvefits$actual_m <- m
mycurvefits$censoring <- c
finaldat <- rbind(finaldat, mycurvefits)

start.time <- Sys.time()
myimps <- mici.impute(ftime = ftime, ftype = ftype, data = u, M = m, scheme = "RSI")
myres <- mici.cuminc(myimps, times = times)
mycurvefits <- data.frame("Incidence" = myres$cuminc, "Time" = myres$time, "SE" = NA,
                          "Lower" = myres$lower, "Upper" = myres$upper)
end.time <- Sys.time()
time.taken <- difftime(end.time, start.time, units="secs")

mycurvefits$time <- time.taken
mycurvefits$truth <- true_dat$cuminc1
mycurvefits$Method <- "KMI_lott"
mycurvefits$n <- n
mycurvefits$m <- "FCI"
mycurvefits$actual_m <- m
mycurvefits$censoring <- c
finaldat <- rbind(finaldat, mycurvefits)

start.time <- Sys.time()
myimps <- mici.impute(ftime = ftime, ftype = ftype, data = u, M = m, scheme = "KMI")
myres <- mici.cuminc(myimps, times = times, int.type = "Bayes", bayes_alpha = shape1, 
                     bayes_beta = shape2, bayes_samps = n_samps)
mycurvefits <- data.frame("Incidence" = myres$cuminc, "Time" = myres$time, "SE" = NA,
                          "Lower" = myres$lower, "Upper" = myres$upper)
end.time <- Sys.time()
time.taken <- difftime(end.time, start.time, units="secs")

mycurvefits$time <- time.taken
mycurvefits$truth <- true_dat$cuminc1
mycurvefits$Method <- "KMI_bayes"
mycurvefits$n <- n
mycurvefits$m <- "FCI"
mycurvefits$actual_m <- m
mycurvefits$censoring <- c
finaldat <- rbind(finaldat, mycurvefits)

start.time <- Sys.time()
source('RGI.R')
end.time <- Sys.time()
time.taken <- difftime(end.time, start.time, units="secs")

mycurvefits$time <- time.taken
mycurvefits$truth <- true_dat$cuminc1
mycurvefits$Method <- "RGI"
mycurvefits$n <- n
mycurvefits$m <- "FCI"
mycurvefits$actual_m <- m
mycurvefits$censoring <- c
finaldat <- rbind(finaldat, mycurvefits)

## For comparison, the Aalen-Johansen estimator (which does not depend on m):
start.time <- Sys.time()
if (length(which(u$ftype==1))==0){
  mycurvefits <- data.frame("Incidence" = 0, "Time" = times, 
                            "SE" = 0, 
                            "Lower" = 0, 
                            "Upper" = 0)
} else{
  fit <- cuminc(u$ftime, u$ftype);
  lcprog <- timepoints(fit, times = new_times);
  
  
  mycurvefits <- data.frame("Incidence" = lcprog$est[1,], "Time" = times, 
                            "SE" = lcprog$var[1,], 
                            "Lower" = lcprog$est[1,]^(exp(-1.96*sqrt(lcprog$var[1,])/(lcprog$est[1,]*log(lcprog$est[1,])))), 
                            "Upper" = lcprog$est[1,]^(exp(1.96*sqrt(lcprog$var[1,])/(lcprog$est[1,]*log(lcprog$est[1,])))))
}
end.time <- Sys.time()
time.taken <- difftime(end.time, start.time, units="secs")

mycurvefits$truth <- true_dat$cuminc1
mycurvefits$time <- time.taken
mycurvefits$Method <- "AalJo_cmp"
mycurvefits$n <- n
mycurvefits$censoring <- c
mycurvefits$m <- NA_character_
mycurvefits$actual_m <- NA_real_

finaldat <- rbind(finaldat, mycurvefits)

if (array_id %% 600 == 0){
  cat(paste0(array_id, "\n"))
}

}

save(list=c("finaldat"), 
     file = paste0(write_to_folder, "cuminc_imp_recs.RData"))