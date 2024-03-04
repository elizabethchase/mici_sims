library(tidyverse);
library(survival);
library(cmprsk);
library(etm);
library(mici);

rm(list=ls())

write_to_folder = "/Users/echase/Library/CloudStorage/OneDrive-RANDCorporation/UMich Research/mici/mici_sims/Results/";
setwd("/Users/echase/Library/CloudStorage/OneDrive-RANDCorporation/UMich Research/mici/mici_sims")
nrep <- 1000

finaldat <- data.frame("Incidence" = NA, 
                       "Time" = NA,
                       "SE" = NA,
                       "Lower" = NA,
                       "Upper" = NA,
                       "Method" = NA,
                       "sim_id" = NA)

for (array_id in 1:nrep){
  set.seed(array_id)
  
  ftime <- sample(c(1:10), size = 100, replace = TRUE, prob = c(0.2, 0.2, 0.2, 0.2, 0.04, 0.04, 0.04, 0.04, 0.02, 0.02))
  ftype <- sample(c(0, 1, 2), size = 100, replace = TRUE, prob = c(0.2, 0.5, 0.3))
  
  u <- data.frame("ftime" = ftime, "ftype" = ftype)
  times <- c(3, 8)
  if (max(u$ftime) < max(times)){
    badinds <- which(times > max(u$ftime))
    new_times <- rep(NA, length(times))
    new_times[badinds] <- times[min(badinds)-1]
    new_times[-badinds] <- times[-badinds]
  } else{
    new_times <- times
  }
  
  m <- 1000
  
  start.time <- Sys.time()
  myimps <- mici.impute(ftime = ftime, ftype = ftype, data = u, M = m, scheme = "RSI")
  myres <- mici.cuminc(myimps, times = times)
  mycurvefits <- data.frame("Incidence" = myres$cuminc, "Time" = myres$time, "SE" = NA,
                            "Lower" = myres$lower, "Upper" = myres$upper)
  end.time <- Sys.time()
  time.taken <- difftime(end.time, start.time, units="secs")
  
  mycurvefits$Method <- "RSI"
  mycurvefits$sim_id <- array_id
  finaldat <- rbind(finaldat, mycurvefits)
  
  start.time <- Sys.time()
  myimps <- mici.impute(ftime = ftime, ftype = ftype, data = u, M = m, scheme = "RSI")
  myres <- mici.cuminc(myimps, times = times)
  mycurvefits <- data.frame("Incidence" = myres$cuminc, "Time" = myres$time, "SE" = NA,
                            "Lower" = myres$lower, "Upper" = myres$upper)
  end.time <- Sys.time()
  time.taken <- difftime(end.time, start.time, units="secs")
  
  mycurvefits$Method <- "KMI"
  mycurvefits$sim_id <- array_id
  finaldat <- rbind(finaldat, mycurvefits)
  
  start.time <- Sys.time()
  source('RGI.R')
  end.time <- Sys.time()
  time.taken <- difftime(end.time, start.time, units="secs")
  
  mycurvefits$Method <- "RGI"
  mycurvefits$sim_id <- array_id
  finaldat <- rbind(finaldat, mycurvefits)
  
  start.time <- Sys.time()
  source('RGI_correct.R')
  end.time <- Sys.time()
  time.taken <- difftime(end.time, start.time, units="secs")
  
  mycurvefits$Method <- "RGI_correct"
  mycurvefits$sim_id <- array_id
  finaldat <- rbind(finaldat, mycurvefits)
  
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
  
  mycurvefits$Method <- "AalJo_cmp"
  mycurvefits$sim_id <- array_id
  
  finaldat <- rbind(finaldat, mycurvefits)
  
}

finaldat <- finaldat[-1,]

save(list=c("finaldat"), 
     file = paste0(write_to_folder, "cuminc_dat_ties.RData"))
