library(tidyverse);
library(survival);
library(mstate);
library(etm);
library(mici);

rm(list=ls())

write_to_folder = "/Users/echase/Library/CloudStorage/OneDrive-RANDCorporation/UMich Research/mici/Code/CumincSimsFinal/";
setwd("/Users/echase/Library/CloudStorage/OneDrive-RANDCorporation/UMich Research/mici/mici_sims")

rho <- 0.2

kappa_r_list <- c(0.5, 2)
kappa_a_list <- c(0.5, 2)
n_list <- c(25, 100, 500)
c_list <- c("low", "med", "high")

fullcombs <- expand.grid(kappa_r_list, kappa_a_list, n_list, c_list)
colnames(fullcombs) <- c("kappar", "kappaa", "n", "c")
fullcombs <- dplyr::filter(fullcombs, (c != "high" | (kappar==0.5 & kappaa==0.5)), 
                    !(kappar==2 & kappaa==2))

fullcombs$setting <- case_when(
  fullcombs$kappar==0.5 & fullcombs$kappaa==0.5 & fullcombs$c=="low" ~ "A",
  fullcombs$kappar==0.5 & fullcombs$kappaa==2 & fullcombs$c=="low" ~ "B",
  fullcombs$kappar==2 & fullcombs$kappaa==0.5 & fullcombs$c=="low" ~ "C",
  fullcombs$kappar==0.5 & fullcombs$kappaa==0.5 & fullcombs$c=="med" ~ "D",
  fullcombs$kappar==0.5 & fullcombs$kappaa==2 & fullcombs$c=="med" ~ "E",
  fullcombs$kappar==2 & fullcombs$kappaa==0.5 & fullcombs$c=="med" ~ "F",
  fullcombs$kappar==0.5 & fullcombs$kappaa==0.5 & fullcombs$c=="high" ~ "G"
)

nrep <- 1000*nrow(fullcombs)

finaldat <- data.frame("Incidence" = NA,
                       "Time" = NA,
                       "SE" = NA,
                       "Lower" = NA,
                       "Upper" = NA,
                       "time" = NA,
                       "truth" = NA,
                       "Method" = NA,
                       "setting" = NA,
                       "n" = NA,
                       "censoring" = NA)

for (array_id in 1:nrep){
  set.seed(array_id)
  
  rowid <- 1 + array_id%%nrow(fullcombs)
  
  kappa_r <- fullcombs[rowid, "kappar"]
  kappa_a <- fullcombs[rowid, "kappaa"]
  n <- fullcombs[rowid, "n"]
  c <- fullcombs[rowid, "c"]
  setting <- fullcombs[rowid, "setting"]
  
  if (c=="low"){
    censor_bound <- 10
  } else if (c=="med"){
    if (kappa_r ==0.5){
      censor_bound <- 2
    } else if (kappa_r==2){
      censor_bound <- 4.5
    }
  } else {
    censor_bound <- 0.2
  }
  
  if (setting=="A" | setting=="B" | setting=="F"){
    times <- c(1, 2)
  } else if (setting=="D" | setting=="E"){
    times <- c(0.1, 0.5)
  } else if (setting=="C"){
    times <- c(2, 4)
  } else if (setting=="G"){
    times <- c(0.03, 0.06)
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
  mycurvefits$setting <- setting
  mycurvefits$n <- n
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
  mycurvefits$setting <- setting
  mycurvefits$n <- n
  mycurvefits$censoring <- c
  finaldat <- rbind(finaldat, mycurvefits)
  
  n_samps <- 1000
  shape1 <- max(0.5, 2*length(which(u$ftype==1))/length(which(u$ftype!=0)))
  shape1 <- min(shape1, 1)
  if (is.nan(shape1)){shape1 <- 1}
  shape2 <- 2-shape1
  
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
  mycurvefits$Method <- "RSI_bayes_custom"
  mycurvefits$setting <- setting
  mycurvefits$n <- n
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
  mycurvefits$Method <- "KMI_bayes_custom"
  mycurvefits$setting <- setting
  mycurvefits$n <- n
  mycurvefits$censoring <- c
  finaldat <- rbind(finaldat, mycurvefits)
  
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
  mycurvefits$Method <- "RSI_bayes_0.8"
  mycurvefits$setting <- setting
  mycurvefits$n <- n
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
  mycurvefits$Method <- "KMI_bayes_0.8"
  mycurvefits$setting <- setting
  mycurvefits$n <- n
  mycurvefits$censoring <- c
  finaldat <- rbind(finaldat, mycurvefits)
  
  shape1 <- 0.5
  shape2 <- 0.5
  
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
  mycurvefits$Method <- "RSI_bayes_0.5"
  mycurvefits$setting <- setting
  mycurvefits$n <- n
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
  mycurvefits$Method <- "KMI_bayes_0.5"
  mycurvefits$setting <- setting
  mycurvefits$n <- n
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
  mycurvefits$setting <- setting
  mycurvefits$n <- n
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
  mycurvefits$setting <- setting
  mycurvefits$n <- n
  mycurvefits$censoring <- c
  finaldat <- rbind(finaldat, mycurvefits)
  
  start.time <- Sys.time()
  source('RGI.R')
  end.time <- Sys.time()
  time.taken <- difftime(end.time, start.time, units="secs")
  
  mycurvefits$time <- time.taken
  mycurvefits$truth <- true_dat$cuminc1
  mycurvefits$Method <- "RGI"
  mycurvefits$setting <- setting
  mycurvefits$n <- n
  mycurvefits$censoring <- c
  finaldat <- rbind(finaldat, mycurvefits)
  
  start.time <- Sys.time()
  if (length(which(u$ftype==1))==0){
    mycurvefits <- data.frame("Incidence" = 0, "Time" = times, 
                              "SE" = 0, 
                              "Lower" = 0, 
                              "Upper" = 0)
  } else{
    fit <- Cuminc(u$ftime, u$ftype);
    lcprog <- summary(fit, times = new_times);
    
    mycurvefits <- data.frame("Incidence" = lcprog$pstate[,2], 
                              "Time" = times, 
                              "SE" = lcprog$std.err[,2], 
                              "Lower" = lcprog$lower[,2], 
                              "Upper" = lcprog$upper[,2])
  }
  end.time <- Sys.time()
  time.taken <- difftime(end.time, start.time, units="secs")
  
  mycurvefits$truth <- true_dat$cuminc1
  mycurvefits$time <- time.taken
  mycurvefits$Method <- "AalJo_cmp"
  mycurvefits$setting <- setting
  mycurvefits$n <- n
  mycurvefits$censoring <- c
  
  finaldat <- rbind(finaldat, mycurvefits)
  
  start.time <- Sys.time()
  if (length(which(u$ftype==1))==0){
    mycurvefits <- data.frame("Incidence" = 0, "Time" = times, 
                              "SE" = 0, 
                              "Lower" = 0, 
                              "Upper" = 0)
  } else{
    u$entry <- 0;
    fit <- etmCIF(Surv(entry, ftime, ftype != 0) ~ 1, data = u, etype = ftype);
    lcprog <- summary(fit, times = new_times, ci.fun = "cloglog")[[1]][[1]];
    mycurvefits <- data.frame("Incidence" = lcprog$P, "Time" = times, 
                              "SE" = sqrt(lcprog$var), "Lower" = lcprog$lower, 
                              "Upper" = lcprog$upper)
  }
  end.time <- Sys.time()
  time.taken <- difftime(end.time, start.time, units="secs")
  
  mycurvefits$truth <- true_dat$cuminc1
  mycurvefits$time <- time.taken
  mycurvefits$Method <- "Gaynor"
  mycurvefits$setting <- setting
  mycurvefits$n <- n
  mycurvefits$censoring <- c
  finaldat <- rbind(finaldat, mycurvefits)
  
  if (array_id %% 600 == 0){
    cat(paste0(array_id, "\n"))
  }
}

finaldat <- finaldat[-1,]

save(list = c("finaldat"), file = "/Users/echase/Desktop/cuminc_dat_final.RData")
