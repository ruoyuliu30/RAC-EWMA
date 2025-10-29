

library(logitnorm)
library(penPHcure)
library(MASS)
library(foreach)
library(doParallel)
library(smoothHR)


### Extract the dataset:
data(whas500)
whas500$t0 <- rep(0, dim(whas500)[1])

# Data from 1997 to 2000 is used as the risk adjustment training set;
# Data from 1999 to 2000 is used as Phase-I data to obtain UCL and LCL.
data_train <- whas500[which(whas500$year != 3),]
data_test <- whas500[which(whas500$year == 3),]
data_phase1 <- whas500[which(whas500$year == 2),]
data_phase2 <- whas500[which(whas500$year == 3),]


### Fit the mixture cure model
mcm <- penPHcure(Surv(time = t0, time2 = lenfol, event = fstat) 
                 ~ age+gender+hr+bmi+chf, 
                 cureform = ~ hr+bmi+sysbp+diasbp+chf, 
                 data = data_train, standardize = FALSE, 
                 ties = "breslow")

# Obtain the fitting results
fit_results_phase1 <- predict(mcm, data_phase1)
survival_prob_phase1 <- fit_results_phase1$SURV
cure_prob_phase1 <- fit_results_phase1$CURE
cumulative_hazard_phase1 <- -log(survival_prob_phase1)
estimate_cure_status_phase1 <- 
  data_phase1$fstat + 
  (1-data_phase1$fstat)*cure_prob_phase1*survival_prob_phase1/(1-cure_prob_phase1+cure_prob_phase1*survival_prob_phase1)

  
fit_results_phase2 <- predict(mcm, data_phase2)
survival_prob_phase2 <- fit_results_phase2$SURV
cure_prob_phase2 <- fit_results_phase2$CURE
cumulative_hazard_phase2 <- -log(survival_prob_phase2)
estimate_cure_status_phase2 <- 
  data_phase2$fstat + 
  (1-data_phase2$fstat)*cure_prob_phase2*survival_prob_phase2/(1-cure_prob_phase2+cure_prob_phase2*survival_prob_phase2)


### RAC-EWMA
RACC <- function(cure_prob, estimate_cure_status, 
                 smooth_parameter, ucl, lcl, 
                 momitoring_direction)
  # upwards when momitoring_direction = 1
  # downwards when momitoring_direction = 2
{
  z <- 0; i <- 1
  while (i > 0){
    nn <- sample(1:(length(cure_prob)), 1)
    # Randomly sample a patient to join the study
    
    
    yn <- (estimate_cure_status[nn]-cure_prob[nn])^2 - 
      cure_prob[nn]*(1-cure_prob[nn])
    z <- (1-smooth_parameter)*z + smooth_parameter*yn
    
    if (momitoring_direction == 1){
      if (z >= ucl){
        return(i)
      } else {
        i <- i+1
        next
      }
    } else {
      if (z <= lcl){
        return(i)
      } else {
        i <- i+1
        next
      }
    }
  }
}

RACL <- function(data, cure_prob, cumulative_hazard, 
                 estimate_cure_status, smooth_parameter, 
                 ucl, lcl, momitoring_direction)
  # upwards when momitoring_direction = 1
  # downwards when momitoring_direction = 2
{
  z <- 0; i <- 1
  while (i > 0){
    nn <- sample(1:(length(cure_prob)), 1)
    # Randomly sample a patient to join the study
    
    yn <- (estimate_cure_status[nn]*
             (data$fstat[nn]-cumulative_hazard[nn]))^2 - 
      estimate_cure_status[nn]*cumulative_hazard[nn]
    z <- (1-smooth_parameter)*z + smooth_parameter*yn

    if (momitoring_direction == 1){
      if (z >= ucl){
        return(i)
      } else {
        i <- i+1
        next
      }
    } else {
      if (z <= lcl){
        return(i)
      } else {
        i <- i+1
        next
      }
    }
  }
}


### Monte Carlo for UCL/LCL
# RACC's UCL
numCores <- detectCores()
registerDoParallel(numCores)
up_rls1 <- foreach(o = 1:10000, .combine = "c", 
                  .packages = c("doParallel", "foreach", 
                                "MASS", "logitnorm")) %dopar% 
  {RACC(cure_prob = cure_prob_phase1, 
        estimate_cure_status = estimate_cure_status_phase1, 
        smooth_parameter = lmd1, ucl = hc1, lcl = hc2, 
        momitoring_direction = 1)}
stopImplicitCluster()

up_arl1 <- mean(up_rls1) # Average run length for upward shitfs

# RACL's UCL
numCores <- detectCores()
registerDoParallel(numCores)
up_rls2 <- foreach(o = 1:10000, .combine = "c", 
                   .packages = c("doParallel", "foreach", 
                                 "MASS", "logitnorm")) %dopar% 
  {RACL(data = data_phase1, cure_prob = cure_prob_phase1, 
        cumulative_hazard = cumulative_hazard_phase1, 
        estimate_cure_status = estimate_cure_status_phase1, 
        smooth_parameter = lmd2, ucl = hl1, lcl = hl2, 
        momitoring_direction = 1)}
stopImplicitCluster()

up_arl2 <- mean(up_rls2) # Average run length for upward shitfs

# If arl is very close to the preset ARL0, 
# then h is the desired UCL.
# Otherwise, adjust h until the above requirements are met.


# RACC's LCL
numCores <- detectCores()
registerDoParallel(numCores)
down_rls1 <- foreach(o = 1:10000, .combine = "c", 
                   .packages = c("doParallel", "foreach", 
                                 "MASS", "logitnorm")) %dopar% 
  {RACC(cure_prob = cure_prob_phase1, 
        estimate_cure_status = estimate_cure_status_phase1, 
        smooth_parameter = lmd1, ucl = hc1, lcl = hc2, 
        momitoring_direction = 2)}
stopImplicitCluster()

down_arl1 <- mean(down_rls1) # Average run length for downward shitfs

# RACL's LCL
numCores <- detectCores()
registerDoParallel(numCores)
down_rls2 <- foreach(o = 1:10000, .combine = "c", 
                   .packages = c("doParallel", "foreach", 
                                 "MASS", "logitnorm")) %dopar% 
  {RACL(data = data_phase1, cure_prob = cure_prob_phase1, 
        cumulative_hazard = cumulative_hazard_phase1, 
        estimate_cure_status = estimate_cure_status_phase1, 
        smooth_parameter = lmd2, ucl = hl1, lcl = hl2, 
        momitoring_direction = 2)}
stopImplicitCluster()

down_arl2 <- mean(down_rls2) # Average run length for downward shitfs

# If arl is very close to the preset ARL0, 
# then h is the desired UCL.
# Otherwise, adjust h until the above requirements are met.


### Monitor Phase-II data
RACCm <- function(cure_prob, estimate_cure_status)
  # upwards when momitoring_direction = 1
  # downwards when momitoring_direction = 2
{
  z <- 0
  zs <- c() # Chart statistic set
  
  for (i in 1:(length(cure_prob))){
    yn <- (estimate_cure_status[i]-cure_prob[i])^2 - 
      cure_prob[i]*(1-cure_prob[i])
    z <- (1-smooth_parameter)*z + smooth_parameter*yn
    zs <- c(zs, z)
  }
  return(zs)
}

RACLm <- function(data, cure_prob, cumulative_hazard, 
                 estimate_cure_status, smooth_parameter)
  # upwards when momitoring_direction = 1
  # downwards when momitoring_direction = 2
{
  z <- 0
  zs <- c() # Chart statistic set
  
  for (i in 1:(length(cure_prob))){
    yn <- (estimate_cure_status[i]*
             (data$fstat[i]-cumulative_hazard[i]))^2 - 
      estimate_cure_status[i]*cumulative_hazard[i]
    z <- (1-smooth_parameter)*z + smooth_parameter*yn
    zs <- c(zs, z)
  }
  return(zs)
}


# Obtain Phase-II's Chart statistic sets
RACC_chart_statistic <- RACCm(
  cure_prob = cure_prob_phase2, 
  estimate_cure_status = estimate_cure_status_phase2, 
  smooth_parameter = lmd1)

RACL_chart_statistic <- RACLm(
  data = data_phase2, cure_prob = cure_prob_phase2, 
  cumulative_hazard = cumulative_hazard_phase2, 
  cure_prob = cure_prob_phase2, 
  estimate_cure_status = estimate_cure_status_phase2, 
  smooth_parameter = lmd2)


# Alarm information
RACC_alarm_up <- which(RACC_chart_statistic >= hc1)
RACC_alarm_down <- which(RACC_chart_statistic <= hc2)
RACL_alarm_up <- which(RACC_chart_statistic >= hl1)
RACL_alarm_down <- which(RACC_chart_statistic <= hl1)

# The serial number of the first sample in the alarm sample list
# is the first alarm position




