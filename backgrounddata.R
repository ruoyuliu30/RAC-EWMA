

### All parameter values in the following code have been 
### reported in the paper, and users can select them 
### according to their own needs.


library(logitnorm)
library(penPHcure)


data_generation <- function(a1 = alpha1, a2 = alpha0, 
                            b = beta, sha = shape_parameter, 
                            l1 = censoring_time, 
                            l2 = scale_parameter, 
                            nn = data_size)
{
  x1 <- rnorm(nn, 0, 1) # Random covariates for the cure rate
  p0 <- invlogit(x1*a1 + a2) # Estimates of the cure rate
  y0 <- rbinom(nn, 1, p0) # Cured status
  x2 <- rnorm(nn, 0, 1) # Random covariates for the hazard
  e0 <- exp(x2*b) # Estimates of the hazard
  t0 <- rep(0, nn) # Start time

  ts <- c()
  ds <- c()
  
  u1 <- runif(nn, 0, 1) # Random variates for time generation
  for (i in 1:nn){
    t <- ifelse(y0[i] == 0, 1000000000,
                (-log(u1[i])/(l2*e0[i]))^(1/sha))
    while (t <= t0[i]) {
      t <- ifelse(y0[i] == 0, 1000000000,
                  (-log(u1[i])/(l2*e0[i]))^(1/sha))
    }
    if (t <= l1){
      ds <- c(ds, 1); ts <- c(ts, t)
    } else {
      ds <- c(ds, 0); ts <- c(ts, tc)
    }
  }
  
  dd <- data.frame(x1, x2, t0, ts, y0, ds)
  dd <- round(dd, 2)
  return(dd)
}


# Generate background data
backgournd_data <- data_generation(a1 = alpha1, a2 = alpha0, 
                                   b = beta, 
                                   sha = shape_parameter, 
                                   l1 = censoring_time, 
                                   l2 = scale_parameter, 
                                   nn = data_size)


# Fit the mixture cure model
mcm <- penPHcure(Surv(time = t0, time2 = ts, event = ds) ~ x2, 
                 cureform = ~ x1, data = backgournd_data)





