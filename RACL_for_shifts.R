

### All parameter values in the following code have been 
### reported in the paper, and users can select them 
### according to their own needs.


library(logitnorm)
library(MASS)
library(foreach)
library(doParallel)
library(penPHcure)


### RACL-EWMA
RACL <- function(a1 = alpha1, a2 = alpha0, b = beta, 
                 sha = shape_parameter, l1 = censoring_time, 
                 l2 = scale_parameter, cure_shift = 0, 
                 hazard_shift = 0, cure_theta = 0, 
                 hazard_theta = 0, backgournddata, 
                 smooth_parameter, ucl, lcl, 
                 momitoring_direction)
  # upwards when momitoring_direction = 1
  # downwards when momitoring_direction = 2
{
  # Fit the mixture cure model for the survival prediction
  mcm <- penPHcure(Surv(time = t0, time2 = ts, event = ds) ~ x2, 
                   cureform = ~ x1, data = backgournddata)
  
  z <- 0; i <- 1
  while (i > 0){
    t0 <- 0
    
    x1 <- rnorm(1, 0, 1)
    p0 <- invlogit(x1*a1 + a2) 
    # the cure rate without shifts
    p1 <- invlogit(x1*a1 + a2 + cure_shift + 
                     cure_theta*rnorm(1, 0, 1))
    # the cure rate with shifts
    y1 <- rbinom(1, 1, p1)
    
    x2 <- rnorm(1, 0, 1)
    e0 <- exp(x2*b)
    # the hazard without shifts
    e1 <- exp(x2*b + hazard_shift + 
                hazard_theta*rnorm(1, 0, 1))
    # the hazard with shifts
    
    u1 <- runif(1, 0, 1)
    ts <- ifelse(y1 == 0, 1000000000,
                 (-log(u1)/(l2*e1))^(1/sha))
    if (ts <= l1){
      ds <- 1
    } else {
      ts <- l1; ds <- 0
    }
    
    dd <- data.frame(x1, x2, t0, ts, y0, ds)
    rr <- predict(ccm, dd)
    pp <- rr$CURE[1]
    ss <- rr$SURV[1]
    yy <- ds + (1-ds) * pp*ss/(1-pp+pp*ss)
    # The estimation of the expectation of the cure status
    HH <- -log(ss) # The estimation of the cumulative_hazard
    
    yn <- (yy*(ds-HH))^2 - yy*HH
    z <- (1-smooth_parameter)*z + smooth_parameter*yn
    # Chart statistic
    
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


### ARL1 for the given shift scale
# Take monitoring the upward shift of random effects as an example
cure_shift <- cure_theta <- hazard_shift <- 0
hazard_theta_set <- c(0.1, 0.2, 0.5, 1, 2, 5, 8, 10)
# shfit scales
arl1 <- c()

numCores <- detectCores()
registerDoParallel(numCores)
for (j in 1:length(cure_theta_set)){
  rls1 <- foreach(o = 1:10000, .combine = "c", 
                  .packages = c("doParallel", "foreach", "MASS", 
                                "logitnorm", "penPHcure")) %dopar% 
    {RACL(a1 = alpha1, a2 = alpha0, b = beta, 
          sha = shape_parameter, l1 = censoring_time, 
          l2 = scale_parameter, cure_shift = 0, 
          hazard_shift = 0, cure_theta = hazard_theta_set[j], 
          hazard_theta = 0, backgournddata = backgournd_data, 
          smooth_parameter = lmd, ucl = h1, lcl = h2, 
          momitoring_direction = 1)}
  
  arl1 <- c(arl1, mean(rls1)) # Record results
}
stopImplicitCluster()





