

### All parameter values in the following code have been 
### reported in the paper, and users can select them 
### according to their own needs.


library(logitnorm)
library(MASS)
library(foreach)
library(doParallel)
library(penPHcure)


### RACC-EWMA
RACC <- function(a1 = alpha1, a2 = alpha0, b = beta, 
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

    yn <- (yy-pp)^2 - pp*(1-pp)
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


### Monte Carlo for UCL/LCL
# For UCL
numCores <- detectCores()
registerDoParallel(numCores)
up_rls <- foreach(o = 1:10000, .combine = "c", 
                  .packages = c("doParallel", "foreach", "MASS", 
                                "logitnorm", "penPHcure")) %dopar% 
  {RACC(a1 = alpha1, a2 = alpha0, b = beta, 
        sha = shape_parameter, l1 = censoring_time, 
        l2 = scale_parameter, cure_shift = 0, 
        hazard_shift = 0, cure_theta = 0, 
        hazard_theta = 0, backgournddata = backgournd_data, 
        smooth_parameter = lmd, ucl = h1, lcl = h2, 
        momitoring_direction = 1)}
stopImplicitCluster()

up_arl <- mean(up_rls) # Average run length for upward shitfs

# If up_arl is very close to the preset ARL0, 
# then h1 is the desired UCL.
# Otherwise, adjust h1 until the above requirements are met.


# For LCL
numCores <- detectCores()
registerDoParallel(numCores)
down_rls <- foreach(o = 1:10000, .combine = "c", 
                    .packages = c("doParallel", "foreach", "MASS", 
                                  "logitnorm", "penPHcure")) %dopar% 
  {RACC(a1 = alpha1, a2 = alpha0, b = beta, 
        sha = shape_parameter, l1 = censoring_time, 
        l2 = scale_parameter, cure_shift = 0, 
        hazard_shift = 0, cure_theta = 0, 
        hazard_theta = 0, backgournddata = backgournd_data, 
        smooth_parameter = lmd, ucl = h1, lcl = h2, 
        momitoring_direction = 2)}
stopImplicitCluster()

down_arl <- mean(down_rls) # Average run length for downward shitfs

# If down_arl is very close to the preset ARL0, 
# then h2 is the desired LCL.
# Otherwise, adjust h2 until the above requirements are met.

