# STAT 433 Final Project
# William Chen and Lucy Wu
# Spring 2018, Prof. Steele

#install.packages('normtest')
library(normtest)

## BASE METROPOLIS-HASTINGS ALGORITHM: ONE ITERATION
onestep_mh <- function(curr_state, prop_state, t_ij, t_ji, mean, sd) {
  
  # Runs one step of the Metropolis-Hastings algorithm.
  #
  # --Parameters--
  # curr_state : float
  #   Current state of the MC.
  # prop_state : float
  #   Proposed state of the MC.
  # t_ji : float
  #   Transition probability from j to i.
  # t_ij : float
  #   Transition probability from i to j.
  # mean : float
  #   Desired mean for the generated Normal distribution.
  # sd : float
  #   Desired stdev for the generated Normal distribution.
  #
  # -- Returns --
  # new_state : float
  #   New state of the MC. 
  
  # Compute elements of the acceptance function
  pi_j <- dnorm(prop_state, mean=mean, sd=sd)
  pi_i <- dnorm(curr_state, mean=mean, sd=sd)
  
  # Compute acceptance function value
  a <- (pi_j*t_ji)/(pi_i*t_ij)
  
  # Move to prop_state if a >= 1
  if (a >= 1) {
    new_state <- prop_state
  }
  
  # Move to prop_state with probability a 
  x <- runif(1, min=0, max=1)
  if (a >= x) {
    new_state <- prop_state
  }
  if (a < x) {
    new_state <- curr_state
  }
  
  return(new_state)
  # (Implied) Stay at curr_state otherwise.
}

## PROPOSAL CHAINS
# Uniform random walk proposal chain on [x-1, x+1]
unif_rw_1 <- function(curr_state) {
  prop_state <- runif(1, min=curr_state-1, max=curr_state+1)
  t_ij <- 1/2
  t_ji <- 1/2
  
  result <- c(prop_state, t_ij, t_ji)
  
  return(result)
}

# Uniform random walk proposal chain on [x-2, x+2]
unif_rw_2 <- function(curr_state) {
  prop_state <- runif(1, min=curr_state-2, max=curr_state+2)
  t_ij <- 1/4
  t_ji <- 1/4
  
  result <- c(prop_state, t_ij, t_ji)
  
  return(result)
}

# Uniform random walk proposal chain w/ symmetric Beta
symmetric_beta_22 <- function(curr_state) {
  prop_state_raw <- rbeta(1, shape1=2, shape2=2)
  prop_state <- curr_state + (prop_state_raw*2 - 1)# Change to sample on [x-1, x+1]
  
  t_ij <- dbeta(prop_state_raw, 2, 2) # Question: is this correct?
  t_ji <- dbeta((curr_state+1-prop_state)/2, 2, 2)
  
  result <- c(prop_state, t_ij, t_ji)
  
  return(result)
}

# Uniform random walk proposal chain w/ asymmetric Beta
asymmetric_beta_31 <- function(curr_state) {
  prop_state_raw <- rbeta(1, shape1=3, shape2=1)
  prop_state <- curr_state + (prop_state_raw*2 - 1)
  
  t_ij <- dbeta(prop_state_raw, 3, 1) # Question: is this correct?
  t_ji <- dbeta((curr_state+1-prop_state)/2, 3, 1)
  
  result <- c(prop_state, t_ij, t_ji)
  
  return(result)
}

asymmetric_beta_21 <- function(curr_state) {
  prop_state_raw <- rbeta(1, shape1=2, shape2=1)
  prop_state <- curr_state + (prop_state_raw*2 - 1)
  
  t_ij <- dbeta(prop_state_raw, 2, 1) # Question: is this correct?
  t_ji <- dbeta((curr_state-prop_state+1)/2, 2, 1)
  
  result <- c(prop_state, t_ij, t_ji)
  
  return(result)
}

# Proposal function

## FULL METROPOLIS-HASTINGS ALGORITHM

full_mh <- function(steps, prop_function, mean, sd) {
  
  curr_state <- 0 # Set some arbitrary initial state # QUESTION: IS this truly arbitrary?
  curr_state_list <- c()
  
  for(i in 1:steps) {
    
    # Get proposal state
    prop_state <- prop_function(curr_state)[1]
    t_ij <- prop_function(curr_state)[2]
    t_ji <- prop_function(curr_state)[3]
    
    # Run one iteration of MH
    new_state <- onestep_mh(curr_state, prop_state, t_ij, t_ji, mean, sd)
    curr_state <- new_state
    
    curr_state_list <- c(curr_state_list, curr_state)
  }
  
  hist(curr_state_list, prob=TRUE)
  curve(dnorm(x, mean=mean, sd=sd), col="darkblue", lwd=2, add=TRUE, yaxt="n")
  
  print(jb.norm.test(curr_state_list))
  
  return(curr_state_list)
}
