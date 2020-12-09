# Sets working directory and loads correct libraries.
setwd("C:/Users/adam/Blakey Cloud/Documents/Education/University of Nottingham/PhD/G14UQN/P1 SEIR Model/G14UQN-Project1/")
library("deSolve")
library("zoo")

# Script with the ODE functions.
source("./ODE-solver.r")

# Time step and time range.
dt      <- 1
t_range <- seq(1, 199, by=dt) 

# Initial conditions.
N       <- 2924
initial <- c(S=N-1, E=0, I=1, R=0, C=0, D=0)

# Parameter step sizes and ranges.
dBeta  <- 0.01
beta_range  <- seq(0.1, 0.8, dBeta)

# Interested time index.
t_index <- 20

# Storing the C values at a single time for all values in the parameter range.
C_values <- vector()
for (i in 1:length(beta_range))
{
  # Current beta.
  beta = beta_range[i]
  
  # Parameters we'll use.
  parameters <- c(sigma=0.25, beta=beta, p=0.25, gamma=1/9)
  
  # Calculates solution and stores C at correct t.
  SEIRCD_solution <- rk4(initial, t_range, SEIRCD, parameters)
  C_values[i] <- SEIRCD_solution[[t_index, 6]]
}

# Approximate derivative of C.
#   We note that the use of 'rollmean' is to get the midpoint of the x-axis
#   variables, where the approximation of the derivative holds best. It also
#   conveniently corrects the variable to the correct size for interpolation.
C_deriv <- approxfun(rollmean(beta_range, 2), diff(C_values)/dBeta)