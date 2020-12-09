setwd("C:/Users/adam/Blakey Cloud/Documents/Education/University of Nottingham/PhD/G14UQN/P1 SEIR Model/G14UQN-Project1/")

# ODE solver.
library("deSolve")
library("nlstools")
library("minpack.lm")
source("./ODE-solver.r")

# Initial population size.
N       <- 2924
#N       <- 10000
#N       <- 1000000
initial <- c(S=N-1, E=0, I=1, R=0, C=0, D=0)

# Initial guess at parameters.
parameters <- c(sigma=0.25, beta=0.5, p=0.25, gamma=0.1)
#parameters <- c(sigma=0.25, Opt_par2)

# Time range.
t_range = seq(1, 199, by=1)

# Reads in the data.
total_n    <- 199
training_n <- 60
data       <- read.csv("./data/Data_groupA.txt", sep="", nrows=training_n)
data_test  <- read.csv("./data/Data_groupA.txt", sep="", skip=training_n, nrows=total_n-training_n)

# Copies data to match nice names.
t           <- t_range
t_data      <- data[,1]
C_data      <- cumsum(data[, 2])
D_data      <- cumsum(data[, 3])
t_data_test <- data_test[,1]
C_data_test <- tail(C_data, n=1)+cumsum(data_test[, 3])
D_data_test <- tail(D_data, n=1)+cumsum(data_test[, 4])

# ODE solver for cases and deaths, taking only the unknown parameters as inputs.
ode_solve <- function(beta, p, gamma)
{
  parm <- c(sigma=parameters[["sigma"]], beta=beta, p=p, gamma=gamma)
  SEIRCD_solution <- rk4(initial, t_range, SEIRCD, parm)
  
  D <- SEIRCD_solution[, "D"]
  
  #return(D)
  return(approxfun(t_range, D)(t_data))
}

#fit=nls(deaths~ode_solve(beta, p, gamma), data=c(cumsum(data[1:10, c(1,2)]), cumsum(data[1:10, c(1,3)])), start=parameters[c(2, 3, 4)], control=c("nprint"=1, "maxiter"=50))
fit=nls(deaths~ode_solve(beta, p, gamma), data=c(cumsum(data[2]), cumsum(data[3])), start=parameters[c(2, 3, 4)], control=c("nprint"=1, "maxiter"=50))
updated_parameters <- c(sigma=parameters[["sigma"]], beta=coef(fit)[["beta"]], p=coef(fit)[["p"]], gamma=coef(fit)[["gamma"]])

# Solves equation with updated parameters.
SIER_solution <- rk4(y=initial, times=t_range, func=SEIRCD, parms=updated_parameters)
C_interp      <- approxfun(SIER_solution[, "time"], SIER_solution[, "C"])
D_interp      <- approxfun(SIER_solution[, "time"], SIER_solution[, "D"])

# Copies outputs to match nice names.
S_mod  <- SIER_solution[, "S"]
E_mod  <- SIER_solution[, "E"]
I_mod  <- SIER_solution[, "I"]
R_mod  <- SIER_solution[, "R"]
#C_mod  <- SIER_solution[, "C"]
C_mod  <- C_interp(t)
#D_mod  <- SIER_solution[, "D"]
D_mod  <- D_interp(t)

# Calculates the 95% confidence interval for the fit.
confidence_interval <- confint.default(fit, level=0.95)

## FROM HERE ON WE JUST PLOT THE BOUNDS ASSOCIATED WITH BETA
lower_parameters    <- c(sigma=parameters[["sigma"]], beta=confidence_interval[1, 1], p=coef(fit)[["p"]], gamma=coef(fit)[["gamma"]])
upper_parameters    <- c(sigma=parameters[["sigma"]], beta=confidence_interval[1, 2], p=coef(fit)[["p"]], gamma=coef(fit)[["gamma"]])
  
solution_beta_lower <- rk4(y=initial, times=t_range, func=SEIRCD, parms=lower_parameters)
solution_beta_upper <- rk4(y=initial, times=t_range, func=SEIRCD, parms=upper_parameters)

C_mod_beta_lower <- solution_beta_lower[, "C"]
D_mod_beta_lower <- solution_beta_lower[, "D"]
C_mod_beta_upper <- solution_beta_upper[, "C"]
D_mod_beta_upper <- solution_beta_upper[, "D"]

# Plots data and solution from parameters.
par(mfrow = c(1, 2))
plot(t_data, C_data, xlab="t", ylab="C", log="", xlim=c(0, max(t)), ylim=c(0, max(max(C_mod), max(C_data_test))))
points(t_data_test, C_data_test, col="green")
lines(t, C_mod, col="red")
lines(t, C_mod_beta_lower, col="blue")
lines(t, C_mod_beta_upper, col="blue")
plot(t_data, D_data, xlab="t", ylab="D", log="", xlim=c(0, max(t)), ylim=c(0, max(max(D_mod), max(D_data_test))))
points(t_data_test, D_data_test, col="green")
lines(t, D_mod, col="red")
lines(t, D_mod_beta_lower, col="blue")
lines(t, D_mod_beta_upper, col="blue")