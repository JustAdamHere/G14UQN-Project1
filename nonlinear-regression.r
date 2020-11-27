setwd("C:/Users/adam/Blakey Cloud/Documents/Education/University of Nottingham/PhD/G14UQN/P1 SEIR Model/G14UQN-Project1/")

# ODE solver.
library("deSolve")
source("./ODE-solver.r")

# Initial guess at parameters.
parameters <- c(sigma=1/3, beta=0.75, p=0.006, gamma=1/8)
initial    <- c(S=10000000, E=20000, I=1, R=0, C=1, D=0)

# Time range.
t_range = seq(0, 0.1, by=0.01)

# Reads in the data.
data <- read.csv("./data/Data_groupA.txt", sep="")

# SOME SORT OF LINEAR APPROXIMATION FOR ODE SOLVING BELOW?

#NLS <- function(ode_system, t_range, initial_conditions, initial_parameters, data)
#{
#  ode_solve <- function(parameters)
#  {
#    parm <- c(sigma=parameters[1], beta=parameters[2], p=parameters[3], gamma=parameters[4])
#    CD_solution <- ode(initial_conditions, t_range, ode_system, parm)
#    
#    C <- CD_solution[, "C"]
#    D <- CD_solution[, "D"]
#    
#    #return(list(C[1], D[1]))
#    return(C - data[[2]])
#  }
#  
#  parm <- c(initial_parameters[["sigma"]], initial_parameters[["beta"]], initial_parameters[["p"]], initial_parameters[["gamma"]])
#  
#  #optimal_parameters = optim(initial_parameters[c(2, 3, 4)], ode_solve, input1=initial_parameters[1])
#  optimal_parameters = optim(parm, ode_solve)
#  
#  return(list(beta=optimal_parameters$par[1], p=optimal_parameters$par[2], gamma=optimal_parameters$par[3]))
#}

#NLS(SEIRCD, t_range, initial, parameters, data)

ode_solve <- function(beta, p, gamma)
{
  parm <- c(sigma=parameters[["sigma"]], beta=beta, p=p, gamma=gamma)
  SEIRCD_solution <- rk4(initial, t_range, SEIRCD, parm)
  
  C <- SEIRCD_solution[, "C"]
  D <- SEIRCD_solution[, "D"]
  
  #return(list(C[1], D[1]))
  return(C)
}

#nls.fit=nls(C~approxfun(ode(initial, t_range, SEIRCD, parameters)[, 5]), data=data["cases"], start=parameters, control=list(maxiter=500))
#nls.fit=nls(C~ode(initial, t_range, SEIRCD, parameters)[, 5], data=list(data[, 2]), start=parameters, control=list(maxiter=500))
nls.fit=nls(cases~approxfun(ode_solve(beta, p, gamma)), data=data, start=parameters)
ode(initial["I"], t_range, D, parameters)