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
parameters <- c(sigma=0.25, beta=0.46129314, p=0.24604643, gamma=0.07520636)

# Time range.
t_range = seq(1, 199, by=1)

# Reads in the data.
total_n    <- 199
training_n <- 199
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

# Solves equation with updated parameters.
SIER_solution <- rk4(y=initial, times=t_range, func=SEIRCD, parms=parameters)
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

additionalPlotsC <- function(updatedParameter)
{
  updatedValues = c(0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2)*parameters[updatedParameter]
  
  for (i in 1:length(updatedValues))
  {
    updatedValue <- updatedValues[i]
    
    newParameters <- parameters
    newParameters[updatedParameter] <- updatedValue
    
    SIER_solution <- rk4(y=initial, times=t_range, func=SEIRCD, parms=newParameters)
    C_interp      <- approxfun(SIER_solution[, "time"], SIER_solution[, "C"])
    
    lines(t, C_interp(t), col=rgb(0, 1-i/length(updatedValues), i/length(updatedValues), 1))
  }
}

additionalPlotsD <- function(updatedParameter)
{
  updatedValues = c(0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2)*parameters[updatedParameter]
  
  for (i in 1:length(updatedValues))
  {
    updatedValue <- updatedValues[i]
    
    newParameters <- parameters
    newParameters[updatedParameter] <- updatedValue
    
    SIER_solution <- rk4(y=initial, times=t_range, func=SEIRCD, parms=newParameters)
    D_interp      <- approxfun(SIER_solution[, "time"], SIER_solution[, "D"])
    
    lines(t, D_interp(t), col=rgb(0, 1-i/length(updatedValues), i/length(updatedValues), 1))
  }
}

# Plots data and solution from parameters.
par(mfrow = c(1, 2))
plot(t_data, C_data, xlab="t", ylab="C", log="", xlim=c(0, max(t)), ylim=c(0, max(max(C_mod), max(C_data_test))))
points(t_data_test, C_data_test, col="green")
lines(t, C_mod, col="red")
additionalPlotsC("sigma")
mtext("Plots for varying sigma", side=3, line=-2, outer=TRUE, cex=1.5)

plot(t_data, D_data, xlab="t", ylab="D", log="", xlim=c(0, max(t)), ylim=c(0, max(max(D_mod), max(D_data_test))))
points(t_data_test, D_data_test, col="green")
lines(t, D_mod, col="red")
additionalPlotsD("sigma")