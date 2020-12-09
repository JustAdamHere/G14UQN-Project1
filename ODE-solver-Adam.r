# Import correct library.
library("deSolve")

# Problem details (taken roughly from first example in Carcione et. al).
parameters <- c(sigma=1/3, beta=0.75, p=0.006, gamma=1/8)
initial    <- c(S=10000000, E=20000, I=1, R=0, C=0, D=0)

# Differential equation system for SEIR+CD.
SEIR <- function(t, initial, parameters)
{
  with(as.list(c(initial, parameters)),
  {
    N  <- S+E+I+R
   
    dS <- -beta*S*I/N
    dE <- beta*S*I/N - sigma*E
    dI <- sigma*E - gamma*I
    dR <- (1-p)*gamma*I
    dC <- sigma*E
    dD <- p*gamma*I
    
    list(c(dS, dE, dI, dR, dC, dD))
  })
}

# Time range for solving.
t <- seq(0, 100, by=0.01)

# Solve the ODEs.
SIER_solution <- rk4(y=initial, times=t, func=SEIR, parms=parameters)
S = SIER_solution[, "S"]
I = SIER_solution[, "I"]
E = SIER_solution[, "E"]
R = SIER_solution[, "R"]
C = SIER_solution[, "C"]
D = SIER_solution[, "D"]

# Gives layout for plotting.
par(mfrow = c(3, 2))

print(SIER_solution[,"S"])

# Plots.
plot(t, S, xlab="t")
plot(t, I, xlab="t")
plot(t, E, xlab="t")
plot(t, R, xlab="t")
plot(t, C, xlab="t", log="xy")
plot(t, D, xlab="t", log="xy")
mtext(outer=TRUE, side=3, paste("SIER Model after", tail(t_range, n=1), "days"), cex=1.5)

# Nonlinear regression.