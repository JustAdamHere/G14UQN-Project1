# Problem details.
parameters <- c(sigma=2.99, beta=0.9, p=0.09, gamma=6.15)
initial    <- c(S=150, E=32, I=1, R=0, C=1, D=0)

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
t_range <- seq(0, 10, by=0.01)

# Solve the ODEs.
SIER_solution <- rk4(y=initial, times=t_range, func=SEIR, parms=parameters)

# Gives outer margins in lines of text for plotting.
par(oma = c(0, 0, 3, 0))

# Plots.
plot(SIER_solution, xlab="t")
mtext(outer=TRUE, side=3, paste("SIER Model after", tail(t_range, n=1), "days"), cex=1.5)