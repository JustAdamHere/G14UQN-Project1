## Modified from https://cran.r-project.org/web/packages/deSolve/vignettes/deSolve.pdf.

# Problem details.
parameters <- c(a = -8/3, b = -10, c = 28)
IC         <- c(X = 1, Y = 1, Z = 1)

# Differential equation system.
Lorenz <- function(t, IC, parameters)
{
  with(as.list(c(IC, parameters)),
  {
    dX <- a*X + Y*Z
    dY <- b*(Y-Z)
    dZ <- -X*Y + c*Y - Z
  
    list(c(dX, dY, dZ))
  })
}

# Time range for solving.
t_range <- seq(0, 100, by=0.01)

# Solve the ODE.
solution <- ode(y=IC, times=t_range, func=Lorenz, parms=parameters)

# Gives outer margins in lines of text for plotting.
par(oma = c(0, 0, 3, 0))

# Plots.
plot(solution, xlab="t")
plot(solution[, "X"], solution[, "Z"], pch=".") # pch adds the plot rather than overwriting.
mtext(outer=TRUE, side=3, "Lorenz model", cex=1.5)