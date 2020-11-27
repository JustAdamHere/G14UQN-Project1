S <- function(t, y, parameters)
{
  with(as.list(c(y, parameters)),
  {
    N  <- S+E+I+R
    dS <- -beta*S*I/N
    
    return(list(dS))
  })
}

E <- function(t, y, parameters)
{
  with(as.list(c(y, parameters)),
  {
    N  <- S+E+I+R
    dE <- beta*S*I/N - sigma*E
    
    return(list(dE))
  })
}

I <- function(t, y, parameters)
{
  with(as.list(c(y, parameters)),
  {
    dI <- sigma*E - gamma*I
    
    return(list(dI))
  })
}

R <- function(t, y, parameters)
{
  with(as.list(c(y, parameters)),
  {
    dR <- (1-p)*gamma*I
    
    return(list(dR))
  })
}

C <- function(t, y, parameters)
{
  with(as.list(c(y, parameters)),
  {
    dC <- sigma*E
    
    return(list(dC))
  })
}

D <- function(t, y, parameters)
{
  with(as.list(c(y, parameters)),
  {
    dD <- p*gamma*I
    
    return(list(dD))
  })
}

SEIRCD <- function(t, y, parameters)
{
  with(as.list(c(y, parameters)),
  {
    N  <- S+E+I+R
    
    dS <- S(t, y, parameters)[[1]]
    dE <- E(t, y, parameters)[[1]]
    dI <- I(t, y, parameters)[[1]]
    dR <- R(t, y, parameters)[[1]]
    dC <- C(t, y, parameters)[[1]]
    dD <- D(t, y, parameters)[[1]]
    
    return(list(c(dS, dE, dI, dR, dC, dD)))
  })
}

CD <- function(t, y, parameters)
{
  with(as.list(c(y, parameters)),
  {
   dC <- C(t, y, parameters)[[1]]
   dD <- D(t, y, parameters)[[1]]
   
   return(list(c(dC, dD)))
  })
}

ODE_solver <- function(t_range, initial_y, parameters)
{
  library("deSolve")
  
  SIER_solution <- ode(y=initial_y, times=t_range, func=SEIRCD, parms=parameters, method="rk4")
  S = SIER_solution[, "S"]
  I = SIER_solution[, "I"]
  E = SIER_solution[, "E"]
  R = SIER_solution[, "R"]
  C = SIER_solution[, "C"]
  D = SIER_solution[, "D"]
  
  return(list(S, I, E, R, C, D))
}