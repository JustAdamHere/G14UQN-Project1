library("deSolve")

# SEIRCD model.
fn<-function(time,state,params){
  with(as.list(c(state, params)),{
    dS=-beta*I*S/(S+E+I+R+D)
    dE=beta*I*S/(S+E+I+R+D)-0.25*E
    dI=0.25*E-gamma*I
    dR=(1-p)*gamma*I
    dC=0.25*E
    dD=p*gamma*I
    return(list(c(dS,dE,dI,dR,dC,dD)))
  })
}

# Data to compare with.
data <- read.csv("./data/Data_groupA.txt", sep="")

# Cumulative cases
Cummulative_cases=c(1)
number_cases_so_far=1
for(i in 1:199){
  Cummulative_cases=append(number_cases_so_far,Cummulative_cases)
  number_cases_so_far=number_cases_so_far+data["cases"][i,]
}
Cummulative_cases=rev(Cummulative_cases) #reverse it
Cummulative_cases=Cummulative_cases+rnorm(length(Cummulative_cases),0,1) #add errors

# Cumulative deaths
Cummulative_deaths=c(0)
number_deaths_so_far=0
for(i in 1:199){
  Cummulative_deaths=append(number_deaths_so_far,Cummulative_deaths)
  number_deaths_so_far=number_deaths_so_far+data["deaths"][i,]
}
Cummulative_deaths=rev(Cummulative_deaths) #reverse it
Cummulative_deaths=Cummulative_deaths+rnorm(length(Cummulative_deaths),0,1) #add errors

# Number of days, and population size.
day<-0:199
N<-2924
init <- c(S = N-1, E = 0,I = 1, R = 0, C = 0, D = 0)

#Estimate residual sum of squares for cases
RSS.SEIRC <- function(parameters)
{
  names(parameters) <- c("beta", "gamma","p")
  out <- rk4(y = init, times = day, func = fn, parms = parameters)
  fit <- out[ ,6]
  RSS <- sum((Cummulative_cases - fit)^2)
  return(RSS)
}

#Estimate residual sum of squares for deaths
RSS.SEIRD <- function(parameters)
{
  names(parameters) <- c("beta", "gamma","p")
  out <- rk4(y = init, times = day, func = fn, parms = parameters)
  fit <- out[ ,7]
  RSS <- sum((Cummulative_deaths - fit)^2)
  return(RSS)
}

# Upper and lower bounds for the guesses.
lower = c(0, 0, 0)
upper = c(2, 0.5, 0.5)

# Initial parameter guesses.
sigma <- 0.25
beta  <- 0.75
p     <- 0.25
gamma <- 0.1

# Optimise over deaths.
Opt2 <- optim(c(beta,gamma,p), RSS.SEIRD, method = "L-BFGS-B", lower = lower, upper = upper)
Opt_par2<-Opt2$par
names(Opt_par2)<-c("beta","gamma","p")

#Solve the ODEs for the 2nd optimal parameter estimates:
solution2<-as.data.frame(rk4(y = init, times = day, func = fn, parms = Opt_par2))
colnames(solution2)<-c("time","S","E","I","R","C","D")
#Plot the resulting C vs the real cummulative cases
plot(Cummulative_deaths)
lines(solution2$D,col="red")
#Plot the resukting D vs real cummulative deaths
plot(Cummulative_cases)
lines(solution2$C,col="red")

caseErrors  <- solution2$C-Cummulative_cases
deathErrors <- solution2$D-Cummulative_deaths

sampleVariance_cases  = sum(caseErrors) ^2/(length(caseErrors) -1)
sampleVariance_deaths = sum(deathErrors)^2/(length(deathErrors)-1)