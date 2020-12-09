setwd("C:/Users/adam/Blakey Cloud/Documents/Education/University of Nottingham/PhD/G14UQN/P1 SEIR Model/G14UQN-Project1/")
library("deSolve")
SEIR <- function(t, initial, parameters)
{
  with(as.list(c(initial, parameters)),
       {
         
         dS <- -beta*S*I/(S+E+I+R+D)
         dE <- beta*S*I/(S+E+I+R+D)- 0.25*E
         dI <- 0.25*E - gamma*I
         dR <- (1-p)*gamma*I
         dC <- 0.25*E
         dD <- p*gamma*I
         
         list(c(dS, dE, dI, dR, dC, dD))
       })
}


data <- read.csv("./data/Data_groupA.txt", sep="")
#Cummulative cases
Cummulative_cases=c(1)
number_cases_so_far=1
for(i in 1:199){
  Cummulative_cases=append(number_cases_so_far,Cummulative_cases)
  number_cases_so_far=number_cases_so_far+data["cases"][i,]
}
Cummulative_cases=rev(Cummulative_cases) #reverse it

#Cummulative deaths
Cummulative_deaths=c(0)
number_deaths_so_far=0
for(i in 1:199){
  Cummulative_deaths=append(number_deaths_so_far,Cummulative_deaths)
  number_deaths_so_far=number_deaths_so_far+data["deaths"][i,]
}
Cummulative_deaths=rev(Cummulative_deaths) #reverse it

N=2924
init <- c(S = N-1, E=0,I = 1, R = 0,C=1,D=0)
day<-0:199
p<-sum(data["deaths"])/sum(data["cases"])
beta=1
gamma=0.1

mLL<-function(parameters){
  names(parameters) <- c("beta", "gamma","p")
  observations<-Cummulative_cases
  predictions <- rk4(init, times=day, func=SEIR, parms=parameters)[,6]
  # returning minus log-likelihood:
  -sum(dnorm(x = observations, mean = predictions, sd = 1, log = TRUE))
}

lower = c(0, 0, 0)
upper = c(2, 0.5, 0.5)

MLE <- optim(c(beta,gamma,p), mLL , method = "L-BFGS-B", lower = lower, upper = upper,hessian=TRUE)
MLE_par<-MLE$par
names(MLE_par)<-c("beta","gamma","p")

#check the fit:
MLE_predictions<-as.data.frame(rk4(y = init, times = day, func = SEIR, parms = MLE_par))
colnames(MLE_predictions)<-c("time","S","E","I","R","C","D")
#Plot the resulting C vs the real cummulative cases
plot(Cummulative_cases)
lines(MLE_predictions$C,col="red")

mLL2<-function(parameters){
  names(parameters) <- c("beta", "gamma","p")
  observations<-Cummulative_deaths
  predictions <- rk4(init, times=day, func=SEIR, parms=parameters)[,7]
  # returning minus log-likelihood:
  -sum(dnorm(x = observations, mean = predictions, sd = 1, log = TRUE))
}

lower = c(0, 0, 0)
upper = c(2, 0.5, 0.5)

MLE2 <- optim(c(beta,gamma,p), mLL2 , method = "L-BFGS-B", lower = lower, upper = upper,hessian=TRUE)
MLE2_par<-MLE2$par
names(MLE2_par)<-c("beta","gamma","p")

#check the fit:
MLE2_predictions<-as.data.frame(rk4(y = init, times = day, func = SEIR, parms = MLE2_par))
colnames(MLE2_predictions)<-c("time","S","E","I","R","C","D")
#Plot the results
par(mfrow=c(1,2))
plot(Cummulative_deaths)
lines(MLE2_predictions$D,col="red")
plot(Cummulative_cases)
lines(MLE2_predictions$C,col="red")


##Try to find CIs using the MLE estimates for the deaths likelihood:
hes2<-MLE2$hessian
sse<-sum((Cummulative_deaths-MLE2_predictions$D)^2)
sigma_m<-sse/197
v<-solve(sigma_m*hes2/200)
upper<-MLE2_par+qt(0.975,197)*sqrt(diag(v))
lower<-MLE2_par-qt(0.975,197)*sqrt(diag(v))
confidence_intervals<-as.data.frame(lower)
confidence_intervals$upper<-upper
confidence_intervals

#Try to find CIs for the predictions of cases: 
#we want to find max and min for each result when the parameters are in this intervals
beta_interval=confidence_intervals[1,]
gamma_interval=confidence_intervals[2,]
p_interval=confidence_intervals[3,]
minimal_predictions<-c()
maximal_predictions<-c()
beta=as.list(MLE2_par)$beta
gamma=as.list(MLE2_par)$gamma
p=as.list(MLE2_par)$p
for(i in 1:nrow(MLE2_predictions)){
  function_min<-function(parameters){
    names(parameters) <- c("beta", "gamma","p")
    prediction<-rk4(init, times=day, func=SEIR, parms=parameters)[,6]
    return(prediction[i])
  }
  function_max<-function(parameters){
    names(parameters) <- c("beta", "gamma","p")
    prediction<-rk4(init, times=day, func=SEIR, parms=parameters)[,6]
    return(-prediction[6])
  }
  lower = c(beta_interval$lower,gamma_interval$lower,p_interval$lower)
  upper = c(beta_interval$upper,gamma_interval$upper,p_interval$upper)
  minimal<-optim(c(beta,gamma,p), function_min , method = "L-BFGS-B", lower = lower, upper = upper)
  minimal_par<-minimal$par
  minimal_prediction<-as.data.frame(rk4(y = init, times = day, func = SEIR, parms = minimal_par))[,6][i]
  minimal_predictions<-append(minimal_prediction,minimal_predictions)
  maximal<-optim(c(beta,gamma,p), function_max , method = "L-BFGS-B", lower = lower, upper = upper)
  maximal_par<-maximal$par
  maximal_prediction<-as.data.frame(rk4(y = init, times = day, func = SEIR, parms = maximal_par))[,6][i]
  maximal_predictions<-append(maximal_prediction,maximal_predictions)
}
print(minimal_par)
print(maximal_par)
minimal_predictions<-rev(minimal_predictions)
maximal_predictions<-rev(maximal_predictions)

##Check on plot:
par(mfrow=c(1,1))
plot(Cummulative_cases)
lines(MLE2_predictions$C,col="red")
lines(minimal_predictions,col="red")
lines(maximal_predictions,col="red")

#Note: The confidence intervals for the parameters are so thin, so
#the results we get when fitting with different parameters in the intervals
#are exactly the same, so  we get that the CI for the predictions has
#the same lower and upper bounds = to the actual predictions. WEIRD! 