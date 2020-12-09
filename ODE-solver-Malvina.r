#Project task 1:
#Try to solve the ODE:
library("deSolve")
seir_ode<-function(t,Y,par){
  S<-Y[1]
  E<-Y[2]
  I<-Y[3]
  R<-Y[4]
  C<-Y[5]
  D<-Y[6]
  
  beta<-par[1]
  sigma<-par[2]
  gamma<-par[3]
  p<-par[4]
  
  dYdt<-vector(length=6)
  dYdt[1]=-beta*I*S/(S+E+I+R)
  dYdt[2]=beta*I*S/(S+E+I+R)-sigma*E
  dYdt[3]=sigma*E-gamma*I
  dYdt[4]=(1-p)*gamma*I
  dYdt[5]=sigma*E
  dYdt[6]=p*gamma*I
  
  return(list(dYdt))
}
#import dataset:
data <- read.csv("C:/Users/adam/Blakey Cloud/Documents/Education/University of Nottingham/PhD/G14UQN/P1 SEIR Model/G14UQN-Project1/data/Data_groupA.txt", sep="")
#set p=Number of deaths/cases:
p=sum(data["deaths"])/sum(data["cases"])
#Set sigma=0.25 as given:
sigma<-0.25 #given
#Set beta=the average number of people each infectious person spreads the disease to each day.
#The total number of infected:
total_number_of_cases=sum(data["cases"])
#Find the proportion of new infected people each day:
proportion=c()
for(i in 1:198){
    if(data["cases"][i+1,]!=0 & data["cases"][i,]!=0){
      proportion=append(proportion,c(data["cases"][i+1,]/data["cases"][i,]))
    }
}
#Find the average of that proportion:
average=sum(proportion)/length(proportion)
beta=average
#Set Gamma = the average time of infectiousness:
gamma<-0.1#Average time of infectiousness
init<-c(1,1,1,0,1,0)
t<-seq(0,nrow(data)-1)
par<-c(beta,sigma,gamma,p)
# Solve system using different methods:
sol1<-lsoda(init,t,seir_ode,par)
sol2<-as.data.frame(rk4(init,t,seir_ode,par))
sol3<-as.data.frame(euler(init,t,seir_ode,par))
S=sol2[,2]
E=sol2[,3]
I=sol2[,4]
R=sol2[,5]
C=sol2[,6]
D=sol2[,7]
sim.data<-as.data.frame(S)
sim.data["E"]<-E
sim.data["I"]<-I
sim.data["R"]<-R
sim.data["C"]<-C
sim.data["D"]<-D
sim.data["t"]<-0:198

#Fit data using nls:
ode1=function(t,e,params){
  s=params[1]
  dd=s*e
return(list(dd))
}
param.start=list(s=0.25)
nls.fit=nls(C~rk(E,t,ode1,c(s=s))[,2],data=sim.data,start=param.start,control = list(maxiter = 500))

##Do the same with the second one:
ode2=function(t,i,params){
  p=params[1]
  gamma=params[2]
  dd=gamma*p*i
  return(list(dd))
}
param.start=list(p=0.25,gamma=0.001)
nls.fit=nls(D~rk(I,t,ode2,c(p=p,gamma=gamma))[,2],data=sim.data,start=param.start,control = list(maxiter = 500))
#Make some calculations using data to find the commulative cases:
C=c()
number_cases_so_far=0
for(i in 1:199){
  C=append(number_cases_so_far,C)
  number_cases_so_far=number_cases_so_far+data["cases"][i,]
}
C=rev(C) #reverse it
D=c()
number_deaths_so_far=0
for(i in 1:199){
  D=append(number_deaths_so_far,D)
  number_deaths_so_far=number_deaths_so_far+data["deaths"][i,]
}
D=rev(D) #reverse it

