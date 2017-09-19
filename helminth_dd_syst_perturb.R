library(deSolve)
library(ggplot2)
library(scales)
library(rootSolve)
library(Deriv)
library(plotly)
library(numDeriv)


############################################################################
############################################################################
#### DENSITY DEPENDENT DISEASE SYSTEM EXPERIMENTS:  ########################
####  REFF'S SENSITIVITY TO DD'S AND SYSTEM PERTURBATIONS  #################
############################################################################
#####################  Author: Luke Strgar  ################################





############################################################################
####### GLOBAL PARAMETER DECLARATIONS #####################################
############################################################################
H <- 300 #Human population size (assumed constant)
mu_N <-  1/60 #Natural snail death rate
mu_I <- ((1/10)-(1/60)) #Increased death rate of infected snails
mu_W <- (1/(3.3*365)) #Natural mature parasite death rate
mu_H <- (1/(60*365)) #Natural human death rate
f_N <- 0.1 #Intrinsic snail fertility rate
m1 <- 140 #Rate of stage 1 larval production per mated female worm
m2 <- 3500 #Rate of stage 2 larval production per infected intermediate host
beta <- 5.6e-03 #Human to snail transmission parameter
sigma <- 1/40 #Exposed to infected snail transition rate
K <- 10000 #Ecosystem snail carrying capacity
nu <- 0.1 #Fraction of miracida that make their way into waterway
omega <- 1 #Degree of crossover between snail habitats and human-water contact sites (site specific)
rho <- .43 #Per-capita human water contact rate
pi <- 5e-07 #Probability of parasite establishment per human water contact for an individual cercariae 
c <- 50 #Number of Intermediate Host per Unit Area



  

  
#####################################################################
###### DD FUNCTION DEFINITIONS ######################################
#####################################################################

##Parasite Mating Function
phi_Wk <- function(W,k=0.08) {
  func <- function(x) {
    a <- ( W / (W + k) )
    c <- ((1-a)^(1+k))/(2*pi)
    return(( c*( 1-cos(x) ) / (( 1 + a*cos(x) )^(1+k)) ))
  }
  val = integrate(func, 0, 2*pi, subdivisions = 10000, 
                  rel.tol = 1e-10, stop.on.error = FALSE)$value
  return(1-val)
}
##Host Immunosuppression
psi_W <- function(W,gamma=2,omega=.1) {
  (1+(gamma*omega*W))/(1+(omega*W))
}
##Parasite Fecundity
f_Wgk <- function(W,gamma=0.000767,k=0.08) {
  (1 + ((W*(1-(exp(-gamma))))/k))^(-k-1)
}
##Parasite Mortality
delta_Wgk <- function(W,gamma=.00001,k=0.08) {
  gamma*(k+1)*(W)/k
}
##Intermediate Host Fertility
theta_NA <- function(N,A) {
  f_N*(1-(N/c*A))
}






############################################################################
####### DIFFERENTIAL EQUATION MODEL ########################################
############################################################################

helminth_ODE <- function(time, state, params) {
  
  ### Intermediate Host States
  S <- state["S"]
  E <- state["E"]
  I <- state["I"]
  ### Total Intermediate Host Population
  N <- S+E+I 
  
  ### Human Infection State
  W <- state["W"]

  ### Positive Density Dependent Function Values
  # Mating Function
  phi <- phi_Wk(W,k=.08)
  # Host Immunosuppression
  Is <- Is_W(W,gamma=5,omega=.1)
  
  ### Negative Density Dependent Function Values
  # Parasite Fecundity
  f <- f_Wgk(W,gamma=.047,k=.08)
  # Parasite Mortality
  delta <- delta_Wgk(W,gamma=.00001,k=.08)
  # Intermediate Host Fertility
  theta <- theta_NA(N,A)
  
  ### Dummy Assignments Used to Isolate Effects of Specific DD's
  #phi <- 1
  #f <- 1
  #delta <- 0
  #theta <- 1
  #Is <- 1
  
  ### Time Dependent Number of Miracidia
  M <- .5*W*H*phi*f*m1*nu*omega
  ### Time Dependent Number of Cercaria
  C <- I*m2
  
  ### Snail to Human Transmission Parameter
  lambda <- rho*omega*pi
  
  ### ODE 
  dS<- theta*(S + E) - mu_N*S - beta*M*S
  dE<- beta*M*S - (sigma + mu_N)*E
  dI<- sigma*E - (mu_N + mu_I)*I
  dW<- lambda*psi*C - (mu_H + mu_W + delta)*W

  return(list(c(dS,dE,dI,dW,dR)))
  
}







############################################################################
####### PLOT MODEL TRAJECTORY ###############################################
############################################################################

initState <- c(S=5000, E=0, I=0, W=0.2, R=0)
time = seq(0,365*10,1)
trajModel <- data.frame(ode(y=initState, times=time, func=helminth_ODE, 
                            parms=c(), method = "ode45"))

trajW <- data.frame(trajModel$W)
trajS <- data.frame(trajModel$S)
trajE <- data.frame(trajModel$E)
trajI <- data.frame(trajModel$I)
trajR <- data.frame(trajModel$R)
ODEtime_data <- data.frame(trajModel$time)
ggplot(trajW, aes(x=ODEtime_data, y=trajW, col=Compartment)) +
  #geom_line(aes(y = trajW, col="W"), size = 1.2) +
  #geom_line(aes(y = trajR, col="R"), size = 1.2) +
  geom_line(aes(y = trajS, col = "S"), size = 1.2) +
  #geom_line(aes(y = trajE, col = "E"), size = 1.2) +
  #geom_line(aes(y = trajI, col = "I"), size = 1.2) +
labs(x = "Time", y = "W")








#####################################################################
##### PLOT REFF PROFILE  #############################################
########################################################################

### Compute Values of Reff(W)
reff<-function(W) {
  ### Density Dependent Function Values
  phi <- phi_Wk(W,k=.08)
  psi <- Is_W(W,gamma=5,omega=.1)
  f <- f_Wgk(W,gamma=.047,k=.08)
  delta <- delta_Wgk(W,gamma=.00001,k=.08)
  theta <- theta_NA(N,A)
  
  ### Dummy Assignments Used to Isolate Effects of Specific DD's
  #phi <- 1
  #f <- 1
  #delta <- 0
  #theta <- 1
  #psi <- 1
  
  x1 <- sigma/(mu_N+mu_I)
  x2 <- beta*M/(sigma+mu_N)
  M <- .5*W*H*phi*f*m1*nu*omega
  S_eq <- (c*A/(1+x2+x1*x2))*(1-(mu_N+beta*M)/(f_N*(1+x2)))
  num <- rho*omega*pi*m2*x1*x2*psi*S_eq
  den <- (mu_W+mu_H+delta)*W
  reff <- as.numeric(num/den)
}


### Vector of W Values to Compute Reff, Reff', and Reff'' Over
#Worm Burden Values
w_vals <- c(seq(from=0, to=10, by = 0.01))
reff_vals <- as.numeric()
grad_vals <- as.numeric()
hess_vals <- as.numeric()
for(i in 1:length(w_vals)) {
  grad_vals[i] <- grad(reff, w_vals[i])
  hess_vals[i] <- hessian(reff, w_vals[i])
  reff_vals[i] <- reff(w_vals[i])
}


### Plot Single Reff / Reff Derivative Profile
reff_data <- data.frame(reff_vals)
grad_data <- data.frame(grad_vals)
w_data <- data.frame(w_vals)
ggplot(grad_data, aes(x=w_data, y=grad_data)) +
  #geom_line(aes(y = reff_data)) +
  #geom_line(aes(y = grad_data)) +
  #geom_line(aes(y = hess_data)) +
  #geom_hline(linetype="longdash", yintercept = 1) +
  #geom_hline(linetype="longdash", yintercept = 0) +
  labs(x=expression('W'), y=expression('R'["eff"]))


### Plot Multiple Reff / Derivative Profiles
w_data <- data.frame(w_vals)
ggplot(data1, aes(x=w_data, y=data1)) +
  #geom_line(aes(y = data1, x=w_data), color="orange", size = .6) +
  #geom_line(aes(y = data2, x=w_data), color="purple", size = .6) +
  #geom_line(aes(y = data3, x=w_data), color="blue", size = .6) +
  #geom_line(aes(y = data4, x=w_data), color="red", size = .6) +
  #geom_line(aes(y = data5, x=w_data), color="green", size = .6) +
  #geom_line(aes(y = data6, x=w_data), color="#FF33EC", size = .6) +
  #geom_line(aes(y = data7, x=w_data), color="brown", size = .6) +
  #geom_line(aes(y = data8, x=w_data), color="cyan", size = .6) +
  geom_hline(linetype="longdash", yintercept = 1)
  #geom_hline(linetype="longdash", yintercept = 0) +
  labs(x=expression('W'), y=expression('R'["eff"]))





  
  
#####################################################################
##### FEATURIZE REFF PROFILE  #############################################
########################################################################

#Vector of Parameter Values to Evaluate Critical Reff Features Over
param_vals <- c(seq(.001, .01, .0001), seq(.01,.1,.001))
#Worm Burden Values
w_vals <- c(seq(from=0, to=5, by = 0.1))
#Empty Vectors For Reff Feature Values, Params of Reff Curve Fits, etc.
eeq_vals <- c()
bp_vals <- c()
Rpeak_vals <- c()
Wpeak_vals <- c()
delta_bp_vals <- c()
delta_eeq_vals <- c()
deltaW_vals <- c()
reffavg_vals <- c()
k_vals <- c()
a_vals <- c()
b_vals <- c()

### Loop Over Parameter Values
for(i in 1:length(param_vals)) {
  
  #Reff / Reff' Profile for Given Param Val
  dreffdw_vals <- as.numeric()
  reff_vals <- as.numeric()
  for(j in 1:length(w_vals)) {
    dreffdw_vals[j] <- grad(reff, w_vals[j])
    reff_vals[j] <- reff(w_vals[j])
  }
  #Fit Reff / Reff' Curve to Exponential Pulse
  params <- nls(dreffdw_vals ~ k*(b*exp(-b*w_vals)-a*exp(-a*w_vals)), start = c(k=1,a=.2,b=15))$m$getAllPars()
  k_vals[i] <- params['k']
  a_vals[i] <- params['a']
  b_vals[i] <- params['b']

  ### Reff Function w/Parmeter Value Substituted. To compute Rpeak, Wpeak
  Reff_peak<-function(W) {
    phi <- phi_Wk(W,k=.08)
    psi <- Is_W(W,gamma=5,omega=.1)
    f <- f_Wgk(W,gamma=param_vals[i],k=.08)
    delta <- delta_Wgk(W,gamma=.00001,k=.08)
    theta <- theta_NA(N,A)
    #phi <- 1
    #f <- 1
    #delta <- 0
    #theta <- 1
    #psi <- 1
    
    x1 <- sigma/(mu_N+mu_I)
    x2 <- beta*M/(sigma+mu_N)
    M <- .5*W*H*phi*f*m1*nu*omega
    S_eq <- (c*A/(1+x2+x1*x2))*(1-(mu_N+beta*M)/(f_N*(1+x2)))
    num <- rho*omega*pi*m2*x1*x2*psi*S_eq
    den <- (mu_W+mu_H+delta)*W
    reff <- as.numeric(num/den)
  }

  ### Compute Rpeak, Wpeak
  Rpeak_vals[i] <- optimize(Reff_peak, interval=c(1e-05, 1e05), maximum=TRUE)$objective
  Wpeak_vals[i] <- optimize(Reff_peak, interval=c(1e-05, 1e05), maximum=TRUE)$maximum
  if (Rpeak_vals[i] <= 1) {
    print("Reff pushed to threshold value of 1")
    cat("  parameter value:", param_vals[i])
    cat("  W peak:", Wpeak_vals[i])
    cat("  R peak:", Rpeak_vals[i])
  }
  
  ### Reff Function w/Param Value Substituted, 1 subtracted from Reff. To Compute Wbp, Weeq
  Reff_roots<-function(W) {
    phi <- phi_Wk(W,k=.08)
    psi <- Is_W(W,gamma=5,omega=.1)
    f <- f_Wgk(W,gamma=param_vals[i],k=.08)
    delta <- delta_Wgk(W,gamma=.00001,k=.08)
    theta <- theta_NA(N,A)
    #phi <- 1
    #f <- 1
    #delta <- 0
    #theta <- 1
    #psi <- 1
    
    x1 <- sigma/(mu_N+mu_I)
    x2 <- beta*M/(sigma+mu_N)
    M <- .5*W*H*phi*f*m1*nu*omega
    S_eq <- (c*A/(1+x2+x1*x2))*(1-(mu_N+beta*M)/(f_N*(1+x2)))
    num <- rho*omega*pi*m2*x1*x2*psi*S_eq
    den <- (mu_W+mu_H+delta)*W
    reff <- as.numeric(num/den - 1)
  }
  
  ### Compute Wbp, Weeq
  result = tryCatch({
    bp_vals[i] <- uniroot(Reff_roots, c(1e-05,.150))$root
    eeq_vals[i] <- uniroot(Reff_roots, c(.150, 10000))$root
  }, error = function(e) {
    print(param_vals[i])
    break
  })
  
  ### Evaluate Various Distances Between Critical W Vals
  delta_bp_vals[i] <- Wpeak_vals[i] - bp_vals[i]
  delta_eeq_vals[i] <- eeq_vals[i] - Wpeak_vals[i]
  deltaW_vals[i] <- delta_bp_vals[i] + delta_eeq_vals[i]

  ### Evaluate Area Between Reff Curve and Reff = 1
  area <- integrate(Reff_roots, bp_vals[i], eeq_vals[i], rel.tol = 1e-10)$value
  reffavg_vals[i] <- area/(eeq_vals[i]-bp_vals[i])

}


### Plot Curves w/ Various Reff Features Over Param Vals
bp_data <- data.frame(bp_vals)
eeq_data <- data.frame(eeq_vals)
rpeak_data <- data.frame(Rpeak_vals)
wpeak_data <- data.frame(Wpeak_vals)
param_data <- data.frame(param_vals)
deltaW_data <- data.frame(deltaW_vals)
delta_bp_data <- data.frame(delta_bp_vals)
delta_eeq_data <- data.frame(delta_eeq_vals)
reffavg_data <- data.frame(reffavg_vals)
b_data <- data.frame(b_vals)
a_data <- data.frame(a_vals)
k_data <- data.frame(k_vals)
param_data <- data.frame(param_vals)
ggplot(b_data, aes(y=b_data,x=param_data)) +
  #geom_line(color="#DC0E0E") +
  #geom_line(aes(y=eeq_data,x=param_data), color="blue") +
  #geom_line(aes(y=rpeak_data, x=param_data), color="brown") +
  #geom_line(aes(y=deltaW_data, x=param_data), color="green") +
  #geom_line(aes(y=b_data/a_data, x=param_data), color="purple") +
  #geom_line(aes(y=delta_eeq_data/deltaW_data, x=param_data), color="orange") +
  #geom_line(aes(y=wpeak_data, x=param_data), color="#FF33EC") +
  #geom_line(aes(y=reffavg_data, x=param_data), color="black") +
  #labs(x=expression(gamma), y=expression(beta/alpha)) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15)) +
  xlim(0,.05)
  # scale_x_log10(breaks=c(.001, .01, .1, 1, 10, 100, 1000),
  #               labels=c(.001, .01, .1, 1, 10, 100, 1000)) +
  # scale_y_log10(breaks=c(.01, .1, 1, 10, 100, 1000),
  #                  labels=c(.01, .1, 1, 10, 100, 1000))


   
   
   
   
   
   
   
######################################################################
######## MISCELLANEOUS TESTING ##########################################
##########################################################################
   
b <- 6
a <- .3
c <- 2
c1 <- .1
c2 <- .5
c3 <- 1
c4 <- 2
c5 <- 3
c6 <- 4
c7 <- 5
c8 <- 6
b1 <- 0
b2 <- .5
b3 <- 1
b4 <- 1.5
b5 <- 2
b6 <- 2.5
b7 <- 3
b8 <- 6
a1 <- .5
a2 <- 1
a3 <- 1.5
a4 <- 2
a5 <- 2.5
a6 <- 3
a7 <- 3.5
a8 <- 4
func1 <- function(x) { c*(exp(-a*x) - exp(-b*x)) }
func2 <- function(x) { c*(b*exp(-b*x) - a2*exp(-a2*x)) }
func3 <- function(x) { c*(b*exp(-b*x) - a3*exp(-a3*x)) }
func4 <- function(x) { c*(b*exp(-b*x) - a4*exp(-a4*x)) }
func5 <- function(x) { c*(b*exp(-b*x) - a5*exp(-a5*x)) }
func6 <- function(x) { c*(b*exp(-b*x) - a6*exp(-a6*x)) }
func7 <- function(x) { c*(b*exp(-b*x) - a7*exp(-a7*x)) }
func8 <- function(x) { c*(b*exp(-b*x) - a8*exp(-a8*x)) }
vect1 <- Vectorize(func1)
vect2 <- Vectorize(func2)
vect3 <- Vectorize(func3)
vect4 <- Vectorize(func4)
vect5 <- Vectorize(func5)
vect6 <- Vectorize(func6)
vect7 <- Vectorize(func7)
vect8 <- Vectorize(func8)
plot(vect1, xlim=c(0,5), ylim=c(0,3))
plot(vect2, xlim=c(0,5), ylim=c(-1,1), add=TRUE)
plot(vect3, xlim=c(0,5), ylim=c(-.5,1), add=TRUE)
plot(vect4, xlim=c(0,5), ylim=c(-1,1), add=TRUE)
plot(vect5, xlim=c(0,5), ylim=c(-1,5), add=TRUE)
plot(vect6, xlim=c(0,5), ylim=c(-1,5), add=TRUE)
plot(vect7, xlim=c(0,5), ylim=c(-1,5), add=TRUE)
plot(vect8, xlim=c(0,5), ylim=c(-1,5), add=TRUE)


# reff_params <- nls(reff_vals ~ k*(exp(-a*w_vals)-exp(-b*w_vals)), start = c(k=1,a=.2,b=12))$m$getAllPars()
# dreff_params <- nls(grad_vals ~ k*(b*exp(-b*w_vals)-a*exp(-a*w_vals)), start = c(k=1,a=.2,b=15))$m$getAllPars()
