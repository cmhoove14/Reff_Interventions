library(deSolve)
library(ggplot2)
library(scales)
library(rootSolve)
library(Deriv)
library(plotly)
library(numDeriv)
library(plot3D)
library(rgl)

source("/Users/lukestrgar/Documents/lab/params.R")
source("/Users/lukestrgar/Documents/lab/funcs.R")
source("/Users/lukestrgar/Documents/lab/ode.R")



############################################################################
####### PLOT MODEL TRAJECTORY ###############################################
############################################################################

initState <- c(S=5000, E=2000, I=500, W=70)
time = seq(0,365*10,1)
trajModel <- data.frame(ode(y=initState, times=time, func=helminth_ODE, 
                            parms=c(), method = "ode45"))

trajW <- data.frame(trajModel$W)
trajS <- data.frame(trajModel$S)
trajE <- data.frame(trajModel$E)
trajI <- data.frame(trajModel$I)
ODEtime_data <- data.frame(time)
ggplot(trajW, aes(x=ODEtime_data, y=trajW, col=Compartment)) +
  geom_line(aes(y = trajW, col="W"), size = 1.2) +
  #geom_line(aes(y = trajS, col = "S"), size = 1.2) +
  #geom_line(aes(y = trajE, col = "E"), size = 1.2) +
  geom_line(aes(y = trajI, col = "I"), size = 1.2) +
  labs(x = "Time", y = "W")




















  
  
#####################################################################
##### FEATURIZE REFF PROFILE  #############################################
########################################################################

#Vector of Parameter Values to Evaluate Critical Reff Features Over
param_vals <- c(.03, .1, .265, .5)
#Worm Burden Values
w_vals <- c(seq(from=0, to=1, by=0.1), seq(from=1, to=120, by=1))
#Empty Vectors For Reff Feature Values, Params of Reff Curve Fits, etc.
eeq_vals <- c()
bp_vals <- c()
wpeak_vals <- c()
rpeak_vals <- c()
m1_vals <- c()
m2_vals <- c()

for(i in 1:length(param_vals)) {
    #nu = params["nu"]*(1-param_vals[i])

    reff_opt<-function(W) {
      phi = phi_Wk(W,k=param_vals[i])
      psi = psi_W(W)
      f = f_Wgk(W,k=param_vals[i])
      delta = delta_Wgk(W,k=param_vals[i])
      
      M <- .5*W*H*phi*f*m1*nu*omega
      x1 <- sigma/(mu_N+mu_I)
      x2 <- beta*M/(sigma+mu_N)
      S_eq <- (s*A/(1+x2+x1*x2))*(1-(mu_N+beta*M)/(f_N*(1+x2)))
      num <- rho*omega*pi_c*m2*x1*x2*psi*S_eq
      den <- (mu_W+mu_H+delta)*W
      reff <- as.numeric(num/den)
    }
    
    rpeak_vals[i] <- optimize(reff_opt, interval=c(1e-05, 1e05), maximum=TRUE)$objective
    wpeak_vals[i] <- optimize(reff_opt, interval=c(1e-05, 1e05), maximum=TRUE)$maximum
    
    ### Reff Function w/1 subtracted from Reff to Compute Wbp, Weeq
    Reff_roots<-function(W) {
      phi = phi_Wk(W,k=param_vals[i])
      psi = psi_W(W)
      f = f_Wgk(W,k=param_vals[i])
      delta = delta_Wgk(W,k=param_vals[i])
      
      M <- .5*W*H*phi*f*m1*nu*omega
      x1 <- sigma/(mu_N+mu_I)
      x2 <- beta*M/(sigma+mu_N)
      S_eq <- (s*A/(1+x2+x1*x2))*(1-(mu_N+beta*M)/(f_N*(1+x2)))
      num <- rho*omega*pi_c*m2*x1*x2*psi*S_eq
      den <- (mu_W+mu_H+delta)*W
      reff <- as.numeric((num/den)-1)
    }
    
    
    result = tryCatch({
      bp_vals[i] <- uniroot(Reff_roots, c(1e-05,6.3))$root
      eeq_vals[i] <- uniroot(Reff_roots, c(6.3, 10000))$root
    }, error = function(e) {
      print(i)
      print("Error")
      break
    })
    
    m1_vals[i] <- 1/bp_vals[i]
    m2_vals[i] <- (rpeak_vals[i]-1)/(wpeak_vals[i]-bp_vals[i])
}


### Plot Curves w/ Various Reff Features Over Param Vals
param_data <- data.frame(param_vals)
bp_data <- data.frame(bp_vals)
eeq_data <- data.frame(eeq_vals)
rpeak_data <- data.frame(rpeak_vals)
wpeak_data <- data.frame(wpeak_vals)
front_slope_data <- data.frame(front_slope_vals)
back_slope_data <- data.frame(back_slope_vals)
ggplot(rpeak_data, aes(y=rpeak_data,x=param_data)) +
  geom_line(color="#DC0E0E") +
  #geom_line(aes(y=, x=param_data), color="blue") +
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
  #labs(x=expression(omega), y=expression('R'['Peak'])) +
  labs(x=expression(omega), y=expression('W')) +
  # scale_x_log10(breaks=c(.001, .01, .1, 1, 10, 100, 1000),
  #               labels=c(.001, .01, .1, 1, 10, 100, 1000)) +
  scale_y_log10(breaks=c(.01, .1, 1, 10, 100, 1000),
                   labels=c(.01, .1, 1, 10, 100, 1000))


