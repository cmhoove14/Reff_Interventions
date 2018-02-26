library(ggplot2)
library(numDeriv)

source("/Users/lukestrgar/Documents/lab/params.R")

###### DD FUNCTION DEFINITIONS ######################################
#####################################################################

phi_Wk <- function(W,k=dd_params["k"]) {
  func <- function(x) {
    a <- ( W / (W + k) )
    b <- ((1-a)^(1+k))/(2*pi)
    return(( b*( 1-cos(x) ) / (( 1 + a*cos(x) )^(1+k)) ))
  }
  val = integrate(func, 0, 2*pi, subdivisions = 10000,
                  rel.tol = 1e-10, stop.on.error = FALSE)$value
  return(1-val)
}
##Host Immunosuppression
psi_W <- function(W,gamma=dd_params["psi_gam"],omega=dd_params["psi_omeg"]) {
  (1+(gamma*omega*W))/(1+(omega*W))
}
##Parasite Fecundity
f_Wgk <- function(W,gamma=dd_params["f_gam"],k=dd_params["k"]) {
  (1 + ((W*(1-(exp(-gamma))))/k))^(-k-1)
}
##Parasite Mortality
delta_Wgk <- function(W,gamma=dd_params["del_gam"],k=dd_params["k"]) {
  gamma*(k+1)*(W)/k
}

###### R0 + REFF ######################################
#####################################################################

r0 <- function() {
  N_eq <- s*A*(1-mu_N/f_N)
  num <- .5*beta*m1*nu*omega*H*N_eq*sigma*rho*omega*pi_c*m2
  denom <- (mu_N + sigma)*(mu_N + mu_I)*(mu_W + mu_H)
  as.numeric(num/denom)
}

reff<-function(W) {
  phi = phi_Wk(W)
  psi = psi_W(W)
  f = f_Wgk(W)
  delta = delta_Wgk(W)
  
  M <- .5*W*H*phi*f*m1*nu*omega
  x1 <- sigma/(mu_N+mu_I)
  x2 <- beta*M/(sigma+mu_N)
  S_eq <- (s*A/(1+x2+x1*x2))*(1-(mu_N+beta*M)/(f_N*(1+x2)))
  num <- rho*omega*pi_c*m2*x1*x2*psi*S_eq
  den <- (mu_W+mu_H+delta)*W
  reff <- as.numeric(num/den)
}


###### TOOLS FOR WORKING WITH REFF ######################################
#####################################################################

w_vals <- c(seq(from=0, to=1, by=0.001), seq(from=1, to=140, by=.01))
reff_prof <- function(w_range = w_vals) {
  reff_vals <- as.numeric()
  grad_vals <- as.numeric()
  
  for(i in 1:length(w_vals)) {
    reff_vals[i] <- reff(w_vals[i])
  }
  reff_vals
}

plot_reff <- function(w_range = w_vals) {
  reff_vals <- reff_prof(w_range)
  reff_data <- data.frame(reff_vals)
  w_data <- data.frame(w_range)
  
  ggplot(reff_data, aes(x=w_data, y=reff_data)) +
    geom_line(aes(y = reff_data)) +
    geom_hline(linetype="longdash", yintercept = 1) +
    labs(x=expression('W'), y=expression('R'["eff"]))
}