source("/Users/lukestrgar/Documents/lab/params.R")
source("/Users/lukestrgar/Documents/lab/model_funcs_tools.R")

schisto_noMDA_ODE <- function(time, state, params) {
  
  H = params["H"] #Human population size (assumed constant)
  mu_N = params["mu_N"] #Natural snail death rate
  mu_I = params["mu_I"] #Increased death rate of infected snails
  mu_W = params["mu_W"] #Natural mature parasite death rate
  mu_H <- params["mu_H"] #Natural human death rate
  f_N <- params["f_N"] #Intrinsic snail fertility rate
  m1 <- params["m1"] #Rate of miracidia production per mated female worm
  m2 <- params["m2"] #Rate of cercariae production per infected intermediate host
  beta <- params["beta"] #Human to snail transmission parameter
  sigma <- params["sigma"] #Exposed to infected snail transition rate
  K <- params["K"] #Ecosystem snail carrying capacity
  nu <- params["nu"] #Fraction of miracida that make their way into waterway
  omega <- params["omega"] #Degree of crossover between snail habitats and human-water contact sites (site specific)
  rho <- params["rho"] #Per-capita human water contact rate
  pi_c <- params["pi_c"] #Probability of parasite establishment per human water contact for an individual cercariae 
  s <- params["s"] #Number of snails per Unit Area
  A <- params["A"] #Habitable units of area
  
  ### Intermediate Host States
  S <- state["S"]
  E <- state["E"]
  I <- state["I"]
  ### Total Intermediate Host Population
  N <- S+E+I 
  
  ### Human Infection State
  W = state['W']
  
  phi = phi_Wk(W)
  psi = psi_W(W)
  f = f_Wgk(W)
  delta = delta_Wgk(W)
  
  ### Time Dependent Number of Miracidia
  M <- .5*W*H*phi*f*m1*nu*omega
  ### Time Dependent Number of Cercaria
  C <- I*m2
  
  ### Snail to Human Transmission Parameter
  lambda <- rho*omega*pi_c
  
  ### ODE 
  dS<- f_N*(1-(N/(s*A)))*(S + E) - mu_N*S - beta*M*S
  dE<- beta*M*S - (sigma + mu_N)*E
  dI<- sigma*E - (mu_N + mu_I)*I
  dW <- lambda*psi*C - (mu_H + mu_W + delta)*W
  
  return(list(c(dS, dE, dI, dW)))
  
}

schisto_MDA_ODE <- function(time, state, params) {
  
  H = params["H"] #Human population size (assumed constant)
  mu_N = params["mu_N"] #Natural snail death rate
  mu_I = params["mu_I"] #Increased death rate of infected snails
  mu_W = params["mu_W"] #Natural mature parasite death rate
  mu_H <- params["mu_H"] #Natural human death rate
  f_N <- params["f_N"] #Intrinsic snail fertility rate
  m1 <- params["m1"] #Rate of miracidia production per mated female worm
  m2 <- params["m2"] #Rate of cercariae production per infected intermediate host
  beta <- params["beta"] #Human to snail transmission parameter
  sigma <- params["sigma"] #Exposed to infected snail transition rate
  K <- params["K"] #Ecosystem snail carrying capacity
  nu <- params["nu"] #Fraction of miracida that make their way into waterway
  omega <- params["omega"] #Degree of crossover between snail habitats and human-water contact sites (site specific)
  rho <- params["rho"] #Per-capita human water contact rate
  pi_c <- params["pi_c"] #Probability of parasite establishment per human water contact for an individual cercariae 
  s <- params["s"] #Number of snails per Unit Area
  A <- params["A"] #Habitable units of area
  
  ### Intermediate Host States
  S <- state["S"]
  E <- state["E"]
  I <- state["I"]
  ### Total Intermediate Host Population
  N <- S+E+I 
  
  ### Human Infection State
  Wt = state['Wt']
  Wu = state['Wu']
  W = cov*Wt+(1-cov)*Wu
  
  phi = phi_Wk(W)
  psi = psi_W(W)
  f = f_Wgk(W)
  delta = delta_Wgk(W)
  
  ### Time Dependent Number of Miracidia
  M <- .5*W*H*phi*f*m1*nu*omega
  ### Time Dependent Number of Cercaria
  C <- I*m2
  
  ### Snail to Human Transmission Parameter
  lambda <- rho*omega*pi_c
  
  ### ODE 
  dS<- f_N*(1-(N/(s*A)))*(S + E) - mu_N*S - beta*M*S
  dE<- beta*M*S - (sigma + mu_N)*E
  dI<- sigma*E - (mu_N + mu_I)*I
  dWt<- lambda*psi*C - (mu_H + mu_W + delta)*Wt
  dWu<- lambda*psi*C - (mu_H + mu_W + delta)*Wu
  
  return(list(c(dS,dE,dI,dWt,dWu)))
}