source("/Users/lukestrgar/Documents/lab/params.R")
source("/Users/lukestrgar/Documents/lab/model_funcs_tools.R")

sfx_mda <- function(x, params, t) {
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
  
  S = x['S']
  E = x['E']
  I = x['I']
  N = S + I + E
  Wt = x['Wt']
  Wu = x['Wu']
  W = cov*Wt+(1-cov)*Wu
  
  phi = phi_Wk(W)
  psi = psi_W(W)
  f = f_Wgk(W)
  delta = delta_Wgk(W)
  
  lambda = rho*omega*pi_c
  C = I*m2
  M = .5*W*H*phi*f*m1*nu*omega
  
  return(c(max(f_N*(1-(N/(s*A)))*(S+E), 0), #Snail birth
           mu_N * S, #Susceptible snail death
           beta*M*S, #Snail exposure
           mu_N*E, #Exposed snail dies
           sigma*E, #Exposed snail becomes infected
           (mu_N + mu_I)*I, #Infectd snail dies
           lambda*psi*C, #Infected snail produces adult worm
           (mu_W + mu_H + delta)*Wt, #Adult worm in treated population dies
           (mu_W + mu_H + delta)*Wu) #Adult worm in untreated population dies
  )
  
}

### Model
sfx_noMDA <- function(x, params, t) {
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
  
  S = x['S']
  E = x['E']
  I = x['I']
  N = S + I + E
  W = x['W']
  
  phi = phi_Wk(W)
  psi = psi_W(W)
  f = f_Wgk(W)
  delta = delta_Wgk(W)
  
  lambda = rho*omega*pi_c
  C = I*m2
  M = .5*W*H*phi*f*m1*nu*omega
  
  return(c(max(f_N*(1-(N/(s*A)))*(S+E), 0), #Snail birth
           mu_N * S, #Susceptible snail death
           beta*M*S, #Snail exposure
           mu_N*E, #Exposed snail dies
           sigma*E, #Exposed snail becomes infected
           (mu_N + mu_I)*I, #Infectd snail dies
           lambda*psi*C, #Infected snail produces adult worm
           (mu_W + mu_H + delta)*W) #Adult worm dies
  )
  
}