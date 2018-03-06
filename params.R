## Non density dependent parameters
params <- c("H" = 300, 
            "mu_N" =  1/60, 
            "mu_I" = ((1/10)-(1/60)), 
            "mu_W" = (1/(3.3*365)),
            "mu_H" = (1/(60*365)), 
            "f_N" = 0.1, 
            "m1" = 140, 
            "m2" = 3500, 
            "beta" = 1.5e-06,
            "sigma" = 1/40, 
            "K" = 10000, 
            "nu" = 0.1, 
            "omega" = 1, 
            "rho" = .43, 
            "pi_c" = 1.5e-08,
            "s" = 50, 
            "A" = 200)

## Density Dependent Parameters
dd_params <- c("k" = 0.2,
               "psi_gam" = 4,
               "psi_omeg" = .1,
               "f_gam" = .01,
               "del_gam" = .000001)

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
