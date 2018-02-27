library(deSolve)
library(adaptivetau)

source("/Users/lukestrgar/Documents/lab/params.R")
source("/Users/lukestrgar/Documents/lab/model_funcs_tools.R")
source("/Users/lukestrgar/Documents/lab/stochastic_models.R")
source("/Users/lukestrgar/Documents/lab/deterministic_models.R")


####### ALTERNATIVE INTERVENTIONS STOCHASTIC SIMULATIONS #####################################
################################################################################################

### Starting State
start = c(S = 5000, # susceptible snails
          E = 2000, # infected snails
          I = 500, # infected snails
          W = 72) # worms in untreated population

### Model Transitions
transitions = list(
  c(S = 1),             #New snail born
  c(S = -1),            #Susceptible snail dies
  c(S = -1, E = 1),     #Susceptible snail becomes exposed
  c(E = -1),            #Exposed snail dies
  c(E = -1, I = 1),     #Exposed snail becomes Infected
  c(I = -1),            #Infected snail dies
  c(W = 1),  #Infected snail emits cercaria that produces an adult worm
  c(W = -1))           #Adult worm in the untreated population dies


#function to simulate transmission over 61 years (1 year transmission spin up, 20 yrs sanitation, 40 yrs recovery)
stoch.sim = function(init, reduc, sim){
  frac = (1-reduc)^(1/20)
  print("Simulation #")
  print(sim)
  eq_traj = as.data.frame(ode(y=init,times=seq(0,200*365,30),
                              func=schisto_noMDA_ODE,parms=params,method="ode23"))
  eq = eq_traj[length(seq(0,200*365,30)), c(2:5)]
  w.traj.nomda = c(w.traj.nomda, eq$W)
  assign('w.traj.nomda', w.traj.nomda, envir = .GlobalEnv)
  init1 = setNames(as.numeric(round(eq)), colnames(eq))
  
  set.seed(sim)
  
  fill[[1]] = ssa.adaptivetau(init1, transitions, sfx_noMDA, params, tf=365)    #simulate 1 year of transmission
  w.traj.nomda = c(w.traj.nomda, as.data.frame(fill[[1]])$W)
  assign('w.traj.nomda', w.traj.nomda, envir = .GlobalEnv)  

for(m in 2:21){    #simulate 20 years of sanitation
  init = setNames(as.numeric(fill[[m-1]][dim(fill[[m-1]])[1],c(2:5)]), 
                  colnames(fill[[m-1]])[c(2:5)]) #reset initial states
  
  params["nu"] = frac*params["nu"] #apply sanitation
  fill[[m]] = ssa.adaptivetau(init, transitions, 
                              sfx_noMDA, params, tf=365) #stochastic sim for a year
  w.traj.nomda = c(w.traj.nomda, as.data.frame(fill[[m]])$W)
  assign('w.traj.nomda', w.traj.nomda, envir = .GlobalEnv)
  
  fill[[m]][,1] = fill[[m]][,1] + (365*(m-1)+(m-1))    #adjust time
}

for(f in 22:years){
  init = setNames(as.numeric(fill[[f-1]][dim(fill[[f-1]])[1],c(2:5)]), 
                  colnames(fill[[f-1]])[c(2:5)]) #reset initial states
  
  #params["nu"] = frac*params["nu"] #no sanitation

  fill[[f]] = ssa.adaptivetau(init, transitions, 
                              sfx_noMDA, params, tf=365) #stochastic sim for a year
  w.traj.nomda = c(w.traj.nomda, as.data.frame(fill[[f]])$W)
  assign('w.traj.nomda', w.traj.nomda, envir = .GlobalEnv)
  
  fill[[f]][,1] = fill[[f]][,1] + (365*(f-1)+(f-1))    #adjust time
}

matfin = do.call(rbind,fill)
fill.test[c(1:nrow(matfin)), , sim] = matfin
assign('fill.test', fill.test, envir = .GlobalEnv)

}


### Number of simulations and parameters to simulate over
par.sims = 1
reduc.range = seq(.5,1,length=par.sims)
year.days = as.numeric()
for(i in 1:20){
  year.days[i] = 365*i + (i-1)
}
stoch.sims = 1
years = 61
fill = list()
fin.vals = c()
w.traj.nomda = c()

### Simulation Loop
for(i in 1:par.sims) {
  print("Iteration #")
  print(length(fin.vals)+1)
  
  fill.test = array(data = NA, dim = c(years*365*3, 5, stoch.sims))    #array to fill with all simulations
  pe1 = as.numeric()        #Vector for number of chains that go extinct
  
  #Run sims for parameter set  
  sapply(c(1:stoch.sims), stoch.sim, init = start, reduc = reduc.range[i], simplify = T)
  #Get probability of elimination (P(e)) as number of chains that lead to extinction out of all chains
  for(p in 1:stoch.sims) {
    print(fill.test[max(which(!is.na(fill.test[ , 5, p]))), , p][3:5])
    if(sum(fill.test[max(which(!is.na(fill.test[ , 5, p]))), , p][3:5]) == 0){ 
      pe1[p] = 1   #if no exposed, infected snails and no adult worms, elimination = 1
    } else {
      pe1[p] = 0   #else, elimination = 0
    }
  }
  
  fin.vals[length(fin.vals)+1] = sum(pe1) / stoch.sims
} ## End of Simulation Loop

