library(deSolve)
library(adaptivetau)
library(pracma)

source("/Users/lukestrgar/Documents/lab/params.R")
source("/Users/lukestrgar/Documents/lab/model_funcs_tools.R")
source("/Users/lukestrgar/Documents/lab/stochastic_models.R")
source("/Users/lukestrgar/Documents/lab/deterministic_models.R")

####### MDA + ALT INT STOCHASTIC SIMULATIONS #####################################
################################################################################################

### MDA Intervention Params
cov = 0.8
eff = 0.99
mda.years = c(2:21)

### Starting State
start = c(S = 5000, # susceptible snails
          E = 2000, # infected snails
          I = 500, # infected snails
          Wt = 72, # worms in treated population
          Wu = 72) # worms in untreated population

### Model Transitions
transitions = list(
  c(S = 1),             #New snail born
  c(S = -1),            #Susceptible snail dies
  c(S = -1, E = 1),     #Susceptible snail becomes exposed
  c(E = -1),            #Exposed snail dies
  c(E = -1, I = 1),     #Exposed snail becomes Infected
  c(I = -1),            #Infected snail dies
  c(Wt = 1, Wu = 1),    #Infected snail emits cercaria that produces an adult worm
  c(Wt = -1),           #Adult worm in the treated population dies
  c(Wu = -1))           #Adult worm in the untreated population dies


#function to simulate transmission over 61 years (1 year transmission spin up, 20 yrs MDA, 40 yrs recovery)
stoch.sim = function(init, cov, reduc, sim){
  
  ## Asign cov., frac to reduce nu by every iteration
  assign('cov', cov, envir = .GlobalEnv)
  frac = (1-reduc)^(1/20)
  
  # Run model to eq. 
  eq_traj = as.data.frame(ode(y=init,times=seq(0,200*365,30),
                              func=schisto_MDA_ODE,parms=params,method="ode23"))
  eq = eq_traj[length(seq(0,200*365,30)), c(2:6)]
  
  ## Save model traj
  w.traj.mdaplusalt = c(w.traj.mdaplusalt, (cov*eq$Wt)+(1-cov)*eq$Wu)
  assign('w.traj.mdaplusalt', w.traj.mdaplusalt, envir = .GlobalEnv)
  
  init1 = setNames(as.numeric(round(eq)), colnames(eq))
  
  set.seed(sim)
  
  fill[[1]] = ssa.adaptivetau(init1, transitions, 
                              sfx_mda, params, tf=365)    #simulate 1 year of transmission
  w.traj.mdaplusalt = c(w.traj.mdaplusalt, (cov*as.data.frame(fill[[1]])$Wt) + (1-cov)*as.data.frame(fill[[1]])$Wu)
  assign('w.traj.mdaplusalt', w.traj.mdaplusalt, envir = .GlobalEnv)
  
  for(m in 2:21){    #simulate 20 years of sanitation + MDA
    init = setNames(as.numeric(fill[[m-1]][dim(fill[[m-1]])[1],c(2:6)]), 
                    colnames(fill[[m-1]])[c(2:6)]) #reset initial states
    
    params["nu"] = frac*params["nu"] #apply sanitation

    fill[[m]] = ssa.adaptivetau(init, transitions, 
                                sfx_mda, params, tf=365) #stochastic sim for a year
    
    ## Save model traj
    w.traj.mdaplusalt = c(w.traj.mdaplusalt, (cov*as.data.frame(fill[[m]])$Wt) + (1-cov)*as.data.frame(fill[[m]])$Wu)
    assign('w.traj.mdaplusalt', w.traj.mdaplusalt, envir = .GlobalEnv)
    
    fill[[m]][,1] = fill[[m]][,1] + (365*(m-1)+(m-1))    #adjust time
  }
  
  for (k in 22:41) {
    init = setNames(as.numeric(fill[[k-1]][dim(fill[[k-1]])[1],c(2:6)]), 
                    colnames(fill[[k-1]])[c(2:6)]) #reset initial states
    
    init[4] = round(init[4]* (1-eff))  #NO MDA
    
    fill[[k]] = ssa.adaptivetau(init, transitions, 
                                sfx_mda, params, tf=365) #stochastic sim for a year
    
    ## Save model traj
    w.traj.mdaplusalt = c(w.traj.mdaplusalt, (cov*as.data.frame(fill[[k]])$Wt) + (1-cov)*as.data.frame(fill[[k]])$Wu)
    assign('w.traj.mdaplusalt', w.traj.mdaplusalt, envir = .GlobalEnv)
    
    fill[[k]][,1] = fill[[k]][,1] + (365*(k-1)+(k-1))    #adjust time
    
  }
  
  fill[[42]] = ssa.adaptivetau(init1, transitions, 
                              sfx_mda, params, tf=365)    #simulate 1 year of transmission
  
  ## Save model traj
  w.traj.mdaplusalt = c(w.traj.mdaplusalt, (cov*as.data.frame(fill[[42]])$Wt) + (1-cov)*as.data.frame(fill[[42]])$Wu)
  assign('w.traj.mdaplusalt', w.traj.mdaplusalt, envir = .GlobalEnv)
  
  for(f in 43:years){
    init = setNames(as.numeric(fill[[f-1]][dim(fill[[f-1]])[1],c(2:6)]), 
                    colnames(fill[[f-1]])[c(2:6)]) #reset initial states
    
    fill[[f]] = ssa.adaptivetau(init, transitions, 
                                sfx_mda, params, tf=365) #stochastic sim for a year
    w.traj.mdaplusalt = c(w.traj.mdaplusalt, (cov*as.data.frame(fill[[f]])$Wt) + (1-cov)*as.data.frame(fill[[f]])$Wu)
    assign('w.traj.mdaplusalt', w.traj.mdaplusalt, envir = .GlobalEnv)
    
    fill[[f]][,1] = fill[[f]][,1] + (365*(f-1)+(f-1))    #adjust time
  }
  
  matfin = do.call(rbind,fill)
  fill.test[c(1:nrow(matfin)), , sim] = cbind(matfin, Wm = cov*matfin[,5] + (1-cov)*matfin[,6])
  assign('fill.test', fill.test, envir = .GlobalEnv)
  
}


### Number of simulations and parameters to simulate over
cov.range = seq(.5,1,length=par.sims)
reduc.range = seq(.7,1,length=par.sims)
par.sims = 1
year.days = as.numeric()
for(i in 1:20){
  year.days[i] = 365*i + (i-1)
}
stoch.sims = 1
years = 62
fill = list()
fin.vals = c()
### vector for tracking W traj. NOTE -- this is a global variable and is setup to be written to for a single simulation
  ### (i.e. par.sims = stoch.sims = 1)
w.traj.mdaplusalt = c()

### Simulation Loop
for(i in 1:par.sims) {
  print("Iteration #")
  print(length(fin.vals)+1)
  
  fill.test = array(data = NA, dim = c(years*365*3, 7, stoch.sims))    #array to fill with all simulations
  pe1 = as.numeric()        #Vector for number of chains that go extinct
  
  #Run sims for parameter set 
  sapply(c(1:stoch.sims), stoch.sim, init = start, cov = cov.range[i], reduc = reduc.range[i], simplify = T)
  #Get probability of elimination (P(e)) as number of chains that lead to extinction out of all chains
  for(p in 1:stoch.sims) {
    print(fill.test[max(which(!is.na(fill.test[ , 7, p]))), , p][3:7])
    if(sum(fill.test[max(which(!is.na(fill.test[ , 7, p]))), , p][3:7]) == 0){ 
      pe1[p] = 1   #if no exposed, infected snails and no adult worms, elimination = 1
    } else {
      pe1[p] = 0   #else, elimination = 0
    }
  }
  
  fin.vals[length(fin.vals)+1] = sum(pe1) / stoch.sims
} ## End of Simulation Loop

