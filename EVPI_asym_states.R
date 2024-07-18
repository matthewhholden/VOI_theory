source(file="VOI_functions.R")

#a vector of states
g = 15
n.states = round( 2 ^ ((2:(2*g-2))/2) )  
ms = length(n.states)

n.actions = 100
n.sims = 100000
EVPI = matrix(NA, ms, n.sims)

print(paste(n.actions, " actions, ", n.states, " states"))

start_time <- Sys.time()
  EVPI = EVPI.asymp.states(m = n.states, n.a =n.actions, nsims = n.sims)
end_time <- Sys.time()
print(end_time - start_time)

save(EVPI,n.sims, n.states, n.actions, ms, file = paste("EVPI_asymp_states_",ms,"_max_states_",
                                             n.actions, "_actions_", n.sims,"_sims.rdata", sep=""))


plot(log10(n.states), apply(EVPI,1,median))
make.plot(EVPI, n.var = n.states, nv=ms, logplot=TRUE, xa=TRUE)
