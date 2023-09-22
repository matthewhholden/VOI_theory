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
EVPI.sim = EVPI.gam.u.states(m = n.states, n.a =n.actions, nsims = n.sims, alpha=1, beta=2)
EVPI = EVPI.sim[[1]]/EVPI.sim[[2]]

end_time <- Sys.time()
print(end_time - start_time)

save(EVPI.sim,n.sims, n.states, n.actions, ms, file = paste("EVPI_ben_gap_gam_states_",ms,"_max_states_",
                                             n.actions, "_actions_", n.sims,"_sims.rdata", sep=""))


plot(log10(n.states), apply(EVPI,1,median))
make.plot(EVPI, n.var = n.states, nv=ms, logplot=TRUE, xa=TRUE)
