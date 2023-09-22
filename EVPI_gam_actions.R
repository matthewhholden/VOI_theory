source(file="VOI_functions.R")

#a vector of actions
g = 15
n.actions = round( 2 ^ ((2:(2*g-2))/2) ) 
ma = length(n.actions)

#other params   
n.states = 100
n.sims = 100000
EVPI = matrix(NA, ma, n.sims)

print(paste(n.actions, " actions, ", n.states, " states"))

start_time <- Sys.time()
EVPI.sim = EVPI.gam.u.actions(m = n.actions, n.s = n.states, nsims = n.sims, alpha = 1, beta = 2)
end_time <- Sys.time()
print(end_time - start_time)

save(EVPI.sim, n.sims, n.actions, ma, n.states, file = paste("EVPI_ben_gap_gam_actions",ma,"_max_actions_",
                                                         n.states, "_states_", n.sims,"_sims.rdata", sep="") )

plot(log10(n.actions), apply(EVPI.sim[[1]]/EVPI.sim[[2]],1,median))

make.plot(EVPI.sim[[1]]/EVPI.sim[[2]], n.var = n.actions, nv=ma, logplot=TRUE, xa=TRUE)
