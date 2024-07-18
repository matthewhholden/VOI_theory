
source(file="VOI_functions.R")

filename =  "EVPI_asymp_multi_plot.pdf"

pdf(file = filename, width = 7, height = 7) # The height of the plot in inches

par(mfrow=c(2,2), oma = c(0,3,1,0), mar = c(4,2,1,.5))

### row 1
c.abc=1.5
x = log10(2.5); y =.46

load("EVPI_asymp_ 27 _max_actions_ 2 _states_ 1e+05 _sims.RData")
make.plot(EVPI, n.var = n.actions, nv=ma, 
          logplot=TRUE, xa=TRUE, ya=TRUE)
text(x=x, y=y, paste("a)  ", n.states, " states", sep=""), cex = c.abc, adj=0)
rm(EVPI, n.states, n.actions, ma)

load("EVPI_asymp_ 27 _max_actions_ 100 _states_ 1e+05 _sims.RData")
make.plot(EVPI, nv = ma, n.var = n.actions,
          xa = TRUE, ya = FALSE, logplot=TRUE)
text(x=x, y=y, paste("b)  ", n.states, " states", sep=""), cex = c.abc, adj=0)
rm(EVPI, n.states, n.actions, ma)


### row 2
load("EVPI_asymp_states_27_max_states_2_actions_1e+05_sims.RData")
make.plot(EVPI, nv = ms, n.var = n.states,
          xa = TRUE, ya = TRUE, logplot=TRUE)
text(x=x, y=y, paste("c)  ", n.actions, " actions", sep=""), cex = c.abc, adj=0)
rm(EVPI, n.states, n.actions, ms)

load("EVPI_asymp_states_27_max_states_100_actions_1e+05_sims.RData")
make.plot(EVPI, nv = ms, n.var = n.states,
          xa = TRUE, ya = FALSE, logplot=TRUE, leg=FALSE)
text(x=x, y=y, paste("d)  ", n.actions, " actions", sep=""), cex = c.abc, adj=0)
rm(EVPI, n.states, n.actions, ms)

mtext(side = 1, "Number of states (log scale)", outer = TRUE, padj = -1.5, adj = .5)
mtext(side = 1, "Number of actions (log scale)", outer = TRUE, padj = -30.5, adj = .5)
mtext(side = 2, " Expected Value of Perfect Information (EVPI) / max utility gap", outer = TRUE, padj = -0.5, adj = .5)

print("g")
dev.off()