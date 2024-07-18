
source(file="VOI_functions.R")

filename =  "EVPI_gamma_multi_plot.pdf"

pdf(file = filename, width = 7, height = 7) # The height of the plot in inches

par(mfrow=c(2,2), oma = c(0,3,1,0), mar = c(4,2,1,.5))

### row 1
cl = 1.35
yl = c(0,0.65) #limits for y-axis of plotting windows
x = log10(2.5); y =.92*yl[2] #location of a),b),c) labels

load("EVPI_ben_gap_gam_actions27_max_actions_2_states_1e+05_sims.rdata")
EVPI = EVPI.sim[[1]]/EVPI.sim[[2]]
make.plot(EVPI = EVPI, n.var = n.actions, nv=ma, 
          logplot=TRUE, xa=TRUE, ya=TRUE, ylimits = yl)
text(x=x, y=y, paste("a)  ", n.states, " states", sep=""), cex = cl, adj=0)
rm(EVPI.sim, EVPI, n.states, n.actions, ma)

load("EVPI_ben_gap_gam_actions27_max_actions_100_states_1e+05_sims.rdata")
EVPI = EVPI.sim[[1]]/EVPI.sim[[2]]
make.plot(EVPI = EVPI, n.var = n.actions, nv=ma, 
          logplot=TRUE, xa=TRUE, ya=FALSE, ylimits = yl)
text(x=x, y=y, paste("b)  ", n.states, " states", sep=""), cex = cl, adj=0)
rm(EVPI.sim, EVPI, n.states, n.actions, ma)


### row 2
load("EVPI_ben_gap_gam_states_27_max_states_2_actions_1e+05_sims.rdata")
EVPI = EVPI.sim[[1]]/EVPI.sim[[2]]
make.plot(EVPI, nv = ms, n.var = n.states,
          xa = TRUE, ya = TRUE, logplot=TRUE, ylimits = yl)
text(x=x, y=y, paste("c)  ", n.actions, " actions", sep=""), cex = cl, adj=0)
rm(EVPI.sim,EVPI, n.states, n.actions, ms)

load("EVPI_ben_gap_gam_states_27_max_states_100_actions_1e+05_sims.rdata")
EVPI = EVPI.sim[[1]]/EVPI.sim[[2]]
make.plot(EVPI, nv = ms, n.var = n.states,
          xa = TRUE, ya = FALSE, logplot=TRUE, leg=FALSE, ylimits = yl)
text(x=x, y=y, paste("d)  ", n.actions, " actions", sep=""), cex = cl, adj=0)

mtext(side = 1, "Number of states (log scale)", outer = TRUE, padj = -1.4, adj = .5, cex = cl)
mtext(side = 1, "Number of actions (log scale)", outer = TRUE, padj = -22.7, adj = .5, cex = cl)
mtext(side = 2, "    Expected value of perfect information (EVPI) / max utility gap", 
      outer = TRUE, padj = -1.3, adj = .5, cex = cl)

dev.off()
