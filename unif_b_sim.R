source(file="VOI_functions.R")

# #a vector of states
# m = 10^seq(-3,3,by=.5) 
# nm = length(m)
# nsims = 1000000
# 
# start_time <- Sys.time()
# dat = EVPI.unif(nsims = nsims, m=m)
# 
# end_time <- Sys.time()
# print(end_time - start_time)
# 
# save(m,nm,dat,nsims, file = paste("unif_", nm,"_b_vals_", nsims,"_sims.rdata", sep=""))
# 

load(file = "unif_13_b_vals_1e+06_sims.rdata")
EVPI.prop = dat$EVPI/dat$max.b.gap
plot(log10(m), apply(dat$EVPI/dat$max.b.gap, 1, median) )

filename =  "EVPI_unif_vary_b.pdf"
pdf(file = filename, width = 4, height = 4) # The height of the plot in inches
cl = 1.35
yl = c(0,0.65) #limits for y-axis of plotting windows
par(mfrow=c(1,1), oma = c(2,2,1,0), mar = c(2,2,1,1))
make.plot(EVPI.prop, n.var = m, nv = nm, logplot=TRUE, 
          xa=TRUE,ya=FALSE)
axis(side = 2, pos = log10(s[1]))
mtext(side = 1, "Upper utility bound", 
      outer = TRUE, padj = .5, adj = .5, cex = cl)
mtext(side = 2, "EVPI / Max benefit gap", 
      outer = TRUE, padj = -.5, adj = .5, cex = cl)
dev.off()