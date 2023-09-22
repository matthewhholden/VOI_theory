source(file="VOI_functions.R")
library(plotrix)
#######################################################################
#### VOI histogram analysis
#######################################################################

#number of simulations
n = 1000000

#max value of utility, for uniform random draws
m = 2
l = 1
mu = 1
sig = sqrt(4/12) #same variance as unif, note variance is 1 (1/l^2) for exp

#initialize variables
u.u = array(NA_real_, c(n, 2 ,2))
u.e = array(NA_real_, c(n, 2 ,2))
u.n = array(NA_real_, c(n, 2 ,2))
p = rep(NA_real_,n)

VOI.u = rep(NA_real_,n)
VOI.e = rep(NA_real_,n)
VOI.n = rep(NA_real_,n)

VOI.u.prop = rep(NA_real_,n)
VOI.e.prop = rep(NA_real_,n)
VOI.n.prop = rep(NA_real_,n)

for(i in 1:n){
  
  #draw u from specified distribution
    u.u[i, , ] = matrix( runif(4, min=0, max=m), 2, 2)
    u.e[i, , ] = matrix( rexp(4, rate=l), 2, 2)
    u.n[i, , ] = matrix( rnorm(4, mean=mu, sd = sig), 2, 2)
    
  #draw p
    p[i] = .5 #runif(1);#rbeta(1, shape1 = 5, shape2 = 5)
  
  #caclulate VOI
    VOI.u[i] = EVPI.func( u=u.u[i,,] , p=c(p[i], 1 - p[i] ) )
    VOI.e[i] = EVPI.func( u=u.e[i,,] , p=c(p[i], 1 - p[i] ) )
    VOI.n[i] = EVPI.func( u=u.n[i,,] , p=c(p[i], 1 - p[i] ) )
  
    VOI.u.prop[i] = VOI.u[i] / ( max(u.u[i,,]) - min(u.u[i,,]) )
    VOI.e.prop[i] = VOI.e[i] / ( max(u.e[i,,]) - min(u.e[i,,]) )
    VOI.n.prop[i] = VOI.n[i] / ( max(u.n[i,,]) - min(u.n[i,,]) )

}


#### VOI histograms
# hist(VOI, #breaks = seq(0,.5*m, by = 1), 
#      freq = FALSE, xlab="EVPI", main ="Histogram of EVPI")
# hist(VOI.prop, #breaks=seq(0,.5, by = .01), 
#      freq = TRUE, xlab="EVPI  /  max{ u }", 
#      main ="Histogram of EVPI as proportion of max utility")

#### non zero VOI histograms
VOI.e.non0 = VOI.e.prop[VOI.e.prop>0]  
VOI.u.non0 = VOI.u.prop[VOI.u.prop>0]  
VOI.n.non0 = VOI.n.prop[VOI.n.prop>0]  


br = seq(0,.5, by = .01) 

h.u =   hist(VOI.u.prop, freq = TRUE, breaks = br )
h.e =   hist(VOI.e.prop, freq = TRUE, breaks = br )
h.n =   hist(VOI.n.prop, freq = TRUE, breaks = br )

h.u.0 = hist(VOI.u.non0, freq = TRUE, breaks=br)
h.e.0 = hist(VOI.e.non0, freq = TRUE, breaks=br)
h.n.0 = hist(VOI.n.non0, freq = TRUE, breaks=br)

plot(h.u)

ylim = c(0, max(h.u.0$counts, h.n.0$counts, h.e.0$counts))
plot(h.u.0, col = rgb(0,0,1,.35), ylim = ylim)
plot(h.e.0, add=TRUE, col = rgb(1,0,0,.35) )
plot(h.n.0, add=TRUE, col = rgb(0,1,0,.35) )

par(mfrow=c(1,3))
plot(h.u.0, col = rgb(0,0,1,.35), ylim = ylim)
plot(h.e.0, col = rgb(1,0,0,.35), ylim = ylim)
plot(h.n.0, col = rgb(0,1,0,.35), ylim = ylim)


#figure distributions
pdf(file = 'EVPI_hists_exp_norm_unif_p_fixed_0_5.pdf', width = 7, height = 2.6)
par(mfrow = c(1,3), mar = c(1,1.2,1,1), oma= c(3,3,0,0))
lw = 1; ca=1.1; cl = 1.2;
plot(h.u.0, col = rgb(0,0,1,.35), ylim = ylim,
     xlab="", ylab="", main = "",
     xaxs='i', yaxs='i', 
     cex.axis = ca, lwd = lw)
text(x=.15, y = .99*ylim[2], 'a) Uniform utility', cex = cl)
plot(h.e.0, col = rgb(1,0,0,.35), ylim = ylim,
     xlab="", ylab="", main = "",
     xaxs='i', yaxs='i', yaxt = 'n', 
     cex.axis = ca, lwd = lw)
text(x=.19, y = .99*ylim[2], 'b) Exponential utility', cex = cl)
plot(h.n.0, col = rgb(0,1,0,.35), ylim = ylim,
     xlab="", ylab="", main = "",
     xaxs='i', yaxs='i', yaxt = 'n',
     cex.axis = ca, lwd = lw)
text(x=.14, y = .99*ylim[2], 'c) Normal utility', cex = cl)
mtext(side=1, "Expected value of perfect information (EVPI) / Max benefit gap", 
      outer = TRUE, padj=1.75, cex = cl)
mtext(side=2, "Frequency", outer = TRUE, padj=-1.5, cex = cl)
dev.off()


print(paste("prop zero, exp: ", sum(VOI.e.prop==0)/n))
print(paste("prop zero, norm: ", sum(VOI.n.prop==0)/n))
print(paste("prop zero, unif: ", sum(VOI.u.prop==0)/n))