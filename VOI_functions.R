#### Functions for VOI analysis


library(LaplacesDemon)

EVPI.func = function(u, p){
  #u = utility matrix with entries a (actions, rows), s (states, columns)
  #p = vector of probabilities
  
  perf.info = sum( p * apply(u, 2, max) );
  no.info = max( apply( t(t(u) * p), 1, sum) );
  
  evpi = perf.info - no.info
  
  return(evpi)
  
}

EVPI.sim.states = function(m, nsims, n.a = 2, max.u=1, min.u=0){
  #m = vector of the number of states from 2 to max states
  #n.a = number of actions
  #nsims = number of simulations
  
  for(j in m){
    for(k in 1:nsims){
      p = as.vector(rdirichlet( 1, rep(1, j) ))
      u = matrix( runif(n.a*j), n.a, j )
      EVPI[j-1, k] = EVPI.func(u,p)
    }
  }
  
  return(EVPI)
}

EVPI.sim.actions = function(m, nsims, n.s = 2, max.u=1, min.u=0){
  #m = vector of the number of actions from 2 to max actions
  #n.a = number of actions
  #nsims = number of simulations
  
  for(j in m){
    for(k in 1:nsims){
      p = as.vector(rdirichlet( 1, rep(1, n.s) ))
      u = matrix( runif(n.s*j), j, n.s )
      EVPI[j-1, k] = EVPI.func(u,p)
    }
  }
  
  return(EVPI)
}

EVPI.sim.asymp.actions = function(m, nsims, n.s = 2, max.u=1, min.u=0){
  #m = vector of the number of actions 
  #n.s = number of states
  #nsims = number of simulations
  
  n=length(m)

  for(j in 1:n){
    print(j);
    for(k in 1:nsims){
      p = as.vector(rdirichlet( 1, rep(1, n.s) ))
      u = matrix( runif(n.s*m[j]), m[j], n.s )
      EVPI[j, k] = EVPI.func(u,p)
    }
  }
  
  return(EVPI)
}

EVPI.asymp.states = function(m, nsims, n.a = 2, max.u=1, min.u=0){
  #m = vector of the number of states from 2 to max states
  #n.a = number of actions
  #nsims = number of simulations
  n = length(m)
  
  for(j in 1:n){
    print(j)
    for(k in 1:nsims){
      p = as.vector(rdirichlet( 1, rep(1, m[j]) ))
      u = matrix( runif(n.a*m[j]), n.a, m[j] )
      EVPI[j, k] = EVPI.func(u,p)
    }
  }
  
  return(EVPI)
}

############### exponetially distributed VOI

EVPI.gam.u.actions = function(m, nsims, n.s = 2, alpha = 1, beta = 1){
  #m = vector of the number of actions 
  #n.s = number of states
  #nsims = number of simulations
  
  n=length(m)
  max.b.gap = matrix(NA, n, nsims)
  
  for(j in 1:n){
    print(j);
    for(k in 1:nsims){
      p = as.vector(rdirichlet( 1, rep(1, n.s) ))
      u = matrix( rgamma(n.s*m[j], shape = alpha, rate = beta ), 
                  m[j], n.s )
      EVPI[j, k] = EVPI.func(u,p)
      max.b.gap[j, k] = max( apply(u, 2, max) - apply(u, 2, min) )
    }
  }
  
  return(list(EVPI, max.b.gap))
}

EVPI.gam.u.states = function(m, nsims, n.a, alpha = 1, beta = 1){
  #m = vector of the number of states from 2 to max states
  #n.a = number of actions
  #nsims = number of simulations
  n = length(m)
  max.b.gap = matrix(NA, n, nsims)
  for(j in 1:n){
    print(j)
    for(k in 1:nsims){
      p = as.vector(rdirichlet( 1, rep(1, m[j]) ))
      u = matrix( rgamma(n.a*m[j], shape = alpha, rate = beta),
                  n.a, m[j] )
      EVPI[j, k] = EVPI.func(u,p)
      max.b.gap[j, k] = max( apply(u, 2, max) - apply(u, 2, min) )
    }
  }
  
  return(list(EVPI, max.b.gap))
}


EVPI.exp = function(nsims, l){
  #l = vector of lambdas (parameter of exponential distribution)
  #nsims = number of simulations
  n = length(l)
  max.b.gap = matrix(NA, n, nsims)
  EVPI = matrix(NA, n, nsims)
  p = matrix( runif(n*nsims), n, nsims )
  u = array( NA_real_, dim = c(n, nsims, 2, 2) )
  
  for(j in 1:n){
    print(j)
    
    u[j,,,] = array( rexp(4*nsims, rate = l[j]), dim=c(nsims, 2, 2 ) )
    
    for(k in 1:nsims){

      EVPI[j, k] = EVPI.func(u[j,k,,], c(p[j,k], 1 - p[j,k]) )
      max.b.gap[j, k] = max( apply(u[j,k,,], 2, max) - apply(u[j,k,,], 2, min) )
      
    }
  }
  
  return(list(EVPI=EVPI, max.b.gap=max.b.gap,u=u,p=p))
  
}

EVPI.unif = function(nsims, m){
  #m = vector of max of uniform (parameter of uniform distribution)
  #nsims = number of simulations
  n = length(m)
  max.b.gap = matrix(NA, n, nsims)
  EVPI = matrix(NA, n, nsims)
  p = matrix( runif(n*nsims), n, nsims )
  u = array( NA_real_, dim = c(n, nsims, 2, 2) )
  
  for(j in 1:n){
    print(j)
    
    u[j,,,] = array( runif(4*nsims, min = 0, max = m[j]), dim=c(nsims, 2, 2 ) )
    
    for(k in 1:nsims){
      
      EVPI[j, k] = EVPI.func(u[j,k,,], c(p[j,k], 1 - p[j,k]) )
      max.b.gap[j, k] = max( apply(u[j,k,,], 2, max) - apply(u[j,k,,], 2, min) )
      
    }
  }
  
  return(list(EVPI=EVPI, max.b.gap=max.b.gap,u=u,p=p))
  
}

EVPI.norm = function(nsims, s, mu=0){
  #m = vector of max of uniform (parameter of uniform distribution)
  #nsims = number of simulations
  n = length(s)
  max.b.gap = matrix(NA, n, nsims)
  EVPI = matrix(NA, n, nsims)
  p = matrix( runif(n*nsims), n, nsims )
  u = array( NA_real_, dim = c(n, nsims, 2, 2) )
  
  for(j in 1:n){
    print(j)
    
    u[j,,,] = array( rnorm(4*nsims, mean = mu, sd = s[j]), dim=c(nsims, 2, 2 ) )
    
    for(k in 1:nsims){
      
      EVPI[j, k] = EVPI.func(u[j,k,,], c(p[j,k], 1 - p[j,k]) )
      max.b.gap[j, k] = max( apply(u[j,k,,], 2, max) - apply(u[j,k,,], 2, min) )
      
    }
  }
  
  return(list(EVPI=EVPI, max.b.gap=max.b.gap,u=u,p=p))
  
}

#### Analytic VOI function
VOIanal = function( b1, b2, p){
  
  voi = NA
  p.star = b2 / (b1 + b2)
  
  if(p > p.star){
    
    voi = (1-p)*b2
    
  } else if(p < p.star){
    
    voi = b1*p
  }
  
  return( voi )
  
}


make.plot = function(EVPI, nv, n.var, 
                     leg = FALSE, labs = FALSE,
                     xa = FALSE, ya = FALSE, logplot=FALSE, ylimits = c(0,.5)){
  #plots EVPI against the variable being varied (actions or states)
  #nv = max number of the variable 
  #n.var = vector of numbers of states or actions
  
  m.EVPI = apply(EVPI, 1, median)
  quants = c(.0005, 0.025, .25, .75, 0.975, .9995)
  q = apply(EVPI, 1, quantile, probs = quants)

  print(m.EVPI)
  
  ## plot on regular scale  
  if(logplot == FALSE){
    plot(n.var, m.EVPI, 
         xlim = c(2,n.var[nv]), 
         ylim = ylimits, 
         type='l', lwd = 2,
         axes = FALSE, yaxs="i", xaxs="i",
         ylab="", xlab="")
    polygon(c(n.var,rev(n.var)),c(t(q[1,]),rev(q[6,])),
            col="gray87", border = NA)
    polygon(c(n.var,rev(n.var)),c(t(q[2,]),rev(q[5,])),
            col="gray75", border = NA)
    polygon(c(n.var,rev(n.var)),c(t(q[3,]),rev(q[4,])),
            col="gray60", border = NA)
    lines(n.var, m.EVPI, lwd = 2)
    box() 
    
    if(xa){ axis(side = 1, pos=0) }
    if(ya){ axis(side = 2, pos=2) }
    
    if(labs){
      fs = 1.2
      mtext(side=1, xlb, adj=0.5, padj=4, cex=fs)
      mtext(side=2, "Expected Value of Perfect Information", adj=0.5, padj=-4, cex=fs)
    }
    
  ## plot on log scale  
  } else if(logplot == TRUE){
    print('good')
    x = log10(n.var)
    if(n.var[1]==2){
      xticks = seq(1, log10(n.var[nv]), by = 1)
    } else {
      xticks = seq(log10(n.var[1]), log10(n.var[nv]), by = 1)
    }

    
    tick.labs = paste( 10^xticks )
    print(x)
    print(xticks)
    print(tick.labs)
    plot(x, m.EVPI, 
         xlim = c(min(x),max(x)), 
         ylim = ylimits, 
         type='l', lwd = 2,
         axes = FALSE, yaxs="i", xaxs="i",
         ylab="", xlab="")
    polygon(c(x,rev(x)),c(t(q[1,]),rev(q[6,])),
            col="gray87", border = NA)
    polygon(c(x,rev(x)),c(t(q[2,]),rev(q[5,])),
            col="gray75", border = NA)
    polygon(c(x,rev(x)),c(t(q[3,]),rev(q[4,])),
            col="gray60", border = NA)
    lines(x, m.EVPI, lwd = 2)
    box() 
    
    if(xa){ axis(at = xticks, labels = tick.labs, side = 1, pos=0) }
    if(ya){ axis(side = 2, pos=log10(2)) }
    
    if(labs){
      fs = 1.2
      mtext(side=1, xlb, adj=0.5, padj=4, cex=fs)
      mtext(side=2, "Expected Value of Perfect Information", adj=0.5, padj=-4, cex=fs)
    }
    
  }

  
  if(leg){
    
    legend(x=log10(30.5), y=.492, lwd=1, legend=c("median"),
           lty=c(1), bty="n",
           x.intersp=c(1.5)
    ) 
    
    legend(x=log10(31.7), y=.15, lwd=1,
           legend=c("50%", "95%", "99.9%"),
           fill=c("grey60", "grey75", "gray87"),
           lty =c(NA,NA,NA),
           bty="n",
           border=c("black", "black", "black"),
           x.intersp=rep(.1,3)
    ) 
  }

}