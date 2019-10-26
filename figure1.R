source("./simulation-models.R")
source("./simulation-functions.R")
setup.A <- c(BETA=0,
             OMEGA=0,
             GAMMA=0,
             alpha=rep(c(0,-1,-1,0,1,-1),2) * log(2),
             q.alpha=c(1,0,-1,0,1,-1) * log(2),
             beta=rep(c(1,-1,0,0,1,-1),2) * log(2),
             q.beta=c(1,-1,0,0,1,-1) * log(2),
             omega=rep(c(0,0,1,-1,1,-1),2) * log(2),
             q.omega=c(0,0,1,-1,1,-1) * log(2),
             gamma=rep(0,12),
             q.gamma=rep(0,12))
setup.B <- c(BETA=-2,
             OMEGA=0,
             GAMMA=0,
             alpha=rep(c(0,-1,-1,0,1,-1),2) * log(2),
             q.alpha=c(1,0,-1,0,1,-1) * log(2),
             beta=rep(c(1,-1,0,0,1,-1),2) * log(2),
             q.beta=c(1,-1,0,0,1,-1) * log(2),
             omega=rep(c(0,0,1,-1,1,-1),2) * log(2),
             q.omega=c(0,0,1,-1,1,-1) * log(2),
             gamma=rep(0,12),
             q.gamma=rep(0,12))

par(mfrow=c(1,2))
set.seed(9)
do.call("plot.model2",list(parms=setup.A,N=10000,timeinterest=8))
mtext("A",cex=2,at=-2,xpd=NA,line=1.5)
set.seed(10)
do.call("plot.model2",list(parms=setup.B,N=10000,timeinterest=8))
mtext("B",cex=2,at=-2,xpd=NA,line=1.5)
