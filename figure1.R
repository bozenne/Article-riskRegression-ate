## load 
source("simulation-functions.R") ## get createFigure1
source("simulation-models.R") ## get model2, setup.A, setup.B


pdf("figures/figure1-sim-setting.pdf")
par(mfrow=c(1,2))
set.seed(9)
do.call("createFigure1",list(model = model2, parms = setup.A, N = 10000, timeinterest = 8))
mtext("A",cex=2,at=-2,xpd=NA,line=1.5)
set.seed(10)
do.call("createFigure1",list(model = model2, parms = setup.B, N = 10000, timeinterest = 8))
mtext("B",cex=2,at=-2,xpd=NA,line=1.5)
dev.off()
