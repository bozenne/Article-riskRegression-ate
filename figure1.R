library(data.table)
library(lava)
library(prodlim)

## * Import
source("./simulation-settings.R") ## import model2, setup.A, setup.B
createFigure1 <- function(model, parms, N, timeinterest){
    require(data.table)
    require(lava)
    require(prodlim)
    
    distribution(model,censtime) <- coxWeibull.lvm(scale=parms["scale.cr"])
    ppp <- parms[!grepl("formula|scale",names(parms))]
    dat <- as.data.table(lava::sim(model,
                                   n = N,
                                   p = ppp))
    cens.percent <- dat[,mean(time<timeinterest & event==0)]
    E1.percent <- dat[,mean(time<timeinterest & event==1)]
    E2.percent <- dat[,mean(time<timeinterest & event==2)]
    treat.percent <- dat[,mean(A==1)]
    dat[,dummy:=1]
    dat[eventtime2<eventtime1,dummy:=2]
    F <- prodlim(Hist(eventtime1, dummy) ~ A,data=dat)

    dat[,dummy:=1]
    dat[eventtime2<eventtime1,dummy:=2]
    
    ## ** RCT: no effect on treatment assignment
    parms0 <- parms
    parms0[grepl("alpha",names(parms0))] <- 0
    dat.RCT <- as.data.table(lava::sim(model,n = N,p = parms0))
    cens.percent.RCT <- dat.RCT[,mean(time<timeinterest & event==0)]
    E1.percent.RCT <- dat.RCT[,mean(time<timeinterest & event==1)]
    E2.percent.RCT <- dat.RCT[,mean(time<timeinterest & event==2)]
    treat.percent.RCT <- dat.RCT[,mean(A==1)]
    dat.RCT[,dummy:=1]
    dat.RCT[eventtime2<eventtime1,dummy:=2]
    F.RCT <- prodlim(Hist(eventtime1, dummy) ~ A,data=dat.RCT)

    ## **  plot the curves
    plot(F,confint=FALSE,xlim=c(0,timeinterest),legend=FALSE,atrisk=FALSE,ylab="Absolute risk of cause 1",col=c("gray66","black"))
    plot(F.RCT,add=TRUE,lty=2,confint=FALSE,xlim=c(0,timeinterest),col=c("gray66","black"))
    legend(col=c("gray66","black"),bty="n",x="topleft",lwd=3,lty=c(1),legend=c("Untreated","Treated"),title="Non-randomized")
    legend(col=c("gray66","black"),bty="n",x="topright",lwd=3,lty=c(2),legend=c("Untreated","Treated"),title="Randomized")
}

## create figure
pdf("figures/figure1-sim-setting.pdf", width = 10, height = 7)
par(mfrow=c(1,2))
set.seed(9)
do.call("createFigure1",list(model = model2, parms = setup.A$parameters, N = 10000, timeinterest = 8))
mtext("A",cex=2,at=-2,xpd=NA,line=1.5)
set.seed(10)
do.call("createFigure1",list(model = model2, parms = setup.B$parameters, N = 10000, timeinterest = 8))
mtext("B",cex=2,at=-2,xpd=NA,line=1.5)
dev.off()
