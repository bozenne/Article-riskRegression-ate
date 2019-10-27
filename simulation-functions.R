### simulation-functions.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Aug 16 2018 (08:47) 
## Version: 
## Last-Updated: okt 27 2019 (17:26) 
##           By: Brice Ozenne
##     Update #: 274
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(lava)
library(riskRegression)
library(data.table)
library(Publish)
library(parallel)
library(CausalGAM)
library(riskRegression)
library(survival)
library(data.table)

## * createFigure1
createFigure1 <- function(model,parms,N,timeinterest){
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
    cat("Observed data world:\n")
    cat(paste0("                     Treated: ",round(100*treat.percent,1),"%\n"))
    cat(paste0("Censored before timeinterest: ",round(100*cens.percent,1),"%\n"))
    cat(paste0(" Event 1 before timeinterest: ",round(100*E1.percent,1),"%\n"))
    cat(paste0(" Event 2 before timeinterest: ",round(100*E2.percent,1),"%\n"))
    cat(paste0("      Effect at timeinterest: ",round(100*diff(predictRisk(F,newdata=data.frame(A=0:1),times=timeinterest,cause=1)),1),"%\n"))
    cat("Atrisk at timeinterest:",paste0("(A=0)",paste(summary(F,cause=1,asMatrix=TRUE,times=timeinterest,newdata=data.frame(A=0:1))$table[,"n.risk"],collapse=" (A=1) ")))
    dat[,dummy:=1]
    dat[eventtime2<eventtime1,dummy:=2]
    ## RCT: no effect on treatment assignment
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
    cat("\nHypothetical RCT data world:\n")
    cat(paste0("                     Treated: ",round(100*treat.percent.RCT,1),"%\n"))
    cat(paste0("Censored before timeinterest: ",round(100*cens.percent.RCT,1),"%\n"))
    cat(paste0(" Event 1 before timeinterest: ",round(100*E1.percent.RCT,1),"%\n"))
    cat(paste0(" Event 2 before timeinterest: ",round(100*E2.percent.RCT,1),"%\n"))
    cat(paste0("      Effect at timeinterest: ",round(100*diff(predictRisk(F.RCT,newdata=data.frame(A=0:1),times=timeinterest,cause=1)),1),"%\n"))
### plot the curves
    plot(F,confint=FALSE,xlim=c(0,timeinterest),legend=FALSE,atrisk=FALSE,ylab="Absolute risk of cause 1",col=c("gray66","black"))
    plot(F.RCT,add=TRUE,lty=2,confint=FALSE,xlim=c(0,timeinterest),col=c("gray66","black"))
    legend(col=c("gray66","black"),bty="n",x="topleft",lwd=3,lty=c(1),legend=c("Untreated","Treated"),title="Non-randomized")
    legend(col=c("gray66","black"),bty="n",x="topright",lwd=3,lty=c(2),legend=c("Untreated","Treated"),title="Randomized")
}

## * rct.model2
# approximate true ATE value based on large uncensored data set
# where the treatment assigment is at random (alpha=0, q.alpha=0)
rct.model2 <- function(parms,N,timeinterest){
    parms0 <- parms
    parms0[grepl("alpha",names(parms0))] <- 0
    dat <- as.data.table(lava::sim(model2,n = N,p = parms0))
    ## uncensored data
    dat[,dummy:=1]
    dat[eventtime2<eventtime1,dummy:=2]
    out <- as.numeric(diff(predictRisk(prodlim(Hist(eventtime1, dummy) ~ A,data=dat),newdata=data.frame(A=0:1),times=timeinterest,cause=1)))
    out
}

## * run.simulation
## simulate from model and apply alternative estimators
run.simulation <- function(model,
                           setup,
                           B=1,
                           N,
                           timeinterest,
                           do.tmle.gam=FALSE,
                           do.tmle=FALSE,
                           do.ateRobust=TRUE,
                           do.known=FALSE,
                           coverage=TRUE,
                           causalgam=FALSE,
                           true.ate,
                           cores=1){
    parms <- setup$parameters
    forms <- setup$formula
    if (missing(true.ate)){
        true.ate <- setup$true.ate
    }
    if (is.null(true.ate)){
        true.ate <- rct.model2(parms,
                               N=50000,
                               timeinterest=timeinterest)
    }
    ## assign parameter values for lava
    assign("timeinterest",timeinterest,pos=-1)
    assign("N",N,pos=-1)
    for (v in names(parms)){assign(v,parms[[v]],pos=-1)}
    for (f in names(forms)){assign(f,forms[[f]],pos=-1)}
    run.one <- function(...,N,timeinterest,coverage=coverage){
        distribution(model,censtime) <- coxWeibull.lvm(scale=scale.cr)
        dat <- as.data.table(lava::sim(model,n = N,p = c(BETA=BETA,OMEGA=OMEGA,GAMMA=GAMMA,alpha1=alpha1,alpha2=alpha2,alpha3=alpha3,alpha4=alpha4,alpha5=alpha5,alpha6=alpha6,alpha7=alpha7,alpha8=alpha8,alpha9=alpha9,alpha10=alpha10,alpha11=alpha11,alpha12=alpha12,q.alpha1=q.alpha1,q.alpha2=q.alpha2,q.alpha3=q.alpha3,q.alpha4=q.alpha4,q.alpha5=q.alpha5,q.alpha6=q.alpha6,beta1=beta1,beta2=beta2,beta3=beta3,beta4=beta4,beta5=beta5,beta6=beta6,beta7=beta7,beta8=beta8,beta9=beta9,beta10=beta10,beta11=beta11,beta12=beta12,q.beta1=q.beta1,q.beta2=q.beta2,q.beta3=q.beta3,q.beta4=q.beta4,q.beta5=q.beta5,q.beta6=q.beta6,omega1=omega1,omega2=omega2,omega3=omega3,omega4=omega4,omega5=omega5,omega6=omega6,omega7=omega7,omega8=omega8,omega9=omega9,omega10=omega10,omega11=omega11,omega12=omega12,q.omega1=q.omega1,q.omega2=q.omega2,q.omega3=q.omega3,q.omega4=q.omega4,q.omega5=q.omega5,q.omega6=q.omega6,gamma1=gamma1,gamma2=gamma2,gamma3=gamma3,gamma4=gamma4,gamma5=gamma5,gamma6=gamma6,gamma7=gamma7,gamma8=gamma8,gamma9=gamma9,gamma10=gamma10,gamma11=gamma11,gamma12=gamma12,q.gamma1=q.gamma1,q.gamma2=q.gamma2,q.gamma3=q.gamma3,q.gamma4=q.gamma4,q.gamma5=q.gamma5,q.gamma6=q.gamma6)))
        while (max(dat$time)<timeinterest){
            dat <- as.data.table(lava::sim(model,n = N,p = c(BETA=BETA,OMEGA=OMEGA,GAMMA=GAMMA,alpha1=alpha1,alpha2=alpha2,alpha3=alpha3,alpha4=alpha4,alpha5=alpha5,alpha6=alpha6,alpha7=alpha7,alpha8=alpha8,alpha9=alpha9,alpha10=alpha10,alpha11=alpha11,alpha12=alpha12,q.alpha1=q.alpha1,q.alpha2=q.alpha2,q.alpha3=q.alpha3,q.alpha4=q.alpha4,q.alpha5=q.alpha5,q.alpha6=q.alpha6,beta1=beta1,beta2=beta2,beta3=beta3,beta4=beta4,beta5=beta5,beta6=beta6,beta7=beta7,beta8=beta8,beta9=beta9,beta10=beta10,beta11=beta11,beta12=beta12,q.beta1=q.beta1,q.beta2=q.beta2,q.beta3=q.beta3,q.beta4=q.beta4,q.beta5=q.beta5,q.beta6=q.beta6,omega1=omega1,omega2=omega2,omega3=omega3,omega4=omega4,omega5=omega5,omega6=omega6,omega7=omega7,omega8=omega8,omega9=omega9,omega10=omega10,omega11=omega11,omega12=omega12,q.omega1=q.omega1,q.omega2=q.omega2,q.omega3=q.omega3,q.omega4=q.omega4,q.omega5=q.omega5,q.omega6=q.omega6,gamma1=gamma1,gamma2=gamma2,gamma3=gamma3,gamma4=gamma4,gamma5=gamma5,gamma6=gamma6,gamma7=gamma7,gamma8=gamma8,gamma9=gamma9,gamma10=gamma10,gamma11=gamma11,gamma12=gamma12,q.gamma1=q.gamma1,q.gamma2=q.gamma2,q.gamma3=q.gamma3,q.gamma4=q.gamma4,q.gamma5=q.gamma5,q.gamma6=q.gamma6)))
        }
        if (causalgam==TRUE){
                                        # only censored data after timeinterest
            dat[,censbefore:=(event==0 & time<timeinterest)]
            dat[censbefore==TRUE,time:=pmin(eventtime2,eventtime1)]                
            dat[censbefore==TRUE,event:=as.numeric(1+1*(eventtime2<eventtime1))]
            dat[,censbefore:=NULL]
        }
        dat[,censtime:=NULL]
        if (do.ateRobust){
            dat[,AF:=as.factor(A)]
            if (do.known==TRUE){
                ateR <- ate(data = dat,estimator = c("Gformula","IPTW","AIPTW"),event = list(formula(paste0("Hist(time, event)~AF+",formula.event)),formula(paste0("Hist(time, event)~A+",formula.cr))),censor = formula(paste0("Surv(time, event==0)~AF+",formula.cens)),treatment = formula(paste0("AF~",formula.treatment)),known.nuisance=FALSE,times = timeinterest,product.limit = FALSE,verbose=0,cause = 1)$riskComparison
                ateR <- ateR[,grep("^diff",names(ateR),value=1),with=0]
                ateR.known <- ate(data = dat,estimator = c("Gformula","IPTW","AIPTW"),event = list(formula(paste0("Hist(time, event)~AF+",formula.event)),formula(paste0("Hist(time, event)~A+",formula.cr))),censor = formula(paste0("Surv(time, event==0)~AF+",formula.cens)),treatment = formula(paste0("AF~",formula.treatment)),known.nuisance=TRUE,times = timeinterest,product.limit = FALSE,verbose=0,cause = 1)$riskComparison
                ateR.known <- ateR.known[,grep("^diff",names(ateR.known),value=1),with=0]
            } else{
                ateR <- ate(data = dat,estimator = c("Gformula","IPTW","AIPTW"),event = list(formula(paste0("Hist(time, event)~AF+",formula.event)),formula(paste0("Hist(time, event)~A+",formula.cr))),censor = formula(paste0("Surv(time, event==0)~AF+",formula.cens)),treatment = formula(paste0("AF~",formula.treatment)),known.nuisance=FALSE,times = timeinterest,product.limit = FALSE,verbose=0,cause = 1)$riskComparison
                ateR <- ateR[,grep("^diff",names(ateR),value=1),with=0]
            }
        }
        km <- predictRisk(prodlim(Hist(time, event) ~ A,data=dat),newdata=data.frame(A=0:1),times=timeinterest,cause=1)
        ## print(km)
        ## if (any(is.na(km))) browser()
        ## gform <- ate(cfit,data=dat,treatment="A", times = timeinterest,se=FALSE,nuisance.iid = FALSE,cause = 1)
        if (do.tmle.gam==TRUE){
            SL.lib <- c("SL.glm","SL.gam")            
            tmle.gam <- survtmle(ftime=dat[["time"]],ftype=dat[["event"]],trt=dat[["A"]],adjustVars=data.frame(dat[,c("W1","W2","W3","W4","W5","W6","W7","W8","W9","W10","W11","W12"),with=FALSE]),SL.ftime=SL.lib,SL.trt=SL.lib,SL.ctime=SL.lib,method="mean",t0=timeinterest,ftypeOfInterest=1,returnModels=FALSE)
        }
        if (do.tmle){
            SL.lib <- c("SL.glm")            
            try(tmle <- survtmle(ftime=dat[["time"]],ftype=dat[["event"]],trt=dat[["A"]],adjustVars=data.frame(dat[,c("W1","W2","W3","W4","W5","W6","W7","W8","W9","W10","W11","W12"),with=FALSE]),glm.trt=formula.treatment,glm.ftime=paste0("trt+",formula.event),glm.ctime=paste0("trt+",formula.cens),method="mean",t0=timeinterest,ftypeOfInterest=1,returnModels=TRUE))
        }
        if (causalgam==TRUE){
            library(CausalGAM)
            dat[,Y:=eventtime1<timeinterest & eventtime1<eventtime2]
            ATE.out <- estimate.ATE(pscore.formula = formula(paste0("A~",formula.treatment)),pscore.family = binomial,outcome.formula.t = formula(paste0("Y~A+",formula.event)),outcome.formula.c = formula(paste0("Y~A+",formula.event)),outcome.family = binomial,treatment.var = "A",data=dat,divby0.action="t",divby0.tol=0.001,var.gam.plot=FALSE,nboot=0)
        }
        out <- c(KM.naive=diff(km))
        if (do.ateRobust) {
            out <- c(out,unlist(ateR[,grep("a$|W$",names(ateR),value=1),with=0]))
            if (do.known==TRUE){
                out <- c(out,unlist(ateR.robust[,grep("a$|W$",names(ateR),value=1),with=0]))
            }
            names(out) <- gsub("diff.","",names(out))
        }
        if (do.tmle) out <- c(out,TMLE=diff(tmle$est))
        if (causalgam==TRUE){
            out <- c(out,c(CausalGam.IPW=ATE.out$ATE.IPW.hat,CausalGam.AIPW=ATE.out$ATE.AIPW.hat,CausalGam.reg=ATE.out$ATE.reg.hat))
        }
        if (do.tmle.gam==TRUE) out <- c(out,TMLE.gam=diff(tmle.gam$est))
        if (coverage){
            ## if (do.ateRobust){
            ## se.rob <- c(ateR$ate.se[3,])
            ## }else se.rob <- NULL
            ## if (do.tmle){
            ## se.tmle <- c("TMLE"=sd(apply(tmle$ic,1,diff))/sqrt(N))
            ## }else se.tmle <- NULL
            ## se <- c(se.rob,se.tmle)
            ## lower <- out[-1]+qnorm(0.025)*se
            ## upper <- out[-1]+qnorm(0.975)*se
            se <- unlist(ateR[,grep(".se$",names(ateR),value=1),with=0])
            names(se) <- gsub("diff.","",names(se))
            if (do.known){
                se.known <- unlist(ateR[,grep(".se$",names(ateR),value=1),with=0])
                names(se.known) <- gsub("diff.","",names(se.known))
                lower.known <- unlist(ateR.known[,grep(".lower$",names(ateR.known),value=1),with=0])
                upper.known <- unlist(ateR.known[,grep(".upper$",names(ateR.known),value=1),with=0])
                cov.known <- lower.known<true.ate & upper.known>true.ate
                names(cov.known) <- gsub(".lower|diff.","",paste0("coverage.",names(cov.known)))
            }
            lower <- unlist(ateR[,grep(".lower$",names(ateR),value=1),with=0])
            upper <- unlist(ateR[,grep(".upper$",names(ateR),value=1),with=0])
            cov <- lower<true.ate & upper>true.ate
            names(cov) <- gsub(".lower|diff.","",paste0("coverage.",names(cov)))
            out <- c(out,true=true.ate,se,cov)
            if (do.known){
                out <- c(out,se.known=se.known,cov.known=cov.known)
            }
        }
        cens.percent <- dat[,mean(time<timeinterest & event==0)]
        out <- c(cens.percent=cens.percent,true.ate=true.ate,out)
        out
    }
    debug(run.one)
    browser()
    lava.options(messages=1)
    lava.options(messages=0)
    out <- sim(run.one,B,mc.cores=cores,timeinterest=timeinterest,N=N,coverage=coverage)
    ## })
    out
}

## * run.simulation1
## simulate from model and apply alternative estimators
run.simulation1 <- function(model,
                            setup,
                            B=1,
                            N,
                            timeinterest,
                            do.ateRobust=TRUE,
                            do.known=FALSE,
                            coverage=TRUE,
                            true.ate,
                            cores=1){
    parms <- setup$parameters
    forms <- setup$formula
    if (missing(true.ate)){
        true.ate <- setup$true.ate
    }
    if (is.null(true.ate)){
        true.ate <- rct.model2(parms,N=50000,timeinterest=timeinterest)
    }
    ## assign parameter values for lava
    assign("timeinterest",timeinterest,pos=-1)
    assign("N",N,pos=-1)
    for (v in names(parms)){assign(v,parms[[v]],pos=-1)}
    for (f in names(forms)){assign(f,forms[[f]],pos=-1)}
    run.one <- function(...,N,timeinterest,coverage=coverage){
        distribution(model,censtime) <- coxWeibull.lvm(scale=scale.cr)
        dat <- as.data.table(lava::sim(model,n = N,p = c(BETA=BETA,OMEGA=OMEGA,GAMMA=GAMMA,alpha1=alpha1,alpha2=alpha2,alpha3=alpha3,alpha4=alpha4,alpha5=alpha5,alpha6=alpha6,alpha7=alpha7,alpha8=alpha8,alpha9=alpha9,alpha10=alpha10,alpha11=alpha11,alpha12=alpha12,q.alpha1=q.alpha1,q.alpha2=q.alpha2,q.alpha3=q.alpha3,q.alpha4=q.alpha4,q.alpha5=q.alpha5,q.alpha6=q.alpha6,beta1=beta1,beta2=beta2,beta3=beta3,beta4=beta4,beta5=beta5,beta6=beta6,beta7=beta7,beta8=beta8,beta9=beta9,beta10=beta10,beta11=beta11,beta12=beta12,q.beta1=q.beta1,q.beta2=q.beta2,q.beta3=q.beta3,q.beta4=q.beta4,q.beta5=q.beta5,q.beta6=q.beta6,omega1=omega1,omega2=omega2,omega3=omega3,omega4=omega4,omega5=omega5,omega6=omega6,omega7=omega7,omega8=omega8,omega9=omega9,omega10=omega10,omega11=omega11,omega12=omega12,q.omega1=q.omega1,q.omega2=q.omega2,q.omega3=q.omega3,q.omega4=q.omega4,q.omega5=q.omega5,q.omega6=q.omega6,gamma1=gamma1,gamma2=gamma2,gamma3=gamma3,gamma4=gamma4,gamma5=gamma5,gamma6=gamma6,gamma7=gamma7,gamma8=gamma8,gamma9=gamma9,gamma10=gamma10,gamma11=gamma11,gamma12=gamma12,q.gamma1=q.gamma1,q.gamma2=q.gamma2,q.gamma3=q.gamma3,q.gamma4=q.gamma4,q.gamma5=q.gamma5,q.gamma6=q.gamma6)))
        while (max(dat$time)<timeinterest){
            dat <- as.data.table(lava::sim(model,n = N,p = c(BETA=BETA,OMEGA=OMEGA,GAMMA=GAMMA,alpha1=alpha1,alpha2=alpha2,alpha3=alpha3,alpha4=alpha4,alpha5=alpha5,alpha6=alpha6,alpha7=alpha7,alpha8=alpha8,alpha9=alpha9,alpha10=alpha10,alpha11=alpha11,alpha12=alpha12,q.alpha1=q.alpha1,q.alpha2=q.alpha2,q.alpha3=q.alpha3,q.alpha4=q.alpha4,q.alpha5=q.alpha5,q.alpha6=q.alpha6,beta1=beta1,beta2=beta2,beta3=beta3,beta4=beta4,beta5=beta5,beta6=beta6,beta7=beta7,beta8=beta8,beta9=beta9,beta10=beta10,beta11=beta11,beta12=beta12,q.beta1=q.beta1,q.beta2=q.beta2,q.beta3=q.beta3,q.beta4=q.beta4,q.beta5=q.beta5,q.beta6=q.beta6,omega1=omega1,omega2=omega2,omega3=omega3,omega4=omega4,omega5=omega5,omega6=omega6,omega7=omega7,omega8=omega8,omega9=omega9,omega10=omega10,omega11=omega11,omega12=omega12,q.omega1=q.omega1,q.omega2=q.omega2,q.omega3=q.omega3,q.omega4=q.omega4,q.omega5=q.omega5,q.omega6=q.omega6,gamma1=gamma1,gamma2=gamma2,gamma3=gamma3,gamma4=gamma4,gamma5=gamma5,gamma6=gamma6,gamma7=gamma7,gamma8=gamma8,gamma9=gamma9,gamma10=gamma10,gamma11=gamma11,gamma12=gamma12,q.gamma1=q.gamma1,q.gamma2=q.gamma2,q.gamma3=q.gamma3,q.gamma4=q.gamma4,q.gamma5=q.gamma5,q.gamma6=q.gamma6)))
        }
        dat[,censtime:=NULL]
        if (do.ateRobust){
            dat[,AF:=as.factor(A)]
            if (do.known==TRUE){
                ateR <- ate(data = dat,
                            estimator = c("Gformula","IPTW","AIPTW"),
                            event = list(formula(paste0("Hist(time, event)~AF+",formula.event)),formula(paste0("Hist(time, event)~A+",formula.cr))),
                            censor = formula(paste0("Surv(time, event==0)~AF+",formula.cens)),
                            treatment = formula(paste0("AF~",formula.treatment)),
                            known.nuisance=FALSE,
                            times = timeinterest,
                            product.limit = FALSE,
                            verbose=0,
                            cause = 1)$riskComparison
                ateR <- ateR[,grep("^diff",names(ateR),value=1),with=0]
                ateR.known <- ate(data = dat,
                                  estimator = c("Gformula","IPTW","AIPTW"),
                                  event = list(formula(paste0("Hist(time, event)~AF+",formula.event)),formula(paste0("Hist(time, event)~A+",formula.cr))),
                                  censor = formula(paste0("Surv(time, event==0)~AF+",formula.cens)),
                                  treatment = formula(paste0("AF~",formula.treatment)),
                                  known.nuisance=TRUE,
                                  times = timeinterest,
                                  product.limit = FALSE,
                                  verbose=0,
                                  cause = 1)$riskComparison
                ateR.known <- ateR.known[,grep("^diff",names(ateR.known),value=1),with=0]
            } else{
                ateR <- ate(data = dat,estimator = c("Gformula","IPTW","AIPTW"),event = list(formula(paste0("Hist(time, event)~AF+",formula.event)),formula(paste0("Hist(time, event)~A+",formula.cr))),censor = formula(paste0("Surv(time, event==0)~AF+",formula.cens)),treatment = formula(paste0("AF~",formula.treatment)),known.nuisance=FALSE,times = timeinterest,product.limit = FALSE,verbose=0,cause = 1)$riskComparison
                ateR <- ateR[,grep("^diff",names(ateR),value=1),with=0]
            }
        }
        km <- predictRisk(prodlim(Hist(time, event) ~ A,data=dat),newdata=data.frame(A=0:1),times=timeinterest,cause=1)
        out <- c(KM.naive=diff(km))
        if (do.ateRobust) {
            out <- c(out,unlist(ateR[,grep("a$|W$",names(ateR),value=1),with=0]))
            if (do.known==TRUE){
                out <- c(out,unlist(ateR.known[,grep("a$|W$",names(ateR),value=1),with=0]))
            }
            names(out) <- gsub("diff.","",names(out))
        }
        if (coverage){
            se <- unlist(ateR[,grep(".se$",names(ateR),value=1),with=0])
            names(se) <- gsub("diff.","",names(se))
            ## print(se)
            if (do.known){
                se.known <- unlist(ateR.known[,grep(".se$",names(ateR.known),value=1),with=0])
                names(se.known) <- gsub("diff.","",names(se.known))
                ## print(se.known)
                lower.known <- unlist(ateR.known[,grep(".lower$",names(ateR.known),value=1),with=0])
                upper.known <- unlist(ateR.known[,grep(".upper$",names(ateR.known),value=1),with=0])
                cov.known <- lower.known<true.ate & upper.known>true.ate
                names(cov.known) <- gsub(".lower|diff.","",paste0("coverage.",names(cov.known)))
            }
            lower <- unlist(ateR[,grep(".lower$",names(ateR),value=1),with=0])
            upper <- unlist(ateR[,grep(".upper$",names(ateR),value=1),with=0])
            cov <- lower<true.ate & upper>true.ate
            names(cov) <- gsub(".lower|diff.","",paste0("coverage.",names(cov)))
            out <- c(out,true=true.ate,se,cov)
            if (do.known){
                out <- c(out,se.known=se.known,cov.known=cov.known)
            }
        }
        cens.percent <- dat[,mean(time<timeinterest & event==0)]
        out <- c(cens.percent=cens.percent,true.ate=true.ate,out)
        out
    }
    ## debug(run.one)
    ## browser()
    lava.options(messages=1)
    lava.options(messages=0)
    out <- sim(run.one,B,mc.cores=cores,timeinterest=timeinterest,N=N,coverage=coverage)
    ## })
    out
}

## * table.simulation
table.simulation <- function(x,true,digits=2,...){
    library(Publish)
    res <- do.call("rbind",lapply(x,function(X){
        sumsim <- summary(X,true=true)
        bias <- sumsim[rownames(sumsim)=="Bias",]
        rmse <- sumsim[rownames(sumsim)=="RMSE",]
        out <- matrix(NA,ncol=NCOL(sumsim),nrow=1)
        for (i in 1:length(bias)){
            ## if (i==2)browser()
            if (is.infinite(bias[i]) || abs(bias[i])>100) 
                out[,i] <- "NA (NA)"
            else
                out[,i] <- sprintf(paste0("%1.",digits,"f (%1.",digits,"f)"),bias[i],rmse[i])
        }
        colnames(out) <- names(bias)
        out}))
    res <- cbind("N"=sub("N:","",names(x)),res)
    res
}

## * boxplot.simulation
boxplot.simulation <- function(x,title,which=c("KM.naive","Gformula","IPWnaive","AIPWnaive","TMLE"),...){
    require(ggplot2)
    require(data.table)
    namesx <- names(x)
    dat <- rbindlist(lapply(1:length(x),function(v){
        u <- x[[v]]
        size <- NROW(u)
        out <- data.table(N=sub("N:","",namesx[[v]]),
                          Method=rep(colnames(u),rep(size,NCOL(u))),
                          Value=c(u))
        out[!is.na(Value) & (abs(Value)<1) & !is.infinite(Value)]
    }))
    dat <- dat[!(Method %in% c("IPWefficient","AIPWefficient"))]
    dat[,Method:=factor(Method,levels=c("KM.naive","Gformula","IPWnaive","AIPWnaive","TMLE"),
                        labels=c("Aalen-Johansen","Cox regression (G-formula)","IPW","AIPW","TMLE"))]
    nn <- gsub("N:","",names(x))
    bp <- ggplot(dat,aes(x=factor(N,levels=nn),y=Value,fill=Method))+geom_boxplot()
    bp+xlab("Sample size")+ylab("Treatment effect")
}

######################################################################
### simulation-functions.R ends here
