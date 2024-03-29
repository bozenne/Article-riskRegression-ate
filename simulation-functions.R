### simulation-functions.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Aug 16 2018 (08:47) 
## Version: 
## Last-Updated: nov 27 2019 (14:10) 
##           By: Brice Ozenne
##     Update #: 350
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(data.table)
library(lava)
library(prodlim)

library(riskRegression)
library(survival)

library(parallel)

## * run.simulation
## simulate from model and apply alternative estimators
run.simulation <- function(model,
                           setup,
                           B=1,
                           N,
                           timeinterest,
                           do.known=FALSE,
                           coverage=TRUE,
                           true.ate,
                           cores=1){

    ## ** prepare
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

    ## ** warper for one simulation
    run.one <- function(...,N,timeinterest,coverage=coverage){

        ## *** generate data
        ## vector of parameter
        vec.p <- c(BETA=BETA,OMEGA=OMEGA,GAMMA=GAMMA,alpha1=alpha1,alpha2=alpha2,alpha3=alpha3,alpha4=alpha4,alpha5=alpha5,alpha6=alpha6,alpha7=alpha7,alpha8=alpha8,alpha9=alpha9,alpha10=alpha10,alpha11=alpha11,alpha12=alpha12,q.alpha1=q.alpha1,q.alpha2=q.alpha2,q.alpha3=q.alpha3,q.alpha4=q.alpha4,q.alpha5=q.alpha5,q.alpha6=q.alpha6,beta1=beta1,beta2=beta2,beta3=beta3,beta4=beta4,beta5=beta5,beta6=beta6,beta7=beta7,beta8=beta8,beta9=beta9,beta10=beta10,beta11=beta11,beta12=beta12,q.beta1=q.beta1,q.beta2=q.beta2,q.beta3=q.beta3,q.beta4=q.beta4,q.beta5=q.beta5,q.beta6=q.beta6,omega1=omega1,omega2=omega2,omega3=omega3,omega4=omega4,omega5=omega5,omega6=omega6,omega7=omega7,omega8=omega8,omega9=omega9,omega10=omega10,omega11=omega11,omega12=omega12,q.omega1=q.omega1,q.omega2=q.omega2,q.omega3=q.omega3,q.omega4=q.omega4,q.omega5=q.omega5,q.omega6=q.omega6,gamma1=gamma1,gamma2=gamma2,gamma3=gamma3,gamma4=gamma4,gamma5=gamma5,gamma6=gamma6,gamma7=gamma7,gamma8=gamma8,gamma9=gamma9,gamma10=gamma10,gamma11=gamma11,gamma12=gamma12,q.gamma1=q.gamma1,q.gamma2=q.gamma2,q.gamma3=q.gamma3,q.gamma4=q.gamma4,q.gamma5=q.gamma5,q.gamma6=q.gamma6)
        ## update scale parameter
        distribution(model,censtime) <- coxWeibull.lvm(scale=scale.cr)

        ## try to simulate until a
        i.max <- 100
        dat <- as.data.table(lava::sim(model, n = N, p = vec.p))

        i <- 1
        
        while (dat[, max(time)<timeinterest] && i < i.max){ ## check if at least one "event" after the prediction time otherwise re-simulate
            dat <- as.data.table(lava::sim(model, n = N, p = vec.p))
            i <- i + 1
        }
        dat[,censtime:=NULL]
        dat[,AF:=as.factor(A)]

        ## naive non-parametric estimator: Kaplan Meier or Aalen Johansen (survival or competing risk setting)
        e.prodlim <- prodlim(Hist(time, event) ~ A, data=dat) ## fit estimator
        pred.prodlim <- predictRisk(e.prodlim, newdata=data.frame(A=0:1), times=timeinterest, cause=1) ## compute risks 
        out <- c(KM.naive = diff(pred.prodlim))

        ## adjusted estimator based estimating the average treatment effect
        ateR <- ate(data = dat,
                    estimator = c("Gformula","IPTW","AIPTW"),
                    event = list(formula(paste0("Hist(time, event)~AF+",formula.event)),formula(paste0("Hist(time, event)~A+",formula.cr))),
                    censor = formula(paste0("Surv(time, event==0)~AF+",formula.cens)),
                    treatment = formula(paste0("AF~",formula.treatment)),
                    known.nuisance = FALSE,
                    times = timeinterest,
                    product.limit = FALSE,
                    se = coverage, iid = FALSE, band = FALSE,
                    verbose=0,
                    cause = 1)$riskComparison
        ateR <- ateR[,grep("^diff",names(ateR),value=1),with=0]
        out <- c(out,unlist(ateR[,grep("a$|W$",names(ateR),value=1),with=0]))

        if (do.known==TRUE){
            ateR.known <- ate(data = dat,
                              estimator = c("Gformula","IPTW","AIPTW"),
                              event = list(formula(paste0("Hist(time, event)~AF+",formula.event)),formula(paste0("Hist(time, event)~A+",formula.cr))),
                              censor = formula(paste0("Surv(time, event==0)~AF+",formula.cens)),
                              treatment = formula(paste0("AF~",formula.treatment)),
                              known.nuisance = TRUE,
                              times = timeinterest,
                              product.limit = FALSE,
                              se = coverage, iid = FALSE, band = FALSE,
                              verbose=0,
                              cause = 1)$riskComparison
            ateR.known <- ateR.known[,grep("^diff",names(ateR.known),value=1),with=0]
            out <- c(out,unlist(ateR.known[,grep("a$|W$",names(ateR),value=1),with=0]))
        }        

        ## adjust name
        names(out) <- gsub("diff.","",names(out))

        ## compute coverage of the CIs
        if (coverage){
            se <- unlist(ateR[,grep(".se$",names(ateR),value=1),with=0])
            names(se) <- gsub("diff.","",names(se))

            lower <- unlist(ateR[,grep(".lower$",names(ateR),value=1),with=0])
            upper <- unlist(ateR[,grep(".upper$",names(ateR),value=1),with=0])
            cov <- lower<true.ate & upper>true.ate
            names(cov) <- gsub(".lower|diff.","",paste0("coverage.",names(cov)))
            out <- c(out,true=true.ate,se,cov)
            
            if (do.known){
                se.known <- unlist(ateR.known[,grep(".se$",names(ateR),value=1),with=0])
                names(se.known) <- gsub("diff.","",names(se.known))

                lower.known <- unlist(ateR.known[,grep(".lower$",names(ateR.known),value=1),with=0])
                upper.known <- unlist(ateR.known[,grep(".upper$",names(ateR.known),value=1),with=0])
                cov.known <- lower.known<true.ate & upper.known>true.ate
                names(cov.known) <- gsub(".lower|diff.","",paste0("coverage.",names(cov.known)))
                out <- c(out,se.known=se.known,cov.known=cov.known)
            }
        }
     
        ## store percentage of censoring in the data and export
        cens.percent <- dat[,mean(time<timeinterest & event==0)]
        out <- c(cens.percent=cens.percent,true.ate=true.ate,out)
        return(out)
    }

    ## ** run simulations
    res <- lava::sim(run.one,
                     B,
                     mc.cores=cores,
                     timeinterest=timeinterest,
                     N=N,
                     coverage=coverage)
    return(res)
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

## * simulate.data
simulate.data <- function(){

    ## group 1
    m.dabigatran <- lvm()
    distribution(m.dabigatran, ~male) <- binomial.lvm(p = 4507/7078, size = 1)
    categorical(m.dabigatran, labels=c("<65","65-74","75-84",">85"), p = c(2697,3505,864,12)/7078) <- "age"
    transform(m.dabigatran, age.num ~ age) <- function(x){
        as.numeric(factor(x[,1], levels = c("<65","65-74","75-84",">85")))
    }
    transform(m.dabigatran, age.num2 ~ age) <- function(x){x^2}
    distribution(m.dabigatran, ~heart.failure) <- binomial.lvm(p = 843/7078, size = 1)
    distribution(m.dabigatran, ~hypertension) <- binomial.lvm(p = 2939/7078, size = 1)
    distribution(m.dabigatran, ~diabetes)  <- binomial.lvm(p = 751/7078, size = 1)
    distribution(m.dabigatran, ~prior.stroke) <- binomial.lvm(p = 833/7078, size = 1)
    distribution(m.dabigatran, ~prior.bleeding) <- binomial.lvm(p = 583/7078, size = 1)

    distribution(m.dabigatran, "eventtime1") <- lava::coxWeibull.lvm(scale = 1/2000, shape = 1.3)
    distribution(m.dabigatran, "eventtime2") <- lava::coxWeibull.lvm(scale = 1/2000)
    distribution(m.dabigatran, "censtime") <- lava::coxWeibull.lvm(scale = 1/5000)
    m.dabigatran <- lava::eventTime(m.dabigatran, time ~ min(eventtime1 = 1, eventtime2 = 2, 
                                                             censtime = 0), "event")
    latent(m.dabigatran) <- ~eventtime1+eventtime2+censtime+age.num+age.num2
    ## regression(m.dabigatran) <- eventtime1 ~ 1*prior.stroke+0.1*heart.failure+0.1*hypertension+0.1*diabetes+0.25*age.num+0.25*age.num2
    ## regression(m.dabigatran) <- eventtime2 ~ 0.1*prior.stroke+0.5*heart.failure+0.2*hypertension+0.2*diabetes

    ## group 2
    m.rivaroxaban <- lvm()
    distribution(m.rivaroxaban, ~male) <- binomial.lvm(p = 3807/6868, size = 1)
    categorical(m.rivaroxaban, labels=c("<65","65-74","75-84",">85"), p = c(1567,2750,1898,653)/6868) <- "age"
    transform(m.rivaroxaban, age.num ~ age) <- function(x){
        as.numeric(factor(x[,1], levels = c("<65","65-74","75-84",">85")))-1
    }
    transform(m.rivaroxaban, age.num2 ~ age) <- function(x){x^2}
    distribution(m.rivaroxaban, ~heart.failure) <- binomial.lvm(p = 821/6868, size = 1)
    distribution(m.rivaroxaban, ~hypertension) <- binomial.lvm(p = 2989/6868, size = 1)
    distribution(m.rivaroxaban, ~diabetes)  <- binomial.lvm(p = 774/6868, size = 1)
    distribution(m.rivaroxaban, ~prior.stroke) <- binomial.lvm(p = 1025/6868, size = 1)
    distribution(m.rivaroxaban, ~prior.bleeding) <- binomial.lvm(p = 656/6868, size = 1)

    distribution(m.rivaroxaban, "eventtime1") <- lava::coxWeibull.lvm(scale = 1/3000, shape = 1.7)
    distribution(m.rivaroxaban, "eventtime2") <- lava::coxWeibull.lvm(scale = 1/3000)
    distribution(m.rivaroxaban, "censtime") <- lava::coxWeibull.lvm(scale = 1/5000)
    m.rivaroxaban <- lava::eventTime(m.rivaroxaban, time ~ min(eventtime1 = 1, eventtime2 = 2, 
                                                             censtime = 0), "event")
    latent(m.rivaroxaban) <- ~eventtime1+eventtime2+censtime+age.num+age.num2
    ## regression(m.rivaroxaban) <- eventtime1 ~ 1*prior.stroke+0.1*heart.failure+0.1*hypertension+0.1*diabetes+0.25*age.num+0.25*age.num2
    ## regression(m.rivaroxaban) <- eventtime2 ~ 0.1*prior.stroke+0.5*heart.failure+0.2*hypertension+0.2*diabetes

    ## simulate data
    df <- rbind(cbind(lava::sim(m.dabigatran, n = 7078, latent = FALSE), treatment = "Dabigatran"),
                cbind(lava::sim(m.rivaroxaban, n = 6868, latent = FALSE), treatment = "Rivaroxaban")
                )

    return(as.data.table(df)[])
}
######################################################################
### simulation-functions.R ends here
