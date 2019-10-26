### simulation-run.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Aug 17 2018 (11:53) 
## Version: 
## Last-Updated: Oct 24 2019 (19:41) 
##           By: Thomas Alexander Gerds
##     Update #: 239
#----------------------------------------------------------------------
## 
### Commentary: 
## Running the empirical studies for the competing-risk-ATE manuscript
##
##  Figure 1. Simulation scenario
##
##  Figure 2. Observational study bias correction.
##
##  Context: there is no treatment effect but naive estimate based on observational data is large
## 
##  - all models correctly specified
##  - outcome model misspecified
##  - treatment model misspecified
##  - competing risk model misspecified
##  - censoring model misspecified
##
##  Figure 3. Coverage AIPTW.AIPCW 
##
##  - increasing sample size
## 
##  Figure A1. effect of IPCW
##  
##  - increasing censoring percentage
##
##  Figure A2. effect of misspecified competing risk model
##  
##  - varying effect on competing risk hazard 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
setwd("~/research/Methods/ATE/dropbox/simulation/")
source("simulation-models.R")
source("simulation-setting.R")
source("simulation-functions.R")

## Figure 1:
if (FALSE){
    setup.A <- c(BETA=0,OMEGA=0,GAMMA=0,alpha=rep(c(0,-1,-1,0,1,-1),2) * log(2),q.alpha=c(1,0,-1,0,1,-1) * log(2),beta=rep(c(1,-1,0,0,1,-1),2) * log(2),q.beta=c(1,-1,0,0,1,-1) * log(2),omega=rep(c(0,0,1,-1,1,-1),2) * log(2),q.omega=c(0,0,1,-1,1,-1) * log(2),gamma=rep(0,12),q.gamma=rep(0,12))
    setup.B <- c(BETA=-2,OMEGA=0,GAMMA=0,alpha=rep(c(0,-1,-1,0,1,-1),2) * log(2),q.alpha=c(1,0,-1,0,1,-1) * log(2),beta=rep(c(1,-1,0,0,1,-1),2) * log(2),q.beta=c(1,-1,0,0,1,-1) * log(2),omega=rep(c(0,0,1,-1,1,-1),2) * log(2),q.omega=c(0,0,1,-1,1,-1) * log(2),gamma=rep(0,12),q.gamma=rep(0,12))
    par(mfrow=c(1,2))
    set.seed(9)
    do.call("plot.model2",list(parms=setup.A,N=10000,timeinterest=8))
    mtext("A",cex=2,at=-2,xpd=NA,line=1.5)
    set.seed(10)
    do.call("plot.model2",list(parms=setup.B,N=10000,timeinterest=8))
    mtext("B",cex=2,at=-2,xpd=NA,line=1.5)
}
## Figure 2: bias under misspecified models (fix sample size, fix censoring %)
if (FALSE){
    set.seed(1)
    cores <- min(detectCores()-1,25)
    B <- 1000
    results.Figure2 <- lapply(names(setup.Figure2),function(name){
        print(name)
        set <- setup.Figure2[[name]]
        ## source("simulation-functions.R")
        run.simulation(model=model2,setup=set,causalgam=FALSE,N=500,timeinterest=10,B=B,cores=cores,true.ate=0,coverage=TRUE)
    })
    names(results.Figure2) <- names(setup.Figure2)
    saveRDS(results.Figure2,file="simulation-results/results-Figure2-Revision.rds")
}


## Figure 3 (coverage increasing sample size)
if (FALSE){
    cores <- min(detectCores()-1,25)
    B <- 1000
    Nseq <- c(100,150,250,350,500,750,1000,1500)
    for (setting in names(setup.Figure2)[c(3,4,1,2)]){
        results.FigureNa <- lapply(Nseq,function(N){
            print(N)
            set <- setup.Figure2[[setting]]
            if (N==2000) cores <- 10 else cores <- min(detectCores()-1,25)
            set.seed(19)
            res <- run.simulation(model=model2,
                                  setup=set,
                                  causalgam=FALSE,
                                  do.tmle=FALSE,
                                  N=N,
                                  timeinterest=8,
                                  B=B,
                                  cores=cores,
                                  true.ate=0,
                                  coverage=TRUE)
            xclass <- class(res)
            res <- cbind(res,N=N)
            class(res) <- xclass
            res
        })
        names(results.FigureNa) <- paste0("N",Nseq)
        saveRDS(results.FigureNa,file=paste0("simulation-results/results-FigureN-Revision-set-",gsub(" ","-",setting),".rds"))
    }
    ## saveRDS(results.FigureNa,file="simulation-results/results-FigureN-Revision.rds")
}

## Figure 3a (se misspecified known vs unknown nuissancemodels increasing sample size)
if (TRUE){
    cores <- min(detectCores()-1,25)
    B <- 5000
    Nseq <- c(100,200,500,750,1000)
    results.FigureSE <- lapply(Nseq,function(N){
        print(N)
        ## set <- setup.Figure2[["Misspecified treatment model"]]
        set <- setup.Figure2[["All models correct"]]
        if (N==2000) cores <- 10 else cores <- min(detectCores()-1,25)
        set.seed(19)
        ## source("simulation-functions.R")
        res <- run.simulation1(model=model2,setup=set,do.known=TRUE,N=N,timeinterest=8,B=B,cores=cores,true.ate=0,coverage=TRUE)
        xclass <- class(res)
        res <- cbind(res,N=N)
        class(res) <- xclass
        res
        ## saveRDS(res,file=paste0("simulation-results/results-FigureSE-ReRevision-mistreat-N-",N,".rds"))
        saveRDS(res,file=paste0("simulation-results/results-FigureSE-ReRevision-correct-N-",N,".rds"))
    })
}
if (TRUE){
    cores <- min(detectCores()-1,25)
    B <- 5000
    Nseq <- c(100,200,500,750,1000)
    results.FigureSE <- lapply(Nseq,function(N){
        print(N)
        ## set <- setup.Figure2[["Misspecified outcome model"]]
        set <- setup.Figure2[["Misspecified censoring model"]]
        if (N==2000) cores <- 10 else cores <- min(detectCores()-1,25)
        set.seed(19)
        ## source("simulation-functions.R")
        res <- run.simulation1(model=model2,setup=set,do.known=TRUE,N=N,timeinterest=8,B=B,cores=cores,true.ate=0,coverage=TRUE)
        xclass <- class(res)
        res <- cbind(res,N=N)
        class(res) <- xclass
        res
        ## saveRDS(res,file=paste0("simulation-results/results-FigureSE-ReRevision-misout-N-",N,".rds"))
        saveRDS(res,file=paste0("simulation-results/results-FigureSE-ReRevision-miscens-N-",N,".rds"))
    })
}

## a <- readRDS("~/research/Methods/ATE/dropbox/simulation/simulation-results/results-FigureSE-ReRevision-N-500.rds")
## a <- readRDS("c:/Users/hpl802/Dropbox/competing-risk-ATE/simulation/simulation-results/results-FigureSE-ReRevision-N-500.rds")

## dt.res <- as.data.table(unclass(a))
## dtRed.res <- dt.res[abs(AIPTW) < 1,] ## rm outliers
## boxplot(dtRed.res$AIPTW)
## dev.new()
## qqtest::qqtest(dtRed.res$AIPTW)
## qqtest::qqtest(atanh(dtRed.res$AIPTW))
## colnames(Ma)
## MASS::fitdistr(dtRed.res$AIPTW, "t")
## dtRed.res[, lower.AIPTW := AIPTW - 2*AIPTW.se]
## dtRed.res[, upper.AIPTW := AIPTW + 2*AIPTW.se]
## dtRed.res[, coverage2.AIPTW := (true>=lower.AIPTW)*(true<=upper.AIPTW)]
## table(dtRed.res$coverage.AIPTW,dtRed.res$coverage2.AIPTW)
## original scale
## dtRed.res[,.(se.truth = sd(AIPTW, na.rm = TRUE),
             ## se.estimatedFull = mean(AIPTW.se, na.rm = TRUE),
             ## se.estimated0 = mean(se.known.AIPTW.se, na.rm = TRUE),
             ## cov.estimatedFull = mean(coverage.AIPTW, na.rm = TRUE),
             ## cov.estimated0 = mean(cov.known.coverage.AIPTW, na.rm = TRUE)
             ## )]
## atanh transformation
## dtRed.res[,.(se.truth = sd(atanh(AIPTW), na.rm = TRUE),
             ## se.estimatedFull = mean(AIPTW.se/(1-AIPTW^2), na.rm = TRUE),
             ## se.estimated0 = mean(se.known.AIPTW.se/(1-AIPTW^2), na.rm = TRUE),
             ## cov.estimatedFull = mean(coverage.AIPTW, na.rm = TRUE),
             ## cov.estimated0 = mean(cov.known.coverage.AIPTW, na.rm = TRUE)
             ## )]

## atanh(-1)
## sd(rt(1e3,df = 5))

## Figure 5 (varying competing risk effect)
if (FALSE){
    set.seed(1)
    cores <- min(detectCores()-1,25)
    B <- 1000
    N <- 500
    results.Figure.omega1 <- lapply(c(1/3,1/2,1,1.5,2,3),function(om){
        print(om)
        set <- list(parameters=parms.miscr,formula=form.miscr)
        set$parameters["omega1"] <- log(om)
        print(set$parameters["omega1"])
        run.simulation(model=model2,
                       setup=set,
                       causalgam=FALSE,
                       N=N,
                       timeinterest=8,
                       B=B,
                       cores=cores,
                       true.ate=0,
                       coverage=FALSE,
                       do.tmle=FALSE)
        res <- run.simulation(model=model2,setup=set,causalgam=FALSE,N=N,timeinterest=8,B=B,cores=cores,true.ate=0,coverage=FALSE,do.tmle=FALSE)
        xclass <- class(res)
        res <- cbind(res,omega1=log(om))
        class(res) <- xclass
        res
    })
    saveRDS(results.Figure.omega1,file="simulation-results/results-Figure.omega1.rds")
}

## Figure 4 (increasing censoring)
if (FALSE){
    set.seed(1)
    cores <- min(detectCores()-1,25)
    B <- 1000
    N <- 500
    results.Figure.cens <- lapply(c(1/1000,1/400,1/100,1/50,1/10,1/5),function(cr){
        print(N)
        set <- setup.Figure2[[1]]
        set$parameters["scale.cr"] <- cr
        run.simulation(model=model2,
                       setup=set,
                       causalgam=FALSE,
                       N=N,
                       timeinterest=5,
                       B=B,
                       cores=cores,
                       true.ate=0,
                       coverage=TRUE)
    })
    saveRDS(results.Figure.cens,file="simulation-results/results-Figure.cens.rds")
}

## Figure R1 (correlated covariates)
if (FALSE){
    set.seed(1)
    cores <- min(detectCores()-1,25)
    B <- 1000
    N <- 500
    for (setting in names(setup.Figure2)){
        ## zero corresponds to Figure 2
        results.Figure.corr <- lapply(c(0,1,1.5,2),function(ccc){ 
            print(ccc)
            set3 <- setup.Figure2[[setting]]
            model3 <- model2
            covariance.effects <- round(rnorm(66,mean=0,sd=ccc),2)
            covariates <- c("one","two","three","four","five","six","seven","eight","nine","ten","eleven","twelve")
            names(covariance.effects) <- unlist(lapply(2:12,function(x){paste0(covariates[x],1:(x-1))}))
            set3$parameters <- c(set3$parameters,covariance.effects)
            # W1 depends on nothing 
            regression(model3,W2~two1*W1)
            regression(model3,W3~three1*W1+three2*W2)
            regression(model3,W4~four1*W1+four2*W2+four3*W3)
            regression(model3,W5~five1*W1+five2*W2+five3*W4+five4*W4)
            regression(model3,W6~six1*W1+six2*W2+six3*W4+six4*W4+six5*W5)
            regression(model3,W7~seven1*W1+seven2*W2+seven3*W4+seven4*W4+seven5*W5+seven6*W6)
            regression(model3,W8~eight1*W1+eight2*W2+eight3*W4+eight4*W4+eight5*W5+eight6*W6+eight7*W7)
            regression(model3,W9~nine1*W1+nine2*W2+nine3*W4+nine4*W4+nine5*W5+nine6*W6+nine7*W7+nine8*W8)
            regression(model3,W10~ten1*W1+ten2*W2+ten3*W4+ten4*W4+ten5*W5+ten6*W6+ten7*W7+ten8*W8+ten9*W9)
            regression(model3,W11~eleven1*W1+eleven2*W2+eleven3*W4+eleven4*W4+eleven5*W5+eleven6*W6+eleven7*W7+eleven8*W8+eleven9*W9+eleven10)
            regression(model3,W12~twelve1*W1+twelve2*W2+twelve3*W4+twelve4*W4+twelve5*W5+twelve6*W6+twelve7*W7+twelve8*W8+twelve9*W9+twelve10*W10+twelve11*W11)
            run.simulation(model=model3,setup=set3,causalgam=FALSE,N=N,timeinterest=5,B=B,cores=cores,true.ate=0,coverage=TRUE)
        })
        names(results.Figure.corr) <- paste0("corr.",c(0,1,1.5,2))
        saveRDS(results.Figure.corr,file=paste0("simulation-results/results-Figure-corr-set-",gsub(" ","-",setting),".rds"))
    }
}

## Figure R2 (near violation of positivity)
if (FALSE){
    set.seed(1)
    cores <- min(detectCores()-1,25)
    B <- 1000
    N <- 500
    results.Figure.positivity <- lapply(c(-5,-4,-3,3,4,5),function(ia){
        print(ia)
        set <- setup.Figure2[[1]]
        intercept(model2,"A") <- ia
        run.simulation(model=model2,setup=set,do.tmle=FALSE,causalgam=FALSE,N=N,timeinterest=5,B=B,cores=cores,true.ate=0,coverage=TRUE)
    })
    names(results.Figure.positivity) <- paste0("probA-intercept-",c(-5,-4,-3,3,4,5))
    saveRDS(results.Figure.positivity,file="simulation-results/results-Figure-positivity.rds")
}




######################################################################
### simulation-run.R ends here
