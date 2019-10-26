### simulation-setting.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Sep 12 2018 (12:16) 
## Version: 
## Last-Updated: apr 29 2019 (18:15) 
##           By: Brice Ozenne
##     Update #: 15
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
## library(drtmle)
##setwd("~/research/Methods/ATE/dropbox/simulation/")
## source("simulation-models.R")
## source("simulation-functions.R")
## B <- 100
## B <- 1000
## N <- c(200,350,500,750,1000)
## N <- 500
## timeinterest <- 10
#----------------------------------------------------------------------
##  A. Observational study bias correction
##     Model 2, Scenario 1: effect in randomized experiment higher than in observed
## Treatment model
##          A ~ alpha * W[1:12]+ q.alpha * (W[1:6])^2
## Outcome of interest
##          T1 ~ BETA * A + beta * W[1:12]+ q.beta * (W[1:6])^2
## Competing risk
##          T2 ~ OMEGA * A + omega * W[1:12]+ q.omega * (W[1:6])^2
## Censoring
##          C ~ GAMMA * A + gamma * W[1:12]+ q.gamma * (W[1:6])^2
#----------------------------------------------------------------------
## Check case without effects
## setup.NULL <- c(BETA=0,OMEGA=0,GAMMA=0,alpha=rep(0,12),q.alpha=rep(0,6),beta=rep(0,12),q.beta=rep(0,6),omega=rep(0,12),q.omega=rep(0,6),gamma=rep(0,12),q.gamma=rep(0,6))
## do.call("plot.model2",list(parms=setup.NULL,N=10000,timeinterest=10))
## Case 0: all models wrong: quadratic effects 
## setup.A <- c(BETA=0,OMEGA=0,GAMMA=0,alpha=rep(c(0,-1,-1,0,1,-1),2) * log(1.5),q.alpha=c(1,0,-1,0,1,-1) * log(1.5),beta=rep(c(1,-1,0,0,1,-1),2) * log(1.5),q.beta=c(1,-1,0,0,1,-1) * log(1.5),omega=rep(c(0,0,1,-1,1,-1),2) * log(1.5),q.omega=c(0,0,1,-1,1,-1) * log(1.5),gamma=rep(c(1,1,1,1,1,-1),2) * log(1.5),q.gamma=c(0,0,1,-1,1,-1) * log(1.5),scale.cr=1/500)
## W6, W12 unmeasured confounders
## W5, W11 only comp.risk
## W4, W10 only cens
## W3, W9 only outcome
## W2, W8 only treatment
## W1, W7 all
setup.A <- list(true.ate=NULL,
                parameters=c(BETA=0,
                             OMEGA=0,
                             GAMMA=0,
                             alpha=rep(c(1,-1,0,0,0,1),2) * log(2),
                             q.alpha=c(1,-1,0,0,0,-1) * log(2),
                             beta=rep(c(1,0,1,0,0,1),2) * log(2),
                             q.beta=c(1,0,-1,0,0,-1) * log(2),
                             omega=rep(c(-1,0,0,0,1,1),2) * log(2),
                             q.omega=c(1,0,0,0,-1,1) * log(2),
                             gamma=rep(c(-1,0,0,1,0,-1),2) * log(2),
                             q.gamma=c(1,0,0,-1,0,-1) * log(2),
                             scale.cr=1/50),
                formula=list(formula.treatment="W1+W2+W3+W4+W5+W6+W7+W8+W9+W10+W11+W12",
                             formula.event="W1+W2+W3+W4+W5+W6+W7+W8+W9+W10+W11+W12",
                             formula.cr="W1+W2+W3+W4+W5+W6+W7+W8+W9+W10+W11+W12",
                             formula.cens="W1+W2+W3+W4+W5+W6+W7+W8+W9+W10+W11+W12"))
## do.call("plot.model2",list(parms=setup.A$parameters,N=20000,timeinterest=10))
## Figure 2: variation of models specification (N=500,30% censored)
## correct formula and no quadratic effects
parms.correct <- setup.A$parameters
parms.correct[grepl("^q.",names(parms.correct))] <- 0
form.correct <- setup.A$formula
## quadratic effects in treatment model
parms.mistreat <- parms.correct
parms.mistreat[grepl("alpha.",names(parms.mistreat))] <- setup.A$parameters[grepl("alpha.",names(parms.mistreat))]
## omitting covariates W6,W12 from treatment model
form.mistreat <- list(formula.treatment="W2+W3+W4+W5+W7+W8+W9+W10+W11",
                      formula.event="W1+W2+W3+W4+W5+W6+W7+W8+W9+W10+W11+W12",
                      formula.cr="W1+W2+W3+W4+W5+W6+W7+W8+W9+W10+W11+W12",
                      formula.cens="W1+W2+W3+W4+W5+W6+W7+W8+W9+W10+W11+W12")
## quadratic effects in event models
parms.misout <- parms.correct
parms.misout[grepl("beta.",names(parms.misout))] <- setup.A$parameters[grepl("beta.",names(parms.misout))]
form.misout <- list(formula.treatment="W1+W2+W3+W4+W5+W6+W7+W8+W9+W10+W11+W12",
                    formula.event="W2+W3+W4+W5+W7+W8+W9+W10+W11",
                    formula.cr="W1+W2+W3+W4+W5+W6+W7+W8+W9+W10+W11+W12",
                    formula.cens="W1+W2+W3+W4+W5+W6+W7+W8+W9+W10+W11+W12")
## quadratic effects in comprisk models
parms.miscr <- parms.correct
parms.miscr[grepl("omega.",names(parms.miscr))] <- setup.A$parameters[grepl("omega.",names(parms.miscr))]
## omitting covariates W6,W12 from outcome model
form.miscr <- list(formula.treatment="W1+W2+W3+W4+W5+W6+W7+W8+W9+W10+W11+W12",
                    formula.event="W1+W2+W3+W4+W5+W6+W7+W8+W9+W10+W11+W12",
                    formula.cr="W2+W3+W4+W5+W7+W8+W9+W10+W11",
                    formula.cens="W1+W2+W3+W4+W5+W6+W7+W8+W9+W10+W11+W12")
## quadratic effects in censoring models
parms.miscens <- parms.correct
parms.miscens[grepl("gamma.",names(parms.miscens))] <- setup.A$parameters[grepl("gamma.",names(parms.miscens))]
## omitting covariates W6,W12 from censoring model
form.miscens <- list(formula.treatment="W1+W2+W3+W4+W5+W6+W7+W8+W9+W10+W11+W12",
                     formula.event="W1+W2+W3+W4+W5+W6+W7+W8+W9+W10+W11+W12",
                     formula.cr="W1+W2+W3+W4+W5+W6+W7+W8+W9+W10+W11+W12",
                     formula.cens="W2+W3+W4+W5+W7+W8+W9+W10+W11")
setup.Figure2 <- list("All models correct"=list(parameters=parms.correct,formula=form.correct),
                      "Misspecified treatment model"=list(parameters=parms.mistreat,formula=form.mistreat),
                      "Misspecified outcome model"=list(parameters=parms.misout,formula=form.misout),
                      ## "misspecified comprisk model"=list(parameters=parms.miscr,formula=form.miscr),
                      "Misspecified censoring model"=list(parameters=parms.miscens,formula=form.miscens))
## sample size & coverage
setup.Figure3 <- list(parameters=parms.miscr,formula=form.miscr)

## censoring & coverage
setup.Figure4 <- list(parameters=parms.miscr,formula=form.miscr)

## comprisk model misspecified
setup.Figure5 <- list(parameters=parms.miscr,formula=form.miscr)

## parms <- setup.Figure2[[1]]$parameters
setup.Figure6 <- list(parameters=c(BETA=-2,
                                   OMEGA=0,
                                   GAMMA=0,
                                   alpha=rep(c(0,-1,-1,0,1,-1),2) * log(2),
                                   q.alpha=c(1,0,-1,0,1,-1) * log(2),
                                   beta=rep(c(1,-1,0,0,1,-1),2) * log(2),
                                   q.beta=c(1,-1,0,0,1,-1) * log(2),
                                   omega=rep(c(0,0,1,-1,1,-1),2) * log(2),
                                   q.omega=c(0,0,1,-1,1,-1) * log(2),
                                   gamma=rep(0,12),
                                   q.gamma=rep(0,12),
                                   scale.cr=1/50),
                      formula=list(formula.treatment="W1+W2+W3+W4+W5+W6+W7+W8+W9+W10+W11+W12",
                                   formula.event="W1+W2+W3+W4+W5+W6+W7+W8+W9+W10+W11+W12",
                                   formula.cr="W1+W2+W3+W4+W5+W6+W7+W8+W9+W10+W11+W12",
                                   formula.cens="W1+W2+W3+W4+W5+W6+W7+W8+W9+W10+W11+W12"))


######################################################################
### simulation-setting.R ends here
