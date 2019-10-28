### simulation-models.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Aug 16 2018 (08:51) 
## Version: 
## Last-Updated: okt 28 2019 (10:59) 
##           By: Brice Ozenne
##     Update #: 25
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(lava)

## * define scenario under the null (A) and the alternative (B)
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

setup.B <- list(parameters=c(BETA=-2,
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


## * Generative model for the data 
model2 <- lvm()
W <- paste("W",seq(12),sep="")
q.vars <- paste0("W",1:6)
eventtime1 <- "eventtime1"
eventtime2 <- "eventtime2"
censtime <- "censtime"
latent(model2) <- c(eventtime1,eventtime2,censtime)
addvar(model2) <- c(eventtime1,eventtime2,censtime,W)
exogenous(model2) <- c(W,q.vars)
## covariate distribution
distribution(model2,"A") <- binomial.lvm()
distribution(model2,W[1:6]) <- normal.lvm()
distribution(model2,W[7:12]) <- binomial.lvm()
## outcome distribution
distribution(model2,eventtime1) <- coxWeibull.lvm(scale=1/100)
distribution(model2,eventtime2) <- coxWeibull.lvm(scale=1/100)
distribution(model2,censtime) <- coxWeibull.lvm(scale=1/500)
## quadratic terms
for (v in q.vars){
    form.v <- formula(paste("quad.",v,"~",v,sep=""))
    transform(model2,form.v) <- function(x)x^2
}
quad.W <- paste0("quad.",q.vars)
## regression effects
regression(model2) <- A~ alpha1*W1+alpha2*W2+alpha3*W3+alpha4*W4+alpha5*W5+alpha6*W6+alpha7*W7+alpha8*W8+alpha9*W9+alpha10*W10+alpha11*W11+alpha12*W12+q.alpha1*quad.W1+q.alpha2*quad.W2+q.alpha3*quad.W3+q.alpha4*quad.W4+q.alpha5*quad.W5+q.alpha6*quad.W6
regression(model2) <- eventtime1~ BETA*A+beta1*W1+beta2*W2+beta3*W3+beta4*W4+beta5*W5+beta6*W6+beta7*W7+beta8*W8+beta9*W9+beta10*W10+beta11*W11+beta12*W12+q.beta1*quad.W1+q.beta2*quad.W2+q.beta3*quad.W3+q.beta4*quad.W4+q.beta5*quad.W5+q.beta6*quad.W6
regression(model2) <- eventtime2~ OMEGA*A+omega1*W1+omega2*W2+omega3*W3+omega4*W4+omega5*W5+omega6*W6+omega7*W7+omega8*W8+omega9*W9+omega10*W10+omega11*W11+omega12*W12+q.omega1*quad.W1+q.omega2*quad.W2+q.omega3*quad.W3+q.omega4*quad.W4+q.omega5*quad.W5+q.omega6*quad.W6
regression(model2) <- censtime~ GAMMA*A+gamma1*W1+gamma2*W2+gamma3*W3+gamma4*W4+gamma5*W5+gamma6*W6+gamma7*W7+gamma8*W8+gamma9*W9+gamma10*W10+gamma11*W11+gamma12*W12+q.gamma1*quad.W1+q.gamma2*quad.W2+q.gamma3*quad.W3+q.gamma4*quad.W4+q.gamma5*quad.W5+q.gamma6*quad.W6
model2 <- eventTime(model2,time~min(eventtime1=1,eventtime2=2,censtime=0),"event")

## * Investigator model used to estimate ATE
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
                      "Misspecified censoring model"=list(parameters=parms.miscens,formula=form.miscens))

######################################################################
### simulation-models.R ends here
