### simulation-tests.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Sep 12 2018 (12:13) 
## Version: 
## Last-Updated: Oct 18 2019 (12:57) 
##           By: Thomas Alexander Gerds
##     Update #: 12
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

setwd("~/research/Methods/ATE/dropbox/simulation/")
source("simulation-functions.R")
source("simulation-setting.R")
source("simulation-models.R")

u <- get(load("~/tmp/u.rda"))
## f <- do.call("CSC",c(u["data"],list(formula=u[["event"]])))
## u$event <- f
u$data[,A:=as.factor(A)]
source("simulation-functions.R");
x <- do.call("ate",u)
as.data.table(x)
## y <- ateRobust(data = u$dat,type = "competing.risks",formula.event = u[["event"]],formula.censor = u[["censor"]],formula.treatment = u[["treatment"]],augment.cens=TRUE,times = 10,product.limit = FALSE,cause = 1)

B <- 1
cores <- 1
name <- names(setup.Figure2)[[2]]
set <- setup.Figure2[[name]]

source("simulation-functions.R")
set.seed(8)
run.simulation1(model=model2,setup=set,N=200,do.known=TRUE,timeinterest=5,B=B,cores=cores,true.ate=0,coverage=TRUE)

## near positivity violations
intercept(model2,"A") <- 3
intercept(model2,"A") <- 0
d <- sim(model2,10000,p=setup.Figure2[[1]]$parameters)
mean(d$A)
f <- glm(A~W1+W2+W3+W4+W5+W6+W7+W8+W9+W10+W11+W12,data=d,family=binomial)
publish(f,intercept=1L)
min(predictRisk(f,newdata=d))




do.call("plot.model2",list(parms=setup.Figure2[[1]]$parameters,N=20000,timeinterest=10))

parms <- setup.Figure2[[1]]$parameters
parms["BETA"] <- -1
do.call("plot.model2",list(parms=parms,N=10000,timeinterest=10))


do.call("plot.model2",list(parms=setup.Figure2[[2]]$parameters,N=20000,timeinterest=10))
do.call("plot.model2",list(parms=setup.Figure2[[3]]$parameters,N=20000,timeinterest=10))
do.call("plot.model2",list(parms=setup.Figure2[[4]]$parameters,N=20000,timeinterest=10))
do.call("plot.model2",list(parms=setup.Figure2[[5]]$parameters,N=20000,timeinterest=10))
do.call("plot.model2",list(parms=setup.Figure5$parameters,N=20000,timeinterest=10))

## testing 
parms <- setup.Figure2[["mispecified comprisk model"]]$parameters
## parms["omega1"] <- -.7
d1 <- sim(model2,p=parms,n=1500)
setDT(d1)
timeinterest <- 10
x <- ateRobust(data = d1,
          type = "competing.risks",
          formula.event = list(formula(paste0("Hist(time, event)~A+","W1+W2+W3+W4+W5+W6+W7+W8+W9+W10+W11+W12")),
                               formula(paste0("Hist(time, event)~A+","W2+W3+W4+W5+W7+W8+W9+W10+W11"))),
                               ## formula(paste0("Hist(time, event)~A+","W4"))),
          formula.censor = formula(paste0("Surv(time, event==0)~A+","W1+W2+W3+W4+W5+W6+W7+W8+W9+W10+W11+W12")),
          formula.treatment = formula(paste0("A~","W1+W2+W3+W4+W5+W6+W7+W8+W9+W10+W11+W12")),
          times = 1,
          product.limit = FALSE,
          cause = 1)
x$ate.se

xx <- survtmle(ftime=d1[["time"]],ftype=d1[["event"]],trt=d1[["A"]],adjustVars=data.frame(d1[,c("W1","W2","W3","W4","W5","W6","W7","W8","W9","W10","W11","W12"),with=FALSE]),glm.trt="W1+W2+W3+W4+W5+W7+W8+W9+W10",glm.ftime="trt+W1+W2+W3+W4+W5+W7+W8+W9+W10",glm.ctime="trt+W1+W2+W3+W4+W5+W7+W8+W9+W10",method="mean",t0=timeinterest,ftypeOfInterest=1,returnModels=TRUE)
sd(apply(xx$ic,1,diff))/sqrt(400)
sqrt(xx$var[1,1]+xx$var[2,2]-2*xx$var[2,1])
sqrt(sum(diag(xx$var))-2*xx$var[2,1])


######################################################################
### simulation-tests.R ends here
