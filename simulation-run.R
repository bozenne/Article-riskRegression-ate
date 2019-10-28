### simulation-run.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Aug 17 2018 (11:53) 
## Version: 
## Last-Updated: okt 28 2019 (13:15) 
##           By: Brice Ozenne
##     Update #: 258
#----------------------------------------------------------------------
## 
### Commentary: 
## Running the empirical studies for the competing-risk-ATE manuscript
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
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
source("simulation-settings.R")
source("simulation-functions.R")

## * Simulation for figure 2: bias under misspecified models (fix sample size, fix censoring %)
n.sim <- 1000 ## number of simulations
n.obs <- 500 ## sample size in each simulation
tau <- 10 ## time at which the ATE is computed
cores <- min(detectCores()-1,32) ## number of CPUs used to run the simulations

set.seed(1)
for(iSetup in names(setup.Figure2)){ ## iSetup <- names(setup.Figure2)[1]
    print(iSetup)
    
    iRes <- run.simulation(model = model2,
                           setup = setup.Figure2[[iSetup]],
                           N = n.obs,
                           timeinterest = tau,
                           B = n.sim,
                           cores = cores,
                           true.ate=0
                           do.known = TRUE,
                           coverage = FALSE)
    
    xclass <- class(iRes)
    iRes <- cbind(iRes, N=n.obs, scenario = iSetup, time = tau)
    class(iRes) <- xclass
    
    saveRDS(iRes,file=paste0("Results/figure2-",gsub(" ","-",iSetup),".rds"))
}



## * Simulation for figure 3: coverage considering increasing sample size
n.sim <- 1000 ## number of simulations
vec.n.obs <- c(100,200,500,750,1000) ## vector of the sample sizes
cores <- min(detectCores()-1,64) ## number of CPUs used to run the simulations
tau <- 8 ## time at which the ATE is computed

set.seed(1)
for(iN.obs in 1:vec.n.obs){ ## iN.obs <- 1
    cat("sample.size: ",vec.n.obs[iN.obs],"\n")

    for(iSetup in names(setup.Figure2)){ ## iSetup <- names(setup.Figure2)[1]
        print(iSetup)
    
        ## set.seed(19)
        iRes <- run.simulation(model = model2,
                               setup = setup.Figure2[[iSetup]],
                               N = vec.n.obs[iN.obs],
                               timeinterest = tau,
                               B = n.sim,
                               cores = cores,
                               true.ate = 0,
                               do.known = TRUE,
                               coverage = TRUE)
        xclass <- class(iRes)
        iRes <- cbind(iRes, N = vec.n.obs[iN.obs], scenario = iSetup, time = tau)
        class(iRes) <- xclass

        saveRDS(iRes, file = paste0("Results/figue2-",gsub(" ","-",iSetup),"-",vec.n.obs[iN.obs],".rds"))
    }
    cat("\n")

}


######################################################################
### simulation-run.R ends here
