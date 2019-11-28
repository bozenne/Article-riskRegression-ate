##                                                                           ##
##                                 WARNING                                   ##
##                                                                           ##
## BECAUSE WE ARE NOT ALLOWED TO SHARE THE DATASET USED IN THE ILLUSTRATION
## THIS SCRIPT DOES NOT REPRODUCE THE FIGURE 4 OF THE ARTICLE
## WE INSTEAD USED SIMULATED DATA AND A SIMILAR ANALYSIS TO RE-PRODUCE FIGURE 4

## source("figure4.R")

library(data.table)
library(riskRegression)
library(ggplot2)
source("simulation-functions.R")

## * Import results
set.seed(10)
dt <- simulate.data() ## n=13946
##dt <- dt[sample.int(.N,1000)] ## subset to make the execution of ate faster.

##table(dt$event)

## * Compute ATE
method <- "read" ## "IF" "bootstrap" "read"
## table(dt$event)
if(method == "IF"){
    tps <- system.time(
        e.ate <- ate(event = Hist(time,event) ~ treatment+male+age+heart.failure+hypertension+prior.stroke+prior.bleeding,
                     treatment = treatment ~ age+male+hypertension,
                     censor = Surv(time,event == 0) ~ heart.failure+hypertension+prior.stroke+prior.bleeding,
                     estimator = c("Gformula","AIPTW"),
                     times = seq(1, 23, by = 1), cause = 1,
                     data = dt)
    )
    ## print(tps)
    ## timing: n= 1000 --> 71 = 1 min
    ## timing: n= 5000 --> 371s = 6 min
    ## timing: n=10000 --> 1521s = 25 min
    ## timing: n=13946 --> 4256 = 1h10 min
	
	## saveRDS(e.ate, file = file.path("Results","figure4-res-ate.rds"))
}else if(method == "bootstap"){
    system.time(
        e.ate <- ate(event = Hist(time,event) ~ treatment+male+age+heart.failure+hypertension+prior.stroke+prior.bleeding,
                     treatment = treatment ~ age+male+hypertension,
                     censor = Surv(time,event == 0) ~ heart.failure+hypertension+prior.stroke+prior.bleeding,
                     estimator = c("Gformula","AIPTW"),
                     ## estimator = c("Gformula"),
                     times = seq(1, 23, by = 1), cause = 1,
                     B = 1e4,
                     data = dt)
    )
}else if(method == "read"){
    e.ate <- readRDS(file.path("Results","figure4-res-ate.rds"))
}
## summary(e.ate, estimator = "AIPTW"),
## summary(e.ate, estimator = "Gformula")

## * figure 4
## 
dt.ate <- as.data.table(e.ate)
dt.ate[,estimator := factor(estimator, levels = c("Gformula","AIPTW"), labels = c("G-formula","AIPTW,AIPCW"))]

figure4 <- ggplot(dt.ate[type == "ate"], aes(x = time, group = level))
figure4 <- figure4 + geom_ribbon(aes(ymin = lower, ymax = upper, fill = level), alpha = 0.3)
figure4 <- figure4 + geom_line(aes(y = value, color = level), size = 1)
figure4 <- figure4 + scale_colour_manual(values = c("black","grey66"))
figure4 <- figure4 + scale_fill_manual(values = c("black","grey66"))
figure4 <- figure4 + theme_minimal() + theme(legend.spacing.x = unit(0.2, 'cm'), legend.position="top")
figure4 <- figure4 + scale_x_continuous(breaks=seq(0,21,3)) + scale_y_continuous(labels = scales::percent)
figure4 <- figure4 + facet_grid(~estimator)
figure4 <- figure4 + xlab("Months since intiation of treatment") + ylab("Absolute risk of major bleeding (%)") + labs(colour="treatment", fill = "treatment")

ggsave(figure4, filename = file.path("Figures","figure4-illustration.pdf"))
