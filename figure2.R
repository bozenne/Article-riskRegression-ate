library(data.table)
library(ggplot2)

## * Import results
ls.dt.fig2 <- list(
    "All models correct" = readRDS(file="./Results/figure2-All-models-correct.rds"),
    "Misspecified treatment model" = readRDS(file="./Results/figure2-Misspecified-treatment-model.rds"),
    "Misspecified outcome model" = readRDS(file="./Results/figure2-Misspecified-outcome-model.rds"),
    "Misspecified censoring model"  = readRDS(file="./Results/figure2-Misspecified-censoring-model.rds")
)

## * Process results
vec.estimator <- c(
    "AJ" = "KM.naive",
    "G-formula" = "Gformula",
    "IPCW,IPTW" = "IPTW",
    "AIPCW,AIPTW" = "AIPTW"
)
vec.scenario <- names(ls.dt.fig2)
keep.correlation <- "corr.1"

dt.fig2 <- do.call(rbind,lapply(vec.scenario, function(iScenario){ ## iScenario <- names(ls.dt.fig2)[1]
    iDT.W <- as.data.table(ls.dt.fig2[[iScenario]][[keep.correlation]])
    iDT.L <- melt(iDT.W, id=c("cens.percent","true.ate"),
                 measure.vars = c("KM.naive","Gformula","IPTW","AIPTW"),
                 value.name = "bias",
                 variable.name = "estimator")
    iDT.L[,scenario := iScenario]
    return(iDT.L)
}))
dt.fig2[, estimator := factor(estimator, levels = vec.estimator, labels = names(vec.estimator))]
dt.fig2[, scenario := factor(scenario, levels = unique(scenario))]

## * Create figure
figure2 <- ggplot(dt.fig2, aes(x = estimator, y = bias, fill = estimator))
figure2 <- figure2 + geom_boxplot(alpha=0.3)
figure2 <- figure2 + geom_abline(intercept=0, slope=0, color=1)
figure2 <- figure2 + facet_wrap(~ scenario)
figure2 <- figure2 + xlab("") + ylab("Average treatment effect") + ylim(c(-.5,.5))
figure2 <- figure2 + theme_minimal() + theme(legend.position="none")
figure2 <- figure2 + scale_fill_manual(values = c("gray1","gray25","gray50","gray90"))

## * Export figure
ggsave(figure2, filename = "figures/figure2-simulation-bias.pdf", width = 8, height = 7)
