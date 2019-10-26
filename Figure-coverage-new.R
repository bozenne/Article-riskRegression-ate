if (system("echo $USER",intern=TRUE)=="tag"){
    setwd("~/research/Methods/ATE/dropbox/manuscript/")
}else{
    setwd("c:/Users/hpl802/Documents/Projects/Gformula-modelChecking/manuscript")
}
source("../simulation/simulation-models.R")
source("../simulation/simulation-functions.R")

## * Import results
path.simulation <- "../simulation/simulation-results"

estimators <- c("G-formula"= "Gformula", "IPCW,IPTW" = "IPTW", "AIPCW,AIPTW (robust CI)" = "AIPTW", "AIPCW,AIPTW (non-robust CI)" = "AIPTW0")

file.correct <- grep("^results-FigureSE-ReRevision-correct-N-",list.files(path.simulation), value = TRUE)
file.outcome <- grep("^results-FigureSE-ReRevision-misout-N-",list.files(path.simulation), value = TRUE)
file.treatment <- grep("^results-FigureSE-ReRevision-mistreat-N-",list.files(path.simulation), value = TRUE)
file.censoring <- grep("^results-FigureSE-ReRevision-miscens-N-",list.files(path.simulation), value = TRUE)

ls.dt.figure3 <- list()
ls.dt.figure3$correct <- do.call(rbind,lapply(file.correct, function(iFile){ ## iFile <- file.correct[1]
    data.table(readRDS(file.path(path.simulation,iFile)), scenario = "All models correct")
}))
ls.dt.figure3$treatment <- do.call(rbind,lapply(file.treatment, function(iFile){ ## iFile <- file.correct[1]
    data.table(readRDS(file.path(path.simulation,iFile)), scenario = "Misspecified treatment model")
}))
ls.dt.figure3$outcome <- do.call(rbind,lapply(file.outcome, function(iFile){ ## iFile <- file.correct[1]
    data.table(readRDS(file.path(path.simulation,iFile)), scenario = "Misspecified outcome model")
}))
ls.dt.figure3$censoring <- do.call(rbind,lapply(file.censoring, function(iFile){ ## iFile <- file.correct[1]
    data.table(readRDS(file.path(path.simulation,iFile)), scenario = "Misspecified censoring model")
}))

dtW.figure3 <- do.call(rbind, ls.dt.figure3)
dtW.figure3[,AIPTW0 := AIPTW]
setnames(dtW.figure3, old = "se.known.AIPTW.se", new = "AIPTW0.se")
setnames(dtW.figure3, old = "cov.known.coverage.AIPTW", new = "coverage.AIPTW0")

## * Process results
dtL.figure3 <- melt(dtW.figure3[!is.na(KM.naive)], id.vars = c("N","scenario","true"),
                    measure.vars = list(estimators,
                                        paste0(estimators,".se"),
                                        paste0("coverage.",estimators)),
                    value.name = c("estimate","se","coverage"), variable.name = "estimator")
dtL.figure3[, estimator := factor(estimator, levels = 1:length(estimators), labels = names(estimators))]
dtL.figure3[, N := factor(N), ]
dtL.figure3[abs(estimate) > 1, c("estimate","se","coverage") := NA]
dtL.figure3[ , scenario := factor(scenario, levels = unique(scenario))]
setkeyv(dtL.figure3, c("scenario","estimator","N"))
dtLS.figure3 <- dtL.figure3[,.(rep = .N, coverage = mean(coverage,na.rm = TRUE)),
                            by = c("N","estimator","scenario")]
dtLS.figure3[scenario == "Misspecified outcome model" & estimator == "G-formula", coverage := NA]
dtLS.figure3[scenario == "Misspecified treatment model" & estimator == "IPCW,IPTW", coverage := NA]
dtLS.figure3[scenario == "Misspecified censoring model" & estimator == "IPCW,IPTW", coverage := NA]

require(ggplot2)
Fig3.cov <- ggplot(dtLS.figure3[estimator %in% c("AIPCW,AIPTW (robust CI)","AIPCW,AIPTW (non-robust CI)")],
                   aes(x=N, y=coverage, color=estimator, group=estimator, shape = estimator))
Fig3.cov <- Fig3.cov + geom_abline(intercept=0.95,slope=0,color=1)
Fig3.cov <- Fig3.cov + facet_wrap(~scenario)
Fig3.cov <- Fig3.cov + geom_line(size = 1.25) + geom_point(size = 2)
Fig3.cov <- Fig3.cov + xlab("Sample size (N)") + ylab("Coverage") ##+ coord_cartesian(ylim = c(0.9,1))
Fig3.cov <- Fig3.cov + theme_minimal() + theme(legend.title=element_blank(), legend.spacing.x = unit(0.2, 'cm'), legend.position="bottom")
Fig3.cov <- Fig3.cov + scale_shape_manual(breaks = names(estimators)[3:4],
                                          values = c(15,17))
Fig3.cov <- Fig3.cov + scale_color_manual(breaks = names(estimators)[3:4],
                                          values = c("black","gray66"))

## getwd()
ggsave(Fig3.cov + coord_cartesian(ylim = c(0.8,1)),
       filename = "./reresubmission/figure-simulation-N.pdf")

