library(data.table)
library(ggplot2)

## * Import results
file.correct <- grep("^figure3-correct-N-",list.files("./Results"), value = TRUE)
file.treatment <- grep("^figure3-mistreat-N-",list.files("./Results"), value = TRUE)
file.outcome <- grep("^figure3-misout-N-",list.files("./Results"), value = TRUE)
file.censoring <- grep("^figure3-miscens-N-",list.files("./Results"), value = TRUE)

ls.dt.fig3 <- list()
ls.dt.fig3$correct <- do.call(rbind,lapply(file.correct, function(iFile){ ## iFile <- file.correct[1]
    data.table(readRDS(file.path("./Results",iFile)), scenario = "All models correct")
}))
ls.dt.fig3$treatment <- do.call(rbind,lapply(file.treatment, function(iFile){ ## iFile <- file.treatment[1]
    data.table(readRDS(file.path("./Results",iFile)), scenario = "Misspecified treatment model")
}))
ls.dt.fig3$outcome <- do.call(rbind,lapply(file.outcome, function(iFile){ ## iFile <- file.outcome[1]
    data.table(readRDS(file.path("./Results",iFile)), scenario = "Misspecified outcome model")
}))
ls.dt.fig3$censoring <- do.call(rbind,lapply(file.censoring, function(iFile){ ## iFile <- file.censoring[1]
    data.table(readRDS(file.path("./Results",iFile)), scenario = "Misspecified censoring model")
}))

## * Process results
vec.estimator <- c("G-formula" = "Gformula",
                   "IPCW,IPTW" = "IPTW",
                   "AIPCW,AIPTW (robust CI)" = "AIPTW",
                   "AIPCW,AIPTW (non-robust CI)" = "AIPTW0")

dtW.fig3 <- do.call(rbind, ls.dt.fig3)
dtW.fig3[,AIPTW0 := AIPTW]
setnames(dtW.fig3, old = "se.known.AIPTW.se", new = "AIPTW0.se")
setnames(dtW.fig3, old = "cov.known.coverage.AIPTW", new = "coverage.AIPTW0")

dtL.fig3 <- melt(dtW.fig3[!is.na(KM.naive)], id.vars = c("N","scenario","true"),
                 measure.vars = list(vec.estimator,
                                     paste0(vec.estimator,".se"),
                                     paste0("coverage.",vec.estimator)),
                 value.name = c("estimate","se","coverage"), variable.name = "estimator")
dtL.fig3[, estimator := factor(estimator, levels = 1:length(vec.estimator), labels = names(vec.estimator))]
dtL.fig3[, N := factor(N), ]
dtL.fig3[abs(estimate) > 1, c("estimate","se","coverage") := NA]
dtL.fig3[ , scenario := factor(scenario, levels = unique(scenario))]
setkeyv(dtL.fig3, c("scenario","estimator","N"))

dtLS.fig3 <- dtL.fig3[,.(rep = .N, coverage = mean(coverage,na.rm = TRUE)),
                      by = c("N","estimator","scenario")]
dtLS.fig3[scenario == "Misspecified outcome model" & estimator == "G-formula", coverage := NA]
dtLS.fig3[scenario == "Misspecified treatment model" & estimator == "IPCW,IPTW", coverage := NA]
dtLS.fig3[scenario == "Misspecified censoring model" & estimator == "IPCW,IPTW", coverage := NA]

## * Create figure
figure3 <- ggplot(dtLS.fig3[estimator %in% c("AIPCW,AIPTW (robust CI)","AIPCW,AIPTW (non-robust CI)")],
                  aes(x=N, y=coverage, color=estimator, group=estimator, shape = estimator))
figure3 <- figure3 + geom_abline(intercept=0.95,slope=0,color=1)
figure3 <- figure3 + facet_wrap(~scenario)
figure3 <- figure3 + geom_line(size = 1.25) + geom_point(size = 2)
figure3 <- figure3 + xlab("Sample size (N)") + ylab("Coverage") ##+ coord_cartesian(ylim = c(0.9,1))
figure3 <- figure3 + theme_minimal() + theme(legend.title=element_blank(), legend.spacing.x = unit(0.2, 'cm'), legend.position="bottom")
figure3 <- figure3 + scale_shape_manual(values = c(15,17))
figure3 <- figure3 + scale_color_manual(values = c("black","gray66"))
figure3 <- figure3 + coord_cartesian(ylim = c(0.8,1))

## * Export figure
ggsave(figure3, filename = "figures/figure3-simulation-coverage.pdf", width = 8, height = 7)

