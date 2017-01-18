source("analysis_functions.r")

data1m <- preprocess1m()


data1m <- data1m %>% filter(condition == 1)
data1m$ll <- as.integer(factor(data1m$uniqueid))
N <- nrow(data1m)
L <- max(data1m$ll)
## patchHorizon <- log2(data1m$patchTrialLeft) / 5
patchHorizonRaw <- data1m$patchTrialLeft
patchLength <- data1m$patchLength
patchTrialNum <- log2(data1m$catTrial)
patchTrialNum <- patchTrialNum - mean(patchTrialNum)
patchTrialNum <- patchTrialNum
value <- data1m$value
expectation <- data1m$expectation / .333
condition <- data1m$condition
y <- data1m$response
ll <- data1m$ll
N <- nrow(data1m)
L <- max(data1m$ll)
Ntrials <- 32
condition <- data1m$condition[1]
standata = c("N", "L", "Ntrials", "y", "ll", "expectation", "patchHorizonRaw", "patchTrialNum", "condition", "value")


fit <- stan(file="stan_models/regression_freehorizon.stan", data=standata, iter=2000)
save(fit, file="model_fits/regression_1m_freehorizon_cond1.RData")
extracted=extract(fit, permuted=FALSE, pars=c("horizon",
                                              "bias_mean", "bias_sd",
                                              "expectation_mean", "expectation_sd",
                                              ## "horizon_mean", "horizon_sd",
                                              "horizon_sd",
                                              "trial_mean", "trial_sd",
                                              "interaction_mean", "interaction_sd"))
monitor(extracted, digits_summary=2)
