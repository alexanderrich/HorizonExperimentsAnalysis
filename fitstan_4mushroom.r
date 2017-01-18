library("rstan");
library("dplyr");
rstan_options(auto_write = TRUE);
options(mc.cores = parallel::detectCores());

data4m <- read.csv('data_exp2_forsupplement.csv');

## fit stan regression model separately to each experiment/condition
for (i in 0:2) {
  dataCondition <- dataCondition %>% filter(condition == i)
  dataCondition$ll <- as.integer(factor(dataCondition$uniqueid))
  N <- nrow(dataCondition)
  L <- max(dataCondition$ll)
  dataCondition$patchType <- 0
  dataCondition[dataCondition$freq > 0.2, ]$patchType = 1
  patchHorizon <- dataCondition$patchType
  patchTrialNum <- log2(dataCondition$catTrial)
  patchTrialNum <- patchTrialNum - mean(patchTrialNum)
  patchTrialNum <- patchTrialNum
  value <- dataCondition$value
  expectation <- dataCondition$expectation / .333
  y <- dataCondition$response
  ll <- dataCondition$ll
  N <- nrow(dataCondition)
  L <- max(dataCondition$ll)
  Ntrials <- max(dataCondition$catTrial)
  condition <- dataCondition$condition[1]
  standata = c("N", "L", "Ntrials", "y", "ll", "expectation", "patchHorizon", "patchTrialNum", "condition", "value")


  fit <- stan(file="stan_models/regression.stan", data=standata, iter=1000)
  save(fit, file=paste0(c("model_fits/regression_4m_cond", i, ".RData")))
  extracted=extract(fit, permuted=FALSE, pars=c("bias_mean", "bias_sd",
                                                "expectation_mean", "expectation_sd",
                                                "horizon_mean", "horizon_sd",
                                                "trial_mean", "trial_sd",
                                                "interaction_mean", "interaction_sd"))
  print(monitor(extracted, digits_summary=2))
}
