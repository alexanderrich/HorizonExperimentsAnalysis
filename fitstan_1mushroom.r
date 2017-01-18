library("rstan");
library("dplyr");
rstan_options(auto_write = TRUE);
options(mc.cores = parallel::detectCores());

data1m<- read.csv('data_exp1_forsupplement.csv');

## fit stan regression model separately to each experiment/condition
for (i in 0:2) {
  dataCondition <- dataCondition %>% filter(condition == i)
  dataCondition$ll <- as.integer(factor(dataCondition$uniqueid))
  N <- nrow(dataCondition)
  L <- max(dataCondition$ll)
  patchHorizon <- log2(dataCondition$patchTrialLeft) / 5
  patchLength <- dataCondition$patchLength
  patchTrialNum <- log2(dataCondition$patchTrial)
  patchTrialNum <- patchTrialNum - mean(patchTrialNum)
  value <- dataCondition$value
  expectation <- dataCondition$expectation / .333
  condition <- dataCondition$condition
  y <- dataCondition$response
  ll <- dataCondition$ll
  N <- nrow(dataCondition)
  L <- max(dataCondition$ll)
  Ntrials <- 32
  condition <- dataCondition$condition[1]
  standata = c("N", "L", "Ntrials", "y", "ll", "expectation", "patchHorizon", "patchTrialNum", "condition", "value")


  fit <- stan(file="stan_models/regression.stan", data=standata, iter=100)
  save(fit, file=paste0(c("model_fits/regression_1m_cond", i, ".RData")))
  extracted=extract(fit, permuted=FALSE, pars=c("bias_mean", "bias_sd",
                                                "expectation_mean", "expectation_sd",
                                                "horizon_mean", "horizon_sd",
                                                "trial_mean", "trial_sd",
                                                "interaction_mean", "interaction_sd"))
  print(monitor(extracted, digits_summary=2))
}
