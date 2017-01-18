data {
  int<lower=0> N; //rows of observations
  int<lower=1> L; //# of subjects
  int<lower=0,upper=1> y[N]; //responses
  int<lower=1, upper=L> ll[N]; //subject producing each response
  int patchHorizonRaw[N]; // horizon of patch
  vector[N] patchTrialNum; // encounter in mushroom species
  vector[N] expectation; // expected reward given past observations
  int value[N]; // 1 for healthy, 0 for poisonous
  int condition; // determines calculations in simulation for generated quantities
}

parameters {
  real bias_mean;
  real<lower=0> bias_sd;

  real expectation_mean;
  real<lower=0> expectation_sd;

  real<lower=0> horizon_sd;

  vector[32] horizon;

  real trial_mean;
  real<lower=0> trial_sd;

  real interaction_mean;
  real<lower=0> interaction_sd;

  vector[L] subj_expectation_raw;
  vector[L] subj_bias_raw;
  vector[L] subj_horizon_raw;
  vector[L] subj_trial_raw;
  vector[L] subj_interaction_raw;
}

transformed parameters {
  vector[L] subj_expectation;
  vector[L] subj_bias;
  vector[L] subj_horizon;
  vector[L] subj_trial;
  vector[L] subj_interaction;

  // for efficiency, we assume random effects have mean 0 and SD 1, then transform them.
  subj_expectation <- subj_expectation_raw * expectation_sd + expectation_mean;
  subj_horizon <- subj_horizon_raw * horizon_sd + 1;
  subj_trial <- subj_trial_raw * trial_sd + trial_mean;
  subj_interaction <- subj_interaction_raw * interaction_sd + interaction_mean;
  subj_bias <- subj_bias_raw * bias_sd + bias_mean;
}

model {
  int subj;
  vector[N] predictor;
  vector[N] p;
  vector[N] b_expectation;
  vector[N] b_horizon;
  vector[N] b_trial;
  vector[N] b_interaction;
  vector[N] b_bias;
  vector[N] patchHorizon;

  horizon ~ normal(0, 10);
  expectation_mean ~ normal(0, 10);
  expectation_sd ~ normal(0, 5);
  horizon_sd ~ normal(0, 5);
  trial_mean ~ normal(0, 5);
  trial_sd ~ normal(0, 5);
  interaction_mean ~ normal(0, 5);
  interaction_sd ~ normal(0, 5);
  bias_mean ~ normal(0, 5);
  bias_sd ~ normal(0, 5);

  subj_expectation_raw ~ normal(0, 1);
  subj_bias_raw ~ normal(0, 1);
  subj_horizon_raw ~ normal(0, 1);
  subj_trial_raw ~ normal(0, 1);
  subj_interaction_raw ~ normal(0, 1);

  for (n in 1:N) {
    subj <- ll[n];
    b_expectation[n] <- subj_expectation[subj];
    b_horizon[n] <- subj_horizon[subj];
    b_trial[n] <- subj_trial[subj];
    b_interaction[n] <- subj_interaction[subj];
    b_bias[n] <- subj_bias[subj];
    patchHorizon[n] <- horizon[patchHorizonRaw[n]];
  }
  predictor <- b_bias + b_expectation .* expectation + b_horizon .* patchHorizon +
    b_trial .* patchTrialNum + b_interaction .* patchHorizon .* patchTrialNum;
  p <- 1 ./ (1 + exp(-predictor));
  y ~ bernoulli(p);
}

// for each sample, generate a set of predictions for each trial given the
// participant's responses up to that trial, and also generate a simulation of
// the model's behavior given the model's past responses
generated quantities {
  vector[N] y_pred;
  vector[N] y_sim;
  {
    int subj;
    vector[N] predictor;
    vector[N] p;
    vector[N] b_expectation;
    vector[N] b_horizon;
    vector[N] b_trial;
    vector[N] b_interaction;
    vector[N] b_bias;
    vector[N] patchHorizon;
    for (n in 1:N) {
      subj <- ll[n];
      b_expectation[n] <- subj_expectation[subj];
      b_horizon[n] <- subj_horizon[subj];
      b_trial[n] <- subj_trial[subj];
      b_interaction[n] <- subj_interaction[subj];
      b_bias[n] <- subj_bias[subj];
      patchHorizon[n] <- horizon[patchHorizonRaw[n]];
    }
    predictor <- b_bias + b_expectation .* expectation + b_horizon .* patchHorizon +
      b_trial .* patchTrialNum + b_interaction .* patchHorizon .* patchTrialNum;
    p <- 1 ./ (1 + exp(-predictor));
    for (n in 1:N) {
      y_pred[n] <- bernoulli_rng(p[n]);
    }
  }
  {
    int subj;
    real predictor;
    real p;
    real b_expectation;
    real b_horizon;
    real b_trial;
    real b_interaction;
    real b_bias;
    real obsgood;
    real obsbad;
    real expectation_sim;
    real probgood;
    real patchHorizon;
    for (n in 1:N) {
      subj <- ll[n];
      if(patchTrialNum[n] == 1) {
        obsgood <- 0;
        obsbad <- 0;
      }
      if (condition != 2) {
        probgood <- pow(2, obsgood-obsbad)/(1+pow(2, obsgood-obsbad));
        expectation_sim <- .33*probgood - .33*(1-probgood);
      } else {
        expectation_sim <- 2 * (obsgood + 1) / (obsgood+obsbad + 2) - 1;
      }
      patchHorizon <- horizon[patchHorizonRaw[n]];
      b_expectation <- subj_expectation[subj];
      b_horizon <- subj_horizon[subj];
      b_trial <- subj_trial[subj];
      b_interaction <- subj_interaction[subj];
      b_bias <- subj_bias[subj];
      predictor <- b_bias + b_expectation .* expectation[n] + b_horizon .* patchHorizon +
        b_trial .* patchTrialNum[n] + b_interaction .* patchHorizon .* patchTrialNum[n];
      p <- 1 ./ (1 + exp(-predictor));
      y_sim[n] <- bernoulli_rng(p);
      if (y_sim[n] == 1) {
        obsgood <- obsgood + value[n];
        obsbad <- obsbad + 1 - value[n];
      } else if (condition == 1) {
        obsgood <- obsgood + value[n];
        obsbad <- obsbad + 1 - value[n];
      }
    }
  }
}
