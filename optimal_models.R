library(boot)

## problem setting: decision making agent has [horizon] opportunities to either
## approach the prospect and receive an outcome, or avoid the prospect. The
## prospect is known to produce outcomes of +1 with probability either [prob] or
## 1-[prob], and to produce outcomes of -1 otherwise. The two prospect types
## have equal prior probability.

## This function calculates relative value of approaching and
## avoiding for each horizon and information set, using dynamic programming.
## Method starts at final trial and solves backward. horizon: total number of
## choices prob:
finite_values <- function (horizon=32, prob=.67) {
  h <- horizon
  values_ap <- matrix(0, h, 2*(h-1) + 1)
  values_av <- matrix(0, h, 2*(h-1) + 1)
  values <- matrix(0, h, 2*(h-1) + 1)
  VoI <- matrix(0, h, 2*(h-1) + 1)
  Vdiff <- matrix(0, h, 2*(h-1) + 1)
  diff <- -(h-1):(h-1)
  p_good <- prob^diff / (prob^diff + (1-prob)^diff)
  p_bad <- 1 - p_good
  for (i in 1:ncol(values)) {
    values_ap[1, i] <- p_good[i] * (prob * 1 + (1-prob) * -1) +
      p_bad[i] * (prob * -1 + (1-prob) * 1)
    values[1, i] <- max(values_ap[1,i], 0)
    Vdiff[1,i] <- values_ap[1,i] - values_av[1,i]
  }
  for (j in 2:nrow(values)) {
    for (i in j:(ncol(values)-j+1)) {
      values_ap[j, i] <- p_good[i] * (prob * (1 + values[j-1, i+1]) +
                                      (1-prob) * (-1 + values[j-1, i-1])) +
        p_bad[i] * (prob * (-1 + values[j-1, i-1]) +
                    ((1-prob) * (1 + values[j-1, i+1])))
      values_av[j,i] <- values[j-1, i]
      values[j,i] <- max(values_av[j,i], values_ap[j,i])
      VoI[j,i] <- values_ap[j,i] - values_av[j,i] - values_ap[1, i]
      Vdiff[j,i] <- values_ap[j,i] - values_av[j,i]
    }
  }
  ind <- matrix(c( 1 + (0:(length(VoI)-1))%%nrow(VoI),1 + (0:(length(VoI)-1))%/%nrow(VoI)), ncol=2)
  VoIdf <- as.data.frame(cbind(ind, VoI[ind], Vdiff[ind]))
  names(VoIdf) <- c("h", "diff", "voi", "vdiff")
  VoIdf$diff <- VoIdf$diff - max(VoIdf$diff) %/% 2 - 1
  list(VoIdf=VoIdf, values_ap=values_ap, values_av=values_av)
}

## problem setting: decision making agent has an infinite number of
## opportunities to either approach the prospect and receive an outcome, or
## avoid the prospect. Future trials are discounted at a rate of [discount]. The
## prospect is known to produce outcomes of +1 with probability either [prob] or
## 1-[prob], and to produce outcomes of -1 otherwise. The two prospect types
## have equal prior probability.

## This function calculates relative value of approaching and avoiding for each
## information set, using the "value iteration" dynamic programming method.
## additional parameters:
## maxdiff=maximum difference between # of positive and negative outcomes considered
## tol=tolerence for convergence of value during value iteration
infinite_values <- function (maxdiff=40, discount=.98, prob=.67, tol=0.000001) {
  n<- maxdiff*2 + 1
  trans_ap <- matrix(0, n, n)
  trans_av <- diag(n)
  rew_ap <- matrix(0, n, n)
  rew_av <- matrix(0, n, n)
  rew_ap[1,1] <- -(prob - (1-prob))
  rew_ap[n,n] <- (prob - (1-prob))
  trans_ap[1,1] <- 1
  trans_ap[n,n] <- 1
  diff <- -maxdiff:maxdiff
  p_good <- prob^diff / (prob^diff + (1-prob)^diff)
  p_bad <- 1 - p_good
  for (i in 2:(n-1)) {
    for (j in 1:n) {
      if (j == i-1) {
        rew_ap[i, j] <- -1
        trans_ap[i, j] <- p_bad[i] * prob + p_good[i] * (1-prob)
      }
      if (j == i+1) {
        rew_ap[i, j] <- 1
        trans_ap[i, j] <- p_good[i] * prob + p_bad[i] * (1-prob)
      }
    }
  }
  values <- rep(0, maxdiff*2 + 1)
  delta <- tol+1
  while ( delta > tol ) {
    old_values <- values
    disc_values <- values * discount
    values_ap <- sweep(rew_ap, MARGIN=2, disc_values, `+`)
    values_ap <- trans_ap * values_ap
    values_ap <- apply(values_ap, MARGIN=1, sum)
    values_av <- sweep(rew_av, MARGIN=2, disc_values, `+`)
    values_av <- trans_av * values_av
    values_av <- apply(values_av, MARGIN=1, sum)
    values <- pmax(values_av, values_ap)
    delta <- max(abs(values - old_values))
  }
  reward_ap <- p_good * (prob - (1-prob)) + p_bad * ((1-prob) - prob)
  voi <- values_ap - reward_ap - values_av
  list(VoIdf=data.frame(diff=diff, values=values, voi=voi),
       values_ap=values_ap,
       values_av=values_av)
}

## given a global discount rate [gamma], calculate the effective discount rate
## for a prospect with frequency of encounter [f].
effective_discount <- function(gamma, f) {
  sum(sapply(1:100, function (x) {f*(1-f)^(x-1)*gamma^x}))
}

## simulate optimal model [nsims] times in a finite-horizon setting with
## [horizon] trials where the prospect's true probability of a positive outcome
## is [p_healthy].
finite_contingent_simulator <- function(nsims=10000, horizon=32, p_healthy=.33) {
  optim <- finite_values()
  values_diff=optim$values_ap - optim$values_av
  diff0 <- ncol(values_diff) %/% 2 + 1
  stimuli <- matrix(rbinom(horizon*nsims, 1, p_healthy), nrow=horizon)
  responses <- matrix(NA, nrow=horizon, ncol=nsims)
  for (s in 1:nsims) {
    value <- stimuli[,s]
    optresp <- rep(NA, horizon)
    diff <- diff0
    for (t in 1:horizon) {
      if (values_diff[horizon-t+1, diff] >  0) {
        if (value[t] == 1) {
          ## a <- a + 1
          diff <-  diff+1
        } else {
          ## b <- b + 1
          diff <-  diff-1
        }
        optresp[t] <- 1
      } else if (values_diff[horizon-t+1, diff] == 0) {
        optresp[t] <- as.numeric(runif(1) > .5)
      } else {
        optresp[t] <- 0
      }
    }
    responses[,s] <- optresp
  }
  meanresp <- rowMeans(responses)
  firstlasttrials <- rep(FALSE, horizon)
  firstlasttrials[1] <- TRUE
  firstlasttrials[length(firstlasttrials)] <- TRUE
  data.frame(type="contingent", horizon=horizon, freq=NA, trial=1:horizon, response=meanresp,
             firstlasttrials=firstlasttrials, phealthy=p_healthy)
}

## simulate optimal model [nsims] times in a infinite-horizon setting with base
## discount rate [base_discount] and prospect frequency [freq] for [ntrials]
## trials, where the prospect's true probability of a positive outcome is
## [p_healthy].
infinite_contingent_simulator <- function(nsims=10000, freq=.4, base_discount=.99, ntrials=32, p_healthy=.33) {
  optim<- infinite_values(discount=effective_discount(base_discount, freq))
  values_diff <- optim$values_ap - optim$values_av
  diff0 <- length(values_diff) %/% 2 + 1
  stimuli <- matrix(rbinom(ntrials*nsims, 1, p_healthy), nrow=ntrials)
  responses <- matrix(NA, nrow=ntrials, ncol=nsims)
  for (s in 1:nsims) {
    value <- stimuli[,s]
    optresp <- rep(NA, ntrials)
    ## a <- 1
    ## b <- 1
    diff <- diff0
    for (t in 1:ntrials) {
      if (values_diff[diff] >  0) {
        if (value[t] == 1) {
          diff <- diff + 1
        } else {
          diff <- diff - 1
        }
        optresp[t] <- 1
      } else if (values_diff[diff] == 0) {
        optresp[t] <- as.numeric(runif(1) > .5)
      } else {
        optresp[t] <- 0
      }
    }
    responses[,s] <- optresp
  }
  meanresp <- rowMeans(responses)
  firstlasttrials <- rep(FALSE, ntrials)
  firstlasttrials[1] <- TRUE
  firstlasttrials[length(firstlasttrials)] <- TRUE
  data.frame(type="contingent", horizon=NA, freq=freq, trial=1:ntrials, response=meanresp,
             firstlasttrials=firstlasttrials, phealthy=p_healthy)
}

## simulate optimal model [nsims] times for [ntrials] with full information,
## where the prospect's value is revealed regardless of the agent's choice. The
## prospect's true probability of a postivie outcoem is [p_healthy]. As the
## optimal strategy is simply myopic reward-maximization, this function works
## for finite and infinite horizons.
full_info_simulator <- function(nsims=10000, ntrials=32, p_healthy=.33) {
  stimuli <- matrix(rbinom(ntrials*nsims, 1, p_healthy), nrow=ntrials)
  responses <- matrix(NA, nrow=ntrials, ncol=nsims)
  for (s in 1:nsims) {
    value <- stimuli[,s]
    optresp <- rep(NA, ntrials)
    diff <- 0
    for (t in 1:ntrials) {
      if (diff > 0) {
        optresp[t] <- 1
      } else if (diff == 0) {
        optresp[t] <- as.numeric(runif(1) > .5)
      }else {
        optresp[t] <- 0
      }
      if (value[t] == 1) {
        diff <- diff + 1
      } else {
        diff <- diff - 1
      }
    }
    responses[,s] <- optresp
  }
  meanresp <- rowMeans(responses)
  firstlasttrials <- rep(FALSE, ntrials)
  firstlasttrials[1] <- TRUE
  firstlasttrials[length(firstlasttrials)] <- TRUE
  data.frame(type="full", horizon=NA, freq=NA, trial=1:ntrials, response=meanresp,
             firstlasttrials=firstlasttrials, phealthy=p_healthy)
}
