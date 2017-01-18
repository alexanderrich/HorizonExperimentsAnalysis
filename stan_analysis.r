library("rstan")
library("ggplot2")
library("dplyr")
library("tidyr")

get_model_df <- function(filename, level, experiment, information) {
  load(filename)
  extracted <- extract(fit)
  fitdf <- data.frame(
    bias=extracted$bias_mean,
    expectation=extracted$expectation_mean,
    horizon=extracted$horizon_mean,
    trial=extracted$trial_mean,
    interaction=extracted$interaction_mean) %>%
    gather(parameter, sample) %>%
    group_by(parameter) %>%
    summarize(low=quantile(sample, .025),
              mean=mean(sample),
              high=quantile(sample, .975))
  fitdf$condition <- level
  ypred <- extract(fit)$y_sim
  extracted=extract(fit, permuted=FALSE, pars=c("bias_mean", "bias_sd",
                                                "expectation_mean", "expectation_sd",
                                                "horizon_mean", "horizon_sd",
                                                "trial_mean", "trial_sd",
                                                "interaction_mean", "interaction_sd"))
  summary <- as.data.frame(monitor(extracted, digits_summary=2, print=FALSE))
  summary$Experiment <- experiment
  summary$Information <- information
  summary$level <- level
  summary$param <- c("bias", "bias", "expectation", "expectation",
                     "horizon", "horizon", "trial", "trial", "interaction", "interaction")
  list(summary, ypred)
}

models <- c(get_model_df("model_fits/regression_1m_cond0.RData", 0, "1b", "contingent"),
            get_model_df("model_fits/regression_1m_cond1.RData", 1, "1b", "full"),
            get_model_df("model_fits/regression_1m_cond2.RData", 2, "1a", "contingent"),
            get_model_df("model_fits/regression_4m_cond0.RData", 3, "2b", "contingent"),
            get_model_df("model_fits/regression_4m_cond1.RData", 4, "2b", "full"),
            get_model_df("model_fits/regression_4m_cond2.RData", 5, "2a", "contingent"))

all_summaries <- rbind.data.frame(models[1][[1]], models[2][[1]], models[3][[1]], models[4][[1]], models[5][[1]], models[6][[1]])
all_summaries$low <- all_summaries[,"2.5%"]
all_summaries$high <- all_summaries[,"97.5%"]
all_summaries <- select_(all_summaries, "param", "level", "Experiment", "Information", "mean", "low", "high")
all_summaries$paramtype <- rep(c(0,1))
all_summaries <- all_summaries %>%
  gather(stat, value, mean, low, high) %>%
  unite(paramstat, paramtype, stat) %>%
  spread(paramstat, value)
all_summaries <- all_summaries[,c(3,4,1,7,6,5,10,9,8)]
all_summaries <- arrange(all_summaries, Experiment, Information)


## print posterior intervals for population regression parameters
library("xtable")
xtable(all_summaries)


## plot posterior intervals of population regression parameters
params <- rbind.data.frame(df0, df1, df2, df3, df4, df5) %>%
  mutate(condition=factor(condition, levels=0:5,
                          labels=c("Exp1a", "Exp1b cont.", "Exp1b full",
                                   "Exp2a", "Exp2b cont.", "Exp2b full")))

ggplot(params, aes(x=parameter, y=mean, ymin=low, ymax=high)) +
  geom_hline(yintercept=0, color="gray", size=1) +
  geom_boxplot(stat="identity", aes(ymin=low, ymax=high, lower=low, upper=high, middle=mean, fill=parameter))+
  scale_x_discrete(labels=c("bias", "expected reward", "horizon", "trial", "horizon x trial")) +
  labs(x="", y="coefficient value") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        plot.margin=unit(c(0.2,0.2,0,0), "cm"),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = "bottom") +
  theme(axis.text.x = element_text(angle = 30, hjust=1)) +
  facet_wrap(~condition, nrow=2) +
  guides(fill=FALSE)
ggsave(filename="figures/Fig3.pdf", width=6.5, height=5, useDingbats=F)

## make graphs of model predictions of behavior to compare to people's actual behavior.

add_model_mean_CI_1m <- function(ypred, df) {
  results <- matrix(, nrow=2000, ncol=126)
  for (i in 1:2000) {
    df$ypred <- ypred[i,]
    dftemp <- df %>%
      group_by(condition, healthyMushroom, patchLength, catTrial) %>%
      summarise(prop = mean(ypred))
    results[i,] <- dftemp$prop
  }
  quantiles <- apply(results, 2,  quantile, probs= c(.025,.975), na.rm=T )
  df$ypred <- colMeans(ypred)
  dftemp <- df %>%
    group_by(condition, healthyMushroom, patchLength, catTrial) %>%
    summarise(prop=mean(ypred))
  dftemp$lower <- quantiles[1,]
  dftemp$upper <- quantiles[2,]
  dftemp
}

data1m<- read.csv('data_exp1_forsupplement.csv');
data1m0 <- data1m %>% filter(condition == 0)
data1m0 <- add_model_mean_CI_1m(models[1][[2]], data1m0)
data1m1 <- data1m %>% filter(condition == 1)
data1m1 <- add_model_mean_CI_1m(models[2][[2]], data1m1)
data1m2 <- data1m %>% filter(condition == 2)
data1m2 <- add_model_mean_CI_1m(models[3][[2]], data1m2)
data1mGraph <- rbind(data1m0, data1m1, data1m2)


data1mGraph <- data1mGraph %>%
  ungroup() %>%
  mutate(patchLength=factor(patchLength, c(32, 16, 8, 4, 2, 1))) %>%
  mutate(healthyMushroom=factor(healthyMushroom, levels=c(TRUE, FALSE), labels=c("Good species", "Bad species")),
         condition=factor(condition, levels=0:2,
                          labels=c("Exp1a",
                                   "Exp1b contingent",
                                   "Exp1b full"
                                   ))) %>%
  filter(catTrial < 17 | catTrial == 32) %>%
  mutate(catTrial=ifelse(catTrial == 32, 20, catTrial)) %>%
  mutate(catTrialDodged=catTrial-(as.numeric(patchLength)-2.5)/10)

ggplot(data1mGraph, aes(x=catTrialDodged, y=prop, color=patchLength, shape=patchLength)) +
  geom_hline(yintercept=.5, color="darkgray") +
  geom_line(data=subset(data1mGraph, catTrial<20), position=position_dodge(width=-0.5)) +
  geom_linerange(data=subset(data1mGraph, catTrial<20), aes(ymin=lower, ymax=upper), position=position_dodge(width=-0.5)) +
  geom_line(position=position_dodge(width=-0.5), linetype="21") +
  ylab("p(approach)") +
  xlab("Trial in species") +
  facet_grid(healthyMushroom~condition) +
  theme_minimal() +
  scale_y_continuous(breaks=c(1, .5, 0)) +
  scale_x_continuous(breaks=c(0, 8, 16, 18, 20), labels=c("0", "8", "16", "", "32")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        strip.text = element_text(size=10),
        plot.margin=unit(c(0,0,0,0), "cm"),
        legend.margin=unit(0, "cm"),
        panel.border = element_rect(colour = "black", fill=NA)) +
  labs(color="Patch\nlength", shape="Patch\nlength") +
  geom_pointrange(data=subset(data1mGraph, (catTrial == 1 | catTrial == 20 | catTrial==as.numeric(as.character(patchLength)))),
                  aes(ymin=lower, ymax=upper),
                  position=position_dodge(width=-0.5))+
  guides(fill=FALSE) +
  theme(legend.position="right")
ggsave(filename="figures/FigS1.pdf", width=9, height=4, useDingbats=F)


add_model_mean_CI_4m <- function(ypred, df) {
  dftemp <- df %>%
    group_by(condition, healthyMushroom, highFreq, catTrial) %>%
    summarise(prop=mean(ypred))
  results <- matrix(, nrow=2000, ncol=nrow(dftemp))
  for (i in 1:2000) {
    df$ypred <- ypred[i,]
    dftemp <- df %>%
      group_by(condition, healthyMushroom, highFreq, catTrial) %>%
      summarise(prop = mean(ypred))
    results[i,] <- dftemp$prop
  }
  quantiles <- apply(results, 2,  quantile, probs= c(.025,.975), na.rm=T )
  df$ypred <- colMeans(ypred)
  dftemp <- df %>%
    group_by(condition, healthyMushroom, highFreq, catTrial) %>%
    summarise(prop=mean(ypred))
  dftemp$lower <- quantiles[1,]
  dftemp$upper <- quantiles[2,]
  dftemp
}

data4m<- read.csv('data_exp2_forsupplement.csv');
data4m0 <- data4m %>% filter(condition == 0)
data4m0 <- add_model_mean_CI_4m(models[4][[2]], data4m0)
data4m1 <- data4m %>% filter(condition == 1)
data4m1 <- add_model_mean_CI_4m(models[5][[2]], data4m1)
data4m2 <- data4m %>% filter(condition == 2)
data4m2 <- add_model_mean_CI_4m(models[6][[2]], data4m2)
data4mGraph <- rbind(data4m0, data4m1, data4m2)

data4mGraph <- data4mGraph %>%
  ungroup() %>%
  mutate(freqlabel = factor(highFreq, levels=c(TRUE, FALSE),
                            labels=c("high\n(f=0.4)", "low\n(f=0.1)"))) %>%
  mutate(healthyMushroom=factor(healthyMushroom, levels=c(TRUE, FALSE),
                                labels=c("Good species", "Bad species")),
         condition=factor(condition,
                          levels=0:2,
                          labels=c("Exp2a",
                                   "Exp2b contingent",
                                   "Exp2b full"))) %>%
  group_by(condition, healthyMushroom) %>%
  filter(catTrial < 21 | catTrial == max(catTrial)) %>%
  mutate(catTrial=ifelse(catTrial > 20, (catTrial-20)/5 + 20, catTrial)) %>%
  mutate(catTrialDodged=catTrial-(as.numeric(freqlabel)-.5)/10)

ggplot(data4mGraph2, aes(x=catTrialDodged, y=prop, color=freqlabel, shape=freqlabel)) +
  geom_hline(yintercept=.5, color="darkgray") +
  geom_line(data=subset(data4mGraph, catTrial<21), position=position_dodge(width=-0.5)) +
  geom_linerange(data=subset(data4mGraph, catTrial<21), aes(ymin=lower, ymax=upper), position=position_dodge(width=-0.5)) +
  geom_line(position=position_dodge(width=-0.5), linetype="21") +
  ylab("p(approach)") +
  xlab("Trial in species") +
  ylim(0,1) +
  labs(color = "Species\nfrequency", shape="Species\nfrequency") +
  facet_grid(healthyMushroom~condition) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        strip.text = element_text(size=10),
        plot.margin=unit(c(0,0,0,0), "cm"),
        legend.margin=unit(0, "cm"),
        legend.key.height=unit(2, "line"),
        panel.border = element_rect(colour = "black", fill=NA)) +
  scale_y_continuous(breaks=c(1, .5, 0)) +
  scale_x_continuous(breaks=c(0, 10, 20, 22, 24, 26), labels=c("0", "10", "20", "", "40", "")) +
  geom_pointrange(data=data4mGraph%>%
               group_by(condition, healthyMushroom, freqlabel)%>%
               filter(catTrial == 1 | catTrial == max(catTrial)),
               aes(ymin=lower, ymax=upper),
               position=position_dodge(width=-.5))+
  guides(fill=FALSE) +
  theme(legend.position="right")
ggsave(filename="figures/FigS3.pdf", width=9, height=4, useDingbats=F)


# get posterior predictive interval for the mean approach proportion on the first trial
get_first_trial_post <- function(data, sim, c) {
  firstTrials <- filter(data, condition == c)$catTrial == 1
  firstsim <- rowMeans(sim[,firstTrials])
  c(lower=quantile(firstsim, .025),
    mean=mean(firstsim),
    upper=quantile(firstsim, .975))
}

get_first_trial_post(data1m, ypred0, 0)
get_first_trial_post(data1m, ypred1, 1)
get_first_trial_post(data1m, ypred2, 2)

get_first_trial_post(data4m, ypred3, 0)
get_first_trial_post(data4m, ypred4, 1)
get_first_trial_post(data4m, ypred5, 2)


## plot individual horizon coefficients from "free horizon" model

gethorizons <- function(extracted, condition) {
  horizondf <- data.frame(extracted$horizon) %>%
    gather(horizon, sample) %>%
    group_by(horizon) %>%
    summarize(low=quantile(sample, .025),
              mean=mean(sample),
              high=quantile(sample, .975)) %>%
    mutate(horizon=as.numeric(horizon))
  horizondf$condition <- condition
  horizondf
}


load("model_fits/regression_1m_freehorizon_cond0.RData")
extracted <- extract(fit, pars=c("horizon"))
horizons0 <- gethorizons(extracted, 0)
load("model_fits/regression_1m_freehorizon_cond1.RData")
extracted <- extract(fit, pars=c("horizon"))
horizons1 <- gethorizons(extracted, 1)

horizons <- rbind(horizons0, horizons1) %>%
  mutate(condition=factor(condition, labels=c("Exp 1a ", "Exp 1b contingent")))

ggplot(horizons, aes(x=horizon, y=mean, ymin=low, ymax=high, group=condition)) +
  geom_line(aes(color=condition)) +
  geom_ribbon(aes(fill=condition), alpha=.2) +
  theme_minimal() +
  ylab("coefficient value") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        plot.margin=unit(c(0.2,0.2,0,0), "cm"),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = "bottom")
ggsave(filename="figures/RichFigS2.pdf", width=6.5, height=5, useDingbats=F)


## print the proportion of individuals in each scenario showing significant horizon sensitivity
get_indiv_summary <- function(exp, condition, stat) {
  load(paste0("model_fits/regression_", exp, "_cond", condition, ".RData"))
  extracted <- extract(fit, pars=c(stat))
  indivdf <- data.frame(extracted[[stat]]) %>%
    gather(indiv, sample) %>%
    group_by(indiv) %>%
    summarize(low=quantile(sample, .025),
              mean=mean(sample),
              median=median(sample),
              high=quantile(sample, .975))
  c(meanover0=mean(indivdf$mean > 0),
    lowover0=mean(indivdf$low > 0))
}

get_indiv_summary("1m", 0, "subj_horizon")
get_indiv_summary("1m", 1, "subj_horizon")
get_indiv_summary("1m", 2, "subj_horizon")
get_indiv_summary("4m", 0, "subj_horizon")
get_indiv_summary("4m", 1, "subj_horizon")
get_indiv_summary("4m", 2, "subj_horizon")
