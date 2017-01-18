library("ggplot2")
library("dplyr")
library("grid")
source("optimal_models.R")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

## number of times to simulate each scenario
nsim=1000

## simulate data for the finite horizon
simulated_data_finite <- rbind.data.frame(
  full_info_simulator(nsims=nsim, ntrials=32),
  full_info_simulator(nsims=nsim, ntrials=32, p_healthy=.67),
  finite_contingent_simulator(nsims=nsim, horizon=1),
  finite_contingent_simulator(nsims=nsim, horizon=2),
  finite_contingent_simulator(nsims=nsim, horizon=4),
  finite_contingent_simulator(nsims=nsim, horizon=8),
  finite_contingent_simulator(nsims=nsim, horizon=16),
  finite_contingent_simulator(nsims=nsim, horizon=32),
  finite_contingent_simulator(nsims=nsim, horizon=1, p_healthy=.67),
  finite_contingent_simulator(nsims=nsim, horizon=2, p_healthy=.67),
  finite_contingent_simulator(nsims=nsim, horizon=4, p_healthy=.67),
  finite_contingent_simulator(nsims=nsim, horizon=8, p_healthy=.67),
  finite_contingent_simulator(nsims=nsim, horizon=16, p_healthy=.67),
  finite_contingent_simulator(nsims=nsim, horizon=32, p_healthy=.67)
)

## simulate data for the infinite horizon
simulated_data_infinite <- rbind.data.frame(
  full_info_simulator(nsims=nsim, ntrials=32),
  full_info_simulator(nsims=nsim, ntrials=32, p_healthy=.67),
  infinite_contingent_simulator(nsims=nsim, freq=.01),
  infinite_contingent_simulator(nsims=nsim, freq=.1),
  infinite_contingent_simulator(nsims=nsim, freq=.4),
  infinite_contingent_simulator(nsims=nsim, freq=.9),
  infinite_contingent_simulator(nsims=nsim, freq=.01, p_healthy=.67),
  infinite_contingent_simulator(nsims=nsim, freq=.1, p_healthy=.67),
  infinite_contingent_simulator(nsims=nsim, freq=.4, p_healthy=.67),
  infinite_contingent_simulator(nsims=nsim, freq=.9, p_healthy=.67))

## graph finite horizon data
simulated_data_finite <- simulated_data_finite %>% mutate(value=factor(ifelse(phealthy>.5, "2/3 positive", "2/3 negative"))) %>%
  mutate(value=factor(value, levels=levels(value)[c(2,1)]),
         horizon=factor(ifelse(is.na(horizon), "full\ninfo", horizon), levels=c(32, 16, 8, 4, 2, 1, "full\ninfo")))

ggplot(simulated_data_finite, aes(x=trial, y=response, linetype=horizon, shape=horizon, color=horizon)) +
  geom_hline(yintercept=.5, color="darkgray") +
  geom_point(data=subset(simulated_data_finite, firstlasttrials)) +
  geom_line() +
  theme_minimal() +
  facet_grid(value~.) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        strip.text = element_text(size=10),
        plot.margin=unit(c(0,0,0,0), "cm"),
        legend.margin=unit(0, "cm"),
        panel.border = element_rect(colour = "black", fill=NA)) +
  labs(y="p(approach)", shape="Horizon", color="Horizon", linetype="Horizon") +
  scale_linetype_manual(values=c(rep(1, 6), 4)) +
  scale_color_manual(values=c(gg_color_hue(6), "#777777")) +
  scale_shape_manual(values=c(16,17,15,3,4,7,8)) +
  ylim(0,1)
ggsave(filename="figures/Fig1.pdf", width=5, height=4, useDingbats=F)

## graph infinite horizon data
simulated_data_infinite <- simulated_data_infinite %>% mutate(value=factor(ifelse(phealthy>.5, "2/3 positive", "2/3 negative"))) %>%
  mutate(value=factor(value, levels=levels(value)[c(2,1)]),
         frequency=factor(ifelse(is.na(freq), "full\ninfo", freq), levels=c(.9, .4, .1, .01, "full\ninfo")))

ggplot(simulated_data_infinite, aes(x=trial, y=response, linetype=frequency, shape=frequency, color=frequency)) +
  geom_hline(yintercept=.5, color="darkgray") +
  geom_point(data=subset(simulated_data_infinite, firstlasttrials)) +
  geom_line() +
  theme_minimal() +
  facet_grid(value~.) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        strip.text = element_text(size=10),
        plot.margin=unit(c(0,0,0,0), "cm"),
        legend.margin=unit(0, "cm"),
        panel.border = element_rect(colour = "black", fill=NA)) +
  labs(x="trial with prospect", y="p(approach)", linetype="Frequency", shape="Frequency", color="Frequency") +
  scale_linetype_manual(values=c(rep(1, 4), 4)) +
  scale_color_manual(values=c(gg_color_hue(4), "#777777")) +
  ylim(0,1)
ggsave(filename="figures/Fig5.pdf", width=5, height=4, useDingbats=F)
