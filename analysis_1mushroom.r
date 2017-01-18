library("dplyr")
library("ggplot2")
library("tidyr")
library("grid")

## import data
data1m<- read.csv('data_exp1_forsupplement.csv') %>% filter(!exclude)

## create a data frame to feed to ggplot for graphing data
data1mGraph <- data1m %>%
  group_by(condition, healthySpecies, patchLength, patchTrial, uniqueid) %>%
  summarise(indivprop=mean(response)) %>%
  summarise(prop=mean(indivprop),
            propsd=sd(indivprop),
            nrow = n()) %>%
  ungroup() %>%
  mutate(lower = prop - propsd/sqrt(nrow),
         upper = prop + propsd/sqrt(nrow)) %>%
  mutate(patchLength=factor(patchLength, c(32, 16, 8, 4, 2, 1))) %>%
  mutate(healthySpecies=factor(healthySpecies, levels=c(TRUE, FALSE),
                               labels=c("Mostly-healthy", "Mostly-poisonous")),
         condition=factor(condition, levels=c(0,1,2),
                          labels=c("Exp1a", "Exp1b contingent", "Exp1b full"))) %>%
  ## remove trials between 17 and 31 for break in X axis
  filter(patchTrial < 17 | patchTrial == 32) %>%
  ## move trial 32 responses to appear where 20 should be on the X axis
  mutate(patchTrial=ifelse(patchTrial == 32, 20, patchTrial)) %>%
  ## manually dodge the x axis a bit
  mutate(patchTrialDodged=patchTrial-(as.numeric(patchLength)-2.5)/10)

## graph data
ggplot(data1mGraph, aes(x=patchTrialDodged, y=prop, color=patchLength, shape=patchLength)) +
  geom_hline(yintercept=.5, color="darkgray") +
  geom_line(data=subset(data1mGraph, patchTrial<20), position=position_dodge(width=-0.5)) +
  geom_linerange(data=subset(data1mGraph, patchTrial<20),
                 aes(ymin=lower, ymax=upper), position=position_dodge(width=-0.5)) +
  ## use dashed line between trial 16 and trial 32
  geom_line(position=position_dodge(width=-0.5), linetype="21") +
  ylab("p(approach)") +
  xlab("Trial in species") +
  facet_grid(healthySpecies~condition) +
  theme_minimal() +
  scale_y_continuous(breaks=c(1, .5, 0)) +
  ## set up axis labels/ticks for break in axis between 16 and 32
  scale_x_continuous(breaks=c(0, 8, 16, 18, 20),
                     labels=c("0", "8", "16", "", "32")) +
  ## set text sizes, style, minimize blank space on edge
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        strip.text = element_text(size=10),
        plot.margin=unit(c(0,0,0,0), "cm"),
        legend.margin=unit(0, "cm"),
        panel.border = element_rect(colour = "black", fill=NA)) +
  labs(color="Patch\nlength", shape="Patch\nlength") +
  ## put points at the beginning and end of each line
  geom_pointrange(data=subset(data1mGraph, (patchTrial == 1 |
                                            patchTrial == 20 |
                                            patchTrial==as.numeric(as.character(patchLength)))),
                  aes(ymin=lower, ymax=upper),
                  position=position_dodge(width=-0.5))+
  guides(fill=FALSE) +
  theme(legend.position="right")
ggsave(filename="figures/Fig2.pdf", width=9, height=4, useDingbats=F)

## extra analysis: is there a bias towards information-seeking with 0 trials
## left? Compare last trials of each patch in the full-info and contingent-info
## conditions to find out.
lastTrials <- data1m %>%
  filter(patchTrialLeft == 1) %>%
  filter(condition > 0)

## Looks like there is a strong effect of condition, controlling for expected reward
summary(glm(response ~ expectation + condition,
    data=lastTrials%>%mutate(condition=1-condition), family="binomial"))
