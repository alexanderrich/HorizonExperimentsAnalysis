library("dplyr")
library("ggplot2")
library("tidyr")
library("grid")

# load data
data4m<- read.csv('data_exp2_forsupplement.csv') %>% filter(!exclude)

## frequency quiz
freqquizdata <- read.csv('data_freqquiz_forsupplement.csv') %>% filter(!exclude)

## mean responses in frequency quiz for high- and low-frequency species
freqquizdata %>%
  mutate(freq=round(freq, 1)) %>%
  group_by(exp, freq, uniqueid) %>%
  summarize(response=mean(response)) %>%
  summarize(meanresp=mean(response),
            sdresp=sd(response))

## people usually got the frequency exactly right
mean(round(freqquizdata$freq*10) == freqquizdata$response + 0)

## create a data frame to feed to ggplot for graphing data
data4mGraph <- data4m %>%
  group_by(condition, healthySpecies, highFreq, catTrial, uniqueid) %>%
  summarise(indivprop=mean(response)) %>%
  summarise(prop=mean(indivprop),
            propsd=sd(indivprop),
            nrow = n()) %>%
  ungroup() %>%
  mutate(lower = prop - propsd/sqrt(nrow),
         upper=prop + propsd/sqrt(nrow)) %>%
  mutate(freqlabel = factor(highFreq, levels=c(TRUE, FALSE),
                            labels=c("high\n(f=0.4)", "low\n(f=0.1)"))) %>%
  mutate(healthySpecies=factor(healthySpecies, levels=c(TRUE, FALSE),
                                labels=c("Mostly-healthy", "Mostly-poisonous")),
         condition=factor(condition,
                          levels=c(0,1,2),
                          labels=c("Exp2a",
                                   "Exp2b contingent",
                                   "Exp2b full"))) %>%
  group_by(condition, healthySpecies) %>%
  ## remove trials between 20 and the final trial for break in X axis
  filter(catTrial < 21 | catTrial == max(catTrial)) %>%
  ## move final trial closer on the X axis
  mutate(catTrial=ifelse(catTrial > 20, (catTrial-20)/5 + 20, catTrial)) %>%
  ## manually dodge the x axis a bit
  mutate(catTrialDodged=catTrial-(as.numeric(freqlabel)-.5)/10) %>%
  ungroup()

## graph data
ggplot(data4mGraph, aes(x=catTrialDodged, y=prop, color=freqlabel, shape=freqlabel)) +
  geom_hline(yintercept=.5, color="darkgray") +
  geom_line(data=subset(data4mGraph, catTrial<21), position=position_dodge(width=-0.5)) +
  geom_linerange(data=subset(data4mGraph, catTrial<21),
                 aes(ymin=lower, ymax=upper), position=position_dodge(width=-0.5)) +
  ## use dashed line between trial 20 and final trail
  geom_line(position=position_dodge(width=-0.5), linetype="21") +
  ylab("p(approach)") +
  xlab("Trial in species") +
  ylim(0,1) +
  labs(color = "Species\nfrequency", shape="Species\nfrequency") +
  facet_grid(healthySpecies~condition) +
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
  ## set up axis labels/ticks for break in axis after trial 20
  scale_x_continuous(breaks=c(0, 10, 20, 22, 24, 26),
                     labels=c("0", "10", "20", "", "40", "")) +
  ## put points at the beginning and end of each line
  geom_pointrange(data=data4mGraph%>%
               group_by(condition, healthySpecies, freqlabel)%>%
               filter(catTrial == 1 | catTrial == max(catTrial)),
               aes(ymin=lower, ymax=upper),
               position=position_dodge(width=-.5))+
  guides(fill=FALSE) +
  theme(legend.position="right")
ggsave(filename="figures/Fig5.pdf", width=9, height=4, useDingbats=F)

## is the apparent effect of frequency driven by the fact that trials with
## high-frequency species tend to come earlier in the habitat? To test, look at
## how response to first encounter is affected by frequency and by trial in
## habitat
firstbyfreq <- data4m %>% filter(condition == 0) %>%
  group_by(uniqueid, species) %>%
  slice(1) %>%
  ungroup() %>%
  select(habitatTrial, response, highFreq)
## looking at first trial for each patch, there's a high effect of frequency, but
## no effect of how late in the habitat the trial occurs
summary(glm(response ~ highFreq + habitatTrial, data=firstbyfreq, family="binomial"))

## repeat for experiment 2b
firstbyfreq <- data4m %>% filter(condition == 1) %>%
  group_by(uniqueid, species) %>%
  slice(1) %>%
  ungroup() %>%
  select(habitatTrial, response, highFreq)
## again, effect of frequency, but not of trial in habitat
summary(glm(response ~ highFreq + habitatTrial, data=firstbyfreq, family="binomial"))

