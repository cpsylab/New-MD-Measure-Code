library(tidyverse)
library(ggstatsplot)
library(ggplot2)
library(ggsci)
library(lmerTest)
library(sjPlot)
library(ggpubr)
library(performance)

data1 <- read_csv("scripts/data_L5.csv")
data2 <- read_csv("scripts/data_L5_fMRI.csv")


data1 <- mutate(data1, l5_REC = participant_l5_max - participant_l5_min)
data1$participant_id <- 1:nrow(data1)

excluded_ids <- c(1, 15, 17) # participant 16 also excluded in supplemental analysis
data1 <- data1[!data1$participant_id %in% excluded_ids, ]


data2 <- mutate(data2, l5_REC = participant_l5_max - participant_l5_min)
data2$participant_id <- 1:nrow(data2)

excluded_ids <- c(41,51,67,68,69) # participant 72 also excluded in supplemental analysis
data2 <- data2[!data2$participant_id %in% excluded_ids, ]


data1 <- data1 %>% mutate(study = 1)
data2 <- data2 %>% mutate(study = 2)
data <- bind_rows(data1, data2)
data <- mutate(data, lmdi = 1 - participant_auc_l5_scale)



# MIXED EFFECTS MODELS TO TEST LDI AND REC ASSOCIATIONS
# lmdi = lambda measure ; l5_REC = delta measure
ldi_m <- lmer(data = data, scale(LDI) ~ scale(lmdi) + scale(l5_REC) + (1|study))
rec_m <- lmer(data = data, scale(REC) ~ scale(lmdi) + scale(l5_REC) + (1|study))
lmdi_lrec_m <- lmer(data = data, scale(lmdi) ~ scale(l5_REC) + (1|study))


tab_model(ldi_m)
tab_model(rec_m)
tab_model(lmdi_lrec_m)


# LDI Figure
#pdf("ldi-fig.pdf", width=5, height=5)
fig_ldi <- ggplot(data, aes(x=lmdi, y=LDI, group=study)) + 
  geom_point(aes(colour=factor(study)),show.legend = FALSE) + 
  geom_smooth(aes(colour=factor(study)), method="lm",show.legend = FALSE) + 
  xlab(expression(lambda)) + 
  ylab("LDI") +
  scale_color_aaas(labels = c("Study 1", "Study 2")) + 
  theme_bw() + 
  labs(colour = "Study") +
  theme(
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.position = c(0.8, 0.15), 
    legend.title = element_text(size=16),
    legend.text = element_text(size=14),
    legend.background = element_rect(
      colour = "black", 
      size = 0.8
    )
  ) # REC Figure
fig_ldi
#dev.off()
#pdf("rec-fig.pdf", width=5, height=5)
fig_rec <- ggplot(data, aes(x=l5_REC, y=REC, group=study)) + 
  geom_point(aes(colour=factor(study)),show.legend = FALSE) + 
  geom_smooth(aes(colour=factor(study)), method="lm",show.legend = FALSE) + 
  xlab(expression(Delta)) + 
  ylab("REC") +
  scale_color_aaas(labels = c("Study 1", "Study 2")) + 
  theme_bw() + 
  labs(colour = "Study") + 
  theme(
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.position = c(0.8, 0.15), 
    legend.title = element_text(size=16),
    legend.text = element_text(size=14),
    legend.background = element_rect(
      colour = "black", 
      size = 1
    )
  )
#dev.off()
#pdf("ldi-rec-fig.pdf", width=5, height=5)
# Figure to show new values are decorrelated
fig_ldi_rec <- ggplot(data, aes(x=l5_REC, y=lmdi, group=study)) + 
  geom_point(aes(colour=factor(study)),show.legend = FALSE) + 
  geom_smooth(aes(colour=factor(study)), method="lm",show.legend = FALSE) + 
  xlab(expression(Delta)) + 
  ylab(expression(lambda)) +
  scale_color_aaas(labels = c("Study 1", "Study 2")) + 
  theme_bw() + 
  labs(colour = "Study") + 
  theme(
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.position = c(0.8, 0.15), 
    legend.title = element_text(size=16),
    legend.text = element_text(size=14),
    legend.background = element_rect(
      colour = "black", 
      size = 1
    )
  )
#dev.off()
fig_ldi_rec

fig_ldi_oldrec <- ggplot(data, aes(x=lmdi, y=REC, group=study)) + 
  geom_point(aes(colour=factor(study))) + 
  geom_smooth(aes(colour=factor(study)), method="lm") + 
  xlab(expression(lambda)) + 
  ylab("REC") +
  scale_color_aaas(labels = c("Study 1", "Study 2")) + 
  theme_bw() + 
  labs(colour = "Study") + 
  theme(
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.position = c(0.8, 0.15), 
    legend.title = element_text(size=16),
    legend.text = element_text(size=14),
    legend.background = element_rect(
      colour = "black", 
      size = 1
    )
  )

fig_ldi_oldrec

fig_oldldi_rec <- ggplot(data, aes(x=l5_REC, y=LDI, group=study)) + 
  geom_point(aes(colour=factor(study)),show.legend = FALSE) + 
  geom_smooth(aes(colour=factor(study)), method="lm",show.legend = FALSE) + 
  xlab(expression(Delta)) + 
  ylab("LDI") +
  scale_color_aaas(labels = c("Study 1", "Study 2")) + 
  theme_bw() + 
  labs(colour = "Study") + 
  theme(
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.position = c(0.8, 0.15), 
    legend.title = element_text(size=16),
    legend.text = element_text(size=14),
    legend.background = element_rect(
      colour = "black", 
      size = 1
    )
  )

fig_oldldi_rec

#pdf("fig.pdf", width=21, height=14)
combined_fig <- ggarrange(fig_ldi, fig_rec, fig_ldi_rec, fig_oldldi_rec,
                          fig_ldi_oldrec, labels = c("A", "B", "C", "D", "E"),
                          font.label = list(size = 20), ncol = 3, nrow = 2)
#dev.off()

combined_fig

