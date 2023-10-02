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

excluded_ids <- c(1, 15, 17)
data1 <- data1[!data1$participant_id %in% excluded_ids, ]

data2 <- mutate(data2, l5_REC = participant_l5_max - participant_l5_min)
data2$participant_id <- 1:nrow(data2)

excluded_ids <- c(41,51,67,68,69)
data2 <- data2[!data2$participant_id %in% excluded_ids, ]

data1 <- data1 %>% mutate(study = 1)
data2 <- data2 %>% mutate(study = 2)

data <- bind_rows(data1, data2)
data <- mutate(data, lmdi = 1 - participant_auc_l5_scale)

# MIXED EFFECTS MODELS TO TEST LDI AND REC ASSOCIATIONS
ldi_m <- lmer(data = data, scale(LDI) ~ scale(lmdi) + scale(l5_REC) + (1|study))
ldi_rs_m <- lmer(data = data, LDI ~ lmdi + l5_REC + (1 + lmdi + l5_REC|study))
rec_m <- lmer(data = data, scale(REC) ~ scale(lmdi) + scale(l5_REC) + (1|study))
rec_lm <- lm(data=data, scale(REC) ~ scale(lmdi) + scale(l5_REC))
lmdi_lrec_m <- lmer(data = data, scale(lmdi) ~ scale(l5_REC) + (1|study))

tab_model(rec_lm)
tab_model(ldi_m)
tab_model(rec_m)
tab_model(ldi_rs_m)
tab_model(lmdi_lrec_m)

performance::icc(rec_m)

performance::r2_(rec_m)

vars <- insight::get_variance(
  rec_m,
  tolerance=1e-5,
  name_fun = "r2()",
  name_full = "r-squared",
  verbose=T
)

vars_ldi <- insight::get_variance(
  ldi_m,
  tolerance=1e-5,
  name_fun = "r2()",
  name_full = "r-squared",
  verbose=T
)


performance::r2_nakagawa(rec_m)


tab_model(ldi_m)
tab_model(rec_m)
tab_model(lmdi_lrec_m)
summary(ldi_m)
summary(rec_m)
summary(lmdi_lrec_m)

# LDI Figure
#pdf("ldi-fig.pdf", width=5, height=5)
fig_ldi <- ggplot(data, aes(x=lmdi, y=LDI, group=study)) + 
  geom_point(aes(colour=factor(study)),show.legend = FALSE) + 
  geom_smooth(aes(colour=factor(study)), method="lm",show.legend = FALSE) + 
  xlab("λ") + 
  ylab("LDI") +
  scale_color_aaas(labels = c("Study 1", "Study 2")) + 
  theme_bw() + 
  labs(colour = "Study") +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
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
  xlab("Δ") + 
  ylab("REC") +
  scale_color_aaas(labels = c("Study 1", "Study 2")) + 
  theme_bw() + 
  labs(colour = "Study") + 
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
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
fig_ldi_rec <- ggplot(data, aes(x=lmdi, y=l5_REC, group=study)) + 
  geom_point(aes(colour=factor(study)),show.legend = FALSE) + 
  geom_smooth(aes(colour=factor(study)), method="lm",show.legend = FALSE) + 
  xlab("λ") + 
  ylab("Δ") +
  scale_color_aaas(labels = c("Study 1", "Study 2")) + 
  theme_bw() + 
  labs(colour = "Study") + 
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
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
  xlab("λ") + 
  ylab("REC") +
  scale_color_aaas(labels = c("Study 1", "Study 2")) + 
  theme_bw() + 
  labs(colour = "Study") + 
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
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

fig_oldldi_rec <- ggplot(data, aes(x=LDI, y=l5_REC, group=study)) + 
  geom_point(aes(colour=factor(study)),show.legend = FALSE) + 
  geom_smooth(aes(colour=factor(study)), method="lm",show.legend = FALSE) + 
  xlab("LDI") + 
  ylab("Δ") +
  scale_color_aaas(labels = c("Study 1", "Study 2")) + 
  theme_bw() + 
  labs(colour = "Study") + 
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
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
                          font.label = list(size = 24), ncol = 3, nrow = 2)
#dev.off()

combined_fig

#cor.test(data$participant_auc_l5, data$LDI)

cor.test(data$lmdi, data$LDI)

cor.test(data$lmdi, data$REC)

cor.test(data$l5_REC, data$REC)

cor.test(data$l5_REC, data$LDI)

ggscatterstats(
  data = data,
  x = lmdi,
  y = LDI,
  xlab = "Area under curve scaled",
  ylab = "LDI",
  xlim = c(0, 1),
  ylim = c(0, 1),
  marginal = FALSE, # remove marginal histograms
  ggtheme = ggplot2::theme_minimal()
)

ggplot(data, aes(x = participant_auc_l5, y = LDI)) +
  geom_point() +
  geom_smooth(method = "lm") +   # Adds a linear regression line
  theme_minimal() +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Area under curve",
       y = "LDI")

ggplot(data, aes(x = lmdi, y = LDI)) +
  geom_point() +
  geom_smooth(method = "lm") +   # Adds a linear regression line
  theme_minimal() +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Area under curve scaled",
       y = "LDI")

ggplot(data, aes(x = lmdi, y = REC)) +
  geom_point() +
  geom_smooth(method = "lm") +   # Adds a linear regression line
  theme_minimal() +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Area under curve scaled",
       y = "REC")

ggplot(data, aes(x = l5_REC, y = REC)) +
  geom_point() +
  geom_smooth(method = "lm") +   # Adds a linear regression line
  theme_minimal() +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Logistic5 Max - Min",
       y = "REC")

ggplot(data, aes(x = l5_REC, y = LDI)) +
  geom_point() +
  geom_smooth(method = "lm") +   # Adds a linear regression line
  theme_minimal() +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Logistic5 Max - Min",
       y = "LDI")

ggplot(data, aes(x = REC, y = LDI)) +
  geom_point() +
  geom_smooth(method = "lm") +   # Adds a linear regression line
  theme_minimal() +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "REC",
       y = "LDI")

