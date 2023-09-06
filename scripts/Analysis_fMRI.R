library(tidyverse)
library(ggstatsplot)
library(lme3)
data1 <- read_csv("data_L5.csv")
data2 <- read_csv("data_L5_fMRI.csv")

data1 <- mutate(data, l5_REC = participant_l5_max - participant_l5_min)
data1$participant_id <- 1:nrow(data1)

excluded_ids <- c(1, 15, 17)
data1 <- data1[!data1$participant_id %in% excluded_ids, ]

data2 <- mutate(data2, l5_REC = participant_l5_max - participant_l5_min)
data2$participant_id <- 1:nrow(data2)

excluded_ids <- c(41,51,67,68,69,72)
data2 <- data2[!data$participant_id %in% excluded_ids, ]

data1 <- data1 %>% mutate(study = 1)
data2 <- data2 %>% mutate(study = 2)

data <- bind_rows(data1, data2)






#cor.test(data$participant_auc_l5, data$LDI)

cor.test(data$participant_auc_l5_scale, data$LDI)

cor.test(data$participant_auc_l5_scale, data$REC)

cor.test(data$l5_REC, data$REC)

cor.test(data$l5_REC, data$LDI)

ggscatterstats(
  data = data,
  x = participant_auc_l5_scale,
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

ggplot(data, aes(x = participant_auc_l5_scale, y = LDI)) +
  geom_point() +
  geom_smooth(method = "lm") +   # Adds a linear regression line
  theme_minimal() +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Area under curve scaled",
       y = "LDI")

ggplot(data, aes(x = participant_auc_l5_scale, y = REC)) +
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
