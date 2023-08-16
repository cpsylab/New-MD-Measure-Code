library(tidyverse)
data <- read_csv("data_L5.csv")

data <- mutate(data, l5_REC = participant_l5_max - participant_l5_min)
data$participant_id <- 1:nrow(data)

excluded_ids <- c(1, 15, 17)
data <- data[!data$participant_id %in% excluded_ids, ]

correlation_result_1 <- cor.test(data$participant_auc_l5, data$LDI)
print(correlation_result_1)

correlation_result_2 <- cor.test(data$participant_auc_l5_scale, data$LDI)
print(correlation_result_2)

correlation_result_3 <- cor.test(data$l5_REC, data$REC)
print(correlation_result_3)

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

ggplot(data, aes(x = l5_REC, y = REC)) +
  geom_point() +
  geom_smooth(method = "lm") +   # Adds a linear regression line
  theme_minimal() +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Logistic5 Max - Min",
       y = "REC")
