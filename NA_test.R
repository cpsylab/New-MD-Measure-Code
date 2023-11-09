library(tidyverse)
library(R.matlab)

mst_fmri_test <- read_csv("Data/mst_fmri_test.csv")
mst_fmri_age <- read_csv("Data/mst_icv.csv")
mst_mat <- readMat("Data/msttDataWithin.mat")

mst_mat <- as.data.frame(mst_mat)

mst_mat <- mst_mat$X1.1

mst_mat_res <- mst_mat$decisionOSN[1,1:4032]
subject_column <- rep(1:21, each = 192)


mst_mat_dataframe <- data.frame(Subject = subject_column, decision = mst_mat_res)


mst_fmri_test <- mst_fmri_test %>%
  mutate(Subject = factor(Subject)) %>%
  mutate(Subject_ID = as.integer(Subject))


num_missing <- mst_fmri_test %>%
  group_by(Subject) %>%
  summarise(missing_count = sum(is.na(TestObj.RESP)))

num_missing_mat <- mst_mat_dataframe %>%
  group_by(Subject) %>%
  summarise(missing_count = sum(is.na(decision)))

mean(num_missing_mat$missing_count) + 1 * sd(num_missing_mat$missing_count)

(43 - mean(num_missing_mat$missing_count)) / sd(num_missing_mat$missing_count)

ggplot(num_missing[order(num_missing$missing_count), ], aes(x = reorder(Subject, missing_count), y = missing_count)) +
  geom_col()

mean(num_missing$missing_count) + 2 * sd(num_missing$missing_count)

(49 - mean(num_missing$missing_count)) / sd(num_missing$missing_count)

mean(num_missing$missing_count) + sd(num_missing$missing_count)

x <- mean(num_missing$missing_count) + 2 * sd(num_missing$missing_count)
# Threshold for number of missing values

# Get a vector of subject IDs with more than x missing responses
subjects_with_missing <- mst_fmri_test %>%
  group_by(Subject) %>%
  summarise(missing_count = sum(is.na(TestObj.RESP))) %>%
  filter(missing_count > x) %>%
  pull(Subject)

subjects_with_missing


compute_p_value <- function(data) {
  contingency_table <- table(data$TestObj.RESP, data$ItemType)
  
  # Check for positive entries
  if(all(contingency_table == 0)) {
    return(NA_real_)
  }
  
  # Perform chi-square test if possible, otherwise return NA
  tryCatch({
    chi_test <- chisq.test(contingency_table)
    return(chi_test$p.value)
  }, error = function(e) {
    return(NA_real_)
  })
}

# Compute the p-values for each subject
p_values <- mst_fmri_test %>%
  group_by(Subject) %>%
  nest() %>%
  mutate(p_value = map_dbl(data, compute_p_value)) %>%
  ungroup()

# Filter subjects with p-value > 0.05
random_responders <- p_values %>%
  filter(is.na(p_value) | p_value > 0.05) %>%
  pull(Subject)

random_responders

setdiff(mst_fmri_test$Subject, mst_fmri_age$Subject) # 27 and 50 lack fmri data

mst_fmri_age

mst_fmri_age <- mst_fmri_age[!mst_fmri_age$Subject %in% random_responders, ]

mean(mst_fmri_age$Age_Years)
sd(mst_fmri_age$Age_Years)



