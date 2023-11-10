library(tidyverse)
library(R.matlab)

mst_fmri_test <- read_csv("Data/mst_fmri_test.csv")
mst_mat <- readMat("Data/msttDataWithin.mat")

mst_mat <- as.data.frame(mst_mat)$X1.1
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


# Lee and Stark dataset:
# Thresholds:
mean(num_missing_mat$missing_count) + 2 * sd(num_missing_mat$missing_count)
mean(num_missing_mat$missing_count) + sd(num_missing_mat$missing_count)

# Case between 1 and 2 sd
(43 - mean(num_missing_mat$missing_count)) / sd(num_missing_mat$missing_count)


#Wahlheim et al dataset:
#Plot cases:
ggplot(num_missing[order(num_missing$missing_count), ], aes(x = reorder(Subject, missing_count), y = missing_count)) +
  geom_col()

# Thresholds:
mean(num_missing$missing_count) + 2 * sd(num_missing$missing_count)
mean(num_missing$missing_count) + sd(num_missing$missing_count)

# Case between 1 and 2 sd
(49 - mean(num_missing$missing_count)) / sd(num_missing$missing_count)
