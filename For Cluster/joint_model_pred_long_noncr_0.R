gender <- "male"
start.age <- 70
end.age <- 90

setwd("/n/home00/jgiardina/stroke_risk_pred")

source("functions_noncr.R")

set.seed(3218)

load("data/02_biolincc_data_pce.RData")

jm_fit <- readRDS(paste("output/jm_fit_noncr_", gender, ".rds", sep = ""))
test_ids <- readRDS(paste("output/training_data_noncr_", gender, ".rds", sep = ""))$test.ids

test_data <- test.data_prep(test_ids, biolincc_data_02, start.age)
saveRDS(test_data, paste("output/test_data_noncr_", gender, "_", start.age, ".rds", sep = ""))

jm_pred <- jm_pred_long_loop(jm_fit, test_data, start.age, end.age)
saveRDS(jm_pred, paste("output/jm_pred_long_noncr_", gender, "_", start.age, "_", end.age, ".rds", sep = ""))
