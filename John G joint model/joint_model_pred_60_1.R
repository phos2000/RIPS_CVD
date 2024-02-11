gender <- "female"
start.age <- 60
end.age <- 90

setwd("/n/home00/jgiardina/stroke_risk_pred")

source("functions.R")

set.seed(3218)

load("data/02_biolincc_data_pce.RData")

jm_fit <- readRDS(paste("output/jm_fit_", gender, ".rds", sep = ""))
test_ids <- readRDS(paste("output/training_data_", gender, ".rds", sep = ""))$test.ids

test_data <- test.data_prep(test_ids, biolincc_data_02, start.age)
saveRDS(test_data, paste("output/test_data_", gender, "_", start.age, ".rds", sep = ""))

jm_pred <- jm_pred_loop(jm_fit, test_data$jm.pred.data, start.age, end.age)
saveRDS(jm_pred, paste("output/jm_pred_", gender, "_", start.age, "_", end.age, ".rds", sep = ""))
