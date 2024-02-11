gender <- "female"

setwd("/n/home00/jgiardina/stroke_risk_pred")

source("functions.R")

set.seed(3218)

load("data/02_biolincc_data_pce.RData")

#Create baseline age variable
biolincc_data_02 <- biolincc_data_02 %>%
  group_by(shield.id) %>%
  mutate(age_baseline = min(age_merged)) %>%
  filter(ttoevent_ischstroke_merged != 0)

biolincc_data_02 <- biolincc_data_02[complete.cases(biolincc_data_02[,c("ttoevent_ischstroke_merged", "event_ischstroke_merged", "ttoevent.death", "event.death", "cursmk", "prior_cvd_merged", "afib_previous_merged", "dm_previous_merged", "htnmed_merged", "sbp")]),] %>%
  mutate(time_elapsed = time_elapsed/365.25,
         ttoevent_ischstroke_merged = ttoevent_ischstroke_merged/365.25,
         ttoevent.death = ttoevent.death/365.25) %>%
  mutate(age_exact = age_baseline + time_elapsed,
         age_event = age_baseline + ttoevent_ischstroke_merged) %>%
  mutate(age_death = age_baseline + ttoevent.death,
         event_cause = case_when(
           event_ischstroke_merged == 1 ~ "isch.stroke",
           event_ischstroke_merged == 0 & event.death == 0 ~ "alive",
           age_death == age_event & event.death == 1 ~ "death",
           age_death >= age_event ~ "alive"
         ))

training.data <- data_prep(biolincc_data_02, gender)
saveRDS(training.data, paste("output/training_data_", gender, ".rds", sep = ""))
print("Finished Preparing Training Data")

models <- model_prep(training.data)
saveRDS(models, paste("output/models_", gender, ".rds", sep = ""))
print("Finished Preparing Models")

jm_fitted <- jm_run(models)
saveRDS(jm_fitted, paste("output/jm_fit_", gender, ".rds", sep = ""))
print("Finished Fitting Joint Model")
