library(dplyr)
library(tidyr)
library(JMbayes2)
library(foreach)
library(survival)
library(ggplot2)
library(ggsci)
library(extrafont)
loadfonts(device = "win")
loadfonts(device = "postscript")
loadfonts(device = "pdf")
library(ggpubr)
library(stringr)
library(bshazard)
library(muhaz)
library(survPen)

#Dufouil C, Beiser A, McLure LA, et al. Revised Framingham Stroke Risk Profile to Reflect Temporal Trends. Circulation. 2017;135(12):1145-1159. doi:10.1161/CIRCULATIONAHA.115.021275
rfsrs <- Vectorize(function(female, age, current.smk, prevalent.cvd, prevalent.af, dm, ahtn, sbp, years = 10){
  age65 <- as.numeric(age >= 65)
  dm.age65 <- dm*age65
  dm.not.age65 <- dm*(1-age65)
  sbp.no.ahtn <- (sbp-120)*(1-ahtn)/10
  sbp.ahtn <- (sbp-120)*(ahtn)/10
  rf.vector <- c(age, current.smk, prevalent.cvd, prevalent.af, age65, dm.age65, dm.not.age65, ahtn, sbp.no.ahtn, sbp.ahtn)
  rf.vector[1] <- rf.vector[1]/10
  if(as.numeric(female) == 1){
    L <- as.numeric(rf.vector %*% c(0.87938, 0.51127, -0.03035, 1.20720, 0.39796, 1.07111, 0.06565, 0.13085, 0.11303, 0.17234))
    M <- 6.6170719
    sb <- c("1" = 0.99901, "2" = 0.99697, "3" = 0.99424, "4" = 0.99074, "5" = 0.98710, "6" = 0.98327, "7" = 0.97890, "8" = 0.97326, "9" = 0.96581, "10" = 0.95911)[as.character(years)]
    return <- 1-sb^exp(L-M)
  }
  if(as.numeric(female) == 0){
    L <- as.numeric(rf.vector %*% c(0.49716, 0.47254, 0.45341, 0.08064, 0.45426, 1.35304, 0.34385, 0.82598, 0.27323, 0.09793))
    M <- 4.4227101
    sb <- c("1" = 0.99598, "2" = 0.99478, "3" = 0.98919, "4" = 0.98531, "5" = 0.98059, "6" = 0.97661, "7" = 0.96983, "8" = 0.96179, "9" = 0.95509, "10" = 0.94451)[as.character(years)]
    return <- 1-sb^exp(L-M)
  }
  return(return)
}, vectorize.args = c("female", "age", "current.smk", "prevalent.cvd", "prevalent.af", "dm", "ahtn", "sbp"))

cloglog <- Vectorize(function(x){
  return(log(-log(1-x)))
})

data_prep <- function(data, gender){
  gender_num <- case_when(gender == "female" ~ 1,
                          gender == "male" ~ 0)
  ids <- unique(data$shield.id[data$female_merged == gender_num])
  training.set <- sample(ids, length(ids)*0.5)
  test.set <- setdiff(ids, training.set)
  training.data_long <- data[data$shield.id %in% training.set,]
  training.data_event <- unique(training.data_long[,c("shield.id", "age_baseline", "age_event", "event_ischstroke_merged")])
  
  return(list(
    gender = gender,
    training.ids = training.set,
    test.ids = test.set,
    training.data_long = training.data_long,
    training.data_surv = training.data_event
  ))
}

model_prep <- function(data.object){
  long_data <- data.object[["training.data_long"]]
  surv_data <- data.object[["training.data_surv"]]
  
  cursmk_model <- lme(fixed = cursmk ~ age_exact, random = ~ age_exact | shield.id, data = long_data)
  
  prior.cvd_model <- lme(fixed = prior_cvd_merged ~ age_exact, random = ~ age_exact | shield.id, data = long_data)
  
  afib_model <- lme(fixed = afib_previous_merged ~ age_exact, random = ~ age_exact | shield.id, data = long_data)
  
  dm_model <- lme(fixed = dm_previous_merged ~ age_exact, random = ~ age_exact | shield.id, data = long_data)
  
  htnmed_model <- lme(fixed = htnmed_merged ~ age_exact, random = ~ age_exact | shield.id, data = long_data)
  
  sbp_model <- lme(fixed = sbp ~ age_exact, random = ~ age_exact | shield.id, data = long_data, na.action = na.exclude)
  
  long_models <- list(cursmk_model, prior.cvd_model, afib_model, dm_model, htnmed_model, sbp_model)
  
  cox_model <- coxph(Surv(time = age_baseline, time2 = age_event, event = event_ischstroke_merged) ~ 1, data = surv_data)
  
  return(list(
    long_models = long_models,
    surv_model = cox_model,
    surv_data = surv_data
  ))
}

jm_run <- function(model.list){
  
  fx.forms <- list(
    "cursmk" = ~ value(cursmk),
    "prior_cvd_merged" = ~ value(prior_cvd_merged),
    "afib_previous_merged" = ~ value(afib_previous_merged),
    "dm_previous_merged" = ~ value(dm_previous_merged),
    "htnmed_merged" = ~ value(htnmed_merged),
    "sbp" = ~ value(sbp)
  )
  
  jm_fit <- jm(model.list$surv_model, model.list$long_models, time_var = "age_exact", functional_forms = fx.forms, data_Surv = model.list$surv_data, control = list(n_burnin = 5000L, n_iter = 25000L, n_chains = 5, Bsplines_degree = 3, cores = 5))
  
  return(jm_fit)
  
}

baseline_run <- function(data.object){
  long_data <- data.object[["training.data_long"]]
  
  #Note that baseline observation is assumed to be first observation with complete data
  baseline_data <- long_data %>%
    group_by(shield.id) %>%
    slice(which.min(age_exact)) %>%
    mutate(age_baseline = min(age_exact))
  
  baseline_fit <- survPen(formula = ~smf(age_event, df = 10) + cursmk + prior_cvd_merged + afib_previous_merged + dm_previous_merged + htnmed_merged + sbp, t1 = age_event, t0 = age_baseline, event = event_ischstroke_merged, data = baseline_data)
  
  return(baseline_fit)
}

time.dep_run <- function(data.object){
  long_data <- data.object[["training.data_long"]]
  
  #Note that baseline observation is assumed to be first observation with complete data
  time.dep_data <- long_data %>%
    group_by(shield.id) %>%
    arrange(visit, .by_group = TRUE) %>%
    mutate(age_baseline = age_exact) %>%
    mutate(age_event = lead(age_exact, default = age_event[which.max(visit)])) %>%
    mutate(event_ischstroke_time.dep = case_when(
      age_event == age_event[which.max(visit)] ~ event_ischstroke_merged,
      TRUE ~ 0
    ))
  
  time.dep_fit <- survPen(formula = ~smf(age_event, df = 10) + cursmk + prior_cvd_merged + afib_previous_merged + dm_previous_merged + htnmed_merged + sbp, t1 = age_event, t0 = age_baseline, event = event_ischstroke_time.dep, data = time.dep_data)
  
  return(time.dep_fit)
}

test.data_prep <- function(test.list, data, start.age){
  
  data <- data %>%
    group_by(shield.id) %>%
    mutate(age_baseline = min(age_merged)) %>%
    filter(ttoevent_ischstroke_merged != 0)
  
  data <- data[complete.cases(data[,c("ttoevent_ischstroke_merged", "event_ischstroke_merged", "ttoevent.death", "event.death", "cursmk", "prior_cvd_merged", "afib_previous_merged", "dm_previous_merged", "htnmed_merged", "sbp")]),] %>%
    mutate(time_elapsed = time_elapsed/365.25,
           ttoevent_ischstroke_merged = ttoevent_ischstroke_merged/365.25) %>%
    mutate(age_exact = age_baseline + time_elapsed,
           age_event = age_baseline + ttoevent_ischstroke_merged)
  
  test.data_pre.start <- data[data$shield.id %in% test.list & data$age_event > start.age & data$age_exact <= start.age,]
  test.data_post.start <- data[data$shield.id %in% test.list & data$age_event > start.age & data$age_exact > start.age,]
  
  test.data_surv <- unique(test.data_pre.start[,c("shield.id", "age_baseline", "age_event", "event_ischstroke_merged")])
  test.data_pre.start$event_ischstroke_merged <- 0
  test.data_pre.start$age_event <- start.age + 1e-7
  
  return(list(
    pre.data = test.data_pre.start,
    post.data = test.data_post.start,
    surv.data = test.data_surv
  ))
}

jm_pred_loop <- function(jm.fit.object, test.data, start.age, end.age){
  test.ids <- unique(test.data$pre.data$shield.id)
  partitions <- ceiling(seq_along(test.ids)/100)
  test.partitions <- split(test.ids, partitions)
  
  jm_pred <- foreach(group=names(test.partitions), .combine = rbind) %do% {
    
    data <- test.data$pre.data[test.data$pre.data$shield.id %in% test.partitions[[group]],]
    
    jm_pred_run <- predict(jm.fit.object, newdata = data, process = "event", times = seq(from = start.age, to = end.age, by = 1), return_newdata = TRUE, return_mcmc = FALSE, cores = 8)
    
    return(jm_pred_run)
  }
  return(jm_pred)
}

jm_pred_long_loop <- function(jm.fit.object, test.data, start.age, end.age){
  test.ids <- unique(test.data$pre.data$shield.id)
  partitions <- ceiling(seq_along(test.ids)/100)
  test.partitions <- split(test.ids, partitions)
  
  jm_pred <- foreach(group=names(test.partitions), .combine = rbind) %do% {
    
    data <- test.data$pre.data[test.data$pre.data$shield.id %in% test.partitions[[group]],]
    
    jm_pred_run <- predict(jm.fit.object, newdata = data, process = "longitudinal", type_pred = "response", times = seq(from = start.age, to = end.age, by = 1), return_newdata = TRUE, return_mcmc = FALSE, cores = 8)$newdata2
    
    return(jm_pred_run)
  }
  return(jm_pred)
}

test_long_data_prep <- function(long.data_raw){
  long.data <- long.data_raw[,c("shield.id",
                                "female_merged",
                                "age_exact",
                                "age_event",
                                "pred_cursmk",
                                "pred_prior_cvd_merged",
                                "pred_afib_previous_merged",
                                "pred_dm_previous_merged",
                                "pred_htnmed_merged",
                                "pred_sbp")] %>%
    dplyr::rename(cursmk = pred_cursmk,
           prior_cvd_merged = pred_prior_cvd_merged,
           afib_previous_merged = pred_afib_previous_merged,
           dm_previous_merged = pred_dm_previous_merged,
           htnmed_merged = pred_htnmed_merged,
           sbp = pred_sbp) %>%
    mutate(age_event = age_exact)
  
  return(long.data)
}

baseline_pred <- function(baseline.fit.object, test.long.data, gender){
  prediction_indv <- data.frame(surv = predict(baseline.fit.object, newdata = as.data.frame(test.long.data))$surv,
                                shield.id = test.long.data$shield.id,
                                age_exact = test.long.data$age_event,
                                gender = gender) %>%
    group_by(shield.id) %>%
    mutate(surv = surv/max(surv))
  
  prediction_avg <- prediction_indv %>%
    ungroup() %>%
    group_by(age_exact) %>%
    summarise(surv = mean(surv)) %>%
    mutate(method = "Baseline Model") %>%
    mutate(gender = gender)
  
  return(list(
    individual.predictions = prediction_indv,
    average.predictions = prediction_avg
  ))
}

time.dep_pred <- function(time.dep.fit.object, test.long.data, gender){
  prediction.object <- baseline_pred(time.dep.fit.object, test.long.data, gender)
  
  prediction.object[["average.predictions"]] <- prediction.object[["average.predictions"]] %>%
    mutate(method = "Time-Dependent Model")
  
  return(prediction.object)
}

rfsrs_pred <- function(test.long.data, gender){
  prediction_indv <- data.frame(shield.id = test.long.data$shield.id,
                                age_exact = test.long.data$age_exact,
                                surv.cond = 1 - rfsrs(female = test.long.data$female_merged,
                                                  age = test.long.data$age_exact,
                                                  current.smk = test.long.data$cursmk,
                                                  prevalent.cvd = test.long.data$prior_cvd_merged,
                                                  prevalent.af = test.long.data$afib_previous_merged,
                                                  dm = test.long.data$dm_previous_merged,
                                                  ahtn = test.long.data$htnmed_merged,
                                                  sbp = test.long.data$sbp,
                                                  years = 1)) %>%
    group_by(shield.id) %>%
    arrange(age_exact, .by_group = TRUE) %>%
    mutate(surv = cumprod(surv.cond)) %>%
    mutate(surv = surv/max(surv)) %>%
    mutate(gender = gender)
  
  prediction_avg <- prediction_indv %>%
    ungroup() %>%
    group_by(age_exact) %>%
    summarise(surv = mean(surv)) %>%
    mutate(method = "Revised Framingham Stroke Risk Score") %>%
    mutate(gender = gender)
  
  return(list(
    individual.predictions = prediction_indv,
    average.predictions = prediction_avg
  ))
}

emp_surv_curve <- function(test.data, end.age, gender){
  surv.data <- test.data$surv.data
  curve.object <- survfit(Surv(time = age_event, event = event_ischstroke_merged) ~ 1, data = surv.data, conf.type = "log-log")
  
  isch.stroke.curve <- data.frame(age_exact = curve.object$time, surv = curve.object$surv, surv.lower = curve.object$lower, surv.upper = curve.object$upper) %>%
    mutate(method = "Observed Stroke Survival") %>%
    mutate(gender = gender) %>%
    filter(age_exact <= end.age)
  
  return(list(
    isch.stroke_surv.curve = isch.stroke.curve[,c("age_exact", "surv", "method", "gender", "surv.lower", "surv.upper")]
  ))
}

emp_hazard_curve <- function(test.data, start.age, end.age, gender){
  surv.data <- test.data$surv.data
  surv.data_isch.stroke <- surv.data
  
  #curve.object_isch.stroke <- bshazard(Surv(time = age_event, event = event) ~ 1, data = surv.data_isch.stroke, nk = 10, degree = 3)
  
  curve.object_isch.stroke <- muhaz(times = surv.data_isch.stroke$age_event, delta = surv.data_isch.stroke$event_ischstroke_merged, min.time = start.age, max.time = end.age, n.est.grid = as.integer(end.age - start.age), bw.method = "local")
  
  # isch.stroke.curve <- data.frame(age_exact = curve.object_isch.stroke$time, hazard = curve.object_isch.stroke$hazard, hazard.lower = curve.object_isch.stroke$lower.ci, hazard.upper = curve.object_isch.stroke$upper.ci) %>%
  #   mutate(method = "Observed Stroke Hazard") %>%
  #   mutate(gender = gender) %>%
  #   filter(age_exact <= end.age & age_exact >= start.age + 1)
  
  isch.stroke.curve <- data.frame(age_exact = curve.object_isch.stroke$est.grid, hazard = curve.object_isch.stroke$haz.est) %>%
    mutate(method = "Observed Stroke Hazard") %>%
    mutate(gender = gender) %>%
    filter(age_exact <= end.age & age_exact >= start.age + 1)

  
  # return(list(
  #   isch.stroke_hazard.curve = isch.stroke.curve[,c("age_exact", "hazard", "method", "gender", "hazard.lower", "hazard.upper")]
  # ))
  
  return(list(
    isch.stroke_hazard.curve = isch.stroke.curve[,c("age_exact", "hazard", "method", "gender")]
  ))
}

prep_jm_surv_data <- function(jm_pred_data, gender){
  
  individual.predictions <- jm_pred_data %>%
    mutate(surv = 1- pred_CIF) %>%
    mutate(method = "Joint Model") %>%
    mutate(gender = gender)
  
  average.predictions <- jm_pred_data %>%
    mutate(surv = 1- pred_CIF) %>%
    group_by(age_exact) %>%
    summarise(surv = mean(surv)) %>%
    mutate(method = "Joint Model") %>%
    mutate(gender = gender)
    
  return(list(
    individual.predictions = individual.predictions[,c("shield.id", "age_exact", "surv", "method", "gender")],
    average.predictions = average.predictions[,c("age_exact", "surv", "method", "gender")]
  ))
}

generate_surv_curves <- function(test.data_female,
                            jm.pred_female,
                            baseline.fit_female,
                            time.dep.fit_female,
                            long.data_female,
                            
                            test.data_male,
                            jm.pred_male,
                            baseline.fit_male,
                            time.dep.fit_male,
                            long.data_male){
  
  long.data_male <- test_long_data_prep(long.data_male)
  long.data_female <- test_long_data_prep(long.data_female)
  
  jm.surv.pred_male <- prep_jm_surv_data(jm.pred_male, "Male")$average.predictions
  jm.surv.pred_female <- prep_jm_surv_data(jm.pred_female, "Female")$average.predictions
  
  baseline.surv.pred_male <- baseline_pred(baseline.fit_male, long.data_male, "Male")$average.predictions
  baseline.surv.pred_female <- baseline_pred(baseline.fit_female, long.data_female, "Female")$average.predictions
  
  time.dep.surv.pred_male <- time.dep_pred(time.dep.fit_male, long.data_male, "Male")$average.predictions
  time.dep.surv.pred_female <- time.dep_pred(time.dep.fit_female, long.data_female, "Female")$average.predictions
  
  rfsrs.surv.pred_male <- rfsrs_pred(long.data_male, "Male")$average.predictions
  rfsrs.surv.pred_female <- rfsrs_pred(long.data_female, "Female")$average.predictions
  
  obs.surv_male <- emp_surv_curve(test.data_male, end.age = 90, gender = "Male")
  obs.surv_female <- emp_surv_curve(test.data_female, end.age = 90, gender = "Female")
  
  isch.stroke.surv_data <- rbind(jm.surv.pred_female[,c("age_exact", "surv", "method", "gender")],
                                 jm.surv.pred_male[,c("age_exact", "surv", "method", "gender")],
                                 
                                 baseline.surv.pred_female[,c("age_exact", "surv", "method", "gender")],
                                 baseline.surv.pred_male[,c("age_exact", "surv", "method", "gender")],
                                 
                                 time.dep.surv.pred_female[,c("age_exact", "surv", "method", "gender")],
                                 time.dep.surv.pred_male[,c("age_exact", "surv", "method", "gender")],
                                 
                                 rfsrs.surv.pred_female[,c("age_exact", "surv", "method", "gender")],
                                 rfsrs.surv.pred_male[,c("age_exact", "surv", "method", "gender")],
                                 
                                 obs.surv_female$isch.stroke_surv.curve[,c("age_exact", "surv", "method", "gender")],
                                 obs.surv_male$isch.stroke_surv.curve[,c("age_exact", "surv", "method", "gender")])
  
  isch.stroke.surv_curve <- ggplot(data = isch.stroke.surv_data, aes(x = age_exact, y = surv)) +
    geom_line(aes(group = method, color = method), size = 1) +
    geom_ribbon(data = rbind(obs.surv_female$isch.stroke_surv.curve, obs.surv_male$isch.stroke_surv.curve), aes(x = age_exact, ymin = surv.lower, ymax = surv.upper), color = "grey", alpha = 0.1) +
    xlab("Age (years)") + ylab("Stroke-Free Probability") +
    theme_bw() + theme(panel.border = element_blank(),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       axis.line = element_line(colour = "black", size = 1, lineend = "square"),
                       axis.title = element_text(size = 13, family = "Franklin Gothic Medium"),
                       axis.text = element_text(size = 12, family = "Franklin Gothic Book"),
                       axis.ticks = element_blank(),
                       legend.text = element_text(size = 14, family = "Franklin Gothic Book"),
                       legend.title = element_text(size = 16, family = "Franklin Gothic Medium"),
                       plot.title = element_text(size = 16, family = "Franklin Gothic Medium"),
                       strip.text = element_text(size = 16, family = "Franklin Gothic Medium"),
                       strip.background = element_blank()) +
    scale_color_manual(name = "Estimation Method", breaks = c("Observed Stroke Survival", "Joint Model", "Baseline Model", "Time-Dependent Model", "Revised Framingham Stroke Risk Score"), values = c("black", "#00BA38", "#F8766D", "#619CFF", "purple"), c("Observed Stroke Survival", "Joint Model", "Baseline Model", "Time-Dependent Model", "Revised Framingham \nStroke Risk Score")) +
    facet_grid(cols = vars(gender))
  
  return(list(
    "isch.stroke.surv_noncr_curve" = isch.stroke.surv_curve
  ))
}

generate_calib_data <- function(test.data_female,
                                  jm.pred_female,
                                  baseline.fit_female,
                                  time.dep.fit_female,
                                  long.data_female,
                                  
                                  test.data_male,
                                  jm.pred_male,
                                  baseline.fit_male,
                                  time.dep.fit_male,
                                  long.data_male){
  
  long.data_male <- test_long_data_prep(long.data_male)
  long.data_female <- test_long_data_prep(long.data_female)
  
  end.ages <- seq(from = min(long.data_male$age_exact), to = max(long.data_male$age_exact) + 1e-7, by = 5)[-1]
  
  jm.surv.pred_male <- prep_jm_surv_data(jm.pred_male, "Male")$individual.predictions
  jm.surv.pred_female <- prep_jm_surv_data(jm.pred_female, "Female")$individual.predictions
  
  baseline.surv.pred_male <- baseline_pred(baseline.fit_male, long.data_male, "Male")$individual.predictions
  baseline.surv.pred_female <- baseline_pred(baseline.fit_female, long.data_female, "Female")$individual.predictions
  
  time.dep.surv.pred_male <- time.dep_pred(time.dep.fit_male, long.data_male, "Male")$individual.predictions
  time.dep.surv.pred_female <- time.dep_pred(time.dep.fit_female, long.data_female, "Female")$individual.predictions
  
  rfsrs.surv.pred_male <- rfsrs_pred(long.data_male, "Male")$individual.predictions
  rfsrs.surv.pred_female <- rfsrs_pred(long.data_female, "Female")$individual.predictions
  
  calib.results <- expand_grid(gender = c("Male", "Female"),
                               method = c("Joint Model", "Baseline Model", "Time-Dependent Model", "Revised Framingham Stroke Risk Score"),
                               age = end.ages) %>%
    mutate(E50 = NA_real_, E90 = NA_real_)
  
  for(gender in c("Male", "Female")){
    obs.data <- if(gender == "Male"){
      test.data_male$surv.data} else if(gender == "Female"){
        test.data_female$surv.data 
      }
    for(method in c("Joint Model", "Baseline Model", "Time-Dependent Model", "Revised Framingham Stroke Risk Score")){
      predictions <- if(method == "Joint Model" & gender == "Male"){
        jm.surv.pred_male
        } else if(method == "Joint Model" & gender == "Female"){
          jm.surv.pred_female
        } else if(method == "Baseline Model" & gender == "Male"){
          baseline.surv.pred_male
        } else if(method == "Baseline Model" & gender == "Female"){
          baseline.surv.pred_female
        } else if(method == "Time-Dependent Model" & gender == "Male"){
          time.dep.surv.pred_male
        } else if(method == "Time-Dependent Model" & gender == "Female"){
          time.dep.surv.pred_female
        } else if(method == "Revised Framingham Stroke Risk Score" & gender == "Male"){
          rfsrs.surv.pred_male
        } else if(method == "Revised Framingham Stroke Risk Score" & gender == "Female"){
          rfsrs.surv.pred_female
        }
      for(age in end.ages){
        
        print(paste("Working on:", gender, method, age))
        
        age.predictions <- predictions %>%
          filter(round(age_exact) == round(age)) %>%
          mutate(prob = 1 - surv)
        
        calib.data <- full_join(age.predictions, obs.data, by = "shield.id")
        bounds <- c(round(quantile(cloglog(calib.data$prob), 0.01), 2), round(quantile(cloglog(calib.data$prob), 0.99), 2))
        
        cox.model <- coxph(Surv(time = age_event, event = event_ischstroke_merged) ~ ns(cloglog(prob), df = 3, Boundary.knots = bounds), data = calib.data)
        
        obs_prob  <- 1 - c(summary(survfit(cox.model, newdata = calib.data), times = age)$surv)
        
        abs.diff <- abs(calib.data$prob - obs_prob)
        
        calib.results[calib.results$gender == gender & calib.results$method == method & calib.results$age == age, c("E50", "E90")] <- t(quantile(abs.diff, c(0.5, 0.9), names = FALSE))
      }
    } 
  }
  
  return(calib.results)
}

generate_calib_curve <- function(calib_results){
  calib_curve <- ggplot(data = calib_results, aes(x = age)) +
    geom_line(aes(y = E90, group = method, color = method, linetype = "90th Percentile Difference"), size = 1) +
    geom_line(aes(y = E50, group = method, color = method, linetype = "Median Difference"), size = 1) +
    xlab("Age (years)") + ylab("Absolute Difference between Estimated and Observed Stroke Probability") + scale_linetype_discrete("Calibration Statistic") +
    theme_bw() + theme(panel.border = element_blank(),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       axis.line = element_line(colour = "black", size = 1, lineend = "square"),
                       axis.title = element_text(size = 12, family = "Franklin Gothic Medium"),
                       axis.text = element_text(size = 12, family = "Franklin Gothic Book"),
                       axis.ticks = element_blank(),
                       legend.text = element_text(size = 14, family = "Franklin Gothic Book"),
                       legend.title = element_text(size = 16, family = "Franklin Gothic Medium"),
                       plot.title = element_text(size = 16, family = "Franklin Gothic Medium"),
                       strip.text = element_text(size = 16, family = "Franklin Gothic Medium"),
                       strip.background = element_blank())  +
    scale_color_manual(name = "Estimation Method", breaks = c("Joint Model", "Baseline Model", "Time-Dependent Model", "Revised Framingham Stroke Risk Score"), values = c("#00BA38", "#F8766D", "#619CFF", "purple"), labels = c("Joint Model", "Baseline Model", "Time-Dependent Model", "Revised Framingham \nStroke Risk Score")) +
    facet_grid(cols = vars(gender))
  
  return(calib_curve)
}
