library(dplyr)
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

data_prep <- function(data, gender){
  gender_num <- case_when(gender == "female" ~ 1,
                          gender == "male" ~ 0)
  ids <- unique(data$shield.id[data$female_merged == gender_num])
  training.set <- sample(ids, length(ids)*0.5)
  test.set <- setdiff(ids, training.set)
  training.data_long <- data[data$shield.id %in% training.set,]
  training.data_event <- unique(training.data_long[,c("shield.id", "age_baseline", "age_event", "event_cause")])
  training.data_cr <- crisk_setup(data = training.data_event, statusVar = "event_cause", censLevel = "alive", nameStrata = "comp.risk")
  
  return(list(
    gender = gender,
    training.ids = training.set,
    test.ids = test.set,
    training.data_long = training.data_long,
    training.data_surv = training.data_cr
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
  
  cox_model <- coxph(Surv(time = age_baseline, time2 = age_event, event = status2) ~ strata(comp.risk), data = surv_data)
  
  return(list(
    long_models = long_models,
    surv_model = cox_model,
    surv_data = surv_data
  ))
}

jm_run <- function(model.list){
  
  fx.forms <- list(
    "cursmk" = ~ value(cursmk):comp.risk,
    "prior_cvd_merged" = ~ value(prior_cvd_merged):comp.risk,
    "afib_previous_merged" = ~ value(afib_previous_merged):comp.risk,
    "dm_previous_merged" = ~ value(dm_previous_merged):comp.risk,
    "htnmed_merged" = ~ value(htnmed_merged):comp.risk,
    "sbp" = ~ value(sbp):comp.risk
  )
  
  jm_fit <- jm(model.list$surv_model, model.list$long_models, time_var = "age_exact", functional_forms = fx.forms, data_Surv = model.list$surv_data, control = list(n_burnin = 5000L, n_iter = 25000L, n_chains = 5, Bsplines_degree = 3, cores = 5))
  
  return(jm_fit)
  
}

test.data_prep <- function(test.list, data, start.age){
  
  data <- data %>%
    group_by(shield.id) %>%
    mutate(age_baseline = min(age_merged)) %>%
    filter(ttoevent_ischstroke_merged != 0)
  
  data <- data[complete.cases(data[,c("ttoevent_ischstroke_merged", "event_ischstroke_merged", "ttoevent.death", "event.death", "cursmk", "prior_cvd_merged", "afib_previous_merged", "dm_previous_merged", "htnmed_merged", "sbp")]),] %>%
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
  
  test.data_pre.start <- data[data$shield.id %in% test.list & data$age_event > start.age & data$age_exact <= start.age,]
  test.data_post.start <- data[data$shield.id %in% test.list & data$age_event > start.age & data$age_exact > start.age,]
  
  test.data_surv <- unique(test.data_pre.start[,c("shield.id", "age_baseline", "age_event", "event_cause")])
  test.data_cr <- crisk_setup(data = test.data_surv, statusVar = "event_cause", censLevel = "alive", nameStrata = "comp.risk")
  test.data_cr$status2 <- 0
  test.data_cr$age_event <- start.age + 1e-7
  
  test.data_surv <- test.data_surv %>%
    mutate(event_cause_factor = case_when(event_cause == "alive" ~ 0,
                                          event_cause == "isch.stroke" ~ 1,
                                          event_cause == "death" ~ 2)) %>%
    mutate(event_cause_factor = factor(event_cause_factor, 0:2, labels = c("alive", "isch.stroke", "death")))
  
  return(list(
    pre.data = test.data_pre.start,
    post.data = test.data_post.start,
    surv.data = test.data_surv,
    jm.pred.data = list(newdataL = test.data_pre.start, newdataE = test.data_cr)
  ))
}

jm_pred_loop <- function(jm.fit.object, test.data, start.age, end.age){
  test.ids <- unique(test.data$newdataE$shield.id)
  partitions <- ceiling(seq_along(test.ids)/100)
  test.partitions <- split(test.ids, partitions)
  
  jm_pred <- foreach(group=names(test.partitions), .combine = rbind) %do% {
    
    data <- list("newdataL" = test.data$newdataL[test.data$newdataL$shield.id %in% test.partitions[[group]],], "newdataE" = test.data$newdataE[test.data$newdataE$shield.id %in% test.partitions[[group]],])
    
    jm_pred_run <- predict(jm.fit.object, newdata = data, process = "event", times = seq(from = start.age, to = end.age, by = 1), return_newdata = TRUE, return_mcmc = FALSE, cores = 8)
    
    return(jm_pred_run)
  }
  return(jm_pred)
}

jm_pred_long_loop <- function(jm.fit.object, test.data, start.age, end.age){
  test.ids <- unique(test.data$newdataE$shield.id)
  partitions <- ceiling(seq_along(test.ids)/100)
  test.partitions <- split(test.ids, partitions)
  
  jm_pred <- foreach(group=names(test.partitions), .combine = rbind) %do% {
    
    data <- list("newdataL" = test.data$newdataL[test.data$newdataL$shield.id %in% test.partitions[[group]],], "newdataE" = test.data$newdataE[test.data$newdataE$shield.id %in% test.partitions[[group]],])
    
    jm_pred_run <- predict(jm.fit.object, newdata = data, process = "longitudinal", type_pred = "response", times = seq(from = start.age, to = end.age, by = 1), return_newdata = TRUE, return_mcmc = FALSE, cores = 8)$newdata2
    
    return(jm_pred_run)
  }
  return(jm_pred)
}

emp_surv_curve_cr <- function(test.data, end.age, gender){
  surv.data <- test.data$surv.data
  curve.object <- survfit(Surv(time = age_event, event = event_cause_factor) ~ 1, data = surv.data, conf.type = "log-log")
  
  isch.stroke.curve <- data.frame(age_exact = curve.object$time, surv = 1 - curve.object$pstate[,2], surv.lower = 1 - curve.object$upper[,2], surv.upper = 1 - curve.object$lower[,2]) %>%
    mutate(method = "Observed Stroke-Free Survival") %>%
    mutate(gender = gender) %>%
    filter(age_exact <= end.age)
  
  death.curve <- data.frame(age_exact = curve.object$time, surv = 1 - curve.object$pstate[,3], surv.lower = 1 - curve.object$upper[,3], surv.upper = 1 - curve.object$lower[,3]) %>%
    mutate(method = "Observed Pre-Stroke Death Survival") %>%
    mutate(gender = gender) %>%
    filter(age_exact <= end.age)
  
  return(list(
    isch.stroke_surv.curve = isch.stroke.curve[,c("age_exact", "surv", "method", "gender", "surv.lower", "surv.upper")],
    death_surv.curve = death.curve[,c("age_exact", "surv", "method", "gender", "surv.lower", "surv.upper")]
  ))
}

emp_hazard_curve_cr <- function(test.data, start.age, end.age, gender){
  surv.data <- test.data$surv.data
  surv.data_isch.stroke <- surv.data %>% mutate(event = case_when(
    event_cause == "alive" | event_cause == "death" ~ 0,
    event_cause == "isch.stroke" ~ 1
  ))
  surv.data_death <- surv.data %>% mutate(event = case_when(
    event_cause == "alive" | event_cause == "isch.stroke" ~ 0,
    event_cause == "death" ~ 1
  ))
  
  #curve.object_isch.stroke <- bshazard(Surv(time = age_event, event = event) ~ 1, data = surv.data_isch.stroke, nk = 10, degree = 3)
  #curve.object_death <- bshazard(Surv(time = age_event, event = event) ~ 1, data = surv.data_death, nk = 10, degree = 3)
  
  curve.object_isch.stroke <- muhaz(times = surv.data_isch.stroke$age_event, delta = surv.data_isch.stroke$event, min.time = start.age, max.time = end.age, n.est.grid = as.integer(end.age - start.age), bw.method = "local")
  curve.object_death <- muhaz(times = surv.data_death$age_event, delta = surv.data_death$event, min.time = start.age, max.time = end.age, n.est.grid = as.integer(end.age - start.age), bw.method = "local")
  
  # isch.stroke.curve <- data.frame(age_exact = curve.object_isch.stroke$time, hazard = curve.object_isch.stroke$hazard, hazard.lower = curve.object_isch.stroke$lower.ci, hazard.upper = curve.object_isch.stroke$upper.ci) %>%
  #   mutate(method = "Observed Stroke Hazard") %>%
  #   mutate(gender = gender) %>%
  #   filter(age_exact <= end.age & age_exact >= start.age + 1)
  
  isch.stroke.curve <- data.frame(age_exact = curve.object_isch.stroke$est.grid, hazard = curve.object_isch.stroke$haz.est) %>%
    mutate(method = "Observed Stroke Hazard") %>%
    mutate(gender = gender) %>%
    filter(age_exact <= end.age & age_exact >= start.age + 1)
  
  # death.curve <- data.frame(age_exact = curve.object_death$time, hazard = curve.object_death$hazard, hazard.lower = curve.object_death$lower.ci, hazard.upper = curve.object_death$upper.ci) %>%
  #   mutate(method = "Observed Mortality Hazard") %>%
  #   mutate(gender = gender) %>%
  #   filter(age_exact <= end.age & age_exact >= start.age + 1)
  
  death.curve <- data.frame(age_exact = curve.object_death$est.grid, hazard = curve.object_death$haz.est) %>%
    mutate(method = "Observed Mortality Hazard") %>%
    mutate(gender = gender) %>%
    filter(age_exact <= end.age & age_exact >= start.age + 1)
  
  # return(list(
  #   isch.stroke_hazard.curve = isch.stroke.curve[,c("age_exact", "hazard", "method", "gender", "hazard.lower", "hazard.upper")],
  #   death_hazard.curve = death.curve[,c("age_exact", "hazard", "method", "gender", "hazard.lower", "hazard.upper")]
  # ))
  
  return(list(
    isch.stroke_hazard.curve = isch.stroke.curve[,c("age_exact", "hazard", "method", "gender")],
    death_hazard.curve = death.curve[,c("age_exact", "hazard", "method", "gender")]
  ))
}

prep_jm_surv_data <- function(jm_pred_data, gender){
  
  isch.stroke.curve <- jm_pred_data %>%
    filter(comp.risk == "isch.stroke") %>%
    mutate(surv = 1- pred_CIF) %>%
    group_by(age_exact) %>%
    summarise(surv = mean(surv)) %>%
    mutate(method = "Joint Model Stroke-Free Survival") %>%
    mutate(gender = gender)
    
  death.curve <- jm_pred_data %>%
    filter(comp.risk == "death") %>%
    mutate(surv = 1- pred_CIF) %>%
    group_by(age_exact) %>%
    summarise(surv = mean(surv)) %>%
    mutate(method = "Joint Model Pre-Stroke Death Survival") %>%
    mutate(gender = gender)
    
  return(list(
    isch.stroke_surv.curve = isch.stroke.curve[,c("age_exact", "surv", "method", "gender")],
    death_surv.curve = death.curve[,c("age_exact", "surv", "method", "gender")]
  ))
}

generate_cr_curves <- function(test.data_female, jm.pred_female, test.data_male, jm.pred_male){
  
  jm.surv.pred_male <- prep_jm_surv_data(jm.pred_male, "Male")
  jm.surv.pred_female <- prep_jm_surv_data(jm.pred_female, "Female")
  
  obs.surv_male <- emp_surv_curve_cr(test.data_male, end.age = 90, gender = "Male")
  obs.surv_female <- emp_surv_curve_cr(test.data_female, end.age = 90, gender = "Female")
  
  isch.stroke.surv_data <- rbind(jm.surv.pred_female$isch.stroke_surv.curve[,c("age_exact", "surv", "method", "gender")], jm.surv.pred_male$isch.stroke_surv.curve[,c("age_exact", "surv", "method", "gender")], obs.surv_female$isch.stroke_surv.curve[,c("age_exact", "surv", "method", "gender")], obs.surv_male$isch.stroke_surv.curve[,c("age_exact", "surv", "method", "gender")])
  
  isch.stroke.surv_cr_curve <- ggplot(data = isch.stroke.surv_data, aes(x = age_exact, y = surv)) +
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
    scale_color_manual(name = "Estimation Method", breaks = c("Observed Stroke-Free Survival", "Joint Model Stroke-Free Survival"), values = c("black", "#00BA38"), c("Observed Stroke-Free Survival", "Joint Model Stroke-Free Survival")) +
    facet_grid(cols = vars(gender))
  
  death.surv_data <- rbind(jm.surv.pred_female$death_surv.curve[,c("age_exact", "surv", "method", "gender")], jm.surv.pred_male$death_surv.curve[,c("age_exact", "surv", "method", "gender")], obs.surv_female$death_surv.curve[,c("age_exact", "surv", "method", "gender")], obs.surv_male$death_surv.curve[,c("age_exact", "surv", "method", "gender")])
  
  death.surv_cr_curve <- ggplot(data = death.surv_data, aes(x = age_exact, y = surv)) +
    geom_line(aes(group = method, color = method), size = 1) +
    geom_ribbon(data = rbind(obs.surv_female$death_surv.curve, obs.surv_male$death_surv.curve), aes(x = age_exact, ymin = surv.lower, ymax = surv.upper), color = "grey", alpha = 0.1) +
    xlab("Age (years)") + ylab("Survival Probability") +
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
    scale_color_manual(name = "Estimation Method", breaks = c("Observed Pre-Stroke Death Survival", "Joint Model Pre-Stroke Death Survival"), values = c("black", "#00BA38"), c("Observed Pre-Stroke Death Survival", "Joint Model Pre-Stroke Death Survival")) +
    facet_grid(cols = vars(gender))
  
  return(list(
    "isch.stroke.surv_cr_curve" = isch.stroke.surv_cr_curve,
    "death.surv_cr_curve" = death.surv_cr_curve
  ))
}

generate_cr_hazard_curves <- function(test.data_female, jm.pred_female, test.data_male, jm.pred_male){
  
  jm.surv.pred_male <- prep_jm_surv_data(jm.pred_male, "Male")
  jm.surv.pred_female <- prep_jm_surv_data(jm.pred_female, "Female")
  
  obs.hazard_male <- emp_hazard_curve_cr(test.data_male, start.age = min(jm.surv.pred_male$isch.stroke_surv.curve$age_exact), end.age = 90, gender = "Male")
  obs.hazard_female <- emp_hazard_curve_cr(test.data_female, start.age = min(jm.surv.pred_female$isch.stroke_surv.curve$age_exact), end.age = 90, gender = "Female")
  
  jm.hazard.pred_male <- lapply(jm.surv.pred_male, function(data){data %>% mutate(cumhaz = -log(surv)) %>% mutate(hazard = cumhaz - lag(cumhaz), method = ifelse(method == "Joint Model Stroke-Free Survival", "Joint Model Stroke Hazard", "Joint Model Mortality Hazard")) %>% filter(!is.na(hazard)) %>% mutate(hazard = (hazard - lag(hazard, default = hazard[1])) + hazard)})
  jm.hazard.pred_female <- lapply(jm.surv.pred_female, function(data){data %>% mutate(cumhaz = -log(surv)) %>% mutate(hazard = cumhaz - lag(cumhaz), method = ifelse(method == "Joint Model Stroke-Free Survival", "Joint Model Stroke Hazard", "Joint Model Mortality Hazard")) %>% filter(!is.na(hazard)) %>% mutate(hazard = (hazard - lag(hazard, default = hazard[1])) + hazard)})
  
  isch.stroke.hazard_data <- rbind(jm.hazard.pred_female$isch.stroke_surv.curve[,c("age_exact", "hazard", "method", "gender")], jm.hazard.pred_male$isch.stroke_surv.curve[,c("age_exact", "hazard", "method", "gender")], obs.hazard_female$isch.stroke_hazard.curve[,c("age_exact", "hazard", "method", "gender")], obs.hazard_male$isch.stroke_hazard.curve[,c("age_exact", "hazard", "method", "gender")])
  
  isch.stroke.hazard_cr_curve <- ggplot(data = isch.stroke.hazard_data, aes(x = age_exact, y = hazard)) +
    geom_line(aes(group = method, color = method), size = 1) +
    #geom_ribbon(data = rbind(obs.hazard_female$isch.stroke_hazard.curve, obs.hazard_male$isch.stroke_hazard.curve), aes(x = age_exact, ymin = hazard.lower, ymax = hazard.upper), color = "grey", alpha = 0.1) +
    xlab("Age (years)") + ylab("Stroke Hazard") +
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
    scale_color_manual(name = "Estimation Method", breaks = c("Observed Stroke Hazard", "Joint Model Stroke Hazard"), values = c("black", "#00BA38"), c("Observed Stroke Hazard", "Joint Model Stroke Hazard")) +
    facet_grid(cols = vars(gender))
  
  isch.stroke.hazard_cr_curve_female <- ggplot(data = filter(isch.stroke.hazard_data, gender == "Female"), aes(x = age_exact, y = hazard)) +
    geom_line(aes(group = method, color = method), size = 1) +
    #geom_ribbon(data = rbind(obs.hazard_female$isch.stroke_hazard.curve, obs.hazard_male$isch.stroke_hazard.curve), aes(x = age_exact, ymin = hazard.lower, ymax = hazard.upper), color = "grey", alpha = 0.1) +
    xlab("Age (years)") + ylab("Stroke Hazard") +
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
    scale_color_manual(name = "Estimation Method", breaks = c("Observed Stroke Hazard", "Joint Model Stroke Hazard"), values = c("black", "#00BA38"), c("Observed Stroke Hazard", "Joint Model Stroke Hazard")) + ggtitle("Female")
  
  death.hazard_data <- rbind(jm.hazard.pred_female$death_surv.curve[,c("age_exact", "hazard", "method", "gender")], jm.hazard.pred_male$death_surv.curve[,c("age_exact", "hazard", "method", "gender")], obs.hazard_female$death_hazard.curve[,c("age_exact", "hazard", "method", "gender")], obs.hazard_male$death_hazard.curve[,c("age_exact", "hazard", "method", "gender")])
  
  death.hazard_cr_curve <- ggplot(data = death.hazard_data, aes(x = age_exact, y = hazard)) +
    geom_line(aes(group = method, color = method), size = 1) +
    #geom_ribbon(data = rbind(obs.hazard_female$death_hazard.curve, obs.hazard_male$death_hazard.curve), aes(x = age_exact, ymin = hazard.lower, ymax = hazard.upper), color = "grey", alpha = 0.1) +
    xlab("Age (years)") + ylab("Stroke Hazard") +
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
    scale_color_manual(name = "Estimation Method", breaks = c("Observed Mortality Hazard", "Joint Model Mortality Hazard"), values = c("black", "#00BA38"), c("Observed Mortality Hazard", "Joint Model Mortality Hazard")) +
    facet_grid(cols = vars(gender))
  
  return(list(
    "isch.stroke.hazard_cr_curve" = isch.stroke.hazard_cr_curve,
    "death.hazard_cr_curve" = death.hazard_cr_curve,
    "isch.stroke.hazard_cr_curve_female" = isch.stroke.hazard_cr_curve_female
  ))
}
