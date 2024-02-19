## title: "Building Sample with NHANES"
## author: "Astrid Yu, Nathaniel Alemayehu, Jinyi Zhu"

## Load packages
  library(tidyverse)
  library(janitor)
  library(dplyr)
  library(readxl)
  library(haven)
  library(mice)

# Read SAS -----------------------------------------------------------------------

## Set working directory for output

## Read SAS files and pull relevant columns into dataframe
  # 2017-2018
  nhanes_data_1 = data.frame(read_xpt("NHANES SAS data/BMX_J.XPT", col_select = c("SEQN","BMXBMI")))
  nhanes_data_2 = data.frame(read_xpt("NHANES SAS data/GHB_J.XPT", col_select = c("SEQN","LBXGH")))
  nhanes_data_3 = data.frame(read_xpt("NHANES SAS data/BPX_J.XPT", col_select = c("SEQN","BPXDI1","BPXSY1")))
  nhanes_data_4 = data.frame(read_xpt("NHANES SAS data/DEMO_J.XPT", col_select = c("SEQN","RIDAGEYR","RIAGENDR","RIDRETH1","WTMEC2YR"))) #WTMEC2YR is the sample weights variable
  nhanes_data_5 = data.frame(read_xpt("NHANES SAS data/DIQ_J.XPT", col_select = c("SEQN","DIQ010")))
  nhanes_data_6 = data.frame(read_xpt("NHANES SAS data/HDL_J.XPT", col_select = c("SEQN","LBDHDD")))
  nhanes_data_7 = data.frame(read_xpt("NHANES SAS data/MCQ_J.XPT", col_select = c("SEQN","MCQ160D","MCQ160F","MCQ160E")))
  nhanes_data_8 = data.frame(read_xpt("NHANES SAS data/SMQ_J.XPT", col_select = c("SEQN","SMQ040")))
  nhanes_data_9 = data.frame(read_xpt("NHANES SAS data/TCHOL_J.XPT", col_select = c("SEQN","LBXTC")))
  nhanes_data_11 = data.frame(read_xpt("NHANES SAS data/BPQ_J.XPT", col_select = c("SEQN", "BPQ050A")))
  nhanes_data_12 = data.frame(read_xpt("NHANES SAS data/TRIGLY_J.XPT", col_select = c("SEQN", "LBDLDL")))
  
  # 2013-2014
  nhanes_data_1A = data.frame(read_xpt("NHANES SAS data/BMX_H.XPT", col_select = c("SEQN","BMXBMI")))
  nhanes_data_2A = data.frame(read_xpt("NHANES SAS data/GHB_H.XPT", col_select = c("SEQN","LBXGH")))
  nhanes_data_3A = data.frame(read_xpt("NHANES SAS data/BPX_H.XPT", col_select = c("SEQN","BPXDI1","BPXSY1")))
  nhanes_data_4A = data.frame(read_xpt("NHANES SAS data/DEMO_H.XPT", col_select = c("SEQN","RIDAGEYR","RIAGENDR","RIDRETH1","WTMEC2YR"))) #WTMEC2YR is the sample weights variable
  nhanes_data_5A = data.frame(read_xpt("NHANES SAS data/DIQ_H.XPT", col_select = c("SEQN","DIQ010")))
  nhanes_data_6A = data.frame(read_xpt("NHANES SAS data/HDL_H.XPT", col_select = c("SEQN","LBDHDD")))
  nhanes_data_7A = data.frame(read_xpt("NHANES SAS data/MCQ_H.XPT", col_select = c("SEQN","MCQ160D","MCQ160F","MCQ160E")))
  nhanes_data_8A = data.frame(read_xpt("NHANES SAS data/SMQ_H.XPT", col_select = c("SEQN","SMQ040")))
  nhanes_data_9A = data.frame(read_xpt("NHANES SAS data/TCHOL_H.XPT", col_select = c("SEQN","LBXTC")))
  nhanes_data_11A = data.frame(read_xpt("NHANES SAS data/BPQ_H.XPT", col_select = c("SEQN", "BPQ050A")))
  nhanes_data_12A = data.frame(read_xpt("NHANES SAS data/TRIGLY_H.XPT", col_select = c("SEQN", "LBDLDL")))
  
  # 2015-2016
  nhanes_data_1B = data.frame(read_xpt("NHANES SAS data/BMX_I.XPT", col_select = c("SEQN","BMXBMI")))
  nhanes_data_2B = data.frame(read_xpt("NHANES SAS data/GHB_I.XPT", col_select = c("SEQN","LBXGH")))
  nhanes_data_3B = data.frame(read_xpt("NHANES SAS data/BPX_I.XPT", col_select = c("SEQN","BPXDI1","BPXSY1")))
  nhanes_data_4B = data.frame(read_xpt("NHANES SAS data/DEMO_I.XPT", col_select = c("SEQN","RIDAGEYR","RIAGENDR","RIDRETH1","WTMEC2YR"))) #WTMEC2YR is the sample weights variable
  nhanes_data_5B = data.frame(read_xpt("NHANES SAS data/DIQ_I.XPT", col_select = c("SEQN","DIQ010")))
  nhanes_data_6B = data.frame(read_xpt("NHANES SAS data/HDL_I.XPT", col_select = c("SEQN","LBDHDD")))
  nhanes_data_7B = data.frame(read_xpt("NHANES SAS data/MCQ_I.XPT", col_select = c("SEQN","MCQ160D","MCQ160F","MCQ160E")))
  nhanes_data_8B = data.frame(read_xpt("NHANES SAS data/SMQ_I.XPT", col_select = c("SEQN","SMQ040")))
  nhanes_data_9B = data.frame(read_xpt("NHANES SAS data/TCHOL_I.XPT", col_select = c("SEQN","LBXTC")))
  nhanes_data_11B = data.frame(read_xpt("NHANES SAS data/BPQ_I.XPT", col_select = c("SEQN", "BPQ050A")))
  nhanes_data_12B = data.frame(read_xpt("NHANES SAS data/TRIGLY_I.XPT", col_select = c("SEQN", "LBDLDL")))
  
## Merge dataframes and rename columns
  nhanes_raw_data_C <- list(nhanes_data_1,nhanes_data_2,nhanes_data_3,nhanes_data_4,
                          nhanes_data_5,nhanes_data_6,nhanes_data_7,nhanes_data_8,nhanes_data_9, nhanes_data_11, nhanes_data_12)  
  nhanes_raw_data_A <- list(nhanes_data_1A,nhanes_data_2A,nhanes_data_3A,nhanes_data_4A,
                          nhanes_data_5A,nhanes_data_6A,nhanes_data_7A,nhanes_data_8A,nhanes_data_9A, nhanes_data_11A, nhanes_data_12A) 
  nhanes_raw_data_B <- list(nhanes_data_1B,nhanes_data_2B,nhanes_data_3B,nhanes_data_4B,
                          nhanes_data_5B,nhanes_data_6B,nhanes_data_7B,nhanes_data_8B,nhanes_data_9B, nhanes_data_11B, nhanes_data_12B) 
  
  nhanes_raw_data_C <- data.frame(reduce(nhanes_raw_data_C, full_join, by='SEQN'))
  nhanes_raw_data_A <- data.frame(reduce(nhanes_raw_data_A, full_join, by='SEQN'))
  nhanes_raw_data_B <- data.frame(reduce(nhanes_raw_data_B, full_join, by='SEQN'))
  
  nhanes_raw_data <- rbind(nhanes_raw_data_C, nhanes_raw_data_A, nhanes_raw_data_B)


# Clearance, Imputation and Filter -------------------------------------------------------

## Cleaning Before Imputing, ifelse won't change NAs
  # matching race to CHS categories
  nhanes_raw_data$race <- ifelse(nhanes_raw_data$RIDRETH1 == 4, 2, 1)
  # compiling vascular history variable
  nhanes_raw_data$vascular_history <- ifelse(nhanes_raw_data$MCQ160D == 1 | nhanes_raw_data$MCQ160F == 1 | nhanes_raw_data$MCQ160E == 1, 1, 0)
  # compiling diabetes rollup variable
  nhanes_raw_data$diabetes_rollup <- ifelse(nhanes_raw_data$DIQ010 == 1 | nhanes_raw_data$LBXGH > 6.5, 1, 0)
  # compile present smoker variable to match CHS categories
  nhanes_raw_data$SMQ040 <- ifelse(nhanes_raw_data$SMQ040 == 1 | nhanes_raw_data$SMQ040 == 2, 1, 0)
  # convert gender to factor variable(originally: male - 1, female - 2), same as CHS
  nhanes_raw_data$RIAGENDR <- nhanes_raw_data$RIAGENDR - 1
  nhanes_raw_data$male = ifelse(nhanes_raw_data$RIAGENDR == 1, 0, 1)
  # nhanes_raw_data$SBPtxt = as.factor(ifelse(nhanes_raw_data$RXDRSC1 == "I10" | nhanes_raw_data$RXDRSC1 == "I10.P", 1, 0))
  nhanes_raw_data$SBPtxt = ifelse(nhanes_raw_data$BPQ050A == 1, 1, 0)
  
  ## Drop variables
  nhanes_raw_data = nhanes_raw_data %>% select(-LBXGH, -DIQ010, -MCQ160D, -MCQ160F, -MCQ160E, -RIDRETH1, -BPQ050A)
  
## Imputing data
  
  ## Calculating missingness in NHANES data
  missingness = (colMeans(is.na(nhanes_raw_data)))*100
  
  init = mice(nhanes_raw_data, maxit=0) 
  predM = init$predictorMatrix
  predM[,"SEQN"]=0
  nhanes_imputed_test <- mice(nhanes_raw_data, predictorMatrix = predM, m = 5, seed = 11)
  nhanes_raw_data <- complete(nhanes_imputed_test)
  
  ## Filter the data to only respondents over 40 years old and < 80 yrs old
  ## NHANES age: 0 to 79, 80 and over
  nhanes_raw_data <- nhanes_raw_data[which(nhanes_raw_data$RIDAGEYR >= 40 
                                           & nhanes_raw_data$RIDAGEYR <= 79
                                           & nhanes_raw_data$LBXTC >= 30 
                                           & nhanes_raw_data$LBXTC <= 500 
                                           & nhanes_raw_data$LBDHDD >= 5 
                                           & nhanes_raw_data$LBDHDD <= 200 
                                           & nhanes_raw_data$BPXSY1 >= 60 
                                           & nhanes_raw_data$BPXSY1 <= 200
                                           & nhanes_raw_data$vascular_history == 0), ]
  
  # check sample size now
  nrow(nhanes_raw_data)
  
# Pooled Cohort Equation --------------------------------------------------
  
  source('pcr.R')
  
  nhanes_raw_data$pce_orig <- sapply(1:nrow(nhanes_raw_data), function(r) {pcr(gender = ifelse(nhanes_raw_data$male[r] == 1, 'M', 'F'),
                   age = nhanes_raw_data$RIDAGEYR[r],
                   race = ifelse(nhanes_raw_data$race[r] == 2, 'African American', 'Other'), 
                   tot_chol = nhanes_raw_data$LBXTC[r],
                   hdl_chol = nhanes_raw_data$LBDHDD[r],
                   systolic_bp = nhanes_raw_data$BPXSY1[r],
                   bp_treatment = (nhanes_raw_data$SBPtxt[r] == 1),
                   smoker = (nhanes_raw_data$SMQ040[r] == 1),
                   diabetic = (nhanes_raw_data$diabetes_rollup[r] == 1),
                   prs_z = 0)[[1]]/100})
  
#convert race back to a 1/0 variable
  nhanes_raw_data$race <- ifelse(nhanes_raw_data$race == 1, 0, 1) 


# Sample Subgroups --------------------------------------------------------

  sub_dat = function(raw_data = nhanes_raw_data, n_pop = 1000000, PCE_lower = 0, PCE_upper = 1, statin_use_change = 1, diabetes_chol_exclude = 1, statin_use_threshold = 0.075){
    # target group not include people who don't need PCE score.
    if (diabetes_chol_exclude == 1){raw_data = subset(raw_data, LBDLDL < 190 & diabetes_rollup != 1)}
    
    # subgroups by PCE
    raw_data = subset(raw_data, !is.na(WTMEC2YR) & pce_orig >= PCE_lower & pce_orig < PCE_upper)
    
    # define statin_use_change
    statement <- "((tmp$pce_prs < statin_use_threshold & tmp$pce_orig >= statin_use_threshold) | (tmp$pce_prs >= statin_use_threshold & tmp$pce_orig < statin_use_threshold))"
    
    # Sample with replacement
    n_sub = 0
    sample = data.frame()
    while(n_sub < n_pop){
      # create a new subset
      tmp = sample_n(raw_data, size = n_pop, replace = TRUE, weight = WTMEC2YR)
      tmp$PRS = rnorm(nrow(tmp), 0, 1)
      tmp$pce_prs = 1 - exp(log(1 - tmp$pce_orig) * 1.73^tmp$PRS)
      
      if (statin_use_change == 1){tmp = subset(tmp, eval(parse(text = statement)))
      } else {tmp = subset(tmp, eval(parse(text = paste0("!", statement))))}
      
      # extend the sample and recalculate the indeicator
      sample = rbind(sample, tmp)
      n_sub = nrow(sample)
    }
    
    sample = sample_n(sample, size = n_pop, replace = FALSE)
    
    return(sample)
  }
  
  n_per_group = 1000000
  nhanes_sample_lower_5_change = sub_dat(n_pop = n_per_group, PCE_lower = 0, PCE_upper = 0.05)
  nhanes_sample_lower_5_unchange = sub_dat(n_pop = n_per_group, PCE_lower = 0, PCE_upper = 0.05, statin_use_change = 0)
  nhanes_sample_5_7.5_change = sub_dat(n_pop = n_per_group, PCE_lower = 0.05, PCE_upper = 0.075)
  nhanes_sample_5_7.5_unchange = sub_dat(n_pop = n_per_group, PCE_lower = 0.05, PCE_upper = 0.075, statin_use_change = 0)
  nhanes_sample_7.5_20_change = sub_dat(n_pop = n_per_group, PCE_lower = 0.075, PCE_upper = 0.20)
  nhanes_sample_7.5_20_unchange = sub_dat(n_pop = n_per_group, PCE_lower = 0.075, PCE_upper = 0.20, statin_use_change = 0)
  nhanes_sample_greater_20_change = sub_dat(n_pop = n_per_group, PCE_lower = 0.20, PCE_upper = 1)
  nhanes_sample_greater_20_unchange = sub_dat(n_pop = n_per_group, PCE_lower = 0.20, PCE_upper = 1, statin_use_change = 0)
  
# Save baseline pop
  save(nhanes_raw_data, 
       nhanes_sample_lower_5_change, 
       nhanes_sample_lower_5_unchange, 
       nhanes_sample_5_7.5_change, 
       nhanes_sample_5_7.5_unchange, 
       nhanes_sample_7.5_20_change, 
       nhanes_sample_7.5_20_unchange, 
       nhanes_sample_greater_20_change, 
       nhanes_sample_greater_20_unchange, 
       file = "basepop.RData")
  