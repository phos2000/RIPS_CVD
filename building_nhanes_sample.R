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
# setwd("./NHANES SAS data")

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
  # selecting only antihypertensive medications from prescription medications file
  nhanes_data_10 = data.frame(read_xpt("NHANES SAS data/RXQ_RX_J.XPT", col_select = c("SEQN","RXDRSC1")))
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
  # selecting only antihypertensive medications from prescription medications file
  nhanes_data_10A = data.frame(read_xpt("NHANES SAS data/RXQ_RX_H.XPT", col_select = c("SEQN","RXDRSC1")))
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
  # selecting only antihypertensive medications from prescription medications file
  nhanes_data_10B = data.frame(read_xpt("NHANES SAS data/RXQ_RX_I.XPT", col_select = c("SEQN","RXDRSC1")))
  nhanes_data_11B = data.frame(read_xpt("NHANES SAS data/BPQ_I.XPT", col_select = c("SEQN", "BPQ050A")))
  nhanes_data_12B = data.frame(read_xpt("NHANES SAS data/TRIGLY_I.XPT", col_select = c("SEQN", "LBDLDL")))
  
## Merge dataframes and rename columns
  nhanes_raw_data_C <- list(nhanes_data_1,nhanes_data_2,nhanes_data_3,nhanes_data_4,
                          nhanes_data_5,nhanes_data_6,nhanes_data_7,nhanes_data_8,nhanes_data_9,nhanes_data_10, nhanes_data_11, nhanes_data_12)  
  nhanes_raw_data_A <- list(nhanes_data_1A,nhanes_data_2A,nhanes_data_3A,nhanes_data_4A,
                          nhanes_data_5A,nhanes_data_6A,nhanes_data_7A,nhanes_data_8A,nhanes_data_9A,nhanes_data_10A, nhanes_data_11A, nhanes_data_12A) 
  nhanes_raw_data_B <- list(nhanes_data_1B,nhanes_data_2B,nhanes_data_3B,nhanes_data_4B,
                          nhanes_data_5B,nhanes_data_6B,nhanes_data_7B,nhanes_data_8,nhanes_data_9B,nhanes_data_10B, nhanes_data_11B, nhanes_data_12B) 
  
  nhanes_raw_data_C <- data.frame(reduce(nhanes_raw_data_C, full_join, by='SEQN'))
  nhanes_raw_data_A <- data.frame(reduce(nhanes_raw_data_A, full_join, by='SEQN'))
  nhanes_raw_data_B <- data.frame(reduce(nhanes_raw_data_B, full_join, by='SEQN'))
  
  nhanes_raw_data <- rbind(nhanes_raw_data_C, nhanes_raw_data_A, nhanes_raw_data_B)


# Filter and Adjust -------------------------------------------------------
  
## Filter the data to only respondents over 40 years old and < 80 yrs old
## NHANES age: 0 to 79, 80 and over
  nhanes_raw_data <- nhanes_raw_data[which(nhanes_raw_data$RIDAGEYR >= 40 
                                           & nhanes_raw_data$RIDAGEYR <= 79
                                           & nhanes_raw_data$LBXTC >= 30 
                                           & nhanes_raw_data$LBXTC <= 500 
                                           & nhanes_raw_data$LBDHDD >= 5 
                                           & nhanes_raw_data$LBDHDD <= 200 
                                           & nhanes_raw_data$BPXSY1 >= 60 
                                           & nhanes_raw_data$BPXSY1 <= 200), ]
  # check sample size now
  nrow(nhanes_raw_data)

## Calculating missingness in NHANES data
  missingness = (colMeans(is.na(nhanes_raw_data)))*100
  
## Imputing data
  nhanes_raw_data$SMQ040 <- ifelse(nhanes_raw_data$SMQ040 == 0 | is.na(nhanes_raw_data$SMQ040), 0, 1)
  nhanes_raw_data$BPQ050A <- ifelse(nhanes_raw_data$BPQ050A == 0 | is.na(nhanes_raw_data$BPQ050A), 0, 1)
  
  init = mice(nhanes_raw_data, maxit=0) 
  predM = init$predictorMatrix
  predM[,"SEQN"]=0
  nhanes_imputed_test <- mice(nhanes_raw_data, predictorMatrix = predM, m = 5, seed = 11)
  nhanes_raw_data <- complete(nhanes_imputed_test)
  # compiling vascular history rollup variable
  nhanes_raw_data$vascular_history <- as.factor(ifelse(nhanes_raw_data$MCQ160D == 1 | nhanes_raw_data$MCQ160F == 1 | nhanes_raw_data$MCQ160E == 1, 1, 0))

## Sample with replacement
  nhanes_sample = sample_n(nhanes_raw_data[!is.na(nhanes_raw_data$WTMEC2YR) & nhanes_raw_data$vascular_history == 0,], size = 100000, replace = TRUE, weight = subset(nhanes_raw_data, !is.na(WTMEC2YR) & vascular_history == 0)$WTMEC2YR)
  
## Data Cleaning
  # matching race to CHS categories
  nhanes_sample$race <- as.factor(ifelse(nhanes_sample$RIDRETH1 == 4, 2, 1))
  # compile present smoker variable to match CHS categories
  nhanes_sample$SMQ040 <- as.factor(ifelse(nhanes_sample$SMQ040 == 1 | nhanes_sample$SMQ040 == 2, 1, 0))
  # convert gender to factor variable(originally: male - 1, female - 2)
  nhanes_sample$RIAGENDR <- as.factor(nhanes_sample$RIAGENDR - 1)
  nhanes_sample$male = ifelse(nhanes_sample$RIAGENDR == 1, 0, 1)
  # convert ICD medicine I10/I10.P to be SBP treatment

#convert race back to a 1/0 variable
  nhanes_sample$race <- as.factor(ifelse(nhanes_sample$race == 1, 0, 1))

# Convert all factor columns in the data frame to numeric
  nhanes_sample[] <-  lapply(nhanes_sample, function(x) {
    if(is.factor(x)) as.numeric(as.character(x)) else x
  })

# Pooled Cohort Equation --------------------------------------------------
  
  source('pcr.R')
  
  nhanes_sample = nhanes_sample %>%
    rowwise() %>%
    mutate(
      PRS = rnorm(1,0,1),
      pce_orig = pcr(gender = ifelse(male == 1, 'M', 'F'),
                     age = RIDAGEYR,
                     race = ifelse(race == 1, 'African American', 'Other'), 
                     tot_chol = LBXTC,
                     hdl_chol = LBDHDD,
                     systolic_bp = BPXSY1,
                     bp_treatment = (BPQ050A == 1),
                     smoker = (SMQ040 == 1),
                     diabetic = (DIQ010 == 1),
                     prs_z = PRS)[[1]]/100,
           pce_prs = 1 - exp(log(1 - pce_orig) * 1.73^PRS))
  
  
# Save baseline pop
  save(nhanes_raw_data, nhanes_sample, file = "basepop.RData")
