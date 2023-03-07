library(tidyverse)
library(readxl)

# get the US life tables from CDC(2020)
# https://www.cdc.gov/nchs/data/nvsr/nvsr71/nvsr71-01.pdf

flt2020_raw = read_xlsx("life table/flt2020.xlsx", skip = 2, col_names = TRUE)
colnames(flt2020_raw)[1] = "age"
flt2020 = flt2020_raw %>%
  drop_na() %>%
  rowwise %>%
  mutate(Age = ifelse(str_detect(age, '–'), 
                      as.numeric(str_sub(age,1,unlist(gregexpr('–', age))[1]-1)), 
                      as.numeric(str_sub(age,1,unlist(gregexpr(' ', age))[1]-1))),
         rate = -log(1-qx)/1)

mlt2020_raw = read_xlsx("life table/flt2020.xlsx", skip = 2, col_names = TRUE)
colnames(mlt2020_raw)[1] = "age"
mlt2020 = mlt2020_raw %>%
  drop_na() %>%
  rowwise %>%
  mutate(Age = ifelse(str_detect(age, '–'), 
                      as.numeric(str_sub(age,1,unlist(gregexpr('–', age))[1]-1)), 
                      as.numeric(str_sub(age,1,unlist(gregexpr(' ', age))[1]-1))),
         rate = -log(1-qx)/1)

# get the CVD-cause death numbers from GBD(2019)
# https://vizhub.healthdata.org/gbd-results/

CVD_death_raw = read.csv("life table/GBD_2019_CVD_death_percent.csv")
CVD_death = CVD_death_raw %>%
  select(sex,age,val,upper,lower)
CVD_death_f = CVD_death %>% filter(sex == "Female") %>% select(-sex)
CVD_death_m = CVD_death %>% filter(sex == "Male") %>% select(-sex)

# change the age group into a constant number and replicate the left parts

CVD_death_f1 = CVD_death_f %>% filter(str_detect(age, "<")) %>%
  mutate(start = 0, 
         reptimes = as.numeric(str_sub(age,
                            unlist(gregexpr('<', age))[1]+1,
                            unlist(gregexpr(' ', age))[1]-1)),
         end = start + reptimes - 1)

CVD_death_f2 = CVD_death_f %>% 
  filter(str_detect(age,"-") & !str_detect(age, " ")) %>%
  rowwise %>%
  mutate(start = as.numeric(str_sub(age,1,unlist(gregexpr('-', age))[1]-1)),
         end = as.numeric(str_sub(age,unlist(gregexpr('-', age))[1]+1,str_length(age))),
         reptimes = end - start + 1)

CVD_death_f3 = CVD_death_f %>% filter(str_detect(age,"-") & str_detect(age, " ")) %>%
  rowwise %>%
  mutate(start = as.numeric(str_sub(age,1,unlist(gregexpr('-', age))[1]-1)),
         end = as.numeric(str_sub(age,unlist(gregexpr('-', age))[1]+1,unlist(gregexpr(' ', age))[1]-1)),
         reptimes = end - start + 1)

CVD_death_female = rbind(
  as.data.frame(lapply(CVD_death_f1, rep, CVD_death_f1$reptimes)) %>%
  group_by(age) %>%
  mutate(Age = start - 1 + 1:n()) %>%
  ungroup(),
  as.data.frame(lapply(CVD_death_f2, rep, CVD_death_f2$reptimes)) %>%
    group_by(age) %>%
    mutate(Age = start - 1 + 1:n()) %>%
    ungroup(),
  as.data.frame(lapply(CVD_death_f3, rep, CVD_death_f3$reptimes)) %>%
    group_by(age) %>%
    mutate(Age = start - 1 + 1:n()) %>%
    ungroup()) %>%
  arrange(Age) %>%
  select(Age, val, upper, lower)

minrow = min(nrow(flt2020), nrow(CVD_death_female))

## dx
flt2020_nonCVD = data.frame(Age = 0:(minrow-1),
                            dx = flt2020[1:minrow, "dx"] * (1-CVD_death_female[1:minrow, "val"]))
flt2020_CVD = data.frame(Age = 0:(minrow-1),
                            dx = flt2020[1:minrow, "dx"] * CVD_death_female[1:minrow, "val"])

## male

CVD_death_m1 = CVD_death_m %>% filter(str_detect(age, "<")) %>%
  mutate(start = 0, 
         reptimes = as.numeric(str_sub(age,
                                       unlist(gregexpr('<', age))[1]+1,
                                       unlist(gregexpr(' ', age))[1]-1)),
         end = start + reptimes - 1)

CVD_death_m2 = CVD_death_m %>% 
  filter(str_detect(age,"-") & !str_detect(age, " ")) %>%
  rowwise %>%
  mutate(start = as.numeric(str_sub(age,1,unlist(gregexpr('-', age))[1]-1)),
         end = as.numeric(str_sub(age,unlist(gregexpr('-', age))[1]+1,str_length(age))),
         reptimes = end - start + 1)

CVD_death_m3 = CVD_death_m %>% filter(str_detect(age,"-") & str_detect(age, " ")) %>%
  rowwise %>%
  mutate(start = as.numeric(str_sub(age,1,unlist(gregexpr('-', age))[1]-1)),
         end = as.numeric(str_sub(age,unlist(gregexpr('-', age))[1]+1,unlist(gregexpr(' ', age))[1]-1)),
         reptimes = end - start + 1)

CVD_death_male = rbind(
  as.data.frame(lapply(CVD_death_m1, rep, CVD_death_m1$reptimes)) %>%
    group_by(age) %>%
    mutate(Age = start - 1 + 1:n()) %>%
    ungroup(),
  as.data.frame(lapply(CVD_death_m2, rep, CVD_death_m2$reptimes)) %>%
    group_by(age) %>%
    mutate(Age = start - 1 + 1:n()) %>%
    ungroup(),
  as.data.frame(lapply(CVD_death_m3, rep, CVD_death_m3$reptimes)) %>%
    group_by(age) %>%
    mutate(Age = start - 1 + 1:n()) %>%
    ungroup()) %>%
  arrange(Age) %>%
  select(Age, val, upper, lower)

minrow = min(nrow(mlt2020), nrow(CVD_death_male))

## dx
mlt2020_nonCVD = data.frame(Age = 0:(minrow-1),
                            dx = mlt2020[1:minrow, "dx"] * (1-CVD_death_male[1:minrow, "val"]))
mlt2020_CVD = data.frame(Age = 0:(minrow-1),
                         dx = mlt2020[1:minrow, "dx"] * CVD_death_male[1:minrow, "val"])

# validation
ggplot(flt2020_nonCVD, aes(x = Age, y = dx)) +
  geom_point(color = "blue")+
  geom_point(data = flt2020_CVD, aes(x = Age, y = dx), color = "red") +
  ggtitle("Female Number dying between ages x and x + 1 (blue: non-CVD cause, red: CVD cause)")

ggplot(flt2020_nonCVD, aes(x = Age, y = cumsum(dx))) +
  geom_point(color = "blue")+
  geom_point(data = flt2020_CVD, aes(x = Age, y = cumsum(dx)), color = "red") +
  ggtitle("Female Number dying between ages 0 and x + 1 (blue: non-CVD cause, red: CVD cause)")


