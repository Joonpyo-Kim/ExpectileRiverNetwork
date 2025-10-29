#######################################
########## Required packages ##########
#######################################

package_list <- c("knitr", "rmarkdown", "expm",
                  "GWmodel", "matrixcalc",
                  "sampling", "lubridate", "mgcv", "spam","grDevices",
                  'ggplot2','ggmap','sf','dplyr', 
                  "fda", "vars", "quadprog", 
                  "maps", "fields")

for(pkg in package_list){
  if(pkg %in% rownames(installed.packages()) == FALSE){
    install.packages(pkg)
  }
}

sapply(package_list, require, character.only = TRUE)

##############################################
######### Loading source codes & data ########
##############################################

source("Miho_code_sources.R")

miho_new <- readRDS("Miho_forpost.RDS")

TweedData <- miho_new$TweedData
TweedPredPoints <- miho_new$TweedPredPoints

realweights <- miho_new$realweights
adjacency <- miho_new$adjacency

##############################################
######## Split data (train & test set) #######
##############################################

TweedData_train <- TweedData %>% 
  mutate(year = year(date), yday = yday(date), decimal_date = ymd_to_decimal(date)) %>% 
  dplyr::filter(year < 2023)
TweedData_test <- TweedData %>% 
  mutate(year = year(date), yday = yday(date), decimal_date = ymd_to_decimal(date)) %>% 
  dplyr::filter(year == 2023)
TweedData_test2 <- TweedData %>% 
  mutate(year = year(date), yday = yday(date), decimal_date = ymd_to_decimal(date)) %>% 
  dplyr::filter(year == 2024)

#############################################
################ CRPS scores ################
#############################################

# Codes from plotting Figure 12 & 14 are required

#############################################
################## Table 1 ##################
#############################################

CRPS_mat <- matrix(c(
  paste0(TweedData_test_crps$CRPSscore %>% mean %>% round(4), " (",
         TweedData_test_crps$CRPSscore %>% sd %>% round(4), ")"), 
  paste0(TweedData_test_crps_bench1$CRPSscore %>% mean %>% round(4), " (",
         TweedData_test_crps_bench1$CRPSscore %>% sd %>% round(4), ")"),
  paste0(TweedData_test_crps_bench2$CRPSscore %>% mean %>% round(4), " (",
         TweedData_test_crps_bench2$CRPSscore %>% sd %>% round(4), ")"),
  paste0(TweedData_test_crps_bench3$CRPSscore %>% mean %>% round(4), " (",
         TweedData_test_crps_bench3$CRPSscore %>% sd %>% round(4), ")"),
  paste0(TweedData_test_crps_2ahead$CRPSscore %>% mean %>% round(4), " (",
         TweedData_test_crps_2ahead$CRPSscore %>% sd %>% round(4), ")"), 
  paste0(TweedData_test_crps_bench1_2ahead$CRPSscore %>% mean %>% round(4), " (",
         TweedData_test_crps_bench1_2ahead$CRPSscore %>% sd %>% round(4), ")"),
  paste0(TweedData_test_crps_bench2_2ahead$CRPSscore %>% mean %>% round(4), " (",
         TweedData_test_crps_bench2_2ahead$CRPSscore %>% sd %>% round(4), ")"),
  paste0(TweedData_test_crps_bench3_2ahead$CRPSscore %>% mean %>% round(4), " (",
         TweedData_test_crps_bench3_2ahead$CRPSscore %>% sd %>% round(4), ")")
  ), ncol = 2)

colnames(CRPS_mat) <- c("1-year-ahead", "2-year-ahead")
rownames(CRPS_mat) <- c("Proposed method", "Benchmark 1", "Benchmark 2", "Benchmark 3")

sink("tables/Table1.txt")
print(CRPS_mat)
