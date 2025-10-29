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

# Following codes from plotting Figure 12 are required;

# expectile_levels = c(.01, .05, c(1:9)/10, .95, .99)
# smnet_fit_tau <- vector('list', length = length(expectile_levels))
# 
# for(i in 1:length(expectile_levels)){
#   tau <- expectile_levels[i]
#   smnet_fit_tau[[i]] <- ExpecSTnet(realweights, adjacency, TweedData_train, TweedPredPoints, 
#                               log.y=TRUE,
#                               plot.fig = FALSE, tau = tau, 
#                               model.type = c("c", "m", "s", "si", "t", "ti", "ts"), 
#                               penalties = NULL, 
#                               station = 87)
# }

##### Proposed method #####

df_crps <- data.frame(location = NA, date = NA, CRPSscore = NA)
station_list <- TweedData_test$location %>% unique()
station_list2 <- TweedData_test2$location %>% unique()
for(station in station_list){
  predict_date_vec <- (TweedData_test %>% 
                         dplyr::filter(location == station))$date
  
  for(predict_date in predict_date_vec){
    predict_date <- as.Date(predict_date)
    
    expectile_estimates <- vector(length = length(expectile_levels))
    
    for(i in 1:length(expectile_levels)){
      tau <- expectile_levels[i]
      forecast.obj <- ExpecSTnet_forecast(station = station, tau = tau, 
                                          n.ahead = 1,
                                          plot.fit = FALSE,
                                          smnet_fit = smnet_fit_tau[[i]])
      
      expectile_estimates[i] <- expectile_forecast(station = station, 
                                                   predict_date = predict_date, 
                                                   n.ahead = 1,
                                                   tau = tau,
                                                   forecast.obj = forecast.obj)
    }
    
    cdf_at_e <- cdf_from_expect(e = expectile_estimates, tau = expectile_levels)
    
    xk <- seq(from = -1, to = 3, by = 0.005)
    cdf_est <- data.frame(xk = xk)
    cdf.estimate = approx(x = cdf_at_e$x, y = cdf_at_e$cdf, xout = xk, 
                          yleft = 0, yright = 1)
    if(is.na(cdf.estimate$y[1])){
      cdf.estimate$y[1:(min(which(!is.na(cdf.estimate$y)))-1)] <- 0 
    }
    if(is.na(cdf.estimate$y[length(cdf.estimate$y)])){
      cdf.estimate$y[(max(which(!is.na(cdf.estimate$y)))+1):length(xk)] <- 1
    }
    cdf_est$cdf.estimate <- cdf.estimate$y

    obs <- TweedData_test %>% 
      dplyr::filter(location == station & date == predict_date) %>%
      mutate(log_nitrate = log(nitrate))
    obs <- obs$log_nitrate %>% mean()
    df_crps <- rbind(df_crps, data.frame(location = station, date = predict_date, CRPSscore = CRPS(cdf_est, obs)))
    
  }
}

df_crps <- df_crps[-1,]
df_crps$date <- as.Date(df_crps$date)

# merge(TweedData_test, df_crps) %>% arrange(location, date) %>% print(n = Inf)

TweedData_test_crps <- merge(TweedData_test, df_crps)

df_crps_2ahead <- data.frame(location = NA, date = NA, CRPSscore = NA)

for(station in station_list2){
  predict_date_vec <- (TweedData_test2 %>% 
                         dplyr::filter(location == station))$date
  
  for(predict_date in predict_date_vec){
    predict_date <- as.Date(predict_date)
    expectile_estimates <- vector(length = length(expectile_levels))
    
    for(i in 1:length(expectile_levels)){
      tau <- expectile_levels[i]
      forecast.obj <- ExpecSTnet_forecast(station = station, tau = tau, 
                                          n.ahead = 2, 
                                          plot.fit = FALSE,
                                          smnet_fit = smnet_fit_tau[[i]])
      
      expectile_estimates[i] <- expectile_forecast(station = station, 
                                                   predict_date = predict_date, 
                                                   n.ahead = 2, 
                                                   tau = tau,
                                                   forecast.obj = forecast.obj)
    }
    
    cdf_at_e <- cdf_from_expect(e = expectile_estimates, tau = expectile_levels)
    
    xk <- seq(from = -1, to = 3, by = 0.005)
    cdf_est <- data.frame(xk = xk)
    cdf.estimate = approx(x = cdf_at_e$x, y = cdf_at_e$cdf, xout = xk, 
                          yleft = 0, yright = 1)
    if(is.na(cdf.estimate$y[1])){
      cdf.estimate$y[1:(min(which(!is.na(cdf.estimate$y)))-1)] <- 0 
    }
    if(is.na(cdf.estimate$y[length(cdf.estimate$y)])){
      cdf.estimate$y[(max(which(!is.na(cdf.estimate$y)))+1):length(xk)] <- 1
    }
    cdf_est$cdf.estimate <- cdf.estimate$y

    obs <- TweedData_test2 %>% 
      dplyr::filter(location == station & date == predict_date) %>%
      mutate(log_nitrate = log(nitrate))
    obs <- obs$log_nitrate %>% mean()
    df_crps_2ahead <- rbind(df_crps_2ahead, data.frame(location = station, date = predict_date, CRPSscore = CRPS(cdf_est, obs)))
    
  }
}

df_crps_2ahead <- df_crps_2ahead[-1,]
df_crps_2ahead$date <- as.Date(df_crps_2ahead$date)

# merge(TweedData_test, df_crps_2ahead) %>% arrange(location, date) %>% print(n = Inf)

TweedData_test_crps_2ahead <- merge(TweedData_test2, df_crps_2ahead)

##### Benchmark 1 #####

df_crps_bench1 <- df_crps_bench1_2ahead <- data.frame(location = NA, date = NA, CRPSscore = NA)

for(station in station_list){
  (TweedData_train %>% dplyr::filter(location == station) %>% 
     summarize(mean = mean(log(nitrate))))$mean -> estimate_mean
  
  xk <- seq(from = -1, to = 3, by = 0.005)
  cdf_est <- data.frame(xk = xk, cdf.estimate = NA)
  
  cdf_est$cdf.estimate <- ifelse(xk <= estimate_mean, 0, 1)
  
  predict_date_vec <- (TweedData_test %>% 
                         dplyr::filter(location == station))$date
  
  for(predict_date in predict_date_vec){
    predict_date <- as.Date(predict_date)
    obs <- TweedData_test %>% 
      dplyr::filter(location == station & date == predict_date) %>% 
      mutate(log_nitrate = log(nitrate))
    obs <- obs$log_nitrate %>% mean()
    df_crps_bench1 <- rbind(df_crps_bench1, data.frame(location = station, date = predict_date, CRPSscore = CRPS(cdf_est, obs)))
    
  }
}

for(station in station_list2){
  (TweedData_train %>% dplyr::filter(location == station) %>% 
     summarize(mean = mean(log(nitrate))))$mean -> estimate_mean
  
  xk <- seq(from = -1, to = 3, by = 0.005)
  cdf_est <- data.frame(xk = xk, cdf.estimate = NA)
  
  cdf_est$cdf.estimate <- ifelse(xk <= estimate_mean, 0, 1)
  
  predict_date_vec <- (TweedData_test2 %>% 
                         dplyr::filter(location == station))$date
  
  for(predict_date in predict_date_vec){
    predict_date <- as.Date(predict_date)
    obs <- TweedData_test2 %>% 
      dplyr::filter(location == station & date == predict_date) %>% 
      mutate(log_nitrate = log(nitrate))
    obs <- obs$log_nitrate %>% mean()
    df_crps_bench1_2ahead <- rbind(df_crps_bench1_2ahead, data.frame(location = station, date = predict_date, CRPSscore = CRPS(cdf_est, obs)))
    
  }
}


df_crps_bench1 <- df_crps_bench1[-1,]
df_crps_bench1$date <- as.Date(df_crps_bench1$date)

df_crps_bench1_2ahead <- df_crps_bench1_2ahead[-1,]
df_crps_bench1_2ahead$date <- as.Date(df_crps_bench1_2ahead$date)

TweedData_test_crps_bench1 <- merge(TweedData_test, df_crps_bench1)
TweedData_test_crps_bench1_2ahead <- merge(TweedData_test2, df_crps_bench1_2ahead)

##### Benchmark 2 #####

df_crps_bench2 <- df_crps_bench2_2ahead <- data.frame(location = NA, date = NA, CRPSscore = NA)
quantile_levels <- expectile_levels

for(station in station_list){
  (TweedData_train %>% dplyr::filter(location == station) %>% 
     reframe(quant = quantile(log(nitrate), probs = quantile_levels))) -> quant_estimates 
  
  xk <- seq(from = -1, to = 3, by = 0.005)
  cdf_est <- data.frame(xk = xk, cdf.estimate = NA)
  
  cdf.estimate <- approx(x = quant_estimates$quant, y = quantile_levels, xout = xk)
  if(is.na(cdf.estimate$y[1])){
    cdf.estimate$y[1:(min(which(!is.na(cdf.estimate$y)))-1)] <- 0 
  }
  if(is.na(cdf.estimate$y[length(cdf.estimate$y)])){
    cdf.estimate$y[(max(which(!is.na(cdf.estimate$y)))+1):length(xk)] <- 1
  }
  cdf_est$cdf.estimate <- cdf.estimate$y
  
  predict_date_vec <- (TweedData_test %>% 
                         dplyr::filter(location == station))$date
  
  for(predict_date in predict_date_vec){
    predict_date <- as.Date(predict_date)
    obs <- TweedData_test %>% 
      dplyr::filter(location == station & date == predict_date) %>% 
      mutate(log_nitrate = log(nitrate))
    obs <- obs$log_nitrate %>% mean()
    df_crps_bench2 <- rbind(df_crps_bench2, data.frame(location = station, date = predict_date, CRPSscore = CRPS(cdf_est, obs)))
    
  }
  
}


for(station in station_list2){
  (TweedData_train %>% dplyr::filter(location == station) %>% 
     reframe(quant = quantile(log(nitrate), probs = quantile_levels))) -> quant_estimates 
  
  xk <- seq(from = -1, to = 3, by = 0.005)
  cdf_est <- data.frame(xk = xk, cdf.estimate = NA)
  
  cdf.estimate <- approx(x = quant_estimates$quant, y = quantile_levels, xout = xk)
  if(is.na(cdf.estimate$y[1])){
    cdf.estimate$y[1:(min(which(!is.na(cdf.estimate$y)))-1)] <- 0 
  }
  if(is.na(cdf.estimate$y[length(cdf.estimate$y)])){
    cdf.estimate$y[(max(which(!is.na(cdf.estimate$y)))+1):length(xk)] <- 1
  }
  cdf_est$cdf.estimate <- cdf.estimate$y
  
  predict_date_vec <- (TweedData_test2 %>% 
                         dplyr::filter(location == station))$date
  
  for(predict_date in predict_date_vec){
    predict_date <- as.Date(predict_date)
    obs <- TweedData_test2 %>% 
      dplyr::filter(location == station & date == predict_date) %>% 
      mutate(log_nitrate = log(nitrate))
    obs <- obs$log_nitrate %>% mean()
    df_crps_bench2_2ahead <- rbind(df_crps_bench2_2ahead, data.frame(location = station, date = predict_date, CRPSscore = CRPS(cdf_est, obs)))
    
  }
}

df_crps_bench2 <- df_crps_bench2[-1,]
df_crps_bench2$date <- as.Date(df_crps_bench2$date)
df_crps_bench2_2ahead <- df_crps_bench2_2ahead[-1,]
df_crps_bench2_2ahead$date <- as.Date(df_crps_bench2_2ahead$date)


TweedData_test_crps_bench2 <- merge(TweedData_test, df_crps_bench2)
TweedData_test_crps_bench2_2ahead <- merge(TweedData_test2, df_crps_bench2_2ahead)

##### Benchmark 3 #####

df_crps_bench3 <- df_crps_bench3_2ahead <- data.frame(location = NA, date = NA, CRPSscore = NA)

for(station in station_list){
  
  
  predict_date_vec <- (TweedData_test %>% 
                         dplyr::filter(location == station))$date
  
  for(predict_date in predict_date_vec){
    predict_date <- as.Date(predict_date)
    
    forecast.obj <- ExpecSTnet_forecast(station = station, tau = .5, 
                                        n.ahead = 1,
                                        plot.fit = FALSE, 
                                        smnet_fit = smnet_fit_tau[[which(expectile_levels == .5)]])
    
    estimate_mean <- expectile_forecast(station = station, 
                                                 predict_date = predict_date, 
                                                 n.ahead = 1,
                                                 tau = .5,
                                                 forecast.obj = forecast.obj) %>%
      as.vector()
    
    
    xk <- seq(from = -1, to = 3, by = 0.005)
    cdf_est <- data.frame(xk = xk, cdf.estimate = NA)
    
    cdf_est$cdf.estimate <- ifelse(xk <= estimate_mean, 0, 1)
    
    
    obs <- TweedData_test %>% 
      dplyr::filter(location == station & date == predict_date) %>% 
      mutate(log_nitrate = log(nitrate))
    obs <- obs$log_nitrate %>% mean()
    df_crps_bench3 <- rbind(df_crps_bench3, data.frame(location = station, date = predict_date, CRPSscore = CRPS(cdf_est, obs)))
    
  }
}

for(station in station_list2){
  predict_date_vec <- (TweedData_test2 %>% 
                         dplyr::filter(location == station))$date
  
  for(predict_date in predict_date_vec){
    predict_date <- as.Date(predict_date)
    
    forecast.obj <- ExpecSTnet_forecast(station = station, tau = .5, 
                                        n.ahead = 2,
                                        plot.fit = FALSE, 
                                        smnet_fit = smnet_fit_tau[[which(expectile_levels == .5)]])
    
    estimate_mean <- expectile_forecast(station = station, 
                                                 predict_date = predict_date, 
                                                 n.ahead = 1,
                                                 tau = .5,
                                                 forecast.obj = forecast.obj) %>%
      as.vector()
    
    xk <- seq(from = -1, to = 3, by = 0.005)
    cdf_est <- data.frame(xk = xk, cdf.estimate = NA)
    
    cdf_est$cdf.estimate <- ifelse(xk <= estimate_mean, 0, 1)
    
    obs <- TweedData_test2 %>% 
      dplyr::filter(location == station & date == predict_date) %>% 
      mutate(log_nitrate = log(nitrate))
    obs <- obs$log_nitrate %>% mean()
    df_crps_bench3_2ahead <- rbind(df_crps_bench3_2ahead, data.frame(location = station, date = predict_date, CRPSscore = CRPS(cdf_est, obs)))
    
  }
}


df_crps_bench3 <- df_crps_bench3[-1,]
df_crps_bench3$date <- as.Date(df_crps_bench3$date)
df_crps_bench3_2ahead <- df_crps_bench3_2ahead[-1,]
df_crps_bench3_2ahead$date <- as.Date(df_crps_bench3_2ahead$date)

TweedData_test_crps_bench3 <- merge(TweedData_test, df_crps_bench3)
TweedData_test_crps_bench3_2ahead <- merge(TweedData_test2, df_crps_bench3_2ahead)

#############################################
################## Table 1 ##################
#############################################

cat("CRPS score of proposed one : ", 
    TweedData_test_crps$CRPSscore %>% mean %>% round(4), "±", 
    TweedData_test_crps$CRPSscore %>% sd %>% round(4), "\n", 
    "CRPS score of benchmark 1 : ", 
    TweedData_test_crps_bench1$CRPSscore %>% mean %>% round(4), "±", 
    TweedData_test_crps_bench1$CRPSscore %>% sd %>% round(4), "\n",
    "CRPS score of benchmark 2 : ", 
    TweedData_test_crps_bench2$CRPSscore %>% mean %>% round(4), "±", 
    TweedData_test_crps_bench2$CRPSscore %>% sd %>% round(4), "\n",
    "CRPS score of benchmark 3 : ",  
    TweedData_test_crps_bench3$CRPSscore %>% mean %>% round(4), "±", 
    TweedData_test_crps_bench3$CRPSscore %>% sd %>% round(4), "\n")

cat("CRPS score of proposed one : ", 
    TweedData_test_crps_2ahead$CRPSscore %>% mean %>% round(4), "±", 
    TweedData_test_crps_2ahead$CRPSscore %>% sd %>% round(4), "\n",
    "CRPS score of benchmark 1 : ", 
    TweedData_test_crps_bench1_2ahead$CRPSscore %>% mean %>% round(4), "±", 
    TweedData_test_crps_bench1_2ahead$CRPSscore %>% sd %>% round(4), "\n",
    "CRPS score of benchmark 2 : ", 
    TweedData_test_crps_bench2_2ahead$CRPSscore %>% mean %>% round(4), "±", 
    TweedData_test_crps_bench2_2ahead$CRPSscore %>% sd %>% round(4), "\n",
    "CRPS score of benchmark 3 : ",  
    TweedData_test_crps_bench3_2ahead$CRPSscore %>% mean %>% round(4), "±", 
    TweedData_test_crps_bench3_2ahead$CRPSscore %>% sd %>% round(4)) 


pdf("plots/Figure14.pdf", width = 10 , height = 12)
  
  par(mfrow = c(2, 1))
  # loc_list <- TweedData_test_crps$location %>% unique %>% sort
  loc_list <- TweedData_test_crps_2ahead$location %>% unique %>% sort
  
  plot(0, 0, type = "n", xlab = "segment", ylab = "Mean CRPS score",
       xlim = range(1, loc_list %>% length +0), 
       ylim = c(0, 0.15), axes = FALSE, 
       main = "Average CRPS score for 1-year-ahead prediction")
  axis(2)
  axis(1, at = loc_list %>% length %>% seq, 
       labels = loc_list, 
       cex.axis = .67)
  box()
  
  TweedData_test_crps %>% dplyr::filter(location %in% loc_list) %>%
    arrange(location) %>% group_by(location) %>% 
    mutate(location = as.factor(location)) %>% 
    summarize(mean_crps = mean(CRPSscore)) %>%
    lines(type = "o", col = "black", lwd = 2, pch = 16)
  
  
  TweedData_test_crps_bench1 %>% dplyr::filter(location %in% loc_list) %>%
    arrange(location) %>% group_by(location) %>% 
    mutate(location = as.factor(location)) %>% 
    summarize(mean_crps = mean(CRPSscore)) %>% 
    lines(type = "o", col = "red", pch = 16)
  
  TweedData_test_crps_bench2 %>% dplyr::filter(location %in% loc_list) %>%
    arrange(location) %>% group_by(location) %>% 
    mutate(location = as.factor(location)) %>% 
    summarize(mean_crps = mean(CRPSscore)) %>% 
    lines(type = "o", col = "blue", pch = 16)
  
  TweedData_test_crps_bench3 %>% dplyr::filter(location %in% loc_list) %>%
    arrange(location) %>% group_by(location) %>% 
    mutate(location = as.factor(location)) %>% 
    summarize(mean_crps = mean(CRPSscore)) %>% 
    lines(type = "o", col = "green3", pch = 16)
  
  legend("topright", c("Proposed method", "Benchmark 1", "Benchmark 2", "Benchmark 3"), 
         lty = 1, lwd = c(2,1,1,1), col = c("black", "red", "blue", "green3"))
  
  plot(0, 0, type = "n", xlab = "segment", ylab = "Mean CRPS score",
       xlim = range(1, loc_list %>% length +0), 
       ylim = c(0, 0.15), axes = FALSE, 
       main = "Average CRPS score for 1-year-ahead prediction")
  axis(2)
  axis(1, at = loc_list %>% length %>% seq, 
       labels = loc_list, 
       cex.axis = .67)
  box()
  
  TweedData_test_crps_2ahead %>% dplyr::filter(location %in% loc_list) %>%
    arrange(location) %>% group_by(location) %>% 
    mutate(location = as.factor(location)) %>% 
    summarize(mean_crps = mean(CRPSscore)) %>%
    lines(type = "o", col = "black", lwd = 2, pch = 16)
  
  
  TweedData_test_crps_bench1_2ahead %>% dplyr::filter(location %in% loc_list) %>%
    arrange(location) %>% group_by(location) %>% 
    mutate(location = as.factor(location)) %>% 
    summarize(mean_crps = mean(CRPSscore)) %>% 
    lines(type = "o", col = "red", pch = 16)
  
  TweedData_test_crps_bench2_2ahead %>% dplyr::filter(location %in% loc_list) %>%
    arrange(location) %>% group_by(location) %>% 
    mutate(location = as.factor(location)) %>% 
    summarize(mean_crps = mean(CRPSscore)) %>% 
    lines(type = "o", col = "blue", pch = 16)
  
  TweedData_test_crps_bench3_2ahead %>% dplyr::filter(location %in% loc_list) %>%
    arrange(location) %>% group_by(location) %>% 
    mutate(location = as.factor(location)) %>% 
    summarize(mean_crps = mean(CRPSscore)) %>% 
    lines(type = "o", col = "green3", pch = 16)
  
  legend("topright", c("Proposed method", "Benchmark 1", "Benchmark 2", "Benchmark 3"), 
         lty = 1, lwd = c(2,1,1,1), col = c("black", "red", "blue", "green3"))
dev.off()
