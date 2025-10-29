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
################# Figure 13 #################
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

pdf("plots/Figure13.pdf", width = 8, height = 5)
  set.seed(250116)
  
  par(mfcol = c(2, 4))
  station <- 87
  
  ### 2023
  
  predict_date_vec <- (TweedData_test %>% 
                         dplyr::filter(location == station))$date
  
  predict_date_vec_selected <- sample(predict_date_vec, 2)
  
  for(predict_date in predict_date_vec_selected){
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
                                                   tau = tau,
                                                   forecast.obj = forecast.obj)
    }
    
    obs <- TweedData_test %>% 
      dplyr::filter(location == station & date == predict_date) %>% 
      mutate(log_nitrate = log(nitrate))
    obs <- obs$log_nitrate
    
    plot(expectile_estimates, expectile_levels, type = "o", 
         xlab = "Expectile", ylab = "Expectile Level",
         main = predict_date)
    abline(v = obs, lty = 2, col = "gray")
    
    predict_date <- as.Date(predict_date)
    cdf_at_e <- cdf_from_expect(e = expectile_estimates, tau = expectile_levels)
    
    xk <- seq(from = -0.4, to = 2.6, length = 601)
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
    obs <- obs$log_nitrate
    print(paste0("CRPS score on ", predict_date, "=", CRPS(cdf_est, obs)))
    
    plot(cdf_est$xk, cdf_est$cdf.estimate, type = "l", 
         xlab = "x", ylab = "Cumulative distribution functions",
         main = 
           # paste0("Estimated distribution function at station ", station, 
           #             " on ", 
           predict_date)
    # )
    abline(v = obs, lty = 2, col = "gray")
  }
  
  
  ### 2024
  
  predict_date_vec <- (TweedData_test2 %>% 
                         dplyr::filter(location == station))$date
  
  predict_date_vec_selected <- sample(predict_date_vec, 2)
  
  for(predict_date in predict_date_vec_selected){
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
    expectile_estimates <- sort(expectile_estimates)
    
    obs <- TweedData_test2 %>% 
      dplyr::filter(location == station & date == predict_date) %>% 
      mutate(log_nitrate = log(nitrate))
    obs <- obs$log_nitrate
    
    plot(expectile_estimates, expectile_levels, type = "o", 
         xlab = "Expectile", ylab = "Expectile Level",
         main = predict_date)
    abline(v = obs, lty = 2, col = "gray")
    
    predict_date <- as.Date(predict_date)
    cdf_at_e <- cdf_from_expect(e = expectile_estimates, tau = expectile_levels)
    
    xk <- seq(from = -0.4, to = 2.6, length = 601)
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
    obs <- obs$log_nitrate
    print(paste0("CRPS score on ", predict_date, "=", CRPS(cdf_est, obs)))
    
    plot(cdf_est$xk, cdf_est$cdf.estimate, type = "l", 
         xlab = "x", ylab = "Cumulative distribution functions",
         main = 
           # paste0("Estimated distribution function at station ", station, 
           #             " on ", 
           predict_date)
    # )
    abline(v = obs, lty = 2, col = "gray")
  }





dev.off()