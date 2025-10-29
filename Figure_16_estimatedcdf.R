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
################# Figure 16 #################
#############################################

# codes from plotting Figure 15 are required;



pdf("plots/Figure16.pdf", width = 8, height = 5)
  
  par(mfcol = c(2, 4))
  station <- 87
  
  options(warn = -1)
  for(predict_date in c("2025-01-30", "2025-04-30", "2025-07-30", "2025-10-30")){
    predict_date <- as.Date(predict_date)
    
    expectile_estimates <- vector(length = length(expectile_levels))
    
    for(i in 1:length(expectile_levels)){
      tau <- expectile_levels[i]
      forecast.obj <- ExpecSTnet_forecast(station = station, tau = tau, 
                                          n.ahead = 1, var.lag = 2,
                                          plot.fit = FALSE,
                                          smnet_fit = smnet_fit_tau_forecast[[i]])
      
      expectile_estimates[i] <- expectile_forecast(station = station, 
                                                   predict_date = predict_date, 
                                                   tau = tau,
                                                   forecast.obj = forecast.obj)
    }
  
    
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
    
  
    plot(cdf_est$xk, cdf_est$cdf.estimate, type = "l", 
         xlab = "x", ylab = "Cumulative distribution functions",
         main = 
           # paste0("Estimated distribution function at station ", station, 
           #             " on ", 
           predict_date)
    # )
  }
  options(warn = 0)

dev.off()