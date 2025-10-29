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
################# Figure 15 #################
#############################################

expectile_levels = c(.01, .05, c(1:9)/10, .95, .99)
smnet_fit_tau_forecast <- vector('list', length = length(expectile_levels))

for(i in 1:length(expectile_levels)){
  tau <- expectile_levels[i]
  smnet_fit_tau_forecast[[i]] <- ExpecSTnet(realweights, adjacency, TweedData, TweedPredPoints, 
                                            log.y=TRUE,
                                            plot.fig = FALSE, tau = tau, 
                                            model.type = c("c", "m", "s", "si", "t", "ti", "ts"), 
                                            penalties = NULL, 
                                            station = 87)
}

pdf("plots/Figure15.pdf", width = 12, height = 6)
  par(mfcol = c(2,4))
  
  col.vec <- c(rep("blue", 3), 
               rep("cyan3", 3), 
               "black", 
               rep("orange", 3), 
               rep("red", 3))
  
  options(warn = -1)
  for(station in c(16, 44, 87, 107)){
    
    plot(0, 0, type = "n", xlim = c(0, 1), 
         ylim = c(-0.5, 2.2),
         xlab = "", ylab = "", 
         main = paste0("Segment ", station), 
         axes = FALSE)
    axis(2)
    axis(1, tick = FALSE, 
         at = c(1:12)/12 - .5/12,
         labels = c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"), 
         cex.axis = .8, 
         line = -1)
    axis(1, labels = FALSE, tick = TRUE, 
         at = c(0:12)/12)
    box()
    m <- 1
    for(i in 1:length(expectile_levels)){
      tau <- expectile_levels[i]
      forecast.obj <- ExpecSTnet_forecast(station = station, tau = tau, 
                                          n.ahead = 1, var.lag = 2,
                                          plot.fit = FALSE,
                                          smnet_fit = smnet_fit_tau_forecast[[i]])
      poo.pca <- forecast.obj$pca.fd
      var2 <- forecast.obj$var.forecast
      forecast_coef <- poo.pca$meanfd$coefs + poo.pca$harmonics$coefs %*% (unlist(var2$fcst) %>% matrix(ncol = 3))[1,] 
      plot(fd(coef = forecast_coef, basisobj = forecast.obj$data.fd$basis), 
           add = TRUE, col = col.vec[m], href = FALSE)
      m <- m + 1
    }
    
    ##### Probabilistic forecast based CI #####
    predict_dates_ci <- seq(from = as.Date("2025-01-05"), to = as.Date("2025-12-31"), by = 10)
    ci_probfore <- data.frame(date = predict_dates_ci, lower = NA, upper = NA)
    
    for(predict_date in predict_dates_ci){
      predict_date <- as.Date(predict_date)
      expectile_estimates <- vector(length = length(expectile_levels))
      
      for(i in 1:length(expectile_levels)){
        expectile_estimates[i] <- expectile_forecast(station = station, 
                                                   predict_date = predict_date, 
                                                   tau = tau, var.lag = 2,
                                                   n.ahead = 1,
                                                   smnet_fit = smnet_fit_tau_forecast[[i]])
      }
      cdf_at_e <- cdf_from_expect(e = expectile_estimates, tau = expectile_levels)
      xk = seq(from = -1, to = 3, by = 0.005)
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
      
      xk_idx_lower <- (cdf_est %>% dplyr::filter(cdf.estimate < 0.025) %>% dim)[1]
      ci_probfore$lower[which(ci_probfore$date == predict_date)] <-
        cdf_est$xk[xk_idx_lower] + 
        diff(cdf_est$xk)[xk_idx_lower] / diff(cdf_est$cdf.estimate)[xk_idx_lower] * 
        (0.025 - cdf_est$cdf.estimate[xk_idx_lower])
      
      xk_idx_upper <- (cdf_est %>% dplyr::filter(cdf.estimate <= 0.975) %>% dim)[1]
      ci_probfore$upper[which(ci_probfore$date == predict_date)] <-
        cdf_est$xk[xk_idx_upper] + 
        diff(cdf_est$xk)[xk_idx_upper] / diff(cdf_est$cdf.estimate)[xk_idx_upper] * 
        (0.975 - cdf_est$cdf.estimate[xk_idx_upper])
      
    }
    
    plot(0, 0, type = "n", xlim = c(0, 1), 
         # ylim = c(min(ci_probfore$lower), max(ci_probfore$upper)),
         ylim = c(-0.5, 2.2),
         xlab = "", ylab = "", 
         main = paste0("Segment ", station), 
         axes = FALSE)
    axis(2)
    axis(1, tick = FALSE, 
         at = c(1:12)/12 - .5/12,
         labels = c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"), 
         cex.axis = .8, 
         line = -1)
    axis(1, labels = FALSE, tick = TRUE, 
         at = c(0:12)/12)
    box()
    
    forecast.obj <- ExpecSTnet_forecast(station = station, tau = tau, 
                                        n.ahead = 1, var.lag = 2,
                                        plot.fit = FALSE,
                                        smnet_fit = smnet_fit_tau_forecast[[which(expectile_levels == .5)]])
    poo.pca <- forecast.obj$pca.fd
    var2 <- forecast.obj$var.forecast
    forecast_coef <- poo.pca$meanfd$coefs + poo.pca$harmonics$coefs %*% (unlist(var2$fcst) %>% matrix(ncol = 3))[1,] 
    plot(fd(coef = forecast_coef, basisobj = forecast.obj$data.fd$basis), 
         add = TRUE, col = "blue")
  
    
    lines(yday(predict_dates_ci)/365, ci_probfore$lower, col = "gray50")
    lines(yday(predict_dates_ci)/365, ci_probfore$upper, col = "gray50")
  }
  options(warn = 0)

dev.off()