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
################# Figure 11 #################
#############################################

smnet_fit_50 <- ExpecSTnet(realweights, adjacency, TweedData_train, TweedPredPoints, 
                           log.y=TRUE,
                           plot.fig = FALSE, tau = 0.5, 
                           model.type = c("c", "m", "s", "si", "t", "ti", "ts"), 
                           penalties = NULL, 
                           station = 87)

smnet_fit_90 <- ExpecSTnet(realweights, adjacency, TweedData_train, TweedPredPoints, 
                           log.y=TRUE,
                           plot.fig = FALSE, tau = 0.9, 
                           model.type = c("c", "m", "s", "si", "t", "ti", "ts"), 
                           penalties = NULL, 
                           station = 87)

pdf("plots/Figure11.pdf", width = 10, height = 5)
  par(mfrow = c(1, 2))
  forecast.obj <- ExpecSTnet_forecast(station = 87, tau = .5, 
                                      plot.fit = TRUE, 
                                      smnet_fit = smnet_fit_50)
  forecast.obj <- ExpecSTnet_forecast(station = 87, tau = .9, 
                                      plot.fit = TRUE, 
                                      smnet_fit = smnet_fit_90)

dev.off()