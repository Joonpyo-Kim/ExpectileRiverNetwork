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

#############################################
################# Figure 3 ##################
#############################################

yday_my <- function(date){
  if(length(date) > 1){
    sapply(date, yday_my)
  }else{
    if(year(date) %% 4 == 0 & month(date) >= 3){
      return(yday(date)-1)
    }else{
      return(yday(date))
    }
  }
  
}

TweedData %>% dplyr::filter(location == 87) %>% 
  dplyr::mutate(year = year(date), yday = yday_my(date)) %>% 
  dplyr::select(year, yday, date) -> dates87

TweedData %>% dplyr::filter(location == 1) %>% 
  dplyr::mutate(year = year(date), yday = yday_my(date)) %>% 
  dplyr::select(year, yday, date) -> dates1

pdf("plots/Figure03.pdf", width = 7, height = 5)
  par(mar = c(.1, .1, 2.1, 1.1))
  plot(0, 0, xlim = c(-15, 400), ylim = c(2012, 2022), xlab = "", ylab = "", axes = FALSE, 
       # main = "Recorded date at sation 1 and 87")
       main = "")
  legend("topright", c("1", "87"), pch = c(17, 16), col = c( "red", "blue"))
  
  abline(v = 365, col = "gray80")
  for(mth in 1:12){
    mth_yday <- yday(paste0("2021-", mth, "-01"))
    abline(v = mth_yday, col = "gray80")
    text(mth_yday + 15, 2012, labels = month.name[mth], cex = .5)
  }
  
  for(yr in 2013:2022){
    text(-15, yr, labels = yr)
    lines(c(1, 365), rep(yr, 2))
    (dates87 %>% dplyr::filter(year == yr))$yday -> dates87_yday
    (dates1  %>% dplyr::filter(year == yr))$yday ->  dates1_yday
    points(dates87_yday, rep(yr, length(dates87_yday)), col = "blue", pch = 16)
    points( dates1_yday, rep(yr, length (dates1_yday)), col = "red",  pch = 17)
  }
dev.off()