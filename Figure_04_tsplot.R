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
################# Figure 4 ##################
#############################################

pdf("plots/Figure04.pdf", width = 10, height = 8)
  par(mfrow = c(2, 2))
  for(station in c(1, 54, 87, 90)){
    
    
    TweedData %>% filter(location == station) %>% 
      mutate(lognitrate = log(nitrate)) %>% 
      dplyr::select(c("date", "lognitrate")) %>% 
      arrange(date) %>%
      plot(type="l", xlab = "Time", ylab = "log Nitrate Level",
           main = paste("Time Series Plot at Segment", station))
    
  }
dev.off()