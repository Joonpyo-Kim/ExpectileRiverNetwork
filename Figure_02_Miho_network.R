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
################# Figure 2 ##################
#############################################

pdf("plots/Figure02.pdf")
  plot(x=TweedPredPoints$Long, y=TweedPredPoints$Lat, 
       xlab = "Longitude", ylab = "Latitude",
       pch=16, cex=.3)
  points((TweedPredPoints %>% dplyr::filter(StreamUnit %in% c(1, 54, 87, 90)))$Long, 
         (TweedPredPoints %>% dplyr::filter(StreamUnit %in% c(1, 54, 87, 90)))$Lat, 
         pch = 16, cex = .6, col = "red")
  points(TweedData$Long, TweedData$Lat, pch=16, cex=1, col="blue")
  miho_streams <- TweedPredPoints %>% group_by(StreamUnit) %>% 
    summarize(mean_long = mean(Longitude), mean_lat = mean(Latitude)) %>% 
    ungroup() %>% 
    dplyr::filter(StreamUnit %in% c(1, 54, 87, 90)) %>% 
    arrange(StreamUnit)
  text(miho_streams$mean_long + c(0.01, 0.01, 0, -0.01), 
       miho_streams$mean_lat + c(-0.01, 0.01, 0.01, 0),
       label = miho_streams$StreamUnit, col = "red")
dev.off()