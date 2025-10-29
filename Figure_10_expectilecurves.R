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
################# Figure 10 #################
#############################################
ExpecSTnet_plot <- function(smnet_fit, 
                            plot_fig = c("spatial", "seasonal", "trend", "estimated"), 
                            station = 87){
  # smnet_fit : result of ExpecSTnet function
  # plot_fig : plot "spatial" : spatial / "seasonal" : seasonal / "trend" : temporal component /
  #                 "estimated" : estimated curves
  # station : only valid for plot_fig = "estimated"; which segment to be plotted 
  tau <- smnet_fit$tau
  model.type <- smnet_fit$model.type 
  parameters <- smnet_fit$beta_hat
  resids <- smnet_fit$resids
  
  dates <- as.character(smnet_fit$TweedData$date)
  dates2 <- ymd(dates); 
  ymd_to_decimal <- function(date){
    if(year(date) %% 4 != 0){
      return(yday(date)/365)
    }else if(year(date) %% 400 != 0){
      return(yday(date)/366)
    }else{
      return(yday(date)/365)
    }
  }
  decimal.day <- sapply(dates2, ymd_to_decimal) 
  
  n.segments<-nrow(smnet_fit$adjacency) 
  p.dims<-n.segments
  knts<-c(10,10) 
  dims<-get.dimension(model.type = model.type, n.segments = n.segments, knts = knts) #1 113 10 1130 10 1130 100 (c, m, s, si, t, ti, ts)
  c.dims<-cumsum(dims)
  n.par<-sum(dims) 
  
  locs<-cumsum(dims) 
  if(length(model.type) == 1){locs<-c(1, locs)}
  fit.point <- parameters[(locs[match("m", model.type)-1]+1):locs[match("m", model.type)]]
  beta_1 <- parameters[1]
  
  comps <- vector("list", length(model.type))
  min.yr <- min(year(dates))
  max.yr <- max(year(dates))
  yrs <- seq(min.yr, max.yr)
  n.yrs <- length(yrs)
  day.seq <- seq(min.yr + 5/365, max.yr + 1, length.out = (365/5)*n.yrs) 
  
  if(plot_fig == "spatial"){
    comps<-c(comps, list(spatial = parameters[(locs[match("m", model.type)-1]+1):locs[match("m", model.type)]]))
    fit.point <-comps$spatial
    n.cols<-112 
    ngrid<-n.cols
    brks <- seq(-1.5,1.5,length = ngrid + 1) ## Should be properly modified! 
    col.nums<-cut(fit.point[as.numeric(TweedPredPoints$StreamUnit)], breaks = brks, labels = FALSE)
    palette<-colorRampPalette(c("cyan", "green", "yellow", "red", "black")) 
    main.cols<-palette(n.cols)
    
    # plot
    par(mar = c(1, 0, 2, 0))
    plot(TweedPredPoints$Longitude,TweedPredPoints$Latitude, 
         pch = 20, col = main.cols[col.nums],
         cex = TweedPredPoints$Weights, bty = "n", xlab = "",
         ylab = "", xaxt = "n", yaxt = "n", ylim = range(TweedPredPoints$Latitude))
    title(paste0("Spatial Component for ", tau*100, "% expectile"))
  }else if(plot_fig == "seasonal"){
    s.par<-parameters[(locs[match("s", model.type) - 1] + 1):locs[match("s", model.type)]]
    seas<-(1:365)/365
    seas.bas<-cSplineDes(1:365/365, knots = 0:knts[1]/knts[1], ord=3)
    seas.comp<-seas.bas %*% s.par
    which.day<-function(n){seas.comp[n]}
    partial.residuals<-resids+apply(as.matrix(yday(dates2)), 1, which.day) 
    day.resids<-yday(dates2)
    
    # plot
    par(mar=c(4,5,2,2))
    plot(day.resids,partial.residuals,  pch = ".", xlab = "",
         ylim = c(-2, 2), ylab = "log Nitrate concentrations (mg/l)", 
         type = "n", main = "", cex.axis = 1, cex.lab = 1, 
         axes = FALSE)
    title(xlab = "Day in year", line = 2)
    axis(2)
    month_length <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    axis(1, tick = FALSE, 
         at = cumsum(month_length) - month_length/2,
         labels = c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"), 
         cex.axis = .8, 
         line = -0.5)
    axis(1, labels = FALSE, tick = TRUE, 
         at = c(0, cumsum(month_length)))
    box() 
    # polygon(c(1:365, 365:1), c( 2*se.s+seas.comp, rev(-2*se.s+seas.comp)),
    #         col = "grey", border = NA)
    points(day.resids,partial.residuals,  col = "brown",pch = ".", cex = 2)
    lines(1:365, seas.comp, xlab = "Day in year", lwd = 3)
    title(paste0("Seasonal component for ", tau, " expectile of response"))
  }else if(plot_fig == "trend"){
    years<-year(dates)
    t.par<-parameters[(locs[match("t", model.type) - 1] + 1):locs[match("t", model.type)]]
    tren<-seq(min(years), max(years))
    tren.bas<-bbase(tren, xl = min(years), xr = max(years), nseg = (knts[2] - 3))
    tren.comp<-tren.bas %*% t.par 
    
    which.yr<-function(n){tren.comp[n]}
    yr.index<-year(dates) %% min(year(dates)) + 1
    partial.residuals<-resids+apply(as.matrix(yr.index), 1, which.yr) 
    
    # plot
    par(mar=c(4,5,2,2))
    plot(year(dates), partial.residuals, pch = ".", xlab = "", 
         ylab = "log Nitrate concentrations (mg/l)", ylim = c(-2, 2), 
         type = "n", main = "", cex.axis = 1, cex.lab = 1)
    title(xlab = "Year of observation", line = 2.5)
    # polygon(c(min(years):max(years),max(years):min(years)),
    #         c(2*se.t + tren.comp, rev( -2*se.t + tren.comp)), 
    #         col = "grey", border = NA)
    points(year(dates), partial.residuals, col = "brown",pch = ".", cex = 0.2)
    lines(tren, tren.comp, lwd = 2)
    title(paste0("Trend component for ", tau, " expectile of response"))
  }else if(plot_fig == "estimated"){
    poo <- ExpecSTnet_plot.estimated(station = station, smnet_obj = smnet_fit, 
                                     plot.fig = TRUE)
  }else{
    stop("Error: Invalid plot_fig")
  }
}


smnet_fit_50 <- ExpecSTnet(realweights, adjacency, TweedData_train, TweedPredPoints, 
                           log.y=TRUE,
                           plot.fig = FALSE, tau = 0.5, 
                           model.type = c("c", "m", "s", "si", "t", "ti", "ts"), 
                           penalties = NULL, 
                           station = 87)

smnet_fit_01 <- ExpecSTnet(realweights, adjacency, TweedData_train, TweedPredPoints, 
                           log.y=TRUE,
                           plot.fig = FALSE, tau = 0.01, 
                           model.type = c("c", "m", "s", "si", "t", "ti", "ts"), 
                           penalties = NULL, 
                           station = 87)

smnet_fit_10 <- ExpecSTnet(realweights, adjacency, TweedData_train, TweedPredPoints, 
                           log.y=TRUE,
                           plot.fig = FALSE, tau = 0.1, 
                           model.type = c("c", "m", "s", "si", "t", "ti", "ts"), 
                           penalties = NULL, 
                           station = 87)

smnet_fit_30 <- ExpecSTnet(realweights, adjacency, TweedData_train, TweedPredPoints, 
                           log.y=TRUE,
                           plot.fig = FALSE, tau = 0.3, 
                           model.type = c("c", "m", "s", "si", "t", "ti", "ts"), 
                           penalties = NULL, 
                           station = 87)

smnet_fit_70 <- ExpecSTnet(realweights, adjacency, TweedData_train, TweedPredPoints, 
                           log.y=TRUE,
                           plot.fig = FALSE, tau = 0.7, 
                           model.type = c("c", "m", "s", "si", "t", "ti", "ts"), 
                           penalties = NULL, 
                           station = 87)

smnet_fit_90 <- ExpecSTnet(realweights, adjacency, TweedData_train, TweedPredPoints, 
                           log.y=TRUE,
                           plot.fig = FALSE, tau = 0.9, 
                           model.type = c("c", "m", "s", "si", "t", "ti", "ts"), 
                           penalties = NULL, 
                           station = 87)

smnet_fit_99 <- ExpecSTnet(realweights, adjacency, TweedData_train, TweedPredPoints, 
                           log.y=TRUE,
                           plot.fig = FALSE, tau = 0.99, 
                           model.type = c("c", "m", "s", "si", "t", "ti", "ts"), 
                           penalties = NULL, 
                           station = 87)

pdf("plots/Figure10.pdf", width = 7, height = 5)
  par(mar = c(5.1, 4.1, 4.1, 8.1), xpd = TRUE)
  smnet_plot_50 <- ExpecSTnet_plot.estimated(smnet_fit_50, plot.fig = FALSE, station = 87)
  obs.ind <-which(smnet_fit_50$TweedData$location == 87)
  dates = as.character(smnet_fit_50$TweedData$date)
  dates.full <-dates[obs.ind]
  dates.stat <-decimal_date(as.Date(dates.full))
  vals.stat <- smnet_fit_50$response[obs.ind]
  min.yr <- min(year(dates))
  max.yr <- max(year(dates))
  yrs <- seq(min.yr, max.yr)
  n.yrs <- length(yrs)
  day.seq <- seq(min.yr + 5/365, max.yr + 1, length.out = (365/5)*n.yrs) 
  
  plot(dates.stat, vals.stat, type = "n", xlab = "Time", ylab = "log(Nitrate concentrations)")
  points(dates.stat, vals.stat)
  
  
  smnet_plot_01 <- ExpecSTnet_plot.estimated(smnet_fit_01, plot.fig = FALSE, station = 87)
  lines(day.seq, c(t(smnet_plot_01$main))+ c(t(smnet_plot_01$inter)), type = "l",  
        lwd=2, lty=1, col = "blue")
  
  smnet_plot_10 <- ExpecSTnet_plot.estimated(smnet_fit_10, plot.fig = FALSE, station = 87)
  lines(day.seq, c(t(smnet_plot_10$main))+ c(t(smnet_plot_10$inter)), type = "l",  
        lwd=2, lty=1, col = "cyan")
  
  smnet_plot_30 <- ExpecSTnet_plot.estimated(smnet_fit_30, plot.fig = FALSE, station = 87)
  lines(day.seq, c(t(smnet_plot_30$main))+ c(t(smnet_plot_30$inter)), type = "l",  
        lwd=2, lty=1, col = "cyan")
  
  smnet_plot_70 <- ExpecSTnet_plot.estimated(smnet_fit_70, plot.fig = FALSE, station = 87)
  lines(day.seq, c(t(smnet_plot_70$main))+ c(t(smnet_plot_70$inter)), type = "l",  
        lwd=2, lty=1, col = "orange")
  
  smnet_plot_90 <- ExpecSTnet_plot.estimated(smnet_fit_90, plot.fig = FALSE, station = 87)
  lines(day.seq, c(t(smnet_plot_90$main))+ c(t(smnet_plot_90$inter)), type = "l",  
        lwd=2, lty=1, col = "orange")
  
  smnet_plot_99 <- ExpecSTnet_plot.estimated(smnet_fit_99, plot.fig = FALSE, station = 87)
  lines(day.seq, c(t(smnet_plot_99$main))+ c(t(smnet_plot_99$inter)), type = "l",  
        lwd=2, lty=1, col = "red")
  
  lines(day.seq, c(t(smnet_plot_50$main))+ c(t(smnet_plot_50$inter)), 
        type = "l",  lwd=2, lty=1) #col = "limegreen",
  
  # title("Segment 87")
  
  legend("bottomright", inset=c(-0.25,0), 
         title = "tau", 
         legend=c("0.99", "0.9", "0.7", "0.5", "0.3", "0.1", "0.01"), 
         lty = 1, 
         col = c("red", "orange", "orange", "black", "cyan", "cyan", "blue"))

dev.off()