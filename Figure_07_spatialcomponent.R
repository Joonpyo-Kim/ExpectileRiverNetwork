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
################# Figure 7 ##################
#############################################

# register_stadiamaps(key = "YOUR-API-KEY")

######## Spatial components (average) #######
poo <- ExpecSTnet(realweights, adjacency, TweedData_train, TweedPredPoints, 
                  log.y=TRUE,
                  plot.fig = FALSE, tau = 0.5, 
                  model.type = c("c", "m", "s", "si", "t", "ti", "ts"), 
                  penalties = NULL, 
                  station = 54)

model.type = c("c", "m", "s", "si", "t", "ti", "ts")
comps <- vector("list", length(model.type))
knts<-c(10,10) #s모형의 knot 개수, t모형의 knot 개수
n.segments<-nrow(adjacency) #강 물줄기 갯수(113개)
beta_1 <- poo$beta_hat[1]
dims<-get.dimension(model.type = model.type, n.segments = n.segments, knts = knts) #1 113 10 1130 10 1130 100 (c, m, s, si, t, ti, ts)
locs<-cumsum(dims)
fit.point <- poo$beta_hat[(locs[match("m", model.type)-1]+1):locs[match("m", model.type)]]
n.cols<-112 
ngrid<-n.cols
brks <- seq(-1.5,1.5,length = ngrid + 1) ## Should be properly modified! 
col.nums<-cut(fit.point[as.numeric(TweedPredPoints$StreamUnit)], breaks = brks, labels = FALSE)
palette<-colorRampPalette(c("cyan", "green", "yellow", "red", "black")) #원래 c("cyan", "green", "yellow", "red", "black")
main.cols<-palette(n.cols)

spatial_component <- TweedPredPoints
spatial_component$color <- main.cols[col.nums]
data_merged_miho_new <- data_merged_miho
data_merged_miho_new$color <- spatial_component$color

city_long <- c(127.38, 127.3, 127.57, 127.18, 127.52, 127.6, 127.43, 127.36, 127.23)
city_lat <- c(36.9, 37.0, 37.0, 36.78, 36.7, 36.63, 36.57, 36.65, 36.61)
city_name <- c("Jincheon", "Anseong", "Eumseong", "Cheonan", 
               "Cheongwon", "Sangdang", "Seowon", "Heungdeok", 
               "Sejong")

bbox2 <- c(127.1, 36.5, 127.7, 37.05)
map2 <-get_stadiamap(bbox2, maptype = "stamen_terrain_background", 
                     color = "bw", force = TRUE) 
p1 <- ggmap(map2) + 
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
           fill = "white", alpha = 0.6)  + 
  geom_sf(data=sigungu_Geum, inherit.aes = FALSE, fill=NA, color="black") + 
  geom_path(data=data_merged_miho_new, size=1, aes(x=long,y=lat,group=group),
            # color = "skyblue") + 
            color=data_merged_miho_new$color) +
  geom_point(data = industrial_data_transform2, size = 2,
             aes(x = long, y = lat, 
                 shape = danji_index, color = danji_index), 
             show.legend = FALSE) + 
  labs(title="Mean spatial component", x = "Longitude", y = "Latitude", 
       color="Type", shape="Type") + 
  annotate("text", x = city_long, y = city_lat, label = city_name, size = 4)
# p1

pdf("plots/Figure07-1.pdf", width = 4.5, height = 4.5)
  p1
dev.off()

######## Spatial components (90% expectile) #######
poo9 <- ExpecSTnet(realweights, adjacency, TweedData_train, TweedPredPoints, 
                   log.y=TRUE,
                   plot.fig = FALSE, tau = 0.9, 
                   model.type = c("c", "m", "s", "si", "t", "ti", "ts"), 
                   penalties = NULL, 
                   station = 54)

model.type = c("c", "m", "s", "si", "t", "ti", "ts")
comps <- vector("list", length(model.type))
knts<-c(10,10) #s모형의 knot 개수, t모형의 knot 개수
n.segments<-nrow(adjacency) #강 물줄기 갯수(113개)
beta_1 <- poo9$beta_hat[1]
dims<-get.dimension(model.type = model.type, n.segments = n.segments, knts = knts) #1 113 10 1130 10 1130 100 (c, m, s, si, t, ti, ts)
locs<-cumsum(dims)
fit.point <- poo9$beta_hat[(locs[match("m", model.type)-1]+1):locs[match("m", model.type)]]
n.cols<-112 
ngrid<-n.cols
brks <- seq(-1.5,1.5,length = ngrid + 1) ## Should be properly modified! 
col.nums<-cut(fit.point[as.numeric(TweedPredPoints$StreamUnit)], breaks = brks, labels = FALSE)
palette<-colorRampPalette(c("cyan", "green", "yellow", "red", "black")) #원래 c("cyan", "green", "yellow", "red", "black")
main.cols<-palette(n.cols)

spatial_component9 <- TweedPredPoints
spatial_component9$color <- main.cols[col.nums]
data_merged_miho_new9 <- data_merged_miho
data_merged_miho_new9$color <- spatial_component9$color

p2 <- ggmap(map2) + 
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
           fill = "white", alpha = 0.6)  + 
  geom_sf(data=sigungu_Geum, inherit.aes = FALSE, fill=NA, color="black") + 
  geom_path(data=data_merged_miho_new9, size=1, aes(x=long,y=lat,group=group),
            # color = "skyblue") + 
            color=data_merged_miho_new9$color) +
  geom_point(data = industrial_data_transform2, size = 2,
             aes(x = long, y = lat, 
                 shape = danji_index, color = danji_index), 
             show.legend = FALSE) + 
  labs(title="90% expectile spatial component", x = "Longitude", y = "Latitude", 
       color="Type", shape="Type") + 
  annotate("text", x = city_long, y = city_lat, label = city_name, size = 4)

pdf("plots/Figure07-2.pdf", width = 4.5, height = 4.5)
  p2
dev.off()
