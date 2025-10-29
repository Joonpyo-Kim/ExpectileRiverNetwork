#######################################
########## Required packages ##########
#######################################

package_list <- c("knitr", "rmarkdown", "expm",
                  "GWmodel", "matrixcalc",
                  "sampling", "lubridate", "mgcv", "spam","grDevices",
                  'ggplot2','ggmap','sf','dplyr', 
                  "fda", "vars", "quadprog", 
                  "maps", "fields", 
                  "shapefiles", "ggmap", "gridExtra")

for(pkg in package_list){
  if(pkg %in% rownames(installed.packages()) == FALSE){
    install.packages(pkg)
  }
}

sapply(package_list, require, character.only = TRUE)


shape_catchment <- st_read("Figure_01_supp/KRF_3.0_Geumgang/KRF_ver3_CATCHMENT_금강수계.shp")
shape <- st_read("Figure_01_supp/KRF_3.0_Geumgang/KRF_ver3_LINE_금강수계.shp")


########################################
##중권역별
########################################
#권역별 출력
region_index <- unique(sapply(as.character(shape$RCH_ID), function(x) substr(x, start=1, stop=4)))[c(1:14)]
#ID의 처음 4자리가 권역별 코드
#3011: 미호천
#3012: 금강공주
#3008: 대청댐
#3007: 보청천
#3005: 초강
#3009: 갑천
#3013: 논산천
#3006: 대청댐상류
#3014: 금강하구언
#3004: 영동천
#3003: 무주남대천
#3001: 용담댐
#3002: 용담댐하류
#3310: 대청댐하류
#3301: 만경강(x)
#3303: 직소천(x)
#3302: 동진강(x)
#3101: 삽교천(x)
#2005: 금강상류부분(필요없음)
#3202: 부남방조제(x)
#1101, 3201: 대호방조제(x)
#3203: 금강서해(x)
#5301, 5302: 주진천(x)

#for(i in 1:length(region_index)){
#  shape_sub <- subset(shape, sapply(as.character(shape@data$RCH_ID), function(x) substr(x, start=1, stop=4)==region_index[i]) )
#  plot(shape, main=region_index[i])
#  plot(shape_sub, col="red", add=T, lwd=3)
#}

name_index <- c("미호천", "금강공주", "대청댐", "보청천", "초강", "갑천", "논산천", "대청댐상류", "금강하구언", "영동천", "무주남대천", "용담댐", "용담댐하류", "대청댐하류")
name_subindex <- list()
name_subindex[[1]] <- c("미호댐", "원남댐", "병천천상류", "병천천하류", "보강천", "미호천중류", "미호천상류", "작천보", "한천", "백곡천", "조천", "석화수위표", "미호천하류", "무심천", "백곡댐" ) #15개
name_subindex[[2]] <- c("세종보", "지천상류", "지천합류후", "어천합류후", "용수천", "백제보", "석성천", "규암수위표", "금천", "논산천합류전", "유구천", "공주보", "대교천", "공주수위표") #14개
name_subindex[[3]] <- c("대청댐", "대청댐조정지", "소옥천상류", "소옥천하류", "대청댐상류") #5개
name_subindex[[4]] <- c("항건천", "보청천중류", "삼가천합류후", "삼가천", "보청천상류", "보청천하류") #6개
name_subindex[[5]] <- c("석천", "초강상류", "초강하류") #3개
name_subindex[[6]] <- c("갑천하류", "유성수위표", "갑천상류", "대전천", "유등천상류", "유등천하류") #6개
name_subindex[[7]] <- c("노성천", "논산천하류", "강경천", "논산천상류", "탑정댐") #5개
name_subindex[[8]] <- c("보청천합류전") #1개
name_subindex[[9]] <- c("입포수위표", "길산천", "금강하구언") #3개
name_subindex[[10]] <- c("영동천", "봉황천하류", "초강합류후", "봉황천상류", "봉황천합류후", "호탄수위표") #6개
name_subindex[[11]] <- c("무주남대천상류", "무주남대천중류", "무주남대천하류") #3개
name_subindex[[12]] <- c("주자천", "용담댐", "정자천", "구량천", "진안천", "장계천합류후", "장계천", "진안천합류후") #8개
name_subindex[[13]] <- c("용담댐하류") #1개
name_subindex[[14]] <- c("미호천합류전", "매포수위표") #2개
region_subindex <- list()

Geum_RCH_ID <- c()
#for(i in 1:length(region_index)){
for(i in 1:1){
  shape_sub <- subset(shape, sapply(as.character(shape$RCH_ID), function(x) substr(x, start=1, stop=4)==region_index[i]) )
  #plot(shape_sub)
  region_subindex[[i]] <-unique(sapply(as.character(shape_sub$RCH_ID), function(x) substr(x, start=1, stop=6)))
  for(j in 1:length(region_subindex[[i]])){
    shape_sub_sub <- subset(shape_sub, sapply(as.character(shape_sub$RCH_ID), function(x) substr(x, start=1, stop=6))==region_subindex[[i]][j])
    #plot(shape_sub_sub, col="red", add=T)
    Geum_RCH_ID_Imsi <- cbind(shape_sub_sub, data.frame(중권역번호=region_index[i], 중권역이름=name_index[i], 소권역번호=region_subindex[[i]][j], 소권역이름=name_subindex[[i]][j]))
    Geum_RCH_ID <- rbind(Geum_RCH_ID, Geum_RCH_ID_Imsi)
  }
}
#Shape_Leng: 리치길이
#RCH_ID: 리치(집수구역) ID

shape <- subset(shape, shape$RCH_ID%in%Geum_RCH_ID$RCH_ID)



# register_stadiamaps(key = "YOUR-API-KEY")

bbox <- c(125.97974, 35.30228, 128.30501, 37.08265)
names(bbox) <- c("left", "bottom", "right", "top")

shp.df <- data.frame(long=double(), lat=double(), group=character(), values=character())
for(i in 1:length(shape$geometry)){
  shp.df <- rbind(shp.df, data.frame(long=shape$geometry[[i]][,1],
                                     lat =shape$geometry[[i]][,2],
                                     group=shape$RCH_ID[i],
                                     values=Geum_RCH_ID$중권역번호[order(Geum_RCH_ID$OBJECTID)][i]))
}

data_fort   <- fortify(shp.df)


data_merged_miho <- data_fort[which(data_fort$values=="3011"),]

########################################
##Plot stations (2019년 6월 19일자)
########################################
places_matrix_new <- read.csv("Figure_01_supp/eventplace.csv",header=T)
places <- data.frame(name= places_matrix_new$X,x=places_matrix_new$long, y=places_matrix_new$lat)
places <- places[which(places_matrix_new$중권역=="미호천"),]



########################################
##Plot stations
########################################
mapkorea <- map

#ggmap(mapkorea)

Geum_catchment <- subset(shape_catchment, substr(shape_catchment$CAT_ID, start=1, stop=4)=="3011")
Geum.id <- substr(Geum_catchment$CAT_ID, start=1, stop=4)
Geum.miho.id <- substr(Geum_catchment$CAT_ID, start=1, stop=6)
#plot(Geum_catchment)
Geum_catchment_border <- st_union(Geum_catchment)
Geum_subcatchment <- c()
for(i in 1:length(unique(Geum.miho.id))){
  p1 <- Geum_catchment[which(Geum.miho.id==unique(Geum.miho.id)[i]),]
  Geum_subcatchment <- rbind(Geum_subcatchment, mutate(st_sf(st_union(p1)), centrum=T))
}


########################################
##산업단지 위치 plotting
########################################
industrial_data <- st_read("Figure_01_supp/Danji/DAM_PDAN.shp")

#sigungu <- readOGR("~/Dropbox/Maps/KOR_SGG(WGS)/TL_SCCO_SIG.shp")
sigungu <- st_read("Figure_01_supp/KOR_SGG(WGS)/TL_SCCO_SIG.shp")

sigungu_index <- c(239:242, #청주
                   168, #세종
                   223, #천안시(동남구)
                   249, #진천
                   251, #음성
                   53, #안성
                   250, #괴산
                   248 #증평
)

sigungu_pop <- c(171129, #청주 상당구
                 216782, #청주 서원구
                 253912, #청주 흥덕구
                 193767, #청주 청원구
                 280100, #세종
                 255878, #천안 동남구
                 73677, #진천
                 97306, #음성
                 182786, #안성
                 39054, #괴산
                 37783 #증평
)


sigungu_Geum <- subset(sigungu, sigungu$SIG_CD%in%sigungu$SIG_CD[sigungu_index])

# compute the central part of the polygon
sigungu_middle_loc <- c()
for(i in 1:length(sigungu_index)){
  if(length(sigungu_Geum$geometry[[i]][1])==1){
    if(i==24){
      sigungu_middle_loc <- rbind(sigungu_middle_loc, colQuantiles(sigungu_Geum$geometry[[i]]@Polygons[[1]]@coords, prob=0.25))
    }else{
      sigungu_middle_loc <- rbind(sigungu_middle_loc, colMeans(sigungu_Geum$geometry[[i]][1][[1]][[1]]))
    }
  }else{
    length.each <- sapply(sigungu_Geum$geometry[[i]]@Polygons, function(x) nrow(x@coords))
    max.index <- which.max(length.each)
    if(i==13){
      #완주군
      sigungu_middle_loc <- rbind(sigungu_middle_loc, colQuantiles(sigungu_Geum$geometry[[i]]@Polygons[[max.index]]@coords, prob=0.75))
    }else{
      sigungu_middle_loc <- rbind(sigungu_middle_loc, colMeans(sigungu_Geum$geometry[[i]]@Polygons[[max.index]]@coords))
    }
  }
}

sigungu_pop_new <- sigungu_pop/100000

industrial_data_t <- st_transform(industrial_data$geometry, 4326)
industrial_data_transform <- st_coordinates(industrial_data_t)
industrial_data_transform <- data.frame(x= industrial_data_transform[,1],
                                        y= industrial_data_transform[,2],
                                        DAN_ID =industrial_data$DAN_ID,
                                        DAN_NAME =industrial_data$DAN_NAME, 
                                        DANJI_TYPE = industrial_data$DANJI_TYPE)
#industrial_data_transform <- spTransform(industrial_data, "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
#industrial_data_transform@coords
#industrial_data_transform@data$DAN_NAME
#industrial_data_transform@data$DANJI_TYPE #


pointinpoly_index <- rep(0, nrow(industrial_data_transform))

for(i in 1:length(sigungu_Geum$geometry)){
    pointinpoly_index  <- pointinpoly_index  + point.in.polygon(industrial_data_transform[,1], industrial_data_transform[,2], st_coordinates(Geum_catchment_border)[,1], st_coordinates(Geum_catchment_border)[,2])
}

danji_index <- as.numeric(industrial_data_transform$DANJI_TYPE[which(pointinpoly_index>0)])

catchment.test <- fortify(Geum_catchment)

catchment_index <- rep(0, nrow(industrial_data_transform))

for(i in 1:length(Geum_catchment$OBJECTID)){
    catchment_index  <- catchment_index  + point.in.polygon(industrial_data_transform[,1], industrial_data_transform[,2], Geum_catchment$geometry[[i]][[1]][[1]][,1], Geum_catchment$geometry[[i]][[1]][[1]][,2])
}


danji_index2 <- as.numeric(industrial_data_transform$DANJI_TYPE[which(catchment_index>0)])
for(k in 1:length(danji_index2)){
  if(danji_index2[k]==1){
    danji_index2[k] <-"National"
  }else if(danji_index2[k]==2){
    danji_index2[k] <- "General"
  }else if(danji_index2[k]==3){
    danji_index2[k] <- "Urban"
  }else if(danji_index2[k]==4){
    danji_index2[k] <- "Rural"
  }
}

danji_data_frame <- data.frame(x=industrial_data_transform[which(catchment_index>0),1], y=industrial_data_transform[which(catchment_index>0),2], Type=danji_index2)


#plot.test.df <- join(plot.test, sigungu_Geum@data)
#
#ggplot(data=plot.test, mapping=aes(x="long", y="lat"))

data_Sigungu <- data.frame()


sigungu_middle_loc2 <- data.frame(long=sigungu_middle_loc[,1],
                                  lat=sigungu_middle_loc[,2],
                                  size=sigungu_pop_new[match(sigungu_Geum$SIG_ENG_NM, sigungu$SIG_ENG_NM[sigungu_index])])

danji_index <- as.numeric(industrial_data_transform$DANJI_TYPE[which(pointinpoly_index>0)])
danji_index[which(danji_index==1)] <- "National"
danji_index[which(danji_index==2)] <- "General"
danji_index[which(danji_index==3)] <- "Urban"
danji_index[which(danji_index==4)] <- "Rural"

industrial_data_transform2 <- data.frame(long=industrial_data_transform[which(pointinpoly_index>0),1],
                                         lat=industrial_data_transform[which(pointinpoly_index>0),2],
                                         danji_index = as.factor(danji_index))


bbox2 <- bbox
bbox2 <- c(127, 36.4, 128.1, 37.2)
map2 <-get_stadiamap(bbox2, maptype = "stamen_terrain_background") 
p1 <- ggmap(map2) + annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "white", alpha = 0.6)  + geom_sf(data=sigungu_Geum, inherit.aes = FALSE, fill=NA, color="black") + geom_path(data=data_merged_miho,size=1, aes(x=long,y=lat,group=group),color="skyblue") + geom_point(data = sigungu_middle_loc2 , aes(x = long, y = lat, size=size)) + labs(title="(a) Population", x = "Longitude", y = "Latitude", size="Unit: 100,000M") + theme(legend.position = "bottom")



bbox3 <- bbox
bbox3 <- c(127.1, 36.4, 127.7, 37.2)
map <-get_stadiamap(bbox3, maptype = "stamen_terrain_background") 
p2 <- ggmap(map) + 
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "white", alpha = 0.6)  + geom_sf(data=Geum_subcatchment, inherit.aes = FALSE, fill=NA, color="black") + 
  geom_path(data=data_merged_miho,size=1, aes(x=long,y=lat,group=group),color="skyblue") + 
  geom_point(data = industrial_data_transform2, aes(x = long, y = lat, shape = danji_index, color = danji_index), size = 3) + 
  labs(title="(b) Industrial Areas", x = "Longitude", y = "Latitude", color="Type", shape="Type") 
#+ theme(legend.position = "bottom")


pdf("plots/Figure01.pdf", width = 12, height = 7)


grid.arrange(p1, p2, nrow = 1, ncol = 2)
dev.off()

