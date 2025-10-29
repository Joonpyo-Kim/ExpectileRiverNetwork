############################################
########## Preliminary functions ###########
############################################

# Some of these functions are from github repository of one of the authors of O'Donnell et al. (2014). 
# Please see github.com/alastairrushworth/smnet. 

# Truncated p-th power function
tpower <- function(x, t, p)  (x - t) ^ p * (x > t)


# Construct B-spline basis
bbase 	<- 	function(x, xl = min(x), xr = max(x), nseg = 10, deg = 3){
  dx 	<- 	(xr - xl) / nseg
  knot<- 	seq(xl - deg * dx, xr + deg * dx, by = dx)
  P 	<- 	outer(x, knot, tpower, deg)
  n 	<- 	dim(P)[2]
  D 	<- 	diff(diag(n), diff = deg + 1) / (gamma(deg + 1) * dx ^ deg)
  B 	<- 	(-1) ^ (deg + 1) * P %*% t(D)
  B}


get.dimension<-function(model.type, n.segments, knts){
  no<-0
  dims<-vector("numeric")
  if("c" %in% model.type){
    dims<-c(1, dims)
  }
  if("m" %in% model.type){
    dims<-c(dims, n.segments)
  }
  if("s" %in% model.type){
    no<-no+1
    dims<-c(dims, knts[no])
    if("si" %in% model.type){
      dims<-c(dims, knts[no]*n.segments)
    }
  }
  if("t" %in% model.type){
    no<-no+1 
    dims<-c(dims, knts[no])
    if("ti" %in% model.type){
      dims<-c(dims, knts[no]*n.segments)
    }
    if("ts" %in% model.type){
      dims<-c(dims, prod(knts))
    }
  }
  dims
}

box3.prod<-function(a, b)
{
  a1<-kronecker.spam(a, as.spam(matrix(1, nrow = 1, ncol = ncol(b))))
  b1<-kronecker.spam(as.spam(matrix(1, nrow = 1, ncol = ncol(a))), b)
  a1*b1
}

get.blockmat<-function(Mat, lengths, i)
{
  cum.dims<-cumsum(lengths)
  list.mats<-as.list(1:length(lengths))
  mat.construct<-function(n, Mat, i)
  {
    if(!n == i)
    {
      mat<-Matrix(0, nrow = lengths[n], ncol = lengths[n])
    }
    if(n == i)
    {
      mat<-Mat
    }
    mat
  }
  
  blocks<-lapply(list.mats, mat.construct, Mat = Mat, i = i)
  bdiag(blocks)
}

crossprodspam<-function(X){
  t(X) %*% X  
}

get.PTP<-function(P, dim.list, i){
  make.sparse<-function(n) as.spam(matrix(0, n, n))        
  Mat.list<-lapply(dim.list, make.sparse)
  Mat.list[[i]]<-t(P) %*% P
  b.diag<-function(L){
    diag.spam<-L[[1]]
    for(i in 2:length(L)){
      diag.spam<-bdiag.spam(diag.spam, L[[i]]) #diag.spam: creates a sparse block-diagonal matrix 
    }
    diag.spam
  }
  b.diag(Mat.list)
}	

b.diag<-function(L){
  diag.spam<-L[[1]]
  for(i in 2:length(L)){
    diag.spam<-bdiag.spam(diag.spam, L[[i]])
  }
  diag.spam
}

############################################
##### Fitting expectile additive model #####
############################################

ExpecSTnet <- function(realweights, adjacency, TweedData, TweedPredPoints, tau = 0.5, 
                       penalties =  NULL, plot.fig=FALSE, station=NULL, 
                       log.y=FALSE, model.type=NULL){
  
  if(log.y==TRUE){
    response <- log(TweedData$nitrate)
  }else{
    response <- TweedData$nitrate
  }
  
  if(is.null(model.type)){
    model.type = c("c", "m", "s", "si", "t", "ti", "ts")
  }
  
  penalty_param_matrix <- matrix(c(
  1,   25,  195,  180,   90,  170,  185,  120,   95,
 15,   60,  175,  180,   40,  175,  180,   10,   50,
 15,   60,  150,  165,   25,  145,  170,   10,   55,
 10,   65,   65,  130,   30,    5,   40,   15,   55,
 15,   80,   35,  110,   40,    5,   25,   10,   85,
 15,   85,   40,   90,   40,    5,   25,   10,   75,
 10,   95,   35,   95,   35,    5,   20,   10,   75,
 15,  105,   35,   95,   60,    5,   25,    5,   85,
 15,  105,   50,   85,   60,    5,   25,   10,   80,
 15,  110,   80,   95,   70,    5,   45,   15,   80,
 20,  110,  145,  155,   70,  125,  175,   25,   65,
 35,   85,  225,  235,   50,  215,  230,   55,   65,
160,  130,  185,  245,  125,  190,  270,  130,  120
  ), ncol = 9, byrow = TRUE)
  
  if(is.null(penalties)){
    if(0.1 <= tau & tau <= 0.9){
      penalties <- penalty_param_matrix[2+floor(tau*10),]
    }else if(tau < 0.05){
      penalties <- penalty_param_matrix[1,]
    }else if(tau < 0.1){
      penalties <- penalty_param_matrix[2,]
    }else if(tau > 0.95){
      penalties <- penalty_param_matrix[13,]
    }else if(tau > 0.9){
      penalties <- penalty_param_matrix[12,]
    }
  }
  
  if(is.null(station)){
    station=38
  }
  
  
  
  response.locs <- TweedData$location
  knts<-c(10,10) 
  dates = as.character(TweedData$date)
  n.segments<-nrow(adjacency) 
  p.dims<-n.segments
  dims<-get.dimension(model.type = model.type, n.segments = n.segments, knts = knts) #1 113 10 1130 10 1130 100 (c, m, s, si, t, ti, ts)
  c.dims<-cumsum(dims)
  n.par<-sum(dims) 
  X.list<-vector("list") 
  
  
  #### -------------- CONSTRUCT INTERCEPT COLUMN OF MODEL MATRIX ("c") -------------- ####
  if("c"%in%model.type){
    X.list<-c(X.list, list(ones = as.spam(matrix(1, nrow = length(response),  ncol = 1))))
    model.mat<-X.list$ones
  }
  
  
  
  #### -------------- CONSTRUCT NETWORK MODEL MATRIX AND NETWORK PENALTY ("m") -------------- ####
  n<-2
  p.n<-0
  if("m"%in%model.type){
    X<-matrix(0, nrow = length(response), ncol = n.segments) 
    for(i in 1:length(response)){X[i, response.locs[i]] <- 1} 
    X.list<-c(X.list, list(spatial = as.spam(X)))
    model.mat<-cbind.spam(model.mat, X.list$spatial) 
    Z.list<-vector("list")
    rm(list = c("X"))
    
    # calculate difference matrix D on all pairs of adjacenct stream segments
    D1<-matrix(0, n.segments, n.segments) 
    D2<-matrix(0, n.segments, n.segments) 
    for(i in 1:n.segments){ 
      a <- which(adjacency[,i] == 1) 
      if(length(a) > 0)
      {
        D1[i, i]<- -realweights[a[1],1]
        D1[i, a[1]]<-realweights[a[1],1] 
      }	
      if(length(a) > 1)
      {
        D2[i, i]<--realweights[a[2],1] 
        D2[i, a[2]]<-realweights[a[2],1] 
      }
    }
    D<-rbind(D1, D2) #((113*2)*113) matrix
    rm.d<-which(rowSums(abs(D)) == 0)
    P_spat<-as.spam(D[-rm.d,]) 
    PP1<-get.PTP(P = P_spat, dim.list = dims, i = match("m", model.type)) 
    
    if(length(dims) == 1) PP1<-crossprod(P_spat) 
    PTP <- PP1*penalties[1] 
  }else{
    PTP <- as.spam(matrix(0, nrow=sum(dims), ncol=sum(dims))) 
  }
  
  
  
  #### -------------- CONSTRUCT SEASONAL COMPONENT AND SEASONAL PENALTY ("s") -------------- ####
  n<-n+1 #n<-3
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
  
  if("s"%in%model.type){
    # cSplineDes: Uses splineDesign to set up the model matrix for a cyclic B-spline basis.
    Bas<-cSplineDes(decimal.day, knots = (0:knts[1])/knts[1], ord=4) #dim(Bas): 6266 10
    X.list<-c(X.list, list(seasonal = as.spam(Bas)))
    model.mat<-cbind.spam(model.mat, X.list$seasonal) #6266*124 matrix
    
    # circular penalty matrix
    S<-diff(diag(dim(Bas)[2]), diff = 2) 
    last.row<-S[nrow(S), ] 
    row1<-c(last.row[ncol(S)], last.row[-ncol(S)]) 
    S<-as.matrix(rbind(S, row1))
    
    last.row<-S[dim(S)[1],] 
    row1<-c(last.row[dim(S)[2]], last.row[-dim(S)[2]])
    P_seasonal<-as.spam(as.matrix(rbind(S, row1)))
    
    P.seasonal<-get.PTP(P = P_seasonal, dim.list = dims, i = match("s", model.type)) #2494*2494 
    PTP<-PTP+P.seasonal*penalties[n-1] 
  }
  
  
  
  #### -------------- CONSTRUCT SEASONAL-NETWORK INTERACTION COMPONENT ("si") -------------- ####
  n<-n+1 #n<-4
  if("si"%in%model.type){
    seas_spat<-box3.prod(X.list$spatial, X.list$seasonal)
    X.list<-c(X.list, list(seas_spat)) 
    model.mat<-cbind.spam(model.mat, seas_spat) #6266*1254 matrix
    
    P_si1<-kronecker.spam(P_spat, diag.spam(ncol(P_seasonal))) #1120*1130 matrix
    P.seas1<-get.PTP(P = P_si1, dim.list = dims, i = match("si", model.type))*penalties[n-1+p.n]
  }
  p.n<-p.n+1
  if("si"%in%model.type){
    P_si2<-kronecker.spam(diag.spam(ncol(P_spat)), P_seasonal)
    P.seas2<-get.PTP(P = P_si2, dim.list = dims, i = match("si", model.type))*penalties[n-1+p.n] 
    PTP<-PTP + P.seas1 + P.seas2
  }
  
  
  
  #### -------------- CONSTRUCT LONG TERM TREND COMPONENT AND PENALTY ("t") -------------- ####
  n<-n+1 #n<-5
  if("t"%in%model.type){
    years <- year(dates) + (decimal_date(as.Date(dates2))-2013)
    X.list<-c(X.list, list(trend = as.spam(bbase(years, nseg = (knts[2] - 3))))) 
    model.mat<-cbind.spam(model.mat, X.list$trend) #model.mat : 6266 1264
    
    P_trend<-as.spam(diff(diag(ncol(X.list$trend)), diff = 1)) #dim(P_trend) : 9 * 10
    P.trend<-get.PTP(P = P_trend, dim.list = dims, i = match("t", model.type))*penalties[n-1+p.n]
    PTP<-PTP + P.trend
  }
  
  
  
  #### -------------- CONSTRUCT NETWORK-LONG TERM TREND INTERACTION AND PENALTY ("ti") -------------- ####
  if("ti"%in%model.type){
    n<-n+1 #n<-6
    tren_spat<-box3.prod(X.list$spatial, X.list$trend)
    X.list<-c(X.list, list(tren_spat))
    model.mat<-cbind.spam(model.mat, tren_spat) #model.mat: 6266*2394 matrix
    
    P_ti1<-kronecker.spam(P_spat, diag.spam(ncol(P_trend)))
    P_ti2<-kronecker.spam(diag.spam(ncol(P_spat)),  P_trend)
    P.ti1<-get.PTP(P = P_ti1, dim.list = dims, i = match("ti", model.type))*penalties[n-1+p.n] 
  }
  p.n<-p.n+1 #p.n<-2
  if("ti"%in%model.type){
    P.ti2<-get.PTP(P = P_ti2, dim.list = dims, i = match("ti", model.type))*penalties[n-1+p.n]
    PTP<-PTP + P.ti1 + P.ti2
  }
  
  
  
  #### -------------- CONSTRUCT TREND-SEASONAL INTERACTION COMPONENT ("ts") -------------- ####
  if("ts"%in%model.type){
    n<-n+1 #n<-7
    tren_seas <- box3.prod(X.list$trend, X.list$seasonal)
    model.mat<-cbind.spam(model.mat, tren_seas) #model.mat : 6266*2494 matrix
    
    X.list<-c(X.list, list(tren_seas))
    P_ts1<-kronecker.spam(P_trend, diag.spam(ncol(P_seasonal)))
    P_ts2<-diag(ncol(P_trend)) %x% P_seasonal
    P.ts1<-get.PTP(P = as.spam(P_ts1), dim.list = dims, i = match("ts", model.type))*penalties[n-1 + p.n]
  }
  p.n<-p.n+1 #p.n<-3
  if("ts"%in%model.type){
    P.ts2<-get.PTP(P = as.spam(P_ts2), dim.list = dims, i = match("ts", model.type))*penalties[n-1 + p.n] 
    PTP<-PTP + P.ts1 + P.ts2
  }
  
  # gc: Garbage Collection
  gc()
  
  
  
  #### -------------- FIT THE MODEL -------------- ####
  ident<-b.diag(lapply(X.list, crossprodspam))
  W <- diag.spam(dim(model.mat)[1])  
  weight <- diag(W)       
  XTX<-t(model.mat) %*% W %*% model.mat    
  info<-XTX + PTP  +  0.001*ident + diag.spam(c(0.0001, rep(0.1, nrow(XTX)-1))) 
  U<-chol.spam(info, pivot = TRUE) 
  poo <- as.matrix(t(model.mat) %*% W) %*% response
  poo <- as.spam(poo)
  
  beta_hat <- spam::backsolve(U, spam::forwardsolve(U, poo)) 
  fit <- model.mat %*% beta_hat 
  
  weight.new <- as.vector(tau * (response >= fit) + (1-tau) * (response < fit)) 
  # W.new <- diag.spam(weight.new)   
  iter <- 1
  
  # expectile regression
  while(sum((weight-weight.new)^2)!=0){    
    weight <- weight.new      
    W <- diag.spam(weight)      
    iter <- iter + 1
    XTX <-t(model.mat) %*% W %*% model.mat
    info<-XTX + PTP  +  0.001*ident + diag.spam(c(0.0001, rep(0.1, nrow(XTX)-1)))
    U<-chol.spam(info, pivot = TRUE)
    poo <- t(as.matrix(model.mat)) %*% W %*% response
    poo <- as.spam(poo)          
    beta_hat<-spam::backsolve(U, spam::forwardsolve(U, poo))
    fit<-model.mat %*% beta_hat
    weight.new <- as.vector(tau * (response >= fit) + (1-tau) * (response < fit))  
    # W.new <- diag.spam(weight.new)
    if(iter == 500){     
      print(paste("procedure does not converge"));      
      break;     
    }    
  }
  
  W.new <- diag.spam(weight.new)
  resids <- response - fit
  
  
  
  #### -------------- AIC -------------- ####
  if(tau == 0.5){
    vec<-spam::forwardsolve(U, t(model.mat))
    pdof<-rowSums(vec^2)
    dof<-sum(pdof)
    sigma.sq<-sum((response - fit)^2)/(length(response) - dof)
    AICc<-log(sigma.sq) + 1 + ((2 + 2*dof)/(length(response) - dof -2))
    
    # expectile version
    w <- as.vector(tau * (response >= fit) + (1-tau) * (response < fit))
    sse <- sum(w * ((as.vector(resids)^2)))
    AIC_exp <- log(sse/length(response)) + 2*dof/length(response)
  }else{
    W_mat <- diag.spam(as.vector(tau*(resids >= 0) + (1-tau)*(resids < 0)))
    XTWX<-t(model.mat) %*% W_mat %*% model.mat
    info<-XTWX + PTP  +  0.001*ident + diag.spam(c(0.0001, rep(0.1, nrow(XTWX)-1)))
    U<-chol.spam(info, pivot = TRUE)
    vec<-spam::forwardsolve(U, t(model.mat))
    pdof<-rowSums(vec^2)
    dof<-sum(pdof)
    sigma.sq<-sum((response - fit)^2)/(length(response) - dof)
    AICc<-log(sigma.sq) + 1 + ((2 + 2*dof)/(length(response) - dof -2))
    
    # expectile version
    w <- as.vector(tau * (response >= fit) + (1-tau) * (response < fit))
    sse <- sum(w * ((as.vector(resids)^2)))
    AIC_exp <- log(sse/length(response)) + 2*dof/length(response)
  }
  
  ##### ----- main & interaction ----- #####
  parameters <- as.vector(beta_hat) 
  locs<-cumsum(dims) 
  if(length(model.type) == 1){locs<-c(1, locs)}
  comps<-vector("list")	
  fit.point <- parameters[(locs[match("m", model.type)-1]+1):locs[match("m", model.type)]]
  beta_1 <- parameters[1]
  
  comps <- vector("list", length(model.type))
  min.yr <- min(year(dates))
  max.yr <- max(year(dates))
  yrs <- seq(min.yr, max.yr)
  n.yrs <- length(yrs)
  day.seq <- seq(min.yr + 5/365, max.yr + 1, length.out = (365/5)*n.yrs) 
  
  # local segment specific means
  if("m"%in%model.type){
    seg.means<-parameters[(locs[match("m", model.type)-1]+1):locs[match("m", model.type)]]
  }else{
    seg.means<-rep(parameters[1], length(station)) 
  }
  
  station.mean<-seg.means[station]
  
  main <- 0; inter <- 0;
  
  # catchment wide constant
  if("c"%in%model.type){
    C<-parameters[1]
    C<-matrix(C, nrow = n.yrs, ncol = 365/5)
    main <- main+C
  }
  
  # spatial component
  if("m"%in%model.type){
    M<-matrix(station.mean, nrow = n.yrs, ncol = 365/5)
    main <- main+M
  }
  
  # seasonal component
  if("s"%in%model.type){
    
    s.bas<-cSplineDes(seq(5,365,by=5)/365, knots = 0:knts[1]/knts[1], ord=4)
    s.par<-parameters[(locs[match("s", model.type) - 1] + 1):locs[match("s", model.type)]] 
    S<-s.par %*% t(s.bas)
    S<-matrix(S, nrow = n.yrs, ncol = 365/5, byrow = T)	
    main <- main+S
  }
  
  # seasonal spatial interaction component
  if("si"%in%model.type){
    all.si.par<-parameters[(locs[match("si", model.type) - 1] + 1):locs[match("si", model.type)]]
    #si.par.stat<-all.si.par[(station*ncol(s.bas) + 1):((station+1)*ncol(s.bas))]
    si.par.stat<-all.si.par[((station*ncol(s.bas) + 1)-10):((station+1)*ncol(s.bas)-10)]
    si.add<-s.bas %*% si.par.stat
    SI<-matrix(si.add, nrow = n.yrs, ncol = 365/5, byrow = T)	
    inter <- inter+SI
  }
  
  # trend additive component
  if("t"%in%model.type){
    yers<-seq(min.yr, max.yr+1, by=5/365)[-1]
    t.bas<-bbase(yers, nseg = (knts[2] - 3))
    t.par<-parameters[(locs[match("t", model.type) - 1] + 1):locs[match("t", model.type)]]
    Tr1<-t.bas %*% t.par 
    Tr<-matrix(Tr1, nrow = n.yrs, ncol = 365/5, byrow = T)
    main <- main+Tr
  }
  
  # trend spatial interaction component
  if("ti"%in%model.type){
    all.ti.par<-parameters[(locs[match("ti", model.type) - 1] + 1):locs[match("ti", model.type)]]
    #ti.par.stat<-all.ti.par[(station*ncol(t.bas) + 1):((station+1)*ncol(t.bas))]
    ti.par.stat<-all.ti.par[((station*ncol(s.bas) + 1)-10):((station+1)*ncol(s.bas)-10)]
    ti.add<-t.bas %*% ti.par.stat
    TI<-matrix(ti.add, nrow = n.yrs, ncol = 365/5, byrow = T)	
    inter <- inter+TI
  }
  
  # trend seasonal interaction
  if("ts"%in%model.type){
    all.ts.par<-parameters[(locs[match("ts", model.type) - 1] + 1):locs[match("ts", model.type)]]
    s.bas.ext<-matrix(1, ncol = 1,  nrow = n.yrs) %x% s.bas
    ts.bas<-box3.prod(as.spam(t.bas), as.spam(s.bas.ext))
    ts.add<-ts.bas %*% all.ts.par
    TS<-matrix(ts.add, nrow = n.yrs, ncol = 365/5, byrow = T)
    inter <- inter+TS
  }
  inter <- inter
  main <- main
  # main <- S + Tr + C + M
  # inter <- TI + SI + TS
  
  if(plot.fig==TRUE){
    
    #### -------------- PLOT THE SPATIAL COMPONENT FOR ALL TIME POINTS, INCLUDING INTERCEPT -------------- ####
    if("m"%in%model.type){
      comps<-c(comps, list(spatial = parameters[(locs[match("m", model.type)-1]+1):locs[match("m", model.type)]]))
      fit.point <-comps$spatial
      main.c<-paste("Spatial component ", " DoF = ", 
                    round(sum(pdof[((locs[match("m", model.type)-1]+1):locs[match("m", model.type)])]), 1), sep="")
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
      # title(main.c)
      title(paste0("Spatial Component for ", tau*100, "% expectile"))
      
    }
    
    
    
    #### -------------- PLOT THE SEASONAL COMPONENT WITH PARTIAL RESIDUALS -------------- ####
    if("s"%in%model.type){
      
      
      s.par<-parameters[(locs[match("s", model.type) - 1] + 1):locs[match("s", model.type)]]
      seas<-(1:365)/365
      seas.bas<-cSplineDes(1:365/365, knots = 0:knts[1]/knts[1], ord=3)
      seas.comp<-seas.bas %*% s.par
      main.s<-"Seasonal component"
      main.s<-paste(main.s, " DoF = ", round(sum(pdof[((locs[match("s", model.type)-1]+1):locs[match("s", model.type)])]), 1), sep="") 
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
    }
    
    
    
    #### -------------- PLOT THE TREND COMPONENT -------------- ####
    if("t"%in%model.type){
      years<-year(dates)
      #n.t<-which(model.type == "t") 일단 안 쓰임
      t.par<-parameters[(locs[match("t", model.type) - 1] + 1):locs[match("t", model.type)]]
      tren<-seq(min(years), max(years))
      tren.bas<-bbase(tren, xl = min(years), xr = max(years), nseg = (knts[2] - 3))
      tren.comp<-tren.bas %*% t.par 
      
      main.t<-"Trend component"
      main.t<-paste(main.t, " DoF = ", round(sum(pdof[((locs[match("t", model.type)-1]+1):locs[match("t", model.type)])]), 1), sep="")
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
    }
    
    
    
    #### -------------- PLOT THE INTERACTION BETWEEN TREND AND SEASONALITY -------------- ####
    if("ts"%in%model.type){
      min.yr<-min(year(dates))
      max.yr<-max(year(dates))
      yers<-seq(min.yr, max.yr+1, length.out = 50)
      n.ts<-which(model.type == "ts")
      ones<-rep(1, length = 50)
      s.bas<-ones %x% cSplineDes(seq(1,365,length.out=50)/365, knots = 0:knts[1]/knts[1], ord=4)
      t.bas<-bbase(yers, nseg = (knts[2] - 3)) %x% ones
      ts.bas2<-box3.prod(t.bas, s.bas)
      ts.par<-parameters[(locs[n.ts - 1] + 1):locs[n.ts]]
      ts.fitted<-matrix(as.vector(ts.bas2 %*% ts.par), nrow = 50, byrow = T)
      
      # plot
      par(mar = c(5, 4, 4, 2))
      filled.contour(ts.fitted, xlab = "Trend", ylab = "Seasonal",axes = T)
    }
    
    
    
    #### -------------- PLOT THE FITTED VALUES WITH TIME AT A SINGLE MONITORING STATION -------------- ####
    
    
    comps <- vector("list", length(model.type))
    min.yr <- min(year(dates))
    max.yr <- max(year(dates))
    yrs <- seq(min.yr, max.yr)
    n.yrs <- length(yrs)
    day.seq <- seq(min.yr + 5/365, max.yr + 1, length.out = (365/5)*n.yrs) 
    
      obs.ind <-which(response.locs == station)
      dates.full <-dates[obs.ind]
      dates.stat <-decimal_date(as.Date(dates.full))
      vals.stat <-response[obs.ind]
      
      # plot
      if(!("t"%in%model.type) & !("ti"%in%model.type) & !("ts"%in%model.type) & ("s"%in%model.type) ){
        par(mar = c(4,4,2,0.5))
        dates.stat2 <- dates.stat - floor(dates.stat)
        plot(dates.stat2, vals.stat, pch = ".", cex = 1, type = "n",ylim = range(vals.stat),
             xlab = "Time", ylab = "log(Nitrate concentrations)",cex.lab = 1, cex.axis=1, main=paste("Segment ", station, '(W/O TREND COMPONENT)'))
        day.seq2 <- seq(0+5/365, 0+1, length.out=365/5)
        lines(day.seq2, c(t(main))[c(1:length(day.seq2))], type = "l", lwd=2) #col = "limegreen", 
        lines(day.seq2, c(t(inter))[c(1:length(day.seq2))]+ c(t(main))[c(1:length(day.seq2))], type = "l", lwd=2, lty=2) #col = "tomato", 
        points(dates.stat2, vals.stat)
        legend("bottomleft", lty=c(1,2), c("Main", "With Interaction"), lwd=c(2,2)) # col=c("limegreen", "tomato"),
      }else{
        par(mar = c(4,4,2,0.5))
        plot(dates.stat, vals.stat, pch = ".", cex = 1, type = "n",ylim = range(vals.stat),
             xlab = "Time", ylab = "log(Nitrate concentrations)",cex.lab = 1, cex.axis=1, main=paste("Segment ", station))
        lines(day.seq, c(t(main)), type = "l",  lwd=2, lty=2) #col = "limegreen",
        lines(day.seq, c(t(inter))+ c(t(main)), type = "l",lwd=2, lty=1) #col = "tomato", 
        points(dates.stat, vals.stat)
        legend("topright", lty=c(1,2), c("With Interaction","Without Interaction"), lwd=c(2,2)) #col=c("limegreen", "tomato"), 
      }
        
    # }
  }
  
  main <- main
  inter <- inter
  

  return(list(AICc=AICc, AIC_exp=AIC_exp, U=U, info = info, W.new = W.new, weight.new = weight.new , 
                beta_hat=beta_hat, fit=fit, resids=resids, model.type = model.type, model.mat=model.mat, XTX=XTX, PTP=PTP, 
                penalties=penalties, realweights=realweights, adjacency=adjacency, TweedData=TweedData, 
                TweedPredPoints=TweedPredPoints,
                response = response, tau = tau,
                main = main, inter = inter))
  
  
}

ExpecSTnet_plot.estimated <- function(station, smnet_obj, plot.fig = FALSE, legend = TRUE){
  response.locs <- smnet_obj$TweedData$location
  knts<-c(10,10) 
  dates = as.character(smnet_obj$TweedData$date)
  n.segments<-nrow(smnet_obj$adjacency) 
  p.dims<-n.segments
  model.type <- smnet_obj$model.type
  dims<-get.dimension(model.type = model.type, 
                      n.segments = n.segments, 
                      knts = knts) 
  c.dims<-cumsum(dims)
  n.par<-sum(dims) 
  X.list<-vector("list") 
  
  parameters <- as.vector(smnet_obj$beta_hat) 
  locs<-cumsum(dims) 
  if(length(model.type) == 1){locs<-c(1, locs)}
  comps<-vector("list")	
  fit.point <- parameters[(locs[match("m", model.type)-1]+1):locs[match("m", model.type)]]
  beta_1 <- parameters[1]
  
  comps <- vector("list", length(model.type))
  min.yr <- min(year(dates))
  max.yr <- max(year(dates))
  yrs <- seq(min.yr, max.yr)
  n.yrs <- length(yrs)
  day.seq <- seq(min.yr + 5/365, max.yr + 1, length.out = (365/5)*n.yrs) 
  
  # local segment specific means
  if("m"%in%model.type){
    seg.means<-parameters[(locs[match("m", model.type)-1]+1):locs[match("m", model.type)]]
  }else{
    seg.means<-rep(parameters[1], length(station)) 
  }
  
  station.mean<-seg.means[station]
  
  main <- 0; inter <- 0;
  
  # catchment wide constant
  if("c"%in%model.type){
    C<-parameters[1]
    C<-matrix(C, nrow = n.yrs, ncol = 365/5)
    main <- main+C
  }
  
  # spatial component
  if("m"%in%model.type){
    M<-matrix(station.mean, nrow = n.yrs, ncol = 365/5)
    main <- main+M
  }
  
  # seasonal component
  if("s"%in%model.type){
    
    s.bas<-cSplineDes(seq(5,365,by=5)/365, knots = 0:knts[1]/knts[1], ord=4)
    s.par<-parameters[(locs[match("s", model.type) - 1] + 1):locs[match("s", model.type)]] 
    S<-s.par %*% t(s.bas)
    S<-matrix(S, nrow = n.yrs, ncol = 365/5, byrow = T)	
    main <- main+S
  }
  
  # seasonal spatial interaction component
  if("si"%in%model.type){
    all.si.par<-parameters[(locs[match("si", model.type) - 1] + 1):locs[match("si", model.type)]]
    #si.par.stat<-all.si.par[(station*ncol(s.bas) + 1):((station+1)*ncol(s.bas))]
    si.par.stat<-all.si.par[((station*ncol(s.bas) + 1)-10):((station+1)*ncol(s.bas)-10)]
    si.add<-s.bas %*% si.par.stat
    SI<-matrix(si.add, nrow = n.yrs, ncol = 365/5, byrow = T)	
    inter <- inter+SI
  }
  
  # trend additive component
  if("t"%in%model.type){
    yers<-seq(min.yr, max.yr+1, by=5/365)[-1]
    t.bas<-bbase(yers, nseg = (knts[2] - 3))
    t.par<-parameters[(locs[match("t", model.type) - 1] + 1):locs[match("t", model.type)]]
    Tr1<-t.bas %*% t.par 
    Tr<-matrix(Tr1, nrow = n.yrs, ncol = 365/5, byrow = T)
    main <- main+Tr
  }
  
  # trend spatial interaction component
  if("ti"%in%model.type){
    all.ti.par<-parameters[(locs[match("ti", model.type) - 1] + 1):locs[match("ti", model.type)]]
    #ti.par.stat<-all.ti.par[(station*ncol(t.bas) + 1):((station+1)*ncol(t.bas))]
    ti.par.stat<-all.ti.par[((station*ncol(s.bas) + 1)-10):((station+1)*ncol(s.bas)-10)]
    ti.add<-t.bas %*% ti.par.stat
    TI<-matrix(ti.add, nrow = n.yrs, ncol = 365/5, byrow = T)	
    inter <- inter+TI
  }
  
  # trend seasonal interaction
  if("ts"%in%model.type){
    all.ts.par<-parameters[(locs[match("ts", model.type) - 1] + 1):locs[match("ts", model.type)]]
    s.bas.ext<-matrix(1, ncol = 1,  nrow = n.yrs) %x% s.bas
    ts.bas<-box3.prod(as.spam(t.bas), as.spam(s.bas.ext))
    ts.add<-ts.bas %*% all.ts.par
    TS<-matrix(ts.add, nrow = n.yrs, ncol = 365/5, byrow = T)
    inter <- inter+TS
  }
  
  return_obj <- smnet_obj
  return_obj$main <- main 
  return_obj$inter <- inter
  
  if(plot.fig){
    obs.ind <-which(response.locs == station)
    dates.full <-dates[obs.ind]
    dates.stat <-decimal_date(as.Date(dates.full))
    vals.stat <- smnet_obj$response[obs.ind]
    
    if(!("t"%in%model.type) & !("ti"%in%model.type) & !("ts"%in%model.type) & ("s"%in%model.type) ){
      par(mar = c(4,4,2,0.5))
      dates.stat2 <- dates.stat - floor(dates.stat)
      plot(dates.stat2, vals.stat, pch = ".", cex = 1, type = "n",ylim = range(vals.stat),
           xlab = "Time", ylab = "log(Nitrate concentrations)",cex.lab = 1, cex.axis=1, main=paste("Segment ", station, '(W/O TREND COMPONENT)'))
      day.seq2 <- seq(0+5/365, 0+1, length.out=365/5)
      lines(day.seq2, c(t(main))[c(1:length(day.seq2))], type = "l", lwd=2) #col = "limegreen", 
      lines(day.seq2, c(t(inter))[c(1:length(day.seq2))]+ c(t(main))[c(1:length(day.seq2))], type = "l", lwd=2, lty=2) #col = "tomato", 
      points(dates.stat2, vals.stat)
      if(legend){
        legend("bottomleft", lty=c(1,2), c("Main", "With Interaction"), lwd=c(2,2)) # col=c("limegreen", "tomato"),
      }
    }else{
      par(mar = c(4,4,2,0.5))
      plot(dates.stat, vals.stat, pch = ".", cex = 1, type = "n",ylim = range(vals.stat),
           xlab = "Time", ylab = "log(Nitrate concentrations)",cex.lab = 1, cex.axis=1, main=paste("Segment ", station))
      lines(day.seq, c(t(main)), type = "l",  lwd=2, lty=2) #col = "limegreen",
      lines(day.seq, c(t(inter))+ c(t(main)), type = "l",lwd=2, lty=1) #col = "tomato", 
      points(dates.stat, vals.stat)
      if(legend){
        legend("topright", lty=c(1,2), c("With Interaction","Without Interaction"), lwd=c(2,2)) #col=c("limegreen", "tomato"), 
      }
    }
  }
  
  return(return_obj)
  
  

}

##### Transform date to decimal yearday ##### 

ymd_to_decimal <- function(date){
  if(length(date) > 1){
    sapply(date, ymd_to_decimal)
  }else{
    if(year(date) %% 4 != 0){
      return(yday(date)/365)
    }else if(year(date) %% 400 != 0){
      return(yday(date)/366)
    }else{
      return(yday(date)/365)
    }
  }
}

##### Expectile curve forecast based on FPCA and VAR #####

ExpecSTnet_forecast <- function(station, tau = 0.5, n.ahead = 1, 
                                  penalties = NULL, 
                                  plot.fit = TRUE, 
                                  log.y = TRUE,
                                  var.lag = c('AIC', 'SC')[1], 
                                  numPCs = 3,
                                  smnet_fit = NULL){
  if(is.null(smnet_fit)){
    poo <- ExpecSTnet(realweights, adjacency, TweedData_train, TweedPredPoints,
                      log.y=log.y,
                      station = station,
                      plot.fig = FALSE, tau = tau, 
                      model.type = c("c", "m", "s", "si", "t", "ti", "ts"), 
                      penalties = penalties)
  }else{
    poo <- ExpecSTnet_plot.estimated(station = station, smnet_obj = smnet_fit, 
                                     plot.fig = FALSE) 
  }
  dates <- as.character(poo$TweedData$date)
  parameters<-as.vector(poo$beta_hat)
  dims <- c(1, 113, 10, 1130, 10, 1130, 100)
  locs <- cumsum(dims) #locs : 1 114 124 1254 1264 2394 2494 (dims : 1 113 10 1130 10 1130 100)
  comps<-vector("list")	
  model.type = c("c", "m", "s", "si", "t", "ti", "ts")
  fit.point <- parameters[(locs[match("m", model.type)-1]+1):locs[match("m", model.type)]]
  beta_1 <- parameters[1]
  knts<-c(10,10) #s모형의 knot 개수, t모형의 knot 개수
  if(log.y){
    response <- log(poo$TweedData$nitrate)
  }else{
    response <- poo$TweedData$nitrate
  }
  response.locs <- poo$TweedData$location
  
  {
  # comps<-vector("list", length(model.type))
  min.yr<-min(year(dates))
  max.yr<-max(year(dates))
  yrs<-seq(min.yr, max.yr)
  n.yrs<-length(yrs)
  day.seq<-seq(min.yr + 5/365, max.yr + 1, length.out = (365/5)*n.yrs)
  }
  main <- poo$main
  inter <- poo$inter
  dates.stat.inyear <- seq(5,365,by=5)/365
  
  # plot
  if(plot.fit){
    par(mar = c(4,4,2,0.5))
    plot(dates.stat.inyear, dates.stat.inyear, pch = ".", cex = 1, type = "n", 
         xlim = c(0, 1), ylim = range(inter + main),
         xlab = "Time", ylab = "log(Nitrate concentrations)",cex.lab = 1, cex.axis=1, 
         main=paste(tau, "th expectile curve at Segment", station), 
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
    colfunc <- colorRampPalette(c("red","yellow","springgreen","royalblue"))
    col.vec <- colfunc(n.yrs)
    for(n in 1:n.yrs){
      # obs.ind <- which(response.locs == station & year(dates) == min.yr + n - 1)
      # dates.year.n <- dates[obs.ind]
      est.year.n <- (inter + main)[n,]
      lines(seq(5,365,by=5)/365, est.year.n, lwd = 2, col = col.vec[n])
    }
    legend("bottomleft", legend = seq(from = min.yr, by = 1, length = n.yrs), 
           lwd = 1, lty = 1, col = col.vec, cex = .5)
  }
  
  ##### FPCA #####
  expect_inyear <- matrix(nrow = 73, ncol = n.yrs)
  for(n in 1:n.yrs){
    expect_inyear[,n] <- (inter + main)[n,]
  }
  
  daybasis10 <- create.bspline.basis(c(0, 1), nbasis = 10, norder = 4)
  day5_expect_fd <- smooth.basis(seq(5,365,by=5)/365, expect_inyear, daybasis10)$fd
  
  day5_expect_pcaobj <- pca.fd(day5_expect_fd, nharm = numPCs)
  
  # par(mfrow = c(1,3))
  # plot.pca.fd(day5_expect_pcaobj)
  
  # dim(day5_expect_pcaobj$scores)
  
  day5_expect_pca_scores <- day5_expect_pcaobj$scores
  day5_expect_pca_scores <- matrix(day5_expect_pca_scores, ncol = numPCs)
  m <- dim(day5_expect_pca_scores)[2]
  colnames(day5_expect_pca_scores) <- paste0("s", 1:m)
  
  ##### VAR forecast #####
  
  var.select <- VARselect(day5_expect_pca_scores, lag.max=5, type="const")
  
  
  
  if(var.lag == 'AIC'){
    criteria <- var.select[["criteria"]]["AIC(n)",]
    var.lag.used <- which.min(criteria[is.finite(criteria)])
  }else if(var.lag == 'SC'){
    criteria <- var.select[["criteria"]]["SC(n)",]
    var.lag.used <- which.min(criteria[is.finite(criteria)])
  }else if(class(var.lag) == 'numeric'){
    var.lag.used <- var.lag
  }
  var <- vars::VAR(day5_expect_pca_scores, p = var.lag.used, type = "const")
  var_forecast <- predict(var, n.ahead = n.ahead, ci = .95)
  
  return(list(data.fd = day5_expect_fd, 
              pca.fd  = day5_expect_pcaobj, 
              var.forecast = var_forecast, 
              var.lag.used = list(var.lag, var.lag.used),
              n.yrs = n.yrs,
              inter = inter, 
              main = main, 
              fit = inter + main))
  
}

##### Expectile forecast #####

expectile_forecast <- function(station = 87, tau = .5, predict_date, 
                               n.ahead = 1, var.lag = c("AIC", "SC")[1],
                               forecast.obj = NULL, 
                               smnet_fit = NULL){
  # only for n.ahead = 1 now 
  if(is.null(forecast.obj)){
    forecast.obj <- ExpecSTnet_forecast(station = station, tau = tau, 
                                        n.ahead = n.ahead,
                                        var.lag = var.lag, 
                                        plot.fit = FALSE, 
                                        smnet_fit = smnet_fit)
  }
  poo.pca <- forecast.obj$pca.fd
  var2 <- forecast.obj$var.forecast
  forecast_coef <- poo.pca$meanfd$coefs + poo.pca$harmonics$coefs %*% (unlist(var2$fcst) %>% matrix(ncol = 3))[n.ahead,] 
  
  predict_decimaldate <- ymd_to_decimal(predict_date)
  
  fd(coef = forecast_coef, basisobj = forecast.obj$data.fd$basis) %>%
    eval.fd(evalarg = predict_decimaldate)
}

##### CDF estimation from expectiles #####

cdf_from_expect <- function (e, tau, 
                             e0 = NA, eR = NA, lambda = 0, var.dat = NA) 
{
  # e : expectile estimates at tau 
  # tau : vector of expectile levels 
  
  epsilon = 1e-05
  max.iter = 20
  step.halfing = 0.5
  p = tau
  e = sort(e)
  K <- length(e)
  if (is.na(var.dat) || var.dat < 0) {
    var.dat <- var(e)
  }
  if (is.null(p) == TRUE) {
    p <- seq(0 + 1/(K + 1), 1 - 1/(K + 1), length = K)
  }
  mat <- apply(matrix(p), 1, all.equal, 0.5)
  k0 <- which(mat == "TRUE")
  k0.i <- length(k0)
  if (is.na(e0) == TRUE) {
    e0 <- min(e) + (min(e) - min(e[-which.min(e)]))
  }
  if (is.na(eR) == TRUE) {
    eR <- max(e) + (max(e) - max(e[-which.max(e)]))
  }
  if (k0.i) 
    mu05 <- e[k0]
  else mu05 <- approx(p, y = e, xout = 0.5)$y
  eg <- c(e0, e)
  step <- eg[-1] - eg[-length(eg)]
  if (any(as.vector(step == 0))) {
    ind <- which(as.vector(step) == 0)
    step[ind] <- 1e-16
  }
  eg <- eg[-1] + eg[-length(eg)]
  eg <- eg/2
  egR <- (e[K] + eR)/2
  P <- diag(1/step[-1], K - 1, K - 1)
  P <- cbind(0, P)
  diag(P) <- -1/step[-K]
  Kmat <- t(P) %*% P
  delta <- rep(1/(K + 1), K)
  loop <- 0
  iter.diff <- 1
  while ((loop < max.iter) & (iter.diff > epsilon)) {
    loop <- loop + 1
    F <- cumsum(delta)
    Fs <- (kronecker(matrix(1:K), matrix(1, 1, K)) >= kronecker(matrix(1, 
                                                                       K, 1), t(matrix(1:K)))) * 1
    G <- cumsum(eg * delta)
    Gs <- kronecker(matrix(1, K, 1), t(matrix(eg))) * Fs
    h <- e - ((1 - p) * G + p * (mu05 - G))/((1 - p) * F + 
                                               p * (1 - F))
    if (k0.i) {
      h.tilde <- h[-k0]
      hk0 <- mu05 - (G[K] + egR * (1 - F[K]))
      h <- c(h.tilde, hk0)
    }
    hs <- -kronecker(matrix((1 - 2 * p)/((1 - p) * F + p * 
                                           (1 - F))), matrix(1, 1, K)) * Gs + kronecker(matrix((((1 - 
                                                                                                    p) * G + p * (mu05 - G)) * (1 - 2 * p))/(((1 - p) * 
                                                                                                                                                F + p * (1 - F))^2)), matrix(1, 1, K)) * Fs
    if (k0.i) {
      hs.tilde <- hs[-k0, ]
      hsk0 <- -eg + egR
      hs <- rbind(hs.tilde, hsk0)
    }
    Ls <- t(hs) %*% h
    Lss1 <- t(hs) %*% hs
    Lss <- Lss1
    dvec <- -Ls
    Dmat <- Lss + (lambda * var.dat^2) * Kmat
    if (qr(Dmat)$rank != K) {
      penalty.term <- 0.001 * var.dat
      Dmat <- Dmat + penalty.term * diag(K)
    }
    Amat <- cbind(diag(K), -matrix(1, K, 1))
    bvec <- matrix(c(-delta, sum(delta) - 1))
    xsi <- solve.QP(Dmat, dvec, Amat, bvec, meq = 0)$solution
    delta <- delta + step.halfing * xsi
    iter.diff <- max(abs(xsi))
  }
  Fdistcum <- cumsum(delta)
  
  result = list(x = e, cdf = Fdistcum)
  return(result)
}

expectile_levels <- c(.01, .05, 
                      seq(from = .1, to = .9, by = .1), 
                      .95, .99)


cdf_forecast <- function(station, predict_date, 
                         expectile_levels = c(.01, .05, c(1:9)/10, .95, .99),
                         xk = seq(from = -0.4, to = 2.6, length = 601), 
                         var.lag = c("AIC", "SC")[1],
                         n.ahead = 1,
                         numPCs = 3){

  
  expectile_estimates <- vector(length = length(expectile_levels))
  
  for(i in 1:length(expectile_levels)){
    tau <- expectile_levels[i]
    smnet_fit_tau <- ExpecSTnet(realweights, adjacency, TweedData_train, TweedPredPoints, 
                                log.y=TRUE,
                                plot.fig = FALSE, tau = tau, 
                                model.type = c("c", "m", "s", "si", "t", "ti", "ts"), 
                                penalties = NULL, 
                                station = 87)
    

    expectile_estimates[i] <- expectile_forecast(station = station, 
                                                 predict_date = predict_date, 
                                                 tau = tau,
                                                 n.ahead = n.ahead,
                                                 smnet_fit = smnet_fit_tau)
  }
  
  cdf_at_e <- cdf_from_expect(e = expectile_estimates, tau = expectile_levels)
  
  return_df <- data.frame(xk = xk)
  cdf.estimate = approx(x = cdf_at_e$x, y = cdf_at_e$cdf, xout = xk, 
                        yleft = 0, yright = 1)
  if(is.na(cdf.estimate$y[1])){
    cdf.estimate$y[1:(min(which(!is.na(cdf.estimate$y)))-1)] <- 0 
  }
  if(is.na(cdf.estimate$y[length(cdf.estimate$y)])){
    cdf.estimate$y[(max(which(!is.na(cdf.estimate$y)))+1):length(xk)] <- 1
  }
  return_df$cdf.estimate <- cdf.estimate$y
  
  return(return_df)
}

##### CRPS score #####

CRPS <- function(cdf_est, obs){
  xk <- cdf_est$xk
  Fk <- cdf_est$cdf.estimate
  return(mean((Fk - (obs <= xk))^2))
}

