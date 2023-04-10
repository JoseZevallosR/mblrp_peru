':=' <- function(lhs, rhs) {
  frame <- parent.frame()
  lhs <- as.list(substitute(lhs))
  if (length(lhs) > 1)
    lhs <- lhs[-1]
  if (length(lhs) == 1) {
    do.call(`=`, list(lhs[[1]], rhs), envir=frame)
    return(invisible(NULL))
  }
  if (is.function(rhs) || is(rhs, 'formula'))
    rhs <- list(rhs)
  if (length(lhs) > length(rhs))
    rhs <- c(rhs, rep(list(NULL), length(lhs) - length(rhs)))
  for (i in 1:length(lhs))
    do.call(`=`, list(lhs[[i]], rhs[[i]]), envir=frame)
  return(invisible(NULL))
}


points_wgs84=function(data){
  sp::coordinates(data) <- ~x+y
  sp::proj4string(data)='+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
  data
}

nonzero=function(obs){
  "cumulative distribution of none zero values"
  idx=obs>0
  #ecdf(obs[idx])
  obs[idx]
}


coords_correction=function(mapa){
  #rotation of trmm to correct positioning
  flip(flip(t(mapa), 1),2)
}

##################################
####Functions for filtering data #
####before regionalization########
##################################

kickOutliers=function(data){
  #delete neighbors based on variance outliers
  data_help=data
  sp::coordinates(data_help) <- ~x+y
  sp::proj4string(data_help)='+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
  mdist <- geosphere::distm(data_help,fun = geosphere::distHaversine)
  neighbors=nearpoints(mdist)
  outlier=c()
  for (station in 1:dim(data)[1]){
    sub=data[c(station,neighbors[[station]]),]
    outvals=which(sub$var24 %in% boxplot(sub$var24)$out)
    outlier=c(outlier,outvals)
  }

  idx=unique(outlier)
  data[-idx,]
}

filter_Neigbors=function(data,min_n=3,radio=60000){
  #Return the statins with at least n neighbors
  #data= location of the stations
  #min_n= number of minumun neighbors required
  n_old=dim(data)[1]
  n_new=n_old*2

  while (n_old-n_new!=0){
    n_old=dim(data)[1]
    data_help=data
    sp::coordinates(data_help) <- ~x+y
    sp::proj4string(data_help)='+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
    mdist <- geosphere::distm(data_help,fun = geosphere::distHaversine)
    neighbors=nearpoints(mdist,radio=radio)
    data=data[lengths(neighbors)>=min_n,]
    n_new=dim(data)[1]
  }
  data
}

clusterIDX=function(data){
  #return cluster index
  #data contains the intial parameter estimation
  #stats is the rainfall statistics
  drop_range=function(x,vector){
    for (element in vector){
      x=x[x!=element]
    }
    x
  }
  #Estaciones cambiantes de intervalos
  range=1:dim(data)[1]


  data_help=data
  sp::coordinates(data_help) <- ~x+y
  sp::proj4string(data_help)='+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'

  mdist <- geosphere::distm(data_help,fun = geosphere::distHaversine)
  vecinos=nearpoints(mdist)

  idx=c()
  while (length(range)!=0){

    station=range[1]
    range=drop_range(range,c(station,vecinos[[station]]))
    idx=c(idx,station)
  }#ends parameter iteration
  idx
}

