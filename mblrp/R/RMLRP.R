####################################
###Model Statistics#################
####################################
#Mean
meanMBLRPM<-function(a,l,v,k,f,mx,h=1) {
  x<-(h*l*mx*v*(1+k/f))/(a-1)
  return(x)
}
#Variance
varMBLRPM<-function(a,l,v,k,f,mx,h=1) {
  A<-(2*l*(1+k/f)*(mx^2)*(v^a))/((f^2)*((f^2)-1)*(a-1)*(a-2)*(a-3))
  B<-(2*(f^2)-2+k*f)*(f^2)*((a-3)*h*(v^(2-a))-(v^(3-a))+((v+h)^(3-a)))
  C<-k*(f*(a-3)*h*(v^(2-a))-(v^(3-a))+((v+f*h)^(3-a)))
  D<-A*(B-C)
  return(D)
}
#Covariance
covarMBLRPM<-function(a,l,v,k,f,mx,h=1,lag=1) {
  A<-(l*(1+k/f)*(mx^2)*(v^a))/((f^2)*((f^2)-1)*(a-1)*(a-2)*(a-3))
  B<-(2*(f^2)-2+k*f)*(f^2)*(((v+(lag+1)*h)^(3-a))-2*((v+lag*h)^(3-a))+((v+(lag-1)*h)^(3-a)))
  C<-k*(((v+(lag+1)*h*f)^(3-a))-(2*((v+h*lag*f)^(3-a)))+((v+(lag-1)*h*f)^(3-a)))
  D<-A*(B-C)
  return(D)
}
#Dry probabilities
#pdrRPBLRPM<-function(a,l,v,k,f,h=1) {
#    mt<-((1+(f*(k+f))-(0.25*f*(k+f)*(k+4*f))+((f/72)*(k+f)*(4*(k^2)+27*k*f+72*(f^2))))*v)/(f*(a-1))
#    G00<-((1-k-f+1.5*k*f+(f^2)+0.5*(k^2))*v)/(f*(a-1))
#    A<-(f+(k*(v/(v+(k+f)*h))^(a-1)))/(f+k)
#    D<-exp(l*(-h-mt+G00*A))
#    return(D)
#}

####################################
#######Optimization Function########
####################################

MBLRPM=function(mean24,var24,cov24lag1,pdr24,var3,cov3lag1,var6,var12,var18,Lmin,Lmax){
  ####################################
  ###Model Statistics#################
  ####################################
  #Mean
  meanMBLRPM<-function(a,l,v,k,f,mx,h=1) {
    x<-(h*l*mx*v*(1+k/f))/(a-1)
    return(x)
  }
  #Variance
  varMBLRPM<-function(a,l,v,k,f,mx,h=1) {
    A<-(2*l*(1+k/f)*(mx^2)*(v^a))/((f^2)*((f^2)-1)*(a-1)*(a-2)*(a-3))
    B<-(2*(f^2)-2+k*f)*(f^2)*((a-3)*h*(v^(2-a))-(v^(3-a))+((v+h)^(3-a)))
    C<-k*(f*(a-3)*h*(v^(2-a))-(v^(3-a))+((v+f*h)^(3-a)))
    D<-A*(B-C)
    return(D)
  }
  #Covariance
  covarMBLRPM<-function(a,l,v,k,f,mx,h=1,lag=1) {
    A<-(l*(1+k/f)*(mx^2)*(v^a))/((f^2)*((f^2)-1)*(a-1)*(a-2)*(a-3))
    B<-(2*(f^2)-2+k*f)*(f^2)*(((v+(lag+1)*h)^(3-a))-2*((v+lag*h)^(3-a))+((v+(lag-1)*h)^(3-a)))
    C<-k*(((v+(lag+1)*h*f)^(3-a))-(2*((v+h*lag*f)^(3-a)))+((v+(lag-1)*h*f)^(3-a)))
    D<-A*(B-C)
    return(D)
  }
  #Dry probabilities
  #pdrRPBLRPM<-function(a,l,v,k,f,h=1) {
  #  mt<-((1+(f*(k+f))-(0.25*f*(k+f)*(k+4*f))+((f/72)*(k+f)*(4*(k^2)+27*k*f+72*(f^2))))*v)/(f*(a-1))
  #  G00<-((1-k-f+1.5*k*f+(f^2)+0.5*(k^2))*v)/(f*(a-1))
  #  A<-(f+(k*(v/(v+(k+f)*h))^(a-1)))/(f+k)
  #  D<-exp(l*(-h-mt+G00*A))
  #  return(D)
  #}
  symvar<- function(a,l,v,k,f,mx,h,var){
    (1-varMBLRPM(a,l,v,k,f,mx,h)/var)^(2)+(1-var/varMBLRPM(a,l,v,k,f,mx,h))^(2)
  }
  symcovar <- function(a,l,v,k,f,mx,h,cov){
    (1-covarMBLRPM(a,l,v,k,f,mx,h)/cov)^(2)+(1-cov/covarMBLRPM(a,l,v,k,f,mx,h))^(2)
  }
  symmean <- function(a,l,v,k,f,mx,h,meann){
    (1-meanMBLRPM(a,l,v,k,f,mx,h)/meann)^(2)+(1-meann/meanMBLRPM(a,l,v,k,f,mx,h))^(2)
  }
  sympdr <- function(a,l,v,k,f,h,pdrr) {
    (1-pdrRPBLRPM(a,l,v,k,f,h)/pdrr)^(2)+(1-pdrr/pdrRPBLRPM(a,l,v,k,f,h))^(2)
  }
  #Objective function
  fopt <- function(x) {
    a<-x[1];l<-x[2];v<-x[3];k<-x[4];f<-x[5];mx<-x[6]


    w1=1;w2=1;w3=1;w4=1;w5=1;w6=1;

    S3<-w2*symvar(a,l,v,k,f,mx,h=3,var3)+w3*symcovar(a,l,v,k,f,mx,h=3,cov3lag1)

    S6<-w2*symvar(a,l,v,k,f,mx,h=6,var6)

    S12<-w2*symvar(a,l,v,k,f,mx,h=12,var12)

    S18<-w2*symvar(a,l,v,k,f,mx,h=18,var18)

    w1=5;w2=5;w3=5;w4=5

    S24 <- w1*symmean(a,l,v,k,f,mx,h=24,mean24)+ w2*symvar(a,l,v,k,f,mx,h=24,var24)+ w3*symcovar(a,l,v,k,f,mx,h=24,cov24lag1)+w4*sympdr(a,l,v,k,f,h=24,pdr24)

    S<-S24+S3+S6+S12+S18



    if(is.infinite(S)) {S<-10^8}
    if(is.na(S)) {S<-10^8}
    return(S)
  }

  # set the interior and exterior parameters bounds
  xmin <- Lmin
  xmax <- Lmax
  xlow <- Lmin
  xup <- c(50,20,runif(1,min = 1.8, max = 2.2),runif(1,min = 0.05, max = 0.09),runif(1,min = 1.8, max = 2.2),runif(1,min = 8, max = 12))

  modecal <- eas(n=6,m=30,xmin,xmax,xlow,xup,fn=fopt,maxeval=5000,ftol=1.e-10,ratio=0.99,pmut=0.95, beta=2,maxclimbs=5)
  modecal
  a<-modecal$bestpar[[1]];
  l<-modecal$bestpar[[2]];
  v<-modecal$bestpar[[3]];
  k<-modecal$bestpar[[4]];
  f<-modecal$bestpar[[5]];
  mx<-modecal$bestpar[[6]]
  # In order to use the derived parameters in the functions of HyetosR
  # as well as in the classic version of Hyetos,please be sure that
  # for parameters mx and sx the length units are millimeters (mm)
  # and for parameters l, v, mx and sx the time units are days (d).
  # For this reason, make the following unit conversions:

  #checar
  #l<-l*24
  #v<-v/24
  #mx<-mx*24

  # parameter set for implementation in HyetosR functions
  par <- c(a=a,l=l,v=v,k=k,f=f,mx=mx)

  par
}


nearpoints=function(mdist,radio=60000){
  #Return the index of the closest gauge station to each station
  #mdist = matrix of distance in meters
  #radio is the maximum distance allowed among station
  distance=list()
  for (i in 1:dim(mdist)[1]){
    values=sort(mdist[i,])
    values=tail(values, -1)

    condition=values[values<radio]
    near=numeric(length(condition))

    for (j in 1:length(condition)){
      near[j]=which(mdist[i,]==condition[j])[1]
    }

    if (length(near)>1){
      distance[[i]]=near
    }else{
      if (is.na(near)){
        distance[[i]]=0
      }else{
        distance[[i]]=near
      }
    }

  }
  distance
}

################################################
#Cross Validations using Inverse Distance Weigth
################################################
idwCV=function(data,parameter='a',power=2){
  x=data
  coordinates(x) <- ~x+y
  proj4string(x)='+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
  mdist <- distm(x) #the answer is in meters

  crossValidated=numeric()
  for (i in 1:dim(data)[1]){
    info=data[[parameter]][-i]
    denominador=sum((1/mdist[i,-i])^power)
    crossValidated[i]=sum(info/mdist[i,-i]^power/denominador)
  }

  rr=cbind(data[c('x','y')],'var1.pred'=crossValidated,'observed'=data[[parameter]])
  rr$residual=rr$observed-rr[['var1.pred']]
  rr
}


####################################
#######Repetitive Cross Validation##
####################################

repetitiveCV=function(times=1,data,Stats,Lmin,Lmax,fun=MBLRPM){
  #data contains the intial parameter estimation
  #stats is the rainfall statistics

  #Estaciones cambiantes de intervalos

  parameters=data[,3:8]
  data_help=data
  coordinates(data_help) <- ~x+y
  proj4string(data_help)='+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'

  mdist <- distm(data_help,fun = distHaversine)
  vecinos=nearpoints(mdist)

  range=clusterIDX(data)
  old_error=rep(100,length(range))
  current_error=numeric(length(range))
  for (iter in 1:times){
    print(paste("Number of cross validation iteration",as.character(iter)))
    k_cluster=1
    for (station in range){
      print(paste('cluster :',as.character(k_cluster),'/',as.character(length(range))))
      mistakes=c()
      for (k in 1:6){#parameters, 6 in total
        x <- idwCV(data[c(station,vecinos[[station]]),],parameter= c('a','l','v','k','f','mx')[k],power=2)
        #Checking region error
        sub=x[c('x','y','var1.pred','observed','residual')]
        sub$porcentaje=abs(sub$residual)*100/sub$observed
        sub$residual=NULL
        sub$location=c(station,vecinos[[station]])

        #median(sub$porcentaje)
        good_neighbors=subset(sub,sub$porcentaje<30)
        wrong_neigbors=subset(sub,sub$porcentaje>=30)

        if (dim(wrong_neigbors)[1]!=0 & dim(good_neighbors)[1]!=0 ){
          for (fix_id in wrong_neigbors$location){
            Lmin[k,fix_id]=min(good_neighbors$var1.pred)
            Lmax[k,fix_id]=max(good_neighbors$var1.pred)

            mistakes=c(mistakes,fix_id) #add the wrong stations
          }
        }
      }
      n_parameters=length(mistakes)
      mistakes=unique(mistakes)

      if (length(mistakes)==0){
        current_error[k_cluster]=0
        print('No parameters to correct')
      }else{
        current_error[k_cluster]=round(n_parameters*100/(6*length(sub$location)),2)
        if(old_error[k_cluster]<current_error[k_cluster]){
          print('Using old parameter with lower error')
          print(paste("Incorrect parameters in the cluster:",as.character(n_parameters),'/',as.character(6*length(sub$location)),' (',as.character(old_error[k_cluster]),'%)'))
          parameters[mistakes,]=old_parameters[mistakes,]
          old_error[k_cluster]=old_error[k_cluster]
        }else{

          print(paste("Incorrect parameters in the cluster:",as.character(n_parameters),'/',as.character(6*length(sub$location)),' (',as.character(current_error[k_cluster]),'%)'))
          n.cores <- parallel::detectCores() - 1
          #create the cluster
          my.cluster <- parallel::makeCluster(
            n.cores,
            type = "PSOCK"
          )
          #register it to be used by %dopar%
          doParallel::registerDoParallel(cl = my.cluster)

          parameters[mistakes,]=t(matrix(foreach(
            i=mistakes,
            .combine = 'c',
            .packages = "HyetosMinute"
          ) %dopar% {
            momentos=Stats[i,]

            mean24 = momentos$mean24
            var24 = momentos$var24
            cov24lag1 =momentos$autocov24
            pdr24=momentos$dryperiod24
            var3=momentos$var3
            cov3lag1=momentos$autocov3
            var6=momentos$var6
            var12=momentos$var12
            var18=momentos$var18

            par=fun(mean24,var24,cov24lag1,pdr24,var3,cov3lag1,var6,var12,var18,Lmin[,i],Lmax[,i])

            return(par)
          },nrow = 6,ncol = length(mistakes)))
          parallel::stopCluster(cl = my.cluster) #closing the cluster


          old_error[k_cluster]=current_error[k_cluster]
        }
      }


      k_cluster=k_cluster+1 #counting clusters

    }
    #parameters=cbind(Stats[,1:2],parameters)
    names(parameters)=c('a','l','v','k','f','mx')#c('x','y','a','l','v','k','f','mx')
    old_parameters=data[,3:8]

    data[,3:8]=parameters#check

    #write.table(data,paste0("D:/Proyectos_GitHub/Bartlet-Lewis_Regionalization/output/CV_parameters/iteraciones/",'parameters_iter_',as.character(iter),'.csv'),sep = ',',row.names = F)
  }

  data
}

##########################################
###### Run Function#######################
##########################################
run=function(rain_stats,path,iterations=5,fun=MBLRPM,FILE_NAME){
  #rain_stats: contains the rainfall statistics
  #path: where to save the results
  #Maskshape: Shape form of the final results


  n=dim(rain_stats)[1]

  #Maximum and minimum search parameters space
  Lmin=matrix(c(0.1,0.001,0.001,0.001,0.0854,1),nrow = 6,ncol = n)
  Lmax=matrix(c(4,0.1,0.1,0.1,0.1,20),nrow=6,ncol=n)


  #print('Calculating the initial parameters ...')
  #Initial parameters
  parameters0=matrix(data=NA,nrow =n,ncol = 6)

  n.cores <- parallel::detectCores() - 1
  #create the cluster
  my.cluster <- parallel::makeCluster(
    n.cores,
    type = "PSOCK"
  )
  #register it to be used by %dopar%
  doParallel::registerDoParallel(cl = my.cluster)

  parameters0=t(matrix(foreach(
    i=1:n,
    .combine = 'c',
    .packages = "HyetosMinute"
  ) %dopar% {

    momentos=rain_stats[i,]

    mean24 = momentos$mean24
    var24 = momentos$var24
    cov24lag1 =momentos$autocov24
    pdr24=momentos$dryperiod24
    var3=momentos$var3
    cov3lag1=momentos$autocov3
    var6=momentos$var6
    var12=momentos$var12
    var18=momentos$var18
    par=fun(mean24,var24,cov24lag1,pdr24,var3,cov3lag1,var6,var12,var18,Lmin[,i],Lmax[,i])

    return(par)


  },nrow = 6,ncol = n))

  parallel::stopCluster(cl = my.cluster) #closing the cluster

  parameters=cbind(rain_stats[,1:2],parameters0)
  names(parameters)=c('x','y','a','l','v','k','f','mx')

  print('Reptitive Cross Validations ...')
  CV_parameters=repetitiveCV(times = iterations,parameters,rain_stats,Lmin = Lmin ,Lmax = Lmax)
  #saving the initial parameters
  write.table(CV_parameters,paste0(path,FILE_NAME),sep = ',',row.names = F)
  CV_parameters
}


##################
SimStats= function(parameters){

  stats=matrix(data=NA,nrow = dim(parameters)[1],ncol = 16)
  for (i in 1:dim(parameters)[1]){
    par=parameters[i,]
    par[2]<-par[2]*24
    par[3]<-par[3]/24
    par[6]=par[6]*24

    est=numeric(16)
    iter=0
    for (j in c(24,3,6,12)){
      m=meanMBLRPM(par[1],par[2],par[3],par[4],par[5],par[6],h=j/24)
      v=varMBLRPM(par[1],par[2],par[3],par[4],par[5],par[6],h=j/24)
      cov=covarMBLRPM(par[1],par[2],par[3],par[4],par[5],par[6],h=j/24)
      pdr=pdrRPBLRPM(par[1],par[2],par[3],par[4],par[5],h=j/24)
      est[(1+iter*4):(4+iter*4)]=c(m,v,cov,pdr)
      iter=iter+1
    }

    stats[i,]=as.numeric(est)
  }
  stats
}

precp_sim=function(par,n,tscale=24){
  #par are the MBLR parameters
  #n is the number of days

  a=par[1]
  l=par[2]*24
  v=par[3]/24
  k=par[4]
  f=par[5]
  mx=par[6]*24

  tt=SequentialSimul(Length=n,BLpar=list(lambda=l,phi=f,kappa=k
                                         ,alpha=a,v=v,mx=mx,sxmx=1),CellIntensityProp=list(Weibull=FALSE,
                                                                                           iota=NA),TimeScale=tscale,ExportSynthData=list(exp=TRUE,FileContent=c("AllDays"),DaysPerSeason=31,
                                                                                                                                          file="SynthRPBLM.txt"),ImportHistData=list(imp=F,file="HistHourlyData.txt",
                                                                                                                                                                                     ImpDataTimeScale=1,na.values="NA",FileContent=c("AllDays"),DaysPerSeason=31,DailyValues=TRUE),
                     PlotTs=FALSE,Statistics=list(print=TRUE,plot=FALSE),RandSeed=NULL )[[1]]
  as.numeric(tt)
}

