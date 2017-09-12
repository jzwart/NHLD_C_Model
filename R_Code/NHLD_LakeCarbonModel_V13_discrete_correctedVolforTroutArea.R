##### lake carbon model for NHLD from VIC / GFLOW input 
# JAZ; 2016-02-10 using code originally from SEJ
# V13 has modeled zmix, kd, and a 2-box lake model and wetland loading and depth integrated GPP and light attenuated by snow depth and updated water temp budget; 
#  2 box model of DOC (labile and recalcitrant) ; splitting stream water into equal proportions of epi and hypo ; wind sheltering coefficient for k600 
# splitting SWout into epi and hypo 
# V13 does better accounting of DIC pool fluxes in the epi (previous versions had problems with discrete solver and too much flux out + instability)
# V13 includes carbonate speciation based on alkalinity, water temperature, and DIC   
# V13 fixing unit issue on gas flux (was dividing by zmix in previous iterations)

# clear workspace
rm(list=ls())

########## load utility functions and packages
source('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/R Code/NHLD_LakeCarbonModel_supporting_V13_discrete.R')
source('/Users/Jake/Desktop/R functions/AveLightClimate.R')
source('/Users/Jake/Desktop/R functions/DOY.r')
require(deSolve)
require(LakeMetabolizer)
require(snow)
require(sp)
require(rgeos)
require(parallel)
require(rgdal)
require(maptools)
require(aspace)
require(AquaEnv)
require(marelac)

# *************************** NEED TO CHANGE DIRECTORIES ONCE RUNNING ON CRC ****************************
dir<-'D:/Jake/My Papers/NHLD Carbon Model/20170815/'

# flux directory (from VIC); don't need 
# fluxDir<-file.path(dir,)
# forcings directory (from VIC)
forceDir<-file.path(dir,'Force')
# daily lake hydrology flux directory
lakeFluxDir<-file.path(dir,'OUT_RESULTS_TroutWatershed_Modified')
# daily lake geomorphology directory 
lakeGeomorphDir<-file.path(dir,'OUT_RESULTS_TroutWatershed_Modified')
# daily lake initial conditions directory 
lakeInitGeomorphDir<-file.path(dir,'OUT_RESULTS_TroutWatershed_Modified')
# lake watershed information 
lakeShedDir<-'/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Data/C model forcing data/'
watersheds<-read.table(file.path(lakeShedDir,"NHLDsheds_20170323.txt"),header=TRUE,sep="\t",stringsAsFactors=FALSE)
watersheds$percentWetland<-watersheds$percentWetland*100
# lake location, area, and perimeter 
lakeLocation<-readShapeSpatial(file.path(dir,'LakeLocations/NHLDandBuffLakes_Original_ZJH_Rev1.shp'))
lakeLocation<-data.frame(lakeLocation)
lakeLocation$Permanent_<-as.character(lakeLocation$Permanent_)
    
# keeling curve for historic CO2 concentrations 
keeling<-read.table('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Data/keelingCurve.txt',sep=' ',header = F,stringsAsFactors = F)
colnames(keeling)<-c('site_code','year','month','day','hour','minute','second','co2','value_unc','nvalue','latitude','longitude','altitude','elevation','intake_height','instrument','qcflag')
keeling$datetime<-as.Date(paste(keeling$year,keeling$month,keeling$day),format='%Y %m %d')
keeling$co2<-ifelse(keeling$co2==-999.99,NA,keeling$co2)
keeling$co2[1]<-keeling$co2[min(which(!is.na(keeling$co2)))]
keeling$co2[length(keeling$site_code)]<-keeling$co2[max(which(!is.na(keeling$co2)))]
keeling$co2<-approx(keeling$datetime,keeling$co2,keeling$datetime)$y #filling in NA's / missing data 
keeling<-keeling[,c('datetime','co2')]

# snow depth on top of ice; which is the same for all lakes 
snowDepth=read.table(file.path(dir,'Lake_Snow_Depth_Land_SWE_Depth_Ratio_Correction.txt'),sep = '\t')
snowDepth$datetime<-as.Date(paste(snowDepth$V1,snowDepth$V2,snowDepth$V3),format = '%Y %m %d')
snowDepth=snowDepth[,c('datetime','V4')]
colnames(snowDepth)[2]<-'snowDepth_m' # m; snow depth on ice 

# coordinate lookup table for forcing data 
forceCoordList<-strsplit(list.files(forceDir),'_')
if(length(grep('format.txt',forceCoordList,fixed=T))>0){
  forceCoordList<-forceCoordList[-grep('format.txt',forceCoordList,fixed=T)]
}
forceCoord<-data.frame()
for(i in 1:length(forceCoordList)){
  cur<-data.frame(Longitude=forceCoordList[[i]][3],Latitude=forceCoordList[[i]][2])
  forceCoord<-rbind(forceCoord,cur)
}
forceCoord$Longitude<-as.numeric(as.character(forceCoord$Longitude))
forceCoord$Latitude<-as.numeric(as.character(forceCoord$Latitude))
sp.forceCoord<-forceCoord
coordinates(sp.forceCoord)<- ~Longitude+Latitude

# starting year/month/day, ending year/mnth/day, & set up force/flux
startYear=1980
startMonth=5
startDay=1

endYear=2010
endMonth=12
endDay=31

# lakes to loop through 
lakes<-list.files(lakeFluxDir)
lakes<-lakes[grep('_DailyLakeFlux.txt',lakes,fixed = T)]
lakes<-as.character(unlist(strsplit(lakes,split = '_DailyLakeFlux.txt')))
lakes<-lakes[lakes%in%as.character(lakeLocation$Permanent_)] # only keeping primary near field lakes for given run  

# constants for biogeochem model 
kH=29.41  # Henry's Law constant for CO2 in water [L atm mol-1]  # temperature sensitivity
b=105 # half saturation constant for GPP model 
umax=0.9 # max growth rate of phytoplantkon for GPP model 
kQP=0.003 # Qmin - minimum P quota for phytoplankton in GPP model 
m=0.05 # mortality of phytoplankton; day-1 
vj=0.04 # uptake rate of phosphorus for phytoplankton; ug P ug C-1 day-1 
hp=12 # half saturation constant for phosphorus
phytoRecycling=0.2 # 1 - recyling rate of phytoplankton; i.e. 0.05 means 95% recycling rate ; Scavia 1979 P recycling in lake Ontraio = 86-93% 
phytoC2P=95#106		#M:M; from Patrick
lossPhyto=0.1 # loss rate of phytoplankton [day-1]
GPPexudeR=0.03  # Hanson et al. 2004; DOC exude from phytoplankton that is raclcitrant 
GPPexudeL=0.07 # DOC exude from phytoplankton that is labile 
Rauto=0.8 	# Hanson et al. 2004; quick respiration of phytoplankton production 
phytoDeath=0.03 # Hanson et al 2004; fraction of phyto that die; day-1
phyto_C2Chl=50 # phytoplankton carbon to chlorophyll ratio, [gC gChl-1]; sort of made up/average of observations in paper (e.g. 
zmix=2 # initial zmix for model 
kdSnow=10 #m-1; Perovich 2007 Journal of Glaciology 53:201-210; ranged from 10-30 m-1 for visible range 
rhoW=1000 #density of water kg m-3
cW=4186 #specific heat of water J kg-1 degC-1
cA=1012 # specific heat of air J kg-1 dec C-1 
#  coefficients from Sparkling in Duan and Bastiaanssen 2015 for Qt equation 
Qa=0.91
Qb=-101.00
Qc=4.49
Gsc=0.082 # MJ m-2 min-1; solar constnat 
conversion=11.6 # unit converter from MJ m-2 day-1 to W m-2 
aSlob=110 # empirical coefficient for Slob's equation; W m-2
#Factors to convert degrees to radians and vice versa
degToRad <- 2*pi/360
radToDeg <- 180/pi
alpha=0.05 # albedo for water 
rSlow=0.0015 # day-1; decay rate of recalcitrant DOC 
rFast=0.08 # day-1; decay rate of labile DOC; Berggren et al. 2010 
fracLabile=0.10 # fraction; fraction of loaded DOC that is labile Berggren et al. 2010
sedBurial_tPOC=0.9 # tPOC fraction that is permanently buried 
sedBurial_phyto=0.25 # phyto fraction that is permanently buried 
alpha_sw=0.07 # albedo for sw from Lenters 2005
alpha_lw=0.03 # albedo for lw from Lenters 2005 
kelvinZero=273.15 
wLoad=2/365 # wetland load per shoreline length; g C m-1 shoreline day-1 (mean from fitted parameter in Hanson et al. 2014) 
sal=0 # salinity set to 0 for all lakes 
leafLoad=300/12 #mol C m-1 shoreline yr-1; autumn leaf fall; check out Likens 1985 and Gasith and Hasler 1976; Hanson et al 2014; France 1996; France and Peters 1996; 300 g C m-1 shoreline yr-1 is roughly the average for all these studies 

length(lakes)
ntlLookUp<-data.frame(lakeID=c('AL','BM','CR','SP','TR','TB','CB'),Permanent_=c(69886156,69886284,69886510,69886444,69886228,69886158,123148117))
lake='SP'
i=grep(ntlLookUp$Permanent_[ntlLookUp$lakeID==lake],lakes) # Big Musky = 69886284, Sparkling = 69886444, Allequash = , Trout = , Crystal = , Crystal Bog = , Trout Bog = 
date<-read.table(file.path(dir,'Date_OUT_ALL.txt'),stringsAsFactors = F,header = F,sep='\t')
datetime<-as.Date(paste(date$V1,date$V2,date$V3),format = '%Y %m %d')

# failing lakes from 20170407 run 
# failLakes=list.files('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Results/C Model Output/20170407/')
# failLakes=unlist(strsplit(failLakes,'_C_model.txt',fixed = T))
# # failLakes=lakes[!lakes%in%failLakes]
# lakesRun=lakes[1:800]
# failLakes=lakesRun[!lakesRun%in%failLakes]
# 
# i=grep(failLakes[12],lakes)
# i=grep('JAZ_201',lakes) # a small lake that causes problems 

# for(i in 1:length(lakes)){
  print(i)
  # load in current lake hydrology output data 
  curLakeHydro<-read.table(file.path(lakeFluxDir,paste(lakes[i],'_DailyLakeFlux.txt',sep='')),header=F,stringsAsFactors = F)
  curLakeHydro<-cbind(datetime,curLakeHydro)
  colnames(curLakeHydro)<-c('datetime','SWin','DirectP','LakeE','GWin','GWout','SWout','IceSnow','LandMelt','IceMelt','SWoutMinusLandMelt','LandMeltEst','Baseflow')
  # load in current lake geomorphology output data 
  curLakeGeomorph<-read.table(file.path(lakeGeomorphDir,paste(lakes[i],'_DailyLakeProp.txt',sep='')),header=F,stringsAsFactors = F)
  curLakeGeomorph<-cbind(datetime,curLakeGeomorph)
  colnames(curLakeGeomorph)<-c('datetime','Area','Radius','Perim','Stage','Elev','Vol')
  curLakeHydro<-merge(curLakeHydro,curLakeGeomorph,by='datetime',all.x=T)
  rm(curLakeGeomorph)
  curLakeHydro<-merge(curLakeHydro,snowDepth,by='datetime',all.x=T)
  # load in current lake initial geopmorphological parameters 
  curLakeInitGeomorph<-read.table(file.path(lakeInitGeomorphDir,paste(lakes[i],'_IniLakeProp.txt',sep='')),header=F,stringsAsFactors = F)
  colnames(curLakeInitGeomorph)<-c('Area0','Vol0','Radius0','Diameter0','Perim0','Stage0','DL','r2h','WA','WALA','Elev0_DEM','Vol_LinRes','Stage_LinRes','Stream_WA')
  # watershed data 
  curWatershed<-watersheds[watersheds$Permanent_==lakes[i],]
  
  # location of lake used for finding correct forcing data 
  curLakeLoc<-lakeLocation[lakeLocation$Permanent_==lakes[i],]
  sp.curLakeLoc<-curLakeLoc
  coordinates(sp.curLakeLoc)<- ~Longitude+Latitude
  # find closest forcing data 
  curForceCoord<-forceCoord[which(gDistance(sp.forceCoord,sp.curLakeLoc,byid=T)==gDistance(sp.forceCoord,sp.curLakeLoc)),]
  # curForce<-read.table(file.path(forceDir,paste('FORCE',curForceCoord$Latitude,curForceCoord$Longitude,sep='_')),header=F,stringsAsFactors = F)
  curForce<-read.table(file.path(forceDir,paste('FORCE',curForceCoord$Latitude,curForceCoord$Longitude,sep='_')),header=F,stringsAsFactors = F)
  colnames(curForce)<-c('YYYY','MM','DD','HH','Precip','AirTemp','Longwave','Shortwave','Density','VPD','Pressure','Wind')
  curForce$datetime<-as.Date(paste(curForce$YYYY,curForce$MM,curForce$DD,curForce$HH),format='%Y %m %d %H')
  # keep forcing data for days over which you want to model 
  curForce<-curForce[curForce$datetime>as.Date(paste(startYear,startMonth,startDay),format='%Y %m %d')&
                       curForce$datetime<as.Date(paste(endYear,endMonth,endDay),format='%Y %m %d'),]
  
  #get forcing data for lake dates 
  curForce<-curForce[curForce$datetime%in%curLakeHydro$datetime,]
  curLakeHydro<-curLakeHydro[curLakeHydro$datetime%in%curForce$datetime,]
  curForce$I0<-sw.to.par(curForce,sw.col = 'Shortwave')$par # converting from shortwave radiation to PAR from LakeAnalyzer 
  curForce$I0<-ifelse(curForce$I0==0,0.001,curForce$I0)
  #merging keeling curve and hydrologic forcings 
  curLakeHydro<-merge(curLakeHydro,keeling,by='datetime',all.x=T)
  curLakeHydro$co2[1]<-curLakeHydro$co2[min(which(!is.na(curLakeHydro$co2)))]
  curLakeHydro$co2[length(curLakeHydro$datetime)]<-curLakeHydro$co2[max(which(!is.na(curLakeHydro$co2)))]
  curLakeHydro$co2<-approx(curLakeHydro$datetime,curLakeHydro$co2,curLakeHydro$datetime)$y
  
  nObs<-1:length(curLakeHydro$datetime)
  nObs_hourly<-rep(nObs,each=24)

  maxWind<-tapply(curForce$Wind,nObs_hourly,FUN=mean)
  #wind sheltering coefficient parameters for Read et al. 2014 equation 
  hs=7.5 # m; average elevation difference between top of surrounding canopy and lake surface; about average value from Read et al. 2014
  # wind sheltering coefficient from Read et al. 2014 
  if(curLakeInitGeomorph$Area0-625*hs^2*pi<0){
    Ws=0.001 # extremely small if small lake 
  }else{
    Ws=(2/pi)*acos(25*hs*sqrt((pi/curLakeInitGeomorph$Area0)))-(50*hs)/(curLakeInitGeomorph$Area0*sqrt(pi))*
      sqrt(curLakeInitGeomorph$Area0-625*hs^2*pi)
  }
  lakeWind=maxWind*Ws # adjusting lake wind by wind sheltering coefficient 
  
  # calculate concentrations of DIC, DOC, TP based on landcover
  # could use Frost et al 2005 for stream DOC; and Johnston et al 2008 ; Jutras etal 2011 
  # could use Frost et al 2009 for stream C:P of seston 
  streamDOC=exp(1.3961+3.245*(curWatershed$percentWetland/100))/12*1000/1000	# from lottig 2012; mol m-3
  streamPOC=3/12*1000/1000		# ~3 mg L-1; buffam 2011; mol m-3
  streamDIC=10.85/12*1000/1000		#10 mg L-1; lottig 2011; mol m-3
  streamP=0.05/31*1000/1000#0.04/31*1000/1000	  # Long inlet has 84 ug/L	#.025 mg L-1 TDP & 0.04 mg L-1 TP; lottig 2011; mol m-3
  gwDOC=13/12*1000/1000#median(ntlGW$doc,na.rm=TRUE)/12*1000/1000	# mol m-3
  gwDIC=0.7025#median(ntlGW$dic,na.rm=TRUE)/12*1000/1000	# mol m-3
  gwP=0.0007742#median(ntlGW$totp,na.rm=TRUE)/31*1000/1e6		# mol m-3
  precipDOC=3.195/12*1000/1000		#mol m-3; from UNDERC rain data; also see Likens et al. 1983 @ Hubbard Brook 1.1 mg L-1 
  # precipDIC=400/1e6*1/kH*1000	# mol C m-3; assumed in equilibrium with atmosphere
  precipDIC=1/12*1000/1000 # mol C m-3; precip DIC from Cardille et al. 2007
  precipP=0.01/31*1000/1000	# mol m-3; Murphy & DOskey 1976 JGLR
  snowDOC=0.6/12*1000/1000 # snow DOC from UNDERC snow; 0.5985 mg C L-1 to mol C m-3; snow DOC from Sebestyen lab
  snowDIC=400/1e6*1/kH*1000 # mol C m-3; assumed in equilibrium with atmosphere 
  snowP=0.01/31*1000/1000 # mol C m-3; Murphy &DOskey 1976 JGLR 
  # alk=-832+889.8*log10(curLakeInitGeomorph$Perim0)-358.1*log10(curLakeInitGeomorph$Area0) # microEq/L; alkalinity based on Hanson et al. 2004 and Cardille et al. 2007 ; see R code for linear model with area and perimeter 
  # alk=1795.634-895.801*log10(curLakeInitGeomorph$Area0)-2.541*log10(curLakeInitGeomorph$Perim0)+177.581*log10(curLakeInitGeomorph$Area0)*log10(curLakeInitGeomorph$Perim0)
  # alk=-512.5+152.81*log10(curLakeInitGeomorph$Area0)
  # if(alk<(100*-1)){
  #   alk=-85 # min alkalinity in Hanson et al. 2007 
  # }
  # initpH=2.66934+0.79234*log10(curLakeInitGeomorph$Area0) # predicting pH from lake area in Hanson et al. 2007; using that to calculate alkalinity and then letting pH vary in model
  # initDIC=(-6.168-6.361*log10(curLakeInitGeomorph$Area0)+12.955*log10(curLakeInitGeomorph$Perim0))/12 # mol m-3; predicting DIC pool from lake area and perimeter 
  # if(initDIC<(0.4/12)){
  #   initDIC=0.4/12 #mol m-3 minimum of hanson et al. 2007 
  # }
  
  # sensativity analysis for wetland loading 
  # curWatershed$percentWetland=100
  
  # initial epi / hypo volumes 
  r1=curLakeInitGeomorph$r2h*curLakeHydro$Stage[1]
  if(r1==0){ # for when stage is zero 
    r1=curLakeInitGeomorph$r2h*1
  }
  r2=curLakeInitGeomorph$r2h*(curLakeHydro$Stage[1]-zmix)
  if(r2<0){
    r2=0
  }
  Vepi=(1/3)*pi*(r1^2+r1*r2+r2^2)*zmix #truncated cone 
  Vhypo=(1/3)*pi*r2^2*(curLakeHydro$Stage[1]-zmix) # cone 
  Rn=0.466
  epiTemp=15 # starting at 15 degree C (in June)
  hypoTemp=7 # always 7 degree C in summer
  k=0.5 # m day-1; gas exchange coefficient 
  
  # use modeled initDIC to estimate an alkalinity that is held constant in the model 
  # alk=aquaenv(S=sal,t = epiTemp,Pa = 960/1000,d = 0,lat = curLakeLoc$Latitude,fCO2atm = 380/1e6,SumCO2 = initDIC/water.density(epiTemp,sal),
  #             pH = initpH)
  # alk=alk$TA[1]*water.density(epiTemp,sal)/1000*1e6
  # if(alk>2000){
  #   alk=2000 # max of hanson et al. 2007 
  # }
  # if(alk<(85*-1)){
  #   alk=-85 # min of hanson et al. 2007
  # }
  
  ## update model with with final hydrology data from Zach 
  if(mean(curLakeHydro$GWin+curLakeHydro$Baseflow)>0){
    alk=10^(3.2815*(1-exp(-1.2765*log10(mean(curLakeHydro$GWin+curLakeHydro$Baseflow+curLakeHydro$SWin))/3.2815))) 
  }else{
    alk=0 # setting alkalinity to zero if no GWin or Baseflow 
  }
  # if(alk>1700){
  #   alk=1700
  # }
  
  # initDIC=aquaenv(S=sal,t = epiTemp,Pa = 960/1000,d = 0,lat = curLakeLoc$Latitude,pH = initpH,fCO2atm = 380/1e6)
  initDIC=1/12
  
  # parameters that are important 
  params=c(gwIn0=curLakeHydro$GWin,gwOut0=curLakeHydro$GWout,precipDIC=precipDIC,precipDOC=precipDOC,precipP=precipP,streamDIC=streamDIC,
           streamDOC=streamDOC,streamPOC=streamPOC,streamP=streamP,gwDIC=gwDIC,snowDOC=snowDOC,snowDIC=snowDIC,snowP=snowP,
           gwDOC=gwDOC,gwP=gwP,kH=kH,r2h=curLakeInitGeomorph$r2h,zmix=zmix,rSlow=rSlow,rFast=rFast,fracLabile=fracLabile)
  
  X=c(phyto=0.4/12*Vepi,
      P_epi=5/1000/31*Vepi,P_hypo=5/1000/31*Vhypo,
      DIC_epi=initDIC*Vepi,DIC_hypo=initDIC*Vhypo,
      DOCr_epi=4/12*Vepi,DOCr_hypo=4/12*Vhypo,
      DOCl_epi=0.1/12*Vepi,DOCl_hypo=0.1/12*Vhypo,
      tPOC_epi=0.8/12*Vepi,tPOC_hypo=0.8/12*Vhypo,
      Emit=0,Sed_tPOC=0,Sed_phyto=0,zmix=zmix,kD=1,GPP=0,Vepi=Vepi,Vhypo=Vhypo,k=k,epiTemp=epiTemp,hypoTemp=hypoTemp,Burial_tPOC=0,Burial_phyto=0,
      DOC_Load=0,DOC_export=0,DOC_Respired=0,DOC_Respired_woExude=0,emergent_d_epi=0,emergent_d_hypo=0,pH=0,fracCO2=0) # pools of constituents in mol 
  # Q = 0.015
  t=nObs
  curForce$datetime<-as.POSIXct(paste(curForce$YYYY,curForce$MM,curForce$DD,curForce$HH),format='%Y %m %d %H',tz='GMT')
  
  # reducing snow depth because it's really high as of 2017-04-07
  # curLakeHydro$snowDepth_m=curLakeHydro$snowDepth_m*0.3
  
  outDir<-'D:/Jake/My Papers/NHLD Carbon Model/Results/20170818/'

  source('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/R Code/NHLD_LakeCarbonModel_supporting_V13_discrete.R')
  
  out<-matrix(NA,nrow=length(t),ncol=length(X))
  #initiallizing states 
  out[1,]<-X
  
  ptm=proc.time()
  out<-tryCatch(timeStep(t=t,out=out,params=params))
  proc.time()-ptm
  beep(sound = 2)

  out<-data.frame(out)
  colnames(out)<-names(X)
  out$time=t
  curLakeHydro$time<-t
  out<-merge(curLakeHydro,out,by='time',all.y=T)
  
  write.table(out,file.path(outDir,paste(lakes[i],"CRWA_V_adj_C_model.txt",sep="_")),row.names=FALSE,sep="\t",quote=F)
# }

windows()
par(mfrow=c(3,3))
plot(out$phyto/out$Vepi*12*1000/phyto_C2Chl,ylab='phyto ug L-1') # chl ug L-1
plot(out$P_epi/out$Vol*31*1000,ylab='P ug L-1') #ug L-1
plot(out$DIC_epi/out$Vepi*12,ylab='DIC mg L-1') #mg L-1
lines(out$co2*1/1e6/kH*1000*12) # mg L-1 atm equil 
plot((out$DOCr_epi+out$DOCl_epi)/out$Vepi*12,ylab='DOC mg L-1',type='l') #mg L-1
plot(out$DOCr_epi/out$Vepi*12,ylab='DOCr mg L-1') # mg L-1 
plot(out$DOCl_epi/out$Vepi*12,ylab='DOCl mg L-1') # mg L-1 
plot(out$tPOC_epi/out$Vepi*12,ylab='tPOC mg L-1') #mg L-1
plot(out$Emit/out$Area*12*1000,ylab='Emit mg C m-2 day-1',type='l') # mg C m-2 day-1
plot(out$Sed_phyto/out$Area*1000*12,ylab='Sed_phyto mg C m-2 day-1',col='green',lwd=2,type='l') # mg C m-2 day-1
lines(out$Sed_tPOC/out$Area*1000*12,ylab='Sed_tPOC mg C m-2 day-1',lwd=2) # mg C m-2 day-1

windows()
par(mfrow=c(2,2))
plot(out$zmix,ylim=c(max(out$zmix),0))
lines(out$Stage,col='red')
plot(out$pH,ylab='pH')
plot(out$fracCO2,ylab='fraction of DIC as CO2')
plot(out$DIC_epi/out$Vepi*12*out$fracCO2,ylab='CO2 mg L-1')

windows()
par(mfrow=c(3,3))
plot(out$kD)
plot(out$P_hypo[out$Vhypo>.1]/out$Vhypo[out$Vhypo>.1]*31*1000~t[out$Vhypo>.1],ylab='P ug L-1') #ug L-1
plot(out$DIC_hypo[out$Vhypo>.1]/out$Vhypo[out$Vhypo>.1]*12~t[out$Vhypo>.1],ylab='DIC mg L-1') #mg L-1
plot((out$DOCr_hypo[out$Vhypo>.1]+out$DOCl_hypo[out$Vhypo>.1])/out$Vhypo[out$Vhypo>.1]*12~t[out$Vhypo>.1],ylab='DOC mg L-1') #mg L-1
plot(out$tPOC_hypo[out$Vhypo>.1]/out$Vhypo[out$Vhypo>.1]*12~t[out$Vhypo>.1],ylab='tPOC mg L-1') #mg L-1
plot(out$Vol)
lines(out$Vepi+out$Vhypo,col='blue') # check to make sure we're conserving hydrologic mass balance 
plot(out$epiTemp,type='l')
plot(out$Vepi,ylim=c(0,max(out$Vepi)),type='l',ylab='Epi / Hypo Vol')
lines(out$Vhypo,col='red')
out$hydroOut<-out$LakeE+out$GWout+out$SWout # for residence time calc
mean(out$Vol)/mean(out$hydroOut)/365
plot(out$Burial_tPOC*1000*12/out$Area,type='l')
lines(out$Burial_phyto*1000*12/out$Area,col='green')
# 

windows()
plot(out$emergent_d_epi)
plot(out$emergent_d_hypo)
plot(out$DOC_Load)
plot(out$DOC_Respired)
points(out$DOC_Respired_woExude,col='red')
plot(out$DOC_export)

mean(1-(mean(out$DOC_export,na.rm = T)/mean(out$DOC_Load,na.rm = T)),na.rm = T) # fraction retained 
mean(mean(out$DOC_Respired,na.rm = T)/mean(out$DOC_Load,na.rm = T),na.rm = T)
mean(.952*out$fracCO2*out$DIC_epi/out$Vepi*1000*29.41,na.rm = T) # pco2
mean(out$k,na.rm = T)

windows()
ccf(out$emergent_d_epi[!is.na(out$emergent_d_epi)],out$DOC_Load[!is.na(out$emergent_d_epi)],lag.max = 20)

plot(out$emergent_d_epi)
plot(out$emergent_d_epi/(1.047^(out$epiTemp-20)))
plot(out$emergent_d_epi/(1.047^(out$epiTemp-20))~out$DOC_Load)
ccf(out$emergent_d_epi[!is.na(out$emergent_d_epi)]/(1.047^(out$epiTemp[!is.na(out$emergent_d_epi)]-20)),out$DOC_Load[!is.na(out$emergent_d_epi)])

par(new=T)
plot(out$phyto/out$Vepi,col='green') 


# windows()
# xlim=c(150,200)
# plot(out$DIC_epi/out$Vepi*12,ylab='DIC mg L-1',xlim=xlim,ylim=c(0,2)) #mg L-1
# par(new=T)
# plot(diff(out$Emit)/out$Area[2:length(out$time)]*12*1000,ylab='Emit mg C m-2 day-1',type='l',xlim=xlim,yaxt='n') # mg C m-2 day-1
# axis(4)
# par(new=T)


# windows()
# par(mfrow=c(3,3))
# plot(out$Vol)
# plot(out$GWout)
# plot(out$GWin)
# plot(out$SWin)
# plot(out$SWout)
# plot(out$DirectP)
# plot(out$LakeE)
# plot(out$IceMelt)
# plot(diff(out$Emit)/out$Area[2:length(out$time)]*12*1000,ylab='Emit mg C m-2 day-1',type='l',ylim=c(0,max(diff(out$Emit)/out$Area[2:length(out$time)]*12*1000))) # mg C m-2 day-1
# abline(0,0,lty=2,lwd=2)

sum(out$Burial_tPOC,out$Burial_phyto,na.rm = T)/sum(out$Sed_tPOC,out$Sed_phyto,na.rm = T)

windows()
plot(out$GPP*12*out$zmix*1000~as.POSIXct(out$datetime),
     type='l',lwd=1,ylab='GPP mg C m-2 day-1')

mean(out$Sed_tPOC/out$Area*1000*12+out$Sed_phyto/out$Area*1000*12,na.rm = T) # mg C m-2 day-1
mean(out$Burial_phyto/out$Area*1000*12+out$Burial_tPOC/out$Area*1000*12,na.rm = T)
mean(out$GPP*12*out$zmix*1000,na.rm = T)
mean(out$Emit/out$Area*12*1000,na.rm = T)
mean((out$DOCr_epi+out$DOCl_epi)/out$Vepi*12,na.rm = T)

# 
# plot(out$GPP/32*12~as.POSIXct(out$datetime),
#      type='l',lwd=1,ylab='GPP mg O2 L-1 day-1')


# LTER NTL validation chem data 
ntl<-read.csv('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Validation Data/NTL LTER Chemistry/chemical_limnology_of_north_temperate_lakes_lter_primary_study_lakes__nutrients_ph_and_carbon.csv',
              stringsAsFactors = F)
ntlChl<-read.csv('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Validation Data/NTL LTER Chl/north_temperate_lakes_lter__chlorophyll_-_trout_lake_area.csv',
                 stringsAsFactors = F)
ntlTemp<-read.csv('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Validation Data/physical_limnology_of_the_north_temperate_lakes_primary_study_lakes.csv',
                  stringsAsFactors = F)
ntl<-ntl[ntl$lakeid%in%ntlLookUp$lakeID,]
ntl$sampledate<-as.POSIXct(ntl$sampledate,format='%Y-%m-%d')
ntlChl<-ntlChl[ntlChl$lakeid%in%ntlLookUp$lakeID,]
ntlChl$sampledate<-as.POSIXct(ntlChl$sampledate, format='%Y-%m-%d')
ntlTemp<-ntlTemp[ntlTemp$lakeid==lake,]
ntlTemp$sampledate<-as.POSIXct(ntlTemp$sampledate, format='%Y-%m-%d')
ntlTemp<-ntlTemp[sort.list(ntlTemp$sampledate),]

############ Validation Figs
windows()
xlim=c(as.POSIXct('1985-01-01 00:00:00'),max(as.POSIXct(out$datetime)))
plot(ntl$doc[ntl$lakeid==lake&ntl$depth==0&!is.na(ntl$doc)&ntl$sampledate<xlim[2]]~
       ntl$sampledate[ntl$lakeid==lake&ntl$depth==0&!is.na(ntl$doc)&ntl$sampledate<xlim[2]],
     col='grey60',xlim=xlim,pch=16,lwd=3,ylab='',xlab='',cex.axis=2,xaxt='n',type='o')
lines((out$DOCr_epi+out$DOCl_epi)/(out$Vepi)*12~as.POSIXct(out$datetime),xlim=xlim,ylab='DOC mg L-1',type='l',lwd=4,xlab='Date',cex.axis=2,xaxt='n')
axis(1,at=c(as.POSIXct('1990-01-01'),as.POSIXct('2000-01-01'),as.POSIXct('2010-01-01')),labels = c(1990,2000,2010),cex.axis=2)
legend('topleft',legend = c('Modeled','Observed'),fill = c('black','grey60'),bty='n')

windows()
xlim=c(as.POSIXct('1985-01-01 00:00:00'),max(as.POSIXct(out$datetime)))
plot(ntl$totpf[ntl$lakeid==lake&ntl$depth==0&!is.na(ntl$totpf)&ntl$sampledate<xlim[2]]~
       ntl$sampledate[ntl$lakeid==lake&ntl$depth==0&!is.na(ntl$totpf)&ntl$sampledate<xlim[2]],ylim=c(0,50),
     col='grey60',xlim=xlim,pch=16,lwd=3,ylab='SRP ug L-1',xlab='',xaxt='n',type='o',cex.axis=2)
lines(out$P_epi/out$Vepi*31*1000~as.POSIXct(out$datetime),xlim=xlim,lwd=4,xlab='Date',cex.axis=2,xaxt='n')
axis(1,at=c(as.POSIXct('1990-01-01'),as.POSIXct('2000-01-01'),as.POSIXct('2010-01-01')),labels = c(1990,2000,2010),cex.axis=2)
legend('topleft',legend = c('Modeled','Observed'),fill = c('black','grey60'),bty='n')

windows()
xlim=c(as.POSIXct('1985-01-01 00:00:00'),max(as.POSIXct(out$datetime)))
plot(ntl$totpuf[ntl$lakeid==lake&ntl$depth==0&!is.na(ntl$totpuf)&ntl$sampledate<xlim[2]]~
       ntl$sampledate[ntl$lakeid==lake&ntl$depth==0&!is.na(ntl$totpuf)&ntl$sampledate<xlim[2]],ylim=c(0,50),
     col='grey60',xlim=xlim,pch=16,lwd=3,ylab='TP ug L-1',xlab='',xaxt='n',type='o',cex.axis=2)
lines((out$P_epi+out$phyto/phytoC2P)/out$Vepi*31*1000~as.POSIXct(out$datetime),xlim=xlim,lwd=4,xlab='Date',cex.axis=2,xaxt='n')
axis(1,at=c(as.POSIXct('1990-01-01'),as.POSIXct('2000-01-01'),as.POSIXct('2010-01-01')),labels = c(1990,2000,2010),cex.axis=2)
legend('topleft',legend = c('Modeled','Observed'),fill = c('black','grey60'),bty='n')

windows()
xlim=c(as.POSIXct('1985-01-01 00:00:00'),max(as.POSIXct(out$datetime)))
plot(ntl$dic[ntl$lakeid==lake&ntl$depth==0&!is.na(ntl$dic)&ntl$sampledate<xlim[2]]~
       ntl$sampledate[ntl$lakeid==lake&ntl$depth==0&!is.na(ntl$dic)&ntl$sampledate<xlim[2]],ylim=c(0,max(ntl$dic[ntl$lakeid==lake&ntl$depth==0&!is.na(ntl$dic)&ntl$sampledate<xlim[2]])),
     col='grey60',xlim=xlim,pch=16,lwd=3,ylab='DIC ug L-1',xlab='',xaxt='n',type='o',cex.axis=2)
lines(out$DIC_epi/out$Vepi*12~as.POSIXct(out$datetime),xlim=xlim,lwd=4,xlab='Date',cex.axis=2,xaxt='n')
axis(1,at=c(as.POSIXct('1990-01-01'),as.POSIXct('2000-01-01'),as.POSIXct('2010-01-01')),labels = c(1990,2000,2010),cex.axis=2)
legend('topleft',legend = c('Modeled','Observed'),fill = c('black','grey60'),bty='n')

windows()
xlim=c(as.POSIXct('1985-01-01 00:00:00'),max(as.POSIXct(out$datetime)))
plot(ntl$ph[ntl$lakeid==lake&ntl$depth==0&!is.na(ntl$dic)&ntl$sampledate<xlim[2]]~
       ntl$sampledate[ntl$lakeid==lake&ntl$depth==0&!is.na(ntl$dic)&ntl$sampledate<xlim[2]],
     ylim=c(min(ntl$ph[ntl$lakeid==lake&ntl$depth==0&!is.na(ntl$dic)&ntl$sampledate<xlim[2]],na.rm = T),max(ntl$ph[ntl$lakeid==lake&ntl$depth==0&!is.na(ntl$dic)&ntl$sampledate<xlim[2]],na.rm = T)),
     col='grey60',xlim=xlim,pch=16,lwd=3,ylab='pH',xlab='',xaxt='n',type='o',cex.axis=2)
lines(out$pH~as.POSIXct(out$datetime),xlim=xlim,lwd=4,xlab='Date',cex.axis=2,xaxt='n')
axis(1,at=c(as.POSIXct('1990-01-01'),as.POSIXct('2000-01-01'),as.POSIXct('2010-01-01')),labels = c(1990,2000,2010),cex.axis=2)
legend('topleft',legend = c('Modeled','Observed'),fill = c('black','grey60'),bty='n')

windows()
xlim=c(as.POSIXct('1985-01-01 00:00:00'),max(as.POSIXct(out$datetime)))
plot(ntl$alk[ntl$lakeid==lake&ntl$depth==0&!is.na(ntl$alk)&ntl$sampledate<xlim[2]]~
       ntl$sampledate[ntl$lakeid==lake&ntl$depth==0&!is.na(ntl$alk)&ntl$sampledate<xlim[2]],
     ylim=c(min(ntl$alk[ntl$lakeid==lake&ntl$depth==0&!is.na(ntl$alk)&ntl$sampledate<xlim[2]],na.rm = T),max(ntl$alk[ntl$lakeid==lake&ntl$depth==0&!is.na(ntl$alk)&ntl$sampledate<xlim[2]],na.rm = T)),
     col='grey60',xlim=xlim,pch=16,lwd=3,ylab='alk',xlab='',xaxt='n',type='o',cex.axis=2)
abline(h=alk)
axis(1,at=c(as.POSIXct('1990-01-01'),as.POSIXct('2000-01-01'),as.POSIXct('2010-01-01')),labels = c(1990,2000,2010),cex.axis=2)
legend('topleft',legend = c('Modeled','Observed'),fill = c('black','grey60'),bty='n')

# 
# windows()
# plot(out$GPP*out$zmix*12*1000~as.POSIXct(out$datetime),
#      type='l',lwd=1,xlim=xlim,ylab='GPP mg C m-2 day-1')

windows()
plot(ntlChl$chlor[ntlChl$lakeid==lake]~ntlChl$sampledate[ntlChl$lakeid==lake],ylim=c(0,max(ntlChl$chlor[ntlChl$lakeid==lake])),type='l')
lines(out$phyto/out$Vepi*12*1000/phyto_C2Chl~as.POSIXct(out$datetime),col='red',lwd=2)

plot(ntlChl$chlor[ntlChl$lakeid==lake]~ntlChl$sampledate[ntlChl$lakeid==lake],ylim=c(0,20),type='l')
lines(out$phyto/out$Vepi*12*1000/phyto_C2Chl~as.POSIXct(out$datetime),col='red',lwd=2)
plot(ntlChl$chlor[ntlChl$lakeid==lake]~ntlChl$sampledate[ntlChl$lakeid==lake],ylim=c(0,40),type='l',xlim=range(as.POSIXct(out$datetime)))
lines(out$phyto/out$Vepi*12*1000/phyto_C2Chl~as.POSIXct(out$datetime),col='red',lwd=2)

windows()
plot(ntlTemp$wtemp[ntlTemp$depth==1]~ntlTemp$sampledate[ntlTemp$depth==1],type='l')
lines(out$epiTemp~as.POSIXct(out$datetime),col='red')
tempMerged<-out[,c('datetime','epiTemp')]
tempMerged$sampledate<-as.Date(tempMerged$datetime)
ntlTemp$sampledate<-as.Date(ntlTemp$sampledate)
tempMerged<-merge(tempMerged,ntlTemp[ntlTemp$depth==1,],by='sampledate',all.x=T)
plot(tempMerged$epiTemp~tempMerged$wtemp)
abline(0,1)
sqrt(mean((tempMerged$epiTemp-tempMerged$wtemp)^2,na.rm = T)) #rmse for epi temp 

# sedimentation validation 
sum(out$Burial_tPOC,out$Burial_phyto,na.rm = T)/sum(out$Sed_tPOC,out$Sed_phyto,na.rm = T) # burial efficiency 
mean(out$Sed_tPOC/out$Area*1000*12+out$Sed_phyto/out$Area*1000*12,na.rm = T) # mg C m-2 day-1
mean(out$Burial_phyto/out$Area*1000*12+out$Burial_tPOC/out$Area*1000*12,na.rm = T)






