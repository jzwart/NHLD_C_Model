##### lake carbon model for NHLD from VIC / GFLOW input 
# JAZ; 2016-02-10 using code originally from SEJ
# V10 has modeled zmix, kd, and a 2-box lake model and wetland loading and depth integrated GPP and light attenuated by snow depth and updated water temp budget; 
#  2 box model of DOC (labile and recalcitrant) ; splitting stream water into equal proportions of epi and hypo ; wind sheltering coefficient for k600 

# clear workspace
rm(list=ls())

########## load utility functions and packages
source('/afs/crc.nd.edu/user/j/jzwart/NHLDLakeCarbonModel/R_code/NHLD_LakeCarbonModel_supporting_V13_discrete.R')
source('/afs/crc.nd.edu/user/j/jzwart/NHLDLakeCarbonModel/R_code/AveLightClimate.R')
source('/afs/crc.nd.edu/user/j/jzwart/NHLDLakeCarbonModel/R_code/DOY.r')
require(deSolve)
require(LakeMetabolizer)
require(snow)
require(sp)
require(rgeos)
require(parallel)
require(rgdal)
require(maptools)
require(AquaEnv)
require(marelac)

dir<-'/afs/crc.nd.edu/user/j/jzwart/NHLDLakeCarbonModel/Data/20170815'

# forcings directory (from VIC)
forceDir<-file.path(dir,'Force')
# daily lake hydrology flux directory
lakeFluxDir<-file.path(dir,'OUT_RESULTS')
# daily lake geomorphology directory 
lakeGeomorphDir<-file.path(dir,'OUT_RESULTS')
# daily lake initial conditions directory 
lakeInitGeomorphDir<-file.path(dir,'OUT_RESULTS')
# lake watershed information 
watersheds<-read.table(file.path(dir,"NHLDsheds_20170323.txt"),header=TRUE,sep="\t",stringsAsFactors=FALSE)
watersheds$percentWetland<-watersheds$percentWetland*100
# lake location, area, and perimeter 
lakeLocation<-readShapeSpatial(file.path(dir,'LakeLocations/NHLDandBuffLakes_Original_ZJH_Rev1.shp'))
lakeLocation<-data.frame(lakeLocation)
lakeLocation$Permanent_<-as.character(lakeLocation$Permanent_)

# keeling curve for historic CO2 concentrations 
keeling<-read.table(file.path(dir,'keelingCurve.txt'),header = F,stringsAsFactors = F)
colnames(keeling)<-c('site_code','year','month','day','hour','minute','second','co2','value_unc','nvalue','latitude','longitude','altitude','elevation','intake_height','instrument','qcflag')
keeling$datetime<-as.Date(paste(keeling$year,keeling$month,keeling$day),format='%Y %m %d')
keeling$co2<-ifelse(keeling$co2==-999.99,NA,keeling$co2)
keeling$co2[1]<-keeling$co2[min(which(!is.na(keeling$co2)))]
keeling$co2[length(keeling$site_code)]<-keeling$co2[max(which(!is.na(keeling$co2)))]
keeling$co2<-approx(keeling$datetime,keeling$co2,keeling$datetime)$y #filling in NA's / missing data 
keeling<-keeling[,c('datetime','co2')]

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

# snow depth on top of ice; which is the same for all lakes 
snowDepth=read.table(file.path(dir,'Lake_Snow_Depth_Land_SWE_Depth_Ratio_Correction.txt'),sep = '\t')
snowDepth$datetime<-as.Date(paste(snowDepth$V1,snowDepth$V2,snowDepth$V3),format = '%Y %m %d')
snowDepth=snowDepth[,c('datetime','V4')]
colnames(snowDepth)[2]<-'snowDepth_m' # m; snow depth on ice 

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
phytoRecycling=0.2 # 1 - recyling rate of phytoplankton; i.e. 0.05 means 95% recycling rate 
phytoC2P=95#106		#M:M; from Patrick
lossPhyto=0.1 # loss rate of phytoplankton [day-1]
GPPexudeR=0.03  # Hanson et al. 2004; DOC exude from phytoplankton that is raclcitrant 
GPPexudeL=0.07 # DOC exude from phytoplankton that is labile 
Rauto=0.8 	# Hanson et al. 2004; quick respiration of phytoplankton production 
phytoDeath=0.03 # Hanson et al 2004; fraction of phyto that die; day-1
phyto_C2Chl=50 # phytoplankton carbon to chlorophyll ratio, [gC gChl-1]; sort of made up/average of observations in paper (e.g. Geider 1987)
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
hs=7.5 # m; average elevation difference between top of surrounding canopy and lake surface; about average value from Read et al. 2014
wLoad=3/365 # wetland load per shoreline length; g C m-1 shoreline year-1 converted to day-1 (mean from fitted parameter in Hanson et al. 2014) 
sal=0 # salinity set to 0 for all lakes 
leafLoad=300/12 #mol C m-1 shoreline yr-1; autumn leaf fall; check out Likens 1985 and Gasith and Hasler 1976; Hanson et al 2014; France 1996; France and Peters 1996; 300 g C m-1 shoreline yr-1 is roughly the average for all these studies 

date<-read.table(file.path(dir,'Date_OUT_ALL.txt'),stringsAsFactors = F,header = F,sep='\t')
datetime<-as.Date(paste(date$V1,date$V2,date$V3),format = '%Y %m %d')

lakeSim<-function(curLakeID){
  # load in current lake hydrology output data 
  curLakeHydro<<-read.table(file.path(lakeFluxDir,paste(curLakeID,'_DailyLakeFlux.txt',sep='')),header=F,stringsAsFactors = F)
  curLakeHydro<<-cbind(datetime,curLakeHydro)
  colnames(curLakeHydro)<<-c('datetime','SWin','DirectP','LakeE','GWin','GWout','SWout','IceSnow','LandMelt','IceMelt','SWoutMinusLandMelt','LandMeltEst','Baseflow')
  # load in current lake geomorphology output data 
  curLakeGeomorph<<-read.table(file.path(lakeGeomorphDir,paste(curLakeID,'_DailyLakeProp.txt',sep='')),header=F,stringsAsFactors = F)
  curLakeGeomorph<<-cbind(datetime,curLakeGeomorph)
  colnames(curLakeGeomorph)<<-c('datetime','Area','Radius','Perim','Stage','Elev','Vol')
  curLakeHydro<<-merge(curLakeHydro,curLakeGeomorph,by='datetime',all.x=T)
  rm(curLakeGeomorph)
  curLakeHydro<<-merge(curLakeHydro,snowDepth,by='datetime',all.x=T)
  # load in current lake initial geopmorphological parameters 
  curLakeInitGeomorph<<-read.table(file.path(lakeInitGeomorphDir,paste(curLakeID,'_IniLakeProp.txt',sep='')),header=F,stringsAsFactors = F)
  colnames(curLakeInitGeomorph)<<-c('Area0','Vol0','Radius0','Diameter0','Perim0','Stage0','DL','r2h','WA','WALA','Elev0_DEM','Vol_LinRes','Stage_LinRes','Stream_WA')
  # watershed data 
  curWatershed<<-watersheds[watersheds$Permanent_==curLakeID,]
  
  # location of lake used for finding correct forcing data 
  curLakeLoc<<-lakeLocation[lakeLocation$Permanent_==curLakeID,]
  sp.curLakeLoc<<-curLakeLoc
  coordinates(sp.curLakeLoc)<<- ~Longitude+Latitude
  # find closest forcing data 
  curForceCoord<<-forceCoord[which(gDistance(sp.forceCoord,sp.curLakeLoc,byid=T)==gDistance(sp.forceCoord,sp.curLakeLoc)),]
  # curForce<-read.table(file.path(forceDir,paste('FORCE',curForceCoord$Latitude,curForceCoord$Longitude,sep='_')),header=F,stringsAsFactors = F)
  curForce<<-read.table(file.path(forceDir,paste('FORCE',curForceCoord$Latitude,curForceCoord$Longitude,sep='_')),header=F,stringsAsFactors = F)
  colnames(curForce)<<-c('YYYY','MM','DD','HH','Precip','AirTemp','Longwave','Shortwave','Density','VPD','Pressure','Wind')
  curForce$datetime<<-as.Date(paste(curForce$YYYY,curForce$MM,curForce$DD,curForce$HH),format='%Y %m %d %H')
  # keep forcing data for days over which you want to model 
  curForce<<-curForce[curForce$datetime>as.Date(paste(startYear,startMonth,startDay),format='%Y %m %d')&
                       curForce$datetime<as.Date(paste(endYear,endMonth,endDay),format='%Y %m %d'),]
  
  #get forcing data for lake dates 
  curForce<<-curForce[curForce$datetime%in%curLakeHydro$datetime,]
  curLakeHydro<<-curLakeHydro[curLakeHydro$datetime%in%curForce$datetime,]
  curForce$I0<<-sw.to.par(curForce,sw.col = 'Shortwave')$par # converting from shortwave radiation to PAR from LakeAnalyzer 
  curForce$I0<<-ifelse(curForce$I0==0,0.001,curForce$I0)
  #merging keeling curve and hydrologic forcings 
  curLakeHydro<<-merge(curLakeHydro,keeling,by='datetime',all.x=T)
  curLakeHydro$co2[1]<<-curLakeHydro$co2[min(which(!is.na(curLakeHydro$co2)))]
  curLakeHydro$co2[length(curLakeHydro$datetime)]<<-curLakeHydro$co2[max(which(!is.na(curLakeHydro$co2)))]
  curLakeHydro$co2<<-approx(curLakeHydro$datetime,curLakeHydro$co2,curLakeHydro$datetime)$y
  
  nObs<<-1:length(curLakeHydro$datetime)
  nObs_hourly<<-rep(nObs,each=24)

  maxWind<<-tapply(curForce$Wind,nObs_hourly,FUN=mean)
  #wind sheltering coefficient parameters for Read et al. 2014 equation 
  # wind sheltering coefficient from Read et al. 2014 
  if(curLakeInitGeomorph$Area0-625*hs^2*pi<0){
    Ws<<-0.001 # extremely small if small lake 
  }else{
    Ws<<-(2/pi)*acos(25*hs*sqrt((pi/curLakeInitGeomorph$Area0)))-(50*hs)/(curLakeInitGeomorph$Area0*sqrt(pi))*
      sqrt(curLakeInitGeomorph$Area0-625*hs^2*pi)
  }
  lakeWind<<-maxWind*Ws # adjusting lake wind by wind sheltering coefficient 
  
  # calculate concentrations of DIC, DOC, TP based on landcover
  # could use Frost et al 2005 for stream DOC; and Johnston et al 2008 ; Jutras etal 2011 
  # could use Frost et al 2009 for stream C:P of seston 
  streamDOC<<-exp(1.3961+3.245*(curWatershed$percentWetland/100))/12*1000/1000	# from lottig 2012; mol m-3
  streamPOC<<-3/12*1000/1000		# ~3 mg L-1; buffam 2011; mol m-3
  streamDIC<<-10.85/12*1000/1000		#10 mg L-1; lottig 2011; mol m-3
  streamP<<-0.05/31*1000/1000#0.04/31*1000/1000	  # Long inlet has 84 ug/L	#.025 mg L-1 TDP & 0.04 mg L-1 TP; lottig 2011; mol m-3
  gwDOC<<-13/12*1000/1000#median(ntlGW$doc,na.rm<<-TRUE)/12*1000/1000	# mol m-3
  gwDIC<<-0.7025#median(ntlGW$dic,na.rm<<-TRUE)/12*1000/1000	# mol m-3
  gwP<<-0.0007742#median(ntlGW$totp,na.rm<<-TRUE)/31*1000/1e6		# mol m-3
  precipDOC<<-3.195/12*1000/1000		#mol m-3; from UNDERC rain data; also see Likens et al. 1983 @ Hubbard Brook 1.1 mg L-1 
  # precipDIC<<-400/1e6*1/kH*1000	# mol C m-3; assumed in equilibrium with atmosphere
  precipDIC<<-1/12*1000/1000 # mol C m-3; precip DIC from Cardille et al. 2007
  precipP<<-0.01/31*1000/1000	# mol m-3; Murphy & DOskey 1976 JGLR
  snowDOC<<-0.6/12*1000/1000 # snow DOC from UNDERC snow; 0.5985 mg C L-1 to mol C m-3; snow DOC from Sebestyen lab
  snowDIC<<-400/1e6*1/kH*1000 # mol C m-3; assumed in equilibrium with atmosphere 
  snowP<<-0.01/31*1000/1000 # mol C m-3; Murphy &DOskey 1976 JGLR 
  if(mean(curLakeHydro$GWin+curLakeHydro$Baseflow)>0){
    alk<<-10^(3.2815*(1-exp(-1.2765*log10(mean(curLakeHydro$GWin+curLakeHydro$Baseflow+curLakeHydro$SWin))/3.2815))) 
  }else{
    alk<<-0 # setting alkalinity to zero if no GWin or Baseflow 
  }
  
  initDIC<<-1/12
  
  # high wetland load so making entire shoreline wetland 
  curWatershed$percentWetland<<-100
  
  # initial epi / hypo volumes 
  r1<<-curLakeInitGeomorph$r2h*curLakeHydro$Stage[1]
  if(r1==0){ # for when stage is zero 
    r1<<-curLakeInitGeomorph$r2h*1
  }
  r2<<-curLakeInitGeomorph$r2h*(curLakeHydro$Stage[1]-zmix)
  if(r2<0){
    r2<<-0
  }
  Vepi<<-(1/3)*pi*(r1^2+r1*r2+r2^2)*zmix #truncated cone 
  Vhypo<<-(1/3)*pi*r2^2*(curLakeHydro$Stage[1]-zmix) # cone 
  Rn<<-0.466
  epiTemp<<-15 # starting at 15 degree C (in June)
  hypoTemp<<-7 # always 7 degree C in summer
  k<<-0.5 # m day-1; gas exchange coefficient 
  
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
  curForce$datetime<<-strftime(strptime(paste(curForce$YYYY,curForce$MM,curForce$DD,curForce$HH),format = '%Y %m %d %H',tz = 'GMT'))
  
  outDir='/afs/crc.nd.edu/user/j/jzwart/NHLDLakeCarbonModel/Results/C_Model_Output/20170829'

  out=matrix(NA,nrow=length(t),ncol=length(X))
  #initiallizing states 
  out[1,]=X

  out=tryCatch(timeStep(t=t,out=out,params=params))

  out=data.frame(out)
  colnames(out)=names(X)
  out$time=t
  curLakeHydro$time=t
  out=merge(curLakeHydro,out,by='time',all.y=T)
  
  write.table(out,file.path(outDir,paste(curLakeID,"C_model.txt",sep="_")),row.names=FALSE,sep="\t",quote=F) 
}

lakes=as.list(lakes[lakes%in%watersheds$Permanent_])
# lakes=as.list(unlist(lakes)[1:800])

mclapply(lakes,lakeSim,mc.cores=48)
