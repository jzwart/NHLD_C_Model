### Function for calculating light at depth (Iz) and average light climate based on incoming light and light attenuation (kD) or DOC and chlorophyll

# support function that does Iz calc.; allows integration in main function
lightAtten<-function(z,I0,kD){
	Iz=I0*exp(-kD*z)
	return(Iz)
}

# calculates Iz and mean light climate (Ibar)
# minimum arguments required are incident light (I0), depth to calculate light at (z), and light attenuation coefficient (kD) (or DOC [mg/L] & chlorophyll [ug/L] concentrations from which a kD is estimated)
# can take a single or vector of depths and incident lights
# returns a list with to elements light at depth(s) z (Iz) and mean light climate in the mixed layer (zmix) if zmix is specified
lightProfileCalc<-function(I0,z,kD=NULL,DOC=NULL,chloro=NULL,zmix=NULL){
	Ibar=NULL
	
	if(is.null(kD)){
		# from GLEON 25 lakes (Solomon et al 2013) and UNDERC lakes (M2M database)
		kD=-0.217+0.0537*chloro+0.186*DOC
	}
	
	if(length(I0)==1){
		Iz=lightAtten(z=z,I0=I0,kD=kD)
		if(!is.null(zmix)){
			Ibar=integrate(lightAtten,0,zmix,I0=I0,kD=kD)$value/zmix
		}
	}else{
		Iz=matrix(NA,length(z),length(I0))
		Ibar=numeric(length(I0))
		for(i in 1:ncol(Iz)){
			Iz[,i]=lightAtten(z,I0[i],kD)
			Ibar[i]=integrate(lightAtten,0,zmix,I0=I0[i],kD=kD)$value/zmix
		}
	}
	
	return(list(Iz=Iz,Ibar=Ibar))
}

require(LakeMetabolizer) 
# time series wrapper for light profile calculation; JAZ 2015-01-16 
ts.lightProfileCalc<-function(ts.data){
  data<-ts.data
  
  if(has.vars(data,'I0')){
    I0=get.vars(data,'I0')
  }else{
    stop('Data must have I0 column\n')
  }
  
  if(has.vars(data,'z')){
    z=get.vars(data,'z')
  }else{
    stop('Data must have z column\n')
  }
  
  if(has.vars(data,'zmix')){
    zmix=get.vars(data,'zmix')
  }else{
    stop('Data must have zmix column\n')
  }
  
  if(has.vars(data,'kD')){
    kD=get.vars(data,'kD')
  }else if(has.vars(data,'DOC')){
    DOC=get.vars(data,'DOC')
    chloro=get.vars(data,'chloro')
  }else{
    stop('Data must have kD or DOC and chloro column\n')
  }
  
  if(has.vars(data,'kD')){
    lightProf<-unlist(Map(f=lightProfileCalc,I0=I0[,2],z=z[,2],kD=kD[,2],zmix=zmix[,2]))
    lightClimate<-as.numeric(lightProf[names(lightProf)=='Ibar'])
    light_z<-as.numeric(lightProf[names(lightProf)=='Iz'])
  }else if(has.vars(data,'DOC')){
    lightProf<-unlist(Map(f=lightProfileCalc,I0=I0[,2],z=z[,2],DOC=DOC[,2],chloro=chloro[,2],zmix=zmix[,2]))
    lightClimate<-as.numeric(lightProf[names(lightProf)=='Ibar'])
    light_z<-as.numeric(lightProf[names(lightProf)=='Iz'])
  }
  
  return(data.frame(datetime=data$datetime,lightClimate=lightClimate,light_z=light_z))
  
}



