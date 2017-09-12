#nhldWatershedModelSupporting.R supporting functions
#light attenuation
lightAtten<-function(z,I0,kD){
  Iz=I0*exp(-kD*z)
  return(Iz)
}
lightAttenTS<-function(I0,kD,zmix){
  avgI=numeric(length(I0))
  for(i in 1:length(avgI)){
    avgI[i]=integrate(lightAtten,0,zmix,I0=I0[i],kD=kD)$value/zmix
  }
  return(avgI)
}

# dailyGPP based on hourly light, etc.
dailyGPP<-function(day,curForce,curForceDOY,chlCur,SRP,DOC,V,kD,zmix,snowDepth){
  hourlyPAR=lightAtten(snowDepth,curForce$I0[curForceDOY==day],kdSnow)	# umol m2 sec; SW to PAR based on Read...
  # adjusting for snow cover during winter 
    
  #### trying different GPP formulation
  PPmax=3.1#2.2		# hr-1; 2.2 hr-1 x Chl comes from Vadeboncoeur et al. 2008 from Guildford et al 1994 
  kSRP=0.3/1000	# mol P m-3; Halmann and Stiller 1974 L&O 19(5): 774-783 (Lake Kinnerrett) use 0.3/1000
  Ik=180 		# light limitation benchmark for GPP, [umol cm-2 s-1]; Vadeboncoeur et al. 2008
  
  avgI=lightAttenTS(hourlyPAR,kD,zmix)		# hourly average light climate in mixed layer [umol cm-2 s-1]
  GPP=sum(chlCur*PPmax*tanh(avgI/Ik)*((SRP/V)/((SRP/V)+kSRP)))/(12*1000) # mol C m-3 day-1
  
  return(GPP)
}

# dailyGPP based on hourly light, etc.
pointGPP<-function(z,PAR,chlCur,SRP,V,kD){
  # hourlyPAR=curForce$I0[curForceDOY==day]	# umol m2 sec; SW to PAR based on Read...
  
  #### trying different GPP formulation
  PPmax=3.1#2.2		# hr-1; 2.2 hr-1 x Chl comes from Vadeboncoeur et al. 2008 from Guildford et al 1994 
  kSRP=0.3/1000	# mol P m-3; Halmann and Stiller 1974 L&O 19(5): 774-783 (Lake Kinnerrett)
  Ik=180 		# light limitation benchmark for GPP, [umol cm-2 s-1]; Vadeboncoeur et al. 2008
  
  Iz=lightAtten(z,PAR,kD) # hourly light at depth z [umol cm-2 s-1]
  GPPpoint=chlCur*PPmax*tanh(Iz/Ik)*((SRP/V)/((SRP/V)+kSRP))/(12*1000) # mol C m-3 day-1
  
  return(GPPpoint)
}

integratedGPP<-function(day,curForce,curForceDOY,chlCur,SRP,DOC,V,kD,zmix,snowDepth){
  hourlyPAR=lightAtten(snowDepth,curForce$I0[curForceDOY==day],kdSnow)	# umol m2 sec; SW to PAR based on Read...
  # adjusting for snow cover during winter 
  intGPP=numeric(length(hourlyPAR))
  for(i in 1:length(intGPP)){
    intGPP[i]=integrate(f = pointGPP,lower = 0,upper = zmix,PAR=hourlyPAR[i],chlCur=chlCur,SRP=SRP,V=V,kD=kD)$value
  }
  GPPout=sum(intGPP)
  return(GPPout)
}

# atmospheric deposition
depWithDist<-function(distShore,precip,maxWind){
  lPOC=10^(0.43+0.0034*precip+0.11*maxWind-1.05*distShore/(48+distShore))/(12*1000)	#mol C m-2 day-1
  sPOC=10^(0.0082+0.0068*precip+0.12*maxWind-0.6*distShore/(49+distShore))/(12*1000) #mol C m-2 day-1
  
  tPOCdep=(lPOC+sPOC) #tPOCdep summed across lake surface	#mol C m-2 day-1
  return(tPOCdep)
}

# daily tPOC deposition based on Preston et al.
dailyTPOCdep<-function(day,Area,curForce,directPrecip,nObs_hourly){
  maxWind=max(curForce$Wind[nObs_hourly==day])*(1/10)^0.15		# convert 10m wind to 1m for preston model tPOC deposition
  
  tPOCdepTotal=integrate(depWithDist,0.01,sqrt(Area/pi),precip=(directPrecip/Area/1000),maxWind=maxWind)$value/sqrt(Area/pi)*Area #tPOCdep summed across lake surface	#mol C day-1
  return(tPOCdepTotal)
}

rFunc<-function(r20,temp){
  r<-r20*1.047^(temp-20) # scaled R based on Holtgrieve et al. 2010 
  return(r)
}

timeStep<-function(t,out,params){
  
    #########################
    #### State variables ####
    #########################
    for(i in 2:length(t)){
    
    # states of previous timestep 
    phyto=out[i-1,1]
    P_epi=out[i-1,2]
    P_hypo=out[i-1,3]
    DIC_epi=out[i-1,4]
    DIC_hypo=out[i-1,5]
    DOCr_epi=out[i-1,6]
    DOCr_hypo=out[i-1,7]
    DOCl_epi=out[i-1,8]
    DOCl_hypo=out[i-1,9]
    tPOC_epi=out[i-1,10]
    tPOC_hypo=out[i-1,11]
    Emit=out[i-1,12]
    Sed_tPOC=out[i-1,13]
    Sed_phyto=out[i-1,14]
    zmix=out[i-1,15]
    kD=out[i-1,16]
    GPP=out[i-1,17]
    Vepi=out[i-1,18]
    Vhypo=out[i-1,19]
    if(Vhypo<=0){
      Vhypo=0.00001
      out[i-1,19]=Vhypo
    }
    kCur=out[i-1,20]
    epiTemp=out[i-1,21]
    hypoTemp=out[i-1,22]
    Burial_tPOC=out[i-1,23]
    Burial_phyto=out[i-1,24]
    
#     Q=X[9]
    
    wind=maxWind[i-1] # m sec-1 
    U10=lakeWind[i-1] # m sec-1
    SWin=curLakeHydro$SWin[i-1] # m3 
    Baseflow=curLakeHydro$Baseflow[i-1] # m3
    SWout=curLakeHydro$SWout[i-1] # m3 
    GWin=curLakeHydro$GWin[i-1] # m3
    GWout=curLakeHydro$GWout[i-1] # m3 
    directPrecip=curLakeHydro$DirectP[i-1] # m3
    lakeE=curLakeHydro$LakeE[i-1] # m3 
    SnowOnIceMelt=curLakeHydro$IceMelt[i-1] # m3 
    vol=curLakeHydro$Vol[i-1] # m3 
    area=curLakeHydro$Area[i-1] # m2
    area0=curLakeInitGeomorph$Area0 # m2
    perim=curLakeHydro$Perim[i-1] # m
    stage=curLakeHydro$Stage[i-1] # m 
    precipTemp=mean(curForce$AirTemp[i-1==nObs_hourly]) # degrees C ; temperature of preciptaiton (set to mean of daily air temp)
    press=mean(curForce$Pressure[i-1==nObs_hourly])*10 # atm pressure; converting from kPa to mb
    meanWind=mean(curForce$Wind[i-1==nObs_hourly]) # m s-1 
    snowDepth=curLakeHydro$snowDepth_m[i-1] # m   
    meanAirT=mean(curForce$AirTemp[i-1==nObs_hourly])
    
    # what to do if lake dries up??? i.e. when area ==0 ; put in a different check later on 
    if(area==0){
      out[i,1]=out[i-1,1]
      out[i,2]=out[i-1,2]
      out[i,3]=out[i-1,3]
      out[i,4]=out[i-1,4]
      out[i,5]=out[i-1,5]
      out[i,6]=out[i-1,6]
      out[i,7]=out[i-1,7]
      out[i,8]=out[i-1,8]
      out[i,9]=out[i-1,9]
      out[i,10]=out[i-1,10]
      out[i,11]=out[i-1,11]
      out[i-1,12]=0
      out[i-1,13]=0
      out[i-1,14]=0
      out[i,15]=out[i-1,15]
      out[i,16]=out[i-1,16]
      out[i-1,17]=0
      out[i,18]=out[i-1,18]
      out[i,19]=out[i-1,19]
      out[i-1,20]=0
      out[i,21]=out[i-1,21]
      out[i,22]=out[i-1,22]
      out[i-1,23]=0
      out[i-1,24]=0
      out[i-1,25]=0
      out[i-1,26]=0
      out[i-1,27]=0
      out[i-1,28]=0
      out[i-1,29]=0
      out[i-1,30]=0
      out[i-1,31]=0
      out[i-1,32]=0
      next
      # area=10
    }
    if(vol==0){
      vol=10
    }
    if(stage==0){
      stage=1
    }
    if(perim==0){
      perim=10
    }
    if(is.na(kD)|is.infinite(kD)){ # quick fix for kD failing 
      kD=1
    }
    if(kD<0){ # quick fix for kD failing 
      kD=0.1 
    }
    if(is.na(zmix)|is.infinite(zmix)){ # quick fix for zmix failing
      zmix=2
    }
    if(zmix>stage){
      zmix=stage
    }
    if(zmix<0.1){
      zmix=0.1
    }
    
    # atmospheric flux turned on or off? basing ice on/off on evap; 1 means flux on, 0 is off 
    fluxDummy<-as.numeric(curLakeHydro$LakeE[i-1]>0)
    
    # C biogeochemistry
    photoOx=0#	******* THIS IS TURNED OFF RIGHT NOW!!!!!! 44/1000/12 	# photooxidation rate constant, [mol c m-2 day-1]; Graneli et al. 1996, L&O, 41: 698-706
    floc=0#0.005		# fraction of DOC that floculates, [day-1]; von Wachenfeldt & Tranvik 2008 via Jones et al. 2012

    # epilimnion vs. hypolimion 
    r1=curLakeInitGeomorph$r2h*stage
    r2=curLakeInitGeomorph$r2h*(stage-zmix)
    Vepi=(1/3)*pi*(r1^2+r1*r2+r2^2)*zmix #truncated cone 
    Vhypo=(1/3)*pi*r2^2*(stage-zmix) # cone 
    if(Vhypo==0){
      Vhypo=0.000001 # making super small to avoid dividing by zero 
    }
    
    chl=phyto/Vepi*12*1000/phyto_C2Chl # current chlorophyll concentrations [mg m-3]  ug L-1
    kD=0.22*((DOCr_epi+DOCl_epi)/Vepi*12)+0.07*(chl/1000)-0.05 # Morris etal 1995
    # kD=0.22*(DOC/vol*12)+50*(chl/1000)-0.05 # Morris etal 1995; slight modifying chl coef to make stronger 
    
    # calculating todays zmix and volume to estimate entrainment 
    if(i<=10){ # checking to see if it's winter; modified on 2016-05-26 to mix 10 days prior to ice and 10 days after ice off
      mixisCheck<-curLakeHydro$LakeE[1:(i+10)]
    }else if(i>=(length(t)-10)){
      mixisCheck<-curLakeHydro$LakeE[(i-10):length(t)]
    }else{
      mixisCheck<-curLakeHydro$LakeE[(i-10):(i+10)]
    }
    if(any(mixisCheck==0)){ # checking to see if it's winter; modified on 2016-05-26 to mix 10 days prior to ice and 10 days after ice off
      zmixToday<-curLakeHydro$Stage[i] # fully mixed water column in winter or within 10 days of ice 
    }else{
      zmixToday<-10^(0.51263+(-0.65701*log10(kD))+(0.13717*log10(curLakeHydro$Area[i]/1000/1000))) # open ice zmix ; converting area in m2 to km2 
    }
    
    if(is.na(kD)|is.infinite(kD)){ # quick fix for kD failing 
      kD=0.5
    }
    if(kD<0){ # quick fix for kD failing 
      kD=0.1 
    }
    
    if(is.na(zmixToday)|is.infinite(zmixToday)){ # quick fix for zmix failing
      zmixToday=2
    }
    if(zmixToday>curLakeHydro$Stage[i]){
      zmixToday=curLakeHydro$Stage[i]
    }
    if(zmixToday<0.1){
      zmixToday=0.1
    }
    
    r1Today=curLakeInitGeomorph$r2h*curLakeHydro$Stage[i]
    r2Today=curLakeInitGeomorph$r2h*(curLakeHydro$Stage[i]-zmixToday)
    VepiToday=(1/3)*pi*(r1Today^2+r1Today*r2Today+r2Today^2)*zmixToday #truncated cone 
    VhypoToday=(1/3)*pi*r2Today^2*(curLakeHydro$Stage[i]-zmixToday) # cone 
    if(VhypoToday==0){
      VhypoToday=0.000001 # making super small to avoid dividing by zero 
    }
    
    
    dVhypo.dt=VhypoToday-Vhypo # use change in hypo volume for entrainment of water into hypo or epi 
    # volToday=curLakeHydro$Vol[i] # todays total lake volume 
    # dVol.dt=volToday-vol
    dVepi.dt=VepiToday-Vepi
    dVol.dt=(VhypoToday+VepiToday)-(Vhypo+Vepi) # making sure that total volume change is same as Vepi + Vhypo change 
    # entrain=dVepi.dt-(SWin+Baseflow)*(Vepi/vol)-directPrecip-SnowOnIceMelt-GWin*(Vepi/vol)+SWout*(Vepi/vol)+GWout*(Vepi/vol)+lakeE # entrainment
    entrain=dVepi.dt-dVol.dt
    # dZmix.dt=zmixToday-zmix # tells direction of entrainment 
    entrainEpi=as.numeric(entrain>0) # 1 for hypo entrain into epi, 0 otherwise 
    entrainHypo=as.numeric(entrain<0) # 1 for epi entrain into hypo, 0 otherwise
    entrain=abs(entrain) # absolute volume of water entrained; previous two line tell direction 
    
    hypoTempToday=7 # hypo is always 7 degrees C ; set to NTL hypo average 
    if(curLakeHydro$LakeE[i-1]>0){ # if open ice period, use simple energy budget for temperature of water 
      # use Lenters et al 2005 instead *********************************
      # see notebook for details of parameters 
  
      Lv=2.501e6-2370*epiTemp # latent heat of vaporization 
      
      A_p=rhoW*cW*(directPrecip/area/24/60/60)*(9.3-0) # W m-2 
      A_gwin=rhoW*cW*(GWin*(Vepi/vol)/area/24/60/60)*(9.3-0) # W m-2; 9.3 is mean GW temp from Lenters et al 2005
      A_swin=rhoW*cW*((Baseflow+SWin)*(Vepi/vol)/area/24/60/60)*(9.3-0) # W m-2; stream water temp same as gw temp?? 
      A_gwout=rhoW*cW*(GWout*(Vepi/vol)/area/24/60/60)*(epiTemp-0)
      A_swout=rhoW*cW*(SWout*(Vepi/vol)/area/24/60/60)*(epiTemp-0)
      A_evap=cW*(lakeE/area/24/60/60*rhoW/Lv)*(epiTemp-hypoTemp)/Lv
      A_net=A_p+A_gwin+A_swin-A_gwout-A_swout-A_evap
      
      # at equation 7 in Lenters 
      R_lwu=0.97*5.67e-8*(epiTemp+kelvinZero)^4 # mean longwave reflected 
      R_swd=mean(curForce$Shortwave[nObs_hourly==i-1]) # mean shortwave downward 
      R_lwd=mean(curForce$Longwave[nObs_hourly==i-1]) # mean longwave downward ; djusted to make colder for now *******************
      R_net=R_swd*(1-alpha_sw)+R_lwd*(1-alpha_lw)-R_lwu
      
      SVP_air=6.11*exp(17.27*meanAirT/(237.3+meanAirT)) # mb ; saturated vapor pressure of air temp; from HeatFluxAnalyzer 
      SVP_epi=6.11*exp(17.27*epiTemp/(237.3+epiTemp)) # mb ; saturated vapor pressure of epi temp
      VPD=mean(curForce$VPD[i-1==nObs_hourly])*10 # mb; vapor pressure deficit 
      RH=100-(VPD/SVP_air*100) # relative humidity in % 
      q_z=0.622*SVP_air/press # specific humidity; kg kg-1
      density=mean(curForce$Density[i-1==nObs_hourly]) # kg m-3 

      hours=length(curForce$AirTemp[i-1==nObs_hourly])
      # LakeMetabolizer 
      zeng=calc.zeng(dateTime = curForce$datetime[i-1==nObs_hourly],
                Ts = rep(epiTemp,hours),airT = curForce$AirTemp[i-1==nObs_hourly],Uz = rep(meanWind,hours),RH = RH,
                atm.press = press,wnd.z = 10,airT.z = 10,RH.z = 10)
      
      E=mean(zeng$alh,na.rm = T)
      H=mean(zeng$ash,na.rm = T)
      if(is.na(E)|is.na(H)){
        E=0
        H=0
      }

      Q_sed=0.24+5.3*cos(2*pi/365*(DOY(curLakeHydro$datetime[i-1])-358))
      
      # small lakes with large SWin gain a lot of heat with SWin and do not loose it until the next day when SWout kicks in
      if(abs(A_net)>abs(R_net)){
        A_sign=ifelse(A_net>0,1,-1)
        A_net=abs(R_net)*A_sign #  A_net can't be greater than the sum of rest of heat gains and losses 
      }
      
      S=R_net+Q_sed+A_net-(E+H) # W m-2 storage change 
      
      depiTemp.dt=S*area/rhoW/cW/out[i-1,18]*24*60*60
      if(zmix<2|(depiTemp.dt+out[i-1,21])>max(curForce$AirTemp)){
        depiTemp.dt=S*area/rhoW/cW/area/2*24*60*60 # forcing small lakes with small zmix to have heat gain across greater depth (2 meters )
      }
      
      if((depiTemp.dt+epiTemp)>max(curForce$AirTemp)){
        depiTemp.dt=max(curForce$AirTemp)-epiTemp
      }
      
      epiTempToday=epiTemp+depiTemp.dt

      if(entrainEpi==1){ # if hypo water is entrained, estimate heat transfer and new epi temp
        epiTempToday=(Vepi*rhoW*epiTempToday+entrainEpi*entrain*rhoW*hypoTemp)/(Vepi*rhoW+entrainEpi*entrain*rhoW)
      }
      if(is.na(epiTempToday)|is.infinite(epiTempToday)){
        epiTempToday=epiTemp
      }
      if(epiTempToday<2){ #quick fix for now  
        epiTempToday=2
      }
    }else{ # if ice cover, then set epiTemp and hypoTemp to 3 degrees C 
      epiTempToday=3
      hypoTempToday=3
    }
    
    DOCr_respired_epi=rFunc(rSlow,epiTemp)
    DOCl_respired_epi=rFunc(rFast,epiTemp)
    if((stage-zmix)==0){
      DOCr_respired_hypo=rFunc(0,hypoTemp) # no hypo volume so making DOC respired = 0 
      DOCl_respired_hypo=rFunc(0,hypoTemp) 
    }else{
      DOCr_respired_hypo=rFunc(rSlow,hypoTemp) 
      DOCl_respired_hypo=rFunc(rFast,hypoTemp)      
    }

    k600=k.vachon.base(U10,area)  #m day-1 from Vachon & Prairie 2013
    kCur=k600.2.kGAS.base(k600,epiTemp,gas = 'CO2')
    
    # atmCO2=curLakeHydro$co2[i-1] # [ppm]
    # atmEquilCO2=1*atmCO2/1e6/kH*1000	# concentration of CO2 in water at equilibrium with atmosphere, [mol C m-3]
    
    atmEquilCO2=gas_satconc(S=sal,t=epiTemp,P = press/1000,species = 'CO2',atm = curLakeHydro$co2[i-1]/1e6)/1000 # [mol C m-3]
    
    deposition=dailyTPOCdep(i-1,area,curForce,directPrecip,nObs_hourly)    
    if(DOY(curForce$datetime[i-1])%in%seq(275,290,1)){ # depositing leaves in autumn during a 2 week period 
      dailyLeafLoad=leafLoad/length(seq(275,290,1))*perim # leaf input in mol C day-1 
    }else{
      dailyLeafLoad=0
    }
    deposition=deposition+dailyLeafLoad
    
    d=5		#particle diameter [um]; Hanson et al. 2004
    # settling_epi=0.0188*(d/2)^2/zmix
    settling_epi=0.01 # settling rate in Hanson et al. 2004
    if((stage-zmix)==0){
      settling_hypo=0 # check if hypo is zero volume or not
    }else{
      # settling_hypo=0.0188*(d/2)^2/(stage-zmix)
      settling_hypo=0.01 # settling rate in Hanson et al. 2004 
    }

    # tells if sediment respiration should go into hypo or epi
    sedRespDummy<-as.numeric(zmix<stage) # 1 tells it should go into hypo, 0 into epi 
    
    GPP=integratedGPP(day=i-1,curForce=curForce,curForceDOY = nObs_hourly,chlCur = chl,SRP = P_epi,
                      V=Vepi,kD=kD, zmix=zmix,snowDepth=snowDepth)
    
    # carbonate speciation 
    # notes: alkalinity is ANC of filtered water; ANC is ANC of unfiltered water - some people assume ANC = Alkalinity 
    # ANC is in microEq/L; convert to mol / kg of solution by dividing by 1e6 multiplying by 1000 (mol/m3) and dividing by water density (kg/m3) to get to mol/kg
    # calculating speciations, fraction of DIC pool that is CO2, and pH by using water temperature and total alkalinity 
    aquaOut=aquaenv(S=sal,t = epiTemp,Pa = press/1000,d = 0,lat = curLakeLoc$Latitude,fCO2atm = curLakeHydro$co2[i-1]/1e6,SumCO2 = DIC_epi/Vepi/water.density(epiTemp,sal),
                    TA = alk/1e6*1000/water.density(epiTemp,sal))
    CO2=aquaOut$CO2[1]*water.density(epiTemp,sal)*Vepi
    HCO3=aquaOut$HCO3[1]*water.density(epiTemp,sal)*Vepi
    CO3=aquaOut$CO3[1]*water.density(epiTemp,sal)*Vepi
    fracCO2=CO2/DIC_epi
    pH=aquaOut$pH[1]
    
    
    ################################
    #### Differential equations ####
    ################################

    #algal biomass
    dphyto.dt=GPP*Vepi-GPP*Vepi*Rauto-settling_epi*phyto-phytoDeath*phyto-SWout*(Vepi/vol)*phyto/Vepi
    
    # SWin split proportionally to epi /hypo, SWout all from Epi; GW split proportionally to epi / hypo; entrainment of water to hypo / epi based on zmix deepening / shallowing 
    dP_epi.dt=(SWin+Baseflow)*(Vepi/vol)*streamP+directPrecip*precipP+SnowOnIceMelt*snowP+GWin*(Vepi/vol)*gwP+entrainEpi*entrain*P_hypo/Vhypo-
      SWout*(Vepi/vol)*P_epi/Vepi-GWout*(Vepi/vol)*P_epi/Vepi-GPP*Vepi*(1-Rauto)/phytoC2P+phytoDeath*phyto/phytoC2P*(1-phytoRecycling)-entrainHypo*entrain*P_epi/Vepi+
      perim*(Vepi/vol)*wLoad/12*streamP/streamDOC*curWatershed$percentWetland/100#*streamP/streamDOC # mol P  adding in wetland loading 2017-02-09; based on stoich of stream P:DOC; constant load rate for now 2017-03-22
    dP_hypo.dt=(SWin+Baseflow)*(Vhypo/vol)*streamP+GWin*(Vhypo/vol)*gwP+entrainHypo*entrain*P_epi/Vepi-
      GWout*(Vhypo/vol)*P_hypo/Vhypo-SWout*(Vhypo/vol)*P_hypo/Vhypo-entrainEpi*entrain*P_hypo/Vhypo+
      perim*(Vhypo/vol)*wLoad/12*streamP/streamDOC*curWatershed$percentWetland/100# mol P
    
    # reducing instability of DIC and emissions from discrete solver by setting limits on how much can be exported downstream after accounting for emissions / decay ;
    #   after accounting for emissions / DOC decay, DIC is exported based on previous concentration up to atm equilibrium and then we account for DIC taken away via GPP  
    dEmit.dt=fluxDummy*kCur*(atmEquilCO2-CO2/Vepi)*area # flux in or out of CO2 pool 
    
    dDIC_epi_biology.dt=photoOx*area+DOCr_respired_epi*DOCr_epi+DOCl_respired_epi*DOCl_epi-GPP*Vepi+GPP*Rauto*Vepi+settling_epi*tPOC_epi*(1-sedRespDummy)*(1-sedBurial_tPOC)+settling_hypo*tPOC_hypo*(1-sedRespDummy)*(1-sedBurial_tPOC)+settling_epi*phyto*(1-sedRespDummy)*(1-sedBurial_phyto)+
      +phytoDeath*phyto*(1-sedRespDummy)*(1-sedBurial_phyto)
    if((-dEmit.dt-dDIC_epi_biology.dt)>(DIC_epi-atmEquilCO2*Vepi)){
      dEmit.dt=-1*(DIC_epi-atmEquilCO2*Vepi-dDIC_epi_biology.dt)
    }
    
    dDIC_epi_hydro.dt=(SWin+Baseflow)*(Vepi/vol)*streamDIC+directPrecip*precipDIC+SnowOnIceMelt*snowDIC+GWin*(Vepi/vol)*gwDIC+entrainEpi*entrain*DIC_hypo/Vhypo-
      SWout*(Vepi/vol)*DIC_epi/Vepi-GWout*(Vepi/vol)*DIC_epi/Vepi-entrainHypo*entrain*DIC_epi/Vepi+
      perim*(Vepi/vol)*wLoad/12*streamDIC/streamDOC*curWatershed$percentWetland/100
    
    dDIC_epi.dt=dEmit.dt+dDIC_epi_biology.dt+dDIC_epi_hydro.dt
    
    if(-dDIC_epi.dt>(DIC_epi-atmEquilCO2*Vepi)){ # setting max that it can export in equilibrium with atm  
      dDIC_epi_hydro_gains.dt=(SWin+Baseflow)*(Vepi/vol)*streamDIC+directPrecip*precipDIC+SnowOnIceMelt*snowDIC+GWin*(Vepi/vol)*gwDIC+
        entrainEpi*entrain*DIC_hypo/Vhypo+perim*(Vepi/vol)*wLoad/12*streamDIC/streamDOC*curWatershed$percentWetland/100
      dDIC_epi_hydro_loss.dt=-1*(DIC_epi-atmEquilCO2*Vepi-dDIC_epi_biology.dt-dEmit.dt-dDIC_epi_hydro_gains.dt)
      dDIC_epi.dt=dDIC_epi_hydro_gains.dt+dDIC_epi_hydro_loss.dt+dDIC_epi_biology.dt+dEmit.dt
    }

    dDIC_hypo.dt=(SWin+Baseflow)*(Vhypo/vol)*streamDIC+GWin*(Vhypo/vol)*gwDIC+entrainHypo*entrain*DIC_epi/Vepi-
      GWout*(Vhypo/vol)*DIC_hypo/Vhypo-SWout*(Vhypo/vol)*DIC_hypo/Vhypo+DOCr_respired_hypo*DOCr_hypo+
      DOCl_respired_hypo*DOCl_hypo-entrainEpi*entrain*DIC_hypo/Vhypo+
      perim*(Vhypo/vol)*wLoad/12*streamDIC/streamDOC*curWatershed$percentWetland/100+
      settling_hypo*tPOC_hypo*(sedRespDummy)*(1-sedBurial_tPOC)+settling_epi*phyto*(sedRespDummy)*(1-sedBurial_phyto)+phytoDeath*phyto*(sedRespDummy)*(1-sedBurial_phyto)
    
    dDOCr_epi.dt=(SWin+Baseflow)*(Vepi/vol)*streamDOC*(1-fracLabile)+directPrecip*precipDOC*(1-fracLabile)+SnowOnIceMelt*snowDOC*(1-fracLabile)+GWin*(Vepi/vol)*gwDOC*(1-fracLabile)+
      entrainEpi*entrain*DOCr_hypo/Vhypo-
      SWout*(Vepi/vol)*DOCr_epi/Vepi-GWout*(Vepi/vol)*DOCr_epi/Vepi-photoOx*area-
      floc*DOCr_epi-DOCr_respired_epi*DOCr_epi+GPP*Vepi*GPPexudeR-entrainHypo*entrain*DOCr_epi/Vepi+
      perim*(Vepi/vol)*wLoad/12*(1-fracLabile)*curWatershed$percentWetland/100 #[mol C]   adding in wetland loading 2017-02-09; values from Hanson et al. 2014; constant for now 2017-03-22 
    dDOCr_hypo.dt=(SWin+Baseflow)*(Vhypo/vol)*streamDOC*(1-fracLabile)+GWin*(Vhypo/vol)*gwDOC*(1-fracLabile)+entrainHypo*entrain*DOCr_epi/Vepi-
      GWout*(Vhypo/vol)*DOCr_hypo/Vhypo-SWout*(Vhypo/vol)*DOCr_hypo/Vhypo-floc*DOCr_hypo-DOCr_respired_hypo*DOCr_hypo-entrainEpi*entrain*DOCr_hypo/Vhypo+
      perim*(Vhypo/vol)*wLoad/12*(1-fracLabile)*curWatershed$percentWetland/100#[mol C]
    
    dDOCl_epi.dt=(SWin+Baseflow)*(Vepi/vol)*streamDOC*fracLabile+directPrecip*precipDOC*fracLabile+SnowOnIceMelt*snowDOC*fracLabile+GWin*(Vepi/vol)*gwDOC*fracLabile+
      entrainEpi*entrain*DOCl_hypo/Vhypo-
      SWout*(Vepi/vol)*DOCl_epi/Vepi-GWout*(Vepi/vol)*DOCl_epi/Vepi-photoOx*area-
      floc*DOCl_epi-DOCl_respired_epi*DOCl_epi+GPP*Vepi*GPPexudeL-entrainHypo*entrain*DOCl_epi/Vepi+
      perim*(Vepi/vol)*wLoad/12*fracLabile*curWatershed$percentWetland/100 #[mol C]   adding in wetland loading 2017-02-09; values from Hanson et al. 2014; constant for now 2017-03-22
    dDOCl_hypo.dt=(SWin+Baseflow)*(Vhypo/vol)*streamDOC*fracLabile+GWin*(Vhypo/vol)*gwDOC*fracLabile+entrainHypo*entrain*DOCl_epi/Vepi-
      GWout*(Vhypo/vol)*DOCl_hypo/Vhypo-SWout*(Vhypo/vol)*DOCl_hypo/Vhypo-floc*DOCl_hypo-DOCl_respired_hypo*DOCl_hypo-entrainEpi*entrain*DOCl_hypo/Vhypo+
      perim*(Vhypo/vol)*wLoad/12*fracLabile*curWatershed$percentWetland/100#[mol C]
    
    dDOC_Load.dt=(SWin+Baseflow)*streamDOC+directPrecip*precipDOC+SnowOnIceMelt*snowDOC+GWin*gwDOC+
      perim*wLoad/12*curWatershed$percentWetland/100 #[mol C]
    
    dDOC_export.dt=GWout*(Vhypo/vol)*DOCl_hypo/Vhypo+SWout*(Vhypo/vol)*DOCl_hypo/Vhypo+
      GWout*(Vhypo/vol)*DOCr_hypo/Vhypo+SWout*(Vhypo/vol)*DOCr_hypo/Vhypo+
      SWout*(Vepi/vol)*DOCl_epi/Vepi+GWout*(Vepi/vol)*DOCl_epi/Vepi+
      SWout*(Vepi/vol)*DOCr_epi/Vepi+GWout*(Vepi/vol)*DOCr_epi/Vepi
      
    dDOC_Respired.dt=DOCr_respired_epi*DOCr_epi+DOCr_respired_hypo*DOCr_hypo+DOCl_respired_epi*DOCl_epi+DOCl_respired_hypo*DOCl_hypo
    
    dDOC_Respired_woExude.dt=DOCr_respired_epi*(DOCr_epi-GPP*Vepi*GPPexudeR)+
      DOCr_respired_hypo*DOCr_hypo+DOCl_respired_epi*(DOCl_epi-GPP*Vepi*GPPexudeL)+DOCl_respired_hypo*DOCl_hypo
    
    emergent_d_epi=DOCr_epi/(DOCr_epi+DOCl_epi)*DOCr_respired_epi+DOCl_epi/(DOCl_epi+DOCr_epi)*DOCl_respired_epi
    
    emergent_d_hypo=DOCr_hypo/(DOCr_hypo+DOCl_hypo)*DOCr_respired_hypo+DOCl_hypo/(DOCl_hypo+DOCr_hypo)*DOCl_respired_hypo
    
    dtPOC_epi.dt=(SWin+Baseflow)*(Vepi/vol)*streamPOC+deposition+floc*DOCr_epi+floc*DOCl_epi+entrainEpi*entrain*tPOC_hypo/Vhypo-
      SWout*(Vepi/vol)*tPOC_epi/Vepi-settling_epi*tPOC_epi-entrainHypo*entrain*tPOC_epi/Vepi  #[mol C]
    dtPOC_hypo.dt=(SWin+Baseflow)*(Vhypo/vol)*streamPOC+floc*DOCr_hypo+floc*DOCl_hypo+entrainHypo*entrain*tPOC_epi/Vepi-
      SWout*(Vhypo/vol)*tPOC_hypo/Vhypo-
      settling_hypo*tPOC_hypo+settling_epi*tPOC_epi*(sedRespDummy)-entrainEpi*entrain*tPOC_hypo/Vhypo  #[mol C]
    
    if(as.logical(sedRespDummy)){ # TRUE means sediment respiration should go into hypo (zmix < stage)
      dSed_tPOC.dt=settling_hypo*tPOC_hypo
      dBurial_tPOC.dt=settling_hypo*tPOC_hypo*sedBurial_tPOC
    }else{
      dSed_tPOC.dt=settling_epi*tPOC_epi
      dBurial_tPOC.dt=settling_epi*tPOC_epi*sedBurial_tPOC
    }
    
    dSed_phyto.dt=settling_epi*phyto+phytoDeath*phyto # or use mortality*phyto?? 
    dBurial_phyto.dt=(phytoDeath*phyto+settling_epi*phyto)*sedBurial_phyto # or use mortality*phyto?? 
    
    dEmit.dt=-dEmit.dt # positive emissions means out of lake 
    
    # quick fix so the state variables don't go negative
    if(-dphyto.dt>phyto|-dphyto.dt>0.9*phyto){
      dphyto.dt=-0.9*phyto
    }
    if(-dP_epi.dt>P_epi){
      dP_epi.dt=-0.9*P_epi
    }
    if(-dP_hypo.dt>P_hypo){
      dP_hypo.dt=-0.9*P_hypo
    }
    if(-dDIC_epi.dt>DIC_epi){
      dDIC_epi.dt=-0.9*DIC_epi
      dDIC_epi.dt=-(DIC_epi-atmEquilCO2*Vepi)+dDIC_epi_biology.dt
    }
    if(-dDIC_hypo.dt>DIC_hypo){
      dDIC_hypo.dt=-0.9*DIC_hypo
    }
    if(-dDOCr_epi.dt>DOCr_epi){
      dDOCr_epi.dt=-0.9*DOCr_epi
      dDOC_export.dt=0.9*(DOCr_epi+DOCr_hypo+DOCl_epi+DOCl_hypo)-dDOC_Respired.dt
    }
    if(-dDOCr_hypo.dt>DOCr_hypo){
      dDOCr_hypo.dt=-0.9*DOCr_hypo
    }
    if(-dDOCl_epi.dt>DOCl_epi){
      dDOCl_epi.dt=-0.9*DOCl_epi
    }
    if(-dDOCl_hypo.dt>DOCl_hypo){
      dDOCl_hypo.dt=-0.9*DOCl_hypo
    }
    if(-dtPOC_epi.dt>tPOC_epi){
      dtPOC_epi.dt=-0.9*tPOC_epi
    }
    if(-dtPOC_hypo.dt>tPOC_hypo){
      dtPOC_hypo.dt=-0.9*tPOC_hypo
    }
    
    ## update states 
    out[i,1]=out[i-1,1]+dphyto.dt
    if(((phyto+dphyto.dt)/Vepi*12*1000/phyto_C2Chl)<0.1){
      out[i,1]=0.1*Vepi/12/1000*phyto_C2Chl
    }
    out[i,2]=out[i-1,2]+dP_epi.dt
    out[i,3]=out[i-1,3]+dP_hypo.dt
    out[i,4]=out[i-1,4]+dDIC_epi.dt
    out[i,5]=out[i-1,5]+dDIC_hypo.dt
    out[i,6]=out[i-1,6]+dDOCr_epi.dt
    out[i,7]=out[i-1,7]+dDOCr_hypo.dt
    out[i,8]=out[i-1,8]+dDOCl_epi.dt
    out[i,9]=out[i-1,9]+dDOCl_hypo.dt
    out[i,10]=out[i-1,10]+dtPOC_epi.dt
    out[i,11]=out[i-1,11]+dtPOC_hypo.dt
    out[i-1,12]=dEmit.dt
    out[i-1,13]=dSed_tPOC.dt
    out[i-1,14]=dSed_phyto.dt
    out[i,15]=zmixToday
    out[i,16]=kD
    out[i-1,17]=GPP
    out[i,18]=VepiToday
    out[i,19]=VhypoToday
    out[i-1,20]=kCur
    out[i,21]=epiTempToday
    out[i,22]=hypoTempToday
    out[i-1,23]=dBurial_tPOC.dt
    out[i-1,24]=dBurial_phyto.dt
    out[i-1,25]=dDOC_Load.dt
    out[i-1,26]=dDOC_export.dt
    out[i-1,27]=dDOC_Respired.dt
    out[i-1,28]=dDOC_Respired_woExude.dt
    out[i-1,29]=emergent_d_epi
    out[i-1,30]=emergent_d_hypo
    out[i-1,31]=pH
    out[i-1,32]=fracCO2
    
    }
  return(out)
}
