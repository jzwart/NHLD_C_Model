# JAZ, 2017-02-06 regional water budget summary and results for NHLD modeling 


dir<-'E:/Jake/My Papers/NHLD Carbon Modeling/20170609/'

files<-list.files(dir)
fluxFiles=list.files(file.path(dir,'OUT_RESULTS'))
dates<-read.table(file.path(dir,'Date_OUT_ALL.txt'),stringsAsFactors = F,header = F,sep='\t')
colnames(dates)<-c('Year','Month','Day') #  dates for model run 
water_flux<-fluxFiles[grep('LakeFlux',fluxFiles)]
daily_prop<-fluxFiles[grep('DailyLakeProp',fluxFiles)]
init_prop<-fluxFiles[grep('IniLakeProp',fluxFiles)]
lakes<-unlist(strsplit(water_flux,'_DailyLakeFlux.txt'))
#checking for duplicates 
any(duplicated(lakes))

# forcing data 
forceDir<-'E:/Jake/My Papers/NHLD Carbon Modeling/20170609/Force/'
forceFiles<-list.files(forceDir)
dir<-'E:/Jake/My Papers/NHLD Carbon Modeling/20170609/OUT_RESULTS/'

# looping through lakes for summary stats
# HRT histogram 
hrt<-data.frame(lake=rep(NA,length(lakes)),HRT=rep(NA,length(lakes)),percent_Evap=rep(NA,length(lakes)),WALA=rep(NA,length(lakes)),LA=rep(NA,length(lakes)),
                CV_stage=rep(NA,length(lakes)),SD_stage=rep(NA,length(lakes)),Range_stage=rep(NA,length(lakes)),
                SWin=rep(NA,length(lakes)),Baseflow=rep(NA,length(lakes)),DirectP=rep(NA,length(lakes)),LakeE=rep(NA,length(lakes)),GWin=rep(NA,length(lakes)),
                GWout=rep(NA,length(lakes)),SWout=rep(NA,length(lakes)),IceMelt=rep(NA,length(lakes)),
                Stage=rep(NA,length(lakes)),Elev=rep(NA,length(lakes)),Vol=rep(NA,length(lakes)))
for(i in 1:length(lakes)){
  if(i%/%100%in%seq(1,length(lakes),by=100)){print(i)}
  curFlux<-read.table(file.path(dir,water_flux[lakes[i]==unlist(strsplit(water_flux,'_DailyLakeFlux.txt'))]),header = F,stringsAsFactors = F)
  colnames(curFlux)<-c('SWin','DirectP','LakeE','GWin','GWout','SWout','IceSnow','LandMelt','IceMelt','SWoutMinusLandMelt','LandMeltEst','Baseflow')
  curProp<-read.table(file.path(dir,daily_prop[lakes[i]==unlist(strsplit(daily_prop,'_DailyLakeProp.txt'))]),header = F,stringsAsFactors = F)
  colnames(curProp)<-c('Area','Radius','Perim','Stage','Elev','Vol')
  curInitGeomorph<-read.table(file.path(dir,init_prop[lakes[i]==unlist(strsplit(init_prop,'_IniLakeProp.txt'))]),header = F,stringsAsFactors = F)
  colnames(curInitGeomorph)<-c('Area0','Vol0','Radius0','Diameter0','Perim0','Stage0','DL','r2h','WA','WALA','Elev0_DEM','Vol_LinRes','Stage_LinRes','Stream_WA')
  
  # mean HRT, average volume divided by average of all inputs 
  curHRT<-mean(curProp$Vol)/mean(curFlux$SWin+curFlux$DirectP+curFlux$GWin+curFlux$IceMelt+curFlux$Baseflow)/365 # years 
  curPercentEvap<-mean(curFlux$LakeE)/mean(curFlux$LakeE+curFlux$GWout+curFlux$SWout) # percent export that is evap 
  curWALA<-curInitGeomorph$WALA
  curLA<-curInitGeomorph$Area0
  curCVstage<-sd(curProp$Stage)/mean(curProp$Stage)
  curSDstage<-sd(curProp$Stage)
  curRangeStage<-range(curProp$Stage)[2]-range(curProp$Stage)[1]
  curSWin=mean(curFlux$SWin)
  curBaseflow=mean(curFlux$Baseflow)
  curDirectP=mean(curFlux$DirectP)
  curLakeE=mean(curFlux$LakeE)
  curGWin=mean(curFlux$GWin)
  curGWout=mean(curFlux$GWout)
  curSWout=mean(curFlux$SWout)
  curIceMelt=mean(curFlux$IceMelt)
  curStage=mean(curProp$Stage)
  curElev=mean(curProp$Elev)
  curVol=mean(curProp$Vol)
  hrt[i,]<-data.frame(lakes[i],curHRT,curPercentEvap,curWALA,curLA,curCVstage,curSDstage,curRangeStage,
                      curSWin,curBaseflow,curDirectP,curLakeE,curGWin,curGWout,curSWout,curIceMelt,curStage,curElev,
                      curVol)
}


load('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/R Data/LakeWaterSummary_20170713.RData')

hrt$DR=(hrt$WALA*hrt$LA-hrt$LA)/hrt$Vol # correcting for wrong WALA where WA includes LA; makes sure to fix this if WALA changes 
hrt$percent_GWin=hrt$GWin/(hrt$GWin+hrt$DirectP+hrt$SWin+hrt$IceMelt+hrt$Baseflow)
hrt$percent_SWin=hrt$SWin/(hrt$GWin+hrt$DirectP+hrt$SWin+hrt$IceMelt+hrt$Baseflow)
hrt$percent_DirectP=hrt$DirectP/(hrt$GWin+hrt$DirectP+hrt$SWin+hrt$IceMelt+hrt$Baseflow)
hrt$percent_IceMelt=hrt$IceMelt/(hrt$GWin+hrt$DirectP+hrt$SWin+hrt$IceMelt+hrt$Baseflow)
hrt$percent_Baseflow=hrt$Baseflow/(hrt$GWin+hrt$DirectP+hrt$SWin+hrt$IceMelt+hrt$Baseflow)
hrt$percent_GWout=hrt$GWout/(hrt$GWout+hrt$LakeE+hrt$SWout)
hrt$percent_SWout=hrt$SWout/(hrt$GWout+hrt$LakeE+hrt$SWout)

cex.axis=2
cex.lab=2

windows()
png('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Figures/regionalHRThist.png',width = 7,
    height = 7,units='in',res = 300)
par(lwd=2,mar=c(5,5,5,5))
hist(hrt$HRT,ylab='Lakes (n)',xlab='HRT (years)',cex.lab=2,cex.axis=2,main='',lwd=2,xlim=c(0,8))
dev.off()

library(LSD)
# png('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Figures/regionalDR.png',width = 21,
#     height=7,units='in',res=300)
# par(lwd=2,mar=c(5,5,5,5))
# plot(hrt$HRT~log10(hrt$DR),pch=21,cex=1.5,col='black',bg=rgb(0,0,0,alpha=.2),cex.axis=cex.axis,cex.lab=cex.lab,
#      ylab='HRT (years)',xlab='WA:LakeVol',xaxt='n')
# axis(1,at = c(log10(0.10),log10(1),log10(10),log10(100),log10(1000)),labels = c(0.1,1,10,100,1000),cex.axis=cex.axis)

plot(hrt$percent_Evap~log10(hrt$DR),pch=21,cex=1.5,col='black',bg=rgb(0,0,0,alpha=.2),cex.axis=cex.axis,cex.lab=cex.lab,
     ylab='Fraction of Export as Evap',xlab='WA:LakeVol',xaxt='n')
axis(1,at = c(log10(0.10),log10(1),log10(10),log10(100),log10(1000)),labels = c(0.1,1,10,100,1000),cex.axis=cex.axis)

dev.off()

# alternative plots with kernal density 
png('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Figures/regionalDR.png',width = 14,
    height=14,units='in',res=300)
par(mar=c(5,6,5,5),mfrow=c(2,2))
cex.axis=3
cex.lab=3
cex=2
col=colorpalette(c('grey70','black'))
heatscatter(log10(hrt$DR),hrt$HRT,pch = 19,colpal = col,cex=cex,main='',cex.axis=cex.axis,cex.lab=cex.lab,
            ylab='HRT (years)',xlab='WA:LakeVol',xaxt='n')
axis(1,at = c(log10(0.10),log10(1),log10(10),log10(100),log10(1000)),labels = c(0.1,1,10,100,1000),cex.axis=cex.axis)

heatscatter(hrt$percent_Evap,hrt$HRT,pch = 19,colpal = col,cex=cex,main='',cex.axis=cex.axis,cex.lab=cex.lab,
            ylab='HRT (years)',xlab='Fraction of Export as Evap')

heatscatter(log10(hrt$DR),(hrt$percent_SWin+hrt$percent_Baseflow),pch = 19,colpal = col,cex=cex,main='',cex.axis=cex.axis,cex.lab=cex.lab,
            ylab='Fraction of Input as SW',xlab='WA:LakeVol',xaxt='n')
axis(1,at = c(log10(0.10),log10(1),log10(10),log10(100),log10(1000)),labels = c(0.1,1,10,100,1000),cex.axis=cex.axis)

heatscatter(hrt$percent_Evap,(hrt$percent_SWin+hrt$percent_Baseflow),pch = 19,colpal = col,cex=cex,main='',cex.axis=cex.axis,cex.lab=cex.lab,
            ylab='Fraction of Input as SW',xlab='Fraction of Export as Evap')
abline(h=0.5,lty=2,lwd=2)

# heatscatter(log10(hrt$DR),hrt$percent_Evap,pch = 19,colpal = col,cex=cex,main='',cex.axis=cex.axis,cex.lab=cex.lab,
#             ylab='Fraction of Input as SW',xlab='WA:LakeVol',xaxt='n')
# axis(1,at = c(log10(0.10),log10(1),log10(10),log10(100),log10(1000)),labels = c(0.1,1,10,100,1000),cex.axis=cex.axis)
# 
# heatscatter(log10(hrt$Vol),hrt$percent_Evap,pch = 19,colpal = col,cex=cex,main='',cex.axis=cex.axis,cex.lab=cex.lab,
#             ylab='Fraction of Export as Evap',xlab='Lake Volume (m3)',xaxt='n')
# axis(1,at = c(log10(10000),log10(1e6),log10(1e8),log10(1e10)),labels = c(10000,1e6,1e8,1e10),cex.axis=cex.axis)

dev.off()

png('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Figures/regionalHRT_vs_hydrMetrics.png',width = 21,
    height=7,units='in',res=300)
par(mar=c(7,6,5,5),mfrow=c(1,3))
cex.axis=3
cex.lab=3
cex=2.5
col=colorpalette(c('grey70','black'))
heatscatter(log10(hrt$LA),hrt$HRT,pch = 19,colpal = col,cex=cex,main='',cex.axis=cex.axis,cex.lab=cex.lab,
            ylab='HRT (years)',xlab='',xaxt='n')
axis(1,at = c(log10(100),log10(10000),log10(1000000),log10(10^8)),labels = c(100,10000,1000000,1e8),cex.axis=cex.axis)
title(xlab=expression(Lake~Area~(m^2)),cex.lab=cex.lab,line = 4.5)

heatscatter(log10(hrt$WALA),hrt$HRT,pch = 19,colpal = col,cex=cex,main='',cex.axis=cex.axis,cex.lab=cex.lab,
            ylab='',xlab='',xaxt='n')
axis(1,at = c(log10(1),log10(10),log10(100),log10(1000),log10(10000)),labels = c(1,10,100,1000,10000),cex.axis=cex.axis)
title(xlab='WA:LA',cex.lab=cex.lab,line = 4.5)

heatscatter(hrt$percent_Evap,hrt$HRT,pch = 19,colpal = col,cex=cex,main='',cex.axis=cex.axis,cex.lab=cex.lab,
            ylab='',xlab='')
title(xlab='Fraction of Export as Evap',cex.lab=cex.lab,line = 4.5)

dev.off()

png('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Figures/fracEvap_WALA.png',width = 7,
    height=7,units='in',res=300)
par(mar=c(6,6,5,5))
cex.axis=2.5
cex.lab=2.5
cex=2
col=colorpalette(c('grey70','black'))
heatscatter(log10(hrt$WALA),hrt$percent_Evap,pch = 19,colpal = col,cex=cex,main='',cex.axis=cex.axis,cex.lab=cex.lab,
            ylab='Fraction of Export as Evap',xlab='',xaxt='n')
axis(1,at = c(log10(1),log10(10),log10(100),log10(1000),log10(10000)),labels = c(1,10,100,1000,10000),cex.axis=cex.axis)
title(xlab='WA:LA',cex.lab=cex.lab,line = 4)
dev.off()


windows()
par(lwd=2,mar=c(5,5,5,5))
hist(hrt$HRT,ylab='Lakes (n)',xlab='HRT (years)',cex.lab=2,cex.axis=2,main='',lwd=2,xlim=c(0,10))

par(lwd=2,mar=c(5,5,5,5))
hist(hrt$HRT*log(hrt$Vol),ylab='Lakes (n)',xlab='HRT (years)',cex.lab=2,cex.axis=2,main='',lwd=2)

# water volume - weighted HRT 
hist(hrt$HRT*log(hrt$Vol))

hrt$waterIn=hrt$GWin+hrt$SWin+hrt$DirectP+hrt$IceMelt
log(hrt$waterIn)*hrt$HRT
allLakesSum<-data.frame(GWin=sum(hrt$GWin),SWin=sum(hrt$SWin),DirectP=sum(hrt$DirectP),IceMelt=sum(hrt$IceMelt),
                        GWout=-sum(hrt$GWout),SWout=-sum(hrt$SWout),Evap=-sum(hrt$LakeE))
sum(allLakesSum$GWin,allLakesSum$SWin,allLakesSum$DirectP,allLakesSum$IceMelt)
sum(allLakesSum$GWout,allLakesSum$SWout,allLakesSum$Evap)

sum(hrt$Vol)/sum(allLakesSum$GWin,allLakesSum$SWin,allLakesSum$DirectP,allLakesSum$IceMelt)/365

windows()
png('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Figures/regionalHRThist.png',width = 14,
    height = 7,units='in',res = 300)
cex.axis=2
cex.lab=2
par(lwd=2,mar=c(5,5,5,5),mfrow=c(1,2))
hist(hrt$HRT,ylab='Lakes (n)',xlab='HRT (years)',cex.lab=2,cex.axis=2,main='',lwd=2,xlim=c(0,10),breaks= seq(0,10,by = 0.5))
abline(v=sum(hrt$Vol)/sum(allLakesSum$GWin,allLakesSum$SWin,allLakesSum$DirectP,allLakesSum$IceMelt)/365,
       lty=2,lwd=2)
par(lwd=1, mar=c(5,5,5,5))
heatscatter(log10(hrt$Vol),hrt$HRT,pch = 19,colpal = col,cex=cex,main='',cex.axis=cex.axis,cex.lab=cex.lab,
            ylab='HRT (years)',xlab=expression(Lake~Volume~(m^3)),xaxt='n')
axis(1,at = c(log10(1000),log10(100000),log10(1e7),log10(1e9),log10(1e10)),labels = c(1000,100000,1e7,1e9,1e10),cex.axis=cex.axis)
dev.off()

png('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Figures/regionalHRThist_onlyHist.png',width = 7,
    height = 7,units='in',res = 300)
cex.axis=2
cex.lab=2
par(lwd=2,mar=c(5,5,5,5))
hist(hrt$HRT,ylab='Lakes (n)',xlab='HRT (years)',cex.lab=2,cex.axis=2,main='',lwd=2,xlim=c(0,10),breaks= seq(0,10,by = 0.5))
abline(v=sum(hrt$Vol)/sum(allLakesSum$GWin,allLakesSum$SWin,allLakesSum$DirectP,allLakesSum$IceMelt)/365,
       lty=2,lwd=2)
dev.off()


heatscatter(hrt$percent_SWin,hrt$percent_Evap,pch = 19,colpal = col,cex=cex,main='',cex.axis=cex.axis,cex.lab=cex.lab,
            ylab='HRT',xlab='Fraction of Input as SW') # shows drainage seepage split and implications for C retained 
abline(h=0.5,lty=2,lwd=2)
abline(v=0.5,lty=2,lwd=2)

heatscatter(hrt$percent_SWin,hrt$percent_DirectP,pch = 19,colpal = col,cex=cex,main='',cex.axis=cex.axis,cex.lab=cex.lab,
            ylab='HRT',xlab='Fraction of Input as SW')
abline(h=0.5,lty=2,lwd=2)
abline(v=0.5,lty=2,lwd=2)

heatscatter(hrt$percent_SWin,hrt$percent_SWout,pch = 19,colpal = col,cex=cex,main='',cex.axis=cex.axis,cex.lab=cex.lab,
            ylab='HRT',xlab='Fraction of Input as SW')
abline(h=0.5,lty=2,lwd=2)
abline(v=0.5,lty=2,lwd=2)

heatscatter(hrt$percent_Evap,hrt$HRT,pch = 19,colpal = col,cex=cex,main='',cex.axis=cex.axis,cex.lab=cex.lab,
            ylab='HRT',xlab='Fraction of Input as SW')
abline(h=0.5,lty=2,lwd=2)
abline(v=0.5,lty=2,lwd=2)

heatscatter(hrt$percent_Evap,hrt$percent_SWin,pch = 19,colpal = col,cex=cex,main='',cex.axis=cex.axis,cex.lab=cex.lab,
            ylab='HRT',xlab='Fraction of Input as SW')
abline(h=0.5,lty=2,lwd=2)
abline(v=0.5,lty=2,lwd=2)

heatscatter(hrt$percent_Evap,hrt$Elev,pch = 19,colpal = col,cex=cex,main='',cex.axis=cex.axis,cex.lab=cex.lab,
            ylab='HRT',xlab='Fraction of Input as SW')
abline(h=0.5,lty=2,lwd=2)
abline(v=0.5,lty=2,lwd=2)

heatscatter(log10(hrt$DR),log10(hrt$LA),pch = 19,colpal = col,cex=cex,main='',cex.axis=cex.axis,cex.lab=cex.lab,
            ylab=expression(Lake~Area~(m^2)),xlab='WA:LakeVol',xaxt='n',yaxt='n')
axis(1,at = c(log10(0.10),log10(1),log10(10),log10(100),log10(1000)),labels = c(0.1,1,10,100,1000),cex.axis=cex.axis)
axis(2,at = c(log10(1000),log10(1e5),log10(1e7)),labels = c(1000,1e5,1e7),cex.axis=cex.axis)


windows()
par(mar=c(5,5,5,5))
plot(hrt$percent_Evap~hrt$HRT,pch=21,cex=1.5,col='black',bg=rgb(0,0,0,alpha=.2),cex.axis=cex.axis,cex.lab=cex.lab,
     ylab='Fraction of Export as Evap',xlab='HRT (years)')

windows()
par(mar=c(5,5,5,5))
plot(hrt$percent_Evap~log10(hrt$WALA),pch=21,cex=1.5,col='black',bg=rgb(0,0,0,alpha=.2),cex.axis=cex.axis,cex.lab=cex.lab,
     ylab='Fraction of Export as Evap',xlab='WALA')

windows()
plot(hrt$HRT~log10(hrt$WALA),pch=21,cex=1.5,col='black',bg=rgb(0,0,0,alpha=.2),cex.axis=cex.axis,cex.lab=cex.lab,
     ylab='HRT (years)',xlab='WA:LA',xaxt='n')
axis(1,at = c(log10(1),log10(10),log10(100),log10(1000)),labels = c(1,10,100,1000),cex.axis=cex.axis)

windows()
plot(hrt$HRT~log10(hrt$LA),pch=21,cex=1.5,col='black',bg=rgb(0,0,0,alpha=.2),cex.axis=cex.axis,cex.lab=cex.lab,
     ylab='HRT (years)',xlab='Lake Area',xaxt='n')
axis(1,at = c(log10(100),log10(10000),log10(100000),log10(1000000)),labels = c(100,10000,100000,1000000),cex.axis=cex.axis)

# windows() 
# hist(log10(hrt$WALA))

windows()
plot(hrt$HRT~log10(hrt$DR),pch=21,cex=1.5,col='black',bg=rgb(0,0,0,alpha=.2),cex.axis=cex.axis,cex.lab=cex.lab,
     ylab='HRT (years)',xlab='WA:LakeVol',xaxt='n')
axis(1,at = c(log10(0.10),log10(1),log10(10),log10(100),log10(1000)),labels = c(0.1,1,10,100,1000),cex.axis=cex.axis)

windows()
plot(hrt$percent_Evap~log10(hrt$DR),pch=21,cex=1.5,col='black',bg=rgb(0,0,0,alpha=.2),cex.axis=cex.axis,cex.lab=cex.lab,
     ylab='Fraction of Export as Evap',xlab='WA:LakeVol',xaxt='n')
axis(1,at = c(log10(0.10),log10(1),log10(10),log10(100),log10(1000)),labels = c(0.1,1,10,100,1000),cex.axis=cex.axis)

windows()
plot(log10(hrt$LA)~log10(hrt$WALA),pch=21,cex=1.5,col='black',bg=rgb(0,0,0,alpha=.2),cex.axis=cex.axis,cex.lab=cex.lab,
     ylab='Lake Area (m2)',xlab='WA:LA',xaxt='n',yaxt='n')
axis(1,at = c(log10(0.10),log10(1),log10(10),log10(100),log10(1000)),labels = c(0.1,1,10,100,1000),cex.axis=cex.axis)
axis(2,at = c(log10(100),log10(10000),log10(100000),log10(1000000)),labels = c(100,10000,100000,1000000),cex.axis=cex.axis)

windows()
plot(log10(hrt$Vol)~log10(hrt$LA),pch=21,cex=1.5,col='black',bg=rgb(0,0,0,alpha=.2),cex.axis=cex.axis,cex.lab=cex.lab,
     ylab='Lake Vol (m3)',xlab='Lake Area (m2)',xaxt='n',yaxt='n')
axis(1,at = c(log10(100),log10(10000),log10(100000),log10(1000000)),labels = c(100,10000,100000,1000000),cex.axis=cex.axis)
axis(2,at = c(log10(100),log10(10000),log10(100000),log10(1000000)),labels = c(100,10000,100000,1000000),cex.axis=cex.axis)



windows()
par(lwd=2)
hist(hrt$percent_GWin,ylab='Lakes (n)',xlab='Fraction of Input as GW',cex.lab=2,cex.axis=2,main='',lwd=2,xlim=c(0,1))

windows()
par(lwd=2)
hist(hrt$percent_SWin,ylab='Lakes (n)',xlab='Fraction of Input as SW',cex.lab=2,cex.axis=2,main='',lwd=2,xlim=c(0,1))

windows()
par(lwd=2)
hist(hrt$percent_DirectP,ylab='Lakes (n)',xlab='Fraction of Input as Precip',cex.lab=2,cex.axis=2,main='',lwd=2,xlim=c(0,1))

windows()
par(lwd=2)
hist(hrt$percent_IceMelt,ylab='Lakes (n)',xlab='Fraction of Input as IceMelt',cex.lab=2,cex.axis=2,main='',lwd=2,xlim=c(0,1))

windows()
par(lwd=2)
hist(hrt$percent_Evap,ylab='Lakes (n)',xlab='Fraction of Export as Evap',cex.lab=2,cex.axis=2,main='',lwd=2,xlim=c(0,1))

windows()
par(lwd=2)
hist(hrt$percent_GWout,ylab='Lakes (n)',xlab='Fraction of Export as GW',cex.lab=2,cex.axis=2,main='',lwd=2,xlim=c(0,1))

windows()
par(lwd=2)
hist(hrt$percent_SWout,ylab='Lakes (n)',xlab='Fraction of Export as SW',cex.lab=2,cex.axis=2,main='',lwd=2,xlim=c(0,1))

allLakesSum<-data.frame(t(data.frame(GWin=sum(hrt$GWin),SWin=sum(hrt$SWin+hrt$Baseflow),DirectP=sum(hrt$DirectP),IceMelt=sum(hrt$IceMelt),
                        GWout=-sum(hrt$GWout),SWout=-sum(hrt$SWout),Evap=-sum(hrt$LakeE))))
colnames(allLakesSum)<-'total'

library(ggplot2)

p <- ggplot(allLakesSum)
p+geom_bar(aes(rownames(allLakesSum),y=total),stat = 'identity')+scale_x_discrete(limits=rownames(allLakesSum))+theme_classic()+
  labs(x='',y='Total Flux (m-3)')+theme(legend.position="none",axis.text=element_text(size=16),
                                        axis.title=element_text(size=16,face="bold"),
                                        axis.line = element_line(colour = "black",size = 1, linetype = "solid"))

allLakesSum<-data.frame(GWin=sum(hrt$GWin),SWin=sum(hrt$SWin),DirectP=sum(hrt$DirectP),IceMelt=sum(hrt$IceMelt),
                                     GWout=-sum(hrt$GWout),SWout=-sum(hrt$SWout),Evap=-sum(hrt$LakeE))
sum(allLakesSum$GWin,allLakesSum$SWin,allLakesSum$DirectP,allLakesSum$IceMelt,allLakesSum$Baseflow)
sum(allLakesSum$GWout,allLakesSum$SWout,allLakesSum$Evap)

sum(hrt$Vol)/sum(allLakesSum$GWin,allLakesSum$SWin,allLakesSum$DirectP,allLakesSum$IceMelt)/365


windows()
plot(hrt$CV_stage~log10(hrt$WALA))

windows()
plot(hrt$CV_stage~hrt$HRT)

windows()
plot(hrt$SD_stage~hrt$HRT)

windows()
plot(hrt$Range_stage~hrt$HRT)

windows() 
plot(hrt$CV_stage~hrt$percent_Evap)


