# NHLD model summary; JAZ 2017-03-26 

rm(list=ls())
library(LSD)
library(maptools)
library(Hmisc)

dir<-'D:/Jake/My Papers/NHLD Carbon Model/Results/HighWload/20170829_HighWLoad/'
files=list.files(dir)
# files=files[780:801]
skip=6*365 # days to skip 

sum<-data.frame() # only open ice  
all<-data.frame()
val<-data.frame()
for(i in 1:length(files)){
  cur<-read.table(file.path(dir,files[i]),stringsAsFactors = F,sep='\t',header=T)
  lake<-strsplit(files[i],'_C_model.txt',fixed = T)[[1]]
  cur<-na.omit(cur)
  if(length(cur$time)<2){
    next
  }
  curVal=cur[skip:nrow(cur),]
  curVal=curVal[curVal$LakeE>0&as.Date(curVal$datetime)>as.Date('2004-01-01')&as.Date(curVal$datetime)<as.Date('2004-12-31'),]
  curVal=curVal[,3:ncol(curVal)]
  cur<-cur[skip:(nrow(cur)-0),3:ncol(cur)] # skipping first X number of days 
  cur<-cur[cur$Vol>0,] # only keeping days when there's actually water 
  curSum=cur[cur$LakeE>0,] # only open ice periods
 
  cur<-data.frame(t(apply(cur,MARGIN = 2,FUN = mean)))
  curSum<-data.frame(t(apply(curSum,MARGIN = 2,FUN = mean)))
  curVal<-data.frame(t(apply(curVal,MARGIN = 2,FUN=mean)))
  cur$Permanent_<-lake
  curSum$Permanent_<-lake
  curVal$Permanent_<-lake
  
  sum<-rbind(sum,curSum)
  all<-rbind(all,cur)
  val<-rbind(val,curVal)
}

sum2=sum

# 
load('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/R Data/LakeCsummary_20170828_highWload.RData')

watersheds<-read.table('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Data/C model forcing data/NHLDsheds_20170323.txt',
                       stringsAsFactors = F,header=T,sep = '\t')

sum=merge(sum,watersheds,by='Permanent_')
all=merge(all,watersheds,by='Permanent_')

all$alk<-ifelse(all$GWin+all$Baseflow>0,10^(3.2815*(1-exp(-1.2765*log10(all$GWin+all$Baseflow+all$SWin)/3.2815))),0)
all<-all[sort.list(all$Permanent_),]
sum<-sum[sort.list(sum$Permanent_),]
sum$alk<-all$alk
# sum$streamDOC=exp(1.3961+3.245*(sum$percentWetland))/12*1000/1000	# from lottig 2012; mol m-3
# sum$gwDOC=13/12*1000/1000#median(ntlGW$doc,na.rm=TRUE)/12*1000/1000	# mol m-3
# sum$precipDOC=3.195/12*1000/1000		#mol m-3; from UNDERC rain data; also see Likens et al. 1983 @ Hubbard Brook 1.1 mg L-1 
# sum$snowDOC=0.6/12*1000/1000 # snow DOC from UNDERC snow; 0.5985 mg C L-1 to mol C m-3; snow DOC from Sebestyen lab

# sum$DOCLoad=sum$SWin*sum$streamDOC+sum$Baseflow*sum$streamDOC+sum$GWin*sum$gwDOC+sum$DirectP*sum$precipDOC # mol C day-1
# sum$DOCLoad=sum$DOCLoad*12*1000 # mg C day-1 
# sum$emergent_d20=sum$DOCl_epi/(sum$DOCl_epi+sum$DOCr_epi)*0.08+sum$DOCr_epi/(sum$DOCl_epi+sum$DOCr_epi)*0.002
sum$HRT<-sum$Vol/(sum$GWin+sum$SWin+sum$DirectP+sum$Baseflow+sum$IceMelt)
sum$HRT_woEvap=sum$Vol/(sum$GWout+sum$SWout)
all$HRT<-all$Vol/(all$GWin+all$SWin+all$DirectP+all$Baseflow+all$IceMelt)
all$HRT_woEvap=all$Vol/(all$GWout+all$SWout)
# sum$DOCRespired=sum$emergent_d20/(sum$emergent_d20+(1/sum$HRT)) # Brett's way of calculating Frac Retained 
# sum$DOCRespired_Evap=sum$emergent_d20/(sum$emergent_d20+(1/sum$HRT_woEvap))
# sum$FracRet=sum$DOC_Respired/sum$DOC_Load # DOC respired includes phyto exude 
# sum$DOC_Respired_woExude=sum$DOC_Respired-(sum$GPP*sum$Vepi*0.07*0.08)-sum$GPP*sum$Vepi*0.07*0.002
# sum$FracRet_woExude=sum$DOC_Respired_woExude/sum$DOC_Load
sum$FracRet=1-(sum$DOC_export/sum$DOC_Load)
all$FracRet=1-(all$DOC_export/all$DOC_Load)

sum$DIC_load=sum$GWin*0.7025+sum$SWin*0.9041667+sum$Baseflow*0.9041667+sum$DirectP*0.0833333+sum$IceMelt*0.01360082 # mol C day-1 
all$DIC_load=all$GWin*0.7025+all$SWin*0.9041667+all$Baseflow*0.9041667+all$DirectP*0.0833333+all$IceMelt*0.01360082 # mol C day-1 

# fraction of DOC to DIC load 
all$DOC_to_DIC_load=all$DOC_Load/all$DIC_load
hist(all$DOC_to_DIC_load)
summary(all$DOC_to_DIC_load)

sum(all$DOC_Load)/sum(all$DIC_load)
sum(all$DOC_export)/sum(all$DIC_export)

sum$sed_resp<-sum$Sed_phyto*0.75+sum$Sed_tPOC*0.1 # 75% of phyto C and 10% of tpoc C was converted to co2 
all$sed_resp<-all$Sed_phyto*0.75+all$Sed_tPOC*0.1 # 75% of phyto C and 10% of tpoc C was converted to co2 

sum$waterIn=sum$GWin+sum$SWin+sum$DirectP+sum$Baseflow+sum$IceMelt
all$waterIn=all$GWin+all$SWin+all$DirectP+all$Baseflow+all$IceMelt

sum$fluvialOut=sum$GWout+sum$SWout # m3 day-1
all$fluvialOut=all$GWout+all$SWout 

# DIC export (update to calculate in model, this is approximated for right now)
sum$DIC_export=(sum$DIC_epi+sum$DIC_hypo)/sum$Vol*sum$fluvialOut # mol C day-1
all$DIC_export=(all$DIC_epi+all$DIC_hypo)/all$Vol*all$fluvialOut # mol C day-1 

# fraction of DOC to DIC fluvial export
all$DOC_to_DIC_export=all$DOC_export/all$DIC_export
summary(all$DOC_to_DIC_export)

a=all$DOC_to_DIC_load/all$DOC_to_DIC_export
summary(a)

sum$dicLoadvResp=sum$DIC_load/(sum$DOC_Respired+sum$sed_resp)
all$dicLoadvResp=all$DIC_load/(all$DOC_Respired+all$sed_resp)

# respiration rates etc.. based on lake size bins from Downing 
lakeSizeBins=c(0,0.01,.1,1,10,100)*1e6 # breaks for max cutoff of lake size from Downing et al. 2006
sum$lakeSizeBins=cut(sum$Area,breaks = lakeSizeBins)
all$lakeSizeBins=cut(all$Area,breaks=lakeSizeBins)

# about a dozen lakes have negative fraction retained, cutting those out 
all<-all[all$FracRet>=0,]
sum<-sum[sum$Permanent_%in%all$Permanent_,]

windows()
plot(sum$FracRet~sum$HRT,ylim=c(0,1))
# points(sum$FracRet_woExude~sum$HRT,col='red')
plot(sum$FracRet~log10(sum$HRT),ylim=c(0,1))
plot(all$FracRet~all$HRT,ylim=c(0,1))

windows()
plot(log10(sum$DOC_Load)~sum$HRT)

# sum$DOC_conc=(sum$DOCr_epi+sum$DOCl_epi+sum$DOCl_hypo+sum$DOCr_hypo)/(sum$Vepi+sum$Vhypo)*12

litFracRet<-read.csv('/Users/Jake/Documents/Jake/Conferences/2015/Gordon Research Conference/GRC/DrainageVseepage_Cretention.csv',
                     stringsAsFactor=F)
sum$percentEvap<-sum$LakeE/(sum$LakeE+sum$GWout+sum$SWout)
all$percentEvap<-all$LakeE/(all$LakeE+all$GWout+all$SWout)

windows()
png('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Figures/Fig3_fracRet_HRT.png',width = 7,
    height=7,units = 'in',res = 300)
par(mar=c(6,6,5,2))
cex.axis=2.5
cex.lab=2.5
cex=2
col=colorpalette(c('grey70','black'))
# heatscatter(log10(sum$HRT),sum$FracRet,pch = 19,colpal = col,cex=cex,main='',cex.axis=cex.axis,cex.lab=cex.lab,ylim=c(0,1),
            # ylab='Fraction Retained',xlab='HRT',xaxt='n')
heatscatter(log10(all$HRT),all$FracRet,pch = 19,colpal = col,cex=cex,main='',cex.axis=cex.axis,cex.lab=cex.lab,ylim=c(0,1),
            ylab='Fraction Retained',xlab='HRT (days)',xaxt='n')
axis(1,at = c(log10(0.10),log10(1),log10(10),log10(100),log10(1000)),labels = c(0.1,1,10,100,1000),cex.axis=cex.axis)
points(log10(litFracRet$HRT_yrs[litFracRet$LakeType=='Drainage']*365),
       litFracRet$C_ret_proportion[litFracRet$LakeType=='Drainage'],pch=19,col='blue',cex=cex)
points(log10(litFracRet$HRT_yrs[litFracRet$LakeType=='Seepage'|litFracRet$LakeType=='UNDERC_Long']*365),
       litFracRet$C_ret_proportion[litFracRet$LakeType=='Seepage'|litFracRet$LakeType=='UNDERC_Long'],pch=19,col='red',cex=cex)
legend('topleft',legend = c('Modeled','Literature Drainage','Literature Seepage'),
      col=c('grey 40','blue','red'),bty = 'n',pch = 19,cex = 1.5,pt.cex = 2)
dev.off()

png('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Figures/Fig3_DOCrespired_HRT.png',width = 7,
    height=7,units = 'in',res = 300)
par(mar=c(6,6,5,2))
cex.axis=2.5
cex.lab=2.5
cex=2
col=colorpalette(c('grey70','black'))
heatscatter(log10(all$HRT),log10(all$DOC_Respired/all$Area*12*365),pch = 19,colpal = col,cex=cex,main='',cex.axis=cex.axis,cex.lab=cex.lab,xaxt='n',
            xlab='HRT (days)',ylab='DOC Respired (g C m-2 yr-1)',yaxt='n')
axis(2,at = c(log10(0.10),log10(1),log10(10),log10(100),log10(1000)),labels = c(0.1,1,10,100,1000),cex.axis=cex.axis)
axis(1,at = c(log10(0.10),log10(1),log10(10),log10(100),log10(1000)),labels = c(0.1,1,10,100,1000),cex.axis=cex.axis)
dev.off()


# windows()
# par(mar=c(6,6,5,2))
# cex.axis=3
# cex.lab=3
# cex=2
# col=colorpalette(c('grey70','black'))
# heatscatter(sum$percentEvap,sum$FracRet,pch = 19,colpal = col,cex=cex,main='',cex.axis=cex.axis,cex.lab=cex.lab,ylim=c(0,1),
#             ylab='Fraction Retained',xlab='Fraction Evap')
# 
# windows()
# par(mar=c(6,6,5,2))
# cex.axis=3
# cex.lab=3
# cex=2
# col=colorpalette(c('grey70','black'))
# heatscatter(sum$percentEvap,log10(sum$DOC_Respired/sum$Area*12*1000),pch = 19,colpal = col,cex=cex,main='',cex.axis=cex.axis,cex.lab=cex.lab,
#             ylab='Fraction Retained',xlab='Fraction Evap')


windows()
png('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Figures/Fig3_d_HRT.png',width = 7,
    height=7,units = 'in',res = 300)
par(mar=c(6,6,5,2))
cex.axis=2.5
cex.lab=2.5
cex=2
col=colorpalette(c('grey70','black'))
heatscatter(log10(all$HRT),all$emergent_d_epi,pch = 19,colpal = col,cex=cex,main='',cex.axis=cex.axis,cex.lab=cex.lab,
            ylab='Turnover Rate (day-1)',xlab='HRT (days)',xaxt='n')
axis(1,at = c(log10(0.10),log10(1),log10(10),log10(100),log10(1000)),labels = c(0.1,1,10,100,1000),cex.axis=cex.axis)
dev.off()

png('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Figures/Fig3_d_fracEvap.png',width = 7,
    height=7,units = 'in',res = 300)
par(mar=c(6,6,5,2))
cex.axis=2.5
cex.lab=2.5
cex=2
col=colorpalette(c('grey70','black'))
heatscatter(log10(all$percentEvap),all$emergent_d_epi,pch = 19,colpal = col,cex=cex,main='',cex.axis=cex.axis,cex.lab=cex.lab,
            ylab='Turnover Rate (day-1)',xlab='Fraction Export as Evap')
dev.off()
summary(lm(all$emergent_d_epi~log10(all$percentEvap)+all$meanDepth))
res=resid(lm(all$emergent_d_epi~log10(all$percentEvap)))
plot(res~all$meanDepth)
summary(lm(res~all$meanDepth))

png('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Figures/Fig3_d_HRT_logged.png',width = 7,
    height=7,units = 'in',res = 300)
par(mar=c(6,6,5,2))
cex.axis=2.5
cex.lab=2.5
cex=2
col=colorpalette(c('grey70','black'))
heatscatter(log(all$HRT/365),log(all$emergent_d_epi*365),pch = 19,colpal = col,cex=cex,main='',cex.axis=cex.axis,cex.lab=cex.lab,
            ylab='log[Turnover Rate (yr-1)]',xlab='log[HRT (years)]',xaxt='n',yaxt='n')
axis(1,at = c(log(.01),log(0.10),log(1),log(5)),labels = c(.01,0.1,1,5),cex.axis=cex.axis)
axis(2,at = c(log(.01),log(0.50),log(1),log(2)),labels = c(.01,0.5,1,2),cex.axis=cex.axis)
abline(lm(log(all$emergent_d_epi*365)~log(all$HRT/365)))
summary(lm(log(all$emergent_d_epi[all$HRT>10]*365)~log(all$HRT[all$HRT>10]/365)))
abline(lm(log(all$emergent_d_epi[all$HRT>10]*365)~log(all$HRT[all$HRT>10]/365)))
dev.off()

heatscatter(log(all$HRT_woEvap/365),log(all$emergent_d_epi*365),pch = 19,colpal = col,cex=cex,main='',cex.axis=cex.axis,cex.lab=cex.lab,
            ylab='log[Turnover Rate (yr-1)]',xlab='log[HRT (years)]',xaxt='n',yaxt='n')
axis(1,at = c(log(.01),log(0.10),log(1),log(5)),labels = c(.01,0.1,1,5),cex.axis=cex.axis)
axis(2,at = c(log(.01),log(0.50),log(1),log(2)),labels = c(.01,0.5,1,2),cex.axis=cex.axis)
summary(lm(log(all$emergent_d_epi*365)~log(all$HRT_woEvap/365)))

dillon=data.frame(meanDepth=c(8.5,8.9,9.2,5,13.3,7.9,14.2),runoff=c(1.5,4.23,5.6,2.66,4.16,2.00,5.44),
                  WALA=c(105.9/52.35,271.8/34.41,521.8/56.74,406.4/93.60,470.7/71.38,95.5/32.14,532.4/57.13),
                  ret=c(.59,.42,.37,.55,.42,.69,.37),settle=c(2.2,3,3.3,3.2,2.9,4.6,3.2))
dillon$modeledRet=(dillon$settle/dillon$meanDepth/365)/((dillon$settle/dillon$meanDepth/365)+dillon$runoff/365/dillon$meanDepth)
plot(dillon$ret~dillon$modeledRet)
abline(0,1)

# see notebook for derivation 
all$emergent_d_all=(all$emergent_d_epi*all$Vepi+all$emergent_d_hypo*all$Vhypo)/(all$Vepi+all$Vhypo)
fracRet=all$emergent_d_all/(all$emergent_d_all+(all$LakeE/all$Area/all$percentEvap-all$LakeE/all$Area)/(all$Vol/all$Area))

fracRet=mean(all$emergent_d_all)/(mean(all$emergent_d_all)+(all$LakeE/all$Area/all$percentEvap-all$LakeE/all$Area)/(all$Vol/all$Area))

plot(fracRet~all$FracRet)
abline(0,1,lwd=2,col='red')
abline(lm(fracRet~all$FracRet))
summary(lm(fracRet~all$FracRet))

# derivation including shape and rate parameters from Vachon et al. 2016. v is about 0.12, which makes alpha about 13.33 
v=.62
alpha=108
decay=(0.85*v)/(all$meanDepth/(all$LakeE/all$Area/all$percentEvap-all$LakeE/all$Area))*
  log(1+(all$meanDepth/(all$LakeE/all$Area/all$percentEvap-all$LakeE/all$Area))/alpha)

fracRet=decay/(decay+(all$LakeE/all$Area/all$percentEvap-all$LakeE/all$Area)/(all$Vol/all$Area))

fracRet=0.85*v*log(1+(all$meanDepth/(all$LakeE/all$Area/all$percentEvap-all$LakeE/all$Area))/alpha)/
  (0.85*v*log(1+(all$meanDepth/(all$LakeE/all$Area/all$percentEvap-all$LakeE/all$Area))/alpha)+1)

plot(fracRet~all$FracRet,ylim=c(0,1),xlim=c(0,1))
abline(0,1)

library(scales)
gamma_est<-function(p,x){
  -sum(dgamma(x,shape=exp(p[1]),scale=exp(p[2]),log=TRUE))
}

guess=log(c(1.052,100))
xout=seq(0.00001,.1,length.out = 100)

initialReactivity=c(rep(0.0015,9),rep(0.08,1))
par=exp(optim(guess,gamma_est,x=initialReactivity)$par)
out=dgamma(xout,shape=par[1],scale=par[2])
out=out/max(out)

windows()
plot(xout,out,type='l')

v=par[2] # shape 
alpha=par[1] # rate 

alpha=12.5 # average lifetime of most reactivit compounds in days 
v=.9 # relative abundance of the most recalcitrant compounds 

fracRet=0.85*v*log(1+(all$meanDepth/(all$LakeE/all$Area/all$percentEvap-all$LakeE/all$Area))/alpha)/
  (0.85*v*log(1+(all$meanDepth/(all$LakeE/all$Area/all$percentEvap-all$LakeE/all$Area))/alpha)+1)

fracRet=0.85*v*log(1+(all$meanDepth/(all$LakeE/all$Area/all$percentEvap-all$LakeE/all$Area))/alpha)/
  (0.85*v*log(1+(all$meanDepth/(all$LakeE/all$Area/all$percentEvap-all$LakeE/all$Area))/alpha)+1)

plot(fracRet~all$FracRet,ylim=c(0,1),xlim=c(0,1))
abline(0,1)


windows()
png('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Figures/Fig3_fracRet_RespRate.png',width = 7,
    height=7,units = 'in',res = 300)
par(mar=c(6,6,5,2))
cex.axis=2.5
cex.lab=2.5
cex=2
col=colorpalette(c('grey70','black'))
heatscatter(all$FracRet,log10(all$DOC_Respired/all$Area*12*365),pch = 19,colpal = col,cex=cex,main='',cex.axis=cex.axis,cex.lab=cex.lab,xlim=c(0,1),
            xlab='Fraction Retained',ylab='DOC Respired (g C m-2 yr-1)',yaxt='n')
axis(2,at = c(log10(0.10),log10(1),log10(10),log10(100),log10(1000)),labels = c(0.1,1,10,100,1000),cex.axis=cex.axis)
dev.off()

png('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Figures/Fig3_fracRet_fracEvap.png',width = 7,
    height=7,units = 'in',res = 300)
par(mar=c(6,6,5,2))
cex.axis=2.5
cex.lab=2.5
cex=2
col=colorpalette(c('grey70','black'))
heatscatter(all$percentEvap,all$FracRet,pch = 19,colpal = col,cex=cex,main='',cex.axis=cex.axis,cex.lab=cex.lab,xlim=c(0,1),
            xlab='Fraction Export as Evap',ylab='Fraction Retained')
abline(0,1)
abline(lm(all$FracRet~all$percentEvap))
all$meanDepth=all$Vol/all$Area
summary(lm(all$FracRet~all$percentEvap))
summary(lm(all$FracRet~all$percentEvap+all$meanDepth))
summary(lm(all$FracRet~all$percentEvap+all$emergent_d_all))


summary(lm(all$FracRet~all$percentEvap*all$meanDepth))
summary(lm(all$FracRet~(all$percentEvap/all$meanDepth)))
res=resid(lm(all$FracRet~all$percentEvap))
plot(res~all$meanDepth)
abline(lm(res~all$meanDepth))
summary(lm(res~all$meanDepth))
dev.off()

fit=nls(FracRet~a*exp(percentEvap)+b,data = all,algorithm = 'port',
        start = c(b=.03,a=-1),lower=c(b=-Inf,a=-Inf),upper=c(b=Inf,a=Inf))
summary(fit)
all$predFracRet=predict(fit)
plot(all$FracRet~all$percentEvap,pch=16,cex=2)
temp=all[sort.list(all$predFracRet),]
lines(temp$predFracRet~temp$percentEvap,col='red')
abline(0,1,col='red',lty=2,lwd=2)

summary(lm(all$FracRet~all$percentEvap))




png('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Figures/Fig3_emit_fracEvap.png',width = 7,
    height=7,units = 'in',res = 300)
par(mar=c(6,6,5,2))
cex.axis=2.5
cex.lab=2.5
cex=2
col=colorpalette(c('grey70','black'))
heatscatter(all$percentEvap,log10(all$Emit/all$Area*12*365),pch = 19,colpal = col,cex=cex,main='',cex.axis=cex.axis,cex.lab=cex.lab,xlim=c(0,1),
            xlab='Fraction Export as Evap',ylab='Emissions (g C m-2 yr-1)',yaxt='n',ylim=c(log10(3),log10(1000)))
axis(2,at = c(log10(1),log10(10),log10(100),log10(1000)),labels = c(1,10,100,1000),cex.axis=cex.axis)
dev.off()

png('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Figures/Fig3_emit_NEP.png',width = 7,
    height=7,units = 'in',res = 300)
par(mar=c(6,6,5,2))
cex.axis=2.5
cex.lab=2.5
cex=2
col=colorpalette(c('grey70','black'))
heatscatter(sum$NEP/sum$Area*12*1000,log10(all$Emit/all$Area*12*365),pch = 19,colpal = col,cex=cex,main='',cex.axis=cex.axis,cex.lab=cex.lab,
            xlab='Fraction Export as Evap',ylab='Emissions (g C m-2 yr-1)',yaxt='n',ylim=c(log10(3),log10(1000)))
axis(2,at = c(log10(1),log10(10),log10(100),log10(1000)),labels = c(1,10,100,1000),cex.axis=cex.axis)
dev.off()

png('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Figures/Fig3_burial_fracEvap.png',width = 7,
    height=7,units = 'in',res = 300)
par(mar=c(6,6,5,2))
cex.axis=2.5
cex.lab=2.5
cex=2
col=colorpalette(c('grey70','black'))
heatscatter(all$percentEvap,log10((all$Burial_tPOC+all$Burial_phyto)/all$Area*12*365),pch = 19,colpal = col,cex=cex,main='',cex.axis=cex.axis,cex.lab=cex.lab,xlim=c(0,1),
            xlab='Fraction Export as Evap',ylab='Burial (g C m-2 yr-1)',yaxt='n',ylim=c(log10(.5),log10(100)))
axis(2,at = c(log10(1),log10(10),log10(100),log10(1000)),labels = c(1,10,100,1000),cex.axis=cex.axis)
dev.off()

summary(lm(log10((all$Burial_tPOC+all$Burial_phyto)/all$Area*12*365)~all$percentEvap+all$meanDepth))

png('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Figures/Fig3_dicLoad_RespRate.png',width = 7,
    height=7,units = 'in',res = 300)
par(mar=c(6,6,5,2))
cex.axis=2.5
cex.lab=2.5
cex=2
col=colorpalette(c('grey70','black'))
heatscatter(log10(all$DIC_load/all$Area*12*365),log10(all$DOC_Respired/all$Area*12*365),pch = 19,colpal = col,cex=cex,main='',cex.axis=cex.axis,cex.lab=cex.lab,
            xlab='DIC Load (g C m-2 yr-1)',ylab='DOC Respired (g C m-2 yr-1)',yaxt='n',xaxt='n',ylim=c(0,log10(1000)),xlim=c(0,log10(50000)))
axis(2,at = c(log10(0.10),log10(1),log10(10),log10(100),log10(1000)),labels = c(0.1,1,10,100,1000),cex.axis=cex.axis)
axis(1,at = c(log10(0.10),log10(1),log10(10),log10(100),log10(1000),log10(10000)),labels = c(0.1,1,10,100,1000,10000),cex.axis=cex.axis)
abline(0,1,lty=2,lwd=2)
dev.off()

png('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Figures/Fig3_dictoResp_fracEvap.png',width = 7,
    height=7,units = 'in',res = 300)
par(mar=c(6,6,5,2))
cex.axis=2.5
cex.lab=2.5
cex=2
col=colorpalette(c('grey70','black'))
heatscatter(all$percentEvap,log10(all$dicLoadvResp),pch = 19,colpal = col,cex=cex,main='',cex.axis=cex.axis,cex.lab=cex.lab,
            xlab='Fraction Export as Evap',ylab='DIC Load : DIC Produced',yaxt='n')
axis(2,at = c(log10(0.10),log10(1),log10(10),log10(100),log10(1000),log10(10000)),labels = c(0.1,1,10,100,1000,10000),cex.axis=cex.axis)
abline(h=log10(1),lty=2,lwd=2)
dev.off()

png('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Figures/Fig3_dictoResp_emit.png',width = 7,
    height=7,units = 'in',res = 300)
par(mar=c(6,6,5,2))
cex.axis=2.5
cex.lab=2.5
cex=2
col=colorpalette(c('grey70','black'))
heatscatter(log10(all$Emit/all$Area*12*365),log10(all$dicLoadvResp),pch = 19,colpal = col,cex=cex,main='',cex.axis=cex.axis,cex.lab=cex.lab,
            xlab='Emissions (g C m-2 yr-1)',ylab='DIC Load : DIC Produced',yaxt='n',xaxt='n',xlim=c(log10(3),log10(1000)))
axis(2,at = c(log10(0.10),log10(1),log10(10),log10(100),log10(1000),log10(10000)),labels = c(0.1,1,10,100,1000,10000),cex.axis=cex.axis)
axis(1,at = c(log10(0.10),log10(1),log10(10),log10(100),log10(1000),log10(10000)),labels = c(0.1,1,10,100,1000,10000),cex.axis=cex.axis)
abline(h=log10(1),lty=2,lwd=2)
points(log10(all$Emit[sum$NEP>0]/all$Area[sum$NEP>0]*12*365),log10(all$dicLoadvResp[sum$NEP>0]),col='green',cex=cex,
       pch=16)
legend('topleft',legend = c('Heterotrophic','Autotrophic'),
       col=c('grey 40','green'),bty = 'n',pch = 19,cex = 1.5,pt.cex = 2)
dev.off()

png('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Figures/Fig3_dictoResp_emitTotal.png',width = 7,
    height=7,units = 'in',res = 300)
par(mar=c(6,6,5,2))
cex.axis=2.5
cex.lab=2.5
cex=2
col=colorpalette(c('grey70','black'))
heatscatter(log10(all$Emit*12*365),log10(all$dicLoadvResp),pch = 19,colpal = col,cex=cex,main='',cex.axis=cex.axis,cex.lab=cex.lab,
            xlab='Emissions (g C yr-1)',ylab='DIC Load : DOC Respired',yaxt='n',xaxt='n')
axis(2,at = c(log10(0.10),log10(1),log10(10),log10(100),log10(1000),log10(10000)),labels = c(0.1,1,10,100,1000,10000),cex.axis=cex.axis)
axis(1,at = c(log10(1e4),log10(1e6),log10(1e8)),labels = c(1e4,1e6,1e8),cex.axis=cex.axis)
abline(h=log10(1),lty=2,lwd=2)
points(log10(all$Emit[sum$NEP>0]*12*365),log10(all$dicLoadvResp[sum$NEP>0]),col='green',cex=cex,
       pch=16)
dev.off()

png('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Figures/Fig3_dictoResp_Area.png',width = 7,
    height=7,units = 'in',res = 300)
par(mar=c(6,6,5,2))
cex.axis=2.5
cex.lab=2.5
cex=2
col=colorpalette(c('grey70','black'))
heatscatter(log10(all$Area),log10(all$dicLoadvResp),pch = 19,colpal = col,cex=cex,main='',cex.axis=cex.axis,cex.lab=cex.lab,
            xlab='Area (m-2)',ylab='DIC Load : DOC Respired',yaxt='n',xaxt='n')
axis(2,at = c(log10(0.10),log10(1),log10(10),log10(100),log10(1000),log10(10000)),labels = c(0.1,1,10,100,1000,10000),cex.axis=cex.axis)
axis(1,at = c(log10(1e4),log10(1e6),log10(1e8)),labels = c(1e4,1e6,1e8),cex.axis=cex.axis)
abline(h=log10(1),lty=2,lwd=2)
points(log10(all$Area[sum$NEP>0]),log10(all$dicLoadvResp[sum$NEP>0]),col='green',cex=cex,
       pch=16)
dev.off()

windows()
png('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Figures/Fig5_Burial_HRT.png',width = 7,
    height=7,units = 'in',res = 300)
par(mar=c(6,6,5,2))
cex.axis=2.5
cex.lab=2.5
cex=2
col=colorpalette(c('grey70','black'))
heatscatter(log10(all$HRT),log10((all$Burial_tPOC+all$Burial_phyto)/all$Area*12*365),pch = 19,colpal = col,cex=cex,main='',cex.axis=cex.axis,cex.lab=cex.lab,
            xlab='HRT (days)',ylab='Sed Burial (g C m-2 yr-1)',yaxt='n',xaxt='n',ylim=c(log10(0.5),log10(100)))
axis(2,at = c(log10(0.10),log10(1),log10(10),log10(100),log10(1000)),labels = c(0.1,1,10,100,1000),cex.axis=cex.axis)
axis(1,at = c(log10(0.10),log10(1),log10(10),log10(100),log10(1000)),labels = c(0.1,1,10,100,1000),cex.axis=cex.axis)
dev.off()

png('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Figures/Fig5_Burial_TP.png',width = 7,
    height=7,units = 'in',res = 300)
par(mar=c(6,6,5,2))
cex.axis=2.5
cex.lab=2.5
cex=2
col=colorpalette(c('grey70','black'))
heatscatter(log10((all$P_epi+all$phyto/95)/all$Vepi*1000*31),log10((all$Burial_tPOC+all$Burial_phyto)/all$Area*12*365),pch = 19,colpal = col,cex=cex,main='',cex.axis=cex.axis,cex.lab=cex.lab,
            xlab='TP ug L-1',ylab='Sed Burial (g C m-2 yr-1)',yaxt='n',xaxt='n',xlim=c(log10(0.5),log10(100)),ylim=c(log10(0.5),log10(100)))
axis(2,at = c(log10(0.10),log10(1),log10(10),log10(100),log10(1000)),labels = c(0.1,1,10,100,1000),cex.axis=cex.axis)
axis(1,at = c(log10(0.10),log10(1),log10(10),log10(100),log10(1000)),labels = c(0.1,1,10,100,1000),cex.axis=cex.axis)
dev.off()


windows()
png('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Figures/Fig6_Emit_Area.png',width = 7,
    height=7,units = 'in',res = 300)
par(mar=c(6,6,5,2))
cex.axis=2.5
cex.lab=2.5
cex=2
col=colorpalette(c('grey70','black'))
heatscatter(log10(all$Area),log10(all$Emit/all$Area*12*365),pch = 19,colpal = col,cex=cex,main='',cex.axis=cex.axis,cex.lab=cex.lab,
            xlab='Lake Area (m2)',ylab='Emissions (g C m-2 yr-1)',yaxt='n',xaxt='n',ylim=c(log10(3),log10(1000)))
axis(2,at = c(log10(1),log10(10),log10(100),log10(1000)),labels = c(1,10,100,1000),cex.axis=cex.axis)
axis(1,at = c(log10(100),log10(10000),log10(1000000)),labels = c(100,10000,1000000),cex.axis=cex.axis)
dev.off()

png('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Figures/Fig6_Burial_Area.png',width = 7,
    height=7,units = 'in',res = 300)
par(mar=c(6,6,5,2))
cex.axis=2.5
cex.lab=2.5
cex=2
col=colorpalette(c('grey70','black'))
heatscatter(log10(all$Area),log10((all$Burial_tPOC+all$Burial_phyto)/all$Area*12*365),pch = 19,colpal = col,cex=cex,main='',cex.axis=cex.axis,cex.lab=cex.lab,
            xlab='Lake Area (m2)',ylab='Burial (g C m-2 yr-1)',yaxt='n',xaxt='n',ylim=c(log10(.5),log10(100)))
axis(2,at = c(log10(1),log10(10),log10(100),log10(1000)),labels = c(1,10,100,1000),cex.axis=cex.axis)
axis(1,at = c(log10(100),log10(10000),log10(1000000)),labels = c(100,10000,1000000),cex.axis=cex.axis)
dev.off()

png('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Figures/Fig6_WALA_Area.png',width = 7,
    height=7,units = 'in',res = 300)
par(mar=c(6,6,5,2))
cex.axis=2.5
cex.lab=2.5
cex=2
col=colorpalette(c('grey70','black'))
heatscatter(log10(all$Area),log10(all$WALA),pch = 19,colpal = col,cex=cex,main='',cex.axis=cex.axis,cex.lab=cex.lab,
            xlab='Lake Area (m2)',ylab='WA:LA',yaxt='n',xaxt='n')
axis(2,at = c(log10(1),log10(10),log10(100),log10(1000),log10(10000)),labels = c(1,10,100,1000,10000),cex.axis=cex.axis)
axis(1,at = c(log10(100),log10(10000),log10(1000000)),labels = c(100,10000,1000000),cex.axis=cex.axis)
dev.off()

png('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Figures/Fig6_HRT_Area.png',width = 7,
    height=7,units = 'in',res = 300)
par(mar=c(6,6,5,2))
cex.axis=2.5
cex.lab=2.5
cex=2
col=colorpalette(c('grey70','black'))
heatscatter(log10(all$Area),all$HRT/365,pch = 19,colpal = col,cex=cex,main='',cex.axis=cex.axis,cex.lab=cex.lab,
            xlab='Lake Area (m2)',ylab='HRT (years)',xaxt='n')
# axis(2,at = c(log10(1),log10(10),log10(100),log10(1000),log10(10000)),labels = c(1,10,100,1000,10000),cex.axis=cex.axis)
axis(1,at = c(log10(100),log10(10000),log10(1000000)),labels = c(100,10000,1000000),cex.axis=cex.axis)
dev.off()

png('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Figures/Fig6_FracEvap_Area.png',width = 7,
    height=7,units = 'in',res = 300)
par(mar=c(6,6,5,2))
cex.axis=2.5
cex.lab=2.5
cex=2
col=colorpalette(c('grey70','black'))
heatscatter(log10(all$Area),all$percentEvap,pch = 19,colpal = col,cex=cex,main='',cex.axis=cex.axis,cex.lab=cex.lab,
            xlab='Lake Area (m2)',ylab='Fraction Export as Evap',xaxt='n')
# axis(2,at = c(log10(1),log10(10),log10(100),log10(1000),log10(10000)),labels = c(1,10,100,1000,10000),cex.axis=cex.axis)
axis(1,at = c(log10(100),log10(10000),log10(1000000)),labels = c(100,10000,1000000),cex.axis=cex.axis)
dev.off()

png('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Figures/Fig6_d_Area.png',width = 7,
    height=7,units = 'in',res = 300)
par(mar=c(6,6,5,2))
cex.axis=2.5
cex.lab=2.5
cex=2
col=colorpalette(c('grey70','black'))
heatscatter(log10(all$Area),all$emergent_d_epi,pch = 19,colpal = col,cex=cex,main='',cex.axis=cex.axis,cex.lab=cex.lab,
            xlab='Lake Area (m2)',ylab='Turnover Rate DOC (day-1)',xaxt='n')
# axis(2,at = c(log10(1),log10(10),log10(100),log10(1000),log10(10000)),labels = c(1,10,100,1000,10000),cex.axis=cex.axis)
axis(1,at = c(log10(100),log10(10000),log10(1000000)),labels = c(100,10000,1000000),cex.axis=cex.axis)
dev.off()


plot(all$percentEvap~all$HRT)

plot(all$percentEvap,all$HRT*all$LakeE/all$Area/all$meanDepth)



windows()
boxplot(sum$Emit/sum$Area*12*365~sum$lakeSizeBins,outline=T)
boxplot((sum$Burial_tPOC+sum$Burial_phyto)/sum$Area*12*365~sum$lakeSizeBins,outline=T)
boxplot(sum$FracRet[sum$FracRet>0]~sum$lakeSizeBins[sum$FracRet>0],outline=F)

boxplot(all$Emit/all$Area*12*365~all$lakeSizeBins,outline=T)
boxplot((all$Burial_tPOC+all$Burial_phyto)/all$Area*12*365~all$lakeSizeBins,outline=T)
boxplot(all$FracRet[all$FracRet>0]~all$lakeSizeBins[all$FracRet>0],outline=F)
boxplot(all$percentEvap[all$FracRet>0]~all$lakeSizeBins[all$FracRet>0],outline=F)

summary(sum$Emit/sum$Area*12*365~sum$lakeSizeBins) # g C m-2 yr-1 
summary(all$Emit/all$Area*12*365~all$lakeSizeBins)

all$allBurial=all$Burial_tPOC+all$Burial_phyto
sum$allBurial=sum$Burial_tPOC+sum$Burial_phyto
summary(all$allBurial/all$Area*12*365~all$lakeSizeBins)
summary(sum$allBurial/sum$Area*12*365~sum$lakeSizeBins)

summary(all$FracRet~all$lakeSizeBins)

summary(all$HRT~all$lakeSizeBins)

summary(sum$emergent_d_epi~sum$lakeSizeBins)

sum(sum$k*sum$Area)/sum(sum$Area)

sum(sum$pH*sum$Area)/sum(sum$Area)

sum(.952*sum$fracCO2*sum$DIC_epi/sum$Vepi*1000*29.41*sum$Area)/sum(sum$Area)

sum(all$allBurial*12*365)/sum(all$Area)
sum(sum$allBurial*12*365)/sum(sum$Area)

sum(all$Emit*12*365)/sum(all$Area)

sum(sum$Emit*12*365)/sum(sum$Area)

sum(all$FracRet*all$Area)/sum(all$Area)

sum(all$HRT*all$Area)/sum(all$Area)

sum(sum$emergent_d_epi*sum$Area)/sum(sum$Area)

sum(all$emergent_d_epi*all$Area)/sum(all$Area)

all$allDOC<-all$DOCr_epi+all$DOCl_epi
summary(all$allDOC/all$Vepi*12~all$lakeSizeBins)

boxplot(sum$DOCr_epi/sum$Vepi*12~sum$lakeSizeBins,outline=T)

boxplot(sum$DIC_epi/sum$Vepi*12~sum$lakeSizeBins,outline=T)

summary(all$dicLoadvResp)

summary(all$dicLoadvResp~all$lakeSizeBins)

summary(all$percentEvap~all$lakeSizeBins)

sum(all$percentEvap*all$Area)/sum(all$Area)

anova(lm(all$dicLoadvResp~all$lakeSizeBins))

anova(lm(all$HRT~all$lakeSizeBins))
summary(lm(all$HRT~log10(all$Area)))

summary(lm(all$waterIn/all$Area~log10(all$Area)))
summary(lm(log10(all$Vol)~log10(all$Area)))
summary(lm(all$percentEvap~log10(all$Area)))

anova(lm(all$FracRet~all$lakeSizeBins))

anova(lm(all$WALA~all$lakeSizeBins))
boxplot(all$WALA~all$lakeSizeBins,outline=F)
plot(log10(all$WALA)~log10(all$Area))
summary(lm(log10(all$WALA)~log10(all$Area)))

summary(lm(all$emergent_d_epi~log10(all$Area)))

anova(lm(all$percentEvap~all$lakeSizeBins))

anova(lm(all$Emit/all$Area~all$lakeSizeBins))

boxplot(all$Emit/all$Area~all$lakeSizeBins)

anova(lm(all$k~all$lakeSizeBins))


png('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Figures/fig8_fracEmitBury_fracLakes.png',width = 14,height = 7,res=300,units='in')
par(mar=c(6,6,5,2),mfrow=c(1,2))
cex.axis=2.5
cex.lab=2.5
lwd=8
pch=16
cex=2
temp=all[sort.list(all$Emit),]
temp$i=seq(1,nrow(temp),1)/nrow(temp)
temp$EmitFrac=temp$Emit/sum(temp$Emit)
temp$EmitFrac=cumsum(temp$EmitFrac)
plot(temp$i~temp$EmitFrac,cex=cex,pch=pch,lwd=lwd,xlab='Fraction of Total Emissions',ylab='Cumulative Fraction',cex.lab=cex.lab,cex.axis=cex.axis)
v=0.1
h=temp$i[min(which(temp$EmitFrac>v))] # fraction of lakes for which 90% of emissions are accounted for  
abline(h=h,lty=2,lwd=3)
abline(v=v,lty=2,lwd=3)

temp=all[sort.list(all$Burial_phyto+all$Burial_tPOC),]
temp$i=seq(1,nrow(temp),1)/nrow(temp)
temp$BuryFrac=(temp$Burial_tPOC+temp$Burial_phyto)/sum(temp$Burial_tPOC+temp$Burial_phyto)
temp$BuryFrac=cumsum(temp$BuryFrac)
plot(temp$i~temp$BuryFrac,cex=cex,pch=pch,lwd=lwd,xlab='Fraction of Total Burial',ylab='',cex.lab=cex.lab,cex.axis=cex.axis)
v=0.1
h=temp$i[min(which(temp$BuryFrac>v))] # fraction of lakes for which 90% of emissions are accounted for  
abline(h=h,lty=2,lwd=3)
abline(v=v,lty=2,lwd=3)
dev.off() 

temp=all[sort.list(all$Burial_phyto+all$Burial_tPOC),]
temp$i=seq(1,nrow(temp),1)/nrow(temp)
temp$BuryFrac=(temp$Burial_tPOC+temp$Burial_phyto)/sum(temp$Burial_tPOC+temp$Burial_phyto)
temp$BuryFrac[temp$Permanent_==ntlLookUp$Permanent_[ntlLookUp$lakeID=='TR']]

temp=all[sort.list(all$Emit),]
temp$i=seq(1,nrow(temp),1)/nrow(temp)
temp$EmitFrac=temp$Emit/sum(temp$Emit)
temp$EmitFrac[temp$Permanent_==ntlLookUp$Permanent_[ntlLookUp$lakeID=='TR']]

temp$Emit[temp$Permanent_==ntlLookUp$Permanent_[ntlLookUp$lakeID=='TR']]/temp$Area[temp$Permanent_==ntlLookUp$Permanent_[ntlLookUp$lakeID=='TR']]*12*365

(temp$Burial_tPOC[temp$Permanent_==ntlLookUp$Permanent_[ntlLookUp$lakeID=='TR']]+temp$Burial_tPOC[temp$Permanent_==ntlLookUp$Permanent_[ntlLookUp$lakeID=='TR']])/
  temp$Area[temp$Permanent_==ntlLookUp$Permanent_[ntlLookUp$lakeID=='TR']]*12*365



sum(temp$Emit[1:37])/sum(temp$Emit)


boxplot(all$percentEvap~all$lakeSizeBins)

boxplot(all$dicLoadvResp~all$lakeSizeBins,outline=F)

plot(log10(all$dicLoadvResp)~log10(all$Area))

plot(log10(all$dicLoadvResp)~log10(all$Emit))

boxplot(sum$fracCO2~sum$lakeSizeBins,outline=T)

all$totalDICcontributed=all$DOC_Respired+all$sed_resp+all$DIC_load
sum(all$DOC_Respired)/sum(all$totalDICcontributed)
sum(all$DIC_load)/sum(all$totalDICcontributed)
sum(all$sed_resp)/sum(all$totalDICcontributed)

sum(all$Emit[all$dicLoadvResp>1])/sum(all$Emit)
sum(all$Emit[all$dicLoadvResp<1])/sum(all$Emit)

# open ice GPP 
sum$NEP=sum$GPP*sum$Vepi*.15-sum$DOC_Respired # mol C day-1 in NPP
summary(sum$NEP)
boxplot(sum$NEP/sum$Area*12*1000,outline=F) # NEP in mg C m-2 day-1
abline(0,0)

sum$DOC_Respired_Epi=sum$emergent_d_epi*(sum$DOCr_epi+sum$DOCl_epi)
sum$NEP=sum$GPP*sum$Vepi*.15-sum$DOC_Respired_Epi

sum[sum$NEP>0,]

summary(all$Emit[sum$NEP>0]/all$Area[sum$NEP>0]*12*365) # g C m-2 yr-1 (summertime rates) 
summary(all$Emit[sum$NEP<0]/all$Area[sum$NEP<0]*12*365) # g C m-2 yr-1 (summertime rates) 
summary(all$Emit/all$Area*12*365)
sum=sum[sort.list(sum$Permanent_),]
all=all[sort.list(all$Permanent_),]

hist(all$HRT)
hist(all$HRT[sum$NEP>0],col='red',add=T)

netAutoEmit=all$Emit[sum$NEP>0]/all$Area[sum$NEP>0]*12*365 # g C m-2 yr-1 (summertime rates) 

hist(log10(all$Emit/all$Area*12*365))
hist(log10(netAutoEmit),add=T,col='red')

gpp=sum$GPP*sum$Vepi/sum$Area*12*1000 

boxplot(sum$fracCO2*sum$DIC_epi/sum$Vepi*12~sum$lakeSizeBins,outline=T)
summary(sum$fracCO2*sum$DIC_epi/sum$Vepi*12~sum$lakeSizeBins)

summary(.952*sum$fracCO2*sum$DIC_epi/sum$Vepi*1000*29.41~sum$lakeSizeBins) # pCO2 
summary(.952*all$fracCO2*all$DIC_epi/all$Vepi*1000*29.41~all$lakeSizeBins) # pCO2 
summary(sum$k~sum$lakeSizeBins)
summary(all$k~all$lakeSizeBins)
summary(sum$pH~sum$lakeSizeBins)
summary(.952*sum$fracCO2*sum$DIC_epi/sum$Vepi*1000*29.41) # pCO2 

hist(.952*sum$fracCO2*sum$DIC_epi/sum$Vepi*1000*29.41)
hist(.952*sum$fracCO2*sum$DIC_epi/sum$Vepi*1000*29.41/390,breaks = c(-1.2,1.2,2,4,8,16,35))



plot(sum$fracCO2*sum$DIC_epi/sum$Vepi*12~sum$DOC_conc)
abline(lm(sum$fracCO2*sum$DIC_epi/sum$Vepi*12~sum$DOC_conc))
summary(lm(sum$fracCO2*sum$DIC_epi/sum$Vepi*12~sum$DOC_conc))

plot(.952*sum$fracCO2*sum$DIC_epi/sum$Vepi*1000*29.41~sum$percentEvap)
plot(sum$Emit/sum$Area*1000*12~sum$percentEvap)
plot(sum$Emit/sum$Area*1000*12~sum$HRT)
plot(sum$Emit/sum$Area*1000*12~log10(sum$Area))
plot(sum$Emit/sum$Area*1000*12~log10(sum$WALA))

summary(sum$Emit[sum$Area<10000]/sum$Area[sum$Area<10000]*12*1000)

windows()
plot(log10(sum$WALA)~log10(sum$Area))

hanson<-read.csv('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Validation Data/NHLD Random Sampling/random_lake_survey_measured_parameters.csv',
                 stringsAsFactor=F)
hansonLoc<-read.csv('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Validation Data/NHLD Random Sampling/random_lake_survey_lakes.csv',
                    stringsAsFactor=F)
hanson<-merge(hanson,hansonLoc,by='wbic',all.x=T)
# hanson<-hanson[hanson$doc>1&hanson$area>10000,] # keeping only greater than 1 ha for now
wbic_perm=readShapeSpatial('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/GIS/NHLD_Lakes_NHD_WBIC_Join.shp')
wbic_perm=data.frame(wbic_perm)
wbic_perm=wbic_perm[,c('Permanent_','WATERBODY_')]
colnames(wbic_perm)[2]='wbic'
hanson=merge(hanson,wbic_perm,by='wbic',all.x=T)
hanson$Permanent_<-as.character(hanson$Permanent_)

windows()
hist(sum$pH)
hist(hanson$pH)
hist(sum$fracCO2)
plot(sum$pH~log10(sum$Area),ylim=range(c(hanson$pH,sum$pH),na.rm = T))
abline(lm(sum$pH~log10(sum$Area)))
points(hanson$pH~log10(hanson$area),col='red',pch=16)
abline(lm(hanson$pH~log10(hanson$area)),col='red')

summary(hanson$pH)
summary(sum$pH)
boxplot(sum$pH,at = c(1),col='red',ylim=range(hanson$pH,sum$pH))
boxplot(hanson$pH,add = T)

hist(log10(0.952*sum$fracCO2*sum$DIC_epi/sum$Vepi/1000*1e6*29.41),xaxt='n') # ppm 
axis(1,c(log10(500),log10(1000),log10(5000)),c(500,1000,5000))

hist(-0.5*(400-sum$fracCO2*sum$DIC_epi/sum$Vepi/1000*1e6*29.41)/sum$zmix)

windows()
hist((sum$DOCr_epi+sum$DOCl_epi)/sum$Vepi*12)
summary((sum$DOCr_epi+sum$DOCl_epi)/sum$Vepi*12)
summary((sum$DOCr_epi+sum$DOCl_epi+sum$DOCl_hypo+sum$DOCr_hypo)/(sum$Vepi+sum$Vhypo)*12)

# cumulative fraction 
png('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Figures/hanson_doc_cumFreq.png',
    width = 7,height = 7,res = 300,units = 'in')
par(mar=c(6,6,5,2))
cex.axis=2.5
cex.lab=2.5
cex=2
lwd=8
sum$doc_conc<-(sum$DOCr_epi+sum$DOCl_epi+sum$DOCr_hypo+sum$DOCl_hypo)/(sum$Vol)*12
sum<-sum[sort.list(sum$doc_conc),]
sum$i<-seq(1,length(sum$SWin),1)
sum$i<-sum$i/length(sum$SWin)
hanson<-hanson[sort.list(hanson$doc),]
hanson$i<-seq(1,length(hanson$wbic),1)/length(hanson$wbic)
temp=hanson[hanson$doc>0,]
temp=temp[sort.list(temp$doc),]
temp$area<-temp$area/sum(temp$area)
temp$area<-cumsum(temp$area)
tempSum=sum[sort.list(sum$doc_conc),]
tempSum$Area=tempSum$Area/sum(tempSum$Area)
tempSum$Area=cumsum(tempSum$Area)
plot(hanson$i[hanson$doc>0]~hanson$doc[hanson$doc>0],type='l',lwd=lwd,xlab='DOC mg L-1',cex.lab=cex.lab,
     cex.axis=cex.axis,ylab='Cumulaive Fraction',col='grey60')
lines(sum$i~sum$doc_conc,lwd=8)
lines(temp$area~temp$doc,lwd=8,lty=3,col='grey60')
lines(tempSum$Area~tempSum$doc_conc,lwd=8,lty=3)
# axis(1,at=c(log10(5),log10(10),log10(30),log10(100)),labels=c(5,10,30,100),cex.axis=2)
legend('topleft',legend = c('Modeled','Random Sampling'),fill = c('black','grey60'),bty='n',cex = 1.5)
dev.off()

png('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Figures/hanson_doc_1to1.png',
    width = 7,height = 7,res = 300,units = 'in')
temp=merge(hanson,sum,by='Permanent_',all.x=T)
par(bty='l',mar=c(5,6,5,5))
plot(temp$doc[temp$doc>0]~temp$doc_conc[temp$doc>0],pch=16,cex=2,xlim=c(0,35),ylim=c(0,35),cex.axis=2,
     cex.lab=2,ylab='Observed (mg C L-1)',xlab='Modeled (mg C L-1)') # see Canham et al 2004 Fig. 2 
abline(0,1,lwd=3,lty=2)
text(x = 30,y = 25,labels = '1:1',cex = 2)
dev.off()
# abline(lm(temp$doc[temp$doc>0]~temp$doc_conc[temp$doc>0]),lwd=2)
summary(lm(temp$doc[temp$doc>0]~temp$doc_conc[temp$doc>0]))
summary(hanson$doc[hanson$doc>0])
summary(sum$doc_conc)

all$doc_conc<-(all$DOCr_epi+all$DOCl_epi+all$DOCr_hypo+all$DOCl_hypo)/(all$Vol)*12
all<-all[sort.list(all$doc_conc),]
all$i<-seq(1,length(all$SWin),1)
all$i<-all$i/length(all$SWin)
hanson<-hanson[sort.list(hanson$doc),]
hanson$i<-seq(1,length(hanson$wbic),1)/length(hanson$wbic)
plot(hanson$i~log10(hanson$doc),type='l',lwd=8,xlim=c(log10(2.5),log10(60)),xlab='DOC mg L-1',xaxt='n',cex.lab=1.5,
     cex.axis=2,ylab='Cumulaive Fraction',col='grey60')
lines(all$i~log10(all$doc_conc),lwd=8)
axis(1,at=c(log10(5),log10(10),log10(30),log10(100)),labels=c(5,10,30,100),cex.axis=2)
legend('topleft',legend = c('Modeled','Random Sampling'),fill = c('black','grey60'),bty='n' )
temp=merge(hanson,all,by='Permanent_',all.x=T)
plot(temp$doc[temp$doc>0]~temp$doc_conc[temp$doc>0],pch=16,cex=2,xlim=c(0,35),ylim=c(0,35),cex.axis=2,
     cex.lab=2,ylab='Observed',xlab='Modeled') # see Canham et al 2004 Fig. 2 
abline(0,1,lwd=2,lty=2)
abline(lm(temp$doc[temp$doc>0]~temp$doc_conc[temp$doc>0]),lwd=2)
summary(lm(temp$doc[temp$doc>0]~temp$doc_conc[temp$doc>0]))
summary(hanson$doc[hanson$doc>0])
summary(sum$doc_conc)

windows()
hist((sum$DIC_epi)/sum$Vepi*12)
png('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Figures/hanson_dic_cumFreq.png',
    width = 7,height = 7,res = 300,units = 'in')
par(mar=c(6,6,5,2))
cex.axis=2.5
cex.lab=2.5
cex=2
lwd=8
hanson<-hanson[sort.list(hanson$dic),]
hanson$i<-seq(1,length(hanson$wbic),1)/length(hanson$wbic)
plot(hanson$i[hanson$dic>0]~hanson$dic[hanson$dic>0],type='l',lwd=8,xlab='DIC mg L-1',cex.lab=cex.lab,
     cex.axis=cex.axis,ylab='Cumulaive Fraction',col='grey60')
sum$dic_conc<-sum$DIC_epi/sum$Vepi*12
sum<-sum[sort.list(sum$dic_conc),]
sum$i=seq(1,length(sum$Area),1)
sum$i<-sum$i/length(sum$Permanent_)
lines(sum$i~sum$dic_conc,lwd=8)
# axis(1,at=c(log10(1),log10(5),log10(10),log10(30),log10(100)),labels=c(1,5,10,30,100),cex.axis=2)
legend(x = 11,y = .85,legend = c('Modeled','Random Sampling'),fill = c('black','grey60'),bty='n',cex = 1.5)
dev.off()

png('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Figures/hanson_dic_1to1.png',
    width = 7,height = 7,res = 300,units = 'in')
par(bty='l',mar=c(5,6,5,5))
temp=merge(hanson,sum,by='Permanent_',all.x=T)
temp=temp[!is.na(temp$dic_conc[temp$dic>0]),]
plot(temp$dic[temp$dic>0]~temp$dic_conc[temp$dic>0],pch=16,cex=2,xlim=c(0,20),ylim=c(0,20),cex.axis=2,
     cex.lab=2,ylab='Observed (mg C L-1)',xlab='Modeled (mg C L-1)')
abline(0,1,lwd=3,lty=2)
text(x = 17,y = 14,labels = '1:1',cex = 2)
dev.off()
# abline(lm(temp$dic[temp$dic>0]~temp$dic_conc[temp$dic>0]))
# res=resid(lm(temp$dic[temp$dic>0]~temp$dic_conc[temp$dic>0]))
# plot(res[temp$dic>0&temp$dic<2]~temp$elevation[temp$dic>0&temp$dic<2],pch=16,cex=2)
summary(temp$dic[temp$dic>0])
summary(temp$dic_conc)



temp=merge(hanson,sum,by='Permanent_',all.x=T)
plot(temp$doc[temp$doc>0]~temp$doc_conc[temp$doc>0],pch=16,cex=2,xlim=c(0,32),ylim=c(0,32),cex.axis=2,
     cex.lab=2,ylab='Observed (mg C L-1)',xlab='Modeled (mg C L-1)') # see Canham et al 2004 Fig. 2 
abline(0,1,lwd=3,lty=2)
dev.off()

windows()
hist((sum$Sed_tPOC+sum$Sed_phyto)/sum$Area*12*1000)

windows()
hist((sum$Burial_tPOC+sum$Burial_phyto)/(sum$Sed_tPOC+sum$Sed_phyto))

windows()
hist(sum$Emit/sum$Area*12*1000)
summary(sum$Emit/sum$Area*12*1000)

windows()
hist(sum$GPP*12*sum$zmix*1000)
summary(sum$GPP*12*sum$zmix*1000)

windows()
hist((sum$DOCl_hypo+sum$DOCr_hypo)/sum$Vhypo*12)
summary((sum$DOCl_hypo+sum$DOCr_hypo)/sum$Vhypo*12)

windows()
hist(sum$epiTemp)

windows()
hist(sum$P_epi/sum$Vepi*31*1000)
hanson<-hanson[sort.list(hanson$totpf),]
hanson$i<-seq(1,length(hanson$wbic),1)/length(hanson$wbic)
plot(hanson$i[hanson$totpf>0]~hanson$totpf[hanson$totpf>0],type='l',lwd=8,xlab='SRP ug L-1',xaxt='n',cex.lab=1.5,xlim=c(0,120),
     cex.axis=2,ylab='Cumulaive Fraction',col='grey60')
axis(1,at=c(log10(1),log10(5),log10(10),log10(30),log10(100)),labels=c(1,5,10,30,100),cex.axis=2)
legend('topleft',legend = c('Modeled','Random Sampling'),fill = c('black','grey60'),bty='n' )
sum$p_conc<-(sum$P_epi)/(sum$Vepi)*31*1000
sum<-sum[sort.list(sum$p_conc),]
sum$i<-seq(1,length(sum$SWin),1)
sum$i<-sum$i/length(sum$SWin)
lines(sum$i~sum$p_conc,lwd=8)
temp=merge(hanson,sum,by='Permanent_',all.x=T)
plot(temp$totpf[temp$totpf>0]~temp$p_conc[temp$totpf>0],pch=16,cex=2,xlim=c(0,70),ylim=c(0,70),cex.axis=2,
     cex.lab=2,ylab='Observed',xlab='Modeled')
abline(0,1,lwd=2,lty=2)

summary(sum$p_conc)
summary(hanson$totpf)

windows()
hist((sum$P_epi+sum$phyto/95)/sum$Vepi*31*1000)
png('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Figures/hanson_tp_cumFreq.png',
    width = 7,height = 7,res = 300,units = 'in')
par(mar=c(6,6,5,2))
cex.axis=2.5
cex.lab=2.5
cex=2
lwd=8
hanson<-hanson[sort.list(hanson$totpuf),]
hanson$i<-seq(1,length(hanson$wbic),1)/length(hanson$wbic)
plot(hanson$i[hanson$totpuf>0]~hanson$totpuf[hanson$totpuf>0],type='l',lwd=8,xlab='TP ug L-1',cex.lab=cex.lab,xlim=c(0,120),
     cex.axis=cex.axis,ylab='Cumulaive Fraction',col='grey60')
# axis(1,at=c(log10(5),log10(10),log10(30),log10(100)),labels=c(5,10,30,100),cex.axis=2)
legend(x=50,y=.85,legend = c('Modeled','Random Sampling'),fill = c('black','grey60'),bty='n' ,cex = 1.5)
sum$p_conc<-(sum$P_epi+sum$phyto/95)/sum$Vepi*31*1000
sum<-sum[sort.list(sum$p_conc),]
sum$i<-seq(1,length(sum$SWin),1)
sum$i<-sum$i/length(sum$SWin)
lines(sum$i~sum$p_conc,lwd=8)
dev.off()

png('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Figures/hanson_tp_1to1.png',
    width = 7,height = 7,res = 300,units = 'in')
par(bty='l',mar=c(5,6,5,5))
temp=merge(hanson,sum,by='Permanent_',all.x=T)
plot(temp$totpuf[temp$totpuf>0]~temp$p_conc[temp$totpuf>0],pch=16,cex=2,xlim=c(0,70),ylim=c(0,70),cex.axis=2,
     cex.lab=2,ylab='Observed (ug P L-1)',xlab='Modeled (ug P L-1)')
abline(0,1,lwd=2,lty=2)
text(x = 60,y = 50,labels = '1:1',cex = 2)
dev.off()

summary(sum$p_conc)
summary(hanson$totpuf[hanson$totpuf>0])

windows()
png('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Figures/hanson_alk_cumFreq.png',
    width = 7,height = 7,res = 300,units = 'in')
par(mar=c(6,6,5,2))
cex.axis=2.5
cex.lab=2.5
cex=2
lwd=8
hanson<-hanson[sort.list(hanson$alk),]
hanson$i<-seq(1,length(hanson$wbic),1)/length(hanson$wbic)
plot(hanson$i~hanson$alk,type='l',lwd=8,xlab='Alkalinity (uEq L-1)',cex.lab=cex.lab,xlim=c(-80,1500),
     cex.axis=cex.axis,ylab='Cumulaive Fraction',col='grey60')
legend(x=600,y=.85,legend = c('Modeled','Random Sampling'),fill = c('black','grey60'),bty='n' ,cex=1.5)
sum<-sum[sort.list(sum$alk),]
sum$i<-seq(1,length(sum$SWin),1)
sum$i<-sum$i/length(sum$SWin)
lines(sum$i~sum$alk,lwd=8)
dev.off()

png('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Figures/hanson_alk_1to1.png',
    width = 7,height = 7,res = 300,units = 'in')
par(bty='l',mar=c(5,6,5,5))
temp=merge(hanson,sum,by='Permanent_',all.x=T)
plot(temp$alk.x~temp$alk.y,pch=16,cex=2,cex.axis=2,ylim=c(-80,1500),xlim=c(-80,1500),
     cex.lab=2,ylab='Observed (uEq L-1)',xlab='Modeled (uEq L-1)')
abline(0,1,lwd=3,lty=2)
text(x = 1200,y = 1000,labels = '1:1',cex = 2)
dev.off()

summary(hanson$alk)
summary(sum$alk)

windows()
hanson<-hanson[sort.list(hanson$pH),]
hanson$i<-seq(1,length(hanson$wbic),1)/length(hanson$wbic)
plot(hanson$i~hanson$pH,type='l',lwd=8,xlab='pH',cex.lab=1.5,
     cex.axis=2,ylab='Cumulaive Fraction',col='grey60')
legend('topleft',legend = c('Modeled','Random Sampling'),fill = c('black','grey60'),bty='n' )
sum<-sum[sort.list(sum$pH),]
sum$i<-seq(1,length(sum$SWin),1)
sum$i<-sum$i/length(sum$SWin)
lines(sum$i~sum$pH,lwd=8)
temp=merge(hanson,sum,by='Permanent_',all.x=T)
plot(temp$pH.x~temp$pH.y,pch=16,cex=2,cex.axis=2,ylim=c(4,9),xlim=c(4,9),
     cex.lab=2,ylab='Observed',xlab='Modeled')
abline(0,1,lwd=2,lty=2)
abline(lm(temp$pH.x~temp$pH.y))
summary(sum$pH)
summary(hanson$pH)

summary(sum$pH~sum$lakeSizeBins)
summary(hanson$pH~hanson$area)

windows()
hist(sum$zmix)

windows()
hist(sum$tPOC_epi/sum$Vepi*12)

windows()
hist(sum$phyto/sum$Vepi*12*1000/50)

windows()
hanson<-hanson[sort.list(hanson$watershed_area),]
hanson$i<-seq(1,length(hanson$wbic),1)/length(hanson$wbic)
plot(hanson$i~log10(hanson$watershed_area),type='l',lwd=8,xlab='WA (m2)',xaxt='n',cex.lab=1.5,
     cex.axis=2,ylab='Cumulaive Fraction',col='grey60')
axis(1,at=c(log10(1e4),log10(1e6),log10(1e8)),labels=c(1e4,1e6,1e8),cex.axis=2)
watersheds<-read.table('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Data/C model forcing data/NHLDsheds_20170323.txt',
                       stringsAsFactors = F,header=T,sep = '\t')
area=sum[,c('Permanent_','Area')]
watersheds<-merge(watersheds,area,all.y = T)
watersheds$WALA<-(watersheds$Area_m2-watersheds$Area)/watersheds$Area
watersheds<-watersheds[sort.list(watersheds$Area_m2),]
watersheds$i<-seq(1,length(watersheds$FID_),1)
watersheds$i<-watersheds$i/length(watersheds$FID_)
lines(watersheds$i~log10(watersheds$Area_m2),lwd=8)
legend('topleft',legend = c('Modeled','Random Sampling'),fill = c('black','grey60'),bty='n' )
summary(watersheds$Area_m2)
summary(hanson$watershed_area)

windows()
hanson<-hanson[sort.list(hanson$area),]
hanson$i<-seq(1,length(hanson$wbic),1)/length(hanson$wbic)
plot(hanson$i~log10(hanson$area),type='l',lwd=8,xlab='LA (m2)',xaxt='n',cex.lab=1.5,
     cex.axis=2,ylab='Cumulaive Fraction',col='grey60')
axis(1,at=c(log10(1e4),log10(1e6),log10(1e8)),labels=c(1e4,1e6,1e8),cex.axis=2)
watersheds<-read.table('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Data/C model forcing data/NHLDsheds_20170323.txt',
                       stringsAsFactors = F,header=T,sep = '\t')
area=sum[,c('Permanent_','Area')]
watersheds<-merge(watersheds,area,all.y = T)
watersheds<-watersheds[sort.list(watersheds$Area),]
watersheds$i<-seq(1,length(watersheds$FID_),1)
watersheds$i<-watersheds$i/length(watersheds$FID_)
lines(watersheds$i~log10(watersheds$Area),lwd=8)
legend('topleft',legend = c('Modeled','Random Sampling'),fill = c('black','grey60'),bty='n' )
summary(watersheds$Area)
summary(hanson$area)

windows()
hanson<-hanson[sort.list(hanson$watershed_area/hanson$area),]
hanson$i<-seq(1,length(hanson$wbic),1)/length(hanson$wbic)
plot(hanson$i~log10(hanson$watershed_area/hanson$area),type='l',lwd=8,xlab='WA:LA',cex.lab=1.5,xlim=c(log10(.1),log10(600)),xaxt='n',
     cex.axis=2,ylab='Cumulaive Fraction',col='grey60')
axis(1,at=c(log10(1),log10(10),log10(100)),labels=c(1,10,100),cex.axis=2)
watersheds<-read.table('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Data/C model forcing data/NHLDsheds_20170323.txt',
                       stringsAsFactors = F,header=T,sep = '\t')
area=sum[,c('Permanent_','Area')]
watersheds<-merge(watersheds,area,all.y = T)
watersheds$WALA<-(watersheds$Area_m2-watersheds$Area)/watersheds$Area
watersheds<-watersheds[sort.list(watersheds$WALA),]
watersheds$i<-seq(1,length(watersheds$FID_),1)
watersheds$i<-watersheds$i/length(watersheds$FID_)
lines(watersheds$i~log10(watersheds$WALA),lwd=8)
legend('topleft',legend = c('Modeled','Random Sampling'),fill = c('black','grey60'),bty='n' )
summary(watersheds$WALA)
summary(hanson$watershed_area/hanson$area)

windows()
hanson<-hanson[sort.list(hanson$elevation),]
hanson$i<-seq(1,length(hanson$wbic),1)/length(hanson$wbic)
plot(hanson$i~hanson$elevation,type='l',lwd=8,xlab='Elevation (m)',cex.lab=1.5,
     cex.axis=2,ylab='Cumulaive Fraction',col='grey60')
sum<-sum[sort.list(sum$Elev),]
sum$i<-seq(1,length(sum$SWin),1)
sum$i<-sum$i/length(sum$SWin)
lines(sum$i~sum$Elev,lwd=8)
legend('topleft',legend = c('Modeled','Random Sampling'),fill = c('black','grey60'),bty='n' )


windows()
hanson<-hanson[sort.list(hanson$depth_at_samplesite),]
hanson$i<-seq(1,length(hanson$wbic),1)/length(hanson$wbic)
plot(hanson$i~hanson$depth_at_samplesite,type='l',lwd=8,xlab='Obs Max Depth (m)',cex.lab=1.5,
     cex.axis=2,ylab='Cumulaive Fraction',col='grey60')
sum<-sum[sort.list(sum$Stage),]
sum$i<-seq(1,length(sum$SWin),1)
sum$i<-sum$i/length(sum$SWin)
lines(sum$i~sum$Stage,lwd=8)
legend('topleft',legend = c('Modeled','Random Sampling'),fill = c('black','grey60'),bty='n' )
temp=merge(hanson,sum,by='Permanent_',all.x=T)
plot(temp$depth_at_samplesite[temp$depth_at_samplesite>0]~temp$Stage[temp$depth_at_samplesite>0],
     pch=16,cex=2,xlim=c(0,25),ylim=c(0,25),cex.axis=2,
     cex.lab=2,ylab='Observed',xlab='Modeled')
abline(0,1,lwd=2,lty=2)
summary(sum$Stage)
summary(hanson$depth_at_samplesite)

windows()
hanson2=hanson
hanson2$top_hypolimnion<-ifelse(hanson2$top_hypolimnion<0,NA,hanson2$top_hypolimnion)
hanson2=hanson2[!is.na(hanson2$top_hypolimnion),]
hanson2<-hanson2[sort.list(hanson2$top_hypolimnion),]
hanson2$i<-seq(1,length(hanson2$wbic),1)/length(hanson2$wbic)
plot(hanson2$i~hanson2$top_hypolimnion,type='l',lwd=8,xlab='Top (m)',cex.lab=1.5,
     cex.axis=2,ylab='Cumulaive Fraction',col='grey60')
sum<-sum[sort.list(sum$zmix),]
sum$i<-seq(1,length(sum$SWin),1)
sum$i<-sum$i/length(sum$SWin)
lines(sum$i~sum$zmix,lwd=8)
legend('topleft',legend = c('Modeled','Random Sampling'),fill = c('black','grey60'),bty='n' )
temp=merge(hanson,sum,by='Permanent_',all.x=T)
plot(temp$top_hypolimnion[temp$top_hypolimnion>0]~temp$zmix[temp$top_hypolimnion>0],pch=16,cex=2,
     xlim=c(0,14),ylim=c(0,14),cex.axis=2,
     cex.lab=2,ylab='Observed',xlab='Modeled')
abline(0,1,lwd=2,lty=2)
summary(hanson2$area)
summary(sum$Area)



# LTER NTL validation chem data 
ntl<-read.csv('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Validation Data/NTL LTER Chemistry/chemical_limnology_of_north_temperate_lakes_lter_primary_study_lakes__nutrients_ph_and_carbon.csv',
              stringsAsFactors = F)
ntlChl<-read.csv('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Validation Data/NTL LTER Chl/north_temperate_lakes_lter__chlorophyll_-_trout_lake_area.csv',
                 stringsAsFactors = F)
ntlTemp<-read.csv('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Validation Data/physical_limnology_of_the_north_temperate_lakes_primary_study_lakes.csv',
                  stringsAsFactors = F)
ntlLookUp<-data.frame(lakeID=c('AL','BM','CR','SP','TR','TB','CB'),Permanent_=c(69886156,69886284,69886510,69886444,69886228,69886158,123148117))
ntl<-ntl[ntl$lakeid%in%ntlLookUp$lakeID,]
ntl$sampledate<-as.POSIXct(ntl$sampledate,format='%Y-%m-%d')
ntlChl<-ntlChl[ntlChl$lakeid%in%ntlLookUp$lakeID,]
ntlChl$sampledate<-as.POSIXct(ntlChl$sampledate, format='%Y-%m-%d')
ntlTemp<-ntlTemp[ntlTemp$lakeid%in%ntlLookUp$lakeID,]
ntlTemp$sampledate<-as.POSIXct(ntlTemp$sampledate, format='%Y-%m-%d')
ntlTemp<-ntlTemp[sort.list(ntlTemp$sampledate),]

dir<-'D:/Jake/My Papers/NHLD Carbon Model/Results/20170818/'
skip=6*365 # days to skip 
for(i in 1:length(ntlLookUp$lakeID)){
    if(ntlLookUp$lakeID[i]=='CB'){
      cur<-read.table(file.path(dir,paste(ntlLookUp$Permanent_[i],'_C_model.txt',sep='')),header=T,sep='\t')
    }else{
      cur<-read.table(file.path(dir,paste(ntlLookUp$Permanent_[i],'_CRWA_V_adj_C_model.txt',sep='')),
                      header=T,sep='\t')
    }
    cur<-cur[skip:nrow(cur),]
    assign(paste(as.character(ntlLookUp$lakeID[i]),'Cout',sep=''),cur)
}


###################################
ntlMod_Obs=data.frame(lakeID=rep(NA,7),DOC_mod=rep(NA,7),DOC_obs=rep(NA,7),DOC_cor=rep(NA,7),SRP_mod=rep(NA,7),SRP_obs=rep(NA,7),
                      SRP_cor=rep(NA,7),TP_mod=rep(NA,7),TP_obs=rep(NA,7),TP_cor=rep(NA,7),pH_mod=rep(NA,7),pH_obs=rep(NA,7),pH_cor=rep(NA,7),DIC_mod=rep(NA,7),
                      DIC_obs=rep(NA,7),DIC_cor=rep(NA,7),Chl_mod=rep(NA,7),Chl_obs=rep(NA,7),Chl_cor=rep(NA,7),Wtr_mod=rep(NA,7),Wtr_obs=rep(NA,7),Wtr_cor=rep(NA,7))
# TR ****************************************************
# windows()
png('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Figures/NTL_validation_all.png',
    width=21,height = 21,res=300,units='in')
xlim=c(as.POSIXct('1987-01-01 00:00:00'),max(as.POSIXct(TRCout$datetime)))
ylim=c(2,7)
docYlim=c(0,22)
srpYlim=c(0,27)
tpYlim=c(0,41)
phYlim=c(3,9)
dicYlim=c(0,14)
chlYlim=c(0,21)
wtrYlim=c(0,25)
par(mfrow=c(7,7),mar=c(5,6,2,2))
cex=2
lwd=2
cex.axis=2
cex.lab=2
# plot(ntl$doc[ntl$lakeid=='TR'&ntl$depth==0&!is.na(ntl$doc)&ntl$doc>ylim[1]&ntl$doc<ylim[2]&ntl$sampledate<xlim[2]]~ntl$sampledate[ntl$lakeid=='TR'&ntl$depth==0&!is.na(ntl$doc)&ntl$doc>ylim[1]&ntl$doc<ylim[2]&ntl$sampledate<xlim[2]],
#      col='grey60',xlim=xlim,ylim=ylim,pch=16,lwd=3,ylab='',xlab='',yaxt='n',xaxt='n',type='o',cex.lab=1.5)
# par(new=T)
# plot((TRCout$DOCr_epi+TRCout$DOCl_epi)/(TRCout$Vepi)*12~as.POSIXct(TRCout$datetime,format='%Y-%m-%d'),cex.lab=1.5,xlim=xlim,ylim=ylim,ylab='DOC mg L-1',type='l',lwd=4,xlab='Date',cex.axis=2,xaxt='n')
# axis(1,at=c(as.POSIXct('1990-01-01'),as.POSIXct('2005-01-01'),as.POSIXct('2010-01-01')),labels = c(1990,2005,2010),cex.axis=2)
# legend('topleft',legend = c('Modeled','Observed'),fill = c('black','grey60'),bty='n')

temp=TRCout
tempntl=ntl[ntl$lakeid=='TR'&ntl$depth<=1&!is.na(ntl$doc)&ntl$doc>ylim[1]&ntl$doc<ylim[2]&ntl$sampledate<xlim[2],]
temp$month<-strftime(strptime(temp$datetime,'%Y-%m-%d'),'%m')
tempntl$month<-strftime(strptime(tempntl$sampledate,'%Y-%m-%d'),'%m')
temp=aggregate((temp$DOCr_epi+temp$DOCl_epi)/temp$Vepi*12,by=list(temp$month),FUN=mean)
tempntl=aggregate(tempntl$doc,by=list(tempntl$month),FUN=mean)
colnames(temp)<-c('month','doc')
colnames(tempntl)<-c('month','doc')
# windows()
ylim=range(temp$doc,tempntl$doc)
plot(temp$doc~temp$month,pch=16,type='o',ylim=docYlim,cex=cex,lwd=lwd,cex.axis=cex.axis,cex.lab=cex.lab,ylab='DOC (mg L-1)',
     xlab='')
points(tempntl$doc~tempntl$month,col='grey60',pch=16,type='o',cex=cex,lwd=lwd)

ntlMod_Obs$lakeID[1]='TR'
ntlMod_Obs$DOC_mod[1]=mean(temp$doc)
ntlMod_Obs$DOC_obs[1]=mean(tempntl$doc)
ntlMod_Obs$DOC_cor[1]=cor(temp$doc[temp$month%in%tempntl$month],tempntl$doc)

# windows()
# xlim=c(as.POSIXct('1987-01-01 00:00:00'),max(as.POSIXct(TRCout$datetime)))
ylim=c(0,20)
# plot(ntl$totpf[ntl$lakeid=='TR'&ntl$depth==0&!is.na(ntl$totpf)&ntl$totpf>ylim[1]&ntl$totpf<ylim[2]&ntl$sampledate<xlim[2]]~
#        ntl$sampledate[ntl$lakeid=='TR'&ntl$depth==0&!is.na(ntl$totpf)&ntl$totpf>ylim[1]&ntl$totpf<ylim[2]&ntl$sampledate<xlim[2]],
#      col='grey60',xlim=xlim,ylim=ylim,pch=16,lwd=3,ylab='',xlab='',yaxt='n',xaxt='n',type='o',cex.lab=1.5)
# par(new=T)
# plot(TRCout$P_epi/(TRCout$Vepi)*31*1000~as.POSIXct(TRCout$datetime,format='%Y-%m-%d'),cex.lab=1.5,xlim=xlim,ylim=ylim,ylab='SRP ug L-1',type='l',lwd=4,xlab='Date',cex.axis=2,xaxt='n')
# axis(1,at=c(as.POSIXct('1990-01-01'),as.POSIXct('2005-01-01'),as.POSIXct('2010-01-01')),labels = c(1990,2005,2010),cex.axis=2)
# legend('topleft',legend = c('Modeled','Observed'),fill = c('black','grey60'),bty='n')

temp=TRCout
tempntl=ntl[ntl$lakeid=='TR'&ntl$depth<=1&!is.na(ntl$totpf)&ntl$totpf>ylim[1]&ntl$totpf<ylim[2]&ntl$sampledate<xlim[2],]
temp$month<-strftime(strptime(temp$datetime,'%Y-%m-%d'),'%m')
tempntl$month<-strftime(strptime(tempntl$sampledate,'%Y-%m-%d'),'%m')
temp=aggregate(temp$P_epi/(temp$Vepi)*31*1000,by=list(temp$month),FUN=mean)
tempntl=aggregate(tempntl$totpf,by=list(tempntl$month),FUN=mean)
colnames(temp)<-c('month','srp')
colnames(tempntl)<-c('month','srp')
# windows()
ylim=range(temp$srp,tempntl$srp)
plot(temp$srp~temp$month,pch=16,type='o',ylim=srpYlim,cex=cex,lwd=lwd,cex.axis=cex.axis,cex.lab=cex.lab,ylab='SRP (ug L-1)',
     xlab='')
points(tempntl$srp~tempntl$month,col='grey60',pch=16,type='o',cex=cex,lwd=lwd)

ntlMod_Obs$SRP_mod[1]=mean(temp$srp)
ntlMod_Obs$SRP_obs[1]=mean(tempntl$srp)
ntlMod_Obs$SRP_cor[1]=cor(temp$srp[temp$month%in%tempntl$month],tempntl$srp)

# windows()
# xlim=c(as.POSIXct('1987-01-01 00:00:00'),max(as.POSIXct(TRCout$datetime)))
ylim=c(0,40)
# plot(ntl$totpuf[ntl$lakeid=='TR'&ntl$depth==0&!is.na(ntl$totpuf)&ntl$totpuf>ylim[1]&ntl$totpuf<ylim[2]&ntl$sampledate<xlim[2]]~
#        ntl$sampledate[ntl$lakeid=='TR'&ntl$depth==0&!is.na(ntl$totpuf)&ntl$totpuf>ylim[1]&ntl$totpuf<ylim[2]&ntl$sampledate<xlim[2]],
#      col='grey60',xlim=xlim,ylim=ylim,pch=16,lwd=3,ylab='',xlab='',yaxt='n',xaxt='n',type='o',cex.lab=1.5)
# par(new=T)
# plot((TRCout$P_epi+TRCout$phyto/95)/(TRCout$Vepi)*31*1000~as.POSIXct(TRCout$datetime,format='%Y-%m-%d'),cex.lab=1.5,xlim=xlim,ylim=ylim,ylab='TP ug L-1',type='l',lwd=4,xlab='Date',cex.axis=2,xaxt='n')
# axis(1,at=c(as.POSIXct('1990-01-01'),as.POSIXct('2005-01-01'),as.POSIXct('2010-01-01')),labels = c(1990,2005,2010),cex.axis=2)
# legend('topleft',legend = c('Modeled','Observed'),fill = c('black','grey60'),bty='n')

temp=TRCout
tempntl=ntl[ntl$lakeid=='TR'&ntl$depth<=1&!is.na(ntl$totpuf)&ntl$totpuf>ylim[1]&ntl$totpuf<ylim[2]&ntl$sampledate<xlim[2],]
temp$month<-strftime(strptime(temp$datetime,'%Y-%m-%d'),'%m')
tempntl$month<-strftime(strptime(tempntl$sampledate,'%Y-%m-%d'),'%m')
temp=aggregate((temp$P_epi+temp$phyto/95)/(temp$Vepi)*31*1000,by=list(temp$month),FUN=mean)
tempntl=aggregate(tempntl$totpuf,by=list(tempntl$month),FUN=mean)
colnames(temp)<-c('month','tp')
colnames(tempntl)<-c('month','tp')
# windows()
ylim=range(temp$tp,tempntl$tp)
plot(temp$tp~temp$month,pch=16,type='o',ylim=tpYlim,cex=cex,lwd=lwd,cex.axis=cex.axis,cex.lab=cex.lab,ylab='TP (ug L-1)',
     xlab='')
points(tempntl$tp~tempntl$month,col='grey60',pch=16,type='o',cex=cex,lwd=lwd)

ntlMod_Obs$TP_mod[1]=mean(temp$tp)
ntlMod_Obs$TP_obs[1]=mean(tempntl$tp)
ntlMod_Obs$TP_cor[1]=cor(temp$tp[temp$month%in%tempntl$month],tempntl$tp)

# windows()
# xlim=c(as.POSIXct('1987-01-01 00:00:00'),max(as.POSIXct(TRCout$datetime)))
ylim=c(6,10)
# plot(ntl$ph[ntl$lakeid=='TR'&ntl$depth==0&!is.na(ntl$ph)&ntl$ph>ylim[1]&ntl$ph<ylim[2]&ntl$sampledate<xlim[2]]~
#        ntl$sampledate[ntl$lakeid=='TR'&ntl$depth==0&!is.na(ntl$ph)&ntl$ph>ylim[1]&ntl$ph<ylim[2]&ntl$sampledate<xlim[2]],
#      col='grey60',xlim=xlim,ylim=ylim,pch=16,lwd=3,ylab='',xlab='',yaxt='n',xaxt='n',type='o',cex.lab=1.5)
# par(new=T)
# plot(TRCout$pH~as.POSIXct(TRCout$datetime,format='%Y-%m-%d'),cex.lab=1.5,xlim=xlim,ylim=ylim,ylab='pH',type='l',lwd=4,xlab='Date',cex.axis=2,xaxt='n')
# axis(1,at=c(as.POSIXct('1990-01-01'),as.POSIXct('2005-01-01'),as.POSIXct('2010-01-01')),labels = c(1990,2005,2010),cex.axis=2)
# legend('topleft',legend = c('Modeled','Observed'),fill = c('black','grey60'),bty='n')

temp=TRCout
tempntl=ntl[ntl$lakeid=='TR'&ntl$depth<=4&!is.na(ntl$ph)&ntl$ph>ylim[1]&ntl$ph<ylim[2]&ntl$sampledate<xlim[2],]
temp$month<-strftime(strptime(temp$datetime,'%Y-%m-%d'),'%m')
tempntl$month<-strftime(strptime(tempntl$sampledate,'%Y-%m-%d'),'%m')
temp=aggregate(temp$pH,by=list(temp$month),FUN=mean,na.rm=T)
tempntl=aggregate(tempntl$ph,by=list(tempntl$month),FUN=mean)
colnames(temp)<-c('month','pH')
colnames(tempntl)<-c('month','pH')
# windows()
ylim=range(temp$pH,tempntl$pH)
plot(temp$pH~temp$month,pch=16,type='o',ylim=phYlim,cex=cex,lwd=lwd,cex.axis=cex.axis,cex.lab=cex.lab,ylab='pH',
     xlab='')
points(tempntl$pH~tempntl$month,col='grey60',pch=16,type='o',cex=cex,lwd=lwd)

ntlMod_Obs$pH_mod[1]=mean(temp$pH)
ntlMod_Obs$pH_obs[1]=mean(tempntl$pH)
ntlMod_Obs$pH_cor[1]=cor(temp$pH[temp$month%in%tempntl$month],tempntl$pH)

# windows()
# xlim=c(as.POSIXct('1987-01-01 00:00:00'),max(as.POSIXct(TRCout$datetime)))
ylim=c(0,20)
# plot(ntl$dic[ntl$lakeid=='TR'&ntl$depth==0&!is.na(ntl$dic)&ntl$dic>ylim[1]&ntl$dic<ylim[2]&ntl$sampledate<xlim[2]]~
#        ntl$sampledate[ntl$lakeid=='TR'&ntl$depth==0&!is.na(ntl$dic)&ntl$dic>ylim[1]&ntl$dic<ylim[2]&ntl$sampledate<xlim[2]],
#      col='grey60',xlim=xlim,ylim=ylim,pch=16,lwd=3,ylab='',xlab='',yaxt='n',xaxt='n',type='o',cex.lab=1.5)
# par(new=T)
# plot(TRCout$DIC_epi/TRCout$Vepi*12~as.POSIXct(TRCout$datetime,format='%Y-%m-%d'),cex.lab=1.5,xlim=xlim,ylim=ylim,ylab='DIC mg L-1',type='l',lwd=4,xlab='Date',cex.axis=2,xaxt='n')
# axis(1,at=c(as.POSIXct('1990-01-01'),as.POSIXct('2005-01-01'),as.POSIXct('2010-01-01')),labels = c(1990,2005,2010),cex.axis=2)
# legend('topleft',legend = c('Modeled','Observed'),fill = c('black','grey60'),bty='n')

temp=TRCout
tempntl=ntl[ntl$lakeid=='TR'&ntl$depth<=4&!is.na(ntl$dic)&ntl$dic>ylim[1]&ntl$dic<ylim[2]&ntl$sampledate<xlim[2],]
temp$month<-strftime(strptime(temp$datetime,'%Y-%m-%d'),'%m')
tempntl$month<-strftime(strptime(tempntl$sampledate,'%Y-%m-%d'),'%m')
temp=aggregate(temp$DIC_epi/temp$Vepi*12,by=list(temp$month),FUN=mean,na.rm=T)
tempntl=aggregate(tempntl$dic,by=list(tempntl$month),FUN=mean)
colnames(temp)<-c('month','dic')
colnames(tempntl)<-c('month','dic')
# windows()
ylim=range(temp$dic,tempntl$dic)
plot(temp$dic~temp$month,pch=16,type='o',ylim=dicYlim,cex=cex,lwd=lwd,cex.axis=cex.axis,cex.lab=cex.lab,ylab='DIC (mg L-1)',
     xlab='')
points(tempntl$dic~tempntl$month,col='grey60',pch=16,type='o',cex=cex,lwd=lwd)

ntlMod_Obs$DIC_mod[1]=mean(temp$dic)
ntlMod_Obs$DIC_obs[1]=mean(tempntl$dic)
ntlMod_Obs$DIC_cor[1]=cor(temp$dic[temp$month%in%tempntl$month],tempntl$dic)

# windows()
# xlim=c(as.POSIXct('1987-01-01 00:00:00'),max(as.POSIXct(TRCout$datetime)))
ylim=c(0,10)
# plot(ntlChl$chlor[ntlChl$lakeid=='TR'&ntlChl$depth==0&!is.na(ntlChl$chlor)&ntlChl$chlor>ylim[1]&ntlChl$chlor<ylim[2]&ntlChl$sampledate<xlim[2]]~
#        ntlChl$sampledate[ntlChl$lakeid=='TR'&ntlChl$depth==0&!is.na(ntlChl$chlor)&ntlChl$chlor>ylim[1]&ntlChl$chlor<ylim[2]&ntlChl$sampledate<xlim[2]],
#      col='grey60',xlim=xlim,ylim=ylim,pch=16,lwd=3,ylab='',xlab='',yaxt='n',xaxt='n',type='o',cex.lab=1.5)
# par(new=T)
# plot((TRCout$phyto/50)/(TRCout$Vepi)*12*1000~as.POSIXct(TRCout$datetime,format='%Y-%m-%d'),cex.lab=1.5,xlim=xlim,ylim=ylim,ylab='Chl-a ug L-1',type='l',lwd=4,xlab='Date',cex.axis=2,xaxt='n')
# axis(1,at=c(as.POSIXct('1990-01-01'),as.POSIXct('2005-01-01'),as.POSIXct('2010-01-01')),labels = c(1990,2005,2010),cex.axis=2)
# legend('topleft',legend = c('Modeled','Observed'),fill = c('black','grey60'),bty='n')

temp=TRCout
tempntl=ntlChl[ntlChl$lakeid=='TR'&ntlChl$depth<=4&!is.na(ntlChl$chlor)&ntlChl$chlor>ylim[1]&ntlChl$chlor<ylim[2]&ntlChl$sampledate<xlim[2],]
temp$month<-strftime(strptime(temp$datetime,'%Y-%m-%d'),'%m')
tempntl$month<-strftime(strptime(tempntl$sampledate,'%Y-%m-%d'),'%m')
temp=aggregate((temp$phyto/50)/(temp$Vepi)*12*1000,by=list(temp$month),FUN=mean)
tempntl=aggregate(tempntl$chlor,by=list(tempntl$month),FUN=mean)
colnames(temp)<-c('month','chl')
colnames(tempntl)<-c('month','chl')
# windows()
ylim=range(temp$chl,tempntl$chl)
plot(temp$chl~temp$month,pch=16,type='o',ylim=chlYlim,cex=cex,lwd=lwd,cex.axis=cex.axis,cex.lab=cex.lab,ylab='Chl (ug L-1)',
     xlab='')
points(tempntl$chl~tempntl$month,col='grey60',pch=16,type='o',cex=cex,lwd=lwd)

ntlMod_Obs$Chl_mod[1]=mean(temp$chl)
ntlMod_Obs$Chl_obs[1]=mean(tempntl$chl)
ntlMod_Obs$Chl_cor[1]=cor(temp$chl[temp$month%in%tempntl$month],tempntl$chl)

# windows()
# xlim=c(as.POSIXct('1987-01-01 00:00:00'),max(as.POSIXct(TRCout$datetime)))
# plot(ntlTemp$wtemp[ntlTemp$depth==1&ntlTemp$lakeid=='TR']~ntlTemp$sampledate[ntlTemp$depth==1&ntlTemp$lakeid=='TR'],type='l',xlim=xlim,
#      ylab='Temp C',xlab='',lwd=2,cex.lab=1.5,cex.axis=2,col='grey60')
# lines(TRCout$epiTemp~as.POSIXct(TRCout$datetime),lwd=2)

temp=TRCout
tempntl=ntlTemp[ntlTemp$depth==1&ntlTemp$lakeid=='TR'&!is.na(ntlTemp$wtemp),]
temp$month<-strftime(strptime(temp$datetime,'%Y-%m-%d'),'%m')
tempntl$month<-strftime(strptime(tempntl$sampledate,'%Y-%m-%d'),'%m')
temp=aggregate(temp$epiTemp,by=list(temp$month),FUN=mean)
tempntl=aggregate(tempntl$wtemp,by=list(tempntl$month),FUN=mean)
colnames(temp)<-c('month','wtr')
colnames(tempntl)<-c('month','wtr')
# windows()
ylim=range(temp$wtr,tempntl$wtr)
plot(temp$wtr~temp$month,pch=16,type='o',ylim=wtrYlim,cex=cex,lwd=lwd,cex.axis=cex.axis,cex.lab=cex.lab,ylab='Epi Temp (deg C)',
     xlab='')
points(tempntl$wtr~tempntl$month,col='grey60',pch=16,type='o',cex=2)

ntlMod_Obs$Wtr_mod[1]=mean(temp$wtr)
ntlMod_Obs$Wtr_obs[1]=mean(tempntl$wtr)
ntlMod_Obs$Wtr_cor[1]=cor(temp$wtr[temp$month%in%tempntl$month],tempntl$wtr)

# dev.off()

# **********************************************************

# CR ****************************************************
# windows()
# png('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Figures/CR_validation.png',width=49,height = 7,res=300,units='in')
# xlim=c(as.POSIXct('1987-01-01 00:00:00'),max(as.POSIXct(CRCout$datetime)))
ylim=c(0,10)
# par(mfrow=c(1,7),mar=c(5,6,2,2))
cex=2
lwd=2
cex.axis=2
cex.lab=2
# plot(ntl$doc[ntl$lakeid=='CR'&ntl$depth==0&!is.na(ntl$doc)&ntl$doc>ylim[1]&ntl$doc<ylim[2]&ntl$sampledate<xlim[2]]~ntl$sampledate[ntl$lakeid=='CR'&ntl$depth==0&!is.na(ntl$doc)&ntl$doc>ylim[1]&ntl$doc<ylim[2]&ntl$sampledate<xlim[2]],
#      col='grey60',xlim=xlim,ylim=ylim,pch=16,lwd=3,ylab='',xlab='',yaxt='n',xaxt='n',type='o',cex.lab=1.5)
# par(new=T)
# plot((CRCout$DOCr_epi+CRCout$DOCl_epi)/(CRCout$Vepi)*12~as.POSIXct(CRCout$datetime,format='%Y-%m-%d'),cex.lab=1.5,xlim=xlim,ylim=ylim,ylab='DOC mg L-1',type='l',lwd=4,xlab='Date',cex.axis=2,xaxt='n')
# axis(1,at=c(as.POSIXct('1990-01-01'),as.POSIXct('2005-01-01'),as.POSIXct('2010-01-01')),labels = c(1990,2005,2010),cex.axis=2)
# legend('topleft',legend = c('Modeled','Observed'),fill = c('black','grey60'),bty='n')

temp=CRCout
tempntl=ntl[ntl$lakeid=='CR'&ntl$depth<=1&!is.na(ntl$doc)&ntl$doc>ylim[1]&ntl$doc<ylim[2]&ntl$sampledate<xlim[2],]
temp$month<-strftime(strptime(temp$datetime,'%Y-%m-%d'),'%m')
tempntl$month<-strftime(strptime(tempntl$sampledate,'%Y-%m-%d'),'%m')
temp=aggregate((temp$DOCr_epi+temp$DOCl_epi)/temp$Vepi*12,by=list(temp$month),FUN=mean)
tempntl=aggregate(tempntl$doc,by=list(tempntl$month),FUN=mean)
colnames(temp)<-c('month','doc')
colnames(tempntl)<-c('month','doc')
# windows()
ylim=range(temp$doc,tempntl$doc)
plot(temp$doc~temp$month,pch=16,type='o',ylim=docYlim,cex=cex,lwd=lwd,cex.axis=cex.axis,cex.lab=cex.lab,ylab='DOC (mg L-1)',
     xlab='')
points(tempntl$doc~tempntl$month,col='grey60',pch=16,type='o',cex=cex,lwd=lwd)

ntlMod_Obs$lakeID[2]='CR'
ntlMod_Obs$DOC_mod[2]=mean(temp$doc)
ntlMod_Obs$DOC_obs[2]=mean(tempntl$doc)
ntlMod_Obs$DOC_cor[2]=cor(temp$doc[temp$month%in%tempntl$month],tempntl$doc)

# windows()
# xlim=c(as.POSIXct('1987-01-01 00:00:00'),max(as.POSIXct(CRCout$datetime)))
ylim=c(0,10)
# plot(ntl$totpf[ntl$lakeid=='CR'&ntl$depth==0&!is.na(ntl$totpf)&ntl$totpf>ylim[1]&ntl$totpf<ylim[2]&ntl$sampledate<xlim[2]]~
#        ntl$sampledate[ntl$lakeid=='CR'&ntl$depth==0&!is.na(ntl$totpf)&ntl$totpf>ylim[1]&ntl$totpf<ylim[2]&ntl$sampledate<xlim[2]],
#      col='grey60',xlim=xlim,ylim=ylim,pch=16,lwd=3,ylab='',xlab='',yaxt='n',xaxt='n',type='o',cex.lab=1.5)
# par(new=T)
# plot(CRCout$P_epi/(CRCout$Vepi)*31*1000~as.POSIXct(CRCout$datetime,format='%Y-%m-%d'),cex.lab=1.5,xlim=xlim,ylim=ylim,ylab='SRP ug L-1',type='l',lwd=4,xlab='Date',cex.axis=2,xaxt='n')
# axis(1,at=c(as.POSIXct('1990-01-01'),as.POSIXct('2005-01-01'),as.POSIXct('2010-01-01')),labels = c(1990,2005,2010),cex.axis=2)
# legend('topleft',legend = c('Modeled','Observed'),fill = c('black','grey60'),bty='n')

temp=CRCout
tempntl=ntl[ntl$lakeid=='CR'&ntl$depth<=1&!is.na(ntl$totpf)&ntl$totpf>ylim[1]&ntl$totpf<ylim[2]&ntl$sampledate<xlim[2],]
temp$month<-strftime(strptime(temp$datetime,'%Y-%m-%d'),'%m')
tempntl$month<-strftime(strptime(tempntl$sampledate,'%Y-%m-%d'),'%m')
temp=aggregate(temp$P_epi/(temp$Vepi)*31*1000,by=list(temp$month),FUN=mean)
tempntl=aggregate(tempntl$totpf,by=list(tempntl$month),FUN=mean)
colnames(temp)<-c('month','srp')
colnames(tempntl)<-c('month','srp')
# windows()
ylim=range(temp$srp,tempntl$srp)
plot(temp$srp~temp$month,pch=16,type='o',ylim=srpYlim,cex=cex,lwd=lwd,cex.axis=cex.axis,cex.lab=cex.lab,ylab='SRP (ug L-1)',
     xlab='')
points(tempntl$srp~tempntl$month,col='grey60',pch=16,type='o',cex=cex,lwd=lwd)

ntlMod_Obs$SRP_mod[2]=mean(temp$srp)
ntlMod_Obs$SRP_obs[2]=mean(tempntl$srp)
ntlMod_Obs$SRP_cor[2]=cor(temp$srp[temp$month%in%tempntl$month],tempntl$srp)

# windows()
# xlim=c(as.POSIXct('1987-01-01 00:00:00'),max(as.POSIXct(CRCout$datetime)))
ylim=c(0,20)
# plot(ntl$totpuf[ntl$lakeid=='CR'&ntl$depth==0&!is.na(ntl$totpuf)&ntl$totpuf>ylim[1]&ntl$totpuf<ylim[2]&ntl$sampledate<xlim[2]]~
#        ntl$sampledate[ntl$lakeid=='CR'&ntl$depth==0&!is.na(ntl$totpuf)&ntl$totpuf>ylim[1]&ntl$totpuf<ylim[2]&ntl$sampledate<xlim[2]],
#      col='grey60',xlim=xlim,ylim=ylim,pch=16,lwd=3,ylab='',xlab='',yaxt='n',xaxt='n',type='o',cex.lab=1.5)
# par(new=T)
# plot((CRCout$P_epi+CRCout$phyto/95)/(CRCout$Vepi)*31*1000~as.POSIXct(CRCout$datetime,format='%Y-%m-%d'),cex.lab=1.5,xlim=xlim,ylim=ylim,ylab='TP ug L-1',type='l',lwd=4,xlab='Date',cex.axis=2,xaxt='n')
# axis(1,at=c(as.POSIXct('1990-01-01'),as.POSIXct('2005-01-01'),as.POSIXct('2010-01-01')),labels = c(1990,2005,2010),cex.axis=2)
# legend('topleft',legend = c('Modeled','Observed'),fill = c('black','grey60'),bty='n')

temp=CRCout
tempntl=ntl[ntl$lakeid=='CR'&ntl$depth<=1&!is.na(ntl$totpuf)&ntl$totpuf>ylim[1]&ntl$totpuf<ylim[2]&ntl$sampledate<xlim[2],]
temp$month<-strftime(strptime(temp$datetime,'%Y-%m-%d'),'%m')
tempntl$month<-strftime(strptime(tempntl$sampledate,'%Y-%m-%d'),'%m')
temp=aggregate((temp$P_epi+temp$phyto/95)/(temp$Vepi)*31*1000,by=list(temp$month),FUN=mean)
tempntl=aggregate(tempntl$totpuf,by=list(tempntl$month),FUN=mean)
colnames(temp)<-c('month','tp')
colnames(tempntl)<-c('month','tp')
# windows()
ylim=range(temp$tp,tempntl$tp)
plot(temp$tp~temp$month,pch=16,type='o',ylim=tpYlim,cex=cex,lwd=lwd,cex.axis=cex.axis,cex.lab=cex.lab,ylab='TP (ug L-1)',
     xlab='')
points(tempntl$tp~tempntl$month,col='grey60',pch=16,type='o',cex=cex,lwd=lwd)

ntlMod_Obs$TP_mod[2]=mean(temp$tp)
ntlMod_Obs$TP_obs[2]=mean(tempntl$tp)
ntlMod_Obs$TP_cor[2]=cor(temp$tp[temp$month%in%tempntl$month],tempntl$tp)

# windows()
# xlim=c(as.POSIXct('1987-01-01 00:00:00'),max(as.POSIXct(CRCout$datetime)))
ylim=c(3,10)
# plot(ntl$ph[ntl$lakeid=='CR'&ntl$depth==0&!is.na(ntl$ph)&ntl$ph>ylim[1]&ntl$ph<ylim[2]&ntl$sampledate<xlim[2]]~
#        ntl$sampledate[ntl$lakeid=='CR'&ntl$depth==0&!is.na(ntl$ph)&ntl$ph>ylim[1]&ntl$ph<ylim[2]&ntl$sampledate<xlim[2]],
#      col='grey60',xlim=xlim,ylim=ylim,pch=16,lwd=3,ylab='',xlab='',yaxt='n',xaxt='n',type='o',cex.lab=1.5)
# par(new=T)
# plot(CRCout$pH~as.POSIXct(CRCout$datetime,format='%Y-%m-%d'),cex.lab=1.5,xlim=xlim,ylim=ylim,ylab='pH',type='l',lwd=4,xlab='Date',cex.axis=2,xaxt='n')
# axis(1,at=c(as.POSIXct('1990-01-01'),as.POSIXct('2005-01-01'),as.POSIXct('2010-01-01')),labels = c(1990,2005,2010),cex.axis=2)
# legend('topleft',legend = c('Modeled','Observed'),fill = c('black','grey60'),bty='n')

temp=CRCout
tempntl=ntl[ntl$lakeid=='CR'&ntl$depth<=4&!is.na(ntl$ph)&ntl$ph>ylim[1]&ntl$ph<ylim[2]&ntl$sampledate<xlim[2],]
temp$month<-strftime(strptime(temp$datetime,'%Y-%m-%d'),'%m')
tempntl$month<-strftime(strptime(tempntl$sampledate,'%Y-%m-%d'),'%m')
temp=aggregate(temp$pH,by=list(temp$month),FUN=mean,na.rm=T)
tempntl=aggregate(tempntl$ph,by=list(tempntl$month),FUN=mean)
colnames(temp)<-c('month','pH')
colnames(tempntl)<-c('month','pH')
# windows()
ylim=range(temp$pH,tempntl$pH)
plot(temp$pH~temp$month,pch=16,type='o',ylim=phYlim,cex=cex,lwd=lwd,cex.axis=cex.axis,cex.lab=cex.lab,ylab='pH',
     xlab='')
points(tempntl$pH~tempntl$month,col='grey60',pch=16,type='o',cex=cex,lwd=lwd)

ntlMod_Obs$pH_mod[2]=mean(temp$pH)
ntlMod_Obs$pH_obs[2]=mean(tempntl$pH)
ntlMod_Obs$pH_cor[2]=cor(temp$pH[temp$month%in%tempntl$month],tempntl$pH)

# windows()
# xlim=c(as.POSIXct('1987-01-01 00:00:00'),max(as.POSIXct(CRCout$datetime)))
ylim=c(0,20)
# plot(ntl$dic[ntl$lakeid=='CR'&ntl$depth==0&!is.na(ntl$dic)&ntl$dic>ylim[1]&ntl$dic<ylim[2]&ntl$sampledate<xlim[2]]~
#        ntl$sampledate[ntl$lakeid=='CR'&ntl$depth==0&!is.na(ntl$dic)&ntl$dic>ylim[1]&ntl$dic<ylim[2]&ntl$sampledate<xlim[2]],
#      col='grey60',xlim=xlim,ylim=ylim,pch=16,lwd=3,ylab='',xlab='',yaxt='n',xaxt='n',type='o',cex.lab=1.5)
# par(new=T)
# plot(CRCout$DIC_epi/CRCout$Vepi*12~as.POSIXct(CRCout$datetime,format='%Y-%m-%d'),cex.lab=1.5,xlim=xlim,ylim=ylim,ylab='DIC mg L-1',type='l',lwd=4,xlab='Date',cex.axis=2,xaxt='n')
# axis(1,at=c(as.POSIXct('1990-01-01'),as.POSIXct('2005-01-01'),as.POSIXct('2010-01-01')),labels = c(1990,2005,2010),cex.axis=2)
# legend('topleft',legend = c('Modeled','Observed'),fill = c('black','grey60'),bty='n')

temp=CRCout
tempntl=ntl[ntl$lakeid=='CR'&ntl$depth<=4&!is.na(ntl$dic)&ntl$dic>ylim[1]&ntl$dic<ylim[2]&ntl$sampledate<xlim[2],]
temp$month<-strftime(strptime(temp$datetime,'%Y-%m-%d'),'%m')
tempntl$month<-strftime(strptime(tempntl$sampledate,'%Y-%m-%d'),'%m')
temp=aggregate(temp$DIC_epi/temp$Vepi*12,by=list(temp$month),FUN=mean,na.rm=T)
tempntl=aggregate(tempntl$dic,by=list(tempntl$month),FUN=mean)
colnames(temp)<-c('month','dic')
colnames(tempntl)<-c('month','dic')
# windows()
ylim=range(temp$dic,tempntl$dic)
plot(temp$dic~temp$month,pch=16,type='o',ylim=dicYlim,cex=cex,lwd=lwd,cex.axis=cex.axis,cex.lab=cex.lab,ylab='DIC (mg L-1)',
     xlab='')
points(tempntl$dic~tempntl$month,col='grey60',pch=16,type='o',cex=cex,lwd=lwd)

ntlMod_Obs$DIC_mod[2]=mean(temp$dic)
ntlMod_Obs$DIC_obs[2]=mean(tempntl$dic)
ntlMod_Obs$DIC_cor[2]=cor(temp$dic[temp$month%in%tempntl$month],tempntl$dic)

# windows()
# xlim=c(as.POSIXct('1987-01-01 00:00:00'),max(as.POSIXct(CRCout$datetime)))
ylim=c(0,30)
# plot(ntlChl$chlor[ntlChl$lakeid=='CR'&ntlChl$depth==0&!is.na(ntlChl$chlor)&ntlChl$chlor>ylim[1]&ntlChl$chlor<ylim[2]&ntlChl$sampledate<xlim[2]]~
#        ntlChl$sampledate[ntlChl$lakeid=='CR'&ntlChl$depth==0&!is.na(ntlChl$chlor)&ntlChl$chlor>ylim[1]&ntlChl$chlor<ylim[2]&ntlChl$sampledate<xlim[2]],
#      col='grey60',xlim=xlim,ylim=ylim,pch=16,lwd=3,ylab='',xlab='',yaxt='n',xaxt='n',type='o',cex.lab=1.5)
# par(new=T)
# plot((CRCout$phyto/50)/(CRCout$Vepi)*12*1000~as.POSIXct(CRCout$datetime,format='%Y-%m-%d'),cex.lab=1.5,xlim=xlim,ylim=ylim,ylab='Chl-a ug L-1',type='l',lwd=4,xlab='Date',cex.axis=2,xaxt='n')
# axis(1,at=c(as.POSIXct('1990-01-01'),as.POSIXct('2005-01-01'),as.POSIXct('2010-01-01')),labels = c(1990,2005,2010),cex.axis=2)
# legend('topleft',legend = c('Modeled','Observed'),fill = c('black','grey60'),bty='n')

temp=CRCout
tempntl=ntlChl[ntlChl$lakeid=='CR'&ntlChl$depth<=2&!is.na(ntlChl$chlor)&ntlChl$chlor>ylim[1]&ntlChl$chlor<ylim[2]&ntlChl$sampledate<xlim[2],]
temp$month<-strftime(strptime(temp$datetime,'%Y-%m-%d'),'%m')
tempntl$month<-strftime(strptime(tempntl$sampledate,'%Y-%m-%d'),'%m')
temp=aggregate((temp$phyto/50)/(temp$Vepi)*12*1000,by=list(temp$month),FUN=mean)
tempntl=aggregate(tempntl$chlor,by=list(tempntl$month),FUN=mean)
colnames(temp)<-c('month','chl')
colnames(tempntl)<-c('month','chl')
# windows()
ylim=range(temp$chl,tempntl$chl)
plot(temp$chl~temp$month,pch=16,type='o',ylim=chlYlim,cex=cex,lwd=lwd,cex.axis=cex.axis,cex.lab=cex.lab,ylab='Chl (ug L-1)',
     xlab='')
points(tempntl$chl~tempntl$month,col='grey60',pch=16,type='o',cex=cex,lwd=lwd)

ntlMod_Obs$Chl_mod[2]=mean(temp$chl)
ntlMod_Obs$Chl_obs[2]=mean(tempntl$chl)
ntlMod_Obs$Chl_cor[2]=cor(temp$chl[temp$month%in%tempntl$month],tempntl$chl)

# windows()
# xlim=c(as.POSIXct('1987-01-01 00:00:00'),max(as.POSIXct(CRCout$datetime)))
# plot(ntlTemp$wtemp[ntlTemp$depth==1&ntlTemp$lakeid=='CR']~ntlTemp$sampledate[ntlTemp$depth==1&ntlTemp$lakeid=='CR'],type='l',xlim=xlim,
#      ylab='Temp C',xlab='',lwd=2,cex.lab=1.5,cex.axis=2,col='grey60')
# lines(CRCout$epiTemp~as.POSIXct(CRCout$datetime),lwd=2)

temp=CRCout
tempntl=ntlTemp[ntlTemp$depth==1&ntlTemp$lakeid=='CR'&!is.na(ntlTemp$wtemp),]
temp$month<-strftime(strptime(temp$datetime,'%Y-%m-%d'),'%m')
tempntl$month<-strftime(strptime(tempntl$sampledate,'%Y-%m-%d'),'%m')
temp=aggregate(temp$epiTemp,by=list(temp$month),FUN=mean)
tempntl=aggregate(tempntl$wtemp,by=list(tempntl$month),FUN=mean)
colnames(temp)<-c('month','wtr')
colnames(tempntl)<-c('month','wtr')
# windows()
ylim=range(temp$wtr,tempntl$wtr)
plot(temp$wtr~temp$month,pch=16,type='o',ylim=wtrYlim,cex=cex,lwd=lwd,cex.axis=cex.axis,cex.lab=cex.lab,ylab='Epi Temp (deg C)',
     xlab='')
points(tempntl$wtr~tempntl$month,col='grey60',pch=16,type='o',cex=2)

ntlMod_Obs$Wtr_mod[2]=mean(temp$wtr)
ntlMod_Obs$Wtr_obs[2]=mean(tempntl$wtr)
ntlMod_Obs$Wtr_cor[2]=cor(temp$wtr[temp$month%in%tempntl$month],tempntl$wtr)

# dev.off()
# **********************************************************


# CB ****************************************************
# windows()
# png('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Figures/CB_validation.png',width=49,height = 7,res=300,units='in')
# xlim=c(as.POSIXct('1987-01-01 00:00:00'),max(as.POSIXct(CBCout$datetime)))
ylim=c(2,30)
# par(mfrow=c(1,7),mar=c(5,6,2,2))
cex=2
lwd=2
cex.axis=2
cex.lab=2
# plot(ntl$doc[ntl$lakeid=='CB'&ntl$depth==0&!is.na(ntl$doc)&ntl$doc>ylim[1]&ntl$doc<ylim[2]&ntl$sampledate<xlim[2]]~ntl$sampledate[ntl$lakeid=='CB'&ntl$depth==0&!is.na(ntl$doc)&ntl$doc>ylim[1]&ntl$doc<ylim[2]&ntl$sampledate<xlim[2]],
#      col='grey60',xlim=xlim,ylim=ylim,pch=16,lwd=3,ylab='',xlab='',yaxt='n',xaxt='n',type='o',cex.lab=1.5)
# par(new=T)
# plot((CBCout$DOCr_epi+CBCout$DOCl_epi)/(CBCout$Vepi)*12~as.POSIXct(CBCout$datetime,format='%Y-%m-%d'),cex.lab=1.5,xlim=xlim,ylim=ylim,ylab='DOC mg L-1',type='l',lwd=4,xlab='Date',cex.axis=2,xaxt='n')
# axis(1,at=c(as.POSIXct('1990-01-01'),as.POSIXct('2005-01-01'),as.POSIXct('2010-01-01')),labels = c(1990,2005,2010),cex.axis=2)
# legend('topleft',legend = c('Modeled','Observed'),fill = c('black','grey60'),bty='n')

temp=CBCout
tempntl=ntl[ntl$lakeid=='CB'&ntl$depth<=1&!is.na(ntl$doc)&ntl$doc>ylim[1]&ntl$doc<ylim[2]&ntl$sampledate<xlim[2],]
temp$month<-strftime(strptime(temp$datetime,'%Y-%m-%d'),'%m')
tempntl$month<-strftime(strptime(tempntl$sampledate,'%Y-%m-%d'),'%m')
temp=aggregate((temp$DOCr_epi+temp$DOCl_epi)/temp$Vepi*12,by=list(temp$month),FUN=mean)
tempntl=aggregate(tempntl$doc,by=list(tempntl$month),FUN=mean)
colnames(temp)<-c('month','doc')
colnames(tempntl)<-c('month','doc')
# windows()
ylim=range(temp$doc,tempntl$doc)
plot(temp$doc~temp$month,pch=16,type='o',ylim=docYlim,cex=cex,lwd=lwd,cex.axis=cex.axis,cex.lab=cex.lab,ylab='DOC (mg L-1)',
     xlab='')
points(tempntl$doc~tempntl$month,col='grey60',pch=16,type='o',cex=cex,lwd=lwd)

ntlMod_Obs$lakeID[3]='CB'
ntlMod_Obs$DOC_mod[3]=mean(temp$doc)
ntlMod_Obs$DOC_obs[3]=mean(tempntl$doc)
ntlMod_Obs$DOC_cor[3]=cor(temp$doc[temp$month%in%tempntl$month],tempntl$doc)

# windows()
# xlim=c(as.POSIXct('1987-01-01 00:00:00'),max(as.POSIXct(CBCout$datetime)))
ylim=c(0,50)
# plot(ntl$totpf[ntl$lakeid=='CB'&ntl$depth==0&!is.na(ntl$totpf)&ntl$totpf>ylim[1]&ntl$totpf<ylim[2]&ntl$sampledate<xlim[2]]~
#        ntl$sampledate[ntl$lakeid=='CB'&ntl$depth==0&!is.na(ntl$totpf)&ntl$totpf>ylim[1]&ntl$totpf<ylim[2]&ntl$sampledate<xlim[2]],
#      col='grey60',xlim=xlim,ylim=ylim,pch=16,lwd=3,ylab='',xlab='',yaxt='n',xaxt='n',type='o',cex.lab=1.5)
# par(new=T)
# plot(CBCout$P_epi/(CBCout$Vepi)*31*1000~as.POSIXct(CBCout$datetime,format='%Y-%m-%d'),cex.lab=1.5,xlim=xlim,ylim=ylim,ylab='SRP ug L-1',type='l',lwd=4,xlab='Date',cex.axis=2,xaxt='n')
# axis(1,at=c(as.POSIXct('1990-01-01'),as.POSIXct('2005-01-01'),as.POSIXct('2010-01-01')),labels = c(1990,2005,2010),cex.axis=2)
# legend('topleft',legend = c('Modeled','Observed'),fill = c('black','grey60'),bty='n')

temp=CBCout
tempntl=ntl[ntl$lakeid=='CB'&ntl$depth<=1&!is.na(ntl$totpf)&ntl$totpf>ylim[1]&ntl$totpf<ylim[2]&ntl$sampledate<xlim[2],]
temp$month<-strftime(strptime(temp$datetime,'%Y-%m-%d'),'%m')
tempntl$month<-strftime(strptime(tempntl$sampledate,'%Y-%m-%d'),'%m')
temp=aggregate(temp$P_epi/(temp$Vepi)*31*1000,by=list(temp$month),FUN=mean)
tempntl=aggregate(tempntl$totpf,by=list(tempntl$month),FUN=mean)
colnames(temp)<-c('month','srp')
colnames(tempntl)<-c('month','srp')
# windows()
ylim=range(temp$srp,tempntl$srp)
plot(temp$srp~temp$month,pch=16,type='o',ylim=srpYlim,cex=cex,lwd=lwd,cex.axis=cex.axis,cex.lab=cex.lab,ylab='SRP (ug L-1)',
     xlab='')
points(tempntl$srp~tempntl$month,col='grey60',pch=16,type='o',cex=cex,lwd=lwd)

ntlMod_Obs$SRP_mod[3]=mean(temp$srp)
ntlMod_Obs$SRP_obs[3]=mean(tempntl$srp)
ntlMod_Obs$SRP_cor[3]=cor(temp$srp[temp$month%in%tempntl$month],tempntl$srp)

# windows()
# xlim=c(as.POSIXct('1987-01-01 00:00:00'),max(as.POSIXct(CBCout$datetime)))
ylim=c(0,50)
# plot(ntl$totpuf[ntl$lakeid=='CB'&ntl$depth==0&!is.na(ntl$totpuf)&ntl$totpuf>ylim[1]&ntl$totpuf<ylim[2]&ntl$sampledate<xlim[2]]~
#        ntl$sampledate[ntl$lakeid=='CB'&ntl$depth==0&!is.na(ntl$totpuf)&ntl$totpuf>ylim[1]&ntl$totpuf<ylim[2]&ntl$sampledate<xlim[2]],
#      col='grey60',xlim=xlim,ylim=ylim,pch=16,lwd=3,ylab='',xlab='',yaxt='n',xaxt='n',type='o',cex.lab=1.5)
# par(new=T)
# plot((CBCout$P_epi+CBCout$phyto/95)/(CBCout$Vepi)*31*1000~as.POSIXct(CBCout$datetime,format='%Y-%m-%d'),cex.lab=1.5,xlim=xlim,ylim=ylim,ylab='TP ug L-1',type='l',lwd=4,xlab='Date',cex.axis=2,xaxt='n')
# axis(1,at=c(as.POSIXct('1990-01-01'),as.POSIXct('2005-01-01'),as.POSIXct('2010-01-01')),labels = c(1990,2005,2010),cex.axis=2)
# legend('topleft',legend = c('Modeled','Observed'),fill = c('black','grey60'),bty='n')

temp=CBCout
tempntl=ntl[ntl$lakeid=='CB'&ntl$depth<=1&!is.na(ntl$totpuf)&ntl$totpuf>ylim[1]&ntl$totpuf<ylim[2]&ntl$sampledate<xlim[2],]
temp$month<-strftime(strptime(temp$datetime,'%Y-%m-%d'),'%m')
tempntl$month<-strftime(strptime(tempntl$sampledate,'%Y-%m-%d'),'%m')
temp=aggregate((temp$P_epi+temp$phyto/95)/(temp$Vepi)*31*1000,by=list(temp$month),FUN=mean)
tempntl=aggregate(tempntl$totpuf,by=list(tempntl$month),FUN=mean)
colnames(temp)<-c('month','tp')
colnames(tempntl)<-c('month','tp')
# windows()
ylim=range(temp$tp,tempntl$tp)
plot(temp$tp~temp$month,pch=16,type='o',ylim=tpYlim,cex=cex,lwd=lwd,cex.axis=cex.axis,cex.lab=cex.lab,ylab='TP (ug L-1)',
     xlab='')
points(tempntl$tp~tempntl$month,col='grey60',pch=16,type='o',cex=cex,lwd=lwd)

ntlMod_Obs$TP_mod[3]=mean(temp$tp)
ntlMod_Obs$TP_obs[3]=mean(tempntl$tp)
ntlMod_Obs$TP_cor[3]=cor(temp$tp[temp$month%in%tempntl$month],tempntl$tp)

# windows()
# xlim=c(as.POSIXct('1987-01-01 00:00:00'),max(as.POSIXct(CBCout$datetime)))
ylim=c(3,10)
# plot(ntl$ph[ntl$lakeid=='CB'&ntl$depth==0&!is.na(ntl$ph)&ntl$ph>ylim[1]&ntl$ph<ylim[2]&ntl$sampledate<xlim[2]]~
#        ntl$sampledate[ntl$lakeid=='CB'&ntl$depth==0&!is.na(ntl$ph)&ntl$ph>ylim[1]&ntl$ph<ylim[2]&ntl$sampledate<xlim[2]],
#      col='grey60',xlim=xlim,ylim=ylim,pch=16,lwd=3,ylab='',xlab='',yaxt='n',xaxt='n',type='o',cex.lab=1.5)
# par(new=T)
# plot(CBCout$pH~as.POSIXct(CBCout$datetime,format='%Y-%m-%d'),cex.lab=1.5,xlim=xlim,ylim=ylim,ylab='pH',type='l',lwd=4,xlab='Date',cex.axis=2,xaxt='n')
# axis(1,at=c(as.POSIXct('1990-01-01'),as.POSIXct('2005-01-01'),as.POSIXct('2010-01-01')),labels = c(1990,2005,2010),cex.axis=2)
# legend('topleft',legend = c('Modeled','Observed'),fill = c('black','grey60'),bty='n')

temp=CBCout
tempntl=ntl[ntl$lakeid=='CB'&ntl$depth<=4&!is.na(ntl$ph)&ntl$ph>ylim[1]&ntl$ph<ylim[2]&ntl$sampledate<xlim[2],]
temp$month<-strftime(strptime(temp$datetime,'%Y-%m-%d'),'%m')
tempntl$month<-strftime(strptime(tempntl$sampledate,'%Y-%m-%d'),'%m')
temp=aggregate(temp$pH,by=list(temp$month),FUN=mean,na.rm=T)
tempntl=aggregate(tempntl$ph,by=list(tempntl$month),FUN=mean)
colnames(temp)<-c('month','pH')
colnames(tempntl)<-c('month','pH')
# windows()
ylim=range(temp$pH,tempntl$pH)
plot(temp$pH~temp$month,pch=16,type='o',ylim=phYlim,cex=cex,lwd=lwd,cex.axis=cex.axis,cex.lab=cex.lab,ylab='pH',
     xlab='')
points(tempntl$pH~tempntl$month,col='grey60',pch=16,type='o',cex=cex,lwd=lwd)

ntlMod_Obs$pH_mod[3]=mean(temp$pH)
ntlMod_Obs$pH_obs[3]=mean(tempntl$pH)
ntlMod_Obs$pH_cor[3]=cor(temp$pH[temp$month%in%tempntl$month],tempntl$pH)

# windows()
# xlim=c(as.POSIXct('1987-01-01 00:00:00'),max(as.POSIXct(CBCout$datetime)))
ylim=c(0,20)
# plot(ntl$dic[ntl$lakeid=='CB'&ntl$depth==0&!is.na(ntl$dic)&ntl$dic>ylim[1]&ntl$dic<ylim[2]&ntl$sampledate<xlim[2]]~
#        ntl$sampledate[ntl$lakeid=='CB'&ntl$depth==0&!is.na(ntl$dic)&ntl$dic>ylim[1]&ntl$dic<ylim[2]&ntl$sampledate<xlim[2]],
#      col='grey60',xlim=xlim,ylim=ylim,pch=16,lwd=3,ylab='',xlab='',yaxt='n',xaxt='n',type='o',cex.lab=1.5)
# par(new=T)
# plot(CBCout$DIC_epi/CBCout$Vepi*12~as.POSIXct(CBCout$datetime,format='%Y-%m-%d'),cex.lab=1.5,xlim=xlim,ylim=ylim,ylab='DIC mg L-1',type='l',lwd=4,xlab='Date',cex.axis=2,xaxt='n')
# axis(1,at=c(as.POSIXct('1990-01-01'),as.POSIXct('2005-01-01'),as.POSIXct('2010-01-01')),labels = c(1990,2005,2010),cex.axis=2)
# legend('topleft',legend = c('Modeled','Observed'),fill = c('black','grey60'),bty='n')

temp=CBCout
tempntl=ntl[ntl$lakeid=='CB'&ntl$depth<=4&!is.na(ntl$dic)&ntl$dic>ylim[1]&ntl$dic<ylim[2]&ntl$sampledate<xlim[2],]
temp$month<-strftime(strptime(temp$datetime,'%Y-%m-%d'),'%m')
tempntl$month<-strftime(strptime(tempntl$sampledate,'%Y-%m-%d'),'%m')
temp=aggregate(temp$DIC_epi/temp$Vepi*12,by=list(temp$month),FUN=mean,na.rm=T)
tempntl=aggregate(tempntl$dic,by=list(tempntl$month),FUN=mean)
colnames(temp)<-c('month','dic')
colnames(tempntl)<-c('month','dic')
# windows()
ylim=range(temp$dic,tempntl$dic)
plot(temp$dic~temp$month,pch=16,type='o',ylim=dicYlim,cex=cex,lwd=lwd,cex.axis=cex.axis,cex.lab=cex.lab,ylab='DIC (mg L-1)',
     xlab='')
points(tempntl$dic~tempntl$month,col='grey60',pch=16,type='o',cex=cex,lwd=lwd)

ntlMod_Obs$DIC_mod[3]=mean(temp$dic)
ntlMod_Obs$DIC_obs[3]=mean(tempntl$dic)
ntlMod_Obs$DIC_cor[3]=cor(temp$dic[temp$month%in%tempntl$month],tempntl$dic)

# windows()
# xlim=c(as.POSIXct('1987-01-01 00:00:00'),max(as.POSIXct(CBCout$datetime)))
ylim=c(0,30)
# plot(ntlChl$chlor[ntlChl$lakeid=='CB'&ntlChl$depth==0&!is.na(ntlChl$chlor)&ntlChl$chlor>ylim[1]&ntlChl$chlor<ylim[2]&ntlChl$sampledate<xlim[2]]~
#        ntlChl$sampledate[ntlChl$lakeid=='CB'&ntlChl$depth==0&!is.na(ntlChl$chlor)&ntlChl$chlor>ylim[1]&ntlChl$chlor<ylim[2]&ntlChl$sampledate<xlim[2]],
#      col='grey60',xlim=xlim,ylim=ylim,pch=16,lwd=3,ylab='',xlab='',yaxt='n',xaxt='n',type='o',cex.lab=1.5)
# par(new=T)
# plot((CBCout$phyto/50)/(CBCout$Vepi)*12*1000~as.POSIXct(CBCout$datetime,format='%Y-%m-%d'),cex.lab=1.5,xlim=xlim,ylim=ylim,ylab='Chl-a ug L-1',type='l',lwd=4,xlab='Date',cex.axis=2,xaxt='n')
# axis(1,at=c(as.POSIXct('1990-01-01'),as.POSIXct('2005-01-01'),as.POSIXct('2010-01-01')),labels = c(1990,2005,2010),cex.axis=2)
# legend('topleft',legend = c('Modeled','Observed'),fill = c('black','grey60'),bty='n')

temp=CBCout
tempntl=ntlChl[ntlChl$lakeid=='CB'&ntlChl$depth<=2&!is.na(ntlChl$chlor)&ntlChl$chlor>ylim[1]&ntlChl$chlor<ylim[2]&ntlChl$sampledate<xlim[2],]
temp$month<-strftime(strptime(temp$datetime,'%Y-%m-%d'),'%m')
tempntl$month<-strftime(strptime(tempntl$sampledate,'%Y-%m-%d'),'%m')
temp=aggregate((temp$phyto/50)/(temp$Vepi)*12*1000,by=list(temp$month),FUN=mean)
tempntl=aggregate(tempntl$chlor,by=list(tempntl$month),FUN=mean)
colnames(temp)<-c('month','chl')
colnames(tempntl)<-c('month','chl')
# windows()
ylim=range(temp$chl,tempntl$chl)
plot(temp$chl~temp$month,pch=16,type='o',ylim=chlYlim,cex=cex,lwd=lwd,cex.axis=cex.axis,cex.lab=cex.lab,ylab='Chl (ug L-1)',
     xlab='')
points(tempntl$chl~tempntl$month,col='grey60',pch=16,type='o',cex=cex,lwd=lwd)

ntlMod_Obs$Chl_mod[3]=mean(temp$chl)
ntlMod_Obs$Chl_obs[3]=mean(tempntl$chl)
ntlMod_Obs$Chl_cor[3]=cor(temp$chl[temp$month%in%tempntl$month],tempntl$chl)

# windows()
# xlim=c(as.POSIXct('1987-01-01 00:00:00'),max(as.POSIXct(CBCout$datetime)))
# plot(ntlTemp$wtemp[ntlTemp$depth==1&ntlTemp$lakeid=='CB']~ntlTemp$sampledate[ntlTemp$depth==1&ntlTemp$lakeid=='CB'],type='l',xlim=xlim,
#      ylab='Temp C',xlab='',lwd=2,cex.lab=1.5,cex.axis=2,col='grey60')
# lines(CBCout$epiTemp~as.POSIXct(CBCout$datetime),lwd=2)

temp=CBCout
tempntl=ntlTemp[ntlTemp$depth==1&ntlTemp$lakeid=='CB'&!is.na(ntlTemp$wtemp),]
temp$month<-strftime(strptime(temp$datetime,'%Y-%m-%d'),'%m')
tempntl$month<-strftime(strptime(tempntl$sampledate,'%Y-%m-%d'),'%m')
temp=aggregate(temp$epiTemp,by=list(temp$month),FUN=mean)
tempntl=aggregate(tempntl$wtemp,by=list(tempntl$month),FUN=mean)
colnames(temp)<-c('month','wtr')
colnames(tempntl)<-c('month','wtr')
# windows()
ylim=range(temp$wtr,tempntl$wtr)
plot(temp$wtr~temp$month,pch=16,type='o',ylim=wtrYlim,cex=cex,lwd=lwd,cex.axis=cex.axis,cex.lab=cex.lab,ylab='Epi Temp (deg C)',
     xlab='')
points(tempntl$wtr~tempntl$month,col='grey60',pch=16,type='o',cex=2)

ntlMod_Obs$Wtr_mod[3]=mean(temp$wtr)
ntlMod_Obs$Wtr_obs[3]=mean(tempntl$wtr)
ntlMod_Obs$Wtr_cor[3]=cor(temp$wtr[temp$month%in%tempntl$month],tempntl$wtr)

# dev.off()
# **********************************************************

# BM ****************************************************
# windows()
# png('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Figures/BM_validation.png',width=49,height = 7,res=300,units='in')
# xlim=c(as.POSIXct('1987-01-01 00:00:00'),max(as.POSIXct(BMCout$datetime)))
ylim=c(2,30)
# par(mfrow=c(1,7),mar=c(5,6,2,2))
cex=2
lwd=2
cex.axis=2
cex.lab=2
# plot(ntl$doc[ntl$lakeid=='BM'&ntl$depth==0&!is.na(ntl$doc)&ntl$doc>ylim[1]&ntl$doc<ylim[2]&ntl$sampledate<xlim[2]]~ntl$sampledate[ntl$lakeid=='BM'&ntl$depth==0&!is.na(ntl$doc)&ntl$doc>ylim[1]&ntl$doc<ylim[2]&ntl$sampledate<xlim[2]],
#      col='grey60',xlim=xlim,ylim=ylim,pch=16,lwd=3,ylab='',xlab='',yaxt='n',xaxt='n',type='o',cex.lab=1.5)
# par(new=T)
# plot((BMCout$DOCr_epi+BMCout$DOCl_epi)/(BMCout$Vepi)*12~as.POSIXct(BMCout$datetime,format='%Y-%m-%d'),cex.lab=1.5,xlim=xlim,ylim=ylim,ylab='DOC mg L-1',type='l',lwd=4,xlab='Date',cex.axis=2,xaxt='n')
# axis(1,at=c(as.POSIXct('1990-01-01'),as.POSIXct('2005-01-01'),as.POSIXct('2010-01-01')),labels = c(1990,2005,2010),cex.axis=2)
# legend('topleft',legend = c('Modeled','Observed'),fill = c('black','grey60'),bty='n')

temp=BMCout
tempntl=ntl[ntl$lakeid=='BM'&ntl$depth<=1&!is.na(ntl$doc)&ntl$doc>ylim[1]&ntl$doc<ylim[2]&ntl$sampledate<xlim[2],]
temp$month<-strftime(strptime(temp$datetime,'%Y-%m-%d'),'%m')
tempntl$month<-strftime(strptime(tempntl$sampledate,'%Y-%m-%d'),'%m')
temp=aggregate((temp$DOCr_epi+temp$DOCl_epi)/temp$Vepi*12,by=list(temp$month),FUN=mean)
tempntl=aggregate(tempntl$doc,by=list(tempntl$month),FUN=mean)
colnames(temp)<-c('month','doc')
colnames(tempntl)<-c('month','doc')
# windows()
ylim=range(temp$doc,tempntl$doc)
plot(temp$doc~temp$month,pch=16,type='o',ylim=docYlim,cex=cex,lwd=lwd,cex.axis=cex.axis,cex.lab=cex.lab,ylab='DOC (mg L-1)',
     xlab='')
points(tempntl$doc~tempntl$month,col='grey60',pch=16,type='o',cex=cex,lwd=lwd)

ntlMod_Obs$lakeID[4]='BM'
ntlMod_Obs$DOC_mod[4]=mean(temp$doc)
ntlMod_Obs$DOC_obs[4]=mean(tempntl$doc)
ntlMod_Obs$DOC_cor[4]=cor(temp$doc[temp$month%in%tempntl$month],tempntl$doc)

# windows()
# xlim=c(as.POSIXct('1987-01-01 00:00:00'),max(as.POSIXct(BMCout$datetime)))
ylim=c(0,50)
# plot(ntl$totpf[ntl$lakeid=='BM'&ntl$depth==0&!is.na(ntl$totpf)&ntl$totpf>ylim[1]&ntl$totpf<ylim[2]&ntl$sampledate<xlim[2]]~
#        ntl$sampledate[ntl$lakeid=='BM'&ntl$depth==0&!is.na(ntl$totpf)&ntl$totpf>ylim[1]&ntl$totpf<ylim[2]&ntl$sampledate<xlim[2]],
#      col='grey60',xlim=xlim,ylim=ylim,pch=16,lwd=3,ylab='',xlab='',yaxt='n',xaxt='n',type='o',cex.lab=1.5)
# par(new=T)
# plot(BMCout$P_epi/(BMCout$Vepi)*31*1000~as.POSIXct(BMCout$datetime,format='%Y-%m-%d'),cex.lab=1.5,xlim=xlim,ylim=ylim,ylab='SRP ug L-1',type='l',lwd=4,xlab='Date',cex.axis=2,xaxt='n')
# axis(1,at=c(as.POSIXct('1990-01-01'),as.POSIXct('2005-01-01'),as.POSIXct('2010-01-01')),labels = c(1990,2005,2010),cex.axis=2)
# legend('topleft',legend = c('Modeled','Observed'),fill = c('black','grey60'),bty='n')

temp=BMCout
tempntl=ntl[ntl$lakeid=='BM'&ntl$depth<=1&!is.na(ntl$totpf)&ntl$totpf>ylim[1]&ntl$totpf<ylim[2]&ntl$sampledate<xlim[2],]
temp$month<-strftime(strptime(temp$datetime,'%Y-%m-%d'),'%m')
tempntl$month<-strftime(strptime(tempntl$sampledate,'%Y-%m-%d'),'%m')
temp=aggregate(temp$P_epi/(temp$Vepi)*31*1000,by=list(temp$month),FUN=mean)
tempntl=aggregate(tempntl$totpf,by=list(tempntl$month),FUN=mean)
colnames(temp)<-c('month','srp')
colnames(tempntl)<-c('month','srp')
# windows()
ylim=range(temp$srp,tempntl$srp)
plot(temp$srp~temp$month,pch=16,type='o',ylim=srpYlim,cex=cex,lwd=lwd,cex.axis=cex.axis,cex.lab=cex.lab,ylab='SRP (ug L-1)',
     xlab='')
points(tempntl$srp~tempntl$month,col='grey60',pch=16,type='o',cex=cex,lwd=lwd)

ntlMod_Obs$SRP_mod[4]=mean(temp$srp)
ntlMod_Obs$SRP_obs[4]=mean(tempntl$srp)
ntlMod_Obs$SRP_cor[4]=cor(temp$srp[temp$month%in%tempntl$month],tempntl$srp)

# windows()
# xlim=c(as.POSIXct('1987-01-01 00:00:00'),max(as.POSIXct(BMCout$datetime)))
ylim=c(0,50)
# plot(ntl$totpuf[ntl$lakeid=='BM'&ntl$depth==0&!is.na(ntl$totpuf)&ntl$totpuf>ylim[1]&ntl$totpuf<ylim[2]&ntl$sampledate<xlim[2]]~
#        ntl$sampledate[ntl$lakeid=='BM'&ntl$depth==0&!is.na(ntl$totpuf)&ntl$totpuf>ylim[1]&ntl$totpuf<ylim[2]&ntl$sampledate<xlim[2]],
#      col='grey60',xlim=xlim,ylim=ylim,pch=16,lwd=3,ylab='',xlab='',yaxt='n',xaxt='n',type='o',cex.lab=1.5)
# par(new=T)
# plot((BMCout$P_epi+BMCout$phyto/95)/(BMCout$Vepi)*31*1000~as.POSIXct(BMCout$datetime,format='%Y-%m-%d'),cex.lab=1.5,xlim=xlim,ylim=ylim,ylab='TP ug L-1',type='l',lwd=4,xlab='Date',cex.axis=2,xaxt='n')
# axis(1,at=c(as.POSIXct('1990-01-01'),as.POSIXct('2005-01-01'),as.POSIXct('2010-01-01')),labels = c(1990,2005,2010),cex.axis=2)
# legend('topleft',legend = c('Modeled','Observed'),fill = c('black','grey60'),bty='n')

temp=BMCout
tempntl=ntl[ntl$lakeid=='BM'&ntl$depth<=1&!is.na(ntl$totpuf)&ntl$totpuf>ylim[1]&ntl$totpuf<ylim[2]&ntl$sampledate<xlim[2],]
temp$month<-strftime(strptime(temp$datetime,'%Y-%m-%d'),'%m')
tempntl$month<-strftime(strptime(tempntl$sampledate,'%Y-%m-%d'),'%m')
temp=aggregate((temp$P_epi+temp$phyto/95)/(temp$Vepi)*31*1000,by=list(temp$month),FUN=mean)
tempntl=aggregate(tempntl$totpuf,by=list(tempntl$month),FUN=mean)
colnames(temp)<-c('month','tp')
colnames(tempntl)<-c('month','tp')
# windows()
ylim=range(temp$tp,tempntl$tp)
plot(temp$tp~temp$month,pch=16,type='o',ylim=tpYlim,cex=cex,lwd=lwd,cex.axis=cex.axis,cex.lab=cex.lab,ylab='TP (ug L-1)',
     xlab='')
points(tempntl$tp~tempntl$month,col='grey60',pch=16,type='o',cex=cex,lwd=lwd)

ntlMod_Obs$TP_mod[4]=mean(temp$tp)
ntlMod_Obs$TP_obs[4]=mean(tempntl$tp)
ntlMod_Obs$TP_cor[4]=cor(temp$tp[temp$month%in%tempntl$month],tempntl$tp)

# windows()
# xlim=c(as.POSIXct('1987-01-01 00:00:00'),max(as.POSIXct(BMCout$datetime)))
ylim=c(3,10)
# plot(ntl$ph[ntl$lakeid=='BM'&ntl$depth==0&!is.na(ntl$ph)&ntl$ph>ylim[1]&ntl$ph<ylim[2]&ntl$sampledate<xlim[2]]~
#        ntl$sampledate[ntl$lakeid=='BM'&ntl$depth==0&!is.na(ntl$ph)&ntl$ph>ylim[1]&ntl$ph<ylim[2]&ntl$sampledate<xlim[2]],
#      col='grey60',xlim=xlim,ylim=ylim,pch=16,lwd=3,ylab='',xlab='',yaxt='n',xaxt='n',type='o',cex.lab=1.5)
# par(new=T)
# plot(BMCout$pH~as.POSIXct(BMCout$datetime,format='%Y-%m-%d'),cex.lab=1.5,xlim=xlim,ylim=ylim,ylab='pH',type='l',lwd=4,xlab='Date',cex.axis=2,xaxt='n')
# axis(1,at=c(as.POSIXct('1990-01-01'),as.POSIXct('2005-01-01'),as.POSIXct('2010-01-01')),labels = c(1990,2005,2010),cex.axis=2)
# legend('topleft',legend = c('Modeled','Observed'),fill = c('black','grey60'),bty='n')

temp=BMCout
tempntl=ntl[ntl$lakeid=='BM'&ntl$depth<=4&!is.na(ntl$ph)&ntl$ph>ylim[1]&ntl$ph<ylim[2]&ntl$sampledate<xlim[2],]
temp$month<-strftime(strptime(temp$datetime,'%Y-%m-%d'),'%m')
tempntl$month<-strftime(strptime(tempntl$sampledate,'%Y-%m-%d'),'%m')
temp=aggregate(temp$pH,by=list(temp$month),FUN=mean,na.rm=T)
tempntl=aggregate(tempntl$ph,by=list(tempntl$month),FUN=mean)
colnames(temp)<-c('month','pH')
colnames(tempntl)<-c('month','pH')
# windows()
ylim=range(temp$pH,tempntl$pH)
plot(temp$pH~temp$month,pch=16,type='o',ylim=phYlim,cex=cex,lwd=lwd,cex.axis=cex.axis,cex.lab=cex.lab,ylab='pH',
     xlab='')
points(tempntl$pH~tempntl$month,col='grey60',pch=16,type='o',cex=cex,lwd=lwd)

ntlMod_Obs$pH_mod[4]=mean(temp$pH)
ntlMod_Obs$pH_obs[4]=mean(tempntl$pH)
ntlMod_Obs$pH_cor[4]=cor(temp$pH[temp$month%in%tempntl$month],tempntl$pH)

# windows()
# xlim=c(as.POSIXct('1987-01-01 00:00:00'),max(as.POSIXct(BMCout$datetime)))
ylim=c(0,20)
# plot(ntl$dic[ntl$lakeid=='BM'&ntl$depth==0&!is.na(ntl$dic)&ntl$dic>ylim[1]&ntl$dic<ylim[2]&ntl$sampledate<xlim[2]]~
#        ntl$sampledate[ntl$lakeid=='BM'&ntl$depth==0&!is.na(ntl$dic)&ntl$dic>ylim[1]&ntl$dic<ylim[2]&ntl$sampledate<xlim[2]],
#      col='grey60',xlim=xlim,ylim=ylim,pch=16,lwd=3,ylab='',xlab='',yaxt='n',xaxt='n',type='o',cex.lab=1.5)
# par(new=T)
# plot(BMCout$DIC_epi/BMCout$Vepi*12~as.POSIXct(BMCout$datetime,format='%Y-%m-%d'),cex.lab=1.5,xlim=xlim,ylim=ylim,ylab='DIC mg L-1',type='l',lwd=4,xlab='Date',cex.axis=2,xaxt='n')
# axis(1,at=c(as.POSIXct('1990-01-01'),as.POSIXct('2005-01-01'),as.POSIXct('2010-01-01')),labels = c(1990,2005,2010),cex.axis=2)
# legend('topleft',legend = c('Modeled','Observed'),fill = c('black','grey60'),bty='n')

temp=BMCout
tempntl=ntl[ntl$lakeid=='BM'&ntl$depth<=4&!is.na(ntl$dic)&ntl$dic>ylim[1]&ntl$dic<ylim[2]&ntl$sampledate<xlim[2],]
temp$month<-strftime(strptime(temp$datetime,'%Y-%m-%d'),'%m')
tempntl$month<-strftime(strptime(tempntl$sampledate,'%Y-%m-%d'),'%m')
temp=aggregate(temp$DIC_epi/temp$Vepi*12,by=list(temp$month),FUN=mean,na.rm=T)
tempntl=aggregate(tempntl$dic,by=list(tempntl$month),FUN=mean)
colnames(temp)<-c('month','dic')
colnames(tempntl)<-c('month','dic')
# windows()
ylim=range(temp$dic,tempntl$dic)
plot(temp$dic~temp$month,pch=16,type='o',ylim=dicYlim,cex=cex,lwd=lwd,cex.axis=cex.axis,cex.lab=cex.lab,ylab='DIC (mg L-1)',
     xlab='')
points(tempntl$dic~tempntl$month,col='grey60',pch=16,type='o',cex=cex,lwd=lwd)

ntlMod_Obs$DIC_mod[4]=mean(temp$dic)
ntlMod_Obs$DIC_obs[4]=mean(tempntl$dic)
ntlMod_Obs$DIC_cor[4]=cor(temp$dic[temp$month%in%tempntl$month],tempntl$dic)

# windows()
# xlim=c(as.POSIXct('1987-01-01 00:00:00'),max(as.POSIXct(BMCout$datetime)))
ylim=c(0,30)
# plot(ntlChl$chlor[ntlChl$lakeid=='BM'&ntlChl$depth==0&!is.na(ntlChl$chlor)&ntlChl$chlor>ylim[1]&ntlChl$chlor<ylim[2]&ntlChl$sampledate<xlim[2]]~
#        ntlChl$sampledate[ntlChl$lakeid=='BM'&ntlChl$depth==0&!is.na(ntlChl$chlor)&ntlChl$chlor>ylim[1]&ntlChl$chlor<ylim[2]&ntlChl$sampledate<xlim[2]],
#      col='grey60',xlim=xlim,ylim=ylim,pch=16,lwd=3,ylab='',xlab='',yaxt='n',xaxt='n',type='o',cex.lab=1.5)
# par(new=T)
# plot((BMCout$phyto/50)/(BMCout$Vepi)*12*1000~as.POSIXct(BMCout$datetime,format='%Y-%m-%d'),cex.lab=1.5,xlim=xlim,ylim=ylim,ylab='Chl-a ug L-1',type='l',lwd=4,xlab='Date',cex.axis=2,xaxt='n')
# axis(1,at=c(as.POSIXct('1990-01-01'),as.POSIXct('2005-01-01'),as.POSIXct('2010-01-01')),labels = c(1990,2005,2010),cex.axis=2)
# legend('topleft',legend = c('Modeled','Observed'),fill = c('black','grey60'),bty='n')

temp=BMCout
tempntl=ntlChl[ntlChl$lakeid=='BM'&ntlChl$depth<=2&!is.na(ntlChl$chlor)&ntlChl$chlor>ylim[1]&ntlChl$chlor<ylim[2]&ntlChl$sampledate<xlim[2],]
temp$month<-strftime(strptime(temp$datetime,'%Y-%m-%d'),'%m')
tempntl$month<-strftime(strptime(tempntl$sampledate,'%Y-%m-%d'),'%m')
temp=aggregate((temp$phyto/50)/(temp$Vepi)*12*1000,by=list(temp$month),FUN=mean)
tempntl=aggregate(tempntl$chlor,by=list(tempntl$month),FUN=mean)
colnames(temp)<-c('month','chl')
colnames(tempntl)<-c('month','chl')
# windows()
ylim=range(temp$chl,tempntl$chl)
plot(temp$chl~temp$month,pch=16,type='o',ylim=chlYlim,cex=cex,lwd=lwd,cex.axis=cex.axis,cex.lab=cex.lab,ylab='Chl (ug L-1)',
     xlab='')
points(tempntl$chl~tempntl$month,col='grey60',pch=16,type='o',cex=cex,lwd=lwd)

ntlMod_Obs$Chl_mod[4]=mean(temp$chl)
ntlMod_Obs$Chl_obs[4]=mean(tempntl$chl)
ntlMod_Obs$Chl_cor[4]=cor(temp$chl[temp$month%in%tempntl$month],tempntl$chl)

# windows()
# xlim=c(as.POSIXct('1987-01-01 00:00:00'),max(as.POSIXct(BMCout$datetime)))
# plot(ntlTemp$wtemp[ntlTemp$depth==1&ntlTemp$lakeid=='BM']~ntlTemp$sampledate[ntlTemp$depth==1&ntlTemp$lakeid=='BM'],type='l',xlim=xlim,
#      ylab='Temp C',xlab='',lwd=2,cex.lab=1.5,cex.axis=2,col='grey60')
# lines(BMCout$epiTemp~as.POSIXct(BMCout$datetime),lwd=2)

temp=BMCout
tempntl=ntlTemp[ntlTemp$depth==1&ntlTemp$lakeid=='BM'&!is.na(ntlTemp$wtemp),]
temp$month<-strftime(strptime(temp$datetime,'%Y-%m-%d'),'%m')
tempntl$month<-strftime(strptime(tempntl$sampledate,'%Y-%m-%d'),'%m')
temp=aggregate(temp$epiTemp,by=list(temp$month),FUN=mean)
tempntl=aggregate(tempntl$wtemp,by=list(tempntl$month),FUN=mean)
colnames(temp)<-c('month','wtr')
colnames(tempntl)<-c('month','wtr')
# windows()
ylim=range(temp$wtr,tempntl$wtr)
plot(temp$wtr~temp$month,pch=16,type='o',ylim=wtrYlim,cex=cex,lwd=lwd,cex.axis=cex.axis,cex.lab=cex.lab,ylab='Epi Temp (deg C)',
     xlab='')
points(tempntl$wtr~tempntl$month,col='grey60',pch=16,type='o',cex=2)

ntlMod_Obs$Wtr_mod[4]=mean(temp$wtr)
ntlMod_Obs$Wtr_obs[4]=mean(tempntl$wtr)
ntlMod_Obs$Wtr_cor[4]=cor(temp$wtr[temp$month%in%tempntl$month],tempntl$wtr)

# dev.off()
# **********************************************************


# SP ****************************************************
# windows()
# png('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Figures/SP_validation.png',width=49,height = 7,res=300,units='in')
# xlim=c(as.POSIXct('1987-01-01 00:00:00'),max(as.POSIXct(SPCout$datetime)))
ylim=c(0,10)
# par(mfrow=c(1,7),mar=c(5,6,2,2))
cex=2
lwd=2
cex.axis=2
cex.lab=2
# plot(ntl$doc[ntl$lakeid=='SP'&ntl$depth==0&!is.na(ntl$doc)&ntl$doc>ylim[1]&ntl$doc<ylim[2]&ntl$sampledate<xlim[2]]~ntl$sampledate[ntl$lakeid=='SP'&ntl$depth==0&!is.na(ntl$doc)&ntl$doc>ylim[1]&ntl$doc<ylim[2]&ntl$sampledate<xlim[2]],
#      col='grey60',xlim=xlim,ylim=ylim,pch=16,lwd=3,ylab='',xlab='',yaxt='n',xaxt='n',type='o',cex.lab=1.5)
# par(new=T)
# plot((SPCout$DOCr_epi+SPCout$DOCl_epi)/(SPCout$Vepi)*12~as.POSIXct(SPCout$datetime,format='%Y-%m-%d'),cex.lab=1.5,xlim=xlim,ylim=ylim,ylab='DOC mg L-1',type='l',lwd=4,xlab='Date',cex.axis=2,xaxt='n')
# axis(1,at=c(as.POSIXct('1990-01-01'),as.POSIXct('2005-01-01'),as.POSIXct('2010-01-01')),labels = c(1990,2005,2010),cex.axis=2)
# legend('topleft',legend = c('Modeled','Observed'),fill = c('black','grey60'),bty='n')

temp=SPCout
tempntl=ntl[ntl$lakeid=='SP'&ntl$depth<=1&!is.na(ntl$doc)&ntl$doc>ylim[1]&ntl$doc<ylim[2]&ntl$sampledate<xlim[2],]
temp$month<-strftime(strptime(temp$datetime,'%Y-%m-%d'),'%m')
tempntl$month<-strftime(strptime(tempntl$sampledate,'%Y-%m-%d'),'%m')
temp=aggregate((temp$DOCr_epi+temp$DOCl_epi)/temp$Vepi*12,by=list(temp$month),FUN=mean)
tempntl=aggregate(tempntl$doc,by=list(tempntl$month),FUN=mean)
colnames(temp)<-c('month','doc')
colnames(tempntl)<-c('month','doc')
# windows()
ylim=range(temp$doc,tempntl$doc)
plot(temp$doc~temp$month,pch=16,type='o',ylim=docYlim,cex=cex,lwd=lwd,cex.axis=cex.axis,cex.lab=cex.lab,ylab='DOC (mg L-1)',
     xlab='')
points(tempntl$doc~tempntl$month,col='grey60',pch=16,type='o',cex=cex,lwd=lwd)

ntlMod_Obs$lakeID[5]='SP'
ntlMod_Obs$DOC_mod[5]=mean(temp$doc)
ntlMod_Obs$DOC_obs[5]=mean(tempntl$doc)
ntlMod_Obs$DOC_cor[5]=cor(temp$doc[temp$month%in%tempntl$month],tempntl$doc)

#windows()
# xlim=c(as.POSIXct('1987-01-01 00:00:00'),max(as.POSIXct(SPCout$datetime)))
ylim=c(0,50)
# plot(ntl$totpf[ntl$lakeid=='SP'&ntl$depth==0&!is.na(ntl$totpf)&ntl$totpf>ylim[1]&ntl$totpf<ylim[2]&ntl$sampledate<xlim[2]]~
#        ntl$sampledate[ntl$lakeid=='SP'&ntl$depth==0&!is.na(ntl$totpf)&ntl$totpf>ylim[1]&ntl$totpf<ylim[2]&ntl$sampledate<xlim[2]],
#      col='grey60',xlim=xlim,ylim=ylim,pch=16,lwd=3,ylab='',xlab='',yaxt='n',xaxt='n',type='o',cex.lab=1.5)
# par(new=T)
# plot(SPCout$P_epi/(SPCout$Vepi)*31*1000~as.POSIXct(SPCout$datetime,format='%Y-%m-%d'),cex.lab=1.5,xlim=xlim,ylim=ylim,ylab='SRP ug L-1',type='l',lwd=4,xlab='Date',cex.axis=2,xaxt='n')
# axis(1,at=c(as.POSIXct('1990-01-01'),as.POSIXct('2005-01-01'),as.POSIXct('2010-01-01')),labels = c(1990,2005,2010),cex.axis=2)
# legend('topleft',legend = c('Modeled','Observed'),fill = c('black','grey60'),bty='n')

temp=SPCout
tempntl=ntl[ntl$lakeid=='SP'&ntl$depth<=1&!is.na(ntl$totpf)&ntl$totpf>ylim[1]&ntl$totpf<ylim[2]&ntl$sampledate<xlim[2],]
temp$month<-strftime(strptime(temp$datetime,'%Y-%m-%d'),'%m')
tempntl$month<-strftime(strptime(tempntl$sampledate,'%Y-%m-%d'),'%m')
temp=aggregate(temp$P_epi/(temp$Vepi)*31*1000,by=list(temp$month),FUN=mean)
tempntl=aggregate(tempntl$totpf,by=list(tempntl$month),FUN=mean)
colnames(temp)<-c('month','srp')
colnames(tempntl)<-c('month','srp')
# windows()
ylim=range(temp$srp,tempntl$srp)
plot(temp$srp~temp$month,pch=16,type='o',ylim=srpYlim,cex=cex,lwd=lwd,cex.axis=cex.axis,cex.lab=cex.lab,ylab='SRP (ug L-1)',
     xlab='')
points(tempntl$srp~tempntl$month,col='grey60',pch=16,type='o',cex=cex,lwd=lwd)

ntlMod_Obs$SRP_mod[5]=mean(temp$srp)
ntlMod_Obs$SRP_obs[5]=mean(tempntl$srp)
ntlMod_Obs$SRP_cor[5]=cor(temp$srp[temp$month%in%tempntl$month],tempntl$srp)

# windows()
# xlim=c(as.POSIXct('1987-01-01 00:00:00'),max(as.POSIXct(SPCout$datetime)))
ylim=c(0,50)
# plot(ntl$totpuf[ntl$lakeid=='SP'&ntl$depth==0&!is.na(ntl$totpuf)&ntl$totpuf>ylim[1]&ntl$totpuf<ylim[2]&ntl$sampledate<xlim[2]]~
#        ntl$sampledate[ntl$lakeid=='SP'&ntl$depth==0&!is.na(ntl$totpuf)&ntl$totpuf>ylim[1]&ntl$totpuf<ylim[2]&ntl$sampledate<xlim[2]],
#      col='grey60',xlim=xlim,ylim=ylim,pch=16,lwd=3,ylab='',xlab='',yaxt='n',xaxt='n',type='o',cex.lab=1.5)
# par(new=T)
# plot((SPCout$P_epi+SPCout$phyto/95)/(SPCout$Vepi)*31*1000~as.POSIXct(SPCout$datetime,format='%Y-%m-%d'),cex.lab=1.5,xlim=xlim,ylim=ylim,ylab='TP ug L-1',type='l',lwd=4,xlab='Date',cex.axis=2,xaxt='n')
# axis(1,at=c(as.POSIXct('1990-01-01'),as.POSIXct('2005-01-01'),as.POSIXct('2010-01-01')),labels = c(1990,2005,2010),cex.axis=2)
# legend('topleft',legend = c('Modeled','Observed'),fill = c('black','grey60'),bty='n')

temp=SPCout
tempntl=ntl[ntl$lakeid=='SP'&ntl$depth<=1&!is.na(ntl$totpuf)&ntl$totpuf>ylim[1]&ntl$totpuf<ylim[2]&ntl$sampledate<xlim[2],]
temp$month<-strftime(strptime(temp$datetime,'%Y-%m-%d'),'%m')
tempntl$month<-strftime(strptime(tempntl$sampledate,'%Y-%m-%d'),'%m')
temp=aggregate((temp$P_epi+temp$phyto/95)/(temp$Vepi)*31*1000,by=list(temp$month),FUN=mean)
tempntl=aggregate(tempntl$totpuf,by=list(tempntl$month),FUN=mean)
colnames(temp)<-c('month','tp')
colnames(tempntl)<-c('month','tp')
# windows()
ylim=range(temp$tp,tempntl$tp)
plot(temp$tp~temp$month,pch=16,type='o',ylim=tpYlim,cex=cex,lwd=lwd,cex.axis=cex.axis,cex.lab=cex.lab,ylab='TP (ug L-1)',
     xlab='')
points(tempntl$tp~tempntl$month,col='grey60',pch=16,type='o',cex=cex,lwd=lwd)

ntlMod_Obs$TP_mod[5]=mean(temp$tp)
ntlMod_Obs$TP_obs[5]=mean(tempntl$tp)
ntlMod_Obs$TP_cor[5]=cor(temp$tp[temp$month%in%tempntl$month],tempntl$tp)

# windows()
# xlim=c(as.POSIXct('1987-01-01 00:00:00'),max(as.POSIXct(SPCout$datetime)))
ylim=c(3,10)
# plot(ntl$ph[ntl$lakeid=='SP'&ntl$depth==0&!is.na(ntl$ph)&ntl$ph>ylim[1]&ntl$ph<ylim[2]&ntl$sampledate<xlim[2]]~
#        ntl$sampledate[ntl$lakeid=='SP'&ntl$depth==0&!is.na(ntl$ph)&ntl$ph>ylim[1]&ntl$ph<ylim[2]&ntl$sampledate<xlim[2]],
#      col='grey60',xlim=xlim,ylim=ylim,pch=16,lwd=3,ylab='',xlab='',yaxt='n',xaxt='n',type='o',cex.lab=1.5)
# par(new=T)
# plot(SPCout$pH~as.POSIXct(SPCout$datetime,format='%Y-%m-%d'),cex.lab=1.5,xlim=xlim,ylim=ylim,ylab='pH',type='l',lwd=4,xlab='Date',cex.axis=2,xaxt='n')
# axis(1,at=c(as.POSIXct('1990-01-01'),as.POSIXct('2005-01-01'),as.POSIXct('2010-01-01')),labels = c(1990,2005,2010),cex.axis=2)
# legend('topleft',legend = c('Modeled','Observed'),fill = c('black','grey60'),bty='n')

temp=SPCout
tempntl=ntl[ntl$lakeid=='SP'&ntl$depth<=4&!is.na(ntl$ph)&ntl$ph>ylim[1]&ntl$ph<ylim[2]&ntl$sampledate<xlim[2],]
temp$month<-strftime(strptime(temp$datetime,'%Y-%m-%d'),'%m')
tempntl$month<-strftime(strptime(tempntl$sampledate,'%Y-%m-%d'),'%m')
temp=aggregate(temp$pH,by=list(temp$month),FUN=mean,na.rm=T)
tempntl=aggregate(tempntl$ph,by=list(tempntl$month),FUN=mean)
colnames(temp)<-c('month','pH')
colnames(tempntl)<-c('month','pH')
# windows()
ylim=range(temp$pH,tempntl$pH)
plot(temp$pH~temp$month,pch=16,type='o',ylim=phYlim,cex=cex,lwd=lwd,cex.axis=cex.axis,cex.lab=cex.lab,ylab='pH',
     xlab='')
points(tempntl$pH~tempntl$month,col='grey60',pch=16,type='o',cex=cex,lwd=lwd)

ntlMod_Obs$pH_mod[5]=mean(temp$pH)
ntlMod_Obs$pH_obs[5]=mean(tempntl$pH)
ntlMod_Obs$pH_cor[5]=cor(temp$pH[temp$month%in%tempntl$month],tempntl$pH)

# windows()
# xlim=c(as.POSIXct('1987-01-01 00:00:00'),max(as.POSIXct(SPCout$datetime)))
ylim=c(0,20)
# plot(ntl$dic[ntl$lakeid=='SP'&ntl$depth==0&!is.na(ntl$dic)&ntl$dic>ylim[1]&ntl$dic<ylim[2]&ntl$sampledate<xlim[2]]~
#        ntl$sampledate[ntl$lakeid=='SP'&ntl$depth==0&!is.na(ntl$dic)&ntl$dic>ylim[1]&ntl$dic<ylim[2]&ntl$sampledate<xlim[2]],
#      col='grey60',xlim=xlim,ylim=ylim,pch=16,lwd=3,ylab='',xlab='',yaxt='n',xaxt='n',type='o',cex.lab=1.5)
# par(new=T)
# plot(SPCout$DIC_epi/SPCout$Vepi*12~as.POSIXct(SPCout$datetime,format='%Y-%m-%d'),cex.lab=1.5,xlim=xlim,ylim=ylim,ylab='DIC mg L-1',type='l',lwd=4,xlab='Date',cex.axis=2,xaxt='n')
# axis(1,at=c(as.POSIXct('1990-01-01'),as.POSIXct('2005-01-01'),as.POSIXct('2010-01-01')),labels = c(1990,2005,2010),cex.axis=2)
# legend('topleft',legend = c('Modeled','Observed'),fill = c('black','grey60'),bty='n')

temp=SPCout
tempntl=ntl[ntl$lakeid=='SP'&ntl$depth<=4&!is.na(ntl$dic)&ntl$dic>ylim[1]&ntl$dic<ylim[2]&ntl$sampledate<xlim[2],]
temp$month<-strftime(strptime(temp$datetime,'%Y-%m-%d'),'%m')
tempntl$month<-strftime(strptime(tempntl$sampledate,'%Y-%m-%d'),'%m')
temp=aggregate(temp$DIC_epi/temp$Vepi*12,by=list(temp$month),FUN=mean,na.rm=T)
tempntl=aggregate(tempntl$dic,by=list(tempntl$month),FUN=mean)
colnames(temp)<-c('month','dic')
colnames(tempntl)<-c('month','dic')
# windows()
ylim=range(temp$dic,tempntl$dic)
plot(temp$dic~temp$month,pch=16,type='o',ylim=dicYlim,cex=cex,lwd=lwd,cex.axis=cex.axis,cex.lab=cex.lab,ylab='DIC (mg L-1)',
     xlab='')
points(tempntl$dic~tempntl$month,col='grey60',pch=16,type='o',cex=cex,lwd=lwd)

ntlMod_Obs$DIC_mod[5]=mean(temp$dic)
ntlMod_Obs$DIC_obs[5]=mean(tempntl$dic)
ntlMod_Obs$DIC_cor[5]=cor(temp$dic[temp$month%in%tempntl$month],tempntl$dic)

# windows()
# xlim=c(as.POSIXct('1987-01-01 00:00:00'),max(as.POSIXct(SPCout$datetime)))
ylim=c(0,10)
# plot(ntlChl$chlor[ntlChl$lakeid=='SP'&ntlChl$depth==0&!is.na(ntlChl$chlor)&ntlChl$chlor>ylim[1]&ntlChl$chlor<ylim[2]&ntlChl$sampledate<xlim[2]]~
#        ntlChl$sampledate[ntlChl$lakeid=='SP'&ntlChl$depth==0&!is.na(ntlChl$chlor)&ntlChl$chlor>ylim[1]&ntlChl$chlor<ylim[2]&ntlChl$sampledate<xlim[2]],
#      col='grey60',xlim=xlim,ylim=ylim,pch=16,lwd=3,ylab='',xlab='',yaxt='n',xaxt='n',type='o',cex.lab=1.5)
# par(new=T)
# plot((SPCout$phyto/50)/(SPCout$Vepi)*12*1000~as.POSIXct(SPCout$datetime,format='%Y-%m-%d'),cex.lab=1.5,xlim=xlim,ylim=ylim,ylab='Chl-a ug L-1',type='l',lwd=4,xlab='Date',cex.axis=2,xaxt='n')
# axis(1,at=c(as.POSIXct('1990-01-01'),as.POSIXct('2005-01-01'),as.POSIXct('2010-01-01')),labels = c(1990,2005,2010),cex.axis=2)
# legend('topleft',legend = c('Modeled','Observed'),fill = c('black','grey60'),bty='n')

temp=SPCout
tempntl=ntlChl[ntlChl$lakeid=='SP'&ntlChl$depth<=2&!is.na(ntlChl$chlor)&ntlChl$chlor>ylim[1]&ntlChl$chlor<ylim[2]&ntlChl$sampledate<xlim[2],]
temp$month<-strftime(strptime(temp$datetime,'%Y-%m-%d'),'%m')
tempntl$month<-strftime(strptime(tempntl$sampledate,'%Y-%m-%d'),'%m')
temp=aggregate((temp$phyto/50)/(temp$Vepi)*12*1000,by=list(temp$month),FUN=mean)
tempntl=aggregate(tempntl$chlor,by=list(tempntl$month),FUN=mean)
colnames(temp)<-c('month','chl')
colnames(tempntl)<-c('month','chl')
# windows()
ylim=range(temp$chl,tempntl$chl)
plot(temp$chl~temp$month,pch=16,type='o',ylim=chlYlim,cex=cex,lwd=lwd,cex.axis=cex.axis,cex.lab=cex.lab,ylab='Chl (ug L-1)',
     xlab='')
points(tempntl$chl~tempntl$month,col='grey60',pch=16,type='o',cex=cex,lwd=lwd)

ntlMod_Obs$Chl_mod[5]=mean(temp$chl)
ntlMod_Obs$Chl_obs[5]=mean(tempntl$chl)
ntlMod_Obs$Chl_cor[5]=cor(temp$chl[temp$month%in%tempntl$month],tempntl$chl)

# windows()
# xlim=c(as.POSIXct('1987-01-01 00:00:00'),max(as.POSIXct(SPCout$datetime)))
# plot(ntlTemp$wtemp[ntlTemp$depth==1&ntlTemp$lakeid=='SP']~ntlTemp$sampledate[ntlTemp$depth==1&ntlTemp$lakeid=='SP'],type='l',xlim=xlim,
#      ylab='Temp C',xlab='',lwd=2,cex.lab=1.5,cex.axis=2,col='grey60')
# lines(SPCout$epiTemp~as.POSIXct(SPCout$datetime),lwd=2)

temp=SPCout
tempntl=ntlTemp[ntlTemp$depth==1&ntlTemp$lakeid=='SP'&!is.na(ntlTemp$wtemp),]
temp$month<-strftime(strptime(temp$datetime,'%Y-%m-%d'),'%m')
tempntl$month<-strftime(strptime(tempntl$sampledate,'%Y-%m-%d'),'%m')
temp=aggregate(temp$epiTemp,by=list(temp$month),FUN=mean)
tempntl=aggregate(tempntl$wtemp,by=list(tempntl$month),FUN=mean)
colnames(temp)<-c('month','wtr')
colnames(tempntl)<-c('month','wtr')
# windows()
ylim=range(temp$wtr,tempntl$wtr)
plot(temp$wtr~temp$month,pch=16,type='o',ylim=wtrYlim,cex=cex,lwd=lwd,cex.axis=cex.axis,cex.lab=cex.lab,ylab='Epi Temp (deg C)',
     xlab='')
points(tempntl$wtr~tempntl$month,col='grey60',pch=16,type='o',cex=2)

ntlMod_Obs$Wtr_mod[5]=mean(temp$wtr)
ntlMod_Obs$Wtr_obs[5]=mean(tempntl$wtr)
ntlMod_Obs$Wtr_cor[5]=cor(temp$wtr[temp$month%in%tempntl$month],tempntl$wtr)

# dev.off()
# **********************************************************

# AL ****************************************************
# windows()
# png('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Figures/AL_validation.png',width=49,height = 7,res=300,units='in')
# xlim=c(as.POSIXct('1987-01-01 00:00:00'),max(as.POSIXct(ALCout$datetime)))
ylim=c(2,30)
# par(mfrow=c(1,7),mar=c(5,6,2,2))
cex=2
lwd=2
cex.axis=2
cex.lab=2
# plot(ntl$doc[ntl$lakeid=='AL'&ntl$depth==0&!is.na(ntl$doc)&ntl$doc>ylim[1]&ntl$doc<ylim[2]&ntl$sampledate<xlim[2]]~ntl$sampledate[ntl$lakeid=='AL'&ntl$depth==0&!is.na(ntl$doc)&ntl$doc>ylim[1]&ntl$doc<ylim[2]&ntl$sampledate<xlim[2]],
#      col='grey60',xlim=xlim,ylim=ylim,pch=16,lwd=3,ylab='',xlab='',yaxt='n',xaxt='n',type='o',cex.lab=1.5)
# par(new=T)
# plot((ALCout$DOCr_epi+ALCout$DOCl_epi)/(ALCout$Vepi)*12~as.POSIXct(ALCout$datetime,format='%Y-%m-%d'),cex.lab=1.5,xlim=xlim,ylim=ylim,ylab='DOC mg L-1',type='l',lwd=4,xlab='Date',cex.axis=2,xaxt='n')
# axis(1,at=c(as.POSIXct('1990-01-01'),as.POSIXct('2005-01-01'),as.POSIXct('2010-01-01')),labels = c(1990,2005,2010),cex.axis=2)
# legend('topleft',legend = c('Modeled','Observed'),fill = c('black','grey60'),bty='n')

temp=ALCout
tempntl=ntl[ntl$lakeid=='AL'&ntl$depth<=1&!is.na(ntl$doc)&ntl$doc>ylim[1]&ntl$doc<ylim[2]&ntl$sampledate<xlim[2],]
temp$month<-strftime(strptime(temp$datetime,'%Y-%m-%d'),'%m')
tempntl$month<-strftime(strptime(tempntl$sampledate,'%Y-%m-%d'),'%m')
temp=aggregate((temp$DOCr_epi+temp$DOCl_epi)/temp$Vepi*12,by=list(temp$month),FUN=mean)
tempntl=aggregate(tempntl$doc,by=list(tempntl$month),FUN=mean)
colnames(temp)<-c('month','doc')
colnames(tempntl)<-c('month','doc')
# windows()
ylim=range(temp$doc,tempntl$doc)
plot(temp$doc~temp$month,pch=16,type='o',ylim=docYlim,cex=cex,lwd=lwd,cex.axis=cex.axis,cex.lab=cex.lab,ylab='DOC (mg L-1)',
     xlab='')
points(tempntl$doc~tempntl$month,col='grey60',pch=16,type='o',cex=cex,lwd=lwd)

ntlMod_Obs$lakeID[6]='AL'
ntlMod_Obs$DOC_mod[6]=mean(temp$doc)
ntlMod_Obs$DOC_obs[6]=mean(tempntl$doc)
ntlMod_Obs$DOC_cor[6]=cor(temp$doc[temp$month%in%tempntl$month],tempntl$doc)

# windows()
# xlim=c(as.POSIXct('1987-01-01 00:00:00'),max(as.POSIXct(ALCout$datetime)))
ylim=c(0,50)
# plot(ntl$totpf[ntl$lakeid=='AL'&ntl$depth==0&!is.na(ntl$totpf)&ntl$totpf>ylim[1]&ntl$totpf<ylim[2]&ntl$sampledate<xlim[2]]~
#        ntl$sampledate[ntl$lakeid=='AL'&ntl$depth==0&!is.na(ntl$totpf)&ntl$totpf>ylim[1]&ntl$totpf<ylim[2]&ntl$sampledate<xlim[2]],
#      col='grey60',xlim=xlim,ylim=ylim,pch=16,lwd=3,ylab='',xlab='',yaxt='n',xaxt='n',type='o',cex.lab=1.5)
# par(new=T)
# plot(ALCout$P_epi/(ALCout$Vepi)*31*1000~as.POSIXct(ALCout$datetime,format='%Y-%m-%d'),cex.lab=1.5,xlim=xlim,ylim=ylim,ylab='SRP ug L-1',type='l',lwd=4,xlab='Date',cex.axis=2,xaxt='n')
# axis(1,at=c(as.POSIXct('1990-01-01'),as.POSIXct('2005-01-01'),as.POSIXct('2010-01-01')),labels = c(1990,2005,2010),cex.axis=2)
# legend('topleft',legend = c('Modeled','Observed'),fill = c('black','grey60'),bty='n')

temp=ALCout
tempntl=ntl[ntl$lakeid=='AL'&ntl$depth<=1&!is.na(ntl$totpf)&ntl$totpf>ylim[1]&ntl$totpf<ylim[2]&ntl$sampledate<xlim[2],]
temp$month<-strftime(strptime(temp$datetime,'%Y-%m-%d'),'%m')
tempntl$month<-strftime(strptime(tempntl$sampledate,'%Y-%m-%d'),'%m')
temp=aggregate(temp$P_epi/(temp$Vepi)*31*1000,by=list(temp$month),FUN=mean)
tempntl=aggregate(tempntl$totpf,by=list(tempntl$month),FUN=mean)
colnames(temp)<-c('month','srp')
colnames(tempntl)<-c('month','srp')
# windows()
ylim=range(temp$srp,tempntl$srp)
plot(temp$srp~temp$month,pch=16,type='o',ylim=srpYlim,cex=cex,lwd=lwd,cex.axis=cex.axis,cex.lab=cex.lab,ylab='SRP (ug L-1)',
     xlab='')
points(tempntl$srp~tempntl$month,col='grey60',pch=16,type='o',cex=cex,lwd=lwd)

ntlMod_Obs$SRP_mod[6]=mean(temp$srp)
ntlMod_Obs$SRP_obs[6]=mean(tempntl$srp)
ntlMod_Obs$SRP_cor[6]=cor(temp$srp[temp$month%in%tempntl$month],tempntl$srp)

# windows()
# xlim=c(as.POSIXct('1987-01-01 00:00:00'),max(as.POSIXct(ALCout$datetime)))
ylim=c(0,50)
# plot(ntl$totpuf[ntl$lakeid=='AL'&ntl$depth==0&!is.na(ntl$totpuf)&ntl$totpuf>ylim[1]&ntl$totpuf<ylim[2]&ntl$sampledate<xlim[2]]~
#        ntl$sampledate[ntl$lakeid=='AL'&ntl$depth==0&!is.na(ntl$totpuf)&ntl$totpuf>ylim[1]&ntl$totpuf<ylim[2]&ntl$sampledate<xlim[2]],
#      col='grey60',xlim=xlim,ylim=ylim,pch=16,lwd=3,ylab='',xlab='',yaxt='n',xaxt='n',type='o',cex.lab=1.5)
# par(new=T)
# plot((ALCout$P_epi+ALCout$phyto/95)/(ALCout$Vepi)*31*1000~as.POSIXct(ALCout$datetime,format='%Y-%m-%d'),cex.lab=1.5,xlim=xlim,ylim=ylim,ylab='TP ug L-1',type='l',lwd=4,xlab='Date',cex.axis=2,xaxt='n')
# axis(1,at=c(as.POSIXct('1990-01-01'),as.POSIXct('2005-01-01'),as.POSIXct('2010-01-01')),labels = c(1990,2005,2010),cex.axis=2)
# legend('topleft',legend = c('Modeled','Observed'),fill = c('black','grey60'),bty='n')

temp=ALCout
tempntl=ntl[ntl$lakeid=='AL'&ntl$depth<=1&!is.na(ntl$totpuf)&ntl$totpuf>ylim[1]&ntl$totpuf<ylim[2]&ntl$sampledate<xlim[2],]
temp$month<-strftime(strptime(temp$datetime,'%Y-%m-%d'),'%m')
tempntl$month<-strftime(strptime(tempntl$sampledate,'%Y-%m-%d'),'%m')
temp=aggregate((temp$P_epi+temp$phyto/95)/(temp$Vepi)*31*1000,by=list(temp$month),FUN=mean)
tempntl=aggregate(tempntl$totpuf,by=list(tempntl$month),FUN=mean)
colnames(temp)<-c('month','tp')
colnames(tempntl)<-c('month','tp')
# windows()
ylim=range(temp$tp,tempntl$tp)
plot(temp$tp~temp$month,pch=16,type='o',ylim=tpYlim,cex=cex,lwd=lwd,cex.axis=cex.axis,cex.lab=cex.lab,ylab='TP (ug L-1)',
     xlab='')
points(tempntl$tp~tempntl$month,col='grey60',pch=16,type='o',cex=cex,lwd=lwd)

ntlMod_Obs$TP_mod[6]=mean(temp$tp)
ntlMod_Obs$TP_obs[6]=mean(tempntl$tp)
ntlMod_Obs$TP_cor[6]=cor(temp$tp[temp$month%in%tempntl$month],tempntl$tp)

# windows()
# xlim=c(as.POSIXct('1987-01-01 00:00:00'),max(as.POSIXct(ALCout$datetime)))
ylim=c(3,10)
# plot(ntl$ph[ntl$lakeid=='AL'&ntl$depth==0&!is.na(ntl$ph)&ntl$ph>ylim[1]&ntl$ph<ylim[2]&ntl$sampledate<xlim[2]]~
#        ntl$sampledate[ntl$lakeid=='AL'&ntl$depth==0&!is.na(ntl$ph)&ntl$ph>ylim[1]&ntl$ph<ylim[2]&ntl$sampledate<xlim[2]],
#      col='grey60',xlim=xlim,ylim=ylim,pch=16,lwd=3,ylab='',xlab='',yaxt='n',xaxt='n',type='o',cex.lab=1.5)
# par(new=T)
# plot(ALCout$pH~as.POSIXct(ALCout$datetime,format='%Y-%m-%d'),cex.lab=1.5,xlim=xlim,ylim=ylim,ylab='pH',type='l',lwd=4,xlab='Date',cex.axis=2,xaxt='n')
# axis(1,at=c(as.POSIXct('1990-01-01'),as.POSIXct('2005-01-01'),as.POSIXct('2010-01-01')),labels = c(1990,2005,2010),cex.axis=2)
# legend('topleft',legend = c('Modeled','Observed'),fill = c('black','grey60'),bty='n')

temp=ALCout
tempntl=ntl[ntl$lakeid=='AL'&ntl$depth<=4&!is.na(ntl$ph)&ntl$ph>ylim[1]&ntl$ph<ylim[2]&ntl$sampledate<xlim[2],]
temp$month<-strftime(strptime(temp$datetime,'%Y-%m-%d'),'%m')
tempntl$month<-strftime(strptime(tempntl$sampledate,'%Y-%m-%d'),'%m')
temp=aggregate(temp$pH,by=list(temp$month),FUN=mean,na.rm=T)
tempntl=aggregate(tempntl$ph,by=list(tempntl$month),FUN=mean)
colnames(temp)<-c('month','pH')
colnames(tempntl)<-c('month','pH')
# windows()
ylim=range(temp$pH,tempntl$pH)
plot(temp$pH~temp$month,pch=16,type='o',ylim=phYlim,cex=cex,lwd=lwd,cex.axis=cex.axis,cex.lab=cex.lab,ylab='pH',
     xlab='')
points(tempntl$pH~tempntl$month,col='grey60',pch=16,type='o',cex=cex,lwd=lwd)

ntlMod_Obs$pH_mod[6]=mean(temp$pH)
ntlMod_Obs$pH_obs[6]=mean(tempntl$pH)
ntlMod_Obs$pH_cor[6]=cor(temp$pH[temp$month%in%tempntl$month],tempntl$pH)

# windows()
# xlim=c(as.POSIXct('1987-01-01 00:00:00'),max(as.POSIXct(ALCout$datetime)))
ylim=c(0,20)
# plot(ntl$dic[ntl$lakeid=='AL'&ntl$depth==0&!is.na(ntl$dic)&ntl$dic>ylim[1]&ntl$dic<ylim[2]&ntl$sampledate<xlim[2]]~
#        ntl$sampledate[ntl$lakeid=='AL'&ntl$depth==0&!is.na(ntl$dic)&ntl$dic>ylim[1]&ntl$dic<ylim[2]&ntl$sampledate<xlim[2]],
#      col='grey60',xlim=xlim,ylim=ylim,pch=16,lwd=3,ylab='',xlab='',yaxt='n',xaxt='n',type='o',cex.lab=1.5)
# par(new=T)
# plot(ALCout$DIC_epi/ALCout$Vepi*12~as.POSIXct(ALCout$datetime,format='%Y-%m-%d'),cex.lab=1.5,xlim=xlim,ylim=ylim,ylab='DIC mg L-1',type='l',lwd=4,xlab='Date',cex.axis=2,xaxt='n')
# axis(1,at=c(as.POSIXct('1990-01-01'),as.POSIXct('2005-01-01'),as.POSIXct('2010-01-01')),labels = c(1990,2005,2010),cex.axis=2)
# legend('topleft',legend = c('Modeled','Observed'),fill = c('black','grey60'),bty='n')

temp=ALCout
tempntl=ntl[ntl$lakeid=='AL'&ntl$depth<=4&!is.na(ntl$dic)&ntl$dic>ylim[1]&ntl$dic<ylim[2]&ntl$sampledate<xlim[2],]
temp$month<-strftime(strptime(temp$datetime,'%Y-%m-%d'),'%m')
tempntl$month<-strftime(strptime(tempntl$sampledate,'%Y-%m-%d'),'%m')
temp=aggregate(temp$DIC_epi/temp$Vepi*12,by=list(temp$month),FUN=mean,na.rm=T)
tempntl=aggregate(tempntl$dic,by=list(tempntl$month),FUN=mean)
colnames(temp)<-c('month','dic')
colnames(tempntl)<-c('month','dic')
# windows()
ylim=range(temp$dic,tempntl$dic)
plot(temp$dic~temp$month,pch=16,type='o',ylim=dicYlim,cex=cex,lwd=lwd,cex.axis=cex.axis,cex.lab=cex.lab,ylab='DIC (mg L-1)',
     xlab='')
points(tempntl$dic~tempntl$month,col='grey60',pch=16,type='o',cex=cex,lwd=lwd)

ntlMod_Obs$DIC_mod[6]=mean(temp$dic)
ntlMod_Obs$DIC_obs[6]=mean(tempntl$dic)
ntlMod_Obs$DIC_cor[6]=cor(temp$dic[temp$month%in%tempntl$month],tempntl$dic)

# windows()
# xlim=c(as.POSIXct('1987-01-01 00:00:00'),max(as.POSIXct(ALCout$datetime)))
ylim=c(0,30)
# plot(ntlChl$chlor[ntlChl$lakeid=='AL'&ntlChl$depth==0&!is.na(ntlChl$chlor)&ntlChl$chlor>ylim[1]&ntlChl$chlor<ylim[2]&ntlChl$sampledate<xlim[2]]~
#        ntlChl$sampledate[ntlChl$lakeid=='AL'&ntlChl$depth==0&!is.na(ntlChl$chlor)&ntlChl$chlor>ylim[1]&ntlChl$chlor<ylim[2]&ntlChl$sampledate<xlim[2]],
#      col='grey60',xlim=xlim,ylim=ylim,pch=16,lwd=3,ylab='',xlab='',yaxt='n',xaxt='n',type='o',cex.lab=1.5)
# par(new=T)
# plot((ALCout$phyto/50)/(ALCout$Vepi)*12*1000~as.POSIXct(ALCout$datetime,format='%Y-%m-%d'),cex.lab=1.5,xlim=xlim,ylim=ylim,ylab='Chl-a ug L-1',type='l',lwd=4,xlab='Date',cex.axis=2,xaxt='n')
# axis(1,at=c(as.POSIXct('1990-01-01'),as.POSIXct('2005-01-01'),as.POSIXct('2010-01-01')),labels = c(1990,2005,2010),cex.axis=2)
# legend('topleft',legend = c('Modeled','Observed'),fill = c('black','grey60'),bty='n')

temp=ALCout
tempntl=ntlChl[ntlChl$lakeid=='AL'&ntlChl$depth<=2&!is.na(ntlChl$chlor)&ntlChl$chlor>ylim[1]&ntlChl$chlor<ylim[2]&ntlChl$sampledate<xlim[2],]
temp$month<-strftime(strptime(temp$datetime,'%Y-%m-%d'),'%m')
tempntl$month<-strftime(strptime(tempntl$sampledate,'%Y-%m-%d'),'%m')
temp=aggregate((temp$phyto/50)/(temp$Vepi)*12*1000,by=list(temp$month),FUN=mean)
tempntl=aggregate(tempntl$chlor,by=list(tempntl$month),FUN=mean)
colnames(temp)<-c('month','chl')
colnames(tempntl)<-c('month','chl')
# windows()
ylim=range(temp$chl,tempntl$chl)
plot(temp$chl~temp$month,pch=16,type='o',ylim=chlYlim,cex=cex,lwd=lwd,cex.axis=cex.axis,cex.lab=cex.lab,ylab='Chl (ug L-1)',
     xlab='')
points(tempntl$chl~tempntl$month,col='grey60',pch=16,type='o',cex=cex,lwd=lwd)

ntlMod_Obs$Chl_mod[6]=mean(temp$chl)
ntlMod_Obs$Chl_obs[6]=mean(tempntl$chl)
ntlMod_Obs$Chl_cor[6]=cor(temp$chl[temp$month%in%tempntl$month],tempntl$chl)

# windows()
# xlim=c(as.POSIXct('1987-01-01 00:00:00'),max(as.POSIXct(ALCout$datetime)))
# plot(ntlTemp$wtemp[ntlTemp$depth==1&ntlTemp$lakeid=='AL']~ntlTemp$sampledate[ntlTemp$depth==1&ntlTemp$lakeid=='AL'],type='l',xlim=xlim,
#      ylab='Temp C',xlab='',lwd=2,cex.lab=1.5,cex.axis=2,col='grey60')
# lines(ALCout$epiTemp~as.POSIXct(ALCout$datetime),lwd=2)

temp=ALCout
tempntl=ntlTemp[ntlTemp$depth==1&ntlTemp$lakeid=='AL'&!is.na(ntlTemp$wtemp),]
temp$month<-strftime(strptime(temp$datetime,'%Y-%m-%d'),'%m')
tempntl$month<-strftime(strptime(tempntl$sampledate,'%Y-%m-%d'),'%m')
temp=aggregate(temp$epiTemp,by=list(temp$month),FUN=mean)
tempntl=aggregate(tempntl$wtemp,by=list(tempntl$month),FUN=mean)
colnames(temp)<-c('month','wtr')
colnames(tempntl)<-c('month','wtr')
# windows()
ylim=range(temp$wtr,tempntl$wtr)
plot(temp$wtr~temp$month,pch=16,type='o',ylim=wtrYlim,cex=cex,lwd=lwd,cex.axis=cex.axis,cex.lab=cex.lab,ylab='Epi Temp (deg C)',
     xlab='')
points(tempntl$wtr~tempntl$month,col='grey60',pch=16,type='o',cex=2)

ntlMod_Obs$Wtr_mod[6]=mean(temp$wtr)
ntlMod_Obs$Wtr_obs[6]=mean(tempntl$wtr)
ntlMod_Obs$Wtr_cor[6]=cor(temp$wtr[temp$month%in%tempntl$month],tempntl$wtr)

# dev.off()
# **********************************************************

# TB ****************************************************
# windows()
# png('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Figures/TB_validation.png',width=49,height = 7,res=300,units='in')
# xlim=c(as.POSIXct('1987-01-01 00:00:00'),max(as.POSIXct(TBCout$datetime)))
ylim=c(2,30)
# par(mfrow=c(1,7),mar=c(5,6,2,2))
cex=2
lwd=2
cex.axis=2
cex.lab=2
# plot(ntl$doc[ntl$lakeid=='TB'&ntl$depth==0&!is.na(ntl$doc)&ntl$doc>ylim[1]&ntl$doc<ylim[2]&ntl$sampledate<xlim[2]]~ntl$sampledate[ntl$lakeid=='TB'&ntl$depth==0&!is.na(ntl$doc)&ntl$doc>ylim[1]&ntl$doc<ylim[2]&ntl$sampledate<xlim[2]],
#      col='grey60',xlim=xlim,ylim=ylim,pch=16,lwd=3,ylab='',xlab='',yaxt='n',xaxt='n',type='o',cex.lab=1.5)
# par(new=T)
# plot((TBCout$DOCr_epi+TBCout$DOCl_epi)/(TBCout$Vepi)*12~as.POSIXct(TBCout$datetime,format='%Y-%m-%d'),cex.lab=1.5,xlim=xlim,ylim=ylim,ylab='DOC mg L-1',type='l',lwd=4,xlab='Date',cex.axis=2,xaxt='n')
# axis(1,at=c(as.POSIXct('1990-01-01'),as.POSIXct('2005-01-01'),as.POSIXct('2010-01-01')),labels = c(1990,2005,2010),cex.axis=2)
# legend('topleft',legend = c('Modeled','Observed'),fill = c('black','grey60'),bty='n')

temp=TBCout
tempntl=ntl[ntl$lakeid=='TB'&ntl$depth<=1&!is.na(ntl$doc)&ntl$doc>ylim[1]&ntl$doc<ylim[2]&ntl$sampledate<xlim[2],]
temp$month<-strftime(strptime(temp$datetime,'%Y-%m-%d'),'%m')
tempntl$month<-strftime(strptime(tempntl$sampledate,'%Y-%m-%d'),'%m')
temp=aggregate((temp$DOCr_epi+temp$DOCl_epi)/temp$Vepi*12,by=list(temp$month),FUN=mean)
tempntl=aggregate(tempntl$doc,by=list(tempntl$month),FUN=mean)
colnames(temp)<-c('month','doc')
colnames(tempntl)<-c('month','doc')
# windows()
ylim=range(temp$doc,tempntl$doc)
plot(temp$doc~temp$month,pch=16,type='o',ylim=docYlim,cex=cex,lwd=lwd,cex.axis=cex.axis,cex.lab=cex.lab,ylab='DOC (mg L-1)',
     xlab='')
points(tempntl$doc~tempntl$month,col='grey60',pch=16,type='o',cex=cex,lwd=lwd)

ntlMod_Obs$lakeID[7]='TB'
ntlMod_Obs$DOC_mod[7]=mean(temp$doc)
ntlMod_Obs$DOC_obs[7]=mean(tempntl$doc)
ntlMod_Obs$DOC_cor[7]=cor(temp$doc[temp$month%in%tempntl$month],tempntl$doc)

# windows()
# xlim=c(as.POSIXct('1987-01-01 00:00:00'),max(as.POSIXct(TBCout$datetime)))
ylim=c(0,50)
# plot(ntl$totpf[ntl$lakeid=='TB'&ntl$depth==0&!is.na(ntl$totpf)&ntl$totpf>ylim[1]&ntl$totpf<ylim[2]&ntl$sampledate<xlim[2]]~
#        ntl$sampledate[ntl$lakeid=='TB'&ntl$depth==0&!is.na(ntl$totpf)&ntl$totpf>ylim[1]&ntl$totpf<ylim[2]&ntl$sampledate<xlim[2]],
#      col='grey60',xlim=xlim,ylim=ylim,pch=16,lwd=3,ylab='',xlab='',yaxt='n',xaxt='n',type='o',cex.lab=1.5)
# par(new=T)
# plot(TBCout$P_epi/(TBCout$Vepi)*31*1000~as.POSIXct(TBCout$datetime,format='%Y-%m-%d'),cex.lab=1.5,xlim=xlim,ylim=ylim,ylab='SRP ug L-1',type='l',lwd=4,xlab='Date',cex.axis=2,xaxt='n')
# axis(1,at=c(as.POSIXct('1990-01-01'),as.POSIXct('2005-01-01'),as.POSIXct('2010-01-01')),labels = c(1990,2005,2010),cex.axis=2)
# legend('topleft',legend = c('Modeled','Observed'),fill = c('black','grey60'),bty='n')

temp=TBCout
tempntl=ntl[ntl$lakeid=='TB'&ntl$depth<=1&!is.na(ntl$totpf)&ntl$totpf>ylim[1]&ntl$totpf<ylim[2]&ntl$sampledate<xlim[2],]
temp$month<-strftime(strptime(temp$datetime,'%Y-%m-%d'),'%m')
tempntl$month<-strftime(strptime(tempntl$sampledate,'%Y-%m-%d'),'%m')
temp=aggregate(temp$P_epi/(temp$Vepi)*31*1000,by=list(temp$month),FUN=mean)
tempntl=aggregate(tempntl$totpf,by=list(tempntl$month),FUN=mean)
colnames(temp)<-c('month','srp')
colnames(tempntl)<-c('month','srp')
# windows()
ylim=range(temp$srp,tempntl$srp)
plot(temp$srp~temp$month,pch=16,type='o',ylim=srpYlim,cex=cex,lwd=lwd,cex.axis=cex.axis,cex.lab=cex.lab,ylab='SRP (ug L-1)',
     xlab='')
points(tempntl$srp~tempntl$month,col='grey60',pch=16,type='o',cex=cex,lwd=lwd)

ntlMod_Obs$SRP_mod[7]=mean(temp$srp)
ntlMod_Obs$SRP_obs[7]=mean(tempntl$srp)
ntlMod_Obs$SRP_cor[7]=cor(temp$srp[temp$month%in%tempntl$month],tempntl$srp)

# windows()
# xlim=c(as.POSIXct('1987-01-01 00:00:00'),max(as.POSIXct(TBCout$datetime)))
ylim=c(0,50)
# plot(ntl$totpuf[ntl$lakeid=='TB'&ntl$depth==0&!is.na(ntl$totpuf)&ntl$totpuf>ylim[1]&ntl$totpuf<ylim[2]&ntl$sampledate<xlim[2]]~
#        ntl$sampledate[ntl$lakeid=='TB'&ntl$depth==0&!is.na(ntl$totpuf)&ntl$totpuf>ylim[1]&ntl$totpuf<ylim[2]&ntl$sampledate<xlim[2]],
#      col='grey60',xlim=xlim,ylim=ylim,pch=16,lwd=3,ylab='',xlab='',yaxt='n',xaxt='n',type='o',cex.lab=1.5)
# par(new=T)
# plot((TBCout$P_epi+TBCout$phyto/95)/(TBCout$Vepi)*31*1000~as.POSIXct(TBCout$datetime,format='%Y-%m-%d'),cex.lab=1.5,xlim=xlim,ylim=ylim,ylab='TP ug L-1',type='l',lwd=4,xlab='Date',cex.axis=2,xaxt='n')
# axis(1,at=c(as.POSIXct('1990-01-01'),as.POSIXct('2005-01-01'),as.POSIXct('2010-01-01')),labels = c(1990,2005,2010),cex.axis=2)
# legend('topleft',legend = c('Modeled','Observed'),fill = c('black','grey60'),bty='n')

temp=TBCout
tempntl=ntl[ntl$lakeid=='TB'&ntl$depth<=1&!is.na(ntl$totpuf)&ntl$totpuf>ylim[1]&ntl$totpuf<ylim[2]&ntl$sampledate<xlim[2],]
temp$month<-strftime(strptime(temp$datetime,'%Y-%m-%d'),'%m')
tempntl$month<-strftime(strptime(tempntl$sampledate,'%Y-%m-%d'),'%m')
temp=aggregate((temp$P_epi+temp$phyto/95)/(temp$Vepi)*31*1000,by=list(temp$month),FUN=mean)
tempntl=aggregate(tempntl$totpuf,by=list(tempntl$month),FUN=mean)
colnames(temp)<-c('month','tp')
colnames(tempntl)<-c('month','tp')
# windows()
ylim=range(temp$tp,tempntl$tp)
plot(temp$tp~temp$month,pch=16,type='o',ylim=tpYlim,cex=cex,lwd=lwd,cex.axis=cex.axis,cex.lab=cex.lab,ylab='TP (ug L-1)',
     xlab='')
points(tempntl$tp~tempntl$month,col='grey60',pch=16,type='o',cex=cex,lwd=lwd)

ntlMod_Obs$TP_mod[7]=mean(temp$tp)
ntlMod_Obs$TP_obs[7]=mean(tempntl$tp)
ntlMod_Obs$TP_cor[7]=cor(temp$tp[temp$month%in%tempntl$month],tempntl$tp)

# windows()
# xlim=c(as.POSIXct('1987-01-01 00:00:00'),max(as.POSIXct(TBCout$datetime)))
ylim=c(3,10)
# plot(ntl$ph[ntl$lakeid=='TB'&ntl$depth==0&!is.na(ntl$ph)&ntl$ph>ylim[1]&ntl$ph<ylim[2]&ntl$sampledate<xlim[2]]~
#        ntl$sampledate[ntl$lakeid=='TB'&ntl$depth==0&!is.na(ntl$ph)&ntl$ph>ylim[1]&ntl$ph<ylim[2]&ntl$sampledate<xlim[2]],
#      col='grey60',xlim=xlim,ylim=ylim,pch=16,lwd=3,ylab='',xlab='',yaxt='n',xaxt='n',type='o',cex.lab=1.5)
# par(new=T)
# plot(TBCout$pH~as.POSIXct(TBCout$datetime,format='%Y-%m-%d'),cex.lab=1.5,xlim=xlim,ylim=ylim,ylab='pH',type='l',lwd=4,xlab='Date',cex.axis=2,xaxt='n')
# axis(1,at=c(as.POSIXct('1990-01-01'),as.POSIXct('2005-01-01'),as.POSIXct('2010-01-01')),labels = c(1990,2005,2010),cex.axis=2)
# legend('topleft',legend = c('Modeled','Observed'),fill = c('black','grey60'),bty='n')

temp=TBCout
tempntl=ntl[ntl$lakeid=='TB'&ntl$depth<=4&!is.na(ntl$ph)&ntl$ph>ylim[1]&ntl$ph<ylim[2]&ntl$sampledate<xlim[2],]
temp$month<-strftime(strptime(temp$datetime,'%Y-%m-%d'),'%m')
tempntl$month<-strftime(strptime(tempntl$sampledate,'%Y-%m-%d'),'%m')
temp=aggregate(temp$pH,by=list(temp$month),FUN=mean,na.rm=T)
tempntl=aggregate(tempntl$ph,by=list(tempntl$month),FUN=mean)
colnames(temp)<-c('month','pH')
colnames(tempntl)<-c('month','pH')
# windows()
ylim=range(temp$pH,tempntl$pH)
plot(temp$pH~temp$month,pch=16,type='o',ylim=phYlim,cex=cex,lwd=lwd,cex.axis=cex.axis,cex.lab=cex.lab,ylab='pH',
     xlab='')
points(tempntl$pH~tempntl$month,col='grey60',pch=16,type='o',cex=cex,lwd=lwd)

ntlMod_Obs$pH_mod[7]=mean(temp$pH)
ntlMod_Obs$pH_obs[7]=mean(tempntl$pH)
ntlMod_Obs$pH_cor[7]=cor(temp$pH[temp$month%in%tempntl$month],tempntl$pH)

# windows()
# xlim=c(as.POSIXct('1987-01-01 00:00:00'),max(as.POSIXct(TBCout$datetime)))
ylim=c(0,20)
# plot(ntl$dic[ntl$lakeid=='TB'&ntl$depth==0&!is.na(ntl$dic)&ntl$dic>ylim[1]&ntl$dic<ylim[2]&ntl$sampledate<xlim[2]]~
#        ntl$sampledate[ntl$lakeid=='TB'&ntl$depth==0&!is.na(ntl$dic)&ntl$dic>ylim[1]&ntl$dic<ylim[2]&ntl$sampledate<xlim[2]],
#      col='grey60',xlim=xlim,ylim=ylim,pch=16,lwd=3,ylab='',xlab='',yaxt='n',xaxt='n',type='o',cex.lab=1.5)
# par(new=T)
# plot(TBCout$DIC_epi/TBCout$Vepi*12~as.POSIXct(TBCout$datetime,format='%Y-%m-%d'),cex.lab=1.5,xlim=xlim,ylim=ylim,ylab='DIC mg L-1',type='l',lwd=4,xlab='Date',cex.axis=2,xaxt='n')
# axis(1,at=c(as.POSIXct('1990-01-01'),as.POSIXct('2005-01-01'),as.POSIXct('2010-01-01')),labels = c(1990,2005,2010),cex.axis=2)
# legend('topleft',legend = c('Modeled','Observed'),fill = c('black','grey60'),bty='n')

temp=TBCout
tempntl=ntl[ntl$lakeid=='TB'&ntl$depth<=4&!is.na(ntl$dic)&ntl$dic>ylim[1]&ntl$dic<ylim[2]&ntl$sampledate<xlim[2],]
temp$month<-strftime(strptime(temp$datetime,'%Y-%m-%d'),'%m')
tempntl$month<-strftime(strptime(tempntl$sampledate,'%Y-%m-%d'),'%m')
temp=aggregate(temp$DIC_epi/temp$Vepi*12,by=list(temp$month),FUN=mean,na.rm=T)
tempntl=aggregate(tempntl$dic,by=list(tempntl$month),FUN=mean)
colnames(temp)<-c('month','dic')
colnames(tempntl)<-c('month','dic')
# windows()
ylim=range(temp$dic,tempntl$dic)
plot(temp$dic~temp$month,pch=16,type='o',ylim=dicYlim,cex=cex,lwd=lwd,cex.axis=cex.axis,cex.lab=cex.lab,ylab='DIC (mg L-1)',
     xlab='')
points(tempntl$dic~tempntl$month,col='grey60',pch=16,type='o',cex=cex,lwd=lwd)

ntlMod_Obs$DIC_mod[7]=mean(temp$dic)
ntlMod_Obs$DIC_obs[7]=mean(tempntl$dic)
ntlMod_Obs$DIC_cor[7]=cor(temp$dic[temp$month%in%tempntl$month],tempntl$dic)

# windows()
# xlim=c(as.POSIXct('1987-01-01 00:00:00'),max(as.POSIXct(TBCout$datetime)))
ylim=c(0,30)
# plot(ntlChl$chlor[ntlChl$lakeid=='TB'&ntlChl$depth==0&!is.na(ntlChl$chlor)&ntlChl$chlor>ylim[1]&ntlChl$chlor<ylim[2]&ntlChl$sampledate<xlim[2]]~
#        ntlChl$sampledate[ntlChl$lakeid=='TB'&ntlChl$depth==0&!is.na(ntlChl$chlor)&ntlChl$chlor>ylim[1]&ntlChl$chlor<ylim[2]&ntlChl$sampledate<xlim[2]],
#      col='grey60',xlim=xlim,ylim=ylim,pch=16,lwd=3,ylab='',xlab='',yaxt='n',xaxt='n',type='o',cex.lab=1.5)
# par(new=T)
# plot((TBCout$phyto/50)/(TBCout$Vepi)*12*1000~as.POSIXct(TBCout$datetime,format='%Y-%m-%d'),cex.lab=1.5,xlim=xlim,ylim=ylim,ylab='Chl-a ug L-1',type='l',lwd=4,xlab='Date',cex.axis=2,xaxt='n')
# axis(1,at=c(as.POSIXct('1990-01-01'),as.POSIXct('2005-01-01'),as.POSIXct('2010-01-01')),labels = c(1990,2005,2010),cex.axis=2)
# legend('topleft',legend = c('Modeled','Observed'),fill = c('black','grey60'),bty='n')

temp=TBCout
tempntl=ntlChl[ntlChl$lakeid=='TB'&ntlChl$depth<=2&!is.na(ntlChl$chlor)&ntlChl$chlor>ylim[1]&ntlChl$chlor<ylim[2]&ntlChl$sampledate<xlim[2],]
temp$month<-strftime(strptime(temp$datetime,'%Y-%m-%d'),'%m')
tempntl$month<-strftime(strptime(tempntl$sampledate,'%Y-%m-%d'),'%m')
temp=aggregate((temp$phyto/50)/(temp$Vepi)*12*1000,by=list(temp$month),FUN=mean)
tempntl=aggregate(tempntl$chlor,by=list(tempntl$month),FUN=mean)
colnames(temp)<-c('month','chl')
colnames(tempntl)<-c('month','chl')
# windows()
ylim=range(temp$chl,tempntl$chl)
plot(temp$chl~temp$month,pch=16,type='o',ylim=chlYlim,cex=cex,lwd=lwd,cex.axis=cex.axis,cex.lab=cex.lab,ylab='Chl (ug L-1)',
     xlab='')
points(tempntl$chl~tempntl$month,col='grey60',pch=16,type='o',cex=cex,lwd=lwd)

ntlMod_Obs$Chl_mod[7]=mean(temp$chl)
ntlMod_Obs$Chl_obs[7]=mean(tempntl$chl)
ntlMod_Obs$Chl_cor[7]=cor(temp$chl[temp$month%in%tempntl$month],tempntl$chl)

# windows()
# xlim=c(as.POSIXct('1987-01-01 00:00:00'),max(as.POSIXct(TBCout$datetime)))
# plot(ntlTemp$wtemp[ntlTemp$depth==1&ntlTemp$lakeid=='TB']~ntlTemp$sampledate[ntlTemp$depth==1&ntlTemp$lakeid=='TB'],type='l',xlim=xlim,
#      ylab='Temp C',xlab='',lwd=2,cex.lab=1.5,cex.axis=2,col='grey60')
# lines(TBCout$epiTemp~as.POSIXct(TBCout$datetime),lwd=2)

temp=TBCout
tempntl=ntlTemp[ntlTemp$depth<=1&ntlTemp$lakeid=='TB'&!is.na(ntlTemp$wtemp),]
temp$month<-strftime(strptime(temp$datetime,'%Y-%m-%d'),'%m')
tempntl$month<-strftime(strptime(tempntl$sampledate,'%Y-%m-%d'),'%m')
temp=aggregate(temp$epiTemp,by=list(temp$month),FUN=mean)
tempntl=aggregate(tempntl$wtemp,by=list(tempntl$month),FUN=mean)
colnames(temp)<-c('month','wtr')
colnames(tempntl)<-c('month','wtr')
# windows()
ylim=range(temp$wtr,tempntl$wtr)
plot(temp$wtr~temp$month,pch=16,type='o',ylim=wtrYlim,cex=cex,lwd=lwd,cex.axis=cex.axis,cex.lab=cex.lab,ylab='Epi Temp (deg C)',
     xlab='')
points(tempntl$wtr~tempntl$month,col='grey60',pch=16,type='o',cex=2)

ntlMod_Obs$Wtr_mod[7]=mean(temp$wtr)
ntlMod_Obs$Wtr_obs[7]=mean(tempntl$wtr)
ntlMod_Obs$Wtr_cor[7]=cor(temp$wtr[temp$month%in%tempntl$month],tempntl$wtr)

# windows()
# ylim=range(ntlMod_Obs$DOC_obs,ntlMod_Obs$DOC_mod)
# plot(ntlMod_Obs$DOC_obs~ntlMod_Obs$DOC_mod,pch=16,cex=cex,cex.axis=cex.axis,cex.lab=cex.lab,ylab='Obs DOC (mg L-1)',
#      xlab='Mod DOC (mg L-1)',ylim=ylim,xlim=ylim)
# abline(0,1,lty=2,lwd=lwd)
# ylim=range(ntlMod_Obs$SRP_mod,ntlMod_Obs$SRP_obs)
# plot(ntlMod_Obs$SRP_obs~ntlMod_Obs$SRP_mod,pch=16,cex=cex,cex.axis=cex.axis,cex.lab=cex.lab,ylab='Obs SRP (ug L-1)',
#      xlab='Mod SRP (ug L-1)',ylim=ylim,xlim=ylim)
# abline(0,1,lty=2,lwd=lwd)
# ylim=range(ntlMod_Obs$TP_mod,ntlMod_Obs$TP_obs)
# plot(ntlMod_Obs$TP_obs~ntlMod_Obs$TP_mod,pch=16,cex=cex,cex.axis=cex.axis,cex.lab=cex.lab,ylab='Obs TP (ug L-1)',
#      xlab='Mod TP (ug L-1)',ylim=ylim,xlim=ylim)
# abline(0,1,lty=2,lwd=lwd)
# ylim=range(ntlMod_Obs$pH_obs,ntlMod_Obs$pH_mod)
# plot(ntlMod_Obs$pH_obs~ntlMod_Obs$pH_mod,pch=16,cex=cex,cex.axis=cex.axis,cex.lab=cex.lab,ylab='Obs pH',
#      xlab='Mod pH',ylim=ylim,xlim=ylim)
# abline(0,1,lty=2,lwd=lwd)
# ylim=range(ntlMod_Obs$DIC_mod,ntlMod_Obs$DIC_obs)
# plot(ntlMod_Obs$DIC_obs~ntlMod_Obs$DIC_mod,pch=16,cex=cex,cex.axis=cex.axis,cex.lab=cex.lab,ylab='Obs DIC (mg L-1)',
#      xlab='Mod DIC (mg L-1)',ylim=ylim,xlim=ylim)
# abline(0,1,lty=2,lwd=lwd)
# ylim=range(ntlMod_Obs$Chl_mod,ntlMod_Obs$Chl_obs)
# plot(ntlMod_Obs$Chl_obs~ntlMod_Obs$Chl_mod,pch=16,cex=cex,cex.axis=cex.axis,cex.lab=cex.lab,ylab='Obs Chl (ug L-1)',
#      xlab='Mod Chl (ug L-1)',ylim=ylim,xlim=ylim)
# abline(0,1,lty=2,lwd=lwd)
# ylim=range(ntlMod_Obs$Wtr_mod,ntlMod_Obs$Wtr_obs)
# plot(ntlMod_Obs$Wtr_obs~ntlMod_Obs$Wtr_mod,pch=16,cex=cex,cex.axis=cex.axis,cex.lab=cex.lab,ylab='Obs Wtr (deg C)',
#      xlab='Mod Wtr (deg C)',ylim=ylim,xlim=ylim)
# abline(0,1,lty=2,lwd=lwd)

dev.off()
# **********************************************************
#########################################################################

summary(lm(ntlMod_Obs$DOC_obs~ntlMod_Obs$DOC_mod))
summary(lm(ntlMod_Obs$SRP_obs~ntlMod_Obs$SRP_mod))
summary(lm(ntlMod_Obs$TP_obs~ntlMod_Obs$TP_mod))
summary(lm(ntlMod_Obs$pH_obs~ntlMod_Obs$pH_mod))
summary(lm(ntlMod_Obs$DIC_obs~ntlMod_Obs$DIC_mod))
summary(lm(ntlMod_Obs$Wtr_obs~ntlMod_Obs$Wtr_mod))


# **********************************************************

# comparing to Cole et al. 2007 estimates 
# 1.1e14 Pg C emitted / year; 4.2 million km2 lake surface area = 26.2 g C m-2 year-1 = 95 mg C m-2 day-1; 1.1e14 Pg C emitted / year; 
# 0.9e14 Pg C emitted / year; 4.2 million km2 lake surface area 
coleEmit=0.11e15/(4200000*1000*1000) # average emit g C m-2 yr-1 => could use Raymond's lake area (using Downing's lake area - not clear what they used in Cole et al. 2007)
coleBury=0.05e15/(4200000*1000*1000) # average burial mg C m-2 day-1

resArea=400000*1000*1000 # km2 converted to m2
coleEmitRes=0.28e15/resArea
coleBuryRes=0.18e15/resArea


raymondEmit=.322e15/(3001009*1000*1000) # average emit from Raymond et al 2013 (see supplementary for lake size bins)

# atmCO2=1*390/1e6/29.41*1000 # atm equil in mol C m-3
# raymond0.1<-1*1303/1e6/29.41*1000
# raymond1<-1*1135/1e6/29.41*1000
# raymond10<-1*966/1e6/29.41*1000
# raymond100<-1*818/1e6/29.41*1000

# emit0.1=0.57*(raymond0.1-atmCO2)*211551*1000*1000*365*12 # I don't think this is right, but it's close 
# emit1=0.80*(raymond1-atmCO2)*212233*1000*1000*365*12 # I don't think this is right, but it's close 
# emit10=0.85*(raymond10-atmCO2)*370262*1000*1000*365*12 # I don't think this is right, but it's close 
# emit100=1.09*(raymond100-atmCO2)*402597*1000*1000*365*12 # I don't think this is right, but it's close 
# 
# (emit100+emit10+emit0.1+emit1)/1e15

raymondEmit0.1=0.027/211551/1000/1000*1e15 # g C m-2 yr-1; dividing Raymonds total emissions by lake surface area by class to get an average emissions per lake class 
raymondEmit1=0.032/212233/1000/1000*1e15
raymondEmit10=0.048/370262/1000/1000*1e15
raymondEmit100=0.050/402597/1000/1000*1e15

raymondEmitNHLD=raymondEmit0.1*sum(sum$Area[sum$Area<=1e5])+raymondEmit1*sum(sum$Area[sum$Area>1e5&sum$Area<=1e6])+
  raymondEmit10*sum(sum$Area[sum$Area>1e6&sum$Area<=1e7])+raymondEmit100*sum(sum$Area[sum$Area>1e7]) # g C yr-1 

sum$epi_CO2=sum$fracCO2*sum$DIC_epi/sum$Vepi # mol C m-3 as CO2 
all$epi_CO2=all$fracCO2*all$DIC_epi/all$Vepi # mol C m-3 as CO2 - don't use this for scaling because super high co2 during winter 

# sum$epi_CO2 # 
NHLDresArea=22281214
NHLDlakeArea=777575138

coleNHLDemit=coleEmit*NHLDlakeArea # C emit / day with Cole estimates 
coleNHLDbury=coleBury*NHLDlakeArea # C burial / day with Cole estimates 
coleNHLDemitRes=coleEmitRes*NHLDresArea
coleNHLDburyRes=coleBuryRes*NHLDresArea
# raymondNHLDemit=raymondEmit*(sum(sum$Area)) 

sum$HRT<-sum$Vol/(sum$GWin+sum$SWin+sum$DirectP+sum$Baseflow)
sum$percentEvap<-sum$LakeE/(sum$LakeE+sum$GWout+sum$SWout)
summary(sum$percentEvap)

ndNHLDemit=sum(sum$Emit*12*365) # g C yr-1 
ndNHLDbury=sum(sum$Burial_phyto*12*365,sum$Burial_tPOC*12*365) # g C yr-1 

ndOpenIceEmit=ndNHLDemit
ndOpenIceBurial=ndNHLDbury


ndNHLDemitSmallLakes=sum(sum$Emit[sum$Area<10000]*12*365)
ndNHLDburySmallLakes=sum(sum$Burial_phyto[sum$Area<10000]*12*365,sum$Burial_tPOC[sum$Area<10000]*12*365)

ndNHLDemitLargeLakes=sum(sum$Emit[sum$Area>10000]*12*365)
ndNHLDburyLargeLakes=sum(sum$Burial_phyto[sum$Area>10000]*12*365,sum$Burial_tPOC[sum$Area>10000]*12*365)

sum(sum$Area[sum$Area>220000])/sum(sum$Area) # 22 ha is average lake area for NHLD 

summary(sum$HRT)

ndNHLDemitlowHRT=sum(sum$Emit[sum$HRT<545]*12*365) # 545 days is median 
ndNHLDemitlongHRT=sum(sum$Emit[sum$HRT>545]*12*365)

ndNHLDemitlowPercE=sum(sum$Emit[sum$percentEvap<.705]*12*365) # 0.705 is median 
ndNHLDemithighPercE=sum(sum$Emit[sum$percentEvap>.705]*12*365)

ndNHLDemitlowHRT/ndNHLDemit
ndNHLDemitlongHRT/ndNHLDemit

ndNHLDburyLargeLakes/ndNHLDbury
ndNHLDemitLargeLakes/ndNHLDemit

ndNHLDemitlowPercE/ndNHLDemit
ndNHLDemithighPercE/ndNHLDemit

ndNHLDemit/coleNHLDemit
ndNHLDbury/coleNHLDbury

ndNHLDemit/raymondEmitNHLD

ndNHLDemit/1e9 # Gigagrams of C per year 
ndNHLDbury/1e9
coleNHLDemit/1e9
coleNHLDbury/1e9
raymondEmitNHLD/1e9

ndNHLDemitSmallLakes/sum(sum$Area[sum$Area<10000]) # average areal rate small lakes 
ndNHLDemitLargeLakes/sum(sum$Area[sum$Area>10000]) # average areal rate large lakes 
ndNHLDemit/sum(sum$Area)

ndNHLDburySmallLakes/sum(sum$Area[sum$Area<10000]) # average areal rate small lakes 
ndNHLDburyLargeLakes/sum(sum$Area[sum$Area>10000]) # average areal rate large lakes 
ndNHLDbury/sum(sum$Area)

ndNHLDemitlowPercE/sum(sum$Area[sum$percentEvap<.63])
ndNHLDemithighPercE/sum(sum$Area[sum$percentEvap>.63])

mean(sum$Emit[sum$Area<10000]/sum$Area[sum$Area<10000]*12*365)
mean(sum$Emit[sum$Area>10000]/sum$Area[sum$Area>10000]*12*365)

# average turnover rate of C 
sum$emergent_d20=sum$DOCl_epi/(sum$DOCl_epi+sum$DOCr_epi)*0.08+sum$DOCr_epi/(sum$DOCl_epi+sum$DOCr_epi)*0.002


windows()
plot(sum$emergent_d20~sum$HRT,pch=16,cex=1,cex.lab=1.5,cex.axis=2,ylab='Emergent d20 (day-1)',xlab='HRT (days)')

windows()
plot(log10((sum$Emit/sum$Area*12*1000))~sum$percentEvap,yaxt='n')
axis(2,at = c(log10(10),log10(100),log10(1000)),labels = c(10,100,1000))
plot(log10(sum$Emit/sum$Area*12*1000)~log10(sum$HRT),yaxt='n')
axis(2,at = c(log10(10),log10(100),log10(1000)),labels = c(10,100,1000))

plot(log10(sum$HRT)~log10(sum$Area),yaxt='n')
axis(2,at = c(log10(10),log10(100),log10(1000)),labels = c(10,100,1000))

plot(sum$pH~log10(sum$Area))
points(hanson$pH~log10(hanson$area),col='red',pch=16)


# comparing if using entire year - including ice cover 


coleNHLDemit=coleEmit*(sum(all$Area)) # C emit / day with Cole estimates 
coleNHLDbury=coleBury*(sum(all$Area)) # C burial / day with Cole estimates 
raymondEmitNHLD=raymondEmit0.1*sum(all$Area[all$Area<=1e5])+raymondEmit1*sum(all$Area[all$Area>1e5&all$Area<=1e6])+
  raymondEmit10*sum(all$Area[all$Area>1e6&all$Area<=1e7])+raymondEmit100*sum(all$Area[all$Area>1e7]) # g C yr-1 

summary(all$percentEvap)

ndNHLDemit=sum(all$Emit*12*365)
ndNHLDbury=sum(all$Burial_phyto*12*365,all$Burial_tPOC*12*365)

ndAnnualEmit=ndNHLDemit
ndAnnualBurial=ndNHLDbury

ndNHLDemitSmallLakes=sum(all$Emit[all$Area<10000]*12*365)
ndNHLDburySmallLakes=sum(all$Burial_phyto[all$Area<10000]*12*365,all$Burial_tPOC[all$Area<10000]*12*365)

ndNHLDemitLargeLakes=sum(all$Emit[all$Area>10000]*12*365)
ndNHLDburyLargeLakes=sum(all$Burial_phyto[all$Area>10000]*12*365,all$Burial_tPOC[all$Area>10000]*12*365)

sum(all$Area[all$Area>220000])/sum(all$Area) # 22 ha is average lake area for NHLD 

summary(all$HRT)

ndNHLDemitlowHRT=sum(all$Emit[all$HRT<400]*12*365)
ndNHLDemitlongHRT=sum(all$Emit[all$HRT>400]*12*365)

ndNHLDemitlowPercE=sum(all$Emit[all$percentEvap<.63]*12*365)
ndNHLDemithighPercE=sum(all$Emit[all$percentEvap>.63]*12*365)

ndNHLDemitlowHRT/ndNHLDemit
ndNHLDemitlongHRT/ndNHLDemit

ndNHLDburyLargeLakes/ndNHLDbury
ndNHLDemitLargeLakes/ndNHLDemit

ndNHLDemitlowPercE/ndNHLDemit
ndNHLDemithighPercE/ndNHLDemit

ndNHLDemit/coleNHLDemit
ndNHLDbury/coleNHLDbury

ndNHLDemit/raymondEmitNHLD

ndNHLDemit/1e9 # Gigagrams of C per year 
ndNHLDbury/1e9
coleNHLDemit/1e9
coleNHLDbury/1e9
raymondEmitNHLD/1e9

(coleNHLDemitRes+coleNHLDemit)/1e9
(coleNHLDburyRes+coleNHLDbury)/1e9

# total fluxes Gg C yr-1 
ndAnnualEmit=ndAnnualEmit/1e9
ndAnnualBurial=ndAnnualBurial/1e9
ndOpenIceEmit=ndOpenIceEmit/1e9
ndOpenIceBurial=ndOpenIceBurial/1e9
coleEmit=21
coleBury=9.55
coleEmitRes=35.96
coleBuryRes=19.28
raymondEmit=106.6

# fraction of carbon retained in NHLD 
sum(all$FracRet*all$Area)/sum(all$Area)

allFlux=c(ndAnnualEmit,ndAnnualBurial,ndOpenIceEmit,ndOpenIceBurial,coleEmit,coleBury,coleEmitRes,
          coleBuryRes,raymondEmit)

windows()
png('/Users/Jake/Documents/Jake/MyPapers/Regional Lake Carbon Model - ECI/Figures/fig7_annualEmit_burial.png',
    width = 7,height = 7,units='in',res = 300)
bp.at=c(0,0,1,0,1,0,1,0,1)
par(mar=c(6,6,5,4))
barplot(allFlux,space = bp.at,col=c('gray80','gray20'),border = 'black',lwd=2,cex.axis = 2,cex.lab=2,ylab='Carbon Flux (Gg C yr-1)')
legend('topleft',legend = c('Emissions','Burial'),col = c('gray80','gray20'),bty = 'n',fill = c('gray80','gray20'),cex=2)
dev.off()


