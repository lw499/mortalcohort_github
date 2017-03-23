setwd("~/Dropbox/phd/Missing data comparisons/pek/simwithoutedu_tonly/mortalcohort_github/demo_data")
library(MASS)
library(nlme)
library(gee)
library(geepack)

logit <- function(term) {
  return( ifelse(!is.na(term),log(term/(1-term)),NA) )
}

EXPIT <- function(term) {
  return( ifelse(!is.na(term),exp(term)/(1+exp(term)),NA) )
}

library(foreign)
tmpdata = read.csv("demodat.csv", sep=",")


#CALCULATING THE PERCENT OF MISSING DATA FOR EACH VARIABLE#
pmiss <- function(data) sapply(data,function(x) data.frame(nmiss=sum(is.na(x)), n=length(x), pmiss=sum(is.na(x))/length(x)))

source("augipw.R")
source("augipwd_strat.R")
source("linearincrements.r")

##Code for missing data probability
pmiss <- function(data) sapply(data,function(x) data.frame(nmiss=sum(is.na(x)), n=length(x), pmiss=sum(is.na(x))/length(x)))

##Code for analyzing data using IEE
fit = geeglm(pek ~  time +  I(time^2) + Female + shiftage + Educyrs + Smoke + time*Female + time*shiftage + time*Educyrs + time*Smoke,id=Case,data=tmpdata,corstr="ind")
##Code for analyzing data using LMM
mixed= lme(fixed = pek ~  time +  I(time^2) + Female + shiftage + Educyrs + Smoke + time*Female + time*shiftage + time*Educyrs + time*Smoke, random=~time|Case, data = tmpdata,na.action=na.omit)


###Codes for analyzing data using IPW (under u-MAR)###
######################################################
dat=tmpdata;
dat = dat[dat$s1>=dat$time,]
dat$yl = c(NA, dat$pek[1:nrow(dat)-1])
dat$yl = ifelse(dat$time==0, NA, dat$yl)
dat$lyl = c(NA, dat$yl[1:nrow(dat)-1])
dat$lyl = ifelse(dat$time<=2, NA, dat$lyl)
dat$llyl = c(NA, dat$lyl[1:nrow(dat)-1])
dat$llyl = ifelse(dat$time<=4, NA, dat$llyl)
dat$lllyl = c(NA, dat$llyl[1:nrow(dat)-1])
dat$lllyl = ifelse(dat$time<=6, NA, dat$lllyl)
dat$Ri = ifelse(!is.na(dat$pek),1,0)
dat$Rim1 = c(NA, dat$Ri[1:nrow(dat)-1])   ##Status of the last observation
dat$Rim1 = ifelse(dat$time==0, NA, dat$Rim1)

ipwdat = dat[,c("Case", "time", "shiftage", "Female", "Educyrs", "Ri","pek","s1","Rim1","yl", "lyl","llyl","lllyl","iadl1","HlthPrv1","mmse1","Smoke")]
ipwdat$pred_obs=NA;ipwdat$pi=NA;

ipwdat22 = ipwdat[ipwdat$time==2 & ipwdat$s1>=2 & ipwdat$Rim1==1,]
ipwdat2 = ipwdat22[,c("yl", "shiftage", "Female","Educyrs","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
logistfit = glm(Ri ~ yl +shiftage+Female+Educyrs +iadl1+HlthPrv1+mmse1, family = binomial(link=logit),data = ipwdat22)
ipwdat22$pred_obs = predict(logistfit, newdata = ipwdat2, type="response")
##########

ipwdat44 = ipwdat[ipwdat$time==4 &  ipwdat$s1>=4 & ipwdat$Rim1==1,]
ipwdat4 = ipwdat44[,c("yl", "lyl","llyl","lllyl", "shiftage", "Female","Educyrs","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
logistfit = glm(Ri ~ yl +shiftage+Female+Educyrs +iadl1+HlthPrv1+mmse1, family = binomial(link=logit),data = ipwdat44)
ipwdat44$pred_obs = predict(logistfit, newdata = ipwdat4, type="response")

###############

ipwdat66 = ipwdat[ipwdat$time==6 &  ipwdat$s1>=6 & ipwdat$Rim1==1,]
ipwdat6 = ipwdat66[,c("yl", "lyl","llyl","lllyl", "shiftage", "Female","Educyrs","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
logistfit = glm(Ri ~ yl + shiftage+Female+Educyrs +iadl1+HlthPrv1+mmse1, family = binomial(link=logit),data = ipwdat66)
ipwdat66$pred_obs = predict(logistfit, newdata = ipwdat6, type="response")

#################

ipwdat88 = ipwdat[ipwdat$time==8  &  ipwdat$s1>=8 & ipwdat$Rim1==1,]
ipwdat8 = ipwdat88[,c("yl", "lyl","llyl","lllyl", "shiftage", "Female","Educyrs","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
logistfit = glm(Ri ~ yl + shiftage+Female+Educyrs +iadl1+HlthPrv1+mmse1, family = binomial(link=logit),data = ipwdat88)
ipwdat88$pred_obs = predict(logistfit, newdata = ipwdat8, type="response")
##PREDICT for waves

ipwdat2 = ipwdat22; ipwdat2$pi = ipwdat2$pred_obs; ipwdat2$pred_obs=NULL;
ipwdat4 = rbind(ipwdat22, ipwdat44)
ipwdat4 = ipwdat4[order(ipwdat4$Case),]
uniqueid =  unique(ipwdat4$Case);
ipwdat4$pi = NA; pi=NULL;
for(j in uniqueid)
{
  tmp = ipwdat4[ipwdat4$Case==j,];pit=NULL;
  pit[1]=ifelse(!is.na(tmp$pek[1]),tmp$pred_obs[1],NA)
  if(!is.na(tmp$pred_obs[2]) & !is.na(tmp$time[2])){pit[2]=tmp$pred_obs[2]*tmp$pred_obs[1]}
  pi = c(pi, pit)
}; ipwdat4$pi = pi;
ipwdat6 = rbind(ipwdat22,ipwdat44, ipwdat66)
ipwdat6 = ipwdat6[order(ipwdat6$Case),] 
ipwdat6$pi=NA; pi=NULL;
uniqueid =  unique(ipwdat6$Case);
for(j in uniqueid)
{
  tmp = ipwdat6[ipwdat6$Case==j,];pit=NULL;
  pit[1]=ifelse(!is.na(tmp$pek[1]),tmp$pred_obs[1],NA)
  if(!is.na(tmp$pred_obs[2]) & !is.na(tmp$time[2])){pit[2]=tmp$pred_obs[2]*tmp$pred_obs[1]}
  if(!is.na(tmp$pred_obs[3]) & !is.na(tmp$time[3])){pit[3]=tmp$pred_obs[3]*tmp$pred_obs[2]*tmp$pred_obs[1]}
  pi = c(pi, pit)
}; ipwdat6$pi = pi;
ipwdat8 = rbind(ipwdat22, ipwdat44, ipwdat66, ipwdat88)
ipwdat8 = ipwdat8[order(ipwdat8$Case),]
ipwdat8$pi=NA;pi=NULL;
uniqueid =  unique(ipwdat8$Case);
for(j in uniqueid)
{
  tmp = ipwdat8[ipwdat8$Case==j,];pit=NULL;
  pit[1]=ifelse(!is.na(tmp$pek[1]),tmp$pred_obs[1],NA)
  if(!is.na(tmp$pred_obs[2]) & !is.na(tmp$time[2])){pit[2]=tmp$pred_obs[2]*tmp$pred_obs[1]}
  if(!is.na(tmp$pred_obs[3]) & !is.na(tmp$time[3])){pit[3]=tmp$pred_obs[3]*tmp$pred_obs[2]*tmp$pred_obs[1]}
  if(!is.na(tmp$pred_obs[4]) & !is.na(tmp$time[4])){pit[4]=tmp$pred_obs[4]*tmp$pred_obs[3]*tmp$pred_obs[2]*tmp$pred_obs[1]}
  pi = c(pi, pit)
}; ipwdat8$pi = pi;
ipwdat4 = ipwdat4[!(ipwdat4$time<4),]
ipwdat6 = ipwdat6[!(ipwdat6$time<6),]
ipwdat8 = ipwdat8[!(ipwdat8$time<8),]

ipwdat0 = ipwdat[ipwdat$time==0,]; ipwdat0$pi=1; 
ipwdat$pred_obs=ipwdat0$pred_obs=ipwdat2$pred_obs=ipwdat2$pred=ipwdat4$pred_obs=ipwdat6$pred_obs=ipwdat8$pred_obs=NULL;
ripwdat = rbind(ipwdat0,ipwdat2,ipwdat4,ipwdat6,ipwdat8,ipwdat[ipwdat$s1>=2 & ipwdat$time==2 & ipwdat$Rim1==0,],ipwdat[ipwdat$s1>=4 & ipwdat$time==4 & ipwdat$Rim1==0,],ipwdat[ipwdat$s1>=6 & ipwdat$time==6 & ipwdat$Rim1==0,],ipwdat[ipwdat$s1>=8 & ipwdat$time==8 & ipwdat$Rim1==0,])

ripwdat = ripwdat[order(ripwdat$Case),]; ripwdat$pi = 1/ripwdat$pi;
attach(ripwdat)
ipw.fit = geeglm(pek ~  time +  I(time^2) + Female + shiftage + Educyrs + Smoke + time*Female + time*shiftage + time*Educyrs + time*Smoke, data = ripwdat, weights = pi, id=Case)


###Codes for analyzing data under IPW (under p-MAR)###
######################################################
dat = tmpdata
dat = dat[dat$s1>=dat$time,]
dat$yl = c(NA, dat$pek[1:nrow(dat)-1])
dat$yl = ifelse(dat$time==0, NA, dat$yl)
dat$lyl = c(NA, dat$yl[1:nrow(dat)-1])
dat$lyl = ifelse(dat$time<=2, NA, dat$lyl)
dat$llyl = c(NA, dat$lyl[1:nrow(dat)-1])
dat$llyl = ifelse(dat$time<=4, NA, dat$llyl)
dat$lllyl = c(NA, dat$llyl[1:nrow(dat)-1])
dat$lllyl = ifelse(dat$time<=6, NA, dat$lllyl)
dat$Ri = ifelse(!is.na(dat$pek),1,0)
dat$Rim1 = c(NA, dat$Ri[1:nrow(dat)-1])   ##Status of the last observation
dat$Rim1 = ifelse(dat$time==0, NA, dat$Rim1)

ipwdat = dat[,c("Case", "time", "shiftage", "Female", "Educyrs", "Ri","pek","s1","Rim1","yl", "lyl","llyl","lllyl","iadl1","HlthPrv1","mmse1","Smoke")]
ipwdat$pred_obs=NA;ipwdat$pi=NA;

ipwdat22 = ipwdat[ipwdat$time==2 & ipwdat$s1>=2 & !is.na(ipwdat$Rim1) & ipwdat$Rim1==1,]; 
ipwdat2 = ipwdat22[,c("yl", "shiftage", "Female","Educyrs","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
logistfit = glm(Ri ~ yl +shiftage+Female+Educyrs +iadl1+HlthPrv1+mmse1, family = binomial(link=logit),data = ipwdat22)
ipwdat22$pred_obs = predict(logistfit, newdata = ipwdat2, type="response")

##########

ipwdat44 = ipwdat[ipwdat$s1>=4 & !is.na(ipwdat$Rim1) & ipwdat$Rim1==1,];
ipwdat442 = ipwdat44[ipwdat44$time==2,]; ipwdat444 = ipwdat44[ipwdat44$time==4,]; 
ipwdat42 = ipwdat442[,c("yl", "lyl","llyl", "shiftage", "Female","Educyrs","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
ipwdat44 = ipwdat444[,c("yl", "lyl","llyl", "shiftage", "Female","Educyrs","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
logistfit2 = glm(Ri ~ yl + shiftage+Female+Educyrs+iadl1+HlthPrv1+mmse1, family = binomial(link=logit),data = ipwdat442)
logistfit4 = glm(Ri ~ yl + shiftage+Female+Educyrs+iadl1+HlthPrv1+mmse1, family = binomial(link=logit),data = ipwdat444)
ipwdat442$pred_obs = predict(logistfit2, newdata = ipwdat42, type="response")
ipwdat444$pred_obs = predict(logistfit4, newdata = ipwdat44, type="response")

###############

ipwdat66 = ipwdat[ipwdat$s1>=6 & !is.na(ipwdat$Rim1) & ipwdat$Rim1==1,];
ipwdat662 = ipwdat66[ipwdat66$time==2,]; ipwdat664 = ipwdat66[ipwdat66$time==4,]; ipwdat666 = ipwdat66[ipwdat66$time==6,];
ipwdat62 = ipwdat662[,c("yl", "lyl","llyl","shiftage", "Female","Educyrs","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
ipwdat64 = ipwdat664[,c("yl", "lyl","llyl","shiftage", "Female","Educyrs","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
ipwdat66 = ipwdat666[,c("yl", "lyl","llyl","shiftage", "Female","Educyrs","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
logistfit2 = glm(Ri ~ yl+shiftage+Female+Educyrs+iadl1+HlthPrv1+mmse1, family = binomial(link=logit),data = ipwdat662)
logistfit4 = glm(Ri ~ yl+shiftage+Female+Educyrs+iadl1+HlthPrv1+mmse1, family = binomial(link=logit),data = ipwdat664)
logistfit6 = glm(Ri ~ yl+shiftage+Female+Educyrs+iadl1+HlthPrv1+mmse1, family = binomial(link=logit),data = ipwdat666)
ipwdat662$pred_obs = predict(logistfit2, newdata = ipwdat62, type="response")
ipwdat664$pred_obs = predict(logistfit4, newdata = ipwdat64, type="response")
ipwdat666$pred_obs = predict(logistfit6, newdata = ipwdat66, type="response")

#################

ipwdat88 = ipwdat[ipwdat$s1==8 & !is.na(ipwdat$Rim1) & ipwdat$Rim1==1,]
ipwdat882 = ipwdat88[ipwdat88$time==2,]; ipwdat884 = ipwdat88[ipwdat88$time==4,]; ipwdat886 = ipwdat88[ipwdat88$time==6,]; ipwdat888 = ipwdat88[ipwdat88$time==8,];
ipwdat82 = ipwdat882[,c("yl", "lyl","llyl","lllyl", "shiftage", "Female","Educyrs","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
ipwdat84 = ipwdat884[,c("yl", "lyl","llyl","lllyl", "shiftage", "Female","Educyrs","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
ipwdat86 = ipwdat886[,c("yl", "lyl","llyl","lllyl", "shiftage", "Female","Educyrs","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
ipwdat88 = ipwdat888[,c("yl", "lyl","llyl","lllyl", "shiftage", "Female","Educyrs","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
logistfit2 = glm(Ri ~ yl + shiftage+Female+Educyrs+iadl1+HlthPrv1+mmse1, family = binomial(link=logit),data = ipwdat882)
logistfit4 = glm(Ri ~ yl + shiftage+Female+Educyrs+iadl1+HlthPrv1+mmse1, family = binomial(link=logit),data = ipwdat884)
logistfit6 = glm(Ri ~ yl + shiftage+Female+Educyrs+iadl1+HlthPrv1+mmse1, family = binomial(link=logit),data = ipwdat886)
logistfit8 = glm(Ri ~ yl + shiftage+Female+Educyrs+iadl1+HlthPrv1+mmse1, family = binomial(link=logit),data = ipwdat888)
ipwdat882$pred_obs = predict(logistfit2, newdata = ipwdat82, type="response")
ipwdat884$pred_obs = predict(logistfit4, newdata = ipwdat84, type="response")
ipwdat886$pred_obs = predict(logistfit6, newdata = ipwdat86, type="response")
ipwdat888$pred_obs = predict(logistfit8, newdata = ipwdat88, type="response")

ipwdat2 = ipwdat22; ipwdat2$pi = ipwdat2$pred_obs; ipwdat2$pred_obs=NULL;
ipwdat4 = rbind(ipwdat442, ipwdat444)
ipwdat4 = ipwdat4[order(ipwdat4$Case),]
uniqueid =  unique(ipwdat4$Case);
ipwdat4$pi = NA; pi=NULL;
for(j in uniqueid)
{
  tmp = ipwdat4[ipwdat4$Case==j,];pit=NULL;
  pit[1]=ifelse(!is.na(tmp$pek[1]),tmp$pred_obs[1],NA)
  if(!is.na(tmp$pred_obs[2]) & !is.na(tmp$time[2])){pit[2]=tmp$pred_obs[2]*tmp$pred_obs[1]}
  pi = c(pi, pit)
}; ipwdat4$pi = pi;
ipwdat6 = rbind(ipwdat662,ipwdat664, ipwdat666)
ipwdat6 = ipwdat6[order(ipwdat6$Case),] 
ipwdat6$pi=NA; pi=NULL;
uniqueid =  unique(ipwdat6$Case);
for(j in uniqueid)
{
  tmp = ipwdat6[ipwdat6$Case==j,];pit=NULL;
  pit[1]=ifelse(!is.na(tmp$pek[1]),tmp$pred_obs[1],NA)
  if(!is.na(tmp$pred_obs[2]) & !is.na(tmp$time[2])){pit[2]=tmp$pred_obs[2]*tmp$pred_obs[1]}
  if(!is.na(tmp$pred_obs[3]) & !is.na(tmp$time[3])){pit[3]=tmp$pred_obs[3]*tmp$pred_obs[2]*tmp$pred_obs[1]}
  pi = c(pi, pit)
}; ipwdat6$pi = pi;
ipwdat8 = rbind(ipwdat882, ipwdat884, ipwdat886, ipwdat888)
ipwdat8 = ipwdat8[order(ipwdat8$Case),]
ipwdat8$pi=NA;pi=NULL;
uniqueid =  unique(ipwdat8$Case);
for(j in uniqueid)
{
  tmp = ipwdat8[ipwdat8$Case==j,];pit=NULL;
  pit[1]=ifelse(!is.na(tmp$pek[1]),tmp$pred_obs[1],NA)
  if(!is.na(tmp$pred_obs[2]) & !is.na(tmp$time[2])){pit[2]=tmp$pred_obs[2]*tmp$pred_obs[1]}
  if(!is.na(tmp$pred_obs[3]) & !is.na(tmp$time[3])){pit[3]=tmp$pred_obs[3]*tmp$pred_obs[2]*tmp$pred_obs[1]}
  if(!is.na(tmp$pred_obs[4]) & !is.na(tmp$time[4])){pit[4]=tmp$pred_obs[4]*tmp$pred_obs[3]*tmp$pred_obs[2]*tmp$pred_obs[1]}
  pi = c(pi, pit)
}; ipwdat8$pi = pi;
ipwdat4 = ipwdat4[!(ipwdat4$time<4),]
ipwdat6 = ipwdat6[!(ipwdat6$time<6),]
ipwdat8 = ipwdat8[!(ipwdat8$time<8),]

ipwdat0 = ipwdat[ipwdat$time==0,]; ipwdat0$pi=1; 
ipwdat$pred_obs=ipwdat0$pred_obs=ipwdat2$pred_obs=ipwdat2$pred=ipwdat4$pred_obs=ipwdat6$pred_obs=ipwdat8$pred_obs=NULL;
ripwdat = rbind(ipwdat0,ipwdat2,ipwdat4,ipwdat6, ipwdat8,ipwdat[ipwdat$time>0 & ipwdat$Rim1==0,])

ripwdat = ripwdat[order(ripwdat$Case),];ripwdat$pi=1/ripwdat$pi;
ipw.fitp = geeglm(pek ~  time +  I(time^2) + Female + shiftage + Educyrs + Smoke + time*Female + time*shiftage + time*Educyrs + time*Smoke, data = ripwdat, weights = pi, id=Case)


###Codes for analyzing data under IPW (f-MAR)###
################################################
dat=tmpdata;
dat = dat[dat$s1>=dat$time,]
dat$yl = c(NA, dat$pek[1:nrow(dat)-1])
dat$yl = ifelse(dat$time==0, NA, dat$yl)
dat$lyl = c(NA, dat$yl[1:nrow(dat)-1])
dat$lyl = ifelse(dat$time<=2, NA, dat$lyl)
dat$llyl = c(NA, dat$lyl[1:nrow(dat)-1])
dat$llyl = ifelse(dat$time<=4, NA, dat$llyl)
dat$lllyl = c(NA, dat$llyl[1:nrow(dat)-1])
dat$lllyl = ifelse(dat$time<=6, NA, dat$lllyl)
dat$Ri = ifelse(!is.na(dat$pek),1,0)
dat$Rim1 = c(NA, dat$Ri[1:nrow(dat)-1])   ##Status of the last observation
dat$Rim1 = ifelse(dat$time==0, NA, dat$Rim1)

###IPW-IEE-D
ipwdat = dat[,c("Case", "time", "shiftage", "Female", "Educyrs", "Ri","pek","s1","Rim1","yl", "lyl","llyl","lllyl","iadl1","HlthPrv1","mmse1","Smoke")]
ipwdat$pred_obs=NA;ipwdat$Rim1=ifelse(ipwdat$time==0, 0, ipwdat$Rim1);

ipwdat22 = ipwdat[ipwdat$s1==2 & ipwdat$Rim1==1,]; 
ipwdat2 = ipwdat22[,c("yl", "shiftage", "Female","Educyrs","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
logistfit = glm(Ri ~ yl +shiftage+Female+Educyrs+iadl1+HlthPrv1+mmse1, family = binomial(link=logit),data = ipwdat22)
ipwdat22$pred_obs = predict(logistfit, newdata = ipwdat2, type="response")

##########

ipwdat44 = ipwdat[ipwdat$s1==4 & ipwdat$Rim1==1,];
ipwdat442 = ipwdat44[ipwdat44$time==2,]; ipwdat444 = ipwdat44[ipwdat44$time==4,]; 
ipwdat42 = ipwdat442[,c("yl", "lyl","llyl","lllyl", "shiftage", "Female","Educyrs","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
ipwdat44 = ipwdat444[,c("yl", "lyl","llyl","lllyl", "shiftage", "Female","Educyrs","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
logistfit2 = glm(Ri ~ yl + shiftage+Female+Educyrs+iadl1+HlthPrv1+mmse1, family = binomial(link=logit),data = ipwdat442)
logistfit4 = glm(Ri ~ yl + shiftage+Female+Educyrs+iadl1+HlthPrv1+mmse1, family = binomial(link=logit),data = ipwdat444)
ipwdat442$pred_obs = predict(logistfit2, newdata = ipwdat42, type="response")
ipwdat444$pred_obs = predict(logistfit4, newdata = ipwdat44, type="response")

###############

ipwdat66 = ipwdat[ipwdat$s1==6 & ipwdat$Rim1==1,];
ipwdat662 = ipwdat66[ipwdat66$time==2,]; ipwdat664 = ipwdat66[ipwdat66$time==4,]; ipwdat666 = ipwdat66[ipwdat66$time==6,];
ipwdat62 = ipwdat662[,c("yl", "lyl","llyl","lllyl", "shiftage", "Female","Educyrs","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
ipwdat64 = ipwdat664[,c("yl", "lyl","llyl","lllyl", "shiftage", "Female","Educyrs","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
ipwdat66 = ipwdat666[,c("yl", "lyl","llyl","lllyl", "shiftage", "Female","Educyrs","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
logistfit2 = glm(Ri ~ yl+shiftage+Female+Educyrs+iadl1+HlthPrv1+mmse1, family = binomial(link=logit),data = ipwdat662)
logistfit4 = glm(Ri ~ yl+shiftage+Female+Educyrs+iadl1+HlthPrv1+mmse1, family = binomial(link=logit),data = ipwdat664)
logistfit6 = glm(Ri ~ yl+shiftage+Female+Educyrs+iadl1+HlthPrv1+mmse1, family = binomial(link=logit),data = ipwdat666)
ipwdat662$pred_obs = predict(logistfit2, newdata = ipwdat62, type="response")
ipwdat664$pred_obs = predict(logistfit4, newdata = ipwdat64, type="response")
ipwdat666$pred_obs = predict(logistfit6, newdata = ipwdat66, type="response")

#################

ipwdat88 = ipwdat[ipwdat$s1==8 & ipwdat$Rim1==1,]
ipwdat882 = ipwdat88[ipwdat88$time==2,]; ipwdat884 = ipwdat88[ipwdat88$time==4,]; ipwdat886 = ipwdat88[ipwdat88$time==6,]; ipwdat888 = ipwdat88[ipwdat88$time==8,];
ipwdat82 = ipwdat882[,c("yl", "lyl","llyl","lllyl", "shiftage", "Female","Educyrs","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
ipwdat84 = ipwdat884[,c("yl", "lyl","llyl","lllyl", "shiftage", "Female","Educyrs","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
ipwdat86 = ipwdat886[,c("yl", "lyl","llyl","lllyl", "shiftage", "Female","Educyrs","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
ipwdat88 = ipwdat888[,c("yl", "lyl","llyl","lllyl", "shiftage", "Female","Educyrs","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
logistfit2 = glm(Ri ~ yl + shiftage+Female+Educyrs+iadl1+HlthPrv1+mmse1, family = binomial(link=logit),data = ipwdat882)
logistfit4 = glm(Ri ~ yl + shiftage+Female+Educyrs+iadl1+HlthPrv1+mmse1, family = binomial(link=logit),data = ipwdat884)
logistfit6 = glm(Ri ~ yl + shiftage+Female+Educyrs+iadl1+HlthPrv1+mmse1, family = binomial(link=logit),data = ipwdat886)
logistfit8 = glm(Ri ~ yl + shiftage+Female+Educyrs+iadl1+HlthPrv1+mmse1, family = binomial(link=logit),data = ipwdat888)
ipwdat882$pred_obs = predict(logistfit2, newdata = ipwdat82, type="response")
ipwdat884$pred_obs = predict(logistfit4, newdata = ipwdat84, type="response")
ipwdat886$pred_obs = predict(logistfit6, newdata = ipwdat86, type="response")
ipwdat888$pred_obs = predict(logistfit8, newdata = ipwdat88, type="response")
##PREDICT for waves

combine = rbind(ipwdat22,ipwdat442,ipwdat444,ipwdat662,ipwdat664,ipwdat666,ipwdat882,ipwdat884,ipwdat886,ipwdat888)
combine = combine[order(combine$Case),]; 
combine$Rim1=combine$Ri=combine$yl=combine$lyl=combine$llyl=combine$lllyl=NULL;
combw = reshape(combine, idvar = "Case", v.names=c("pek","pred_obs"),timevar="time", direction="wide")
combw$pi.2=1/combw$pred_obs.2; combw$pi.4=1/(combw$pred_obs.2*combw$pred_obs.4);
combw$pi.6=1/(combw$pred_obs.2*combw$pred_obs.4*combw$pred_obs.6); combw$pi.8=1/(combw$pred_obs.2*combw$pred_obs.4*combw$pred_obs.6*combw$pred_obs.8);
combw$pred_obs.2 = combw$pred_obs.4 = combw$pred_obs.6 = combw$pred_obs.8 = NULL;

long = reshape(combw, varying=list(c(10:13),c(14:17)), v.names = c("pek","pi"), new.row.names=1:1000000, times=c(2,4,6,8), direction="long")
long$id=NULL;
ipwdat0 = ipwdat[ipwdat$time==0,]; ipwdat0$pi=1; 
ipwdat0$yl=ipwdat0$lyl=ipwdat0$llyl=ipwdat0$lllyl=ipwdat0$pred_obs=ipwdat0$Ri=ipwdat0$Rim1=NULL;
ripwdatd = rbind(ipwdat0,long)

ripwdatd = ripwdatd[order(ripwdatd$Case),];
ripwdatd = ripwdatd[ripwdatd$time<= ripwdatd$s1,]; ripwdatd$Ri=ifelse(!is.na(ripwdatd$pek),1,0);

attach(ripwdatd)
ipw.fitd = geeglm(pek ~ time +  I(time^2) + Female + shiftage + Educyrs + Smoke + time*Female + time*shiftage + time*Educyrs + time*Smoke, data = ripwdatd, weights = pi, id=Case)

###Codes for analyzing data under MI (without d)###
###################################################
w = reshape(tmpdata, idvar = "Case", v.names=c("pek"),timevar="time", direction="wide")
midat = w[, c("Female","Educyrs","shiftage","pek.0","pek.2","pek.4","pek.6","pek.8","iadl1","Smoke","HlthPrv1","mmse1")]

imdata = as.matrix(midat)

nsim = 30
PE = VAR = matrix(0,nsim,11)
rngseed(112789)
for(i in 1:nsim)
{
  s = prelim.norm(imdata)
  thetahat = em.norm(s,showits=F)
  ximp = imp.norm(s, thetahat, imdata)
  tmp = as.data.frame(ximp); tmp$s1=w$s1;
  
  long = reshape(tmp, varying=list(c(4:8)), v.names = c("pek"), new.row.names=1:1000000, times=c(0,2,4,6,8), direction="long")
  longsort = long[order(long$id),];
  midat = longsort[longsort$time<= longsort$s1,]
  
  fitmi = geeglm(pek ~  time +  I(time^2) + Female + shiftage + Educyrs + Smoke + time*Female + time*shiftage + time*Educyrs + time*Smoke, data = midat, id=id,na.action=na.omit)
  PE[i,]=(fitmi)$coef
  VAR[i,] = summary(fitmi)$coef[,2]^2
}

COEF = apply(PE, 2, mean)
U1 = apply(VAR, 2, mean)
B1 = apply(PE, 2, var)
SE = (U1+(1+1/nsim)*B1)^0.5
RES = cbind(COEF, SE)

Wald = RES[,1]^2/RES[,2]^2
pvalues <- pchisq(Wald, 1, lower=F)
cbind(COEF, SE, pvalues)

###Coes to analyze data using MI(strat D)###
############################################

wd = reshape(tmpdata, idvar = "Case", v.names=c("pek"),timevar="time", direction="wide")
wd1=wd[wd$s1==0,];wd2 = wd[wd$s1==2,];wd3 = wd[wd$s1==4,];
wd4 = wd[wd$s1==6,];wd5 = wd[wd$s1==8,];
midatd1 = wd1[, c("Female","Educyrs","shiftage","pek.0","pek.2","pek.4","pek.6","pek.8","s1","iadl1","Smoke","HlthPrv1","mmse1")]
midatd2 = cbind(wd2[, c("Female","Educyrs","shiftage","pek.0","pek.2","iadl1","Smoke","HlthPrv1","mmse1")])
midatd3 = cbind(wd3[, c("Female","Educyrs","shiftage","pek.0","pek.2","pek.4","iadl1","Smoke","HlthPrv1","mmse1")]) 
midatd4 = cbind(wd4[, c("Female","Educyrs","shiftage","pek.0","pek.2","pek.4","pek.6","iadl1","Smoke","HlthPrv1","mmse1")])
midatd5 = cbind(wd5[, c("Female","Educyrs","shiftage","pek.0","pek.2","pek.4","pek.6","pek.8","iadl1","Smoke","HlthPrv1","mmse1")])

imdatad2 = as.matrix(midatd2)
imdatad3 = as.matrix(midatd3); 
imdatad4 = as.matrix(midatd4); 
imdatad5 = as.matrix(midatd5);  

nsim = 30
PE = VAR = matrix(0,nsim,11)
rngseed(112789)
for(i in 1:nsim)
{
  s = prelim.norm(imdatad2)
  thetahat = em.norm(s,showits=F)
  ximp = imp.norm(s, thetahat, imdatad2)
  tmp = as.data.frame(ximp); tmp$pek.4=tmp$pek.6=tmp$pek.8=NA; tmp$s1=wd2$s1;
  
  mydat = rbind(midatd1, tmp)
  
  s = prelim.norm(imdatad3)
  thetahat = em.norm(s,showits=F)
  ximp = imp.norm(s, thetahat, imdatad3)
  tmp = as.data.frame(ximp); tmp$pek.6=tmp$pek.8=NA; tmp$s1=wd3$s1;
  
  mydat = rbind(mydat, tmp)
  
  s = prelim.norm(imdatad4)
  thetahat = em.norm(s,showits=F)
  ximp = imp.norm(s, thetahat, imdatad4)
  tmp = as.data.frame(ximp); tmp$pek.8=NA; tmp$s1=wd4$s1;
  
  mydat = rbind(mydat, tmp)
  
  s = prelim.norm(imdatad5)
  thetahat = em.norm(s,showits=F)
  ximp = imp.norm(s, thetahat, imdatad5)
  tmp = as.data.frame(ximp); tmp$s1=wd5$s1;
  
  mydat = rbind(mydat, tmp)
  
  long = reshape(mydat, varying=list(c(4:8)), v.names = c("pek"), new.row.names=1:1000000, times=c(0,2,4,6,8), direction="long")
  longsort = long[order(long$id),];
  midatd = longsort[longsort$time<= longsort$s1,]
  
  fitmi = geeglm(pek ~  time +  I(time^2) + Female + shiftage + Educyrs + Smoke + time*Female + time*shiftage + time*Educyrs + time*Smoke, data = midatd, id=id,na.action=na.omit)
  PE[i,]=(fitmi)$coef
  VAR[i,] = summary(fitmi)$coef[,2]^2
  
}

COEFd = apply(PE, 2, mean)
U1 = apply(VAR, 2, mean)
B1 = apply(PE, 2, var)
SEd = (U1+(1+1/nsim)*B1)^0.5
RESd = cbind(COEFd, SEd)

Wald = RESd[,1]^2/RESd[,2]^2
pvalues <- pchisq(Wald, 1, lower=F)
cbind(COEFd, SEd, pvalues)


###Codes to analyze data using AIPW (not conditioning on D)###
##############################################################
dat = tmpdata
dat$yl = c(NA, dat$pek[1:nrow(dat)-1])
dat$yl = ifelse(dat$time==0, NA, dat$yl)
dat$lyl = c(NA, dat$yl[1:nrow(dat)-1])
dat$lyl = ifelse(dat$time<=2, NA, dat$lyl)
dat$llyl = c(NA, dat$lyl[1:nrow(dat)-1])
dat$llyl = ifelse(dat$time<=4, NA, dat$llyl)
dat$lllyl = c(NA, dat$llyl[1:nrow(dat)-1])
dat$lllyl = ifelse(dat$time<=6, NA, dat$lllyl)
dat$Ri = ifelse(!is.na(dat$pek),1,0)
dat$Rim1 = c(NA, dat$Ri[1:nrow(dat)-1])   ##Status of the last observation
dat$Rim1 = ifelse(dat$time==0, NA, dat$Rim1)

midat12 = dat[dat$time==2 & dat$Rim1==1,]
mifit1 = lm(pek  ~ yl + shiftage + Educyrs + Female +iadl1+HlthPrv1+mmse1+Smoke, data = midat12)
mi12 = midat12[,c("yl", "shiftage", "Female","Educyrs","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
midat12$pred_obs = predict(mifit1, newdata = mi12, type="response")
midat12$pred_obs = ifelse(!is.na(midat12$pek), midat12$pek, midat12$pred_obs)
midat12$pek = midat12$pred_obs; 

midat23 = dat[dat$time==4 & dat$Rim1==1,]
mifit2 = lm(pek  ~ yl + lyl + shiftage + Educyrs+Female +iadl1+HlthPrv1+mmse1+Smoke, data = midat23)
mi23 = midat23[,c("yl", "lyl","shiftage", "Female","Educyrs","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
midat23$pred_obs = predict(mifit2, newdata = mi23, type="response")
midat23$pred_obs = ifelse(!is.na(midat23$pek), midat23$pek, midat23$pred_obs)

midat34 = dat[dat$time==6 & dat$Rim1==1,]
mifit3 = lm(pek  ~ yl + lyl + llyl + shiftage + Educyrs+Female +iadl1+HlthPrv1+mmse1+Smoke, data = midat34)
mi34 = midat34[,c("yl", "lyl","llyl", "shiftage", "Female","Educyrs","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
midat34$pred_obs = predict(mifit3, newdata = mi34, type="response")
midat34$pred_obs = ifelse(!is.na(midat34$pek), midat34$pek, midat34$pred_obs)

midat45 = dat[dat$time==8 & dat$Rim1==1,]
mifit4 = lm(pek  ~ yl + lyl + llyl + lllyl + shiftage+Educyrs+Female +iadl1+HlthPrv1+mmse1+Smoke, data = midat45)
mi45 = midat45[,c("yl", "lyl","llyl","lllyl", "shiftage", "Female","Educyrs","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
midat45$pred_obs = predict(mifit4, newdata = mi45, type="response")
midat45$pred_obs = ifelse(!is.na(midat45$pek), midat45$pek, midat45$pred_obs)

###########################Fitted values are not covariates but outcomes
tmpmidat23 = dat[dat$time==4 & is.na(dat$yl),]
midat23$pek = midat23$pred_obs; midat23$pred_obs=NULL; midat23 = rbind(tmpmidat23,midat23);
midat23$Ri = ifelse(!is.na(midat23$pek),1,0)
mifit5 = lm(pek  ~  lyl + shiftage + Educyrs+Female +iadl1+HlthPrv1+mmse1+Smoke, data = midat23)
mi23 = midat23[,c("lyl", "shiftage", "Female","Educyrs","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
midat23$pred_obs = predict(mifit5, newdata = mi23, type="response")
midat23$pred_obs = ifelse(!is.na(midat23$pek), midat23$pek, midat23$pred_obs)
midat23$pek = midat23$pred_obs; 

tmpmidat34 = dat[dat$time==6 & is.na(dat$yl) & !is.na(dat$lyl),]
midat34$pek = midat34$pred_obs; midat34$pred_obs=NULL; midat34 = rbind(tmpmidat34,midat34);
midat34$Ri = ifelse(!is.na(midat34$pek),1,0)
mifit6 = lm(pek  ~ lyl + llyl + shiftage + Educyrs+Female +iadl1+HlthPrv1+mmse1+Smoke, data = midat34)
mi34 = midat34[,c("lyl","llyl", "shiftage", "Female","Educyrs","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
midat34$pred_obs = predict(mifit6, newdata = mi34, type="response")
midat34$pred_obs = ifelse(!is.na(midat34$pek), midat34$pek, midat34$pred_obs)

tmpmidat45 = dat[dat$time==8 & is.na(dat$yl) & !is.na(dat$lyl),]
midat45$pek = midat45$pred_obs; midat45$pred_obs=NULL; midat45 = rbind(tmpmidat45,midat45);
midat45$Ri = ifelse(!is.na(midat45$pek),1,0)
mifit7 = lm(pek  ~ lyl + llyl + lllyl+ shiftage + Educyrs+Female +iadl1+HlthPrv1+mmse1+Smoke, data = midat45)
mi45 = midat45[,c("lyl","llyl","lllyl" , "shiftage", "Female","Educyrs","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
midat45$pred_obs = predict(mifit7, newdata = mi45, type="response")
midat45$pred_obs = ifelse(!is.na(midat45$pek), midat45$pek, midat45$pred_obs)

######
tmpmidat34 = dat[dat$time==6 & is.na(dat$yl) & is.na(dat$lyl),]
midat34$pek = midat34$pred_obs; midat34$pred_obs=NULL; midat34 = rbind(tmpmidat34,midat34);
midat34$Ri = ifelse(!is.na(midat34$pek),1,0)
mifit8 = lm(pek  ~ llyl + shiftage + Educyrs+Female +iadl1+HlthPrv1+mmse1+Smoke, data = midat34)
mi34 = midat34[,c("llyl", "shiftage", "Female","Educyrs","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
midat34$pred_obs = predict(mifit8, newdata = mi34, type="response")
midat34$pred_obs = ifelse(!is.na(midat34$pek), midat34$pek, midat34$pred_obs)
midat34$pek = midat34$pred_obs; 

tmpmidat45 = dat[dat$time==8 & is.na(dat$yl) & is.na(dat$lyl) & !is.na(dat$llyl),]
midat45$pek = midat45$pred_obs; midat45$pred_obs=NULL; midat45 = rbind(tmpmidat45,midat45);
midat45$Ri = ifelse(!is.na(midat45$pek),1,0)
mifit9 = lm(pek  ~  llyl + lllyl + shiftage + Educyrs+Female +iadl1+HlthPrv1+mmse1+Smoke, data = midat45)
mi45 = midat45[,c("llyl","lllyl", "shiftage", "Female","Educyrs","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
midat45$pred_obs = predict(mifit9, newdata = mi45, type="response")
midat45$pred_obs = ifelse(!is.na(midat45$pek), midat45$pek, midat45$pred_obs)

######
tmpmidat45 = dat[dat$time==8 & is.na(dat$yl) & is.na(dat$llyl),]
midat45$pek = midat45$pred_obs; midat45$pred_obs=NULL; midat45 = rbind(tmpmidat45,midat45);
midat45$Ri = ifelse(!is.na(midat45$pek),1,0)
mifit10 = lm(pek  ~ lllyl + shiftage + Educyrs+Female +iadl1+HlthPrv1+mmse1+Smoke, data = midat45)
mi45 = midat45[,c("lllyl", "shiftage", "Female","Educyrs","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
midat45$pred_obs = predict(mifit10, newdata = mi45, type="response")
midat45$pred_obs = ifelse(!is.na(midat45$pek), midat45$pek, midat45$pred_obs)
midat45$pek = midat45$pred_obs; 

#####
midat = rbind(midat12, midat23, midat34,midat45); midat$pred_obs = NULL;
midat = rbind(dat[dat$time==0,],midat)
midatsort = midat[order(midat$Case),]
midatsort = midatsort[midatsort$time<=midatsort$s1,]

########## complete data set
midatsort$pi=1/ripwdat$pi
midatsort$Rim1 = midatsort$yl = midatsort$lyl = midatsort$llyl = midatsort$lllyl = NULL;
midatsort$Ri = ripwdat$Ri; midatsort$digreal = ripwdat$pek;
##NR##
beta = summary(ipw.fit)$coef[,1];
mybeta = NR(crit=0.00001,beta=beta)



###Codes to analyze data using AIPW (stratify on D)###
######################################################
dat = tmpdata
dat = dat[dat$s1>=dat$time,]
dat$yl = c(NA, dat$pek[1:nrow(dat)-1])
dat$yl = ifelse(dat$time==0, NA, dat$yl)
dat$lyl = c(NA, dat$yl[1:nrow(dat)-1])
dat$lyl = ifelse(dat$time<=2, NA, dat$lyl)
dat$llyl = c(NA, dat$lyl[1:nrow(dat)-1])
dat$llyl = ifelse(dat$time<=4, NA, dat$llyl)
dat$lllyl = c(NA, dat$llyl[1:nrow(dat)-1])
dat$lllyl = ifelse(dat$time<=6, NA, dat$lllyl)
dat$Ri = ifelse(!is.na(dat$pek),1,0)
dat$Rim1 = c(NA, dat$Ri[1:nrow(dat)-1])   ##Status of the last observation
dat$Rim1 = ifelse(dat$time==0, NA, dat$Rim1)

midat12 = dat[dat$time==2 & dat$Rim1==1 & dat$s1==2,]
mifit1 = lm(pek  ~ yl + shiftage + Educyrs+Female+iadl1+HlthPrv1+mmse1+Smoke, data = midat12)
mi12 = midat12[,c("yl", "shiftage", "Educyrs","Female","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
midat12$pred_obs = predict(mifit1, newdata = mi12, type="response")
midat12$pred_obs = ifelse(!is.na(midat12$pek), midat12$pek, midat12$pred_obs)
midat12$pek = midat12$pred_obs; 

midat122 = dat[dat$time==2 & dat$Rim1==1 & dat$s1==4,]
mifit12 = lm(pek  ~ yl + shiftage + Educyrs+Female+iadl1+HlthPrv1+mmse1+Smoke, data = midat122)
mi12 = midat122[,c("yl", "shiftage","Educyrs","Female","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
midat122$pred_obs = predict(mifit12, newdata = mi12, type="response")
midat122$pred_obs = ifelse(!is.na(midat122$pek), midat122$pek, midat122$pred_obs)
midat122$pek = midat122$pred_obs; 

midat123 = dat[dat$time==2 & dat$Rim1==1 & dat$s1==6,]
mifit13 = lm(pek  ~ yl + shiftage + Educyrs+Female+iadl1+HlthPrv1+mmse1+Smoke, data = midat123)
mi12 = midat123[,c("yl", "shiftage", "Educyrs","Female","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
midat123$pred_obs = predict(mifit13, newdata = mi12, type="response")
midat123$pred_obs = ifelse(!is.na(midat123$pek), midat123$pek, midat123$pred_obs)
midat123$pek = midat123$pred_obs; 

midat124 = dat[dat$time==2 & dat$Rim1==1 & dat$s1==8,]
mifit14 = lm(pek  ~ yl + shiftage + Educyrs+Female+iadl1+HlthPrv1+mmse1+Smoke, data = midat124)
mi12 = midat124[,c("yl", "shiftage","Educyrs","Female","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
midat124$pred_obs = predict(mifit14, newdata = mi12, type="response")
midat124$pred_obs = ifelse(!is.na(midat124$pek), midat124$pek, midat124$pred_obs)
midat124$pek = midat124$pred_obs; 

midat23 = dat[dat$time==4 & dat$Rim1==1 & dat$s1==4,]
mifit2 = lm(pek  ~ yl + lyl + shiftage + Educyrs+Female+iadl1+HlthPrv1+mmse1+Smoke, data = midat23)
mi23 = midat23[,c("yl", "lyl","shiftage", "Educyrs","Female","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
midat23$pred_obs = predict(mifit2, newdata = mi23, type="response")
midat23$pred_obs = ifelse(!is.na(midat23$pek), midat23$pek, midat23$pred_obs)

midat232 = dat[dat$time==4 & dat$Rim1==1 & dat$s1==6,]
mifit22 = lm(pek  ~ yl + lyl + shiftage + Educyrs+Female+iadl1+HlthPrv1+mmse1+Smoke, data = midat232)
mi23 = midat232[,c("yl", "lyl","shiftage","Educyrs","Female","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
midat232$pred_obs = predict(mifit22, newdata = mi23, type="response")
midat232$pred_obs = ifelse(!is.na(midat232$pek), midat232$pek, midat232$pred_obs)

midat233 = dat[dat$time==4 & dat$Rim1==1 & dat$s1==8,]
mifit23 = lm(pek  ~ yl + lyl + shiftage + Educyrs+Female+iadl1+HlthPrv1+mmse1+Smoke, data = midat233)
mi23 = midat233[,c("yl", "lyl","shiftage", "Educyrs","Female","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
midat233$pred_obs = predict(mifit23, newdata = mi23, type="response")
midat233$pred_obs = ifelse(!is.na(midat233$pek), midat233$pek, midat233$pred_obs)

midat34 = dat[dat$time==6 & dat$Rim1==1 & dat$s1==6,]
mifit3 = lm(pek  ~ yl + lyl + llyl + shiftage + Educyrs+Female+iadl1+HlthPrv1+mmse1+Smoke, data = midat34)
mi34 = midat34[,c("yl", "lyl","llyl", "shiftage","Educyrs","Female","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
midat34$pred_obs = predict(mifit3, newdata = mi34, type="response")
midat34$pred_obs = ifelse(!is.na(midat34$pek), midat34$pek, midat34$pred_obs)

midat342 = dat[dat$time==6 & dat$Rim1==1 & dat$s1==8,]
mifit32 = lm(pek  ~ yl + lyl + llyl + shiftage + Educyrs+Female+iadl1+HlthPrv1+mmse1+Smoke, data = midat342)
mi34 = midat342[,c("yl", "lyl","llyl", "shiftage","Educyrs","Female","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
midat342$pred_obs = predict(mifit32, newdata = mi34, type="response")
midat342$pred_obs = ifelse(!is.na(midat342$pek), midat342$pek, midat342$pred_obs)

midat45 = dat[dat$time==8 & dat$Rim1==1,]
mifit4 = lm(pek  ~ yl + lyl + llyl + lllyl + shiftage + Educyrs+Female+iadl1+HlthPrv1+mmse1+Smoke, data = midat45)
mi45 = midat45[,c("yl", "lyl","llyl","lllyl", "shiftage","Educyrs","Female","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
midat45$pred_obs = predict(mifit4, newdata = mi45, type="response")
midat45$pred_obs = ifelse(!is.na(midat45$pek), midat45$pek, midat45$pred_obs)

###########################Fitted values are not covariates but outcomes
tmpmidat23 = dat[dat$time==4 & is.na(dat$yl) & dat$s1==4,]
midat23$pek = midat23$pred_obs; midat23$pred_obs=NULL; midat23 = rbind(tmpmidat23,midat23);
midat23$Ri = ifelse(!is.na(midat23$pek),1,0)
mifit5 = lm(pek  ~  lyl + shiftage + Educyrs+Female+iadl1+HlthPrv1+mmse1+Smoke, data = midat23)
mi23 = midat23[,c("lyl", "shiftage","Educyrs","Female","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
midat23$pred_obs = predict(mifit5, newdata = mi23, type="response")
midat23$pred_obs = ifelse(!is.na(midat23$pek), midat23$pek, midat23$pred_obs)
midat23$pek = midat23$pred_obs; 

tmpmidat23 = dat[dat$time==4 & is.na(dat$yl) & dat$s1==6,]
midat232$pek = midat232$pred_obs; midat232$pred_obs=NULL; midat232 = rbind(tmpmidat23,midat232);
midat232$Ri = ifelse(!is.na(midat232$pek),1,0)
mifit52 = lm(pek  ~  lyl + shiftage + Educyrs+Female+iadl1+HlthPrv1+mmse1+Smoke, data = midat232)
mi23 = midat232[,c("lyl", "shiftage","Educyrs","Female","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
midat232$pred_obs = predict(mifit52, newdata = mi23, type="response")
midat232$pred_obs = ifelse(!is.na(midat232$pek), midat232$pek, midat232$pred_obs)
midat232$pek = midat232$pred_obs;

tmpmidat23 = dat[dat$time==4 & is.na(dat$yl) & dat$s1==8,]
midat233$pek = midat233$pred_obs; midat233$pred_obs=NULL; midat233 = rbind(tmpmidat23,midat233);
midat233$Ri = ifelse(!is.na(midat233$pek),1,0)
mifit53 = lm(pek  ~  lyl + shiftage + Educyrs+Female+iadl1+HlthPrv1+mmse1+Smoke, data = midat233)
mi23 = midat233[,c("lyl", "shiftage", "Educyrs","Female","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
midat233$pred_obs = predict(mifit53, newdata = mi23, type="response")
midat233$pred_obs = ifelse(!is.na(midat233$pek), midat233$pek, midat233$pred_obs)
midat233$pek = midat233$pred_obs; 

tmpmidat34 = dat[dat$time==6 & is.na(dat$yl) & !is.na(dat$lyl) & dat$s1==6,]
midat34$pek = midat34$pred_obs; midat34$pred_obs=NULL; midat34 = rbind(tmpmidat34,midat34);
midat34$Ri = ifelse(!is.na(midat34$pek),1,0)
mifit6 = lm(pek  ~ lyl + llyl + shiftage + Educyrs+Female+iadl1+HlthPrv1+mmse1+Smoke, data = midat34)
mi34 = midat34[,c("lyl","llyl", "shiftage", "Educyrs","Female","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
midat34$pred_obs = predict(mifit6, newdata = mi34, type="response")
midat34$pred_obs = ifelse(!is.na(midat34$pek), midat34$pek, midat34$pred_obs)

tmpmidat342 = dat[dat$time==6 & is.na(dat$yl) & !is.na(dat$lyl) & dat$s1==8,]
midat342$pek = midat342$pred_obs; midat342$pred_obs=NULL; midat342 = rbind(tmpmidat342,midat342);
midat342$Ri = ifelse(!is.na(midat342$pek),1,0)
mifit62 = lm(pek  ~ lyl + llyl + shiftage + Educyrs+Female+iadl1+HlthPrv1+mmse1+Smoke, data = midat342)
mi34 = midat342[,c("lyl","llyl", "shiftage","Educyrs","Female","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
midat342$pred_obs = predict(mifit62, newdata = mi34, type="response")
midat342$pred_obs = ifelse(!is.na(midat342$pek), midat342$pek, midat342$pred_obs)

tmpmidat45 = dat[dat$time==8 & is.na(dat$yl) & !is.na(dat$lyl),]
midat45$pek = midat45$pred_obs; midat45$pred_obs=NULL; midat45 = rbind(tmpmidat45,midat45);
midat45$Ri = ifelse(!is.na(midat45$pek),1,0)
mifit7 = lm(pek  ~ lyl + llyl + lllyl+ shiftage + Educyrs+Female+iadl1+HlthPrv1+mmse1+Smoke, data = midat45)
mi45 = midat45[,c("lyl","llyl","lllyl" , "shiftage", "Educyrs","Female","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
midat45$pred_obs = predict(mifit7, newdata = mi45, type="response")
midat45$pred_obs = ifelse(!is.na(midat45$pek), midat45$pek, midat45$pred_obs)

######
tmpmidat34 = dat[dat$time==6 & is.na(dat$yl) & is.na(dat$lyl)  & !is.na(dat$llyl) & dat$s1==6,]
midat34$pek = midat34$pred_obs; midat34$pred_obs=NULL; midat34 = rbind(tmpmidat34,midat34);
midat34$Ri = ifelse(!is.na(midat34$pek),1,0)
mifit8 = lm(pek  ~ llyl + shiftage + Educyrs+Female+iadl1+HlthPrv1+mmse1+Smoke, data = midat34)
mi34 = midat34[,c("llyl", "shiftage","Educyrs","Female","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
midat34$pred_obs = predict(mifit8, newdata = mi34, type="response")
midat34$pred_obs = ifelse(!is.na(midat34$pek), midat34$pek, midat34$pred_obs)
midat34$pek = midat34$pred_obs; 

tmpmidat342 = dat[dat$time==6 & is.na(dat$yl) & is.na(dat$lyl)  & !is.na(dat$llyl) & dat$s1==8,]
midat342$pek = midat342$pred_obs; midat342$pred_obs=NULL; midat342 = rbind(tmpmidat342,midat342);
midat342$Ri = ifelse(!is.na(midat342$pek),1,0)
mifit82 = lm(pek  ~ llyl + shiftage + Educyrs+Female+iadl1+HlthPrv1+mmse1+Smoke, data = midat342)
mi34 = midat342[,c("llyl", "shiftage","Educyrs","Female","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
midat342$pred_obs = predict(mifit82, newdata = mi34, type="response")
midat342$pred_obs = ifelse(!is.na(midat342$pek), midat342$pek, midat342$pred_obs)
midat342$pek = midat342$pred_obs; 

tmpmidat45 = dat[dat$time==8 & is.na(dat$yl) & is.na(dat$lyl) & !is.na(dat$llyl),]
midat45$pek = midat45$pred_obs; midat45$pred_obs=NULL; midat45 = rbind(tmpmidat45,midat45);
midat45$Ri = ifelse(!is.na(midat45$pek),1,0)
mifit9 = lm(pek  ~  llyl + lllyl + shiftage + Educyrs+Female+iadl1+HlthPrv1+mmse1+Smoke, data = midat45)
mi45 = midat45[,c("llyl","lllyl", "shiftage","Educyrs","Female","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
midat45$pred_obs = predict(mifit9, newdata = mi45, type="response")
midat45$pred_obs = ifelse(!is.na(midat45$pek), midat45$pek, midat45$pred_obs)

######
tmpmidat45 = dat[dat$time==8 & is.na(dat$yl) & is.na(dat$llyl),]
midat45$pek = midat45$pred_obs; midat45$pred_obs=NULL; midat45 = rbind(tmpmidat45,midat45);
midat45$Ri = ifelse(!is.na(midat45$pek),1,0)
mifit10 = lm(pek  ~ lllyl + shiftage + Educyrs+Female+iadl1+HlthPrv1+mmse1+Smoke, data = midat45)
mi45 = midat45[,c("lllyl", "shiftage","Educyrs","Female","time","s1","iadl1","HlthPrv1","mmse1","Smoke")]
midat45$pred_obs = predict(mifit10, newdata = mi45, type="response")
midat45$pred_obs = ifelse(!is.na(midat45$pek), midat45$pek, midat45$pred_obs)
midat45$pek = midat45$pred_obs; 


#####
midat = rbind(midat12, midat122, midat123, midat124, midat23, midat232, midat233, midat34, midat342, midat45); 
midat$pred_obs = NULL;
midat = rbind(dat[dat$time==0,],midat)
midatsortd = midat[order(midat$Case),]
midatsortd = midatsortd[midatsortd$time<=midatsortd$s1,]

########## complete data set
midatsortd$pi=1/ripwdatd$pi
midatsortd$Rim1 = midatsortd$yl = midatsortd$lyl = midatsortd$llyl = midatsortd$lllyl=NULL;
midatsortd$Ri = ripwdatd$Ri; midatsortd$digreal = ripwdatd$pek;
##NR##

beta = summary(ipw.fitd)$coef[,1];
mybetaD = NRd(crit=0.00001,beta=beta)
