library(doParallel)
library(foreach)

# Calculate the number of cores
getDoParWorkers()                                          
detectCores()                                                      
cl=makeCluster(detectCores()-2)                                              
registerDoParallel(cl)                                                                          
getDoParWorkers()   

#for(m in 1:sim)
myfunc = function(m)
{
  library(geepack)
  source("simdatDAG3LMD.R")
  set.seed(1127)
  seeds = floor(runif(1000)*10^8);
  
  EXPIT <- function(term) {
    return( exp(term)/(1+exp(term)) )
  }
  
  set.seed(seeds[m])
tmpdata = gendata(N=500)
tmpdata$ager = ifelse(tmpdata$shiftage<=4,1,0);
###IPWIEE data
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

###IPW-IEE-D
ipwdat = dat[,c("Case", "time", "shiftage", "sex", "edu", "Ri","pek","s1","Rim1","yl", "lyl","llyl","lllyl","ager")]
ipwdat$pred_obs=NA;ipwdat$pi=NA;

ipwdat22 = ipwdat[ipwdat$time==2 & ipwdat$s1>=2 & !is.na(ipwdat$Rim1) & ipwdat$Rim1==1,]; 
ipwdat2 = ipwdat22[,c("yl", "shiftage", "sex","edu","time","s1")]
logistfit = glm(Ri ~ yl +sex, family = binomial(link=logit),data = ipwdat22)
ipwdat22$pred_obs = predict(logistfit, newdata = ipwdat2, type="response")

##########

ipwdat44 = ipwdat[ipwdat$s1>=4 & !is.na(ipwdat$Rim1) & ipwdat$Rim1==1,];
ipwdat442 = ipwdat44[ipwdat44$time==2,]; ipwdat444 = ipwdat44[ipwdat44$time==4,]; 
ipwdat42 = ipwdat442[,c("yl", "lyl","llyl", "shiftage", "sex","edu","time","s1")]
ipwdat44 = ipwdat444[,c("yl", "lyl","llyl", "shiftage", "sex","edu","time","s1")]
logistfit2 = glm(Ri ~ yl + sex, family = binomial(link=logit),data = ipwdat442)
logistfit4 = glm(Ri ~ yl + sex, family = binomial(link=logit),data = ipwdat444)
ipwdat442$pred_obs = predict(logistfit2, newdata = ipwdat42, type="response")
ipwdat444$pred_obs = predict(logistfit4, newdata = ipwdat44, type="response")

###############

ipwdat66 = ipwdat[ipwdat$s1>=6 & !is.na(ipwdat$Rim1) & ipwdat$Rim1==1,];
ipwdat662 = ipwdat66[ipwdat66$time==2,]; ipwdat664 = ipwdat66[ipwdat66$time==4,]; ipwdat666 = ipwdat66[ipwdat66$time==6,];
ipwdat62 = ipwdat662[,c("yl", "lyl","llyl","shiftage", "sex","edu","time","s1")]
ipwdat64 = ipwdat664[,c("yl", "lyl","llyl","shiftage", "sex","edu","time","s1")]
ipwdat66 = ipwdat666[,c("yl", "lyl","llyl","shiftage", "sex","edu","time","s1")]
logistfit2 = glm(Ri ~ yl+sex, family = binomial(link=logit),data = ipwdat662)
logistfit4 = glm(Ri ~ yl+sex, family = binomial(link=logit),data = ipwdat664)
logistfit6 = glm(Ri ~ yl+sex, family = binomial(link=logit),data = ipwdat666)
ipwdat662$pred_obs = predict(logistfit2, newdata = ipwdat62, type="response")
ipwdat664$pred_obs = predict(logistfit4, newdata = ipwdat64, type="response")
ipwdat666$pred_obs = predict(logistfit6, newdata = ipwdat66, type="response")

#################

ipwdat88 = ipwdat[ipwdat$s1==8 & !is.na(ipwdat$Rim1) & ipwdat$Rim1==1,]
ipwdat882 = ipwdat88[ipwdat88$time==2,]; ipwdat884 = ipwdat88[ipwdat88$time==4,]; ipwdat886 = ipwdat88[ipwdat88$time==6,]; ipwdat888 = ipwdat88[ipwdat88$time==8,];
ipwdat82 = ipwdat882[,c("yl", "lyl","llyl","lllyl", "shiftage", "sex","edu","time","s1")]
ipwdat84 = ipwdat884[,c("yl", "lyl","llyl","lllyl", "shiftage", "sex","edu","time","s1")]
ipwdat86 = ipwdat886[,c("yl", "lyl","llyl","lllyl", "shiftage", "sex","edu","time","s1")]
ipwdat88 = ipwdat888[,c("yl", "lyl","llyl","lllyl", "shiftage", "sex","edu","time","s1")]
logistfit2 = glm(Ri ~ yl + sex, family = binomial(link=logit),data = ipwdat882)
logistfit4 = glm(Ri ~ yl + sex, family = binomial(link=logit),data = ipwdat884)
logistfit6 = glm(Ri ~ yl + sex, family = binomial(link=logit),data = ipwdat886)
logistfit8 = glm(Ri ~ yl + sex, family = binomial(link=logit),data = ipwdat888)
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
ipw.fitp = geeglm(pek ~  factor(time) , data = ripwdat, weights = pi, id=Case)

outtmp = c(
  summary(ipw.fitp)$coefficient[,1])
outtmp
}
test = foreach(m=1:1000) %dopar% myfunc(m)
test2 = do.call("rbind", test)

write.csv(test2,"ipwp.csv")

