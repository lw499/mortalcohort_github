setwd("~/Dropbox/phd/Missing data comparisons/pek/simwithoutedu_tonly/mortalcohort_github/demo_data")

tmpdata = read.csv("demodat.csv", sep=",")
fixdatw = reshape(tmpdata, idvar = "Case", v.names=c("pek"),timevar="time", direction="wide")

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
  set.seed(30)
  seeds = floor(runif(1000)*10^8);
  
  EXPIT <- function(term) {
    return( exp(term)/(1+exp(term)) )}
  ###Bootstrap
  set.seed(seeds[m])
  random = fixdatw[sample(1:nrow(fixdatw), 500, replace=T),]
  longdat = reshape(random, varying=list(c(10:14)), v.names = c("pek"), new.row.names=1:100000, times=c(0,2,4,6,8), direction="long")
  dat = longdat[order(longdat$id),]; id=dat$id; dat$id=NULL;
  shiftage = dat$shiftage; s1=dat$s1; pek=dat$pek; dat$Case=id;
  dat$shiftage=NULL; dat$s1=NULL; dat$pek=NULL;
  dat$pek=pek; dat$shiftage=shiftage; dat$s1=s1;
  
  #dat=tmpdata;
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

  return(c(summary(ipw.fit)$coef[,1]))
}
test = foreach(m=1:1000) %dopar% myfunc(m)
test2 = do.call("rbind", test)

write.csv(test2,"listipw.csv")



