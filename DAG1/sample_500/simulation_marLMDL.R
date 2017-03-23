## Results from IPW methods, MI, LI, AIPW (without stratification)

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
  source("linearincrements.r")
  source("simdatDAG1LMD.R")
  set.seed(60)
  seeds = floor(runif(1000)*10^8);
  
  EXPIT <- function(term) {
    return( exp(term)/(1+exp(term)) )
  }
  
  set.seed(seeds[m])
  tmpdata = gendata(N=500)
  tmpdata$ager = ifelse(tmpdata$shiftage<=4,1,0);  
  
  fit = geeglm(pek ~  factor(time)+ factor(time)*sex,id=Case,data=tmpdata,corstr="ind")
  
  

###IPWIEE
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

ipwdat = dat[,c("Case", "time", "shiftage", "sex", "edu", "Ri","pek","s1","Rim1","yl", "lyl","llyl","lllyl","ager")]
ipwdat$pred_obs=NA;ipwdat$pi=NA;

ipwdat22 = ipwdat[ipwdat$time==2 & ipwdat$s1>=2 & ipwdat$Rim1==1,]
ipwdat2 = ipwdat22[,c("yl", "shiftage", "sex","edu","time","s1")]
logistfit = glm(Ri ~ yl +sex, family = binomial(link=logit),data = ipwdat22)
ipwdat22$pred_obs = predict(logistfit, newdata = ipwdat2, type="response")

##########

ipwdat44 = ipwdat[ipwdat$time==4 &  ipwdat$s1>=4 & ipwdat$Rim1==1,]
ipwdat4 = ipwdat44[,c("yl", "lyl","llyl","lllyl", "shiftage", "sex","edu","time","s1")]
logistfit = glm(Ri ~ yl +sex, family = binomial(link=logit),data = ipwdat44)
ipwdat44$pred_obs = predict(logistfit, newdata = ipwdat4, type="response")

###############

ipwdat66 = ipwdat[ipwdat$time==6 &  ipwdat$s1>=6 & ipwdat$Rim1==1,]
ipwdat6 = ipwdat66[,c("yl", "lyl","llyl","lllyl", "shiftage", "sex","edu","time","s1")]
logistfit = glm(Ri ~ yl + sex, family = binomial(link=logit),data = ipwdat66)
ipwdat66$pred_obs = predict(logistfit, newdata = ipwdat6, type="response")

#################

ipwdat88 = ipwdat[ipwdat$time==8  &  ipwdat$s1>=8 & ipwdat$Rim1==1,]
ipwdat8 = ipwdat88[,c("yl", "lyl","llyl","lllyl", "shiftage", "sex","edu","time","s1")]
logistfit = glm(Ri ~ yl + sex, family = binomial(link=logit),data = ipwdat88)
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

ripwdat = ripwdat[order(ripwdat$Case),];ripwdat$pi=1/ripwdat$pi;
ipw.fit = geeglm(pek ~  factor(time) + factor(time)*sex, data = ripwdat, weights = pi, id=Case)

###IPW-IEE-D
ipwdat = dat[,c("Case", "time", "shiftage", "sex", "edu", "Ri","pek","s1","Rim1","yl", "lyl","llyl","lllyl","ager")]
ipwdat$pred_obs=NA;ipwdat$pi=NA;

ipwdat22 = ipwdat[ipwdat$time==2 & ipwdat$Rim1==1,]
ipwdat2 = ipwdat22[,c("yl", "shiftage", "sex","edu","time","s1")]
logistfit = glm(Ri ~ yl +sex+s1, family = binomial(link=logit),data = ipwdat22)
ipwdat22$pred_obs = predict(logistfit, newdata = ipwdat2, type="response")

##########

ipwdat44 = ipwdat[ipwdat$time==4 & ipwdat$Rim1==1,]
ipwdat4 = ipwdat44[,c("yl", "lyl","llyl","lllyl", "shiftage", "sex","edu","time","s1")]
logistfit = glm(Ri ~ yl + sex+s1, family = binomial(link=logit),data = ipwdat44)
ipwdat44$pred_obs = predict(logistfit, newdata = ipwdat4, type="response")

###############

ipwdat66 = ipwdat[ipwdat$time==6 & ipwdat$Rim1==1,]
ipwdat6 = ipwdat66[,c("yl", "lyl","llyl","lllyl", "shiftage", "sex","edu","time","s1")]
logistfit = glm(Ri ~ yl+sex+s1, family = binomial(link=logit),data = ipwdat66)
ipwdat66$pred_obs = predict(logistfit, newdata = ipwdat6, type="response")

#################

ipwdat88 = ipwdat[ipwdat$time==8 & ipwdat$Rim1==1,]
ipwdat8 = ipwdat88[,c("yl", "lyl","llyl","lllyl", "shiftage", "sex","edu","time","s1")]
logistfit = glm(Ri ~ yl + sex, family = binomial(link=logit),data = ipwdat88)
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
ripwdatd = rbind(ipwdat0,ipwdat2,ipwdat4,ipwdat6,ipwdat8,ipwdat[ipwdat$time==2 & ipwdat$Rim1==0,],ipwdat[ipwdat$time==4 & ipwdat$Rim1==0,],ipwdat[ipwdat$time==6 & ipwdat$Rim1==0,],ipwdat[ipwdat$time==8 & ipwdat$Rim1==0,])

ripwdatd = ripwdatd[order(ripwdatd$Case),]; ripwdatd$pi = 1/ripwdatd$pi;
ipw.fitd = geeglm(pek ~ factor(time)+factor(time)*sex, data = ripwdatd, weights = pi, id=Case)



## Augmented IPW-D ###
#### Imputation models (10)

midat12 = dat[dat$time==2 & dat$Rim1==1,]
mifit1 = lm(pek  ~ yl + sex + s1, data = midat12)
mi12 = midat12[,c("yl", "shiftage", "sex","edu","time","s1")]
midat12$pred_obs = predict(mifit1, newdata = mi12, type="response")
midat12$pred_obs = ifelse(!is.na(midat12$pek), midat12$pek, midat12$pred_obs)
midat12$pek = midat12$pred_obs; 

midat23 = dat[dat$time==4 & dat$Rim1==1,]
mifit2 = lm(pek  ~ yl + sex +  s1, data = midat23)
mi23 = midat23[,c("yl", "lyl","shiftage", "sex","edu","time","s1")]
midat23$pred_obs = predict(mifit2, newdata = mi23, type="response")
midat23$pred_obs = ifelse(!is.na(midat23$pek), midat23$pek, midat23$pred_obs)

midat34 = dat[dat$time==6 & dat$Rim1==1,]
mifit3 = lm(pek  ~ yl + sex + s1, data = midat34)
mi34 = midat34[,c("yl", "lyl","llyl", "shiftage", "sex","edu","time","s1")]
midat34$pred_obs = predict(mifit3, newdata = mi34, type="response")
midat34$pred_obs = ifelse(!is.na(midat34$pek), midat34$pek, midat34$pred_obs)

midat45 = dat[dat$time==8 & dat$Rim1==1,]
mifit4 = lm(pek  ~ yl + sex, data = midat45)
mi45 = midat45[,c("yl", "lyl","llyl","lllyl", "shiftage", "sex","edu","time","s1")]
midat45$pred_obs = predict(mifit4, newdata = mi45, type="response")
midat45$pred_obs = ifelse(!is.na(midat45$pek), midat45$pek, midat45$pred_obs)

###########################Fitted values are not covariates but outcomes
tmpmidat23 = dat[dat$time==4 & is.na(dat$yl),]
midat23$pek = midat23$pred_obs; midat23$pred_obs=NULL; midat23 = rbind(tmpmidat23,midat23);
midat23$Ri = ifelse(!is.na(midat23$pek),1,0)
mifit5 = lm(pek  ~  lyl + sex +s1 , data = midat23)
mi23 = midat23[,c("lyl", "shiftage", "sex","edu","time","s1")]
midat23$pred_obs = predict(mifit5, newdata = mi23, type="response")
midat23$pred_obs = ifelse(!is.na(midat23$pek), midat23$pek, midat23$pred_obs)
midat23$pek = midat23$pred_obs; 

tmpmidat34 = dat[dat$time==6 & is.na(dat$yl) & !is.na(dat$lyl),]
midat34$pek = midat34$pred_obs; midat34$pred_obs=NULL; midat34 = rbind(tmpmidat34,midat34);
midat34$Ri = ifelse(!is.na(midat34$pek),1,0)
mifit6 = lm(pek  ~ lyl + sex +s1, data = midat34)
mi34 = midat34[,c("lyl","llyl", "shiftage", "sex","edu","time","s1")]
midat34$pred_obs = predict(mifit6, newdata = mi34, type="response")
midat34$pred_obs = ifelse(!is.na(midat34$pek), midat34$pek, midat34$pred_obs)

tmpmidat45 = dat[dat$time==8 & is.na(dat$yl) & !is.na(dat$lyl),]
midat45$pek = midat45$pred_obs; midat45$pred_obs=NULL; midat45 = rbind(tmpmidat45,midat45);
midat45$Ri = ifelse(!is.na(midat45$pek),1,0)
mifit7 = lm(pek  ~ lyl + sex, data = midat45)
mi45 = midat45[,c("lyl","llyl","lllyl" , "shiftage", "sex","edu","time","s1")]
midat45$pred_obs = predict(mifit7, newdata = mi45, type="response")
midat45$pred_obs = ifelse(!is.na(midat45$pek), midat45$pek, midat45$pred_obs)

######
tmpmidat34 = dat[dat$time==6 & is.na(dat$yl) & is.na(dat$lyl),]
midat34$pek = midat34$pred_obs; midat34$pred_obs=NULL; midat34 = rbind(tmpmidat34,midat34);
midat34$Ri = ifelse(!is.na(midat34$pek),1,0)
mifit8 = lm(pek  ~ llyl + sex +s1, data = midat34)
mi34 = midat34[,c("llyl", "shiftage", "sex","edu","time","s1")]
midat34$pred_obs = predict(mifit8, newdata = mi34, type="response")
midat34$pred_obs = ifelse(!is.na(midat34$pek), midat34$pek, midat34$pred_obs)
midat34$pek = midat34$pred_obs; 

tmpmidat45 = dat[dat$time==8 & is.na(dat$yl) & is.na(dat$lyl) & !is.na(dat$llyl),]
midat45$pek = midat45$pred_obs; midat45$pred_obs=NULL; midat45 = rbind(tmpmidat45,midat45);
midat45$Ri = ifelse(!is.na(midat45$pek),1,0)
mifit9 = lm(pek  ~  llyl + sex, data = midat45)
mi45 = midat45[,c("llyl","lllyl", "shiftage", "sex","edu","time","s1")]
midat45$pred_obs = predict(mifit9, newdata = mi45, type="response")
midat45$pred_obs = ifelse(!is.na(midat45$pek), midat45$pek, midat45$pred_obs)

######
tmpmidat45 = dat[dat$time==8 & is.na(dat$yl) & is.na(dat$llyl),]
midat45$pek = midat45$pred_obs; midat45$pred_obs=NULL; midat45 = rbind(tmpmidat45,midat45);
midat45$Ri = ifelse(!is.na(midat45$pek),1,0)
mifit10 = lm(pek  ~ lllyl + sex, data = midat45)
mi45 = midat45[,c("lllyl", "shiftage", "sex","edu","time","s1")]
midat45$pred_obs = predict(mifit10, newdata = mi45, type="response")
midat45$pred_obs = ifelse(!is.na(midat45$pek), midat45$pek, midat45$pred_obs)
midat45$pek = midat45$pred_obs; 

#####
midat = rbind(midat12, midat23, midat34,midat45); midat$pred_obs = NULL;
midat = rbind(dat[dat$time==0,],midat)
midatsortd = midat[order(midat$Case),]
midatsortd = midatsortd[midatsortd$time<=midatsortd$s1,]

########## complete data set
midatsortd$pi=1/ripwdatd$pi
midatsortd$Rim1 = midatsortd$yl = midatsortd$lyl = midatsortd$llyl = midatsortd$lllyl = NULL;
midatsortd$Ri = ripwdatd$Ri; midatsortd$digreal = ripwdatd$pek;
##NR##

beta = c(18,0,0,0,0,0,0,0,0,0);
diff=10; crit=0.00001
while(diff>=crit)
{
  phi = phic = 1;
  U = 0; dU=0;
  auniqueid = unique(midatsortd$Case);
  for(i in auniqueid)
  {
    tmp5 = midatsortd[midatsortd$Case==i,];       
    n = nrow(tmp5);
    tmpdiag = diag(5); tmpdiag = tmpdiag[1:n,!tmpdiag[,1]]; tmpdiag2 = tmp5$sex[1]*tmpdiag;
    x = matrix(c(rep(1, n),as.vector(tmpdiag), tmp5$sex,as.vector(tmpdiag2)),nrow=n) #nix4
    last = tmp5[n,]; last$time=10; last2 = do.call("rbind",replicate(5-n, last, simplify = FALSE));tmp5 = rbind(tmp5,last2);
    
    ##j=1
    xtmp = as.data.frame(matrix(c(rep(tmp5$digreal[1],4), tmp5$sex[1], tmp5$edu[1],tmp5$s1[1]),nrow=1))
    colnames(xtmp) = c("yl","lyl","llyl","lllyl","sex","edu","s1")
    
    if(!is.na(tmp5$digreal[1])) {W1=c(tmp5$digreal[1],predict(mifit1, newdata = xtmp, type="response"), predict(mifit5, newdata = xtmp, type="response"),predict(mifit8, newdata = xtmp, type="response"),predict(mifit10, newdata = xtmp, type="response"))[1:n]}else{W1=NULL}
    
    #j=2
    xtmp2l = data.frame(cbind(tmp5$digreal[2],tmp5$sex[1],tmp5$s1[1]));
    xtmp2 = data.frame(cbind(tmp5$digreal[2],tmp5$sex[1]));
    
    if(!is.na(tmp5$digreal[2])){W2=c(tmp5$digreal[1], tmp5$digreal[2],mifit2$coef[1]+sum(mifit2$coef[2:length(mifit2$coef)]*xtmp2l),mifit6$coef[1]+sum(mifit6$coef[2:length(mifit6$coef)]*xtmp2l),mifit9$coef[1]+sum(mifit9$coef[2:length(mifit9$coef)]*xtmp2))[1:n]}else{W2=NULL}
    
    #j=3
    xtmp3l = data.frame(cbind(tmp5$digreal[3],tmp5$sex[1],tmp5$s1[1]));
    xtmp3 = data.frame(cbind(tmp5$digreal[3],tmp5$sex[1]));
    
    if(!is.na(tmp5$digreal[3])){W3=c(tmp5$digreal[1], tmp5$digreal[2],tmp5$digreal[3],mifit3$coef[1]+sum(mifit3$coef[2:length(mifit3$coef)]*xtmp3l),mifit7$coef[1]+sum(mifit7$coef[2:length(mifit7$coef)]*xtmp3))[1:n]}else{W3=NULL}
    
    #j=4
    xtmp4 = data.frame(cbind(tmp5$digreal[4],tmp5$sex[1]));
    
    if(!is.na(tmp5$digreal[3])){W4=c(tmp5$digreal[1], tmp5$digreal[2],tmp5$digreal[3],tmp5$digreal[4],mifit4$coef[1]+sum(mifit4$coef[2:length(mifit4$coef)]*xtmp4))[1:n]}else{W4=NULL}
    
    fitted = x %*% t(t(beta)); tmp5 = tmp5[tmp5$s1>=tmp5$time,];
    
    ### For the left hand side of AIPW (complete side)
    var = diag(phi,n); varc = diag(phic,n)
    if(!is.na(tmp5$pi[n]) & !is.na(tmp5$digreal[n])) {uleft=(1/tmp5$pi[n])*(t(x) %*% var %*% (tmp5$digreal - x%*%t(t(beta)) ) )} else{uleft=matrix(c(rep(0,ncol(x))),ncol=1)}
    if(!is.na(tmp5$pi[n]) & !is.na(tmp5$digreal[n])) {duleft = (1/tmp5$pi[n])*(t(x) %*% varc %*% x)} else{duleft=diag(0,ncol(x))}
    
    ### For the right hand side
    if(!is.null(W1) & !is.na(tmp5$pi[2])) {c1=as.matrix((tmp5$Ri[1]/(tmp5$pi[1])-(tmp5$Ri[2])/(tmp5$pi[2]))*(W1-fitted)); p1=(tmp5$Ri[1]/(tmp5$pi[1])-(tmp5$Ri[2])/(tmp5$pi[2]))} else{c1=matrix(c(rep(0,n)),ncol=1);p1=0}
    if(!is.null(W2)& !is.na(tmp5$pi[3])) {c2=as.matrix((tmp5$Ri[2]/(tmp5$pi[2])-(tmp5$Ri[3])/(tmp5$pi[3]))*(W2-fitted));p2=(tmp5$Ri[2]/(tmp5$pi[2])-(tmp5$Ri[3])/(tmp5$pi[3]))} else{c2=matrix(c(rep(0,n)),ncol=1); p2=0;}
    if(!is.null(W3)& !is.na(tmp5$pi[4])) {c3=as.matrix((tmp5$Ri[3]/(tmp5$pi[3])-(tmp5$Ri[4])/(tmp5$pi[4]))*(W3-fitted)); p3=(tmp5$Ri[3]/(tmp5$pi[3])-(tmp5$Ri[4])/(tmp5$pi[4]))} else{c3=matrix(c(rep(0,n)),ncol=1); p3=0;}
    if(!is.null(W4)& !is.na(tmp5$pi[5])) {c4=as.matrix((tmp5$Ri[4]/(tmp5$pi[4])-(tmp5$Ri[5])/(tmp5$pi[5]))*(W4-fitted));p4=(tmp5$Ri[4]/(tmp5$pi[4])-(tmp5$Ri[5])/(tmp5$pi[5]))} else{c4=matrix(c(rep(0,n)),ncol=1);p4=0;}
    uright = t(x) %*% var %*% (c1+c2+c3+c4)
    duright = t(x) %*% var %*% x *(p1+p2+p3+p4)
    
    if(tmp5$Ri[n]==1) {Ubeta = uleft+uright; dUbeta = duleft+duright} else{Ubeta = uright; dUbeta = duright;}
    
    U = U + Ubeta
    dU = dU + dUbeta
  }
  diff = max(solve(-dU) %*% U)
  beta = beta - solve(-dU) %*% U
  }
mybetad = beta

###AIPWIEE (not conditioning on D)
dat = tmpdata
#dat = dat[dat$s1>=dat$time,]
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
mifit1 = lm(pek  ~ yl + sex, data = midat12)
mi12 = midat12[,c("yl", "shiftage", "sex","edu","time","s1")]
midat12$pred_obs = predict(mifit1, newdata = mi12, type="response")
midat12$pred_obs = ifelse(!is.na(midat12$pek), midat12$pek, midat12$pred_obs)
midat12$pek = midat12$pred_obs; 

midat23 = dat[dat$time==4 & dat$Rim1==1,]
mifit2 = lm(pek  ~ yl + sex, data = midat23)
mi23 = midat23[,c("yl", "lyl","shiftage", "sex","edu","time","s1")]
midat23$pred_obs = predict(mifit2, newdata = mi23, type="response")
midat23$pred_obs = ifelse(!is.na(midat23$pek), midat23$pek, midat23$pred_obs)

midat34 = dat[dat$time==6 & dat$Rim1==1,]
mifit3 = lm(pek  ~ yl + sex, data = midat34)
mi34 = midat34[,c("yl", "lyl","llyl", "shiftage", "sex","edu","time","s1")]
midat34$pred_obs = predict(mifit3, newdata = mi34, type="response")
midat34$pred_obs = ifelse(!is.na(midat34$pek), midat34$pek, midat34$pred_obs)

midat45 = dat[dat$time==8 & dat$Rim1==1,]
mifit4 = lm(pek  ~ yl + sex, data = midat45)
mi45 = midat45[,c("yl", "lyl","llyl","lllyl", "shiftage", "sex","edu","time","s1")]
midat45$pred_obs = predict(mifit4, newdata = mi45, type="response")
midat45$pred_obs = ifelse(!is.na(midat45$pek), midat45$pek, midat45$pred_obs)

###########################Fitted values are not covariates but outcomes
tmpmidat23 = dat[dat$time==4 & is.na(dat$yl),]
midat23$pek = midat23$pred_obs; midat23$pred_obs=NULL; midat23 = rbind(tmpmidat23,midat23);
midat23$Ri = ifelse(!is.na(midat23$pek),1,0)
mifit5 = lm(pek  ~  lyl + sex, data = midat23)
mi23 = midat23[,c("lyl", "shiftage", "sex","edu","time","s1")]
midat23$pred_obs = predict(mifit5, newdata = mi23, type="response")
midat23$pred_obs = ifelse(!is.na(midat23$pek), midat23$pek, midat23$pred_obs)
midat23$pek = midat23$pred_obs; 

tmpmidat34 = dat[dat$time==6 & is.na(dat$yl) & !is.na(dat$lyl),]
midat34$pek = midat34$pred_obs; midat34$pred_obs=NULL; midat34 = rbind(tmpmidat34,midat34);
midat34$Ri = ifelse(!is.na(midat34$pek),1,0)
mifit6 = lm(pek  ~ lyl + sex, data = midat34)
mi34 = midat34[,c("lyl","llyl", "shiftage", "sex","edu","time","s1")]
midat34$pred_obs = predict(mifit6, newdata = mi34, type="response")
midat34$pred_obs = ifelse(!is.na(midat34$pek), midat34$pek, midat34$pred_obs)

tmpmidat45 = dat[dat$time==8 & is.na(dat$yl) & !is.na(dat$lyl),]
midat45$pek = midat45$pred_obs; midat45$pred_obs=NULL; midat45 = rbind(tmpmidat45,midat45);
midat45$Ri = ifelse(!is.na(midat45$pek),1,0)
mifit7 = lm(pek  ~ lyl + sex, data = midat45)
mi45 = midat45[,c("lyl","llyl","lllyl" , "shiftage", "sex","edu","time","s1")]
midat45$pred_obs = predict(mifit7, newdata = mi45, type="response")
midat45$pred_obs = ifelse(!is.na(midat45$pek), midat45$pek, midat45$pred_obs)

######
tmpmidat34 = dat[dat$time==6 & is.na(dat$yl) & is.na(dat$lyl),]
midat34$pek = midat34$pred_obs; midat34$pred_obs=NULL; midat34 = rbind(tmpmidat34,midat34);
midat34$Ri = ifelse(!is.na(midat34$pek),1,0)
mifit8 = lm(pek  ~ llyl + sex, data = midat34)
mi34 = midat34[,c("llyl", "shiftage", "sex","edu","time","s1")]
midat34$pred_obs = predict(mifit8, newdata = mi34, type="response")
midat34$pred_obs = ifelse(!is.na(midat34$pek), midat34$pek, midat34$pred_obs)
midat34$pek = midat34$pred_obs; 

tmpmidat45 = dat[dat$time==8 & is.na(dat$yl) & is.na(dat$lyl) & !is.na(dat$llyl),]
midat45$pek = midat45$pred_obs; midat45$pred_obs=NULL; midat45 = rbind(tmpmidat45,midat45);
midat45$Ri = ifelse(!is.na(midat45$pek),1,0)
mifit9 = lm(pek  ~  llyl + sex, data = midat45)
mi45 = midat45[,c("llyl","lllyl", "shiftage", "sex","edu","time","s1")]
midat45$pred_obs = predict(mifit9, newdata = mi45, type="response")
midat45$pred_obs = ifelse(!is.na(midat45$pek), midat45$pek, midat45$pred_obs)

######
tmpmidat45 = dat[dat$time==8 & is.na(dat$yl) & is.na(dat$llyl),]
midat45$pek = midat45$pred_obs; midat45$pred_obs=NULL; midat45 = rbind(tmpmidat45,midat45);
midat45$Ri = ifelse(!is.na(midat45$pek),1,0)
mifit10 = lm(pek  ~ lllyl + sex, data = midat45)
mi45 = midat45[,c("lllyl", "shiftage", "sex","edu","time","s1")]
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
beta = c(18,0,0,0,0,0,0,0,0,0);
diff=10;crit=0.00001;
while(diff>=crit)
{
  phi = phic = 1;
  U = 0; dU=0;
  auniqueid = unique(midatsort$Case);
  for(i in auniqueid)
  {
    tmp5 = midatsort[midatsort$Case==i,];       
    n = nrow(tmp5);
    tmpdiag = diag(5); tmpdiag = tmpdiag[1:n,!tmpdiag[,1]]; tmpdiag2 = tmp5$sex[1]*tmpdiag;
    x = matrix(c(rep(1, n),as.vector(tmpdiag), tmp5$sex,as.vector(tmpdiag2)),nrow=n) #nix4
    last = tmp5[n,]; last$time=10; last2 = do.call("rbind",replicate(5-n, last, simplify = FALSE));tmp5 = rbind(tmp5,last2);
    
    ##j=1
    xtmp = as.data.frame(matrix(c(rep(tmp5$digreal[1],4), tmp5$sex[1], tmp5$edu[1]),nrow=1))
    colnames(xtmp) = c("yl","lyl","llyl","lllyl","sex","edu")
    
    if(!is.na(tmp5$digreal[1])) {W1=c(tmp5$digreal[1],predict(mifit1, newdata = xtmp, type="response"), predict(mifit5, newdata = xtmp, type="response"),predict(mifit8, newdata = xtmp, type="response"),predict(mifit10, newdata = xtmp, type="response"))[1:n]}else{W1=NULL}
    
    #j=2
    xtmp2 = data.frame(cbind(tmp5$digreal[2],tmp5$sex[1]));
    
    if(!is.na(tmp5$digreal[2])){W2=c(tmp5$digreal[1], tmp5$digreal[2],mifit2$coef[1]+sum(mifit2$coef[2:length(mifit2$coef)]*xtmp2),mifit6$coef[1]+sum(mifit6$coef[2:length(mifit6$coef)]*xtmp2),mifit9$coef[1]+sum(mifit9$coef[2:length(mifit9$coef)]*xtmp2))[1:n]}else{W2=NULL}
    
    #j=3
    xtmp3 = data.frame(cbind(tmp5$digreal[3],tmp5$sex[1]));
    
    if(!is.na(tmp5$digreal[3])){W3=c(tmp5$digreal[1], tmp5$digreal[2],tmp5$digreal[3],mifit3$coef[1]+sum(mifit3$coef[2:length(mifit3$coef)]*xtmp3),mifit7$coef[1]+sum(mifit7$coef[2:length(mifit7$coef)]*xtmp3))[1:n]}else{W3=NULL}
    
    #j=4
    xtmp4 = data.frame(cbind(tmp5$digreal[4],tmp5$sex[1]));
    
    if(!is.na(tmp5$digreal[4])){W4=c(tmp5$digreal[1], tmp5$digreal[2],tmp5$digreal[3],tmp5$digreal[4],mifit4$coef[1]+sum(mifit4$coef[2:length(mifit4$coef)]*xtmp4))[1:n]}else{W4=NULL}
    
    fitted = x %*% t(t(beta)); tmp5 = tmp5[tmp5$s1>=tmp5$time,];
    
    ### For the left hand side of AIPW (complete side)
    var = diag(phi,n); varc = diag(phic,n)
    if(!is.na(tmp5$pi[n]) & !is.na(tmp5$digreal[n])) {uleft=(1/tmp5$pi[n])*(t(x) %*% varc %*% (tmp5$digreal - x%*%t(t(beta)) ) )} else{uleft=matrix(c(rep(0,ncol(x))),ncol=1)}
    if(!is.na(tmp5$pi[n]) & !is.na(tmp5$digreal[n])) {duleft = (1/tmp5$pi[n])*(t(x) %*% varc %*% x)} else{duleft=diag(0,ncol(x))}
    
    
    ### For the right hand side
    if(!is.null(W1) & !is.na(tmp5$pi[2])) {c1=as.matrix((tmp5$Ri[1]/(tmp5$pi[1])-(tmp5$Ri[2])/(tmp5$pi[2]))*(W1-fitted)); p1=(tmp5$Ri[1]/(tmp5$pi[1])-(tmp5$Ri[2])/(tmp5$pi[2]))} else{c1=matrix(c(rep(0,n)),ncol=1);p1=0}
    if(!is.null(W2)& !is.na(tmp5$pi[3])) {c2=as.matrix((tmp5$Ri[2]/(tmp5$pi[2])-(tmp5$Ri[3])/(tmp5$pi[3]))*(W2-fitted));p2=(tmp5$Ri[2]/(tmp5$pi[2])-(tmp5$Ri[3])/(tmp5$pi[3]))} else{c2=matrix(c(rep(0,n)),ncol=1); p2=0;}
    if(!is.null(W3)& !is.na(tmp5$pi[4])) {c3=as.matrix((tmp5$Ri[3]/(tmp5$pi[3])-(tmp5$Ri[4])/(tmp5$pi[4]))*(W3-fitted)); p3=(tmp5$Ri[3]/(tmp5$pi[3])-(tmp5$Ri[4])/(tmp5$pi[4]))} else{c3=matrix(c(rep(0,n)),ncol=1); p3=0;}
    if(!is.null(W4)& !is.na(tmp5$pi[5])) {c4=as.matrix((tmp5$Ri[4]/(tmp5$pi[4])-(tmp5$Ri[5])/(tmp5$pi[5]))*(W4-fitted));p4=(tmp5$Ri[4]/(tmp5$pi[4])-(tmp5$Ri[5])/(tmp5$pi[5]))} else{c4=matrix(c(rep(0,n)),ncol=1);p4=0;}
    uright = t(x) %*% var %*% (c1+c2+c3+c4)
    duright = t(x) %*% var %*% x *(p1+p2+p3+p4)
    
    if(tmp5$Ri[n]==1) {Ubeta = uleft+uright; dUbeta = duleft+duright} else{Ubeta = uright; dUbeta = duright;}
    
    U = U + Ubeta
    dU = dU + dUbeta
  }
  diff = max(solve(-dU) %*% U)
  beta = beta - solve(-dU) %*% U
}
mybeta=beta

###MI
wd = reshape(tmpdata, idvar = "Case", v.names=c("pek"),timevar="time", direction="wide")
midatd = wd[, c("sex","s1","pek.0","pek.2","pek.4","pek.6","pek.8")]

imdatad = as.matrix(midatd)

nsim = 30
PE = matrix(0,nsim,10)
for(i in 1:nsim)
{
  s = prelim.norm(imdatad)
  thetahat = em.norm(s,showits=F)
  rngseed(i)
  ximp = imp.norm(s, thetahat, imdatad)
  tmp = as.data.frame(ximp)
  
  long = reshape(tmp, varying=list(c(3:7)), v.names = c("pek"), new.row.names=1:1000000, times=c(0,2,4,6,8), direction="long")
  longsort = long[order(long$id),];
  midatd = longsort[longsort$time<= longsort$s1,]
  
  fitmi = geeglm(pek ~  factor(time)+factor(time)*sex, data = midatd, id=id,na.action=na.omit)
  PE[i,]=summary(fitmi)$coef[,1]
}

COEFd = apply(PE, 2, mean)



###MI (without d)
w = reshape(tmpdata, idvar = "Case", v.names=c("pek"),timevar="time", direction="wide")
midat = w[, c("sex","pek.0","pek.2","pek.4","pek.6","pek.8")]

imdata = as.matrix(midat)
rngseed(112789)

nsim = 30
PE = matrix(0,nsim,10)
for(i in 1:nsim)
{
  s = prelim.norm(imdata)
  thetahat = em.norm(s,showits=F)
  #theta1 = da.norm(s,thetahat,steps=5)
  ximp = imp.norm(s, thetahat, imdata)
  tmp = as.data.frame(ximp); tmp$s1=w$s1;
  
  long = reshape(tmp, varying=list(c(2:6)), v.names = c("pek"), new.row.names=1:1000000, times=c(0,2,4,6,8), direction="long")
  longsort = long[order(long$id),];
  midat = longsort[longsort$time<= longsort$s1,]
  
  fitmi = geeglm(pek ~  factor(time)+factor(time)*sex, data = midat, id=id,na.action=na.omit)
  PE[i,]=summary(fitmi)$coef[,1]
}

COEF = apply(PE, 2, mean)


###LI-D uMVN
wd$alive1=1;wd$alive2=ifelse(wd$s1>=2,1,0);wd$alive3=ifelse(wd$s1>=4,1,0);wd$alive4=ifelse(wd$s1>=6,1,0);wd$alive5 = ifelse(wd$s1>=8,1,0);
wdm = as.matrix(wd)
tmpli = linearincrements(y=wdm[,7:11],x=wdm[,c(3,5)],alive=wdm[,12:16],method="uMVN",returnimps=T)
tmpliout = tmpli$yimp
newdatlid = cbind(wd[,1:6],tmpliout)
longlid = reshape(newdatlid, varying=list(c(7:11)), v.names = c("pek"), new.row.names=1:1000000, times=c(0,2,4,6,8), direction="long")
lidatd = longlid[order(longlid$id),];
lidatd = lidatd[lidatd$time<= lidatd$s1,]
fitlid = geeglm(pek ~  factor(time)+factor(time)*sex, id=Case, data = lidatd)
PElid=summary(fitlid)$coef[,1];

###LI uMVN
w$alive1=1;w$alive2=ifelse(w$s1>=2,1,0);w$alive3=ifelse(w$s1>=4,1,0);w$alive4=ifelse(w$s1>=6,1,0);w$alive5 = ifelse(w$s1>=8,1,0);
wm = as.matrix(wd)
tmpli = linearincrements(y=wm[,7:11],x=wm[,c(3)],alive=wm[,12:16],method="uMVN",returnimps=T)
tmpliout = tmpli$yimp
newdatli = cbind(w[,1:6],tmpliout)
longli = reshape(newdatli, varying=list(c(7:11)), v.names = c("pek"), new.row.names=1:1000000, times=c(0,2,4,6,8), direction="long")
lidat = longli[order(longli$id),];
lidat = lidat[lidat$time<= lidat$s1,]
fitli = geeglm(pek ~ factor(time)+factor(time)*sex, id=Case, data = lidat)
PEli=summary(fitli)$coef[,1];


outtmp = c(summary(fit)$coefficient[,1],
           summary(ipw.fit)$coefficient[,1],
           summary(ipw.fitd)$coefficient[,1],
           mybeta,
           mybetad,
           PEli,
           PElid,
           COEF,
           COEFd)
outtmp
}
test = foreach(m=1:1000) %dopar% myfunc(m)
test2 = do.call("rbind", test)

write.csv(test2,"miss_methods.csv")

