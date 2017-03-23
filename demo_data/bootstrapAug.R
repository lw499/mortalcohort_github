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
  mytmpdat = dat

dat = mytmpdat[mytmpdat$s1>=mytmpdat$time,] 
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

###AIPWIEE (not conditioning on D)

dat = mytmpdat
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

diff=10;crit=0.00001;r=1;
while(diff>=crit)
{
  phi = phic = 1;
  U = 0; dU=0;
  auniqueid = unique(midatsort$Case);
  for(i in auniqueid)
  {
    tmp5 = midatsort[midatsort$Case==i,];       
    n = nrow(tmp5);
    if(tmp5$Smoke[1]==0){smokestat1 = rep(0,n)} else if(tmp5$Smoke[1]==1){smokestat1 = rep(1,n);}
    x = matrix(c(rep(1, n),tmp5$time, tmp5$time^2, tmp5$Female, tmp5$shiftage, tmp5$Educyrs, smokestat1, tmp5$time*tmp5$Female,tmp5$time*tmp5$shiftage,tmp5$time*tmp5$Educyrs, tmp5$time*smokestat1),nrow=n) #nix4
    
    last = tmp5[n,]; last$time=10; last2 = do.call("rbind",replicate(5-n, last, simplify = FALSE));tmp5 = rbind(tmp5,last2);
    
    ##j=1
    xtmp = as.data.frame(matrix(c(rep(tmp5$digreal[1],4), tmp5$shiftage[1], tmp5$Educyrs[1], tmp5$Female[1],tmp5$mmse[1],tmp5$Smoke[1],tmp5$iadl1[1],tmp5$HlthPrv1[1]),nrow=1))
    colnames(xtmp) = c("yl","lyl","llyl","lllyl","shiftage","Educyrs","Female","mmse1","Smoke","iadl1","HlthPrv1")
    
    if(!is.na(tmp5$digreal[1])) {W1=c(tmp5$digreal[1],predict(mifit1, newdata = xtmp, type="response"), predict(mifit5, newdata = xtmp, type="response"),predict(mifit8, newdata = xtmp, type="response"),predict(mifit10, newdata = xtmp, type="response"))[1:n]}else{W1=NULL}
    
    #j=2
    xtmp2 = data.frame(cbind(tmp5$digreal[2],tmp5$digreal[1],tmp5$shiftage[1],tmp5$Educyrs[1],tmp5$Female[1],tmp5$iadl1[1],tmp5$HlthPrv1[1],tmp5$mmse1[1],tmp5$Smoke[1]));
    
    if(!is.na(tmp5$digreal[2])){W2=c(tmp5$digreal[1], tmp5$digreal[2],mifit2$coef[1]+sum(mifit2$coef[2:length(mifit2$coef)]*xtmp2),mifit6$coef[1]+sum(mifit6$coef[2:length(mifit6$coef)]*xtmp2),mifit9$coef[1]+sum(mifit9$coef[2:length(mifit9$coef)]*xtmp2))[1:n]}else{W2=NULL}
    
    #j=3
    xtmp3 = data.frame(cbind(tmp5$digreal[3],tmp5$digreal[2],tmp5$digreal[1],tmp5$shiftage[1],tmp5$Educyrs[1],tmp5$Female[1],tmp5$iadl1[1],tmp5$HlthPrv1[1],tmp5$mmse1[1],tmp5$Smoke[1]));
    
    if(!is.na(tmp5$digreal[3])){W3=c(tmp5$digreal[1], tmp5$digreal[2],tmp5$digreal[3],mifit3$coef[1]+sum(mifit3$coef[2:length(mifit3$coef)]*xtmp3),mifit7$coef[1]+sum(mifit7$coef[2:length(mifit7$coef)]*xtmp3))[1:n]}else{W3=NULL}
    
    #j=4
    xtmp4 = data.frame(cbind(tmp5$digreal[4],tmp5$digreal[3],tmp5$digreal[2],tmp5$digreal[1],tmp5$shiftage[1],tmp5$Educyrs[1],tmp5$Female[1],tmp5$iadl1[1],tmp5$HlthPrv1[1],tmp5$mmse1[1],tmp5$Smoke[1]));
    
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
  #r=r+1
  #cat(r, "\n")
  }
return(c(beta))
}
test = foreach(m=1:1000) %dopar% myfunc(m)
test2 = do.call("rbind", test)

write.csv(test2,"listaug.csv")