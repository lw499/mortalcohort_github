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
  set.seed(10)
  seeds = floor(runif(1000)*10^8);
  
  EXPIT <- function(term) {
    return( exp(term)/(1+exp(term)) )
  }
  
  set.seed(seeds[m])
N=1000000
shiftagesd=2.613984 ; 
edusd = 2.336318;
shiftage = rnorm(N, 3.576063, shiftagesd); 
sex = rbinom(N,1,0.5); 
edu = rnorm(N,5.19222,edusd);   # Age, sex, Edu

y1 = rnorm(N, 17.63033,3.132368*0.75)

b01=5.29500; b11=0.73491; b21=-0.51502; e1=1.85;
b02=6.53720; b22=0.64448; b32=-0.71516; e2=1.818;
b03=4.18131; b33=0.74205; b43=-0.77404; e3=1.787;
b04=4.20122; b44=0.71654; b54=-0.02033; e4=1.971;

g01=1.502150*-0.75; g11=0.043*4.5; g21=0.35048;
g02=0.13903*-42.5; g22=0.09607*4.5; g32=0.40396;
g03=0.67365*-4.75;  g33=0.05374*4.5; g43=0.98483; 
g04=0.30571*-10.5; g44=0.05598*4.5; g54=0.62457;

y2 =  b01 + y1*b11 + sex*b21 + rnorm(N,0,e1)
y3 =  b02 + y2*b22 + sex*b32 + rnorm(N,0,e2)
y4 =  b03 + y3*b33 + sex*b43 + rnorm(N,0,e3)
y5 =  b04 + y4*b44 + sex*b54 + rnorm(N,0,e4)
edu = ifelse(edu<0, 0, edu)
edu = ifelse(edu<0, 0, edu)

Case = seq(1,N,1)
dat = as.data.frame(cbind(Case, edu, sex, shiftage, y1, y2, y3, y4, y5 ))

##Simulating missingness and survival probabilities##
dat$D2 = EXPIT(g01+g11*dat$y1+g21*sex)
dat$D3 = EXPIT(g02+g22*dat$y2+g32*sex)
dat$D4 = EXPIT(g03+g33*dat$y3+g43*sex)
dat$D5 = EXPIT(g04+g44*dat$y4+g54*sex)

### Di>=0 and Di>=2 ###
dat$s1=ifelse(runif(nrow(dat))<=dat$D2,999,0); dat0 = dat[dat$s1==0,]; dattmp = dat[dat$s1==999,];
dattmp$s1=ifelse(runif(nrow(dattmp))<=dattmp$D3, 999, 2); dat2 = dattmp[dattmp$s1==2,]; dattmp = dattmp[dattmp$s1==999,];
dattmp$s1=ifelse(runif(nrow(dattmp))<=dattmp$D4, 999, 4); dat4 = dattmp[dattmp$s1==4,]; dattmp = dattmp[dattmp$s1==999,];
dattmp$s1=ifelse(runif(nrow(dattmp))<=dattmp$D5, 999, 6); dat6 = dattmp[dattmp$s1==6,]; dat8 = dattmp[dattmp$s1==999,];


dat8$s1=8;
datfin = rbind(dat0,dat2,dat4,dat6,dat8)
datfin$D2=datfin$D3=datfin$D4=datfin$D5=NULL;

dat=datfin
#dat$ager = ifelse(dat$shiftage<=4,1,0);
#dat$sex = ifelse(dat$sex<0,0,1);

tmp1 = dat[dat$s1>=0 & dat$sex==0,]; m1=mean(tmp1$y1);
tmp2 = dat[dat$s1>=2 & dat$sex==0,]; m2=mean(tmp2$y2);
tmp3 = dat[dat$s1>=4 & dat$sex==0,]; m3=mean(tmp3$y3);
tmp4 = dat[dat$s1>=6 & dat$sex==0,]; m4=mean(tmp4$y4);
tmp5 = dat[dat$s1==8 & dat$sex==0,]; m5=mean(tmp5$y5);

tmp1 = dat[dat$s1>=0 & dat$sex==1,]; m11=mean(tmp1$y1);
tmp2 = dat[dat$s1>=2 & dat$sex==1,]; m21=mean(tmp2$y2);
tmp3 = dat[dat$s1>=4 & dat$sex==1,]; m31=mean(tmp3$y3);
tmp4 = dat[dat$s1>=6 & dat$sex==1,]; m41=mean(tmp4$y4);
tmp5 = dat[dat$s1==8 & dat$sex==1,]; m51=mean(tmp5$y5);

b0=m1; b1=m2-b0; b2=m3-b0; b3=m4-b0; b4=m5-b0;
b5 = m11-b0; b6=m21-m2-b5;
b7=m31-m3-b5; b8 = m41-m4-b5; b9=m51-m5-b5 
out = c(b0,b1,b2,b3,b4,b5,b6,b7,b8,b9)
out
}
test = foreach(m=1:1000) %dopar% myfunc(m)
test2 = do.call("rbind", test)



write.csv(test2,"true_mean.csv")
