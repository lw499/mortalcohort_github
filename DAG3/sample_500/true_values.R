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
t = c(0,2,4,6,8);
s = rnorm(N, 7.492475, 3.714435)
s1=ifelse(s<2,0,8);s1=ifelse(s>=2 & s<4, 2, s1);
s1=ifelse(s>=4 & s<6, 4, s1); s1=ifelse(s>=6 & s<8, 6, s1);

y1 = rnorm(N, 17.63033,3.132368*0.75)

#b01=4.88947; b11=0.73191; b21=-0.53041; b41=0.67732; b51=0.31168; b61=0.54708; e1=1.849;
#b02=5.73597; b22=0.64988; b32=-0.80995; b52=0.65343; b62=0.89493; e2=1.807;
#b03=4.14446; b33=0.74143; b43=-0.77977; b63=0.06189; e3=1.793; 
#b04=4.20122; b44=0.71654; b54=-0.02033; e4=1.971;

b01=5.05421; b11=0.73159; b21=-0.56623; b41=0.05197; e1=1.85;
b02=5.16183; b22=0.64555; b32=-0.83485; b52=0.19805; e2=1.804;
b03=3.95880; b33=0.74143; b43=-0.77977; b63=0.03094; e3=1.793; 
b04=4.20122; b44=0.71654; b54=-0.02033; e4=1.971;

y2 =  b01 + y1*b11 + sex*b21 + s1*b41 + rnorm(N,0,e1)
y3 =  b02 + y2*b22 + sex*b32 + s1*b52 + rnorm(N,0,e2)
y4 =  b03 + y3*b33 + sex*b43 + s1*b63 + rnorm(N,0,e3)
y5 =  b04 + y4*b44 + sex*b54 + rnorm(N,0,e4)

#y2 =  b01 + y1*b11 + sex*b21 + I(s1==4)*b41 + I(s1==6)*b51 + I(s1==8)*b61 + rnorm(N,0,e1)
#y3 =  b02 + y2*b22 + sex*b32 + I(s1==6)*b52 + I(s1==8)*b62 + rnorm(N,0,e2)
#y4 =  b03 + y3*b33 + sex*b43 + I(s1==8)*b63 + rnorm(N,0,e3)
#y5 =  b04 + y4*b44 + sex*b54 + rnorm(N,0,e4)
edu = ifelse(edu<0, 0, edu)
edu = ifelse(edu<0, 0, edu)

Case = seq(1,N,1)
dat = as.data.frame(cbind(Case, edu, sex, shiftage, y1, y2, y3, y4, y5,s))

### Di>=0 and Di>=2 ###
dat$s1=ifelse(dat$s<2,0,8);dat$s1=ifelse(dat$s>=2 & dat$s<4, 2, dat$s1);
dat$s1=ifelse(dat$s>=4 & dat$s<6, 4, dat$s1); dat$s1=ifelse(dat$s>=6 & dat$s<8, 6, dat$s1);

#tmp1 = dat[dat$s1>=0 & dat$sex==0,]; m1=mean(tmp1$y1);
#tmp2 = dat[dat$s1>=2 & dat$sex==0,]; m2=mean(tmp2$y2);
#tmp3 = dat[dat$s1>=4 & dat$sex==0,]; m3=mean(tmp3$y3);
#tmp4 = dat[dat$s1>=6 & dat$sex==0,]; m4=mean(tmp4$y4);
#tmp5 = dat[dat$s1>=8 & dat$sex==0,]; m5=mean(tmp5$y5);

#tmp1 = dat[dat$s1>=0 & dat$sex==1,]; m11=mean(tmp1$y1);
#tmp2 = dat[dat$s1>=2 & dat$sex==1,]; m21=mean(tmp2$y2);
#tmp3 = dat[dat$s1>=4 & dat$sex==1,]; m31=mean(tmp3$y3);
#tmp4 = dat[dat$s1>=6 & dat$sex==1,]; m41=mean(tmp4$y4);
#tmp5 = dat[dat$s1>=8 & dat$sex==1,]; m51=mean(tmp5$y5);

#b0=m1; b1=m2-b0; b2=m3-b0; b3=m4-b0; b4=m5-b0;
#b5 = m11-b0; b6=m21-m2-b5;
#b7=m31-m3-b5; b8 = m41-m4-b5; b9=m51-m5-b5 
#out = c(b0,b1,b2,b3,b4,b5,b6,b7,b8,b9)
tmp1 = dat[dat$s1>=0,]; m1=mean(tmp1$y1);
tmp2 = dat[dat$s1>=2,]; m2=mean(tmp2$y2);
tmp3 = dat[dat$s1>=4,]; m3=mean(tmp3$y3);
tmp4 = dat[dat$s1>=6,]; m4=mean(tmp4$y4);
tmp5 = dat[dat$s1>=8,]; m5=mean(tmp5$y5);

b0=m1; b1=m2-b0; b2=m3-b0; b3=m4-b0; b4=m5-b0;
out = c(b0,b1,b2,b3,b4)
out
}
test = foreach(m=1:1000) %dopar% myfunc(m)
test2 = do.call("rbind", test)



write.csv(test2,"true_mean.csv")
