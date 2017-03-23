NRd = function(crit,beta)
{
  beta = beta
  diff=10;r=1;
  while(diff>=crit)
  {
    phi = phic = 1;
    U = 0; dU=0;
    auniqueid = unique(midatsortd$Case);
    for(i in auniqueid)
    {
      tmp5 = midatsortd[midatsortd$Case==i,];       
      n = nrow(tmp5);
      if(tmp5$Smoke[1]==0){smokestat1 = rep(0,n)} else if(tmp5$Smoke[1]==1){smokestat1 = rep(1,n);}
      x = matrix(c(rep(1, n),tmp5$time, tmp5$time^2, tmp5$Female, tmp5$shiftage, tmp5$Educyrs, smokestat1, tmp5$time*tmp5$Female,tmp5$time*tmp5$shiftage,tmp5$time*tmp5$Educyrs, tmp5$time*smokestat1),nrow=n)
      last = tmp5[n,]; last$time=10; last2 = do.call("rbind",replicate(5-n, last, simplify = FALSE));tmp5 = rbind(tmp5,last2);
      
      ##j=1
      xtmp = as.data.frame(matrix(c(rep(tmp5$digreal[1],4), tmp5$shiftage[1], tmp5$Educyrs[1],tmp5$Female[1],tmp5$s1[1],tmp5$mmse1[1],tmp5$Smoke[1],tmp5$height[1],tmp5$iadl1[1],tmp5$HlthPrv1[1]),nrow=1))
      colnames(xtmp) = c("yl","lyl","llyl","lllyl","shiftage","Educyrs","Female","s1","mmse1","Smoke","iadl1","HlthPrv1")
      
      W1 = NULL;
      if(!is.na(tmp5$digreal[1]) & tmp5$s1[1]==0) {W1=c(tmp5$digreal[1],predict(mifit1, newdata = xtmp, type="response"), predict(mifit5, newdata = xtmp, type="response"),predict(mifit8, newdata = xtmp, type="response"),predict(mifit10, newdata = xtmp, type="response"))[1:n]}else{W1=W1}
      
      if(!is.na(tmp5$digreal[1]) & tmp5$s1[1]==2) {W1=c(tmp5$digreal[1],predict(mifit1, newdata = xtmp, type="response"), predict(mifit5, newdata = xtmp, type="response"),predict(mifit8, newdata = xtmp, type="response"),predict(mifit10, newdata = xtmp, type="response"))[1:n]}else{W1=W1}
      
      if(!is.na(tmp5$digreal[1]) & tmp5$s1[1]==4) {W1=c(tmp5$digreal[1],predict(mifit12, newdata = xtmp, type="response"), predict(mifit5, newdata = xtmp, type="response"),predict(mifit8, newdata = xtmp, type="response"),predict(mifit10, newdata = xtmp, type="response"))[1:n]}else{W1=W1}
      
      if(!is.na(tmp5$digreal[1]) & tmp5$s1[1]==6) {W1=c(tmp5$digreal[1],predict(mifit13, newdata = xtmp, type="response"), predict(mifit52, newdata = xtmp, type="response"),predict(mifit8, newdata = xtmp, type="response"),predict(mifit10, newdata = xtmp, type="response"))[1:n]}else{W1=W1}
      
      if(!is.na(tmp5$digreal[1]) & tmp5$s1[1]==8) {W1=c(tmp5$digreal[1],predict(mifit14, newdata = xtmp, type="response"), predict(mifit53, newdata = xtmp, type="response"),predict(mifit82, newdata = xtmp, type="response"),predict(mifit10, newdata = xtmp, type="response"))[1:n]}else{W1=W1}
      
      #j=2
      xtmp2l = data.frame(cbind(tmp5$digreal[2],tmp5$digreal[1],tmp5$shiftage[1],tmp5$Educyrs[1],tmp5$Female[1],tmp5$s1[1],tmp5$mmse1[1],tmp5$Smoke[1],tmp5$height[1],tmp5$iadl1[1],tmp5$HlthPrv1[1]));
      colnames(xtmp2l) = c("yl","lyl","shiftage","Educyrs","Female","s1","mmse1","Smoke","iadl1","HlthPrv1")
      xtmp2 = data.frame(cbind(tmp5$digreal[2],tmp5$digreal[1],tmp5$shiftage[1],tmp5$Educyrs[1],tmp5$Female[1],tmp5$s1[1],tmp5$mmse1[1],tmp5$Smoke[1],tmp5$height[1],tmp5$iadl1[1],tmp5$HlthPrv1[1]));
      colnames(xtmp2) = c("lyl","llyl","shiftage","Educyrs","Female","s1","mmse1","Smoke","iadl1","HlthPrv1")
      xtmp2ll = data.frame(cbind(tmp5$digreal[2],tmp5$digreal[1],tmp5$shiftage[1],tmp5$Educyrs[1],tmp5$Female[1],tmp5$mmse1[1],tmp5$Smoke[1],tmp5$height[1],tmp5$iadl1[1],tmp5$HlthPrv1[1]));
      colnames(xtmp2ll) = c("llyl","lllyl","shiftage","Educyrs","Female","mmse1","Smoke","iadl1","HlthPrv1")
      
      W2=NULL
      if(!is.na(tmp5$digreal[2]) & tmp5$s1[1]==2){W2=c(tmp5$digreal[1], tmp5$digreal[2],predict(mifit2, newdata=xtmp2l,type="response"),predict(mifit6, newdata=xtmp2,type="response"),predict(mifit9, newdata=xtmp2ll,type="response"))[1:n]}else{W2=W2}
      
      if(!is.na(tmp5$digreal[2]) & tmp5$s1[1]==4){W2=c(tmp5$digreal[1], tmp5$digreal[2],predict(mifit2, newdata=xtmp2l,type="response"),predict(mifit6, newdata=xtmp2,type="response"),predict(mifit9, newdata=xtmp2ll,type="response"))[1:n]}else{W2=W2}
      
      if(!is.na(tmp5$digreal[2]) & tmp5$s1[1]==6){W2=c(tmp5$digreal[1], tmp5$digreal[2],predict(mifit22, newdata=xtmp2l,type="response"),predict(mifit6, newdata=xtmp2,type="response"),predict(mifit9, newdata=xtmp2ll,type="response"))[1:n]}else{W2=W2}
      
      if(!is.na(tmp5$digreal[2]) & tmp5$s1[1]==8){W2=c(tmp5$digreal[1], tmp5$digreal[2],predict(mifit23, newdata=xtmp2l,type="response"),predict(mifit62, newdata=xtmp2,type="response"),predict(mifit9, newdata=xtmp2ll,type="response"))[1:n]}else{W2=W2}
      
      #j=3
      xtmp3 = data.frame(cbind(tmp5$digreal[3],tmp5$digreal[2],tmp5$digreal[1],tmp5$shiftage[1],tmp5$Educyrs[1],tmp5$Female[1],tmp5$s1[1],tmp5$mmse1[1],tmp5$Smoke[1],tmp5$height[1],tmp5$iadl1[1],tmp5$HlthPrv1[1]));
      colnames(xtmp3) = c("yl","lyl","llyl","shiftage","Educyrs","Female","s1","mmse1","Smoke","iadl1","HlthPrv1")
      xtmp3l = data.frame(cbind(tmp5$digreal[3],tmp5$digreal[2],tmp5$digreal[1],tmp5$shiftage[1],tmp5$Educyrs[1],tmp5$Female[1],tmp5$mmse1[1],tmp5$Smoke[1],tmp5$height[1],tmp5$iadl1[1],tmp5$HlthPrv1[1]));
      colnames(xtmp3l) = c("lyl","llyl","lllyl","shiftage","Educyrs","Female","mmse1","Smoke","iadl1","HlthPrv1")
      
      W3=NULL
      if(!is.na(tmp5$digreal[3]) & tmp5$s1[1]==6){W3=c(tmp5$digreal[1], tmp5$digreal[2],tmp5$digreal[3],predict(mifit3, newdata=xtmp3,type="response"),predict(mifit7, newdata=xtmp3l,type="response"))[1:n]}else{W3=W3}
      
      if(!is.na(tmp5$digreal[3]) & tmp5$s1[1]==8){W3=c(tmp5$digreal[1], tmp5$digreal[2],tmp5$digreal[3],predict(mifit32, newdata=xtmp3,type="response"),predict(mifit7, newdata=xtmp3l,type="response"))[1:n]}else{W3=W3}
      
      #j=4
      xtmp4 = data.frame(cbind(tmp5$digreal[4],tmp5$digreal[3],tmp5$digreal[2],tmp5$digreal[1],tmp5$shiftage[1],tmp5$Educyrs[1],tmp5$Female[1],tmp5$mmse1[1],tmp5$Smoke[1],tmp5$height[1],tmp5$iadl1[1],tmp5$HlthPrv1[1]));
      colnames(xtmp4) = c("yl","lyl","llyl","lllyl","shiftage","Educyrs","Female","mmse1","Smoke","iadl1","HlthPrv1")
      
      W4=NULL
      if(!is.na(tmp5$digreal[3]) & tmp5$s1[1]==8){W4=c(tmp5$digreal[1], tmp5$digreal[2],tmp5$digreal[3],tmp5$digreal[4],predict(mifit4, newdata=xtmp4,type="response"))[1:n]}else{W4=W4}
      
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
    r=r+1
    #if(r==10){diff=0.0000001}
    #cat(r, "\n")
  }
  return(beta)
}