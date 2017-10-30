library(MASS)
library(glmnet)
library(ggplot2)
library(ivmodel)
#z=instruments, c=confounders, x=exposure, y=outcome

#one confounder

#covariance matrix for instruments
dta1<-data.frame(x=1:8)
for (j in 1:8){
  list1<-list()
  for (i in 1:8){
    tl<-0.5^abs(i-j)
    list1[[length(list1)+1]]<-tl
  }
  dta1<-cbind(dta1,unlist(list1))
}
s<-matrix(c(dta1[,2],dta1[,3],dta1[,4],dta1[,5],dta1[,6],dta1[,7],dta1[,8],dta1[,9]),ncol=8,nrow=8)

#OLS method
list12<-list()
list13<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+3*c+2*rnorm(20))
  y<-as.vector(x+c+2*rnorm(20,0,2))
  fit1<-lm(y~x)
  pp2<-summary(fit1)
  pp3<-pp2$coefficients[2,4]
  list12[[length(list12)+1]]<-pp3
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+3*c+2*rnorm(20))
  y<-as.vector(x+c+2*rnorm(20,0,2))
  ypre<-as.numeric(predict(fit1,as.data.frame(x)))
  list13[[length(list13)+1]]<-mean((y-ypre)^2)
}
ap6<-unlist(list12)

#TSL method
list2<-list()
list4<-list()
list6<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-b%*%t(z)+3*c+2*rnorm(20)
  y<-as.vector(x+c+2*rnorm(20,0,2))
  crossval <-  cv.glmnet(y=x,x=z)
  penalty <- crossval$lambda.min
  fit3 <-glmnet( y=x,x=z, alpha = 1, lambda = penalty )
  xhat<-predict(crossval, newx = z, s = penalty)
  xhat<-as.vector(xhat[,1])
  fit4<-lm(y~xhat)
  pp2<-summary(fit4)
  pp3<-1
  if(nrow(coef(summary(fit4)))==2) pp3 = coef(summary(fit4))[2,4]
  tt2<-deviance(fit4)
  list2[[length(list2)+1]]<-tt2
  list4[[length(list4)+1]]<-pp3
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-b%*%t(z)+3*c+2*rnorm(20)
  y<-as.vector(x+c+2*rnorm(20,0,2))
  
  xhat<-as.numeric(predict(fit3 ,z, s=fit3$lambda.min))
  ypre<-as.numeric(predict(fit4,as.data.frame(as.numeric(xhat))))
  list6[[length(list6)+1]]<-mean((y-ypre)^2)
}
ap<-unlist(list4)

#TSLS method
list1<-list()
list3<-list()
list5<-list()
for (i in 1:100){
  set.seed(i)
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+3*c+2*rnorm(20))
  y<-as.vector(x+c+2*rnorm(20,0,2))
  z<-cbind(z,x)
  fit1<-lm(x~.,z)
  xhat<-fit1$fitted.values
  fit2<-lm(y~xhat)
  pp<-summary(fit2)
  pp1<-pp$coefficients[2,4]
  tt1<-deviance(fit2)
  mse1<-mean(pp$residuals^2)
  list1[[length(list1)+1]]<-tt1
  list3[[length(list3)+1]]<-pp1
  
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-b%*%t(z)+3*c+2*rnorm(20)
  y<-as.vector(x+c+2*rnorm(20,0,2))
  z<-cbind(z,x)
  xhat<-as.numeric(predict(fit1,as.data.frame(z[,1:8])))
  ypre<-as.numeric(predict(fit2,as.data.frame(xhat)))
  list5[[length(list5)+1]]<-mean((y-ypre)^2)
}

ap2<-unlist(list3)

#TSPL method
list8<-list()
list9<-list()
list10<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-b%*%t(z)+3*c+2*rnorm(20)
  y<-as.vector(x+c+2*rnorm(20,0,2))
  x = x-mean(x)
  z = t(t(z)-colMeans(z))
  fit5 <-  cv.glmnet(z,x,keep=TRUE,nfolds = length(x),intercept=FALSE)
  lambda = fit5$lambda.min
  xhat = as.numeric(fit5$fit.preval[,which(fit5$lambda==lambda)[1]])
  fit6 = lm(y~xhat)
  rp<-fit6$coefficients[2]
  list10[[length(list10)+1]]<-rp
  pp4 = 1
  if(nrow(coef(summary(fit6)))==2) pp4 = coef(summary(fit6))[2,4]
  list8[[length(list8)+1]]<-pp4
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-b%*%t(z)+3*c+2*rnorm(20)
  y<-as.vector(x+c+2*rnorm(20,0,2))
  x = x-mean(x)
  z = t(t(z)-colMeans(z))
  
  xhat<-as.numeric(predict(fit5 ,z, s=fit5$lambda.min))
  ppre<-as.numeric(predict(fit6,data.frame(xhat)))
  list9[[length(list9)+1]]<-mean((y-ppre)^2)
}

ap4<-unlist(list8)

#MSE plot when there is association btw x and y
lasso.mse<-as.numeric(unlist(list6))
ols.mse<-as.numeric(unlist(list5))
pv.mse<-as.numeric(unlist(list9))
ol.mse<-as.numeric(unlist(list13))
v<-c(rep("TSLS",100),rep("TSL",100),rep("TSPL",100),rep("OLS",100))
mse<-as.data.frame(cbind(v,c(ols.mse,lasso.mse,pv.mse,ol.mse)))
mse[,2]<- as.numeric(as.character(mse[,2]))
ggplot(mse,aes(x=mse$v,y=mse[,2],group=mse$v))+scale_y_continuous(limits = c(0, 500))+geom_boxplot(fill=c("palegreen","darkblue","skyblue","darkgreen"))+labs(x="Model",y="MSE",title="Association between X and Y")+theme(plot.title = element_text(hjust = 0.5))

#true pos and false neg
fnn<-list4>0.05
fn<-length(fnn[fnn=="TRUE"])
tp<-100-fn
fnn1<-list3>0.05
fn1<-length(fnn1[fnn1=="TRUE"])
tp1<-100-fn1
fnn2<-list8>0.05
fn2<-length(fnn2[fnn2=="TRUE"])
tp2<-100-fn2
fnn3<-list12>0.05
fn3<-length(fnn3[fnn3=="TRUE"])
tp3<-100-fn3

##no association
dta1<-data.frame(x=1:8)
for (j in 1:8){
  list1<-list()
  for (i in 1:8){
    tl<-0.5^abs(i-j)
    list1[[length(list1)+1]]<-tl
  }
  dta1<-cbind(dta1,unlist(list1))
}
s<-matrix(c(dta1[,2],dta1[,3],dta1[,4],dta1[,5],dta1[,6],dta1[,7],dta1[,8],dta1[,9]),ncol=8,nrow=8)

list12<-list()
list13<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+3*c+2*rnorm(20))
  y<-as.vector(c+2*rnorm(20,0,2))
  fit1<-lm(y~x)
  pp2<-summary(fit1)
  pp3<-pp2$coefficients[2,4]
  list12[[length(list12)+1]]<-pp3
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+3*c+2*rnorm(20))
  y<-as.vector(c+2*rnorm(20,0,2))
  ypre<-as.numeric(predict(fit1,as.data.frame(x)))
  list13[[length(list13)+1]]<-mean((y-ypre)^2)
}
ap7<-unlist(list12)

list2<-list()
list4<-list()
list6<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20)
  x<-b%*%t(z)+3*c+2*rnorm(20)
  y<-as.vector(c+2*rnorm(20,0,2))
  crossval <-  cv.glmnet(y=x,x=z)
  penalty <- crossval$lambda.min
  fit3 <-glmnet( y=x,x=z, alpha = 1, lambda = penalty )
  xhat<-predict(crossval, newx = z, s = penalty)
  xhat<-as.vector(xhat[,1])
  fit4<-lm(y~xhat)
  pp2<-summary(fit4)
  pp3<-1
  if(nrow(coef(summary(fit4)))==2) pp3 = coef(summary(fit4))[2,4]
  tt2<-deviance(fit4)
  list2[[length(list2)+1]]<-tt2
  list4[[length(list4)+1]]<-pp3
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-b%*%t(z)+3*c+2*rnorm(20)
  y<-as.vector(c+2*rnorm(20,0,2))
  
  xhat<-as.numeric(predict(fit3 ,z, s=fit3$lambda.min))
  ypre<-as.numeric(predict(fit4,as.data.frame(as.numeric(xhat))))
  list6[[length(list6)+1]]<-mean((y-ypre)^2)
}
ap1<-unlist(list4)


list1<-list()
list3<-list()
list5<-list()
for (i in 1:100){
  set.seed(i)
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+3*c+2*rnorm(20))
  y<-as.vector(c+2*rnorm(20,0,2))
  z<-cbind(z,x)
  fit1<-lm(x~.,z)
  xhat<-fit1$fitted.values
  fit2<-lm(y~xhat)
  pp<-summary(fit2)
  pp1<-pp$coefficients[2,4]
  tt1<-deviance(fit2)
  mse1<-mean(pp$residuals^2)
  list1[[length(list1)+1]]<-tt1
  list3[[length(list3)+1]]<-pp1
  
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-b%*%t(z)+3*c+2*rnorm(20)
  y<-as.vector(c+2*rnorm(20,0,2))
  z<-cbind(z,x)
  xhat<-as.numeric(predict(fit1,as.data.frame(z[,1:8])))
  ypre<-as.numeric(predict(fit2,as.data.frame(xhat)))
  list5[[length(list5)+1]]<-mean((y-ypre)^2)
}

ap3<-unlist(list3)



list8<-list()
list9<-list()
list10<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-b%*%t(z)+3*c+2*rnorm(20)
  y<-as.vector(c+2*rnorm(20,0,2))
  x = x-mean(x)
  z = t(t(z)-colMeans(z))
  fit5 <-  cv.glmnet(z,x,keep=TRUE,nfolds = length(x),intercept=FALSE)
  lambda = fit5$lambda.min
  xhat = as.numeric(fit5$fit.preval[,which(fit5$lambda==lambda)[1]])
  fit6 = lm(y~xhat)
  rp<-fit6$coefficients[2]
  list10[[length(list10)+1]]<-rp
  pp4 = 1
  if(nrow(coef(summary(fit6)))==2) pp4 = coef(summary(fit6))[2,4]
  list8[[length(list8)+1]]<-pp4
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-b%*%t(z)+3*c+2*rnorm(20)
  y<-as.vector(c+2*rnorm(20,0,2))
  x = x-mean(x)
  z = t(t(z)-colMeans(z))
  
  xhat<-as.numeric(predict(fit5 ,z, s=fit5$lambda.min))
  ppre<-as.numeric(predict(fit6,data.frame(xhat)))
  list9[[length(list9)+1]]<-mean((y-ppre)^2)
}

ap5<-unlist(list8)

#MSE plot when x and y are not associated
lasso.mse<-as.numeric(unlist(list6))
ols.mse<-as.numeric(unlist(list5))
pv.mse<-as.numeric(unlist(list9))
ol.mse<-as.numeric(unlist(list13))
v<-c(rep("TSLS",100),rep("TSL",100),rep("TSPL",100),rep("OLS",100))
mse<-as.data.frame(cbind(v,c(ols.mse,lasso.mse,pv.mse,ol.mse)))
mse[,2]<- as.numeric(as.character(mse[,2]))
ggplot(mse,aes(x=mse$v,y=mse[,2],group=mse$v))+scale_y_continuous(limits = c(0, 500))+geom_boxplot(fill=c("palegreen","darkblue","skyblue","darkgreen"))+labs(x="Model",y="MSE",title="No Association between X and Y")+theme(plot.title = element_text(hjust = 0.5))

tnn<-list4>0.05
tn<-length(tnn[tnn=="TRUE"])
fp<-100-tn
tnn1<-list3>0.05
tn1<-length(tnn1[tnn1=="TRUE"])
fp1<-100-tn1
tnn2<-list8>0.05
tn2<-length(tnn2[tnn2=="TRUE"])
fp2<-100-tn2
tnn3<-list12>0.05
tn3<-length(tnn3[tnn3=="TRUE"])
fp3<-100-tn3

#make MCC function
mcc<-function(a,b,c,d){
  (a*b-c*d)/sqrt((a+c)*(b+d)*(c+b)*(a+d))
}
#a=tp b=tn c=fp d=fn

#Matthews coefficient correlation
mcc(tp,tn,fp,fn)
mcc(tp1,tn1,fp1,fn1)
mcc(tp2,tn2,fp2,fn2)
mcc(tp3,tn3,fp3,fn3)

#box plot of p-values
df<-as.data.frame(cbind(as.numeric(c(ap,ap1,ap2,ap3,ap4,ap5)),as.factor(c(rep("Ass",100),rep("No Ass",100))),(c(rep("TSL",200),rep("TSLS",200),rep("TSPL",200)))))
colnames(df)<-c("p.value","ass","method")
df[,1]<- as.numeric(as.character(df[,1]))
df<-df[c(1:100,200:300,400:500),]
ggplot(df,aes(x=df$method,y=df$p.value,group=df$method))+labs(x="Method",y="p-value",title="Range of p-values from the Simulations")+theme(plot.title = element_text(hjust = 0.5))+geom_boxplot(fill=c("darkblue","skyblue","darkgreen"))

#ROC plot
app<-c(ap,ap1)
class = rep(c(1,0),c(length(ap1),length(ap)))
ord = order(app)
plot(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),type = 'l',main = 'ROC curve',xlab = 'FPR',ylab = 'TPR',col="darkblue",lwd=4)
abline(0,1)
p = c(ap2,ap3)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="skyblue",type = 'l',lwd=4)
p = c(ap4,ap5)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="darkgreen",type = 'l',lwd=4)
p = c(ap6,ap7)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="palegreen",type = 'l',lwd=4)
legend('bottomright',c('TSL','TSLS','TSPL','OLS'),col = c("darkblue","skyblue","darkgreen","palegreen"),lty = 1,cex=1,lwd=3)


##no confounders 
dta1<-data.frame(x=1:8)
for (j in 1:8){
  list1<-list()
  for (i in 1:8){
    tl<-0.5^abs(i-j)
    list1[[length(list1)+1]]<-tl
  }
  dta1<-cbind(dta1,unlist(list1))
}
s<-matrix(c(dta1[,2],dta1[,3],dta1[,4],dta1[,5],dta1[,6],dta1[,7],dta1[,8],dta1[,9]),ncol=8,nrow=8)

list12<-list()
list13<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  x<-as.vector(b%*%t(z)+2*rnorm(20,0,4))
  y<-as.vector(x+2*rnorm(20,0,6))
  fit1<-lm(y~x)
  pp2<-summary(fit1)
  pp3<-pp2$coefficients[2,4]
  list12[[length(list12)+1]]<-pp3
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  x<-as.vector(b%*%t(z)+2*rnorm(20,0,4))
  y<-as.vector(x+2*rnorm(20,0,6))
  ypre<-as.numeric(predict(fit1,as.data.frame(x)))
  list13[[length(list13)+1]]<-mean((y-ypre)^2)
}
ap6<-unlist(list12)

list2<-list()
list4<-list()
list6<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  x<-as.vector(b%*%t(z)+2*rnorm(20,0,4))
  y<-as.vector(x+2*rnorm(20,0,6))
  crossval <-  cv.glmnet(y=x,x=z)
  penalty <- crossval$lambda.min
  fit3 <-glmnet( y=x,x=z, alpha = 1, lambda = penalty )
  xhat<-predict(crossval, newx = z, s = penalty)
  xhat<-as.vector(xhat[,1])
  fit4<-lm(y~xhat)
  pp2<-summary(fit4)
  pp3 = 1
  if(nrow(coef(summary(fit4)))==2) pp3 = coef(summary(fit4))[2,4]
  tt2<-deviance(fit4)
  list2[[length(list2)+1]]<-tt2
  list4[[length(list4)+1]]<-pp3
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  x<-as.vector(b%*%t(z)+2*rnorm(20,0,4))
  y<-as.vector(x+2*rnorm(20,0,6))
  
  xhat<-as.numeric(predict(fit3 ,z, s=fit3$lambda.min))
  ypre<-as.numeric(predict(fit4,as.data.frame(as.numeric(xhat))))
  list6[[length(list6)+1]]<-mean((y-ypre)^2)
}
ap<-unlist(list4)

list1<-list()
list3<-list()
list5<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  x<-as.vector(b%*%t(z)+2*rnorm(20,0,4))
  y<-as.vector(x+2*rnorm(20,0,6))
  z<-as.data.frame(cbind(z,x))
  fit1<-lm(x~.,z)
  xhat<-fit1$fitted.values
  fit2<-lm(y~xhat)
  pp<-summary(fit2)
  pp1<-pp$coefficients[2,4]
  tt1<-deviance(fit2)
  mse1<-mean(pp$residuals^2)
  list1[[length(list1)+1]]<-tt1
  list3[[length(list3)+1]]<-pp1
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  x<-as.vector(b%*%t(z)+2*rnorm(20,0,4))
  y<-as.vector(x+2*rnorm(20,0,6))
  z<-cbind(z,x)
  xhat<-as.numeric(predict(fit1,as.data.frame(z[,1:8])))
  ypre<-as.numeric(predict(fit2,as.data.frame(xhat)))
  list5[[length(list5)+1]]<-mean((y-ypre)^2)
}
ap2<-unlist(list3)


list8<-list()
list9<-list()
list10<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  x<-as.vector(b%*%t(z)+2*rnorm(20,0,4))
  y<-as.vector(x+2*rnorm(20,0,6))
  x = x-mean(x)
  z = t(t(z)-colMeans(z))
  fit5 <-  cv.glmnet(z,x,keep=TRUE,nfolds = length(x),intercept=FALSE)
  lambda = fit5$lambda.min
  xhat = as.numeric(fit5$fit.preval[,which(fit5$lambda==lambda)[1]])
  fit6 = lm(y~xhat)
  rp<-fit6$coefficients[2]
  list10[[length(list10)+1]]<-rp
  pp4 = 1
  if(nrow(coef(summary(fit6)))==2) pp4 = coef(summary(fit6))[2,4]
  list8[[length(list8)+1]]<-pp4
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  x<-as.vector(b%*%t(z)+2*rnorm(20,0,4))
  y<-as.vector(x+2*rnorm(20,0,6))
  x = x-mean(x)
  z = t(t(z)-colMeans(z))
  
  xhat<-as.numeric(predict(fit5 ,z, s=fit5$lambda.min))
  ppre<-as.numeric(predict(fit6,data.frame(xhat)))
  list9[[length(list9)+1]]<-mean((y-ppre)^2)
}

ap4<-unlist(list8)

lasso.mse<-as.numeric(unlist(list6))
ols.mse<-as.numeric(unlist(list5))
pv.mse<-as.numeric(unlist(list9))
ol.mse<-as.numeric(unlist(list13))
v<-c(rep("TSLS",100),rep("TSL",100),rep("TSPL",100),rep("OLS",100))
mse<-as.data.frame(cbind(v,c(ols.mse,lasso.mse,pv.mse,ol.mse)))
mse[,2]<- as.numeric(as.character(mse[,2]))
ggplot(mse,aes(x=mse$v,y=mse[,2],group=mse$v))+scale_y_continuous(limits = c(0, 90))+geom_boxplot(fill=c("palegreen","darkblue","skyblue","darkgreen"))+labs(x="Model",y="MSE",title="Association between X and Y")+theme(plot.title = element_text(hjust = 0.5))

fnn<-list4>0.05
fn<-length(fnn[fnn=="TRUE"])
tp<-100-fn
fnn1<-list3>0.05
fn1<-length(fnn1[fnn1=="TRUE"])
tp1<-100-fn1
fnn2<-list8>0.05
fn2<-length(fnn2[fnn2=="TRUE"])
tp2<-100-fn2
fnn3<-list12>0.05
fn3<-length(fnn3[fnn3=="TRUE"])
tp3<-100-fn3


dta1<-data.frame(x=1:8)
for (j in 1:8){
  list1<-list()
  for (i in 1:8){
    tl<-0.5^abs(i-j)
    list1[[length(list1)+1]]<-tl
  }
  dta1<-cbind(dta1,unlist(list1))
}
s<-matrix(c(dta1[,2],dta1[,3],dta1[,4],dta1[,5],dta1[,6],dta1[,7],dta1[,8],dta1[,9]),ncol=8,nrow=8)

list12<-list()
list13<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  x<-as.vector(b%*%t(z)+2*rnorm(20,0,4))
  y<-as.vector(2*rnorm(20,0,6))
  fit1<-lm(y~x)
  pp2<-summary(fit1)
  pp3<-pp2$coefficients[2,4]
  list12[[length(list12)+1]]<-pp3
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  x<-as.vector(b%*%t(z)+2*rnorm(20,0,4))
  y<-as.vector(2*rnorm(20,0,6))
  ypre<-as.numeric(predict(fit1,as.data.frame(x)))
  list13[[length(list13)+1]]<-mean((y-ypre)^2)
}
ap7<-unlist(list12)


list2<-list()
list4<-list()
list6<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  x<-as.vector(b%*%t(z)+2*rnorm(20,0,4))
  y<-as.vector(2*rnorm(20,0,6))
  crossval <-  cv.glmnet(y=x,x=z)
  penalty <- crossval$lambda.min
  fit3 <-glmnet( y=x,x=z, alpha = 1, lambda = penalty )
  xhat<-predict(crossval, newx = z, s = penalty)
  xhat<-as.vector(xhat[,1])
  fit4<-lm(y~xhat)
  pp2<-summary(fit4)
  pp3<-1
  if(nrow(coef(summary(fit4)))==2) pp3 = coef(summary(fit4))[2,4]
  tt2<-deviance(fit4)
  list2[[length(list2)+1]]<-tt2
  list4[[length(list4)+1]]<-pp3
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  x<-as.vector(b%*%t(z)+2*rnorm(20,0,4))
  y<-as.vector(2*rnorm(20,0,6))
  
  xhat<-as.numeric(predict(fit3 ,z, s=fit3$lambda.min))
  ypre<-as.numeric(predict(fit4,as.data.frame(as.numeric(xhat))))
  list6[[length(list6)+1]]<-mean((y-ypre)^2)
}
ap1<-unlist(list4)

list1<-list()
list3<-list()
list5<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  x<-as.vector(b%*%t(z)+2*rnorm(20,0,4))
  y<-as.vector(2*rnorm(20,0,2))
  z<-as.data.frame(cbind(z,x))
  fit1<-lm(x~.,z)
  xhat<-fit1$fitted.values
  fit2<-lm(y~xhat)
  pp<-summary(fit2)
  pp1<-pp$coefficients[2,4]
  tt1<-deviance(fit2)
  mse1<-mean(pp$residuals^2)
  list1[[length(list1)+1]]<-tt1
  list3[[length(list3)+1]]<-pp1
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  x<-as.vector(b%*%t(z)+2*rnorm(20,0,4))
  y<-as.vector(2*rnorm(20,0,6))
  z<-cbind(z,x)
  xhat<-as.numeric(predict(fit1,as.data.frame(z[,1:8])))
  ypre<-as.numeric(predict(fit2,as.data.frame(xhat)))
  list5[[length(list5)+1]]<-mean((y-ypre)^2)
}
ap3<-unlist(list3)


list8<-list()
list9<-list()
list10<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  x<-as.vector(b%*%t(z)+2*rnorm(20,0,4))
  y<-as.vector(2*rnorm(20,0,6))
  x = x-mean(x)
  z = t(t(z)-colMeans(z))
  fit5 <-  cv.glmnet(z,x,keep=TRUE,nfolds = length(x),intercept=FALSE)
  lambda = fit5$lambda.min
  xhat = as.numeric(fit5$fit.preval[,which(fit5$lambda==lambda)[1]])
  fit6 = lm(y~xhat)
  rp<-fit6$coefficients[2]
  list10[[length(list10)+1]]<-rp
  pp4 = 1
  if(nrow(coef(summary(fit6)))==2) pp4 = coef(summary(fit6))[2,4]
  list8[[length(list8)+1]]<-pp4
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  x<-as.vector(b%*%t(z)+2*rnorm(20,0,4))
  y<-as.vector(2*rnorm(20,0,6))
  x = x-mean(x)
  z = t(t(z)-colMeans(z))
  
  xhat<-as.numeric(predict(fit5 ,z, s=fit5$lambda.min))
  ppre<-as.numeric(predict(fit6,data.frame(xhat)))
  list9[[length(list9)+1]]<-mean((y-ppre)^2)
}

ap5<-unlist(list8)

lasso.mse<-as.numeric(unlist(list6))
ols.mse<-as.numeric(unlist(list5))
pv.mse<-as.numeric(unlist(list9))
ol.mse<-as.numeric(unlist(list13))
v<-c(rep("TSLS",100),rep("TSL",100),rep("TSPL",100),rep("OLS",100))
mse<-as.data.frame(cbind(v,c(ols.mse,lasso.mse,pv.mse,ol.mse)))
mse[,2]<- as.numeric(as.character(mse[,2]))
ggplot(mse,aes(x=mse$v,y=mse[,2],group=mse$v))+scale_y_continuous(limits = c(0, 90))+geom_boxplot(fill=c("palegreen","darkblue","skyblue","darkgreen"))+labs(x="Model",y="MSE",title="No Association between X and Y")+theme(plot.title = element_text(hjust = 0.5))

tnn<-list4>0.05
tn<-length(tnn[tnn=="TRUE"])
fp<-100-tn
tnn1<-list3>0.05
tn1<-length(tnn1[tnn1=="TRUE"])
fp1<-100-tn1
tnn2<-list8>0.05
tn2<-length(tnn2[tnn2=="TRUE"])
fp2<-100-tn2
tnn3<-list12>0.05
tn3<-length(tnn3[tnn3=="TRUE"])
fp3<-100-tn3

mcc(tp,tn,fp,fn)
mcc(tp1,tn1,fp1,fn1)
mcc(tp2,tn2,fp2,fn2)
mcc(tp3,tn3,fp3,fn3)


#ROC plot
app<-c(ap,ap1)
class = rep(c(1,0),c(length(ap1),length(ap)))
ord = order(app)
plot(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),type = 'l',main = 'ROC curve',xlab = 'FPR',ylab = 'TPR',col="darkblue",lwd=4)
abline(0,1)
p = c(ap2,ap3)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="skyblue",type = 'l',lwd=4)
p = c(ap4,ap5)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="darkgreen",type = 'l',lwd=4)
p = c(ap6,ap7)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="palegreen",type = 'l',lwd=4)
legend('bottomright',c('TSL','TSLS','TSPL','OLS'),col = c("darkblue","skyblue","darkgreen","palegreen"),lty = 1,cex=1,lwd=3)

app<-c(ap,ap1)
class = rep(c(1,0),c(length(ap1),length(ap)))
ord = order(app)
plot(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),type = 'l',main = 'ROC curve',xlab = 'FPR',ylab = 'TPR',col="darkblue",lwd=4)
abline(0,1)
p = c(ap2,ap3)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="skyblue",type = 'l',lwd=4)
p = c(ap4,ap5)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="darkgreen",type = 'l',lwd=4)
legend('bottomright',c('TSL','TSLS','TSPL'),col = c("darkblue","skyblue","darkgreen"),lty = 1,cex=1,lwd=3)

df<-as.data.frame(cbind(as.numeric(c(ap,ap1,ap2,ap3,ap4,ap5)),as.factor(c(rep("Ass",100),rep("No Ass",100))),(c(rep("TSL",200),rep("TSLS",200),rep("TSPL",200)))))
colnames(df)<-c("p.value","ass","method")
df[,1]<- as.numeric(as.character(df[,1]))
df<-df[c(1:100,200:300,400:500),]
ggplot(df,aes(x=df$method,y=df$p.value,group=df$method))+labs(x="Method",y="p-value",title="Range of p-values from the Simulations")+theme(plot.title = element_text(hjust = 0.5))+geom_boxplot(fill=c("darkblue","skyblue","darkgreen"))


##extra noise
dta1<-data.frame(x=1:8)
for (j in 1:8){
  list1<-list()
  for (i in 1:8){
    tl<-0.5^abs(i-j)
    list1[[length(list1)+1]]<-tl
  }
  dta1<-cbind(dta1,unlist(list1))
}
s<-matrix(c(dta1[,2],dta1[,3],dta1[,4],dta1[,5],dta1[,6],dta1[,7],dta1[,8],dta1[,9]),ncol=8,nrow=8)


list2<-list()
list4<-list()
list6<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+3*c+8*rnorm(20))
  y<-as.vector(x+c+8*rnorm(20,0,2))
  crossval <-  cv.glmnet(y=x,x=z)
  penalty <- crossval$lambda.min
  fit3 <-glmnet( y=x,x=z, alpha = 1, lambda = penalty )
  xhat<-predict(crossval, newx = z, s = penalty)
  xhat<-as.vector(xhat[,1])
  fit4<-lm(y~xhat)
  pp2<-summary(fit4)
  pp3<-1
  if(nrow(coef(summary(fit4)))==2) pp3 = coef(summary(fit4))[2,4]
  tt2<-deviance(fit4)
  list2[[length(list2)+1]]<-tt2
  list4[[length(list4)+1]]<-pp3
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+3*c+8*rnorm(20))
  y<-as.vector(x+c+8*rnorm(20,0,2))
  xhat<-as.numeric(predict(fit3 ,z, s=fit3$lambda.min))
  ypre<-as.numeric(predict(fit4,as.data.frame(as.numeric(xhat))))
  list6[[length(list6)+1]]<-mean((y-ypre)^2)
}
ap<-unlist(list4)

list1<-list()
list3<-list()
list5<-list()
for (i in 1:100){
  set.seed(i)
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+3*c+8*rnorm(20))
  y<-as.vector(x+c+8*rnorm(20,0,2))
  z<-cbind(z,x)
  fit1<-lm(x~.,z)
  xhat<-fit1$fitted.values
  fit2<-lm(y~xhat)
  pp<-summary(fit2)
  pp1<-pp$coefficients[2,4]
  tt1<-deviance(fit2)
  mse1<-mean(pp$residuals^2)
  list1[[length(list1)+1]]<-tt1
  list3[[length(list3)+1]]<-pp1
  
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+3*c+8*rnorm(20))
  y<-as.vector(x+c+8*rnorm(20,0,2))
  z<-cbind(z,x)
  xhat<-as.numeric(predict(fit1,as.data.frame(z[,1:8])))
  ypre<-as.numeric(predict(fit2,as.data.frame(xhat)))
  list5[[length(list5)+1]]<-mean((y-ypre)^2)
}
ap2<-unlist(list3)


list8<-list()
list9<-list()
list10<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+3*c+8*rnorm(20))
  y<-as.vector(x+c+8*rnorm(20,0,2))
  x = x-mean(x)
  z = t(t(z)-colMeans(z))
  fit5 <-  cv.glmnet(z,x,keep=TRUE,nfolds = length(x),intercept=FALSE)
  lambda = fit5$lambda.min
  xhat = as.numeric(fit5$fit.preval[,which(fit5$lambda==lambda)[1]])
  fit6 = lm(y~xhat)
  rp<-fit6$coefficients[2]
  list10[[length(list10)+1]]<-rp
  pp4 = 1
  if(nrow(coef(summary(fit6)))==2) pp4 = coef(summary(fit6))[2,4]
  list8[[length(list8)+1]]<-pp4
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+3*c+8*rnorm(20))
  y<-as.vector(x+c+8*rnorm(20,0,2))
  x = x-mean(x)
  z = t(t(z)-colMeans(z))
  
  xhat<-as.numeric(predict(fit5 ,z, s=fit5$lambda.min))
  ppre<-as.numeric(predict(fit6,data.frame(xhat)))
  list9[[length(list9)+1]]<-mean((y-ppre)^2)
}

ap4<-unlist(list8)

lasso.mse<-as.numeric(unlist(list6))
ols.mse<-as.numeric(unlist(list5))
pv.mse<-as.numeric(unlist(list9))
v<-c(rep("TSLS",100),rep("TSL",100),rep("TSPL",100))
mse<-as.data.frame(cbind(v,c(ols.mse,lasso.mse,pv.mse)))
mse[,2]<- as.numeric(as.character(mse[,2]))
ggplot(mse,aes(x=mse$v,y=mse[,2],group=mse$v))+scale_y_continuous(limits = c(0, 1250))+geom_boxplot(fill=c("darkblue","skyblue","darkgreen"))+labs(x="Model",y="MSE",title="Association between X and Y")+theme(plot.title = element_text(hjust = 0.5))

fnn<-list4>0.05
fn<-length(fnn[fnn=="TRUE"])
tp<-100-fn
fnn1<-list3>0.05
fn1<-length(fnn1[fnn1=="TRUE"])
tp1<-100-fn1
fnn2<-list8>0.05
fn2<-length(fnn2[fnn2=="TRUE"])
tp2<-100-fn2

dta1<-data.frame(x=1:8)
for (j in 1:8){
  list1<-list()
  for (i in 1:8){
    tl<-0.5^abs(i-j)
    list1[[length(list1)+1]]<-tl
  }
  dta1<-cbind(dta1,unlist(list1))
}
s<-matrix(c(dta1[,2],dta1[,3],dta1[,4],dta1[,5],dta1[,6],dta1[,7],dta1[,8],dta1[,9]),ncol=8,nrow=8)


list2<-list()
list4<-list()
list6<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+3*c+8*rnorm(20))
  y<-as.vector(c+8*rnorm(20,0,2))
  crossval <-  cv.glmnet(y=x,x=z)
  penalty <- crossval$lambda.min
  fit3 <-glmnet( y=x,x=z, alpha = 1, lambda = penalty )
  xhat<-predict(crossval, newx = z, s = penalty)
  xhat<-as.vector(xhat[,1])
  fit4<-lm(y~xhat)
  pp2<-summary(fit4)
  pp3<-1
  if(nrow(coef(summary(fit4)))==2) pp3 = coef(summary(fit4))[2,4]
  tt2<-deviance(fit4)
  list2[[length(list2)+1]]<-tt2
  list4[[length(list4)+1]]<-pp3
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+3*c+8*rnorm(20))
  y<-as.vector(c+8*rnorm(20,0,2))
  xhat<-as.numeric(predict(fit3 ,z, s=fit3$lambda.min))
  ypre<-as.numeric(predict(fit4,as.data.frame(as.numeric(xhat))))
  list6[[length(list6)+1]]<-mean((y-ypre)^2)
}
ap1<-unlist(list4)

list1<-list()
list3<-list()
list5<-list()
for (i in 1:100){
  set.seed(i)
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+3*c+8*rnorm(20))
  y<-as.vector(c+8*rnorm(20,0,2))
  z<-cbind(z,x)
  fit1<-lm(x~.,z)
  xhat<-fit1$fitted.values
  fit2<-lm(y~xhat)
  pp<-summary(fit2)
  pp1<-pp$coefficients[2,4]
  tt1<-deviance(fit2)
  mse1<-mean(pp$residuals^2)
  list1[[length(list1)+1]]<-tt1
  list3[[length(list3)+1]]<-pp1
  
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+3*c+8*rnorm(20))
  y<-as.vector(c+8*rnorm(20,0,2))
  z<-cbind(z,x)
  xhat<-as.numeric(predict(fit1,as.data.frame(z[,1:8])))
  ypre<-as.numeric(predict(fit2,as.data.frame(xhat)))
  list5[[length(list5)+1]]<-mean((y-ypre)^2)
}
ap3<-unlist(list3)



list8<-list()
list9<-list()
list10<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+3*c+8*rnorm(20))
  y<-as.vector(c+8*rnorm(20,0,2))
  x = x-mean(x)
  z = t(t(z)-colMeans(z))
  fit5 <-  cv.glmnet(z,x,keep=TRUE,nfolds = length(x),intercept=FALSE)
  lambda = fit5$lambda.min
  xhat = as.numeric(fit5$fit.preval[,which(fit5$lambda==lambda)[1]])
  fit6 = lm(y~xhat)
  rp<-fit6$coefficients[2]
  list10[[length(list10)+1]]<-rp
  pp4 = 1
  if(nrow(coef(summary(fit6)))==2) pp4 = coef(summary(fit6))[2,4]
  list8[[length(list8)+1]]<-pp4
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+3*c+8*rnorm(20))
  y<-as.vector(c+8*rnorm(20,0,2))
  x = x-mean(x)
  z = t(t(z)-colMeans(z))
  
  xhat<-as.numeric(predict(fit5 ,z, s=fit5$lambda.min))
  ppre<-as.numeric(predict(fit6,data.frame(xhat)))
  list9[[length(list9)+1]]<-mean((y-ppre)^2)
}

ap5<-unlist(list8)

lasso.mse<-as.numeric(unlist(list6))
ols.mse<-as.numeric(unlist(list5))
pv.mse<-as.numeric(unlist(list9))
v<-c(rep("TSLS",100),rep("TSL",100),rep("TSPL",100))
mse<-as.data.frame(cbind(v,c(ols.mse,lasso.mse,pv.mse)))
mse[,2]<- as.numeric(as.character(mse[,2]))
ggplot(mse,aes(x=mse$v,y=mse[,2],group=mse$v))+scale_y_continuous(limits = c(0, 1250))+geom_boxplot(fill=c("darkblue","skyblue","darkgreen"))+labs(x="Model",y="MSE",title="No Association between X and Y")+theme(plot.title = element_text(hjust = 0.5))

tnn<-list4>0.05
tn<-length(tnn[tnn=="TRUE"])
fp<-100-tn
tnn1<-list3>0.05
tn1<-length(tnn1[tnn1=="TRUE"])
fp1<-100-tn1
tnn2<-list8>0.05
tn2<-length(tnn2[tnn2=="TRUE"])
fp2<-100-tn2

mcc(tp,tn,fp,fn)
mcc(tp1,tn1,fp1,fn1)
mcc(tp2,tn2,fp2,fn2)

#ROC plot
app<-c(ap,ap1)
class = rep(c(1,0),c(length(ap1),length(ap)))
ord = order(app)
plot(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),type = 'l',main = 'ROC curve',xlab = 'FPR',ylab = 'TPR',col="darkblue",lwd=4)
abline(0,1)
p = c(ap2,ap3)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="skyblue",type = 'l',lwd=4)
p = c(ap4,ap5)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="darkgreen",type = 'l',lwd=4)
legend('bottomright',c('TSL','TSLS','TSPL'),col = c("darkblue","skyblue","darkgreen"),lty = 1,cex=1,lwd=3)

df<-as.data.frame(cbind(as.numeric(c(ap,ap1,ap2,ap3,ap4,ap5)),as.factor(c(rep("Ass",100),rep("No Ass",100))),(c(rep("TSL",200),rep("TSLS",200),rep("TSPL",200)))))
colnames(df)<-c("p.value","ass","method")
df[,1]<- as.numeric(as.character(df[,1]))
df<-df[c(1:100,200:300,400:500),]
ggplot(df,aes(x=df$method,y=df$p.value,group=df$method))+labs(x="Method",y="p-value",title="Range of p-values from the Simulations")+theme(plot.title = element_text(hjust = 0.5))+geom_boxplot(fill=c("darkblue","skyblue","darkgreen"))


#More Ivs
dta1<-data.frame(x=1:40)
for (j in 1:40){
  list1<-list()
  for (i in 1:40){
    tl<-0.5^abs(i-j)
    list1[[length(list1)+1]]<-tl
  }
  dta1<-cbind(dta1,unlist(list1))
}
s<-matrix(c(dta1[,2],dta1[,3],dta1[,4],dta1[,5],dta1[,6],dta1[,7],dta1[,8],dta1[,9],dta1[,10],dta1[,11],dta1[,12],dta1[,13],dta1[,14],dta1[,15],dta1[,16],dta1[,17],dta1[,18],dta1[,19],dta1[,20],dta1[,21],dta1[,22],dta1[,23],dta1[,24],dta1[,25],dta1[,26],dta1[,27],dta1[,28],dta1[,29],dta1[,30],dta1[,31],dta1[,32],dta1[,33],dta1[,34],dta1[,35],dta1[,36],dta1[,37],dta1[,38],dta1[,39],dta1[,40],dta1[,41]),ncol=40,nrow=40)


list2<-list()
list4<-list()
list6<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(100,rep(2,40),s)
  b<-(c(3,1.5,0,0,2,0,0,0,2,0,0,1,2,5,0,0,0,0,1,1,0,0,0,0,1,2,0,0,4,2,1,2,0,0,0,7,0,0,0,1))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(x+c+2*rnorm(100,0,6))
  crossval <-  cv.glmnet(y=x,x=z)
  penalty <- crossval$lambda.min
  fit3 <-glmnet( y=x,x=z, alpha = 1, lambda = penalty )
  xhat<-predict(crossval, newx = z, s = penalty)
  xhat<-as.vector(xhat[,1])
  fit4<-lm(y~xhat)
  pp2<-summary(fit4)
  pp3<-1
  if(nrow(coef(summary(fit4)))==2) pp3 = coef(summary(fit4))[2,4]
  tt2<-deviance(fit4)
  list2[[length(list2)+1]]<-tt2
  list4[[length(list4)+1]]<-pp3
  
  z<-mvrnorm(100,rep(2,40),s)
  b<-(c(3,1.5,0,0,2,0,0,0,2,0,0,1,2,5,0,0,0,0,1,1,0,0,0,0,1,2,0,0,4,2,1,2,0,0,0,7,0,0,0,1))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(x+c+2*rnorm(100,0,6))
  
  xhat<-as.numeric(predict(fit3 ,z, s=fit3$lambda.min))
  ypre<-as.numeric(predict(fit4,as.data.frame(as.numeric(xhat))))
  list6[[length(list6)+1]]<-mean((y-ypre)^2)
}
ap<-unlist(list4)

list1<-list()
list3<-list()
list5<-list()
for (i in 1:100){
  set.seed(i)
  z<-as.data.frame(mvrnorm(100,rep(2,40),s))
  b<-(c(3,1.5,0,0,2,0,0,0,2,0,0,1,2,5,0,0,0,0,1,1,0,0,0,0,1,2,0,0,4,2,1,2,0,0,0,7,0,0,0,1))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(x+c+2*rnorm(100,0,6))
  z<-cbind(z,x)
  fit1<-lm(x~.,z)
  xhat<-fit1$fitted.values
  fit2<-lm(y~xhat)
  pp<-summary(fit2)
  pp1<-pp$coefficients[2,4]
  tt1<-deviance(fit2)
  mse1<-mean(pp$residuals^2)
  list1[[length(list1)+1]]<-tt1
  list3[[length(list3)+1]]<-pp1
  
  z<-as.data.frame(mvrnorm(100,rep(2,40),s))
  b<-(c(3,1.5,0,0,2,0,0,0,2,0,0,1,2,5,0,0,0,0,1,1,0,0,0,0,1,2,0,0,4,2,1,2,0,0,0,7,0,0,0,1))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(x+c+2*rnorm(100,0,6))
  z<-cbind(z,x)
  xhat<-as.numeric(predict(fit1,as.data.frame(z[,1:40])))
  ypre<-as.numeric(predict(fit2,as.data.frame(xhat)))
  list5[[length(list5)+1]]<-mean((y-ypre)^2)
}
ap2<-unlist(list3)


list8<-list()
list9<-list()
list10<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(100,rep(2,40),s)
  b<-(c(3,1.5,0,0,2,0,0,0,2,0,0,1,2,5,0,0,0,0,1,1,0,0,0,0,1,2,0,0,4,2,1,2,0,0,0,7,0,0,0,1))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(x+c+2*rnorm(100,0,6))
  x = x-mean(x)
  z = t(t(z)-colMeans(z))
  fit5 <-  cv.glmnet(z,x,keep=TRUE,nfolds = length(x),intercept=FALSE)
  lambda = fit5$lambda.min
  xhat = as.numeric(fit5$fit.preval[,which(fit5$lambda==lambda)[1]])
  fit6 = lm(y~xhat)
  rp<-fit6$coefficients[2]
  list10[[length(list10)+1]]<-rp
  pp4 = 1
  if(nrow(coef(summary(fit6)))==2) pp4 = coef(summary(fit6))[2,4]
  list8[[length(list8)+1]]<-pp4
  
  z<-mvrnorm(100,rep(2,40),s)
  b<-(c(3,1.5,0,0,2,0,0,0,2,0,0,1,2,5,0,0,0,0,1,1,0,0,0,0,1,2,0,0,4,2,1,2,0,0,0,7,0,0,0,1))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(x+c+2*rnorm(100,0,6))
  x = x-mean(x)
  z = t(t(z)-colMeans(z))
  
  xhat<-as.numeric(predict(fit5 ,z, s=fit5$lambda.min))
  ppre<-as.numeric(predict(fit6,data.frame(xhat)))
  list9[[length(list9)+1]]<-mean((y-ppre)^2)
}

ap4<-unlist(list8)

lasso.mse<-as.numeric(unlist(list6))
ols.mse<-as.numeric(unlist(list5))
pv.mse<-as.numeric(unlist(list9))
v<-c(rep("TSLS",100),rep("TSL",100),rep("TSPL",100))
mse<-as.data.frame(cbind(v,c(ols.mse,lasso.mse,pv.mse)))
mse[,2]<- as.numeric(as.character(mse[,2]))
ggplot(mse,aes(x=mse$v,y=mse[,2],group=mse$v))+scale_y_continuous(limits = c(0, 25000))+geom_boxplot(fill=c("darkblue","skyblue","darkgreen"))+labs(x="Model",y="MSE",title="Association between X and Y")+theme(plot.title = element_text(hjust = 0.5))

fnn<-list4>0.05
fn<-length(fnn[fnn=="TRUE"])
tp<-100-fn
fnn1<-list3>0.05
fn1<-length(fnn1[fnn1=="TRUE"])
tp1<-100-fn1
fnn2<-list8>0.05
fn2<-length(fnn2[fnn2=="TRUE"])
tp2<-100-fn2


dta1<-data.frame(x=1:40)
for (j in 1:40){
  list1<-list()
  for (i in 1:40){
    tl<-0.5^abs(i-j)
    list1[[length(list1)+1]]<-tl
  }
  dta1<-cbind(dta1,unlist(list1))
}
s<-matrix(c(dta1[,2],dta1[,3],dta1[,4],dta1[,5],dta1[,6],dta1[,7],dta1[,8],dta1[,9],dta1[,10],dta1[,11],dta1[,12],dta1[,13],dta1[,14],dta1[,15],dta1[,16],dta1[,17],dta1[,18],dta1[,19],dta1[,20],dta1[,21],dta1[,22],dta1[,23],dta1[,24],dta1[,25],dta1[,26],dta1[,27],dta1[,28],dta1[,29],dta1[,30],dta1[,31],dta1[,32],dta1[,33],dta1[,34],dta1[,35],dta1[,36],dta1[,37],dta1[,38],dta1[,39],dta1[,40],dta1[,41]),ncol=40,nrow=40)


list2<-list()
list4<-list()
list6<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(100,rep(2,40),s)
  b<-(c(3,1.5,0,0,2,0,0,0,2,0,0,1,2,5,0,0,0,0,1,1,0,0,0,0,1,2,0,0,4,2,1,2,0,0,0,7,0,0,0,1))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(c+2*rnorm(100,0,6))
  crossval <-  cv.glmnet(y=x,x=z)
  penalty <- crossval$lambda.min
  fit3 <-glmnet( y=x,x=z, alpha = 1, lambda = penalty )
  xhat<-predict(crossval, newx = z, s = penalty)
  xhat<-as.vector(xhat[,1])
  fit4<-lm(y~xhat)
  pp2<-summary(fit4)
  pp3<-1
  if(nrow(coef(summary(fit4)))==2) pp3 = coef(summary(fit4))[2,4]
  tt2<-deviance(fit4)
  list2[[length(list2)+1]]<-tt2
  list4[[length(list4)+1]]<-pp3
  
  z<-mvrnorm(100,rep(2,40),s)
  b<-(c(3,1.5,0,0,2,0,0,0,2,0,0,1,2,5,0,0,0,0,1,1,0,0,0,0,1,2,0,0,4,2,1,2,0,0,0,7,0,0,0,1))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(c+2*rnorm(100,0,6))
  
  xhat<-as.numeric(predict(fit3 ,z, s=fit3$lambda.min))
  ypre<-as.numeric(predict(fit4,as.data.frame(as.numeric(xhat))))
  list6[[length(list6)+1]]<-mean((y-ypre)^2)
}
ap1<-unlist(list4)

list1<-list()
list3<-list()
list5<-list()
for (i in 1:100){
  set.seed(i)
  z<-as.data.frame(mvrnorm(100,rep(2,40),s))
  b<-(c(3,1.5,0,0,2,0,0,0,2,0,0,1,2,5,0,0,0,0,1,1,0,0,0,0,1,2,0,0,4,2,1,2,0,0,0,7,0,0,0,1))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(c+2*rnorm(100,0,6))
  z<-cbind(z,x)
  fit1<-lm(x~.,z)
  xhat<-fit1$fitted.values
  fit2<-lm(y~xhat)
  pp<-summary(fit2)
  pp1<-pp$coefficients[2,4]
  tt1<-deviance(fit2)
  mse1<-mean(pp$residuals^2)
  list1[[length(list1)+1]]<-tt1
  list3[[length(list3)+1]]<-pp1
  
  z<-as.data.frame(mvrnorm(100,rep(2,40),s))
  b<-(c(3,1.5,0,0,2,0,0,0,2,0,0,1,2,5,0,0,0,0,1,1,0,0,0,0,1,2,0,0,4,2,1,2,0,0,0,7,0,0,0,1))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(c+2*rnorm(100,0,6))
  z<-cbind(z,x)
  xhat<-as.numeric(predict(fit1,as.data.frame(z[,1:40])))
  ypre<-as.numeric(predict(fit2,as.data.frame(xhat)))
  list5[[length(list5)+1]]<-mean((y-ypre)^2)
}
ap3<-unlist(list3)


list8<-list()
list9<-list()
list10<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(100,rep(2,40),s)
  b<-(c(3,1.5,0,0,2,0,0,0,2,0,0,1,2,5,0,0,0,0,1,1,0,0,0,0,1,2,0,0,4,2,1,2,0,0,0,7,0,0,0,1))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(c+2*rnorm(100,0,6))
  x = x-mean(x)
  z = t(t(z)-colMeans(z))
  fit5 <-  cv.glmnet(z,x,keep=TRUE,nfolds = length(x),intercept=FALSE)
  lambda = fit5$lambda.min
  xhat = as.numeric(fit5$fit.preval[,which(fit5$lambda==lambda)[1]])
  fit6 = lm(y~xhat)
  rp<-fit6$coefficients[2]
  list10[[length(list10)+1]]<-rp
  pp4 = 1
  if(nrow(coef(summary(fit6)))==2) pp4 = coef(summary(fit6))[2,4]
  list8[[length(list8)+1]]<-pp4
  
  z<-mvrnorm(100,rep(2,40),s)
  b<-(c(3,1.5,0,0,2,0,0,0,2,0,0,1,2,5,0,0,0,0,1,1,0,0,0,0,1,2,0,0,4,2,1,2,0,0,0,7,0,0,0,1))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(c+2*rnorm(100,0,6))
  x = x-mean(x)
  z = t(t(z)-colMeans(z))
  
  xhat<-as.numeric(predict(fit5 ,z, s=fit5$lambda.min))
  ppre<-as.numeric(predict(fit6,data.frame(xhat)))
  list9[[length(list9)+1]]<-mean((y-ppre)^2)
}

ap5<-unlist(list8)

lasso.mse<-as.numeric(unlist(list6))
ols.mse<-as.numeric(unlist(list5))
pv.mse<-as.numeric(unlist(list9))
v<-c(rep("TSLS",100),rep("TSL",100),rep("TSPL",100))
mse<-as.data.frame(cbind(v,c(ols.mse,lasso.mse,pv.mse)))
mse[,2]<- as.numeric(as.character(mse[,2]))
ggplot(mse,aes(x=mse$v,y=mse[,2],group=mse$v))+scale_y_continuous(limits = c(0, 25000))+geom_boxplot(fill=c("darkblue","skyblue","darkgreen"))+labs(x="Model",y="MSE",title="No Association between X and Y")+theme(plot.title = element_text(hjust = 0.5))

tnn<-list4>0.05
tn<-length(tnn[tnn=="TRUE"])
fp<-100-tn
tnn1<-list3>0.05
tn1<-length(tnn1[tnn1=="TRUE"])
fp1<-100-tn1
tnn2<-list8>0.05
tn2<-length(tnn2[tnn2=="TRUE"])
fp2<-100-tn2

mcc(tp,tn,fp,fn)
mcc(tp1,tn1,fp1,fn1)
mcc(tp2,tn2,fp2,fn2)

#ROC plot
app<-c(ap,ap1)
class = rep(c(1,0),c(length(ap1),length(ap)))
ord = order(app)
plot(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),type = 'l',main = 'ROC curve',xlab = 'FPR',ylab = 'TPR',col="darkblue",lwd=4)
abline(0,1)
p = c(ap2,ap3)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="skyblue",type = 'l',lwd=4)
p = c(ap4,ap5)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="darkgreen",type = 'l',lwd=4)
legend('bottomright',c('TSL','TSLS','TSPL'),col = c("darkblue","skyblue","darkgreen"),lty = 1,cex=1,lwd=3)

df<-as.data.frame(cbind(as.numeric(c(ap,ap1,ap2,ap3,ap4,ap5)),as.factor(c(rep("Ass",100),rep("No Ass",100))),(c(rep("TSL",200),rep("TSLS",200),rep("TSPL",200)))))
colnames(df)<-c("p.value","ass","method")
df[,1]<- as.numeric(as.character(df[,1]))
df<-df[c(1:100,200:300,400:500),]
ggplot(df,aes(x=df$method,y=df$p.value,group=df$method))+labs(x="Method",y="p-value",title="Range of p-values from the Simulations")+theme(plot.title = element_text(hjust = 0.5))+geom_boxplot(fill=c("darkblue","skyblue","darkgreen"))

##only more uninformative IVs added

dta1<-data.frame(x=1:40)
for (j in 1:40){
  list1<-list()
  for (i in 1:40){
    tl<-0.5^abs(i-j)
    list1[[length(list1)+1]]<-tl
  }
  dta1<-cbind(dta1,unlist(list1))
}
s<-matrix(c(dta1[,2],dta1[,3],dta1[,4],dta1[,5],dta1[,6],dta1[,7],dta1[,8],dta1[,9],dta1[,10],dta1[,11],dta1[,12],dta1[,13],dta1[,14],dta1[,15],dta1[,16],dta1[,17],dta1[,18],dta1[,19],dta1[,20],dta1[,21],dta1[,22],dta1[,23],dta1[,24],dta1[,25],dta1[,26],dta1[,27],dta1[,28],dta1[,29],dta1[,30],dta1[,31],dta1[,32],dta1[,33],dta1[,34],dta1[,35],dta1[,36],dta1[,37],dta1[,38],dta1[,39],dta1[,40],dta1[,41]),ncol=40,nrow=40)


list2<-list()
list4<-list()
list6<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(100,rep(2,40),s)
  b<-(c(rep(0,12),3,1.5,0,0,2,0,0,0,rep(0,20)))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(x+c+2*rnorm(100,0,6))
  crossval <-  cv.glmnet(y=x,x=z)
  penalty <- crossval$lambda.min
  fit3 <-glmnet( y=x,x=z, alpha = 1, lambda = penalty )
  xhat<-predict(crossval, newx = z, s = penalty)
  xhat<-as.vector(xhat[,1])
  fit4<-lm(y~xhat)
  pp2<-summary(fit4)
  pp3<-1
  if(nrow(coef(summary(fit4)))==2) pp3 = coef(summary(fit4))[2,4]
  tt2<-deviance(fit4)
  list2[[length(list2)+1]]<-tt2
  list4[[length(list4)+1]]<-pp3
  
  z<-mvrnorm(100,rep(2,40),s)
  b<-(c(rep(0,12),3,1.5,0,0,2,0,0,0,rep(0,20)))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(x+c+2*rnorm(100,0,6))
  
  xhat<-as.numeric(predict(fit3 ,z, s=fit3$lambda.min))
  ypre<-as.numeric(predict(fit4,as.data.frame(as.numeric(xhat))))
  list6[[length(list6)+1]]<-mean((y-ypre)^2)
}
ap<-unlist(list4)

list1<-list()
list3<-list()
list5<-list()
for (i in 1:100){
  set.seed(i)
  z<-as.data.frame(mvrnorm(100,rep(2,40),s))
  b<-(c(rep(0,12),3,1.5,0,0,2,0,0,0,rep(0,20)))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(x+c+2*rnorm(100,0,6))
  z<-cbind(z,x)
  fit1<-lm(x~.,z)
  xhat<-fit1$fitted.values
  fit2<-lm(y~xhat)
  pp<-summary(fit2)
  pp1<-pp$coefficients[2,4]
  tt1<-deviance(fit2)
  mse1<-mean(pp$residuals^2)
  list1[[length(list1)+1]]<-tt1
  list3[[length(list3)+1]]<-pp1
  
  z<-as.data.frame(mvrnorm(100,rep(2,40),s))
  b<-(c(rep(0,12),3,1.5,0,0,2,0,0,0,rep(0,20)))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(x+c+2*rnorm(100,0,6))
  z<-cbind(z,x)
  xhat<-as.numeric(predict(fit1,as.data.frame(z[,1:40])))
  ypre<-as.numeric(predict(fit2,as.data.frame(xhat)))
  list5[[length(list5)+1]]<-mean((y-ypre)^2)
}
ap2<-unlist(list3)

list8<-list()
list9<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(100,rep(2,40),s)
  b<-(c(rep(0,12),3,1.5,0,0,2,0,0,0,rep(0,20)))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(x+c+2*rnorm(100,0,6))
  fit5 <-  cv.glmnet(y=x,x=z,keep=TRUE)
  xhat = as.numeric(fit5$fit.preval[,which(fit5$lambda==fit5$lambda.min)[1]])
  fit6 = lm(y~xhat)
  pp4 = 1
  if(nrow(coef(summary(fit6)))==2) pp4 = coef(summary(fit6))[2,4]
  list8[[length(list8)+1]]<-pp4
  
  z<-mvrnorm(100,rep(2,40),s)
  b<-(c(rep(0,12),3,1.5,0,0,2,0,0,0,rep(0,20)))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(x+c+2*rnorm(100,0,6))
  
  xhat<-as.numeric(predict(fit5 ,z, s=fit5$lambda.min))
  ppre<-as.numeric(predict(fit6,data.frame(xhat)))
  list9[[length(list9)+1]]<-mean((y-ppre)^2)
}
ap4<-unlist(list8)

lasso.mse<-as.numeric(unlist(list6))
ols.mse<-as.numeric(unlist(list5))
pv.mse<-as.numeric(unlist(list9))
v<-c(rep("TSLS",100),rep("TSL",100),rep("TSPL",100))
mse<-as.data.frame(cbind(v,c(ols.mse,lasso.mse,pv.mse)))
mse[,2]<- as.numeric(as.character(mse[,2]))
ggplot(mse,aes(x=mse$v,y=mse[,2],group=mse$v))+scale_y_continuous(limits = c(0, 25000))+geom_boxplot(fill=c("darkblue","skyblue","darkgreen"))+labs(x="Model",y="MSE",title="Association between X and Y")+theme(plot.title = element_text(hjust = 0.5))

fnn<-list4>0.05
fn<-length(fnn[fnn=="TRUE"])
tp<-100-fn
fnn1<-list3>0.05
fn1<-length(fnn1[fnn1=="TRUE"])
tp1<-100-fn1
fnn2<-list8>0.05
fn2<-length(fnn2[fnn2=="TRUE"])
tp2<-100-fn2


dta1<-data.frame(x=1:40)
for (j in 1:40){
  list1<-list()
  for (i in 1:40){
    tl<-0.5^abs(i-j)
    list1[[length(list1)+1]]<-tl
  }
  dta1<-cbind(dta1,unlist(list1))
}
s<-matrix(c(dta1[,2],dta1[,3],dta1[,4],dta1[,5],dta1[,6],dta1[,7],dta1[,8],dta1[,9],dta1[,10],dta1[,11],dta1[,12],dta1[,13],dta1[,14],dta1[,15],dta1[,16],dta1[,17],dta1[,18],dta1[,19],dta1[,20],dta1[,21],dta1[,22],dta1[,23],dta1[,24],dta1[,25],dta1[,26],dta1[,27],dta1[,28],dta1[,29],dta1[,30],dta1[,31],dta1[,32],dta1[,33],dta1[,34],dta1[,35],dta1[,36],dta1[,37],dta1[,38],dta1[,39],dta1[,40],dta1[,41]),ncol=40,nrow=40)


list2<-list()
list4<-list()
list6<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(100,rep(2,40),s)
  b<-(c(rep(0,12),3,1.5,0,0,2,0,0,0,rep(0,20)))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(c+2*rnorm(100,0,6))
  crossval <-  cv.glmnet(y=x,x=z)
  penalty <- crossval$lambda.min
  fit3 <-glmnet( y=x,x=z, alpha = 1, lambda = penalty )
  xhat<-predict(crossval, newx = z, s = penalty)
  xhat<-as.vector(xhat[,1])
  fit4<-lm(y~xhat)
  pp2<-summary(fit4)
  pp3<-1
  if(nrow(coef(summary(fit4)))==2) pp3 = coef(summary(fit4))[2,4]
  tt2<-deviance(fit4)
  list2[[length(list2)+1]]<-tt2
  list4[[length(list4)+1]]<-pp3
  
  z<-mvrnorm(100,rep(2,40),s)
  b<-(c(rep(0,12),3,1.5,0,0,2,0,0,0,rep(0,20)))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(c+2*rnorm(100,0,6))
  
  xhat<-as.numeric(predict(fit3 ,z, s=fit3$lambda.min))
  ypre<-as.numeric(predict(fit4,as.data.frame(as.numeric(xhat))))
  list6[[length(list6)+1]]<-mean((y-ypre)^2)
}
ap1<-unlist(list4)

list1<-list()
list3<-list()
list5<-list()
for (i in 1:100){
  set.seed(i)
  z<-as.data.frame(mvrnorm(100,rep(2,40),s))
  b<-(c(rep(0,12),3,1.5,0,0,2,0,0,0,rep(0,20)))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(c+2*rnorm(100,0,6))
  z<-cbind(z,x)
  fit1<-lm(x~.,z)
  xhat<-fit1$fitted.values
  fit2<-lm(y~xhat)
  pp<-summary(fit2)
  pp1<-pp$coefficients[2,4]
  tt1<-deviance(fit2)
  mse1<-mean(pp$residuals^2)
  list1[[length(list1)+1]]<-tt1
  list3[[length(list3)+1]]<-pp1
  
  z<-as.data.frame(mvrnorm(100,rep(2,40),s))
  b<-(c(rep(0,12),3,1.5,0,0,2,0,0,0,rep(0,20)))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(c+2*rnorm(100,0,6))
  z<-cbind(z,x)
  xhat<-as.numeric(predict(fit1,as.data.frame(z[,1:40])))
  ypre<-as.numeric(predict(fit2,as.data.frame(xhat)))
  list5[[length(list5)+1]]<-mean((y-ypre)^2)
}
ap3<-unlist(list3)

list8<-list()
list9<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(100,rep(2,40),s)
  b<-(c(rep(0,12),3,1.5,0,0,2,0,0,0,rep(0,20)))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(c+2*rnorm(100,0,6))
  fit5 <-  cv.glmnet(y=x,x=z,keep=TRUE)
  xhat = as.numeric(fit5$fit.preval[,which(fit5$lambda==fit5$lambda.min)[1]])
  fit6 = lm(y~xhat)
  pp4 = 1
  if(nrow(coef(summary(fit6)))==2) pp4 = coef(summary(fit6))[2,4]
  list8[[length(list8)+1]]<-pp4
  
  z<-mvrnorm(100,rep(2,40),s)
  b<-(c(rep(0,12),3,1.5,0,0,2,0,0,0,rep(0,20)))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(c+2*rnorm(100,0,6))
  
  xhat<-as.numeric(predict(fit5 ,z, s=fit5$lambda.min))
  ppre<-as.numeric(predict(fit6,data.frame(xhat)))
  list9[[length(list9)+1]]<-mean((y-ppre)^2)
}
ap5<-unlist(list8)

lasso.mse<-as.numeric(unlist(list6))
ols.mse<-as.numeric(unlist(list5))
pv.mse<-as.numeric(unlist(list9))
v<-c(rep("TSLS",100),rep("TSL",100),rep("TSPL",100))
mse<-as.data.frame(cbind(v,c(ols.mse,lasso.mse,pv.mse)))
mse[,2]<- as.numeric(as.character(mse[,2]))
ggplot(mse,aes(x=mse$v,y=mse[,2],group=mse$v))+scale_y_continuous(limits = c(0, 25000))+geom_boxplot(fill=c("darkblue","skyblue","darkgreen"))+labs(x="Model",y="MSE",title="No Association between X and Y")+theme(plot.title = element_text(hjust = 0.5))

tnn<-list4>0.05
tn<-length(tnn[tnn=="TRUE"])
fp<-100-tn
tnn1<-list3>0.05
tn1<-length(tnn1[tnn1=="TRUE"])
fp1<-100-tn1
tnn2<-list8>0.05
tn2<-length(tnn2[tnn2=="TRUE"])
fp2<-100-tn2

mcc(tp,tn,fp,fn)
mcc(tp1,tn1,fp1,fn1)
mcc(tp2,tn2,fp2,fn2)

#ROC plot
app<-c(ap,ap1)
class = rep(c(1,0),c(length(ap1),length(ap)))
ord = order(app)
plot(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),type = 'l',main = 'ROC curve',xlab = 'FPR',ylab = 'TPR',col="darkblue",lwd=4)
abline(0,1)
p = c(ap2,ap3)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="skyblue",type = 'l',lwd=4)
p = c(ap4,ap5)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="darkgreen",type = 'l',lwd=4)
legend('bottomright',c('TSL','TSLS','TSPL'),col = c("darkblue","skyblue","darkgreen"),lty = 1,cex=1,lwd=3)


#increase correlation between IVs
s<-diag(1,nrow=40,ncol=40)
s[1,2]<-0.5
s[2,1]<-0.5
s[10,5]<-0.10
s[5,10]<-0.10
s[25,18]<-0.4
s[18,25]<-0.4
s[4,33]<-0.8
s[33,4]<-0.8
s[14,15]<-0.7
s[15,14]<-0.7
s[16,12]<-0.3
s[12,16]<-0.3


list2<-list()
list4<-list()
list6<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(100,rep(2,40),s)
  b<-(c(rep(0,12),3,1.5,0,0,2,0,0,0,rep(0,20)))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(x+c+2*rnorm(100,0,6))
  crossval <-  cv.glmnet(y=x,x=z)
  penalty <- crossval$lambda.min
  fit3 <-glmnet( y=x,x=z, alpha = 1, lambda = penalty )
  xhat<-predict(crossval, newx = z, s = penalty)
  xhat<-as.vector(xhat[,1])
  fit4<-lm(y~xhat)
  pp2<-summary(fit4)
  pp3<-1
  if(nrow(coef(summary(fit4)))==2) pp3 = coef(summary(fit4))[2,4]
  tt2<-deviance(fit4)
  list2[[length(list2)+1]]<-tt2
  list4[[length(list4)+1]]<-pp3
  
  z<-mvrnorm(100,rep(2,40),s)
  b<-(c(rep(0,12),3,1.5,0,0,2,0,0,0,rep(0,20)))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(x+c+2*rnorm(100,0,6))
  
  xhat<-as.numeric(predict(fit3 ,z, s=fit3$lambda.min))
  ypre<-as.numeric(predict(fit4,as.data.frame(as.numeric(xhat))))
  list6[[length(list6)+1]]<-mean((y-ypre)^2)
}
ap<-unlist(list4)

list1<-list()
list3<-list()
list5<-list()
for (i in 1:100){
  set.seed(i)
  z<-as.data.frame(mvrnorm(100,rep(2,40),s))
  b<-(c(rep(0,12),3,1.5,0,0,2,0,0,0,rep(0,20)))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(x+c+2*rnorm(100,0,6))
  z<-cbind(z,x)
  fit1<-lm(x~.,z)
  xhat<-fit1$fitted.values
  fit2<-lm(y~xhat)
  pp<-summary(fit2)
  pp1<-pp$coefficients[2,4]
  tt1<-deviance(fit2)
  mse1<-mean(pp$residuals^2)
  list1[[length(list1)+1]]<-tt1
  list3[[length(list3)+1]]<-pp1
  
  z<-as.data.frame(mvrnorm(100,rep(2,40),s))
  b<-(c(rep(0,12),3,1.5,0,0,2,0,0,0,rep(0,20)))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(x+c+2*rnorm(100,0,6))
  z<-cbind(z,x)
  xhat<-as.numeric(predict(fit1,as.data.frame(z[,1:40])))
  ypre<-as.numeric(predict(fit2,as.data.frame(xhat)))
  list5[[length(list5)+1]]<-mean((y-ypre)^2)
}
ap2<-unlist(list3)

list8<-list()
list9<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(100,rep(2,40),s)
  b<-(c(rep(0,12),3,1.5,0,0,2,0,0,0,rep(0,20)))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(x+c+2*rnorm(100,0,6))
  fit5 <-  cv.glmnet(y=x,x=z,keep=TRUE)
  xhat = as.numeric(fit5$fit.preval[,which(fit5$lambda==fit5$lambda.min)[1]])
  fit6 = lm(y~xhat)
  pp4 = 1
  if(nrow(coef(summary(fit6)))==2) pp4 = coef(summary(fit6))[2,4]
  list8[[length(list8)+1]]<-pp4
  
  z<-mvrnorm(100,rep(2,40),s)
  b<-(c(rep(0,12),3,1.5,0,0,2,0,0,0,rep(0,20)))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(x+c+2*rnorm(100,0,6))
  
  xhat<-as.numeric(predict(fit5 ,z, s=fit5$lambda.min))
  ppre<-as.numeric(predict(fit6,data.frame(xhat)))
  list9[[length(list9)+1]]<-mean((y-ppre)^2)
}
ap4<-unlist(list8)

lasso.mse<-as.numeric(unlist(list6))
ols.mse<-as.numeric(unlist(list5))
pv.mse<-as.numeric(unlist(list9))
v<-c(rep("TSLS",100),rep("TSL",100),rep("TSPL",100))
mse<-as.data.frame(cbind(v,c(ols.mse,lasso.mse,pv.mse)))
mse[,2]<- as.numeric(as.character(mse[,2]))
ggplot(mse,aes(x=mse$v,y=mse[,2],group=mse$v))+scale_y_continuous(limits = c(0, 1000))+geom_boxplot(fill=c("darkblue","skyblue","darkgreen"))+labs(x="Model",y="MSE",title="Association between X and Y")+theme(plot.title = element_text(hjust = 0.5))

fnn<-list4>0.05
fn<-length(fnn[fnn=="TRUE"])
tp<-100-fn
fnn1<-list3>0.05
fn1<-length(fnn1[fnn1=="TRUE"])
tp1<-100-fn1
fnn2<-list8>0.05
fn2<-length(fnn2[fnn2=="TRUE"])
tp2<-100-fn2



list2<-list()
list4<-list()
list6<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(100,rep(2,40),s)
  b<-(c(rep(0,12),3,1.5,0,0,2,0,0,0,rep(0,20)))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(c+2*rnorm(100,0,6))
  crossval <-  cv.glmnet(y=x,x=z)
  penalty <- crossval$lambda.min
  fit3 <-glmnet( y=x,x=z, alpha = 1, lambda = penalty )
  xhat<-predict(crossval, newx = z, s = penalty)
  xhat<-as.vector(xhat[,1])
  fit4<-lm(y~xhat)
  pp2<-summary(fit4)
  pp3<-1
  if(nrow(coef(summary(fit4)))==2) pp3 = coef(summary(fit4))[2,4]
  tt2<-deviance(fit4)
  list2[[length(list2)+1]]<-tt2
  list4[[length(list4)+1]]<-pp3
  
  z<-mvrnorm(100,rep(2,40),s)
  b<-(c(rep(0,12),3,1.5,0,0,2,0,0,0,rep(0,20)))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(c+2*rnorm(100,0,6))
  
  xhat<-as.numeric(predict(fit3 ,z, s=fit3$lambda.min))
  ypre<-as.numeric(predict(fit4,as.data.frame(as.numeric(xhat))))
  list6[[length(list6)+1]]<-mean((y-ypre)^2)
}
ap1<-unlist(list4)

list1<-list()
list3<-list()
list5<-list()
for (i in 1:100){
  set.seed(i)
  z<-as.data.frame(mvrnorm(100,rep(2,40),s))
  b<-(c(rep(0,12),3,1.5,0,0,2,0,0,0,rep(0,20)))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(c+2*rnorm(100,0,6))
  z<-cbind(z,x)
  fit1<-lm(x~.,z)
  xhat<-fit1$fitted.values
  fit2<-lm(y~xhat)
  pp<-summary(fit2)
  pp1<-pp$coefficients[2,4]
  tt1<-deviance(fit2)
  mse1<-mean(pp$residuals^2)
  list1[[length(list1)+1]]<-tt1
  list3[[length(list3)+1]]<-pp1
  
  z<-as.data.frame(mvrnorm(100,rep(2,40),s))
  b<-(c(rep(0,12),3,1.5,0,0,2,0,0,0,rep(0,20)))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(c+2*rnorm(100,0,6))
  z<-cbind(z,x)
  xhat<-as.numeric(predict(fit1,as.data.frame(z[,1:40])))
  ypre<-as.numeric(predict(fit2,as.data.frame(xhat)))
  list5[[length(list5)+1]]<-mean((y-ypre)^2)
}
ap3<-unlist(list3)

list8<-list()
list9<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(100,rep(2,40),s)
  b<-(c(rep(0,12),3,1.5,0,0,2,0,0,0,rep(0,20)))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(c+2*rnorm(100,0,6))
  fit5 <-  cv.glmnet(y=x,x=z,keep=TRUE)
  xhat = as.numeric(fit5$fit.preval[,which(fit5$lambda==fit5$lambda.min)[1]])
  fit6 = lm(y~xhat)
  pp4 = 1
  if(nrow(coef(summary(fit6)))==2) pp4 = coef(summary(fit6))[2,4]
  list8[[length(list8)+1]]<-pp4
  
  z<-mvrnorm(100,rep(2,40),s)
  b<-(c(rep(0,12),3,1.5,0,0,2,0,0,0,rep(0,20)))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(c+2*rnorm(100,0,6))
  
  xhat<-as.numeric(predict(fit5 ,z, s=fit5$lambda.min))
  ppre<-as.numeric(predict(fit6,data.frame(xhat)))
  list9[[length(list9)+1]]<-mean((y-ppre)^2)
}
ap5<-unlist(list8)

lasso.mse<-as.numeric(unlist(list6))
ols.mse<-as.numeric(unlist(list5))
pv.mse<-as.numeric(unlist(list9))
v<-c(rep("TSLS",100),rep("TSL",100),rep("TSPL",100))
mse<-as.data.frame(cbind(v,c(ols.mse,lasso.mse,pv.mse)))
mse[,2]<- as.numeric(as.character(mse[,2]))
ggplot(mse,aes(x=mse$v,y=mse[,2],group=mse$v))+scale_y_continuous(limits = c(0, 1000))+geom_boxplot(fill=c("darkblue","skyblue","darkgreen"))+labs(x="Model",y="MSE",title="No Association between X and Y")+theme(plot.title = element_text(hjust = 0.5))

tnn<-list4>0.05
tn<-length(tnn[tnn=="TRUE"])
fp<-100-tn
tnn1<-list3>0.05
tn1<-length(tnn1[tnn1=="TRUE"])
fp1<-100-tn1
tnn2<-list8>0.05
tn2<-length(tnn2[tnn2=="TRUE"])
fp2<-100-tn2

#make MCC function
mcc<-function(a,b,c,d){
  (a*b-c*d)/sqrt((a+c)*(b+d)*(c+b)*(a+d))
}

mcc(tp,tn,fp,fn)
mcc(tp1,tn1,fp1,fn1)
mcc(tp2,tn2,fp2,fn2)

#ROC plot
app<-c(ap,ap1)
class = rep(c(1,0),c(length(ap1),length(ap)))
ord = order(app)
plot(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),type = 'l',main = 'ROC curve',xlab = 'FPR',ylab = 'TPR',col="darkblue",lwd=4)
abline(0,1)
p = c(ap2,ap3)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="skyblue",type = 'l',lwd=4)
p = c(ap4,ap5)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="darkgreen",type = 'l',lwd=4)
legend('bottomright',c('TSL','TSLS','TSPL'),col = c("darkblue","skyblue","darkgreen"),lty = 1,cex=1,lwd=3)



#neg associatio between some IVs and x
s<-diag(1,nrow=40,ncol=40)
s[1,2]<-0.5
s[2,1]<-0.5
s[10,5]<-0.10
s[5,10]<-0.10
s[25,18]<-0.4
s[18,25]<-0.4
s[4,33]<-0.8
s[33,4]<-0.8
s[14,15]<-0.7
s[15,14]<-0.7
s[16,12]<-0.3
s[12,16]<-0.3
b<-(c(rep(0,12),-3,1.5,0,0,-2,0,0,0,rep(0,20)))


list2<-list()
list4<-list()
list6<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(100,rep(2,40),s)
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(x+c+2*rnorm(100,0,6))
  crossval <-  cv.glmnet(y=x,x=z)
  penalty <- crossval$lambda.min
  fit3 <-glmnet( y=x,x=z, alpha = 1, lambda = penalty )
  xhat<-predict(crossval, newx = z, s = penalty)
  xhat<-as.vector(xhat[,1])
  fit4<-lm(y~xhat)
  pp2<-summary(fit4)
  pp3<-1
  if(nrow(coef(summary(fit4)))==2) pp3 = coef(summary(fit4))[2,4]
  tt2<-deviance(fit4)
  list2[[length(list2)+1]]<-tt2
  list4[[length(list4)+1]]<-pp3
  
  z<-mvrnorm(100,rep(2,40),s)
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(x+c+2*rnorm(100,0,6))
  
  xhat<-as.numeric(predict(fit3 ,z, s=fit3$lambda.min))
  ypre<-as.numeric(predict(fit4,as.data.frame(as.numeric(xhat))))
  list6[[length(list6)+1]]<-mean((y-ypre)^2)
}
ap<-unlist(list4)

list1<-list()
list3<-list()
list5<-list()
for (i in 1:100){
  set.seed(i)
  z<-as.data.frame(mvrnorm(100,rep(2,40),s))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(x+c+2*rnorm(100,0,6))
  z<-cbind(z,x)
  fit1<-lm(x~.,z)
  xhat<-fit1$fitted.values
  fit2<-lm(y~xhat)
  pp<-summary(fit2)
  pp1<-pp$coefficients[2,4]
  tt1<-deviance(fit2)
  mse1<-mean(pp$residuals^2)
  list1[[length(list1)+1]]<-tt1
  list3[[length(list3)+1]]<-pp1
  
  z<-as.data.frame(mvrnorm(100,rep(2,40),s))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(x+c+2*rnorm(100,0,6))
  z<-cbind(z,x)
  xhat<-as.numeric(predict(fit1,as.data.frame(z[,1:40])))
  ypre<-as.numeric(predict(fit2,as.data.frame(xhat)))
  list5[[length(list5)+1]]<-mean((y-ypre)^2)
}
ap2<-unlist(list3)

list8<-list()
list9<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(100,rep(2,40),s)
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(x+c+2*rnorm(100,0,6))
  fit5 <-  cv.glmnet(y=x,x=z,keep=TRUE)
  xhat = as.numeric(fit5$fit.preval[,which(fit5$lambda==fit5$lambda.min)[1]])
  fit6 = lm(y~xhat)
  pp4 = 1
  if(nrow(coef(summary(fit6)))==2) pp4 = coef(summary(fit6))[2,4]
  list8[[length(list8)+1]]<-pp4
  
  z<-mvrnorm(100,rep(2,40),s)
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(x+c+2*rnorm(100,0,6))
  
  xhat<-as.numeric(predict(fit5 ,z, s=fit5$lambda.min))
  ppre<-as.numeric(predict(fit6,data.frame(xhat)))
  list9[[length(list9)+1]]<-mean((y-ppre)^2)
}
ap4<-unlist(list8)

lasso.mse<-as.numeric(unlist(list6))
ols.mse<-as.numeric(unlist(list5))
pv.mse<-as.numeric(unlist(list9))
v<-c(rep("TSLS",100),rep("TSL",100),rep("TSPL",100))
mse<-as.data.frame(cbind(v,c(ols.mse,lasso.mse,pv.mse)))
mse[,2]<- as.numeric(as.character(mse[,2]))
ggplot(mse,aes(x=mse$v,y=mse[,2],group=mse$v))+scale_y_continuous(limits = c(0, 1000))+geom_boxplot(fill=c("darkblue","skyblue","darkgreen"))+labs(x="Model",y="MSE",title="Association between X and Y")+theme(plot.title = element_text(hjust = 0.5))

fnn<-list4>0.05
fn<-length(fnn[fnn=="TRUE"])
tp<-100-fn
fnn1<-list3>0.05
fn1<-length(fnn1[fnn1=="TRUE"])
tp1<-100-fn1
fnn2<-list8>0.05
fn2<-length(fnn2[fnn2=="TRUE"])
tp2<-100-fn2



list2<-list()
list4<-list()
list6<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(100,rep(2,40),s)
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(c+2*rnorm(100,0,6))
  crossval <-  cv.glmnet(y=x,x=z)
  penalty <- crossval$lambda.min
  fit3 <-glmnet( y=x,x=z, alpha = 1, lambda = penalty )
  xhat<-predict(crossval, newx = z, s = penalty)
  xhat<-as.vector(xhat[,1])
  fit4<-lm(y~xhat)
  pp2<-summary(fit4)
  pp3<-1
  if(nrow(coef(summary(fit4)))==2) pp3 = coef(summary(fit4))[2,4]
  tt2<-deviance(fit4)
  list2[[length(list2)+1]]<-tt2
  list4[[length(list4)+1]]<-pp3
  
  z<-mvrnorm(100,rep(2,40),s)
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(c+2*rnorm(100,0,6))
  
  xhat<-as.numeric(predict(fit3 ,z, s=fit3$lambda.min))
  ypre<-as.numeric(predict(fit4,as.data.frame(as.numeric(xhat))))
  list6[[length(list6)+1]]<-mean((y-ypre)^2)
}
ap1<-unlist(list4)

list1<-list()
list3<-list()
list5<-list()
for (i in 1:100){
  set.seed(i)
  z<-as.data.frame(mvrnorm(100,rep(2,40),s))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(c+2*rnorm(100,0,6))
  z<-cbind(z,x)
  fit1<-lm(x~.,z)
  xhat<-fit1$fitted.values
  fit2<-lm(y~xhat)
  pp<-summary(fit2)
  pp1<-pp$coefficients[2,4]
  tt1<-deviance(fit2)
  mse1<-mean(pp$residuals^2)
  list1[[length(list1)+1]]<-tt1
  list3[[length(list3)+1]]<-pp1
  
  z<-as.data.frame(mvrnorm(100,rep(2,40),s))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(c+2*rnorm(100,0,6))
  z<-cbind(z,x)
  xhat<-as.numeric(predict(fit1,as.data.frame(z[,1:40])))
  ypre<-as.numeric(predict(fit2,as.data.frame(xhat)))
  list5[[length(list5)+1]]<-mean((y-ypre)^2)
}
ap3<-unlist(list3)

list8<-list()
list9<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(100,rep(2,40),s)
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(c+2*rnorm(100,0,6))
  fit5 <-  cv.glmnet(y=x,x=z,keep=TRUE)
  xhat = as.numeric(fit5$fit.preval[,which(fit5$lambda==fit5$lambda.min)[1]])
  fit6 = lm(y~xhat)
  pp4 = 1
  if(nrow(coef(summary(fit6)))==2) pp4 = coef(summary(fit6))[2,4]
  list8[[length(list8)+1]]<-pp4
  
  z<-mvrnorm(100,rep(2,40),s)
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(c+2*rnorm(100,0,6))
  
  xhat<-as.numeric(predict(fit5 ,z, s=fit5$lambda.min))
  ppre<-as.numeric(predict(fit6,data.frame(xhat)))
  list9[[length(list9)+1]]<-mean((y-ppre)^2)
}
ap5<-unlist(list8)

lasso.mse<-as.numeric(unlist(list6))
ols.mse<-as.numeric(unlist(list5))
pv.mse<-as.numeric(unlist(list9))
v<-c(rep("TSLS",100),rep("TSL",100),rep("TSPL",100))
mse<-as.data.frame(cbind(v,c(ols.mse,lasso.mse,pv.mse)))
mse[,2]<- as.numeric(as.character(mse[,2]))
ggplot(mse,aes(x=mse$v,y=mse[,2],group=mse$v))+scale_y_continuous(limits = c(0, 1000))+geom_boxplot(fill=c("darkblue","skyblue","darkgreen"))+labs(x="Model",y="MSE",title="No Association between X and Y")+theme(plot.title = element_text(hjust = 0.5))

tnn<-list4>0.05
tn<-length(tnn[tnn=="TRUE"])
fp<-100-tn
tnn1<-list3>0.05
tn1<-length(tnn1[tnn1=="TRUE"])
fp1<-100-tn1
tnn2<-list8>0.05
tn2<-length(tnn2[tnn2=="TRUE"])
fp2<-100-tn2

#make MCC function
mcc<-function(a,b,c,d){
  (a*b-c*d)/sqrt((a+c)*(b+d)*(c+b)*(a+d))
}

mcc(tp,tn,fp,fn)
mcc(tp1,tn1,fp1,fn1)
mcc(tp2,tn2,fp2,fn2)

#ROC plot
app<-c(ap,ap1)
class = rep(c(1,0),c(length(ap1),length(ap)))
ord = order(app)
plot(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),type = 'l',main = 'ROC curve',xlab = 'FPR',ylab = 'TPR',col="darkblue",lwd=4)
abline(0,1)
p = c(ap2,ap3)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="skyblue",type = 'l',lwd=4)
p = c(ap4,ap5)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="darkgreen",type = 'l',lwd=4)
legend('bottomright',c('TSL','TSLS','TSPL'),col = c("darkblue","skyblue","darkgreen"),lty = 1,cex=1,lwd=3)

##increasing strength of confounders and more noise
s<-diag(1,nrow=40,ncol=40)
s[1,2]<-0.5
s[2,1]<-0.5
s[10,5]<-0.10
s[5,10]<-0.10
s[25,18]<-0.4
s[18,25]<-0.4
s[4,33]<-0.8
s[33,4]<-0.8
s[14,15]<-0.7
s[15,14]<-0.7
s[16,12]<-0.3
s[12,16]<-0.3


list2<-list()
list4<-list()
list6<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(100,rep(2,40),s)
  b<-(c(rep(0,12),3,1.5,0,0,2,0,0,0,rep(0,20)))
  c<-rnorm(100,0,4)
  x<-as.vector(b%*%t(z)+8*c+10*rnorm(100,0,20))
  y<-as.vector(x+3*c+5*rnorm(100,0,17))
  crossval <-  cv.glmnet(y=x,x=z)
  penalty <- crossval$lambda.min
  fit3 <-glmnet( y=x,x=z, alpha = 1, lambda = penalty )
  xhat<-predict(crossval, newx = z, s = penalty)
  xhat<-as.vector(xhat[,1])
  fit4<-lm(y~xhat)
  pp2<-summary(fit4)
  pp3<-1
  if(nrow(coef(summary(fit4)))==2) pp3 = coef(summary(fit4))[2,4]
  tt2<-deviance(fit4)
  list2[[length(list2)+1]]<-tt2
  list4[[length(list4)+1]]<-pp3
  
  z<-mvrnorm(100,rep(2,40),s)
  b<-(c(rep(0,12),3,1.5,0,0,2,0,0,0,rep(0,20)))
  c<-rnorm(100,0,4)
  x<-as.vector(b%*%t(z)+8*c+10*rnorm(100,0,20))
  y<-as.vector(x+3*c+5*rnorm(100,0,17))
  
  xhat<-as.numeric(predict(fit3 ,z, s=fit3$lambda.min))
  ypre<-as.numeric(predict(fit4,as.data.frame(as.numeric(xhat))))
  list6[[length(list6)+1]]<-mean((y-ypre)^2)
}
ap<-unlist(list4)

list1<-list()
list3<-list()
list5<-list()
for (i in 1:100){
  set.seed(i)
  z<-as.data.frame(mvrnorm(100,rep(2,40),s))
  b<-(c(rep(0,12),3,1.5,0,0,2,0,0,0,rep(0,20)))
  c<-rnorm(100,0,4)
  x<-as.vector(b%*%t(z)+8*c+10*rnorm(100,0,20))
  y<-as.vector(x+3*c+5*rnorm(100,0,17))
  z<-cbind(z,x)
  fit1<-lm(x~.,z)
  xhat<-fit1$fitted.values
  fit2<-lm(y~xhat)
  pp<-summary(fit2)
  pp1<-pp$coefficients[2,4]
  tt1<-deviance(fit2)
  mse1<-mean(pp$residuals^2)
  list1[[length(list1)+1]]<-tt1
  list3[[length(list3)+1]]<-pp1
  
  z<-as.data.frame(mvrnorm(100,rep(2,40),s))
  b<-(c(rep(0,12),3,1.5,0,0,2,0,0,0,rep(0,20)))
  c<-rnorm(100,0,4)
  x<-as.vector(b%*%t(z)+8*c+10*rnorm(100,0,20))
  y<-as.vector(x+3*c+5*rnorm(100,0,17))
  z<-cbind(z,x)
  xhat<-as.numeric(predict(fit1,as.data.frame(z[,1:40])))
  ypre<-as.numeric(predict(fit2,as.data.frame(xhat)))
  list5[[length(list5)+1]]<-mean((y-ypre)^2)
}
ap2<-unlist(list3)

list8<-list()
list9<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(100,rep(2,40),s)
  b<-(c(rep(0,12),3,1.5,0,0,2,0,0,0,rep(0,20)))
  c<-rnorm(100,0,4)
  x<-as.vector(b%*%t(z)+8*c+10*rnorm(100,0,20))
  y<-as.vector(x+3*c+5*rnorm(100,0,17))
  fit5 <-  cv.glmnet(y=x,x=z,keep=TRUE)
  xhat = as.numeric(fit5$fit.preval[,which(fit5$lambda==fit5$lambda.min)[1]])
  fit6 = lm(y~xhat)
  pp4 = 1
  if(nrow(coef(summary(fit6)))==2) pp4 = coef(summary(fit6))[2,4]
  list8[[length(list8)+1]]<-pp4
  
  z<-mvrnorm(100,rep(2,40),s)
  b<-(c(rep(0,12),3,1.5,0,0,2,0,0,0,rep(0,20)))
  c<-rnorm(100,0,4)
  x<-as.vector(b%*%t(z)+8*c+10*rnorm(100,0,20))
  y<-as.vector(x+3*c+5*rnorm(100,0,17))
  
  xhat<-as.numeric(predict(fit5 ,z, s=fit5$lambda.min))
  ppre<-as.numeric(predict(fit6,data.frame(xhat)))
  list9[[length(list9)+1]]<-mean((y-ppre)^2)
}
ap4<-unlist(list8)

lasso.mse<-as.numeric(unlist(list6))
ols.mse<-as.numeric(unlist(list5))
pv.mse<-as.numeric(unlist(list9))
v<-c(rep("TSLS",100),rep("TSL",100),rep("TSPL",100))
mse<-as.data.frame(cbind(v,c(ols.mse,lasso.mse,pv.mse)))
mse[,2]<- as.numeric(as.character(mse[,2]))
ggplot(mse,aes(x=mse$v,y=mse[,2],group=mse$v))+scale_y_continuous(limits = c(0, 1000))+geom_boxplot(fill=c("darkblue","skyblue","darkgreen"))+labs(x="Model",y="MSE",title="Association between X and Y")+theme(plot.title = element_text(hjust = 0.5))

fnn<-list4>0.05
fn<-length(fnn[fnn=="TRUE"])
tp<-100-fn
fnn1<-list3>0.05
fn1<-length(fnn1[fnn1=="TRUE"])
tp1<-100-fn1
fnn2<-list8>0.05
fn2<-length(fnn2[fnn2=="TRUE"])
tp2<-100-fn2



list2<-list()
list4<-list()
list6<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(100,rep(2,40),s)
  b<-(c(rep(0,12),3,1.5,0,0,2,0,0,0,rep(0,20)))
  c<-rnorm(100,0,4)
  x<-as.vector(b%*%t(z)+8*c+10*rnorm(100,0,20))
  y<-as.vector(3*c+5*rnorm(100,0,17))
  crossval <-  cv.glmnet(y=x,x=z)
  penalty <- crossval$lambda.min
  fit3 <-glmnet( y=x,x=z, alpha = 1, lambda = penalty )
  xhat<-predict(crossval, newx = z, s = penalty)
  xhat<-as.vector(xhat[,1])
  fit4<-lm(y~xhat)
  pp2<-summary(fit4)
  pp3<-1
  if(nrow(coef(summary(fit4)))==2) pp3 = coef(summary(fit4))[2,4]
  tt2<-deviance(fit4)
  list2[[length(list2)+1]]<-tt2
  list4[[length(list4)+1]]<-pp3
  
  z<-mvrnorm(100,rep(2,40),s)
  b<-(c(rep(0,12),3,1.5,0,0,2,0,0,0,rep(0,20)))
  c<-rnorm(100,0,4)
  x<-as.vector(b%*%t(z)+8*c+10*rnorm(100,0,20))
  y<-as.vector(3*c+5*rnorm(100,0,17))
  
  xhat<-as.numeric(predict(fit3 ,z, s=fit3$lambda.min))
  ypre<-as.numeric(predict(fit4,as.data.frame(as.numeric(xhat))))
  list6[[length(list6)+1]]<-mean((y-ypre)^2)
}
ap1<-unlist(list4)

list1<-list()
list3<-list()
list5<-list()
for (i in 1:100){
  set.seed(i)
  z<-as.data.frame(mvrnorm(100,rep(2,40),s))
  b<-(c(rep(0,12),3,1.5,0,0,2,0,0,0,rep(0,20)))
  c<-rnorm(100,0,4)
  x<-as.vector(b%*%t(z)+8*c+10*rnorm(100,0,20))
  y<-as.vector(3*c+5*rnorm(100,0,17))
  z<-cbind(z,x)
  fit1<-lm(x~.,z)
  xhat<-fit1$fitted.values
  fit2<-lm(y~xhat)
  pp<-summary(fit2)
  pp1<-pp$coefficients[2,4]
  tt1<-deviance(fit2)
  mse1<-mean(pp$residuals^2)
  list1[[length(list1)+1]]<-tt1
  list3[[length(list3)+1]]<-pp1
  
  z<-as.data.frame(mvrnorm(100,rep(2,40),s))
  b<-(c(rep(0,12),3,1.5,0,0,2,0,0,0,rep(0,20)))
  c<-rnorm(100,0,4)
  x<-as.vector(b%*%t(z)+8*c+10*rnorm(100,0,20))
  y<-as.vector(3*c+5*rnorm(100,0,17))
  z<-cbind(z,x)
  xhat<-as.numeric(predict(fit1,as.data.frame(z[,1:40])))
  ypre<-as.numeric(predict(fit2,as.data.frame(xhat)))
  list5[[length(list5)+1]]<-mean((y-ypre)^2)
}
ap3<-unlist(list3)

list8<-list()
list9<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(100,rep(2,40),s)
  b<-(c(rep(0,12),3,1.5,0,0,2,0,0,0,rep(0,20)))
  c<-rnorm(100,0,4)
  x<-as.vector(b%*%t(z)+8*c+10*rnorm(100,0,20))
  y<-as.vector(3*c+5*rnorm(100,0,17))
  fit5 <-  cv.glmnet(y=x,x=z,keep=TRUE)
  xhat = as.numeric(fit5$fit.preval[,which(fit5$lambda==fit5$lambda.min)[1]])
  fit6 = lm(y~xhat)
  pp4 = 1
  if(nrow(coef(summary(fit6)))==2) pp4 = coef(summary(fit6))[2,4]
  list8[[length(list8)+1]]<-pp4
  
  z<-mvrnorm(100,rep(2,40),s)
  b<-(c(rep(0,12),3,1.5,0,0,2,0,0,0,rep(0,20)))
  c<-rnorm(100,0,4)
  x<-as.vector(b%*%t(z)+8*c+10*rnorm(100,0,20))
  y<-as.vector(3*c+5*rnorm(100,0,17))
  
  xhat<-as.numeric(predict(fit5 ,z, s=fit5$lambda.min))
  ppre<-as.numeric(predict(fit6,data.frame(xhat)))
  list9[[length(list9)+1]]<-mean((y-ppre)^2)
}
ap5<-unlist(list8)

lasso.mse<-as.numeric(unlist(list6))
ols.mse<-as.numeric(unlist(list5))
pv.mse<-as.numeric(unlist(list9))
v<-c(rep("TSLS",100),rep("TSL",100),rep("TSPL",100))
mse<-as.data.frame(cbind(v,c(ols.mse,lasso.mse,pv.mse)))
mse[,2]<- as.numeric(as.character(mse[,2]))
ggplot(mse,aes(x=mse$v,y=mse[,2],group=mse$v))+scale_y_continuous(limits = c(0, 1000))+geom_boxplot(fill=c("darkblue","skyblue","darkgreen"))+labs(x="Model",y="MSE",title="No Association between X and Y")+theme(plot.title = element_text(hjust = 0.5))

tnn<-list4>0.05
tn<-length(tnn[tnn=="TRUE"])
fp<-100-tn
tnn1<-list3>0.05
tn1<-length(tnn1[tnn1=="TRUE"])
fp1<-100-tn1
tnn2<-list8>0.05
tn2<-length(tnn2[tnn2=="TRUE"])
fp2<-100-tn2

#make MCC function
mcc<-function(a,b,c,d){
  (a*b-c*d)/sqrt((a+c)*(b+d)*(c+b)*(a+d))
}

mcc(tp,tn,fp,fn)
mcc(tp1,tn1,fp1,fn1)
mcc(tp2,tn2,fp2,fn2)

#ROC plot
app<-c(ap,ap1)
class = rep(c(1,0),c(length(ap1),length(ap)))
ord = order(app)
plot(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),type = 'l',main = 'ROC curve',xlab = 'FPR',ylab = 'TPR',col="darkblue",lwd=4)
abline(0,1)
p = c(ap2,ap3)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="skyblue",type = 'l',lwd=4)
p = c(ap4,ap5)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="darkgreen",type = 'l',lwd=4)
legend('bottomright',c('TSL','TSLS','TSPL'),col = c("darkblue","skyblue","darkgreen"),lty = 1,cex=1,lwd=3)



#extra confounders sim1-showing lasso better
dta1<-data.frame(x=1:8)
for (j in 1:8){
  list1<-list()
  for (i in 1:8){
    tl<-0.5^abs(i-j)
    list1[[length(list1)+1]]<-tl
  }
  dta1<-cbind(dta1,unlist(list1))
}
s<-matrix(c(dta1[,2],dta1[,3],dta1[,4],dta1[,5],dta1[,6],dta1[,7],dta1[,8],dta1[,9]),ncol=8,nrow=8)


list2<-list()
list4<-list()
list6<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-b%*%t(z)+10*c+2*rnorm(20)
  y<-as.vector(x+5*c+2*rnorm(20,0,2))
  crossval <-  cv.glmnet(y=x,x=z)
  penalty <- crossval$lambda.min
  fit3 <-glmnet( y=x,x=z, alpha = 1, lambda = penalty )
  xhat<-predict(crossval, newx = z, s = penalty)
  xhat<-as.vector(xhat[,1])
  fit4<-lm(y~xhat)
  pp2<-summary(fit4)
  pp3<-1
  if(nrow(coef(summary(fit4)))==2) pp3 = coef(summary(fit4))[2,4]
  tt2<-deviance(fit4)
  list2[[length(list2)+1]]<-tt2
  list4[[length(list4)+1]]<-pp3
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-b%*%t(z)+10*c+2*rnorm(20)
  y<-as.vector(x+5*c+2*rnorm(20,0,2))
  
  xhat<-as.numeric(predict(fit3 ,z, s=fit3$lambda.min))
  ypre<-as.numeric(predict(fit4,as.data.frame(as.numeric(xhat))))
  list6[[length(list6)+1]]<-mean((y-ypre)^2)
}

ap<-unlist(list4)

list1<-list()
list3<-list()
list5<-list()
for (i in 1:100){
  set.seed(i)
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+10*c+2*rnorm(20))
  y<-as.vector(x+5*c+2*rnorm(20,0,2))
  z<-cbind(z,x)
  fit1<-lm(x~.,z)
  xhat<-fit1$fitted.values
  fit2<-lm(y~xhat)
  pp<-summary(fit2)
  pp1<-pp$coefficients[2,4]
  tt1<-deviance(fit2)
  mse1<-mean(pp$residuals^2)
  list1[[length(list1)+1]]<-tt1
  list3[[length(list3)+1]]<-pp1
  
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-b%*%t(z)+10*c+2*rnorm(20)
  y<-as.vector(x+5*c+2*rnorm(20,0,2))
  z<-cbind(z,x)
  xhat<-as.numeric(predict(fit1,as.data.frame(z[,1:8])))
  ypre<-as.numeric(predict(fit2,as.data.frame(xhat)))
  list5[[length(list5)+1]]<-mean((y-ypre)^2)
}
ap2<-unlist(list3)


list8<-list()
list9<-list()
list10<-list()
for (i in 1:100){
  set.seed(i)
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+10*c+2*rnorm(20))
  y<-as.vector(x+5*c+2*rnorm(20,0,2))
  x = x-mean(x)
  z = t(t(z)-colMeans(z))
  fit5 <-  cv.glmnet(z,x,keep=TRUE,nfolds = length(x),intercept=FALSE)
  lambda = fit5$lambda.min
  xhat = as.numeric(fit5$fit.preval[,which(fit5$lambda==lambda)[1]])
  fit6 = lm(y~xhat)
  rp<-fit6$coefficients[2]
  list10[[length(list10)+1]]<-rp
  pp4 = 1
  if(nrow(coef(summary(fit6)))==2) pp4 = coef(summary(fit6))[2,4]
  list8[[length(list8)+1]]<-pp4
  
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+10*c+2*rnorm(20))
  y<-as.vector(x+5*c+2*rnorm(20,0,2))
  x = x-mean(x)
  z = t(t(z)-colMeans(z))
  
  xhat<-as.numeric(predict(fit5 ,z, s=fit5$lambda.min))
  ppre<-as.numeric(predict(fit6,data.frame(xhat)))
  list9[[length(list9)+1]]<-mean((y-ppre)^2)
}

ap4<-unlist(list8)

lasso.mse<-as.numeric(unlist(list6))
ols.mse<-as.numeric(unlist(list5))
pv.mse<-as.numeric(unlist(list9))
v<-c(rep("TSLS",100),rep("TSL",100),rep("TSPL",100))
mse<-as.data.frame(cbind(v,c(ols.mse,lasso.mse,pv.mse)))
mse[,2]<- as.numeric(as.character(mse[,2]))
ggplot(mse,aes(x=mse$v,y=mse[,2],group=mse$v))+scale_y_continuous(limits = c(0, 6500))+geom_boxplot(fill=c("darkblue","skyblue","darkgreen"))+labs(x="Model",y="MSE",title="Association between X and Y")+theme(plot.title = element_text(hjust = 0.5))

fnn<-list4>0.05
fn<-length(fnn[fnn=="TRUE"])
tp<-100-fn
fnn1<-list3>0.05
fn1<-length(fnn1[fnn1=="TRUE"])
tp1<-100-fn1
fnn2<-list8>0.05
fn2<-length(fnn2[fnn2=="TRUE"])
tp2<-100-fn2

##no association
set.seed(5)
dta1<-data.frame(x=1:8)
for (j in 1:8){
  list1<-list()
  for (i in 1:8){
    tl<-0.5^abs(i-j)
    list1[[length(list1)+1]]<-tl
  }
  dta1<-cbind(dta1,unlist(list1))
}
s<-matrix(c(dta1[,2],dta1[,3],dta1[,4],dta1[,5],dta1[,6],dta1[,7],dta1[,8],dta1[,9]),ncol=8,nrow=8)


list2<-list()
list4<-list()
list6<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20)
  x<-b%*%t(z)+10*c+2*rnorm(20)
  y<-as.vector(5*c+2*rnorm(20,0,2))
  crossval <-  cv.glmnet(y=x,x=z)
  penalty <- crossval$lambda.min
  fit3 <-glmnet( y=x,x=z, alpha = 1, lambda = penalty )
  xhat<-predict(crossval, newx = z, s = penalty)
  xhat<-as.vector(xhat[,1])
  fit4<-lm(y~xhat)
  pp2<-summary(fit4)
  pp3<-1
  if(nrow(coef(summary(fit4)))==2) pp3 = coef(summary(fit4))[2,4]
  tt2<-deviance(fit4)
  list2[[length(list2)+1]]<-tt2
  list4[[length(list4)+1]]<-pp3
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-b%*%t(z)+10*c+2*rnorm(20)
  y<-as.vector(5*c+2*rnorm(20,0,2))
  
  xhat<-as.numeric(predict(fit3 ,z, s=fit3$lambda.min))
  ypre<-as.numeric(predict(fit4,as.data.frame(as.numeric(xhat))))
  list6[[length(list6)+1]]<-mean((y-ypre)^2)
}
ap1<-unlist(list4)

list1<-list()
list3<-list()
list5<-list()
for (i in 1:100){
  set.seed(i)
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+10*c+2*rnorm(20))
  y<-as.vector(5*c+2*rnorm(20,0,2))
  z<-cbind(z,x)
  fit1<-lm(x~.,z)
  xhat<-fit1$fitted.values
  fit2<-lm(y~xhat)
  pp<-summary(fit2)
  pp1<-pp$coefficients[2,4]
  tt1<-deviance(fit2)
  mse1<-mean(pp$residuals^2)
  list1[[length(list1)+1]]<-tt1
  list3[[length(list3)+1]]<-pp1
  
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-b%*%t(z)+10*c+2*rnorm(20)
  y<-as.vector(5*c+2*rnorm(20,0,2))
  z<-cbind(z,x)
  xhat<-as.numeric(predict(fit1,as.data.frame(z[,1:8])))
  ypre<-as.numeric(predict(fit2,as.data.frame(xhat)))
  list5[[length(list5)+1]]<-mean((y-ypre)^2)
}
ap3<-unlist(list3)


list8<-list()
list9<-list()
list10<-list()
for (i in 1:100){
  set.seed(i)
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+10*c+2*rnorm(20))
  y<-as.vector(5*c+2*rnorm(20,0,2))
  x = x-mean(x)
  z = t(t(z)-colMeans(z))
  fit5 <-  cv.glmnet(z,x,keep=TRUE,nfolds = length(x),intercept=FALSE)
  lambda = fit5$lambda.min
  xhat = as.numeric(fit5$fit.preval[,which(fit5$lambda==lambda)[1]])
  fit6 = lm(y~xhat)
  rp<-fit6$coefficients[2]
  list10[[length(list10)+1]]<-rp
  pp4 = 1
  if(nrow(coef(summary(fit6)))==2) pp4 = coef(summary(fit6))[2,4]
  list8[[length(list8)+1]]<-pp4
  
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+10*c+2*rnorm(20))
  y<-as.vector(5*c+2*rnorm(20,0,2))
  x = x-mean(x)
  z = t(t(z)-colMeans(z))
  
  xhat<-as.numeric(predict(fit5 ,z, s=fit5$lambda.min))
  ppre<-as.numeric(predict(fit6,data.frame(xhat)))
  list9[[length(list9)+1]]<-mean((y-ppre)^2)
}

ap5<-unlist(list8)

lasso.mse<-as.numeric(unlist(list6))
ols.mse<-as.numeric(unlist(list5))
pv.mse<-as.numeric(unlist(list9))
v<-c(rep("TSLS",100),rep("TSL",100),rep("TSPL",100))
mse<-as.data.frame(cbind(v,c(ols.mse,lasso.mse,pv.mse)))
mse[,2]<- as.numeric(as.character(mse[,2]))
ggplot(mse,aes(x=mse$v,y=mse[,2],group=mse$v))+scale_y_continuous(limits = c(0, 6500))+geom_boxplot(fill=c("darkblue","skyblue","darkgreen"))+labs(x="Model",y="MSE",title="No Association between X and Y")+theme(plot.title = element_text(hjust = 0.5))

tnn<-list4>0.05
tn<-length(tnn[tnn=="TRUE"])
fp<-100-tn
tnn1<-list3>0.05
tn1<-length(tnn1[tnn1=="TRUE"])
fp1<-100-tn1
tnn2<-list8>0.05
tn2<-length(tnn2[tnn2=="TRUE"])
fp2<-100-tn2

mcc(tp,tn,fp,fn)
mcc(tp1,tn1,fp1,fn1)
mcc(tp2,tn2,fp2,fn2)

app<-c(ap,ap1)
class = rep(c(1,0),c(length(ap1),length(ap)))
ord = order(app)
plot(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),type = 'l',main = 'ROC curve',xlab = 'FPR',ylab = 'TPR',col="darkblue",lwd=4)
abline(0,1)
p = c(ap2,ap3)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="skyblue",type = 'l',lwd=4)
p = c(ap4,ap5)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="darkgreen",type = 'l',lwd=4)
legend('bottomright',c('TSL','TSLS','TSPL'),col = c("darkblue","skyblue","darkgreen"),lty = 1,cex=1,lwd=3)

df<-as.data.frame(cbind(as.numeric(c(ap,ap1,ap2,ap3,ap4,ap5)),as.factor(c(rep("Ass",100),rep("No Ass",100))),(c(rep("TSL",200),rep("TSLS",200),rep("TSPL",200)))))
colnames(df)<-c("p.value","ass","method")
df[,1]<- as.numeric(as.character(df[,1]))
df<-df[c(1:100,200:300,400:500),]
ggplot(df,aes(x=df$method,y=df$p.value,group=df$method))+labs(x="Method",y="p-value",title="Range of p-values from the Simulations")+theme(plot.title = element_text(hjust = 0.5))+geom_boxplot(fill=c("darkblue","skyblue","darkgreen"))


#neg association x and y - preval better
set.seed(5)
dta1<-data.frame(x=1:8)
for (j in 1:8){
  list1<-list()
  for (i in 1:8){
    tl<-0.5^abs(i-j)
    list1[[length(list1)+1]]<-tl
  }
  dta1<-cbind(dta1,unlist(list1))
}
s<-matrix(c(dta1[,2],dta1[,3],dta1[,4],dta1[,5],dta1[,6],dta1[,7],dta1[,8],dta1[,9]),ncol=8,nrow=8)


list2<-list()
list4<-list()
list6<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(1,1)
  x<-g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20)
  y<-as.vector(-0.25*x+c+3*rnorm(20))
  crossval <-  cv.glmnet(y=x,x=z)
  penalty <- crossval$lambda.min
  fit3 <-glmnet( y=x,x=z, alpha = 1, lambda = penalty )
  xhat<-predict(crossval, newx = z, s = penalty)
  xhat<-as.vector(xhat[,1])
  fit4<-lm(y~xhat)
  pp2<-summary(fit4)
  pp3<-1
  if(nrow(coef(summary(fit4)))==2) pp3 = coef(summary(fit4))[2,4]
  tt2<-deviance(fit4)
  list2[[length(list2)+1]]<-tt2
  list4[[length(list4)+1]]<-pp3
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(1,1)
  x<-g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20)
  y<-as.vector(-0.25*x+c+3*rnorm(20))
  
  xhat<-as.numeric(predict(fit3 ,z, s=fit3$lambda.min))
  ypre<-as.numeric(predict(fit4,as.data.frame(as.numeric(xhat))))
  list6[[length(list6)+1]]<-mean((y-ypre)^2)
}
ap<-unlist(list4)

list1<-list()
list3<-list()
list5<-list()
for (i in 1:100){
  set.seed(i)
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(1,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(-0.25*x+c+3*rnorm(20))
  z<-cbind(z,x)
  fit1<-lm(x~.,z)
  xhat<-fit1$fitted.values
  fit2<-lm(y~xhat)
  pp<-summary(fit2)
  pp1<-pp$coefficients[2,4]
  tt1<-deviance(fit2)
  mse1<-mean(pp$residuals^2)
  list1[[length(list1)+1]]<-tt1
  list3[[length(list3)+1]]<-pp1
  
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(1,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(-0.25*x+c+3*rnorm(20))
  z<-cbind(z,x)
  xhat<-as.numeric(predict(fit1,as.data.frame(z[,1:8])))
  ypre<-as.numeric(predict(fit2,as.data.frame(xhat)))
  list5[[length(list5)+1]]<-mean((y-ypre)^2)
}
ap2<-unlist(list3)


list8<-list()
list9<-list()
list10<-list()
for (i in 1:100){
  set.seed(i)
  z<-(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(1,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(-0.25*x+c+3*rnorm(20))
  x = x-mean(x)
  z = t(t(z)-colMeans(z))
  fit5 <-  cv.glmnet(z,x,keep=TRUE,nfolds = length(x),intercept=FALSE)
  lambda = fit5$lambda.min
  xhat = as.numeric(fit5$fit.preval[,which(fit5$lambda==lambda)[1]])
  fit6 = lm(y~xhat)
  rp<-fit6$coefficients[2]
  list10[[length(list10)+1]]<-rp
  pp4 = 1
  if(nrow(coef(summary(fit6)))==2) pp4 = coef(summary(fit6))[2,4]
  list8[[length(list8)+1]]<-pp4
  
  z<-(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(1,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(-0.25*x+c+3*rnorm(20))
  x = x-mean(x)
  z = t(t(z)-colMeans(z))
  
  xhat<-as.numeric(predict(fit5 ,z, s=fit5$lambda.min))
  ppre<-as.numeric(predict(fit6,data.frame(xhat)))
  list9[[length(list9)+1]]<-mean((y-ppre)^2)
}

ap4<-unlist(list8)

lasso.mse<-as.numeric(unlist(list6))
ols.mse<-as.numeric(unlist(list5))
pv.mse<-as.numeric(unlist(list9))
v<-c(rep("TSLS",100),rep("TSL",100),rep("TSPL",100))
mse<-as.data.frame(cbind(v,c(ols.mse,lasso.mse,pv.mse)))
mse[,2]<- as.numeric(as.character(mse[,2]))
ggplot(mse,aes(x=mse$v,y=mse[,2],group=mse$v))+scale_y_continuous(limits = c(0, 100))+geom_boxplot(fill=c("darkblue","skyblue","darkgreen"))+labs(x="Model",y="MSE",title="Association between X and Y")+theme(plot.title = element_text(hjust = 0.5))

fnn<-list4>0.05
fn<-length(fnn[fnn=="TRUE"])
tp<-100-fn
fnn1<-list3>0.05
fn1<-length(fnn1[fnn1=="TRUE"])
tp1<-100-fn1
fnn2<-list8>0.05
fn2<-length(fnn2[fnn2=="TRUE"])
tp2<-100-fn2

##no association
set.seed(5)
dta1<-data.frame(x=1:8)
for (j in 1:8){
  list1<-list()
  for (i in 1:8){
    tl<-0.5^abs(i-j)
    list1[[length(list1)+1]]<-tl
  }
  dta1<-cbind(dta1,unlist(list1))
}
s<-matrix(c(dta1[,2],dta1[,3],dta1[,4],dta1[,5],dta1[,6],dta1[,7],dta1[,8],dta1[,9]),ncol=8,nrow=8)


list2<-list()
list4<-list()
list6<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(1,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(c+3*rnorm(20))
  crossval <-  cv.glmnet(y=x,x=z)
  penalty <- crossval$lambda.min
  fit3 <-glmnet( y=x,x=z, alpha = 1, lambda = penalty )
  xhat<-predict(crossval, newx = z, s = penalty)
  xhat<-as.vector(xhat[,1])
  fit4<-lm(y~xhat)
  pp2<-summary(fit4)
  pp3<-1
  if(nrow(coef(summary(fit4)))==2) pp3 = coef(summary(fit4))[2,4]
  tt2<-deviance(fit4)
  list2[[length(list2)+1]]<-tt2
  list4[[length(list4)+1]]<-pp3
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(1,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(c+3*rnorm(20))
  
  xhat<-as.numeric(predict(fit3 ,z, s=fit3$lambda.min))
  ypre<-as.numeric(predict(fit4,as.data.frame(as.numeric(xhat))))
  list6[[length(list6)+1]]<-mean((y-ypre)^2)
}
ap1<-unlist(list4)

list1<-list()
list3<-list()
list5<-list()
for (i in 1:100){
  set.seed(i)
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(1,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(c+3*rnorm(20))
  z<-cbind(z,x)
  fit1<-lm(x~.,z)
  xhat<-fit1$fitted.values
  fit2<-lm(y~xhat)
  pp<-summary(fit2)
  pp1<-pp$coefficients[2,4]
  tt1<-deviance(fit2)
  mse1<-mean(pp$residuals^2)
  list1[[length(list1)+1]]<-tt1
  list3[[length(list3)+1]]<-pp1
  
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(1,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(c+3*rnorm(20))
  z<-cbind(z,x)
  xhat<-as.numeric(predict(fit1,as.data.frame(z[,1:8])))
  ypre<-as.numeric(predict(fit2,as.data.frame(xhat)))
  list5[[length(list5)+1]]<-mean((y-ypre)^2)
}

ap3<-unlist(list3)



list8<-list()
list9<-list()
list10<-list()
for (i in 1:100){
  set.seed(i)
  z<-(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(1,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(c+3*rnorm(20))
  x = x-mean(x)
  z = t(t(z)-colMeans(z))
  fit5 <-  cv.glmnet(z,x,keep=TRUE,nfolds = length(x),intercept=FALSE)
  lambda = fit5$lambda.min
  xhat = as.numeric(fit5$fit.preval[,which(fit5$lambda==lambda)[1]])
  fit6 = lm(y~xhat)
  rp<-fit6$coefficients[2]
  list10[[length(list10)+1]]<-rp
  pp4 = 1
  if(nrow(coef(summary(fit6)))==2) pp4 = coef(summary(fit6))[2,4]
  list8[[length(list8)+1]]<-pp4
  
  z<-(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(1,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(c+3*rnorm(20))
  x = x-mean(x)
  z = t(t(z)-colMeans(z))
  
  xhat<-as.numeric(predict(fit5 ,z, s=fit5$lambda.min))
  ppre<-as.numeric(predict(fit6,data.frame(xhat)))
  list9[[length(list9)+1]]<-mean((y-ppre)^2)
}

ap5<-unlist(list8)

lasso.mse<-as.numeric(unlist(list6))
ols.mse<-as.numeric(unlist(list5))
pv.mse<-as.numeric(unlist(list9))
v<-c(rep("TSLS",100),rep("TSL",100),rep("TSPL",100))
mse<-as.data.frame(cbind(v,c(ols.mse,lasso.mse,pv.mse)))
mse[,2]<- as.numeric(as.character(mse[,2]))
ggplot(mse,aes(x=mse$v,y=mse[,2],group=mse$v))+scale_y_continuous(limits = c(0, 100))+geom_boxplot(fill=c("darkblue","skyblue","darkgreen"))+labs(x="Model",y="MSE",title="No Association between X and Y")+theme(plot.title = element_text(hjust = 0.5))

tnn<-list4>0.05
tn<-length(tnn[tnn=="TRUE"])
fp<-100-tn
tnn1<-list3>0.05
tn1<-length(tnn1[tnn1=="TRUE"])
fp1<-100-tn1
tnn2<-list8>0.05
tn2<-length(tnn2[tnn2=="TRUE"])
fp2<-100-tn2

mcc(tp,tn,fp,fn)
mcc(tp1,tn1,fp1,fn1)
mcc(tp2,tn2,fp2,fn2)

app<-c(ap,ap1)
class = rep(c(1,0),c(length(ap1),length(ap)))
ord = order(app)
plot(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),type = 'l',main = 'ROC curve',xlab = 'FPR',ylab = 'TPR',col="darkblue",lwd=4)
abline(0,1)
p = c(ap2,ap3)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="skyblue",type = 'l',lwd=4)
p = c(ap4,ap5)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="darkgreen",type = 'l',lwd=4)
legend('bottomright',c('TSL','TSLS','TSPL'),col = c("darkblue","skyblue","darkgreen"),lty = 1,cex=1,lwd=3)

df<-as.data.frame(cbind(as.numeric(c(ap,ap1,ap2,ap3,ap4,ap5)),as.factor(c(rep("Ass",100),rep("No Ass",100))),(c(rep("TSL",200),rep("TSLS",200),rep("TSPL",200)))))
colnames(df)<-c("p.value","ass","method")
df[,1]<- as.numeric(as.character(df[,1]))
df<-df[c(1:100,200:300,400:500),]
ggplot(df,aes(x=df$method,y=df$p.value,group=df$method))+labs(x="Method",y="p-value",title="Range of p-values from the Simulations")+theme(plot.title = element_text(hjust = 0.5))+geom_boxplot(fill=c("darkblue","skyblue","darkgreen"))


##two confounders

dta1<-data.frame(x=1:8)
for (j in 1:8){
  list1<-list()
  for (i in 1:8){
    tl<-0.5^abs(i-j)
    list1[[length(list1)+1]]<-tl
  }
  dta1<-cbind(dta1,unlist(list1))
}
s<-matrix(c(dta1[,2],dta1[,3],dta1[,4],dta1[,5],dta1[,6],dta1[,7],dta1[,8],dta1[,9]),ncol=8,nrow=8)


list2<-list()
list4<-list()
list6<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  c1<-rnorm(20,0,3)
  g=c(1,1,2)
  x<-g%*%t(cbind(as.numeric(b%*%t(z)),c,c1))+3*rnorm(20)
  y<-as.vector(x+-2*c+-5*c1+3*rnorm(20))
  crossval <-  cv.glmnet(y=x,x=z)
  penalty <- crossval$lambda.min
  fit3 <-glmnet( y=x,x=z, alpha = 1, lambda = penalty )
  xhat<-predict(crossval, newx = z, s = penalty)
  xhat<-as.vector(xhat[,1])
  fit4<-lm(y~xhat)
  pp2<-summary(fit4)
  pp3<-1
  if(nrow(coef(summary(fit4)))==2) pp3 = coef(summary(fit4))[2,4]
  tt2<-deviance(fit4)
  list2[[length(list2)+1]]<-tt2
  list4[[length(list4)+1]]<-pp3
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  c1<-rnorm(20,0,3)
  g=c(1,1,2)
  x<-g%*%t(cbind(as.numeric(b%*%t(z)),c,c1))+3*rnorm(20)
  y<-as.vector(x+-2*c+-5*c1+3*rnorm(20))
  
  xhat<-as.numeric(predict(fit3 ,z, s=fit3$lambda.min))
  ypre<-as.numeric(predict(fit4,as.data.frame(as.numeric(xhat))))
  list6[[length(list6)+1]]<-mean((y-ypre)^2)
}
ap<-unlist(list4)

list1<-list()
list3<-list()
list5<-list()
for (i in 1:100){
  set.seed(i)
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  c1<-rnorm(20,0,3)
  g=c(1,1,2)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c,c1))+3*rnorm(20))
  y<-as.vector(x+-2*c+-5*c1+3*rnorm(20))
  z<-cbind(z,x)
  fit1<-lm(x~.,z)
  xhat<-fit1$fitted.values
  fit2<-lm(y~xhat)
  pp<-summary(fit2)
  pp1<-pp$coefficients[2,4]
  tt1<-deviance(fit2)
  mse1<-mean(pp$residuals^2)
  list1[[length(list1)+1]]<-tt1
  list3[[length(list3)+1]]<-pp1
  
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  c1<-rnorm(20,0,3)
  g=c(1,1,2)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c,c1))+3*rnorm(20))
  y<-as.vector(x+-2*c+-5*c1+3*rnorm(20))
  z<-cbind(z,x)
  xhat<-as.numeric(predict(fit1,as.data.frame(z[,1:8])))
  ypre<-as.numeric(predict(fit2,as.data.frame(xhat)))
  list5[[length(list5)+1]]<-mean((y-ypre)^2)
}
ap2<-unlist(list3)


list8<-list()
list9<-list()
list10<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  c1<-rnorm(20,0,3)
  g=c(1,1,2)
  x<-g%*%t(cbind(as.numeric(b%*%t(z)),c,c1))+3*rnorm(20)
  y<-as.vector(x+-2*c+-5*c1+3*rnorm(20))
  x = x-mean(x)
  z = t(t(z)-colMeans(z))
  fit5 <-  cv.glmnet(z,x,keep=TRUE,nfolds = length(x),intercept=FALSE)
  lambda = fit5$lambda.min
  xhat = as.numeric(fit5$fit.preval[,which(fit5$lambda==lambda)[1]])
  fit6 = lm(y~xhat)
  rp<-fit6$coefficients[2]
  list10[[length(list10)+1]]<-rp
  pp4 = 1
  if(nrow(coef(summary(fit6)))==2) pp4 = coef(summary(fit6))[2,4]
  list8[[length(list8)+1]]<-pp4
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  c1<-rnorm(20,0,3)
  g=c(1,1,2)
  x<-g%*%t(cbind(as.numeric(b%*%t(z)),c,c1))+3*rnorm(20)
  y<-as.vector(x+-2*c+-5*c1+3*rnorm(20))
  x = x-mean(x)
  z = t(t(z)-colMeans(z))
  
  xhat<-as.numeric(predict(fit5 ,z, s=fit5$lambda.min))
  ppre<-as.numeric(predict(fit6,data.frame(xhat)))
  list9[[length(list9)+1]]<-mean((y-ppre)^2)
}

ap4<-unlist(list8)

lasso.mse<-as.numeric(unlist(list6))
ols.mse<-as.numeric(unlist(list5))
pv.mse<-as.numeric(unlist(list9))
v<-c(rep("TSLS",100),rep("TSL",100),rep("TSPL",100))
mse<-as.data.frame(cbind(v,c(ols.mse,lasso.mse,pv.mse)))
mse[,2]<- as.numeric(as.character(mse[,2]))
ggplot(mse,aes(x=mse$v,y=mse[,2],group=mse$v))+scale_y_continuous(limits = c(0, 1000))+geom_boxplot(fill=c("darkblue","skyblue","darkgreen"))+labs(x="Model",y="MSE",title="Association between X and Y")+theme(plot.title = element_text(hjust = 0.5))

fnn<-list4>0.05
fn<-length(fnn[fnn=="TRUE"])
tp<-100-fn
fnn1<-list3>0.05
fn1<-length(fnn1[fnn1=="TRUE"])
tp1<-100-fn1
fnn2<-list8>0.05
fn2<-length(fnn2[fnn2=="TRUE"])
tp2<-100-fn2

##no association
set.seed(5)
dta1<-data.frame(x=1:8)
for (j in 1:8){
  list1<-list()
  for (i in 1:8){
    tl<-0.5^abs(i-j)
    list1[[length(list1)+1]]<-tl
  }
  dta1<-cbind(dta1,unlist(list1))
}
s<-matrix(c(dta1[,2],dta1[,3],dta1[,4],dta1[,5],dta1[,6],dta1[,7],dta1[,8],dta1[,9]),ncol=8,nrow=8)


list2<-list()
list4<-list()
list6<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  c1<-rnorm(20,0,3)
  g=c(1,1,2)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c,c1))+3*rnorm(20))
  y<-as.vector(-2*c+-5*c1+3*rnorm(20))
  crossval <-  cv.glmnet(y=x,x=z)
  penalty <- crossval$lambda.min
  fit3 <-glmnet( y=x,x=z, alpha = 1, lambda = penalty )
  xhat<-predict(crossval, newx = z, s = penalty)
  xhat<-as.vector(xhat[,1])
  fit4<-lm(y~xhat)
  pp2<-summary(fit4)
  pp3<-1
  if(nrow(coef(summary(fit4)))==2) pp3 = coef(summary(fit4))[2,4]
  tt2<-deviance(fit4)
  list2[[length(list2)+1]]<-tt2
  list4[[length(list4)+1]]<-pp3
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  c1<-rnorm(20,0,3)
  g=c(1,1,2)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c,c1))+3*rnorm(20))
  y<-as.vector(-2*c+-5*c1+3*rnorm(20))
  
  xhat<-as.numeric(predict(fit3 ,z, s=fit3$lambda.min))
  ypre<-as.numeric(predict(fit4,as.data.frame(as.numeric(xhat))))
  list6[[length(list6)+1]]<-mean((y-ypre)^2)
}
ap1<-unlist(list4)

list1<-list()
list3<-list()
list5<-list()
for (i in 1:100){
  set.seed(i)
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  c1<-rnorm(20,0,3)
  g=c(1,1,2)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c,c1))+3*rnorm(20))
  y<-as.vector(-2*c+-5*c1+3*rnorm(20))
  z<-cbind(z,x)
  fit1<-lm(x~.,z)
  xhat<-fit1$fitted.values
  fit2<-lm(y~xhat)
  pp<-summary(fit2)
  pp1<-pp$coefficients[2,4]
  tt1<-deviance(fit2)
  mse1<-mean(pp$residuals^2)
  list1[[length(list1)+1]]<-tt1
  list3[[length(list3)+1]]<-pp1
  
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  c1<-rnorm(20,0,3)
  g=c(1,1,2)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c,c1))+3*rnorm(20))
  y<-as.vector(-2*c+-5*c1+3*rnorm(20))
  z<-cbind(z,x)
  xhat<-as.numeric(predict(fit1,as.data.frame(z[,1:8])))
  ypre<-as.numeric(predict(fit2,as.data.frame(xhat)))
  list5[[length(list5)+1]]<-mean((y-ypre)^2)
}

ap3<-unlist(list3)



list8<-list()
list9<-list()
list10<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  c1<-rnorm(20,0,3)
  g=c(1,1,2)
  x<-g%*%t(cbind(as.numeric(b%*%t(z)),c,c1))+3*rnorm(20)
  y<-as.vector(-2*c+-5*c1+3*rnorm(20))
  x = x-mean(x)
  z = t(t(z)-colMeans(z))
  fit5 <-  cv.glmnet(z,x,keep=TRUE,nfolds = length(x),intercept=FALSE)
  lambda = fit5$lambda.min
  xhat = as.numeric(fit5$fit.preval[,which(fit5$lambda==lambda)[1]])
  fit6 = lm(y~xhat)
  rp<-fit6$coefficients[2]
  list10[[length(list10)+1]]<-rp
  pp4 = 1
  if(nrow(coef(summary(fit6)))==2) pp4 = coef(summary(fit6))[2,4]
  list8[[length(list8)+1]]<-pp4
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  c1<-rnorm(20,0,3)
  g=c(1,1,2)
  x<-g%*%t(cbind(as.numeric(b%*%t(z)),c,c1))+3*rnorm(20)
  y<-as.vector(-2*c+-5*c1+3*rnorm(20))
  x = x-mean(x)
  z = t(t(z)-colMeans(z))
  
  xhat<-as.numeric(predict(fit5 ,z, s=fit5$lambda.min))
  ppre<-as.numeric(predict(fit6,data.frame(xhat)))
  list9[[length(list9)+1]]<-mean((y-ppre)^2)
}

ap5<-unlist(list8)

lasso.mse<-as.numeric(unlist(list6))
ols.mse<-as.numeric(unlist(list5))
pv.mse<-as.numeric(unlist(list9))
v<-c(rep("TSLS",100),rep("TSL",100),rep("TSPL",100))
mse<-as.data.frame(cbind(v,c(ols.mse,lasso.mse,pv.mse)))
mse[,2]<- as.numeric(as.character(mse[,2]))
ggplot(mse,aes(x=mse$v,y=mse[,2],group=mse$v))+geom_boxplot(fill=c("darkblue","skyblue","darkgreen"))+scale_y_continuous(limits = c(0, 1000))+labs(x="Model",y="MSE",title="No Association between X and Y")+theme(plot.title = element_text(hjust = 0.5))

tnn<-list4>0.05
tn<-length(tnn[tnn=="TRUE"])
fp<-100-tn
tnn1<-list3>0.05
tn1<-length(tnn1[tnn1=="TRUE"])
fp1<-100-tn1
tnn2<-list8>0.05
tn2<-length(tnn2[tnn2=="TRUE"])
fp2<-100-tn2

mcc(tp,tn,fp,fn)
mcc(tp1,tn1,fp1,fn1)
mcc(tp2,tn2,fp2,fn2)

app<-c(ap,ap1)
class = rep(c(1,0),c(length(ap1),length(ap)))
ord = order(app)
plot(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),type = 'l',main = 'ROC curve',xlab = 'FPR',ylab = 'TPR',col="darkblue",lwd=4)
abline(0,1)
p = c(ap2,ap3)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="skyblue",type = 'l',lwd=4)
p = c(ap4,ap5)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="darkgreen",type = 'l',lwd=4)
legend('bottomright',c('TSL','TSLS','TSPL'),col = c("darkblue","skyblue","darkgreen"),lty = 1,cex=1,lwd=3)

df<-as.data.frame(cbind(as.numeric(c(ap,ap1,ap2,ap3,ap4,ap5)),as.factor(c(rep("Ass",100),rep("No Ass",100))),(c(rep("TSL",200),rep("TSLS",200),rep("TSPL",200)))))
colnames(df)<-c("p.value","ass","method")
df[,1]<- as.numeric(as.character(df[,1]))
df<-df[c(1:100,200:300,400:500),]
ggplot(df,aes(x=df$method,y=df$p.value,group=df$method))+labs(x="Method",y="p-value",title="Range of p-values from the Simulations")+theme(plot.title = element_text(hjust = 0.5))+geom_boxplot(fill=c("darkblue","skyblue","darkgreen"))





##no informative IVs
dta1<-data.frame(x=1:8)
for (j in 1:8){
  list1<-list()
  for (i in 1:8){
    tl<-0.5^abs(i-j)
    list1[[length(list1)+1]]<-tl
  }
  dta1<-cbind(dta1,unlist(list1))
}
s<-matrix(c(dta1[,2],dta1[,3],dta1[,4],dta1[,5],dta1[,6],dta1[,7],dta1[,8],dta1[,9]),ncol=8,nrow=8)


list2<-list()
list4<-list()
list6<-list()
list12<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(.5,1)
  x<-g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20)
  y<-as.vector(x+c+3*rnorm(20))
  crossval <-  cv.glmnet(y=x,x=z)
  penalty <- crossval$lambda.min
  fit3 <-glmnet( y=x,x=z, alpha = 1, lambda = penalty )
  xhat<-predict(crossval, newx = z, s = penalty)
  xhat<-as.vector(xhat[,1])
  fit4<-lm(y~xhat)
  rl<-fit4$coefficients[2]
  list12[[length(list12)+1]]<-rl
  pp2<-summary(fit4)
  pp3<-1
  if(nrow(coef(summary(fit4)))==2) pp3 = coef(summary(fit4))[2,4]
  tt2<-deviance(fit4)
  list2[[length(list2)+1]]<-tt2
  list4[[length(list4)+1]]<-pp3
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(.5,1)
  x<-g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20)
  y<-as.vector(x+c+3*rnorm(20))
  
  xhat<-as.numeric(predict(fit3 ,z, s=fit3$lambda.min))
  ypre<-as.numeric(predict(fit4,as.data.frame(as.numeric(xhat))))
  list6[[length(list6)+1]]<-mean((y-ypre)^2)
}
ap<-unlist(list4)

list1<-list()
list3<-list()
list5<-list()
list11<-list()
for (i in 1:100){
  set.seed(i)
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(.5,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(x+c+3*rnorm(20))
  z<-cbind(z,x)
  fit1<-lm(x~.,z)
  xhat<-fit1$fitted.values
  fit2<-lm(y~xhat)
  ro<-fit2$coefficients[2]
  list11[[length(list11)+1]]<-ro
  pp<-summary(fit2)
  pp1<-pp$coefficients[2,4]
  tt1<-deviance(fit2)
  mse1<-mean(pp$residuals^2)
  list1[[length(list1)+1]]<-tt1
  list3[[length(list3)+1]]<-pp1
  
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(.5,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(x+c+3*rnorm(20))
  z<-cbind(z,x)
  xhat<-as.numeric(predict(fit1,as.data.frame(z[,1:8])))
  ypre<-as.numeric(predict(fit2,as.data.frame(xhat)))
  list5[[length(list5)+1]]<-mean((y-ypre)^2)
}
ap2<-unlist(list3)

list8<-list()
list9<-list()
list10<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(.5,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(x+c+3*rnorm(20))
  fit5 <-  cv.glmnet(y=x,x=z,keep=TRUE)
  xhat = as.numeric(fit5$fit.preval[,which(fit5$lambda==fit5$lambda.min)[1]])
  fit6 = lm(y~xhat)
  rp<-fit6$coefficients[2]
  list10[[length(list10)+1]]<-rp
  pp4 = 1
  if(nrow(coef(summary(fit6)))==2) pp4 = coef(summary(fit6))[2,4]
  list8[[length(list8)+1]]<-pp4
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(.5,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(x+c+3*rnorm(20))
  
  xhat<-as.numeric(predict(fit5 ,z, s=fit5$lambda.min))
  ppre<-as.numeric(predict(fit6,data.frame(xhat)))
  list9[[length(list9)+1]]<-mean((y-ppre)^2)
}

ap4<-unlist(list8)
cd1<-unlist(list12)
cd2<-unlist(list11)
cd3<-unlist(list10)

lasso.mse<-as.numeric(unlist(list6))
ols.mse<-as.numeric(unlist(list5))
pv.mse<-as.numeric(unlist(list9))
v<-c(rep("TSLS",100),rep("TSL",100),rep("TSPL",100))
mse<-as.data.frame(cbind(v,c(ols.mse,lasso.mse,pv.mse)))
mse[,2]<- as.numeric(as.character(mse[,2]))
ggplot(mse,aes(x=mse$v,y=mse[,2],group=mse$v))+geom_boxplot(fill=c("darkblue","skyblue","darkgreen"))+labs(x="Model",y="MSE",title="Informative Instruments")+theme(plot.title = element_text(hjust = 0.5))

fnn<-list4>0.05
fn<-length(fnn[fnn=="TRUE"])
tp<-100-fn
fnn1<-list3>0.05
fn1<-length(fnn1[fnn1=="TRUE"])
tp1<-100-fn1
fnn2<-list8>0.05
fn2<-length(fnn2[fnn2=="TRUE"])
tp2<-100-fn2

##no association
set.seed(5)
dta1<-data.frame(x=1:8)
for (j in 1:8){
  list1<-list()
  for (i in 1:8){
    tl<-0.5^abs(i-j)
    list1[[length(list1)+1]]<-tl
  }
  dta1<-cbind(dta1,unlist(list1))
}
s<-matrix(c(dta1[,2],dta1[,3],dta1[,4],dta1[,5],dta1[,6],dta1[,7],dta1[,8],dta1[,9]),ncol=8,nrow=8)


list2<-list()
list4<-list()
list6<-list()
list12<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(0,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(x+c+3*rnorm(20))
  crossval <-  cv.glmnet(y=x,x=z)
  penalty <- crossval$lambda.min
  fit3 <-glmnet( y=x,x=z, alpha = 1, lambda = penalty )
  xhat<-predict(crossval, newx = z, s = penalty)
  xhat<-as.vector(xhat[,1])
  fit4<-lm(y~xhat)
  rl<-fit4$coefficients[2]
  list12[[length(list12)+1]]<-rl
  pp2<-summary(fit4)
  pp3<-1
  if(nrow(coef(summary(fit4)))==2) pp3 = coef(summary(fit4))[2,4]
  tt2<-deviance(fit4)
  list2[[length(list2)+1]]<-tt2
  list4[[length(list4)+1]]<-pp3
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(0,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(x+c+3*rnorm(20))
  
  xhat<-as.numeric(predict(fit3 ,z, s=fit3$lambda.min))
  ypre<-as.numeric(predict(fit4,as.data.frame(as.numeric(xhat))))
  list6[[length(list6)+1]]<-mean((y-ypre)^2)
}
ap1<-unlist(list4)

list1<-list()
list3<-list()
list5<-list()
list11<-list()
for (i in 1:100){
  set.seed(i)
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(0,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(x+c+3*rnorm(20))
  z<-cbind(z,x)
  fit1<-lm(x~.,z)
  xhat<-fit1$fitted.values
  fit2<-lm(y~xhat)
  ro<-fit2$coefficients[2]
  list11[[length(list11)+1]]<-ro
  pp<-summary(fit2)
  pp1<-pp$coefficients[2,4]
  tt1<-deviance(fit2)
  mse1<-mean(pp$residuals^2)
  list1[[length(list1)+1]]<-tt1
  list3[[length(list3)+1]]<-pp1
  
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(0,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(x+c+3*rnorm(20))
  z<-cbind(z,x)
  xhat<-as.numeric(predict(fit1,as.data.frame(z[,1:8])))
  ypre<-as.numeric(predict(fit2,as.data.frame(xhat)))
  list5[[length(list5)+1]]<-mean((y-ypre)^2)
}

ap3<-unlist(list3)


list8<-list()
list9<-list()
list10<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(0,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(x+c+3*rnorm(20))
  fit5 <-  cv.glmnet(y=x,x=z,keep=TRUE)
  xhat = as.numeric(fit5$fit.preval[,which(fit5$lambda==fit5$lambda.min)[1]])
  fit6 = lm(y~xhat)
  rp<-fit6$coefficients[2]
  list10[[length(list10)+1]]<-rp
  pp4 = 1
  if(nrow(coef(summary(fit6)))==2) pp4 = coef(summary(fit6))[2,4]
  list8[[length(list8)+1]]<-pp4
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(0,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(x+c+3*rnorm(20))
  
  xhat<-as.numeric(predict(fit5 ,z, s=fit5$lambda.min))
  ppre<-as.numeric(predict(fit6,data.frame(xhat)))
  list9[[length(list9)+1]]<-mean((y-ppre)^2)
}
ap5<-unlist(list8)

cd4<-unlist(list12)
cd5<-unlist(list11)
cd6<-unlist(list10)

lasso.mse<-as.numeric(unlist(list6))
ols.mse<-as.numeric(unlist(list5))
pv.mse<-as.numeric(unlist(list9))
v<-c(rep("TSLS",100),rep("TSL",100),rep("TSPL",100))
mse<-as.data.frame(cbind(v,c(ols.mse,lasso.mse,pv.mse)))
mse[,2]<- as.numeric(as.character(mse[,2]))
ggplot(mse,aes(x=mse$v,y=mse[,2],group=mse$v))+geom_boxplot(fill=c("darkblue","skyblue","darkgreen"))+labs(x="Model",y="MSE",title="No Informative Instruments")+theme(plot.title = element_text(hjust = 0.5))

tnn<-list4>0.05
tn<-length(tnn[tnn=="TRUE"])
fp<-100-tn
tnn1<-list3>0.05
tn1<-length(tnn1[tnn1=="TRUE"])
fp1<-100-tn1
tnn2<-list8>0.05
tn2<-length(tnn2[tnn2=="TRUE"])
fp2<-100-tn2



#make MCC function
mcc<-function(a,b,c,d){
  (a*b-c*d)/sqrt((a+c)*(b+d)*(c+b)*(a+d))
}
#a=tp b=tn c=fp d=fn

mcc(tp,tn,fp,fn)
mcc(tp1,tn1,fp1,fn1)
mcc(tp2,tn2,fp2,fn2)

app<-c(ap,ap1)
class = rep(c(1,0),c(length(ap1),length(ap)))
ord = order(app)
plot(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),type = 'l',main = 'ROC curve',xlab = 'FPR',ylab = 'TPR',col="darkblue",lwd=4)
abline(0,1)
p = c(ap2,ap3)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="skyblue",type = 'l',lwd=4)
p = c(ap4,ap5)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="darkgreen",type = 'l',lwd=4)
legend('bottomright',c('TSL','TSLS','TSPL'),col = c("darkblue","skyblue","darkgreen"),lty = 1,cex=1,lwd=3)

df<-as.data.frame(cbind(as.numeric(c(ap,ap1,ap2,ap3,ap4,ap5)),as.factor(c(rep("Ass",100),rep("No Ass",100))),(c(rep("TSL",200),rep("TSLS",200),rep("TSPL",200)))))
colnames(df)<-c("p.value","ass","method")
df[,1]<- as.numeric(as.character(df[,1]))
df<-df[c(1:100,200:300,400:500),]
ggplot(df,aes(x=df$method,y=df$p.value,group=df$method))+labs(x="Method",y="p-value",title="Range of p-values from the Simulations")+theme(plot.title = element_text(hjust = 0.5))+geom_boxplot(fill=c("darkblue","skyblue","darkgreen"))


##Preval small k-folds
##no informative IVs
dta1<-data.frame(x=1:8)
for (j in 1:8){
  list1<-list()
  for (i in 1:8){
    tl<-0.5^abs(i-j)
    list1[[length(list1)+1]]<-tl
  }
  dta1<-cbind(dta1,unlist(list1))
}
s<-matrix(c(dta1[,2],dta1[,3],dta1[,4],dta1[,5],dta1[,6],dta1[,7],dta1[,8],dta1[,9]),ncol=8,nrow=8)


list2<-list()
list4<-list()
list6<-list()
list12<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(4,3)
  x<-g%*%t(cbind(as.numeric(b%*%t(z)),c))+2*rnorm(20)
  y<-as.vector(.5*x+3*c+3*rnorm(20))
  x = x-mean(x)
  z = t(t(z)-colMeans(z))
  crossval <-  cv.glmnet(y=x,x=z)
  penalty <- crossval$lambda.min
  fit3 <-glmnet( y=x,x=z, alpha = 1, lambda = penalty )
  xhat<-predict(crossval, newx = z, s = penalty)
  xhat<-as.vector(xhat[,1])
  fit4<-lm(y~xhat)
  rl<-fit4$coefficients[2]
  list12[[length(list12)+1]]<-rl
  pp2<-summary(fit4)
  pp3<-1
  if(nrow(coef(summary(fit4)))==2) pp3 = coef(summary(fit4))[2,4]
  tt2<-deviance(fit4)
  list2[[length(list2)+1]]<-tt2
  list4[[length(list4)+1]]<-pp3
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(4,3)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+2*rnorm(20))
  y<-as.vector(.5*x+3*c+3*rnorm(20))
  x = x-mean(x)
  z = t(t(z)-colMeans(z))
  
  xhat<-as.numeric(predict(fit3 ,z, s=fit3$lambda.min))
  ypre<-as.numeric(predict(fit4,as.data.frame(as.numeric(xhat))))
  list6[[length(list6)+1]]<-mean((y-ypre)^2)
}
ap<-unlist(list4)

list1<-list()
list3<-list()
list5<-list()
list11<-list()
for (i in 1:100){
  set.seed(i)
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(4,3)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+2*rnorm(20))
  y<-as.vector(.5*x+3*c+3*rnorm(20))
  x = x-mean(x)
  z = as.data.frame(t(t(z)-colMeans(z)))
  z<-cbind(z,x)
  fit1<-lm(x~.,z)
  xhat<-fit1$fitted.values
  fit2<-lm(y~xhat)
  ro<-fit2$coefficients[2]
  list11[[length(list11)+1]]<-ro
  pp<-summary(fit2)
  pp1<-pp$coefficients[2,4]
  tt1<-deviance(fit2)
  mse1<-mean(pp$residuals^2)
  list1[[length(list1)+1]]<-tt1
  list3[[length(list3)+1]]<-pp1
  
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(4,3)
  y<-as.vector(.5*x+3*c+3*rnorm(20))
  x = x-mean(x)
  z = as.data.frame(t(t(z)-colMeans(z)))
  z<-cbind(z,x)
  xhat<-as.numeric(predict(fit1,as.data.frame(z[,1:8])))
  ypre<-as.numeric(predict(fit2,as.data.frame(xhat)))
  list5[[length(list5)+1]]<-mean((y-ypre)^2)
}
ap2<-unlist(list3)

list8<-list()
list9<-list()
list10<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(4,3)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+2*rnorm(20))
  y<-as.vector(.5*x+3*c+3*rnorm(20))
  x = x-mean(x)
  z = t(t(z)-colMeans(z))
  fit5 <-  cv.glmnet(z,x,keep=TRUE,nfolds = length(x),intercept=FALSE)
  lambda = fit5$lambda.min
  xhat = as.numeric(fit5$fit.preval[,which(fit5$lambda==lambda)[1]])
  fit6 = lm(y~xhat)
  rp<-fit6$coefficients[2]
  list10[[length(list10)+1]]<-rp
  pp4 = 1
  if(nrow(coef(summary(fit6)))==2) pp4 = coef(summary(fit6))[2,4]
  list8[[length(list8)+1]]<-pp4
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(4,3)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+2*rnorm(20))
  y<-as.vector(.5*x+3*c+3*rnorm(20))
  x = x-mean(x)
  z = t(t(z)-colMeans(z))
  
  xhat<-as.numeric(predict(fit5 ,z, s=fit5$lambda.min))
  ppre<-as.numeric(predict(fit6,data.frame(xhat)))
  list9[[length(list9)+1]]<-mean((y-ppre)^2)
}

ap4<-unlist(list8)
cd1<-unlist(list12)
cd2<-unlist(list11)
cd3<-unlist(list10)

lasso.mse<-as.numeric(unlist(list6))
ols.mse<-as.numeric(unlist(list5))
pv.mse<-as.numeric(unlist(list9))
v<-c(rep("TSLS",100),rep("TSL",100),rep("TSPL",100))
mse<-as.data.frame(cbind(v,c(ols.mse,lasso.mse,pv.mse)))
mse[,2]<- as.numeric(as.character(mse[,2]))
ggplot(mse,aes(x=mse$v,y=mse[,2],group=mse$v))+geom_boxplot(fill=c("darkblue","skyblue","darkgreen"))+labs(x="Model",y="MSE",title="Informative Instruments")+theme(plot.title = element_text(hjust = 0.5))

fnn<-list4>0.05
fn<-length(fnn[fnn=="TRUE"])
tp<-100-fn
fnn1<-list3>0.05
fn1<-length(fnn1[fnn1=="TRUE"])
tp1<-100-fn1
fnn2<-list8>0.05
fn2<-length(fnn2[fnn2=="TRUE"])
tp2<-100-fn2

##no association
set.seed(5)
dta1<-data.frame(x=1:8)
for (j in 1:8){
  list1<-list()
  for (i in 1:8){
    tl<-0.5^abs(i-j)
    list1[[length(list1)+1]]<-tl
  }
  dta1<-cbind(dta1,unlist(list1))
}
s<-matrix(c(dta1[,2],dta1[,3],dta1[,4],dta1[,5],dta1[,6],dta1[,7],dta1[,8],dta1[,9]),ncol=8,nrow=8)


list2<-list()
list4<-list()
list6<-list()
list12<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(0,3)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+2*rnorm(20))
  y<-as.vector(.5*x+3*c+3*rnorm(20))
  x = x-mean(x)
  z = t(t(z)-colMeans(z))
  crossval <-  cv.glmnet(y=x,x=z)
  penalty <- crossval$lambda.min
  fit3 <-glmnet( y=x,x=z, alpha = 1, lambda = penalty )
  xhat<-predict(crossval, newx = z, s = penalty)
  xhat<-as.vector(xhat[,1])
  fit4<-lm(y~xhat)
  rl<-fit4$coefficients[2]
  list12[[length(list12)+1]]<-rl
  pp2<-summary(fit4)
  pp3<-1
  if(nrow(coef(summary(fit4)))==2) pp3 = coef(summary(fit4))[2,4]
  tt2<-deviance(fit4)
  list2[[length(list2)+1]]<-tt2
  list4[[length(list4)+1]]<-pp3
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(0,3)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+2*rnorm(20))
  y<-as.vector(.5*x+3*c+3*rnorm(20))
  x = x-mean(x)
  z = t(t(z)-colMeans(z))
  
  xhat<-as.numeric(predict(fit3 ,z, s=fit3$lambda.min))
  ypre<-as.numeric(predict(fit4,as.data.frame(as.numeric(xhat))))
  list6[[length(list6)+1]]<-mean((y-ypre)^2)
}
ap1<-unlist(list4)

list1<-list()
list3<-list()
list5<-list()
list11<-list()
for (i in 1:100){
  set.seed(i)
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(0,3)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+2*rnorm(20))
  y<-as.vector(.5*x+3*c+3*rnorm(20))
  x = x-mean(x)
  z = as.data.frame(t(t(z)-colMeans(z)))
  z<-(cbind(x,z))
  fit1<-lm(x~.,z)
  xhat<-fit1$fitted.values
  fit2<-lm(y~xhat)
  ro<-fit2$coefficients[2]
  list11[[length(list11)+1]]<-ro
  pp<-summary(fit2)
  pp1<-pp$coefficients[2,4]
  tt1<-deviance(fit2)
  mse1<-mean(pp$residuals^2)
  list1[[length(list1)+1]]<-tt1
  list3[[length(list3)+1]]<-pp1
  
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(0,3)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+2*rnorm(20))
  y<-as.vector(.5*x+3*c+3*rnorm(20))
  x = x-mean(x)
  z = as.data.frame(t(t(z)-colMeans(z)))
  z<-cbind(z,x)
  xhat<-as.numeric(predict(fit1,as.data.frame(z[,1:8])))
  ypre<-as.numeric(predict(fit2,as.data.frame(xhat)))
  list5[[length(list5)+1]]<-mean((y-ypre)^2)
}

ap3<-unlist(list3)


list8<-list()
list9<-list()
list10<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(0,3)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+2*rnorm(20))
  y<-as.vector(.5*x+3*c+3*rnorm(20))
  x = x-mean(x)
  z = t(t(z)-colMeans(z))
  fit5 <-  cv.glmnet(z,x,keep=TRUE,nfolds = length(x),intercept=FALSE)
  lambda = fit5$lambda.min
  xhat = as.numeric(fit5$fit.preval[,which(fit5$lambda==lambda)[1]])
  fit6 = lm(y~xhat)
  rp<-fit6$coefficients[2]
  list10[[length(list10)+1]]<-rp
  pp4 = 1
  if(nrow(coef(summary(fit6)))==2) pp4 = coef(summary(fit6))[2,4]
  list8[[length(list8)+1]]<-pp4
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(0,3)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+2*rnorm(20))
  y<-as.vector(.5*x+3*c+3*rnorm(20))
  x = x-mean(x)
  z = t(t(z)-colMeans(z))
  
  xhat<-as.numeric(predict(fit5 ,z, s=fit5$lambda.min))
  ppre<-as.numeric(predict(fit6,data.frame(xhat)))
  list9[[length(list9)+1]]<-mean((y-ppre)^2)
}
ap5<-unlist(list8)

cd4<-unlist(list12)
cd5<-unlist(list11)
cd6<-unlist(list10)

lasso.mse<-as.numeric(unlist(list6))
ols.mse<-as.numeric(unlist(list5))
pv.mse<-as.numeric(unlist(list9))
v<-c(rep("TSLS",100),rep("TSL",100),rep("TSPL",100))
mse<-as.data.frame(cbind(v,c(ols.mse,lasso.mse,pv.mse)))
mse[,2]<- as.numeric(as.character(mse[,2]))
ggplot(mse,aes(x=mse$v,y=mse[,2],group=mse$v))+geom_boxplot(fill=c("darkblue","skyblue","darkgreen"))+labs(x="Model",y="MSE",title="No Informative Instruments")+theme(plot.title = element_text(hjust = 0.5))


tnn<-list4>0.05
tn<-length(tnn[tnn=="TRUE"])
fp<-100-tn
tnn1<-list3>0.05
tn1<-length(tnn1[tnn1=="TRUE"])
fp1<-100-tn1
tnn2<-list8>0.05
tn2<-length(tnn2[tnn2=="TRUE"])
fp2<-100-tn2



#make MCC function
mcc<-function(a,b,c,d){
  (a*b-c*d)/sqrt((a+c)*(b+d)*(c+b)*(a+d))
}
#a=tp b=tn c=fp d=fn

mcc(tp,tn,fp,fn)
mcc(tp1,tn1,fp1,fn1)
mcc(tp2,tn2,fp2,fn2)

app<-c(ap,ap1)
class = rep(c(1,0),c(length(ap1),length(ap)))
ord = order(app)
plot(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),type = 'l',main = 'ROC curve',xlab = 'FPR',ylab = 'TPR',col="darkblue",lwd=4)
abline(0,1)
p = c(ap2,ap3)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="skyblue",type = 'l',lwd=4)
p = c(ap4,ap5)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="darkgreen",type = 'l',lwd=4)
legend('bottomright',c('TSL','TSLS','TSPL'),col = c("darkblue","skyblue","darkgreen"),lty = 1,cex=1,lwd=3)

df<-as.data.frame(cbind(as.numeric(c(ap,ap1,ap2,ap3,ap4,ap5)),as.factor(c(rep("Ass",100),rep("No Ass",100))),(c(rep("TSL",200),rep("TSLS",200),rep("TSPL",200)))))
colnames(df)<-c("p.value","ass","method")
df[,1]<- as.numeric(as.character(df[,1]))
df<-df[c(1:100,200:300,400:500),]
ggplot(df,aes(x=df$method,y=df$p.value,group=df$method))+labs(x="Method",y="p-value",title="Range of p-values from the Simulations")+theme(plot.title = element_text(hjust = 0.5))+geom_boxplot(fill=c("darkblue","skyblue","darkgreen"))

##bigsim
b = c(3,1.5,0,0,2,0,0,0,rep(0,40),1,4,6,rep(0,100),1,7,rep(0,100),2,rep(0,120))
i = matrix(1:length(b),length(b),length(b))
s = 0.5^(abs(i-t(i)))
list12<-list()
list13<-list()
for (i in 1:100){
  set.seed(i)
  z<-as.data.frame(mvrnorm(1000,rep(2,374),s))
  c<-rnorm(1000,0,6)
  x<-as.vector(b%*%t(z)+20*c+2*rnorm(1000,0,40))
  y<-as.vector(x+6*c+2*rnorm(1000,0,30))
  z<-as.data.frame(cbind(z,x))
  fit1<-lm(y~x)
  pp2<-summary(fit1)
  pp3<-pp2$coefficients[2,4]
  list12[[length(list12)+1]]<-pp3
  
  z<-mvrnorm(1000,rep(2,374),s)
  c<-rnorm(1000,0,6)
  x<-as.vector(b%*%t(z)+20*c+2*rnorm(1000,0,40))
  y<-as.vector(x+6*c+2*rnorm(1000,0,30))
  z<-as.data.frame(cbind(z,x))
  ypre<-as.numeric(predict(fit1,as.data.frame(x)))
  list13[[length(list13)+1]]<-mean((y-ypre)^2)
}
ap6<-unlist(list12)


list2<-list()
list4<-list()
list6<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(1000,rep(2,374),s)
  c<-rnorm(1000,0,6)
  x<-as.vector(b%*%t(z)+20*c+2*rnorm(1000,0,40))
  y<-as.vector(x+6*c+2*rnorm(1000,0,30))
  crossval <-  cv.glmnet(y=x,x=z)
  penalty <- crossval$lambda.min
  fit3 <-glmnet( y=x,x=z, alpha = 1, lambda = penalty )
  xhat<-predict(crossval, newx = z, s = penalty)
  xhat<-as.vector(xhat[,1])
  fit4<-lm(y~xhat)
  pp2<-summary(fit4)
  pp3<-1
  if(nrow(coef(summary(fit4)))==2) pp3 = coef(summary(fit4))[2,4]
  tt2<-deviance(fit4)
  list2[[length(list2)+1]]<-tt2
  list4[[length(list4)+1]]<-pp3
  
  z<-mvrnorm(1000,rep(2,374),s)
  c<-rnorm(1000,0,6)
  x<-as.vector(b%*%t(z)+20*c+2*rnorm(1000,0,40))
  y<-as.vector(x+6*c+2*rnorm(1000,0,30))
  
  xhat<-as.numeric(predict(fit3 ,z, s=fit3$lambda.min))
  ypre<-as.numeric(predict(fit4,as.data.frame(as.numeric(xhat))))
  list6[[length(list6)+1]]<-mean((y-ypre)^2)
}
ap<-unlist(list4)

list1<-list()
list3<-list()
list5<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(1000,rep(2,374),s)
  c<-rnorm(1000,0,6)
  x<-as.vector(b%*%t(z)+20*c+2*rnorm(1000,0,40))
  y<-as.vector(x+6*c+2*rnorm(1000,0,30))
  z<-as.data.frame(cbind(z,x))
  fit1<-lm(x~.,z)
  xhat<-fit1$fitted.values
  fit2<-lm(y~xhat)
  pp<-summary(fit2)
  pp1<-pp$coefficients[2,4]
  tt1<-deviance(fit2)
  mse1<-mean(pp$residuals^2)
  list1[[length(list1)+1]]<-tt1
  list3[[length(list3)+1]]<-pp1
  
  z<-mvrnorm(1000,rep(2,374),s)
  c<-rnorm(1000,0,6)
  x<-as.vector(b%*%t(z)+20*c+2*rnorm(1000,0,40))
  y<-as.vector(x+6*c+2*rnorm(1000,0,30))
  z<-cbind(z,x)
  xhat<-as.numeric(predict(fit1,as.data.frame(z[,1:374])))
  ypre<-as.numeric(predict(fit2,as.data.frame(xhat)))
  list5[[length(list5)+1]]<-mean((y-ypre)^2)
}
ap2<-unlist(list3)

list8<-list()
list9<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(1000,rep(2,374),s)
  c<-rnorm(1000,0,6)
  x<-as.vector(b%*%t(z)+20*c+2*rnorm(1000,0,40))
  y<-as.vector(x+6*c+2*rnorm(1000,0,30))
  fit5 <-  cv.glmnet(y=x,x=z,keep=TRUE)
  xhat = as.numeric(fit5$fit.preval[,which(fit5$lambda==fit5$lambda.min)[1]])
  fit6 = lm(y~xhat)
  pp4 = 1
  if(nrow(coef(summary(fit6)))==2) pp4 = coef(summary(fit6))[2,4]
  list8[[length(list8)+1]]<-pp4
  
  z<-mvrnorm(1000,rep(2,374),s)
  c<-rnorm(1000,0,6)
  x<-as.vector(b%*%t(z)+20*c+2*rnorm(1000,0,40))
  y<-as.vector(x+6*c+2*rnorm(1000,0,30))
  xhat<-as.numeric(predict(fit5 ,z, s=fit5$lambda.min))
  ppre<-as.numeric(predict(fit6,data.frame(xhat)))
  list9[[length(list9)+1]]<-mean((y-ppre)^2)
}

ap4<-unlist(list8)

lasso.mse<-as.numeric(unlist(list6))
ols.mse<-as.numeric(unlist(list5))
pv.mse<-as.numeric(unlist(list9))
ol.mse<-as.numeric(unlist(list13))
v<-c(rep("TSLS",100),rep("TSL",100),rep("TSPL",100),rep("OLS",100))
mse<-as.data.frame(cbind(v,c(ols.mse,lasso.mse,pv.mse,ol.mse)))
mse[,2]<- as.numeric(as.character(mse[,2]))
ggplot(mse,aes(x=mse$v,y=mse[,2],group=mse$v))+geom_boxplot(fill=c("green","yellow","blue","red"))+labs(x="Model",y="MSE",title="Association between X and Y")+theme(plot.title = element_text(hjust = 0.5))

fnn<-list4>0.05
fn<-length(fnn[fnn=="TRUE"])
tp<-100-fn
fnn1<-list3>0.05
fn1<-length(fnn1[fnn1=="TRUE"])
tp1<-100-fn1
fnn2<-list8>0.05
fn2<-length(fnn2[fnn2=="TRUE"])
tp2<-100-fn2
fnn3<-list12>0.05
fn3<-length(fnn3[fnn3=="TRUE"])
tp3<-100-fn3

##no association
list12<-list()
list13<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(1000,rep(2,374),s)
  c<-rnorm(1000,0,6)
  x<-as.vector(b%*%t(z)+20*c+2*rnorm(1000,0,40))
  y<-as.vector(6*c+2*rnorm(1000,0,30))
  z<-as.data.frame(cbind(z,x))
  fit1<-lm(y~x)
  pp2<-summary(fit1)
  pp3<-pp2$coefficients[2,4]
  list12[[length(list12)+1]]<-pp3
  
  z<-mvrnorm(1000,rep(2,374),s)
  c<-rnorm(1000,0,6)
  x<-as.vector(b%*%t(z)+20*c+2*rnorm(1000,0,40))
  y<-as.vector(6*c+2*rnorm(1000,0,30))
  ypre<-as.numeric(predict(fit1,as.data.frame(x)))
  list13[[length(list13)+1]]<-mean((y-ypre)^2)
}
ap7<-unlist(list12)

list2<-list()
list4<-list()
list6<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(1000,rep(2,374),s)
  c<-rnorm(1000,0,6)
  x<-as.vector(b%*%t(z)+20*c+2*rnorm(1000,0,40))
  y<-as.vector(6*c+2*rnorm(1000,0,30))
  crossval <-  cv.glmnet(y=x,x=z)
  penalty <- crossval$lambda.min
  fit3 <-glmnet( y=x,x=z, alpha = 1, lambda = penalty )
  xhat<-predict(crossval, newx = z, s = penalty)
  xhat<-as.vector(xhat[,1])
  fit4<-lm(y~xhat)
  pp2<-summary(fit4)
  pp3<-1
  if(nrow(coef(summary(fit4)))==2) pp3 = coef(summary(fit4))[2,4]
  tt2<-deviance(fit4)
  list2[[length(list2)+1]]<-tt2
  list4[[length(list4)+1]]<-pp3
  
  z<-mvrnorm(1000,rep(2,374),s)
  c<-rnorm(1000,0,6)
  x<-as.vector(b%*%t(z)+20*c+2*rnorm(1000,0,40))
  y<-as.vector(6*c+2*rnorm(1000,0,30))
  
  xhat<-as.numeric(predict(fit3 ,z, s=fit3$lambda.min))
  ypre<-as.numeric(predict(fit4,as.data.frame(as.numeric(xhat))))
  list6[[length(list6)+1]]<-mean((y-ypre)^2)
}
ap1<-unlist(list4)

list1<-list()
list3<-list()
list5<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(1000,rep(2,374),s)
  c<-rnorm(1000,0,6)
  c<-rnorm(1000,0,6)
  x<-as.vector(b%*%t(z)+20*c+2*rnorm(1000,0,40))
  y<-as.vector(6*c+2*rnorm(1000,0,30))
  z<-as.data.frame(cbind(z,x))
  fit1<-lm(x~.,z)
  xhat<-fit1$fitted.values
  fit2<-lm(y~xhat)
  pp<-summary(fit2)
  pp1<-pp$coefficients[2,4]
  tt1<-deviance(fit2)
  mse1<-mean(pp$residuals^2)
  list1[[length(list1)+1]]<-tt1
  list3[[length(list3)+1]]<-pp1
  
  z<-mvrnorm(1000,rep(2,374),s)
  c<-rnorm(1000,0,6)
  x<-as.vector(b%*%t(z)+20*c+2*rnorm(1000,0,40))
  y<-as.vector(6*c+2*rnorm(1000,0,30))
  z<-cbind(z,x)
  xhat<-as.numeric(predict(fit1,as.data.frame(z[,1:374])))
  ypre<-as.numeric(predict(fit2,as.data.frame(xhat)))
  list5[[length(list5)+1]]<-mean((y-ypre)^2)
}

ap3<-unlist(list3)


list8<-list()
list9<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(1000,rep(2,374),s)
  c<-rnorm(1000,0,6)
  x<-as.vector(b%*%t(z)+20*c+2*rnorm(1000,0,40))
  y<-as.vector(6*c+2*rnorm(1000,0,30))
  fit5 <-  cv.glmnet(y=x,x=z,keep=TRUE)
  xhat = as.numeric(fit5$fit.preval[,which(fit5$lambda==fit5$lambda.min)[1]])
  fit6 = lm(y~xhat)
  pp4 = 1
  if(nrow(coef(summary(fit6)))==2) pp4 = coef(summary(fit6))[2,4]
  list8[[length(list8)+1]]<-pp4
  
  z<-mvrnorm(1000,rep(2,374),s)
  c<-rnorm(1000,0,6)
  x<-as.vector(b%*%t(z)+20*c+2*rnorm(1000,0,40))
  y<-as.vector(6*c+2*rnorm(1000,0,30))
  
  xhat<-as.numeric(predict(fit5 ,z, s=fit5$lambda.min))
  ppre<-as.numeric(predict(fit6,data.frame(xhat)))
  list9[[length(list9)+1]]<-mean((y-ppre)^2)
}
ap5<-unlist(list8)

lasso.mse<-as.numeric(unlist(list6))
ols.mse<-as.numeric(unlist(list5))
pv.mse<-as.numeric(unlist(list9))
ol.mse<-as.numeric(unlist(list13))
v<-c(rep("TSLS",100),rep("TSL",100),rep("TSPL",100),rep("OLS",100))
mse<-as.data.frame(cbind(v,c(ols.mse,lasso.mse,pv.mse,ol.mse)))
mse[,2]<- as.numeric(as.character(mse[,2]))
ggplot(mse,aes(x=mse$v,y=mse[,2],group=mse$v))+geom_boxplot(fill=c("green","yellow","blue","red"))+labs(x="Model",y="MSE",title="No Association between X and Y")+theme(plot.title = element_text(hjust = 0.5))

tnn<-list4>0.05
tn<-length(tnn[tnn=="TRUE"])
fp<-100-tn
tnn1<-list3>0.05
tn1<-length(tnn1[tnn1=="TRUE"])
fp1<-100-tn1
tnn2<-list8>0.05
tn2<-length(tnn2[tnn2=="TRUE"])
fp2<-100-tn2
tnn3<-list12>0.05
tn3<-length(tnn3[tnn3=="TRUE"])
fp3<-100-tn3

#make MCC function
mcc<-function(a,b,c,d){
  (a*b-c*d)/sqrt((a+c)*(b+d)*(c+b)*(a+d))
}
#a=tp b=tn c=fp d=fn

mcc(tp,tn,fp,fn)
mcc(tp1,tn1,fp1,fn1)
mcc(tp2,tn2,fp2,fn2)
mcc(tp3,tn3,fp3,fn3)


#ROC plot
app<-c(ap,ap1)
class = rep(c(1,0),c(length(ap1),length(ap)))
ord = order(app)
plot(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),type = 'l',main = 'ROC curve',xlab = 'FPR',ylab = 'TPR',col="darkblue",lwd=4)
abline(0,1)
p = c(ap2,ap3)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="skyblue",type = 'l',lwd=4)
p = c(ap4,ap5)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="darkgreen",type = 'l',lwd=4)
legend('bottomright',c('TSL','TSLS','TSPL'),col = c("darkblue","skyblue","darkgreen"),lty = 1,cex=1,lwd=3)



##binary outcome

#base
dta1<-data.frame(x=1:8)
for (j in 1:8){
  list1<-list()
  for (i in 1:8){
    tl<-0.5^abs(i-j)
    list1[[length(list1)+1]]<-tl
  }
  dta1<-cbind(dta1,unlist(list1))
}
s<-matrix(c(dta1[,2],dta1[,3],dta1[,4],dta1[,5],dta1[,6],dta1[,7],dta1[,8],dta1[,9]),ncol=8,nrow=8)


list2<-list()
list4<-list()
list6<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-b%*%t(z)+3*c+2*rnorm(20)
  y<-as.vector(x+c+2*rnorm(20,0,2))
  y <- ifelse(y >= median(y), 1, 0)
  crossval <-  cv.glmnet(y=x,x=z)
  penalty <- crossval$lambda.min
  fit3 <-glmnet( y=x,x=z, alpha = 1, lambda = penalty )
  xhat<-predict(crossval, newx = z, s = penalty)
  xhat<-as.vector(xhat[,1])
  fit4<-glm(y~xhat)
  pp2<-summary(fit4)
  pp3<-1
  if(nrow(coef(summary(fit4)))==2) pp3 = coef(summary(fit4))[2,4]
  tt2<-deviance(fit4)
  list2[[length(list2)+1]]<-tt2
  list4[[length(list4)+1]]<-pp3
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-b%*%t(z)+3*c+2*rnorm(20)
  y<-as.vector(x+c+2*rnorm(20,0,2))
  y <- ifelse(y >= median(y), 1, 0)
  
  xhat<-as.numeric(predict(fit3 ,z, s=fit3$lambda.min))
  ypre<-as.numeric(predict(fit4,as.data.frame(as.numeric(xhat))))
  list6[[length(list6)+1]]<-mean((y-ypre)^2)
}
ap<-unlist(list4)

list1<-list()
list3<-list()
list5<-list()
for (i in 1:100){
  set.seed(i)
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+3*c+2*rnorm(20))
  y<-as.vector(x+c+2*rnorm(20,0,2))
  y <- ifelse(y >= median(y), 1, 0)
  z<-cbind(z,x)
  fit1<-lm(x~.,z)
  xhat<-fit1$fitted.values
  fit2<-lm(y~xhat)
  pp<-summary(fit2)
  pp1<-pp$coefficients[2,4]
  tt1<-deviance(fit2)
  mse1<-mean(pp$residuals^2)
  list1[[length(list1)+1]]<-tt1
  list3[[length(list3)+1]]<-pp1
  
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-b%*%t(z)+3*c+2*rnorm(20)
  y<-as.vector(x+c+2*rnorm(20,0,2))
  y <- ifelse(y >= median(y), 1, 0)
  z<-cbind(z,x)
  xhat<-as.numeric(predict(fit1,as.data.frame(z[,1:8])))
  ypre<-as.numeric(predict(fit2,as.data.frame(xhat)))
  list5[[length(list5)+1]]<-mean((y-ypre)^2)
}
ap2<-unlist(list3)

list8<-list()
list9<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-b%*%t(z)+3*c+2*rnorm(20)
  y<-as.vector(x+c+2*rnorm(20,0,2))
  y <- ifelse(y >= median(y), 1, 0)
  fit5 <-  cv.glmnet(y=x,x=z,keep=TRUE)
  xhat = as.numeric(fit5$fit.preval[,which(fit5$lambda==fit5$lambda.min)[1]])
  fit6 = lm(y~xhat)
  pp4 = 1
  if(nrow(coef(summary(fit6)))==2) pp4 = coef(summary(fit6))[2,4]
  list8[[length(list8)+1]]<-pp4
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-b%*%t(z)+3*c+2*rnorm(20)
  y<-as.vector(x+c+2*rnorm(20,0,2))
  y <- ifelse(y >= median(y), 1, 0)
  
  xhat<-as.numeric(predict(fit5 ,z, s=fit5$lambda.min))
  ppre<-as.numeric(predict(fit6,data.frame(xhat)))
  list9[[length(list9)+1]]<-mean((y-ppre)^2)
}

ap4<-unlist(list8)

lasso.mse<-as.numeric(unlist(list6))
ols.mse<-as.numeric(unlist(list5))
pv.mse<-as.numeric(unlist(list9))
v<-c(rep("OLS",100),rep("Lasso",100),rep("PV",100))
mse<-as.data.frame(cbind(v,c(ols.mse,lasso.mse,pv.mse)))
mse[,2]<- as.numeric(as.character(mse[,2]))
ggplot(mse,aes(x=mse$v,y=mse[,2],group=mse$v))+geom_boxplot(fill=c("yellow","blue","red"))+labs(x="Model",y="MSE",title="Association between X and Y")+theme(plot.title = element_text(hjust = 0.5))

fnn<-list4>0.05
fn<-length(fnn[fnn=="TRUE"])
tp<-100-fn
fnn1<-list3>0.05
fn1<-length(fnn1[fnn1=="TRUE"])
tp1<-100-fn1
fnn2<-list8>0.05
fn2<-length(fnn2[fnn2=="TRUE"])
tp2<-100-fn2

##no association
set.seed(5)
dta1<-data.frame(x=1:8)
for (j in 1:8){
  list1<-list()
  for (i in 1:8){
    tl<-0.5^abs(i-j)
    list1[[length(list1)+1]]<-tl
  }
  dta1<-cbind(dta1,unlist(list1))
}
s<-matrix(c(dta1[,2],dta1[,3],dta1[,4],dta1[,5],dta1[,6],dta1[,7],dta1[,8],dta1[,9]),ncol=8,nrow=8)


list2<-list()
list4<-list()
list6<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20)
  x<-b%*%t(z)+3*c+2*rnorm(20)
  y<-as.vector(c+2*rnorm(20,0,2))
  y <- ifelse(y >= median(y), 1, 0)
  crossval <-  cv.glmnet(y=x,x=z)
  penalty <- crossval$lambda.min
  fit3 <-glmnet( y=x,x=z, alpha = 1, lambda = penalty )
  xhat<-predict(crossval, newx = z, s = penalty)
  xhat<-as.vector(xhat[,1])
  fit4<-lm(y~xhat)
  pp2<-summary(fit4)
  pp3<-1
  if(nrow(coef(summary(fit4)))==2) pp3 = coef(summary(fit4))[2,4]
  tt2<-deviance(fit4)
  list2[[length(list2)+1]]<-tt2
  list4[[length(list4)+1]]<-pp3
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-b%*%t(z)+3*c+2*rnorm(20)
  y<-as.vector(c+2*rnorm(20,0,2))
  y <- ifelse(y >= median(y), 1, 0)
  
  xhat<-as.numeric(predict(fit3 ,z, s=fit3$lambda.min))
  ypre<-as.numeric(predict(fit4,as.data.frame(as.numeric(xhat))))
  list6[[length(list6)+1]]<-mean((y-ypre)^2)
}
ap1<-unlist(list4)

list1<-list()
list3<-list()
list5<-list()
for (i in 1:100){
  set.seed(i)
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+3*c+2*rnorm(20))
  y<-as.vector(c+2*rnorm(20,0,2))
  y <- ifelse(y >= median(y), 1, 0)
  z<-cbind(z,x)
  fit1<-lm(x~.,z)
  xhat<-fit1$fitted.values
  fit2<-lm(y~xhat)
  pp<-summary(fit2)
  pp1<-pp$coefficients[2,4]
  tt1<-deviance(fit2)
  mse1<-mean(pp$residuals^2)
  list1[[length(list1)+1]]<-tt1
  list3[[length(list3)+1]]<-pp1
  
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-b%*%t(z)+3*c+2*rnorm(20)
  y<-as.vector(c+2*rnorm(20,0,2))
  y <- ifelse(y >= median(y), 1, 0)
  z<-cbind(z,x)
  xhat<-as.numeric(predict(fit1,as.data.frame(z[,1:8])))
  ypre<-as.numeric(predict(fit2,as.data.frame(xhat)))
  list5[[length(list5)+1]]<-mean((y-ypre)^2)
}

ap3<-unlist(list3)


list8<-list()
list9<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20)
  x<-b%*%t(z)+3*c+2*rnorm(20)
  y<-as.vector(c+2*rnorm(20,0,2))
  y <- ifelse(y >= median(y), 1, 0)
  fit5 <-  cv.glmnet(y=x,x=z,keep=TRUE)
  xhat = as.numeric(fit5$fit.preval[,which(fit5$lambda==fit5$lambda.min)[1]])
  fit6 = lm(y~xhat)
  pp4 = 1
  if(nrow(coef(summary(fit6)))==2) pp4 = coef(summary(fit6))[2,4]
  list8[[length(list8)+1]]<-pp4
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-b%*%t(z)+3*c+2*rnorm(20)
  y<-as.vector(c+2*rnorm(20,0,2))
  y <- ifelse(y >= median(y), 1, 0)
  
  xhat<-as.numeric(predict(fit5 ,z, s=fit5$lambda.min))
  ppre<-as.numeric(predict(fit6,data.frame(xhat)))
  list9[[length(list9)+1]]<-mean((y-ppre)^2)
}
ap5<-unlist(list8)

lasso.mse<-as.numeric(unlist(list6))
ols.mse<-as.numeric(unlist(list5))
pv.mse<-as.numeric(unlist(list9))
v<-c(rep("OLS",100),rep("Lasso",100),rep("PV",100))
mse<-as.data.frame(cbind(v,c(ols.mse,lasso.mse,pv.mse)))
mse[,2]<- as.numeric(as.character(mse[,2]))
ggplot(mse,aes(x=mse$v,y=mse[,2],group=mse$v))+geom_boxplot(fill=c("yellow","blue","red"))+labs(x="Model",y="MSE",title="No Association between X and Y")+theme(plot.title = element_text(hjust = 0.5))

tnn<-list4>0.05
tn<-length(tnn[tnn=="TRUE"])
fp<-100-tn
tnn1<-list3>0.05
tn1<-length(tnn1[tnn1=="TRUE"])
fp1<-100-tn1
tnn2<-list8>0.05
tn2<-length(tnn2[tnn2=="TRUE"])
fp2<-100-tn2



#make MCC function
mcc<-function(a,b,c,d){
  (a*b-c*d)/sqrt((a+c)*(b+d)*(c+b)*(a+d))
}
#a=tp b=tn c=fp d=fn

mcc(tp,tn,fp,fn)
mcc(tp1,tn1,fp1,fn1)
mcc(tp2,tn2,fp2,fn2)

app<-c(ap,ap1)
class = rep(c(1,0),c(length(ap1),length(ap)))
ord = order(app)
plot(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),type = 'l',main = 'ROC curve',xlab = 'FPR',ylab = 'TPR',col="darkblue",lwd=4)
abline(0,1)
p = c(ap2,ap3)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="skyblue",type = 'l',lwd=4)
p = c(ap4,ap5)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="darkgreen",type = 'l',lwd=4)
legend('bottomright',c('TSL','TSLS','TSPL'),col = c("darkblue","skyblue","darkgreen"),lty = 1,cex=1,lwd=3)


#exnoise

dta1<-data.frame(x=1:8)
for (j in 1:8){
  list1<-list()
  for (i in 1:8){
    tl<-0.5^abs(i-j)
    list1[[length(list1)+1]]<-tl
  }
  dta1<-cbind(dta1,unlist(list1))
}
s<-matrix(c(dta1[,2],dta1[,3],dta1[,4],dta1[,5],dta1[,6],dta1[,7],dta1[,8],dta1[,9]),ncol=8,nrow=8)


list2<-list()
list4<-list()
list6<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+3*c+8*rnorm(20))
  y<-as.vector(x+c+8*rnorm(20,0,2))
  y <- ifelse(y >= median(y), 1, 0)
  crossval <-  cv.glmnet(y=x,x=z)
  penalty <- crossval$lambda.min
  fit3 <-glmnet( y=x,x=z, alpha = 1, lambda = penalty )
  xhat<-predict(crossval, newx = z, s = penalty)
  xhat<-as.vector(xhat[,1])
  fit4<-glm(y~xhat)
  pp2<-summary(fit4)
  pp3<-1
  if(nrow(coef(summary(fit4)))==2) pp3 = coef(summary(fit4))[2,4]
  tt2<-deviance(fit4)
  list2[[length(list2)+1]]<-tt2
  list4[[length(list4)+1]]<-pp3
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+3*c+8*rnorm(20))
  y<-as.vector(x+c+8*rnorm(20,0,2))
  y <- ifelse(y >= median(y), 1, 0)
  
  xhat<-as.numeric(predict(fit3 ,z, s=fit3$lambda.min))
  ypre<-as.numeric(predict(fit4,as.data.frame(as.numeric(xhat))))
  list6[[length(list6)+1]]<-mean((y-ypre)^2)
}
ap<-unlist(list4)

list1<-list()
list3<-list()
list5<-list()
for (i in 1:100){
  set.seed(i)
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+3*c+8*rnorm(20))
  y<-as.vector(x+c+8*rnorm(20,0,2))
  y <- ifelse(y >= median(y), 1, 0)
  z<-cbind(z,x)
  fit1<-lm(x~.,z)
  xhat<-fit1$fitted.values
  fit2<-lm(y~xhat)
  pp<-summary(fit2)
  pp1<-pp$coefficients[2,4]
  tt1<-deviance(fit2)
  mse1<-mean(pp$residuals^2)
  list1[[length(list1)+1]]<-tt1
  list3[[length(list3)+1]]<-pp1
  
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+3*c+8*rnorm(20))
  y<-as.vector(x+c+8*rnorm(20,0,2))
  y <- ifelse(y >= median(y), 1, 0)
  z<-cbind(z,x)
  xhat<-as.numeric(predict(fit1,as.data.frame(z[,1:8])))
  ypre<-as.numeric(predict(fit2,as.data.frame(xhat)))
  list5[[length(list5)+1]]<-mean((y-ypre)^2)
}
ap2<-unlist(list3)

list8<-list()
list9<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+3*c+8*rnorm(20))
  y<-as.vector(x+c+8*rnorm(20,0,2))
  y <- ifelse(y >= median(y), 1, 0)
  fit5 <-  cv.glmnet(y=x,x=z,keep=TRUE)
  xhat = as.numeric(fit5$fit.preval[,which(fit5$lambda==fit5$lambda.min)[1]])
  fit6 = lm(y~xhat)
  pp4 = 1
  if(nrow(coef(summary(fit6)))==2) pp4 = coef(summary(fit6))[2,4]
  list8[[length(list8)+1]]<-pp4
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+3*c+8*rnorm(20))
  y<-as.vector(x+c+8*rnorm(20,0,2))
  y <- ifelse(y >= median(y), 1, 0)
  
  xhat<-as.numeric(predict(fit5 ,z, s=fit5$lambda.min))
  ppre<-as.numeric(predict(fit6,data.frame(xhat)))
  list9[[length(list9)+1]]<-mean((y-ppre)^2)
}

ap4<-unlist(list8)

fnn<-list4>0.05
fn<-length(fnn[fnn=="TRUE"])
tp<-100-fn
fnn1<-list3>0.05
fn1<-length(fnn1[fnn1=="TRUE"])
tp1<-100-fn1
fnn2<-list8>0.05
fn2<-length(fnn2[fnn2=="TRUE"])
tp2<-100-fn2

##no association
set.seed(5)
dta1<-data.frame(x=1:8)
for (j in 1:8){
  list1<-list()
  for (i in 1:8){
    tl<-0.5^abs(i-j)
    list1[[length(list1)+1]]<-tl
  }
  dta1<-cbind(dta1,unlist(list1))
}
s<-matrix(c(dta1[,2],dta1[,3],dta1[,4],dta1[,5],dta1[,6],dta1[,7],dta1[,8],dta1[,9]),ncol=8,nrow=8)


list2<-list()
list4<-list()
list6<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20)
  x<-as.vector(b%*%t(z)+3*c+8*rnorm(20))
  y<-as.vector(c+8*rnorm(20,0,2))
  y <- ifelse(y >= median(y), 1, 0)
  crossval <-  cv.glmnet(y=x,x=z)
  penalty <- crossval$lambda.min
  fit3 <-glmnet( y=x,x=z, alpha = 1, lambda = penalty )
  xhat<-predict(crossval, newx = z, s = penalty)
  xhat<-as.vector(xhat[,1])
  fit4<-lm(y~xhat)
  pp2<-summary(fit4)
  pp3<-1
  if(nrow(coef(summary(fit4)))==2) pp3 = coef(summary(fit4))[2,4]
  tt2<-deviance(fit4)
  list2[[length(list2)+1]]<-tt2
  list4[[length(list4)+1]]<-pp3
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+3*c+8*rnorm(20))
  y<-as.vector(c+8*rnorm(20,0,2))
  y <- ifelse(y >= median(y), 1, 0)
  
  xhat<-as.numeric(predict(fit3 ,z, s=fit3$lambda.min))
  ypre<-as.numeric(predict(fit4,as.data.frame(as.numeric(xhat))))
  list6[[length(list6)+1]]<-mean((y-ypre)^2)
}
ap1<-unlist(list4)

list1<-list()
list3<-list()
list5<-list()
for (i in 1:100){
  set.seed(i)
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+3*c+8*rnorm(20))
  y<-as.vector(c+8*rnorm(20,0,2))
  y <- ifelse(y >= median(y), 1, 0)
  z<-cbind(z,x)
  fit1<-lm(x~.,z)
  xhat<-fit1$fitted.values
  fit2<-lm(y~xhat)
  pp<-summary(fit2)
  pp1<-pp$coefficients[2,4]
  tt1<-deviance(fit2)
  mse1<-mean(pp$residuals^2)
  list1[[length(list1)+1]]<-tt1
  list3[[length(list3)+1]]<-pp1
  
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+3*c+8*rnorm(20))
  y<-as.vector(c+8*rnorm(20,0,2))
  y <- ifelse(y >= median(y), 1, 0)
  z<-cbind(z,x)
  xhat<-as.numeric(predict(fit1,as.data.frame(z[,1:8])))
  ypre<-as.numeric(predict(fit2,as.data.frame(xhat)))
  list5[[length(list5)+1]]<-mean((y-ypre)^2)
}

ap3<-unlist(list3)


list8<-list()
list9<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20)
  x<-as.vector(b%*%t(z)+3*c+8*rnorm(20))
  y<-as.vector(c+8*rnorm(20,0,2))
  y <- ifelse(y >= median(y), 1, 0)
  fit5 <-  cv.glmnet(y=x,x=z,keep=TRUE)
  xhat = as.numeric(fit5$fit.preval[,which(fit5$lambda==fit5$lambda.min)[1]])
  fit6 = lm(y~xhat)
  pp4 = 1
  if(nrow(coef(summary(fit6)))==2) pp4 = coef(summary(fit6))[2,4]
  list8[[length(list8)+1]]<-pp4
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+3*c+8*rnorm(20))
  y<-as.vector(c+8*rnorm(20,0,2))
  y <- ifelse(y >= median(y), 1, 0)
  
  xhat<-as.numeric(predict(fit5 ,z, s=fit5$lambda.min))
  ppre<-as.numeric(predict(fit6,data.frame(xhat)))
  list9[[length(list9)+1]]<-mean((y-ppre)^2)
}
ap5<-unlist(list8)

tnn<-list4>0.05
tn<-length(tnn[tnn=="TRUE"])
fp<-100-tn
tnn1<-list3>0.05
tn1<-length(tnn1[tnn1=="TRUE"])
fp1<-100-tn1
tnn2<-list8>0.05
tn2<-length(tnn2[tnn2=="TRUE"])
fp2<-100-tn2



#make MCC function
mcc<-function(a,b,c,d){
  (a*b-c*d)/sqrt((a+c)*(b+d)*(c+b)*(a+d))
}
#a=tp b=tn c=fp d=fn

mcc(tp,tn,fp,fn)
mcc(tp1,tn1,fp1,fn1)
mcc(tp2,tn2,fp2,fn2)

app<-c(ap,ap1)
class = rep(c(1,0),c(length(ap1),length(ap)))
ord = order(app)
plot(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),type = 'l',main = 'ROC curve',xlab = 'FPR',ylab = 'TPR',col="darkblue",lwd=4)
abline(0,1)
p = c(ap2,ap3)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="skyblue",type = 'l',lwd=4)
p = c(ap4,ap5)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="darkgreen",type = 'l',lwd=4)
legend('bottomright',c('TSL','TSLS','TSPL'),col = c("darkblue","skyblue","darkgreen"),lty = 1,cex=1,lwd=3)


##increased strength conf lasso

dta1<-data.frame(x=1:8)
for (j in 1:8){
  list1<-list()
  for (i in 1:8){
    tl<-0.5^abs(i-j)
    list1[[length(list1)+1]]<-tl
  }
  dta1<-cbind(dta1,unlist(list1))
}
s<-matrix(c(dta1[,2],dta1[,3],dta1[,4],dta1[,5],dta1[,6],dta1[,7],dta1[,8],dta1[,9]),ncol=8,nrow=8)


list2<-list()
list4<-list()
list6<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-b%*%t(z)+10*c+2*rnorm(20)
  y<-as.vector(x+5*c+2*rnorm(20,0,2))
  y <- ifelse(y >= median(y), 1, 0)
  crossval <-  cv.glmnet(y=x,x=z)
  penalty <- crossval$lambda.min
  fit3 <-glmnet( y=x,x=z, alpha = 1, lambda = penalty )
  xhat<-predict(crossval, newx = z, s = penalty)
  xhat<-as.vector(xhat[,1])
  fit4<-glm(y~xhat)
  pp2<-summary(fit4)
  pp3<-1
  if(nrow(coef(summary(fit4)))==2) pp3 = coef(summary(fit4))[2,4]
  tt2<-deviance(fit4)
  list2[[length(list2)+1]]<-tt2
  list4[[length(list4)+1]]<-pp3
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-b%*%t(z)+10*c+2*rnorm(20)
  y<-as.vector(x+5*c+2*rnorm(20,0,2))
  y <- ifelse(y >= median(y), 1, 0)
  
  xhat<-as.numeric(predict(fit3 ,z, s=fit3$lambda.min))
  ypre<-as.numeric(predict(fit4,as.data.frame(as.numeric(xhat))))
  list6[[length(list6)+1]]<-mean((y-ypre)^2)
}
ap<-unlist(list4)

list1<-list()
list3<-list()
list5<-list()
for (i in 1:100){
  set.seed(i)
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+10*c+2*rnorm(20))
  y<-as.vector(x+5*c+2*rnorm(20,0,2))
  y <- ifelse(y >= median(y), 1, 0)
  z<-cbind(z,x)
  fit1<-lm(x~.,z)
  xhat<-fit1$fitted.values
  fit2<-lm(y~xhat)
  pp<-summary(fit2)
  pp1<-pp$coefficients[2,4]
  tt1<-deviance(fit2)
  mse1<-mean(pp$residuals^2)
  list1[[length(list1)+1]]<-tt1
  list3[[length(list3)+1]]<-pp1
  
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-b%*%t(z)+10*c+2*rnorm(20)
  y<-as.vector(x+5*c+2*rnorm(20,0,2))
  y <- ifelse(y >= median(y), 1, 0)
  z<-cbind(z,x)
  xhat<-as.numeric(predict(fit1,as.data.frame(z[,1:8])))
  ypre<-as.numeric(predict(fit2,as.data.frame(xhat)))
  list5[[length(list5)+1]]<-mean((y-ypre)^2)
}
ap2<-unlist(list3)

list8<-list()
list9<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-b%*%t(z)+10*c+2*rnorm(20)
  y<-as.vector(x+5*c+2*rnorm(20,0,2))
  y <- ifelse(y >= median(y), 1, 0)
  fit5 <-  cv.glmnet(y=x,x=z,keep=TRUE)
  xhat = as.numeric(fit5$fit.preval[,which(fit5$lambda==fit5$lambda.min)[1]])
  fit6 = lm(y~xhat)
  pp4 = 1
  if(nrow(coef(summary(fit6)))==2) pp4 = coef(summary(fit6))[2,4]
  list8[[length(list8)+1]]<-pp4
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-b%*%t(z)+10*c+2*rnorm(20)
  y<-as.vector(x+5*c+2*rnorm(20,0,2))
  y <- ifelse(y >= median(y), 1, 0)
  
  xhat<-as.numeric(predict(fit5 ,z, s=fit5$lambda.min))
  ppre<-as.numeric(predict(fit6,data.frame(xhat)))
  list9[[length(list9)+1]]<-mean((y-ppre)^2)
}

ap4<-unlist(list8)

fnn<-list4>0.05
fn<-length(fnn[fnn=="TRUE"])
tp<-100-fn
fnn1<-list3>0.05
fn1<-length(fnn1[fnn1=="TRUE"])
tp1<-100-fn1
fnn2<-list8>0.05
fn2<-length(fnn2[fnn2=="TRUE"])
tp2<-100-fn2

##no association
set.seed(5)
dta1<-data.frame(x=1:8)
for (j in 1:8){
  list1<-list()
  for (i in 1:8){
    tl<-0.5^abs(i-j)
    list1[[length(list1)+1]]<-tl
  }
  dta1<-cbind(dta1,unlist(list1))
}
s<-matrix(c(dta1[,2],dta1[,3],dta1[,4],dta1[,5],dta1[,6],dta1[,7],dta1[,8],dta1[,9]),ncol=8,nrow=8)


list2<-list()
list4<-list()
list6<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20)
  x<-b%*%t(z)+10*c+2*rnorm(20)
  y<-as.vector(5*c+2*rnorm(20,0,2))
  y <- ifelse(y >= median(y), 1, 0)
  crossval <-  cv.glmnet(y=x,x=z)
  penalty <- crossval$lambda.min
  fit3 <-glmnet( y=x,x=z, alpha = 1, lambda = penalty )
  xhat<-predict(crossval, newx = z, s = penalty)
  xhat<-as.vector(xhat[,1])
  fit4<-lm(y~xhat)
  pp2<-summary(fit4)
  pp3<-1
  if(nrow(coef(summary(fit4)))==2) pp3 = coef(summary(fit4))[2,4]
  tt2<-deviance(fit4)
  list2[[length(list2)+1]]<-tt2
  list4[[length(list4)+1]]<-pp3
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-b%*%t(z)+10*c+2*rnorm(20)
  y<-as.vector(5*c+2*rnorm(20,0,2))
  y <- ifelse(y >= median(y), 1, 0)
  
  xhat<-as.numeric(predict(fit3 ,z, s=fit3$lambda.min))
  ypre<-as.numeric(predict(fit4,as.data.frame(as.numeric(xhat))))
  list6[[length(list6)+1]]<-mean((y-ypre)^2)
}
ap1<-unlist(list4)

list1<-list()
list3<-list()
list5<-list()
for (i in 1:100){
  set.seed(i)
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+10*c+2*rnorm(20))
  y<-as.vector(5*c+2*rnorm(20,0,2))
  y <- ifelse(y >= median(y), 1, 0)
  z<-cbind(z,x)
  fit1<-lm(x~.,z)
  xhat<-fit1$fitted.values
  fit2<-lm(y~xhat)
  pp<-summary(fit2)
  pp1<-pp$coefficients[2,4]
  tt1<-deviance(fit2)
  mse1<-mean(pp$residuals^2)
  list1[[length(list1)+1]]<-tt1
  list3[[length(list3)+1]]<-pp1
  
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-b%*%t(z)+10*c+2*rnorm(20)
  y<-as.vector(5*c+2*rnorm(20,0,2))
  y <- ifelse(y >= median(y), 1, 0)
  z<-cbind(z,x)
  xhat<-as.numeric(predict(fit1,as.data.frame(z[,1:8])))
  ypre<-as.numeric(predict(fit2,as.data.frame(xhat)))
  list5[[length(list5)+1]]<-mean((y-ypre)^2)
}

ap3<-unlist(list3)


list8<-list()
list9<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20)
  x<-b%*%t(z)+10*c+2*rnorm(20)
  y<-as.vector(5*c+2*rnorm(20,0,2))
  y <- ifelse(y >= median(y), 1, 0)
  fit5 <-  cv.glmnet(y=x,x=z,keep=TRUE)
  xhat = as.numeric(fit5$fit.preval[,which(fit5$lambda==fit5$lambda.min)[1]])
  fit6 = lm(y~xhat)
  pp4 = 1
  if(nrow(coef(summary(fit6)))==2) pp4 = coef(summary(fit6))[2,4]
  list8[[length(list8)+1]]<-pp4
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-b%*%t(z)+10*c+2*rnorm(20)
  y<-as.vector(5*c+2*rnorm(20,0,2))
  y <- ifelse(y >= median(y), 1, 0)
  
  xhat<-as.numeric(predict(fit5 ,z, s=fit5$lambda.min))
  ppre<-as.numeric(predict(fit6,data.frame(xhat)))
  list9[[length(list9)+1]]<-mean((y-ppre)^2)
}
ap5<-unlist(list8)

tnn<-list4>0.05
tn<-length(tnn[tnn=="TRUE"])
fp<-100-tn
tnn1<-list3>0.05
tn1<-length(tnn1[tnn1=="TRUE"])
fp1<-100-tn1
tnn2<-list8>0.05
tn2<-length(tnn2[tnn2=="TRUE"])
fp2<-100-tn2



#make MCC function
mcc<-function(a,b,c,d){
  (a*b-c*d)/sqrt((a+c)*(b+d)*(c+b)*(a+d))
}
#a=tp b=tn c=fp d=fn

mcc(tp,tn,fp,fn)
mcc(tp1,tn1,fp1,fn1)
mcc(tp2,tn2,fp2,fn2)

app<-c(ap,ap1)
class = rep(c(1,0),c(length(ap1),length(ap)))
ord = order(app)
plot(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),type = 'l',main = 'ROC curve',xlab = 'FPR',ylab = 'TPR',col="darkblue",lwd=4)
abline(0,1)
p = c(ap2,ap3)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="skyblue",type = 'l',lwd=4)
p = c(ap4,ap5)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="darkgreen",type = 'l',lwd=4)
legend('bottomright',c('TSL','TSLS','TSPL'),col = c("darkblue","skyblue","darkgreen"),lty = 1,cex=1,lwd=3)


#negass
dta1<-data.frame(x=1:8)
for (j in 1:8){
  list1<-list()
  for (i in 1:8){
    tl<-0.5^abs(i-j)
    list1[[length(list1)+1]]<-tl
  }
  dta1<-cbind(dta1,unlist(list1))
}
s<-matrix(c(dta1[,2],dta1[,3],dta1[,4],dta1[,5],dta1[,6],dta1[,7],dta1[,8],dta1[,9]),ncol=8,nrow=8)


list2<-list()
list4<-list()
list6<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(1,1)
  x<-g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20)
  y<-as.vector(-0.25*x+c+3*rnorm(20))
  y <- ifelse(y >= median(y), 1, 0)
  crossval <-  cv.glmnet(y=x,x=z)
  penalty <- crossval$lambda.min
  fit3 <-glmnet( y=x,x=z, alpha = 1, lambda = penalty )
  xhat<-predict(crossval, newx = z, s = penalty)
  xhat<-as.vector(xhat[,1])
  fit4<-lm(y~xhat)
  pp2<-summary(fit4)
  pp3<-1
  if(nrow(coef(summary(fit4)))==2) pp3 = coef(summary(fit4))[2,4]
  tt2<-deviance(fit4)
  list2[[length(list2)+1]]<-tt2
  list4[[length(list4)+1]]<-pp3
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(1,1)
  x<-g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20)
  y<-as.vector(-0.25*x+c+3*rnorm(20))
  y <- ifelse(y >= median(y), 1, 0)
  
  xhat<-as.numeric(predict(fit3 ,z, s=fit3$lambda.min))
  ypre<-as.numeric(predict(fit4,as.data.frame(as.numeric(xhat))))
  list6[[length(list6)+1]]<-mean((y-ypre)^2)
}
ap<-unlist(list4)

list1<-list()
list3<-list()
list5<-list()
for (i in 1:100){
  set.seed(i)
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(1,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(-0.25*x+c+3*rnorm(20))
  y <- ifelse(y >= median(y), 1, 0)
  z<-cbind(z,x)
  fit1<-lm(x~.,z)
  xhat<-fit1$fitted.values
  fit2<-lm(y~xhat)
  pp<-summary(fit2)
  pp1<-pp$coefficients[2,4]
  tt1<-deviance(fit2)
  mse1<-mean(pp$residuals^2)
  list1[[length(list1)+1]]<-tt1
  list3[[length(list3)+1]]<-pp1
  
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(1,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(-0.25*x+c+3*rnorm(20))
  y <- ifelse(y >= median(y), 1, 0)
  z<-cbind(z,x)
  xhat<-as.numeric(predict(fit1,as.data.frame(z[,1:8])))
  ypre<-as.numeric(predict(fit2,as.data.frame(xhat)))
  list5[[length(list5)+1]]<-mean((y-ypre)^2)
}
ap2<-unlist(list3)

list8<-list()
list9<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(1,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(-0.25*x+c+3*rnorm(20))
  y <- ifelse(y >= median(y), 1, 0)
  fit5 <-  cv.glmnet(y=x,x=z,keep=TRUE)
  xhat = as.numeric(fit5$fit.preval[,which(fit5$lambda==fit5$lambda.min)[1]])
  fit6 = lm(y~xhat)
  pp4 = 1
  if(nrow(coef(summary(fit6)))==2) pp4 = coef(summary(fit6))[2,4]
  list8[[length(list8)+1]]<-pp4
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(1,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(-0.25*x+c+3*rnorm(20))
  y <- ifelse(y >= median(y), 1, 0)
  
  xhat<-as.numeric(predict(fit5 ,z, s=fit5$lambda.min))
  ppre<-as.numeric(predict(fit6,data.frame(xhat)))
  list9[[length(list9)+1]]<-mean((y-ppre)^2)
}

ap4<-unlist(list8)

lasso.mse<-as.numeric(unlist(list6))
ols.mse<-as.numeric(unlist(list5))
pv.mse<-as.numeric(unlist(list9))
v<-c(rep("OLS",100),rep("Lasso",100),rep("PV",100))
mse<-as.data.frame(cbind(v,c(ols.mse,lasso.mse,pv.mse)))
mse[,2]<- as.numeric(as.character(mse[,2]))
ggplot(mse,aes(x=mse$v,y=mse[,2],group=mse$v))+geom_boxplot(fill=c("yellow","blue","red"))+labs(x="Model",y="MSE",title="Association between X and Y")+theme(plot.title = element_text(hjust = 0.5))

fnn<-list4>0.05
fn<-length(fnn[fnn=="TRUE"])
tp<-100-fn
fnn1<-list3>0.05
fn1<-length(fnn1[fnn1=="TRUE"])
tp1<-100-fn1
fnn2<-list8>0.05
fn2<-length(fnn2[fnn2=="TRUE"])
tp2<-100-fn2

##no association
set.seed(5)
dta1<-data.frame(x=1:8)
for (j in 1:8){
  list1<-list()
  for (i in 1:8){
    tl<-0.5^abs(i-j)
    list1[[length(list1)+1]]<-tl
  }
  dta1<-cbind(dta1,unlist(list1))
}
s<-matrix(c(dta1[,2],dta1[,3],dta1[,4],dta1[,5],dta1[,6],dta1[,7],dta1[,8],dta1[,9]),ncol=8,nrow=8)


list2<-list()
list4<-list()
list6<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(1,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(c+3*rnorm(20))
  y <- ifelse(y >= median(y), 1, 0)
  crossval <-  cv.glmnet(y=x,x=z)
  penalty <- crossval$lambda.min
  fit3 <-glmnet( y=x,x=z, alpha = 1, lambda = penalty )
  xhat<-predict(crossval, newx = z, s = penalty)
  xhat<-as.vector(xhat[,1])
  fit4<-lm(y~xhat)
  pp2<-summary(fit4)
  pp3<-1
  if(nrow(coef(summary(fit4)))==2) pp3 = coef(summary(fit4))[2,4]
  tt2<-deviance(fit4)
  list2[[length(list2)+1]]<-tt2
  list4[[length(list4)+1]]<-pp3
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(1,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(c+3*rnorm(20))
  y <- ifelse(y >= median(y), 1, 0)
  
  xhat<-as.numeric(predict(fit3 ,z, s=fit3$lambda.min))
  ypre<-as.numeric(predict(fit4,as.data.frame(as.numeric(xhat))))
  list6[[length(list6)+1]]<-mean((y-ypre)^2)
}
ap1<-unlist(list4)

list1<-list()
list3<-list()
list5<-list()
for (i in 1:100){
  set.seed(i)
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(1,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(c+3*rnorm(20))
  y <- ifelse(y >= median(y), 1, 0)
  z<-cbind(z,x)
  fit1<-lm(x~.,z)
  xhat<-fit1$fitted.values
  fit2<-lm(y~xhat)
  pp<-summary(fit2)
  pp1<-pp$coefficients[2,4]
  tt1<-deviance(fit2)
  mse1<-mean(pp$residuals^2)
  list1[[length(list1)+1]]<-tt1
  list3[[length(list3)+1]]<-pp1
  
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(1,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(c+3*rnorm(20))
  y <- ifelse(y >= median(y), 1, 0)
  z<-cbind(z,x)
  xhat<-as.numeric(predict(fit1,as.data.frame(z[,1:8])))
  ypre<-as.numeric(predict(fit2,as.data.frame(xhat)))
  list5[[length(list5)+1]]<-mean((y-ypre)^2)
}

ap3<-unlist(list3)


list8<-list()
list9<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(1,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(c+3*rnorm(20))
  y <- ifelse(y >= median(y), 1, 0)
  fit5 <-  cv.glmnet(y=x,x=z,keep=TRUE)
  xhat = as.numeric(fit5$fit.preval[,which(fit5$lambda==fit5$lambda.min)[1]])
  fit6 = lm(y~xhat)
  pp4 = 1
  if(nrow(coef(summary(fit6)))==2) pp4 = coef(summary(fit6))[2,4]
  list8[[length(list8)+1]]<-pp4
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(1,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(c+3*rnorm(20))
  y <- ifelse(y >= median(y), 1, 0)
  
  xhat<-as.numeric(predict(fit5 ,z, s=fit5$lambda.min))
  ppre<-as.numeric(predict(fit6,data.frame(xhat)))
  list9[[length(list9)+1]]<-mean((y-ppre)^2)
}
ap5<-unlist(list8)

lasso.mse<-as.numeric(unlist(list6))
ols.mse<-as.numeric(unlist(list5))
pv.mse<-as.numeric(unlist(list9))
v<-c(rep("OLS",100),rep("Lasso",100),rep("PV",100))
mse<-as.data.frame(cbind(v,c(ols.mse,lasso.mse,pv.mse)))
mse[,2]<- as.numeric(as.character(mse[,2]))
ggplot(mse,aes(x=mse$v,y=mse[,2],group=mse$v))+geom_boxplot(fill=c("yellow","blue","red"))+labs(x="Model",y="MSE",title="No Association between X and Y")+theme(plot.title = element_text(hjust = 0.5))

tnn<-list4>0.05
tn<-length(tnn[tnn=="TRUE"])
fp<-100-tn
tnn1<-list3>0.05
tn1<-length(tnn1[tnn1=="TRUE"])
fp1<-100-tn1
tnn2<-list8>0.05
tn2<-length(tnn2[tnn2=="TRUE"])
fp2<-100-tn2



#make MCC function
mcc<-function(a,b,c,d){
  (a*b-c*d)/sqrt((a+c)*(b+d)*(c+b)*(a+d))
}
#a=tp b=tn c=fp d=fn

mcc(tp,tn,fp,fn)
mcc(tp1,tn1,fp1,fn1)
mcc(tp2,tn2,fp2,fn2)

app<-c(ap,ap1)
class = rep(c(1,0),c(length(ap1),length(ap)))
ord = order(app)
plot(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),type = 'l',main = 'ROC curve',xlab = 'FPR',ylab = 'TPR',col="darkblue",lwd=4)
abline(0,1)
p = c(ap2,ap3)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="skyblue",type = 'l',lwd=4)
p = c(ap4,ap5)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="darkgreen",type = 'l',lwd=4)
legend('bottomright',c('TSL','TSLS','TSPL'),col = c("darkblue","skyblue","darkgreen"),lty = 1,cex=1,lwd=3)


#liml sims

dta1<-data.frame(x=1:8)
for (j in 1:8){
  list1<-list()
  for (i in 1:8){
    tl<-0.5^abs(i-j)
    list1[[length(list1)+1]]<-tl
  }
  dta1<-cbind(dta1,unlist(list1))
}
s<-matrix(c(dta1[,2],dta1[,3],dta1[,4],dta1[,5],dta1[,6],dta1[,7],dta1[,8],dta1[,9]),ncol=8,nrow=8)

list1<-list()
list3<-list()
for(i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+3*c+2*rnorm(20))
  y<-as.vector(x+c+2*rnorm(20,0,2))
  
  p<-ivmodel(Y=y,D=x,Z=z)
  ttp<-LIML(p)
  limlp<-ttp$p.value
  list1[[length(list1)+1]]<-limlp
  tts<-p$kClass
  tts1<-tts$p.value[2]
  list3[[length(list3)+1]]<-tts1
}
ap8<-unlist(list1)
ap10<-unlist(list3)
fnn3<-list1>0.05
fn3<-length(fnn3[fnn3=="TRUE"])
tp3<-100-fn3

list2<-list()
list4<-list()
list6<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-b%*%t(z)+3*c+2*rnorm(20)
  y<-as.vector(x+c+2*rnorm(20,0,2))
  crossval <-  cv.glmnet(y=x,x=z)
  penalty <- crossval$lambda.min
  fit3 <-glmnet( y=x,x=z, alpha = 1, lambda = penalty )
  xhat<-predict(crossval, newx = z, s = penalty)
  xhat<-as.vector(xhat[,1])
  fit4<-lm(y~xhat)
  pp2<-summary(fit4)
  pp3<-1
  if(nrow(coef(summary(fit4)))==2) pp3 = coef(summary(fit4))[2,4]
  tt2<-deviance(fit4)
  list2[[length(list2)+1]]<-tt2
  list4[[length(list4)+1]]<-pp3
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-b%*%t(z)+3*c+2*rnorm(20)
  y<-as.vector(x+c+2*rnorm(20,0,2))
  
  xhat<-as.numeric(predict(fit3 ,z, s=fit3$lambda.min))
  ypre<-as.numeric(predict(fit4,as.data.frame(as.numeric(xhat))))
  list6[[length(list6)+1]]<-mean((y-ypre)^2)
}
ap<-unlist(list4)

list1<-list()
list3<-list()
list5<-list()
for (i in 1:100){
  set.seed(i)
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+3*c+2*rnorm(20))
  y<-as.vector(x+c+2*rnorm(20,0,2))
  z<-cbind(z,x)
  fit1<-lm(x~.,z)
  xhat<-fit1$fitted.values
  fit2<-lm(y~xhat)
  pp<-summary(fit2)
  pp1<-pp$coefficients[2,4]
  tt1<-deviance(fit2)
  mse1<-mean(pp$residuals^2)
  list1[[length(list1)+1]]<-tt1
  list3[[length(list3)+1]]<-pp1
  
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-b%*%t(z)+3*c+2*rnorm(20)
  y<-as.vector(x+c+2*rnorm(20,0,2))
  z<-cbind(z,x)
  xhat<-as.numeric(predict(fit1,as.data.frame(z[,1:8])))
  ypre<-as.numeric(predict(fit2,as.data.frame(xhat)))
  list5[[length(list5)+1]]<-mean((y-ypre)^2)
}
ap2<-unlist(list3)

list8<-list()
list9<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-b%*%t(z)+3*c+2*rnorm(20)
  y<-as.vector(x+c+2*rnorm(20,0,2))
  fit5 <-  cv.glmnet(y=x,x=z,keep=TRUE)
  xhat = as.numeric(fit5$fit.preval[,which(fit5$lambda==fit5$lambda.min)[1]])
  fit6 = lm(y~xhat)
  pp4 = 1
  if(nrow(coef(summary(fit6)))==2) pp4 = coef(summary(fit6))[2,4]
  list8[[length(list8)+1]]<-pp4
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-b%*%t(z)+3*c+2*rnorm(20)
  y<-as.vector(x+c+2*rnorm(20,0,2))
  
  xhat<-as.numeric(predict(fit5 ,z, s=fit5$lambda.min))
  ppre<-as.numeric(predict(fit6,data.frame(xhat)))
  list9[[length(list9)+1]]<-mean((y-ppre)^2)
}

ap4<-unlist(list8)

lasso.mse<-as.numeric(unlist(list6))
ols.mse<-as.numeric(unlist(list5))
pv.mse<-as.numeric(unlist(list9))
v<-c(rep("TSLS",100),rep("TSL",100),rep("TSPL",100))
mse<-as.data.frame(cbind(v,c(ols.mse,lasso.mse,pv.mse)))
mse[,2]<- as.numeric(as.character(mse[,2]))
ggplot(mse,aes(x=mse$v,y=mse[,2],group=mse$v))+geom_boxplot(fill=c("yellow","blue","red"))+labs(x="Model",y="MSE",title="Association between X and Y")+theme(plot.title = element_text(hjust = 0.5))

fnn<-list4>0.05
fn<-length(fnn[fnn=="TRUE"])
tp<-100-fn
fnn1<-list3>0.05
fn1<-length(fnn1[fnn1=="TRUE"])
tp1<-100-fn1
fnn2<-list8>0.05
fn2<-length(fnn2[fnn2=="TRUE"])
tp2<-100-fn2

##no association
dta1<-data.frame(x=1:8)
for (j in 1:8){
  list1<-list()
  for (i in 1:8){
    tl<-0.5^abs(i-j)
    list1[[length(list1)+1]]<-tl
  }
  dta1<-cbind(dta1,unlist(list1))
}
s<-matrix(c(dta1[,2],dta1[,3],dta1[,4],dta1[,5],dta1[,6],dta1[,7],dta1[,8],dta1[,9]),ncol=8,nrow=8)


list2<-list()
list4<-list()
list6<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20)
  x<-b%*%t(z)+3*c+2*rnorm(20)
  y<-as.vector(c+2*rnorm(20,0,2))
  crossval <-  cv.glmnet(y=x,x=z)
  penalty <- crossval$lambda.min
  fit3 <-glmnet( y=x,x=z, alpha = 1, lambda = penalty )
  xhat<-predict(crossval, newx = z, s = penalty)
  xhat<-as.vector(xhat[,1])
  fit4<-lm(y~xhat)
  pp2<-summary(fit4)
  pp3<-1
  if(nrow(coef(summary(fit4)))==2) pp3 = coef(summary(fit4))[2,4]
  tt2<-deviance(fit4)
  list2[[length(list2)+1]]<-tt2
  list4[[length(list4)+1]]<-pp3
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-b%*%t(z)+3*c+2*rnorm(20)
  y<-as.vector(c+2*rnorm(20,0,2))
  
  xhat<-as.numeric(predict(fit3 ,z, s=fit3$lambda.min))
  ypre<-as.numeric(predict(fit4,as.data.frame(as.numeric(xhat))))
  list6[[length(list6)+1]]<-mean((y-ypre)^2)
}
ap1<-unlist(list4)

list1<-list()
list3<-list()
list5<-list()
for (i in 1:100){
  set.seed(i)
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+3*c+2*rnorm(20))
  y<-as.vector(c+2*rnorm(20,0,2))
  z<-cbind(z,x)
  fit1<-lm(x~.,z)
  xhat<-fit1$fitted.values
  fit2<-lm(y~xhat)
  pp<-summary(fit2)
  pp1<-pp$coefficients[2,4]
  tt1<-deviance(fit2)
  mse1<-mean(pp$residuals^2)
  list1[[length(list1)+1]]<-tt1
  list3[[length(list3)+1]]<-pp1
  
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-b%*%t(z)+3*c+2*rnorm(20)
  y<-as.vector(c+2*rnorm(20,0,2))
  z<-cbind(z,x)
  xhat<-as.numeric(predict(fit1,as.data.frame(z[,1:8])))
  ypre<-as.numeric(predict(fit2,as.data.frame(xhat)))
  list5[[length(list5)+1]]<-mean((y-ypre)^2)
}

ap3<-unlist(list3)


list8<-list()
list9<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20)
  x<-b%*%t(z)+3*c+2*rnorm(20)
  y<-as.vector(c+2*rnorm(20,0,2))
  fit5 <-  cv.glmnet(y=x,x=z,keep=TRUE)
  xhat = as.numeric(fit5$fit.preval[,which(fit5$lambda==fit5$lambda.min)[1]])
  fit6 = lm(y~xhat)
  pp4 = 1
  if(nrow(coef(summary(fit6)))==2) pp4 = coef(summary(fit6))[2,4]
  list8[[length(list8)+1]]<-pp4
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-b%*%t(z)+3*c+2*rnorm(20)
  y<-as.vector(c+2*rnorm(20,0,2))
  
  xhat<-as.numeric(predict(fit5 ,z, s=fit5$lambda.min))
  ppre<-as.numeric(predict(fit6,data.frame(xhat)))
  list9[[length(list9)+1]]<-mean((y-ppre)^2)
}
ap5<-unlist(list8)


list2<-list()
list4<-list()
for(i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+3*c+2*rnorm(20))
  y<-as.vector(c+2*rnorm(20,0,2))
  
  p<-ivmodel(Y=y,D=x,Z=z)
  ttp<-LIML(p)
  limlp<-ttp$p.value
  list2[[length(list2)+1]]<-limlp
  tts3<-p$kClass
  tts2<-tts3$p.value[2]
  list4[[length(list4)+1]]<-tts2
}
ap9<-unlist(list2)
ap11<-unlist(list4)

tnn3<-list2>0.05
tn3<-length(tnn3[tnn3=="TRUE"])
fp3<-100-tn3

tnn<-list4>0.05
tn<-length(tnn[tnn=="TRUE"])
fp<-100-tn
tnn1<-list3>0.05
tn1<-length(tnn1[tnn1=="TRUE"])
fp1<-100-tn1
tnn2<-list8>0.05
tn2<-length(tnn2[tnn2=="TRUE"])
fp2<-100-tn2



#make MCC function
mcc<-function(a,b,c,d){
  (a*b-c*d)/sqrt((a+c)*(b+d)*(c+b)*(a+d))
}
#a=tp b=tn c=fp d=fn

mcc(tp,tn,fp,fn)
mcc(tp1,tn1,fp1,fn1)
mcc(tp2,tn2,fp2,fn2)
mcc(tp3,tn3,fp3,fn3)

#ROC plot
app<-c(ap,ap1)
class = rep(c(1,0),c(length(ap1),length(ap)))
ord = order(app)
plot(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),type = 'l',main = 'ROC curve',xlab = 'FPR',ylab = 'TPR',col="darkblue",lwd=4)
abline(0,1)
p = c(ap2,ap3)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="skyblue",type = 'l',lwd=4)
p = c(ap4,ap5)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="darkgreen",type = 'l',lwd=4)
p = c(ap8,ap9)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="palegreen",type = 'l',lwd=4)
legend('bottomright',c('TSL','TSLS','TSPL','LIML'),col = c("darkblue","skyblue","darkgreen","palegreen"),lty = 1,cex=1,lwd=3)


##extra noise
dta1<-data.frame(x=1:8)
for (j in 1:8){
  list1<-list()
  for (i in 1:8){
    tl<-0.5^abs(i-j)
    list1[[length(list1)+1]]<-tl
  }
  dta1<-cbind(dta1,unlist(list1))
}
s<-matrix(c(dta1[,2],dta1[,3],dta1[,4],dta1[,5],dta1[,6],dta1[,7],dta1[,8],dta1[,9]),ncol=8,nrow=8)

list1<-list()
for(i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+3*c+8*rnorm(20))
  y<-as.vector(x+c+8*rnorm(20,0,2))
  
  p<-ivmodel(Y=y,D=x,Z=z)
  ttp<-LIML(p)
  limlp<-ttp$p.value
  list1[[length(list1)+1]]<-limlp
}
ap8<-unlist(list1)
fnn3<-list1>0.05
fn3<-length(fnn3[fnn3=="TRUE"])
tp3<-100-fn3

list2<-list()
for(i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+3*c+8*rnorm(20))
  y<-as.vector(c+8*rnorm(20,0,2))
  
  p<-ivmodel(Y=y,D=x,Z=z)
  ttp<-LIML(p)
  limlp<-ttp$p.value
  list2[[length(list2)+1]]<-limlp
}
ap9<-unlist(list2)

tnn3<-list2>0.05
tn3<-length(tnn3[tnn3=="TRUE"])
fp3<-100-tn3


list2<-list()
list4<-list()
list6<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+3*c+8*rnorm(20))
  y<-as.vector(x+c+8*rnorm(20,0,2))
  crossval <-  cv.glmnet(y=x,x=z)
  penalty <- crossval$lambda.min
  fit3 <-glmnet( y=x,x=z, alpha = 1, lambda = penalty )
  xhat<-predict(crossval, newx = z, s = penalty)
  xhat<-as.vector(xhat[,1])
  fit4<-lm(y~xhat)
  pp2<-summary(fit4)
  pp3<-1
  if(nrow(coef(summary(fit4)))==2) pp3 = coef(summary(fit4))[2,4]
  tt2<-deviance(fit4)
  list2[[length(list2)+1]]<-tt2
  list4[[length(list4)+1]]<-pp3
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+3*c+8*rnorm(20))
  y<-as.vector(x+c+8*rnorm(20,0,2))
  xhat<-as.numeric(predict(fit3 ,z, s=fit3$lambda.min))
  ypre<-as.numeric(predict(fit4,as.data.frame(as.numeric(xhat))))
  list6[[length(list6)+1]]<-mean((y-ypre)^2)
}
ap<-unlist(list4)

list1<-list()
list3<-list()
list5<-list()
for (i in 1:100){
  set.seed(i)
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+3*c+8*rnorm(20))
  y<-as.vector(x+c+8*rnorm(20,0,2))
  z<-cbind(z,x)
  fit1<-lm(x~.,z)
  xhat<-fit1$fitted.values
  fit2<-lm(y~xhat)
  pp<-summary(fit2)
  pp1<-pp$coefficients[2,4]
  tt1<-deviance(fit2)
  mse1<-mean(pp$residuals^2)
  list1[[length(list1)+1]]<-tt1
  list3[[length(list3)+1]]<-pp1
  
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+3*c+8*rnorm(20))
  y<-as.vector(x+c+8*rnorm(20,0,2))
  z<-cbind(z,x)
  xhat<-as.numeric(predict(fit1,as.data.frame(z[,1:8])))
  ypre<-as.numeric(predict(fit2,as.data.frame(xhat)))
  list5[[length(list5)+1]]<-mean((y-ypre)^2)
}
ap2<-unlist(list3)

list8<-list()
list9<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+3*c+8*rnorm(20))
  y<-as.vector(x+c+8*rnorm(20,0,2))
  fit5 <-  cv.glmnet(y=x,x=z,keep=TRUE)
  xhat = as.numeric(fit5$fit.preval[,which(fit5$lambda==fit5$lambda.min)[1]])
  fit6 = lm(y~xhat)
  pp4 = 1
  if(nrow(coef(summary(fit6)))==2) pp4 = coef(summary(fit6))[2,4]
  list8[[length(list8)+1]]<-pp4
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+3*c+8*rnorm(20))
  y<-as.vector(x+c+8*rnorm(20,0,2))
  
  xhat<-as.numeric(predict(fit5 ,z, s=fit5$lambda.min))
  ppre<-as.numeric(predict(fit6,data.frame(xhat)))
  list9[[length(list9)+1]]<-mean((y-ppre)^2)
}
ap4<-unlist(list8)

fnn<-list4>0.05
fn<-length(fnn[fnn=="TRUE"])
tp<-100-fn
fnn1<-list3>0.05
fn1<-length(fnn1[fnn1=="TRUE"])
tp1<-100-fn1
fnn2<-list8>0.05
fn2<-length(fnn2[fnn2=="TRUE"])
tp2<-100-fn2

dta1<-data.frame(x=1:8)
for (j in 1:8){
  list1<-list()
  for (i in 1:8){
    tl<-0.5^abs(i-j)
    list1[[length(list1)+1]]<-tl
  }
  dta1<-cbind(dta1,unlist(list1))
}
s<-matrix(c(dta1[,2],dta1[,3],dta1[,4],dta1[,5],dta1[,6],dta1[,7],dta1[,8],dta1[,9]),ncol=8,nrow=8)


list2<-list()
list4<-list()
list6<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+3*c+8*rnorm(20))
  y<-as.vector(c+8*rnorm(20,0,2))
  crossval <-  cv.glmnet(y=x,x=z)
  penalty <- crossval$lambda.min
  fit3 <-glmnet( y=x,x=z, alpha = 1, lambda = penalty )
  xhat<-predict(crossval, newx = z, s = penalty)
  xhat<-as.vector(xhat[,1])
  fit4<-lm(y~xhat)
  pp2<-summary(fit4)
  pp3<-1
  if(nrow(coef(summary(fit4)))==2) pp3 = coef(summary(fit4))[2,4]
  tt2<-deviance(fit4)
  list2[[length(list2)+1]]<-tt2
  list4[[length(list4)+1]]<-pp3
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+3*c+8*rnorm(20))
  y<-as.vector(c+8*rnorm(20,0,2))
  xhat<-as.numeric(predict(fit3 ,z, s=fit3$lambda.min))
  ypre<-as.numeric(predict(fit4,as.data.frame(as.numeric(xhat))))
  list6[[length(list6)+1]]<-mean((y-ypre)^2)
}
ap1<-unlist(list4)

list1<-list()
list3<-list()
list5<-list()
for (i in 1:100){
  set.seed(i)
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+3*c+8*rnorm(20))
  y<-as.vector(c+8*rnorm(20,0,2))
  z<-cbind(z,x)
  fit1<-lm(x~.,z)
  xhat<-fit1$fitted.values
  fit2<-lm(y~xhat)
  pp<-summary(fit2)
  pp1<-pp$coefficients[2,4]
  tt1<-deviance(fit2)
  mse1<-mean(pp$residuals^2)
  list1[[length(list1)+1]]<-tt1
  list3[[length(list3)+1]]<-pp1
  
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+3*c+8*rnorm(20))
  y<-as.vector(c+8*rnorm(20,0,2))
  z<-cbind(z,x)
  xhat<-as.numeric(predict(fit1,as.data.frame(z[,1:8])))
  ypre<-as.numeric(predict(fit2,as.data.frame(xhat)))
  list5[[length(list5)+1]]<-mean((y-ypre)^2)
}
ap3<-unlist(list3)


list8<-list()
list9<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+3*c+8*rnorm(20))
  y<-as.vector(c+8*rnorm(20,0,2))
  fit5 <-  cv.glmnet(y=x,x=z,keep=TRUE)
  xhat = as.numeric(fit5$fit.preval[,which(fit5$lambda==fit5$lambda.min)[1]])
  fit6 = lm(y~xhat)
  pp4 = 1
  if(nrow(coef(summary(fit6)))==2) pp4 = coef(summary(fit6))[2,4]
  list8[[length(list8)+1]]<-pp4
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+3*c+8*rnorm(20))
  y<-as.vector(c+8*rnorm(20,0,2))
  
  xhat<-as.numeric(predict(fit5 ,z, s=fit5$lambda.min))
  ppre<-as.numeric(predict(fit6,data.frame(xhat)))
  list9[[length(list9)+1]]<-mean((y-ppre)^2)
}
ap5<-unlist(list8)

tnn<-list4>0.05
tn<-length(tnn[tnn=="TRUE"])
fp<-100-tn
tnn1<-list3>0.05
tn1<-length(tnn1[tnn1=="TRUE"])
fp1<-100-tn1
tnn2<-list8>0.05
tn2<-length(tnn2[tnn2=="TRUE"])
fp2<-100-tn2

mcc(tp,tn,fp,fn)
mcc(tp1,tn1,fp1,fn1)
mcc(tp2,tn2,fp2,fn2)

#ROC plot
#ROC plot
app<-c(ap,ap1)
class = rep(c(1,0),c(length(ap1),length(ap)))
ord = order(app)
plot(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),type = 'l',main = 'ROC curve',xlab = 'FPR',ylab = 'TPR',col="darkblue",lwd=4)
abline(0,1)
p = c(ap2,ap3)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="skyblue",type = 'l',lwd=4)
p = c(ap4,ap5)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="darkgreen",type = 'l',lwd=4)
p = c(ap8,ap9)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="palegreen",type = 'l',lwd=4)
legend('bottomright',c('TSL','TSLS','TSPL','LIML'),col = c("darkblue","skyblue","darkgreen","palegreen"),lty = 1,cex=1,lwd=3)


#More Ivs
dta1<-data.frame(x=1:40)
for (j in 1:40){
  list1<-list()
  for (i in 1:40){
    tl<-0.5^abs(i-j)
    list1[[length(list1)+1]]<-tl
  }
  dta1<-cbind(dta1,unlist(list1))
}
s<-matrix(c(dta1[,2],dta1[,3],dta1[,4],dta1[,5],dta1[,6],dta1[,7],dta1[,8],dta1[,9],dta1[,10],dta1[,11],dta1[,12],dta1[,13],dta1[,14],dta1[,15],dta1[,16],dta1[,17],dta1[,18],dta1[,19],dta1[,20],dta1[,21],dta1[,22],dta1[,23],dta1[,24],dta1[,25],dta1[,26],dta1[,27],dta1[,28],dta1[,29],dta1[,30],dta1[,31],dta1[,32],dta1[,33],dta1[,34],dta1[,35],dta1[,36],dta1[,37],dta1[,38],dta1[,39],dta1[,40],dta1[,41]),ncol=40,nrow=40)

list1<-list()
for(i in 1:100){
  set.seed(i)
  z<-mvrnorm(100,rep(2,40),s)
  b<-(c(3,1.5,0,0,2,0,0,0,2,0,0,1,2,5,0,0,0,0,1,1,0,0,0,0,1,2,0,0,4,2,1,2,0,0,0,7,0,0,0,1))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(x+c+2*rnorm(100,0,6))
  
  p<-ivmodel(Y=y,D=x,Z=z)
  ttp<-LIML(p)
  limlp<-ttp$p.value
  list1[[length(list1)+1]]<-limlp
}
ap8<-unlist(list1)
fnn3<-list1>0.05
fn3<-length(fnn3[fnn3=="TRUE"])
tp3<-100-fn3

list2<-list()
for(i in 1:100){
  set.seed(i)
  z<-mvrnorm(100,rep(2,40),s)
  b<-(c(3,1.5,0,0,2,0,0,0,2,0,0,1,2,5,0,0,0,0,1,1,0,0,0,0,1,2,0,0,4,2,1,2,0,0,0,7,0,0,0,1))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(c+2*rnorm(100,0,6))
  
  p<-ivmodel(Y=y,D=x,Z=z)
  ttp<-LIML(p)
  limlp<-ttp$p.value
  list2[[length(list2)+1]]<-limlp
}
ap9<-unlist(list2)

tnn3<-list2>0.05
tn3<-length(tnn3[tnn3=="TRUE"])
fp3<-100-tn3

list2<-list()
list4<-list()
list6<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(100,rep(2,40),s)
  b<-(c(3,1.5,0,0,2,0,0,0,2,0,0,1,2,5,0,0,0,0,1,1,0,0,0,0,1,2,0,0,4,2,1,2,0,0,0,7,0,0,0,1))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(x+c+2*rnorm(100,0,6))
  crossval <-  cv.glmnet(y=x,x=z)
  penalty <- crossval$lambda.min
  fit3 <-glmnet( y=x,x=z, alpha = 1, lambda = penalty )
  xhat<-predict(crossval, newx = z, s = penalty)
  xhat<-as.vector(xhat[,1])
  fit4<-lm(y~xhat)
  pp2<-summary(fit4)
  pp3<-1
  if(nrow(coef(summary(fit4)))==2) pp3 = coef(summary(fit4))[2,4]
  tt2<-deviance(fit4)
  list2[[length(list2)+1]]<-tt2
  list4[[length(list4)+1]]<-pp3
  
  z<-mvrnorm(100,rep(2,40),s)
  b<-(c(3,1.5,0,0,2,0,0,0,2,0,0,1,2,5,0,0,0,0,1,1,0,0,0,0,1,2,0,0,4,2,1,2,0,0,0,7,0,0,0,1))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(x+c+2*rnorm(100,0,6))
  
  xhat<-as.numeric(predict(fit3 ,z, s=fit3$lambda.min))
  ypre<-as.numeric(predict(fit4,as.data.frame(as.numeric(xhat))))
  list6[[length(list6)+1]]<-mean((y-ypre)^2)
}
ap<-unlist(list4)

list1<-list()
list3<-list()
list5<-list()
for (i in 1:100){
  set.seed(i)
  z<-as.data.frame(mvrnorm(100,rep(2,40),s))
  b<-(c(3,1.5,0,0,2,0,0,0,2,0,0,1,2,5,0,0,0,0,1,1,0,0,0,0,1,2,0,0,4,2,1,2,0,0,0,7,0,0,0,1))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(x+c+2*rnorm(100,0,6))
  z<-cbind(z,x)
  fit1<-lm(x~.,z)
  xhat<-fit1$fitted.values
  fit2<-lm(y~xhat)
  pp<-summary(fit2)
  pp1<-pp$coefficients[2,4]
  tt1<-deviance(fit2)
  mse1<-mean(pp$residuals^2)
  list1[[length(list1)+1]]<-tt1
  list3[[length(list3)+1]]<-pp1
  
  z<-as.data.frame(mvrnorm(100,rep(2,40),s))
  b<-(c(3,1.5,0,0,2,0,0,0,2,0,0,1,2,5,0,0,0,0,1,1,0,0,0,0,1,2,0,0,4,2,1,2,0,0,0,7,0,0,0,1))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(x+c+2*rnorm(100,0,6))
  z<-cbind(z,x)
  xhat<-as.numeric(predict(fit1,as.data.frame(z[,1:40])))
  ypre<-as.numeric(predict(fit2,as.data.frame(xhat)))
  list5[[length(list5)+1]]<-mean((y-ypre)^2)
}
ap2<-unlist(list3)

list8<-list()
list9<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(100,rep(2,40),s)
  b<-(c(3,1.5,0,0,2,0,0,0,2,0,0,1,2,5,0,0,0,0,1,1,0,0,0,0,1,2,0,0,4,2,1,2,0,0,0,7,0,0,0,1))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(x+c+2*rnorm(100,0,6))
  fit5 <-  cv.glmnet(y=x,x=z,keep=TRUE)
  xhat = as.numeric(fit5$fit.preval[,which(fit5$lambda==fit5$lambda.min)[1]])
  fit6 = lm(y~xhat)
  pp4 = 1
  if(nrow(coef(summary(fit6)))==2) pp4 = coef(summary(fit6))[2,4]
  list8[[length(list8)+1]]<-pp4
  
  z<-mvrnorm(100,rep(2,40),s)
  b<-(c(3,1.5,0,0,2,0,0,0,2,0,0,1,2,5,0,0,0,0,1,1,0,0,0,0,1,2,0,0,4,2,1,2,0,0,0,7,0,0,0,1))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(x+c+2*rnorm(100,0,6))
  
  xhat<-as.numeric(predict(fit5 ,z, s=fit5$lambda.min))
  ppre<-as.numeric(predict(fit6,data.frame(xhat)))
  list9[[length(list9)+1]]<-mean((y-ppre)^2)
}
ap4<-unlist(list8)

fnn<-list4>0.05
fn<-length(fnn[fnn=="TRUE"])
tp<-100-fn
fnn1<-list3>0.05
fn1<-length(fnn1[fnn1=="TRUE"])
tp1<-100-fn1
fnn2<-list8>0.05
fn2<-length(fnn2[fnn2=="TRUE"])
tp2<-100-fn2


dta1<-data.frame(x=1:40)
for (j in 1:40){
  list1<-list()
  for (i in 1:40){
    tl<-0.5^abs(i-j)
    list1[[length(list1)+1]]<-tl
  }
  dta1<-cbind(dta1,unlist(list1))
}
s<-matrix(c(dta1[,2],dta1[,3],dta1[,4],dta1[,5],dta1[,6],dta1[,7],dta1[,8],dta1[,9],dta1[,10],dta1[,11],dta1[,12],dta1[,13],dta1[,14],dta1[,15],dta1[,16],dta1[,17],dta1[,18],dta1[,19],dta1[,20],dta1[,21],dta1[,22],dta1[,23],dta1[,24],dta1[,25],dta1[,26],dta1[,27],dta1[,28],dta1[,29],dta1[,30],dta1[,31],dta1[,32],dta1[,33],dta1[,34],dta1[,35],dta1[,36],dta1[,37],dta1[,38],dta1[,39],dta1[,40],dta1[,41]),ncol=40,nrow=40)


list2<-list()
list4<-list()
list6<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(100,rep(2,40),s)
  b<-(c(3,1.5,0,0,2,0,0,0,2,0,0,1,2,5,0,0,0,0,1,1,0,0,0,0,1,2,0,0,4,2,1,2,0,0,0,7,0,0,0,1))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(c+2*rnorm(100,0,6))
  crossval <-  cv.glmnet(y=x,x=z)
  penalty <- crossval$lambda.min
  fit3 <-glmnet( y=x,x=z, alpha = 1, lambda = penalty )
  xhat<-predict(crossval, newx = z, s = penalty)
  xhat<-as.vector(xhat[,1])
  fit4<-lm(y~xhat)
  pp2<-summary(fit4)
  pp3<-1
  if(nrow(coef(summary(fit4)))==2) pp3 = coef(summary(fit4))[2,4]
  tt2<-deviance(fit4)
  list2[[length(list2)+1]]<-tt2
  list4[[length(list4)+1]]<-pp3
  
  z<-mvrnorm(100,rep(2,40),s)
  b<-(c(3,1.5,0,0,2,0,0,0,2,0,0,1,2,5,0,0,0,0,1,1,0,0,0,0,1,2,0,0,4,2,1,2,0,0,0,7,0,0,0,1))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(c+2*rnorm(100,0,6))
  
  xhat<-as.numeric(predict(fit3 ,z, s=fit3$lambda.min))
  ypre<-as.numeric(predict(fit4,as.data.frame(as.numeric(xhat))))
  list6[[length(list6)+1]]<-mean((y-ypre)^2)
}
ap1<-unlist(list4)

list1<-list()
list3<-list()
list5<-list()
for (i in 1:100){
  set.seed(i)
  z<-as.data.frame(mvrnorm(100,rep(2,40),s))
  b<-(c(3,1.5,0,0,2,0,0,0,2,0,0,1,2,5,0,0,0,0,1,1,0,0,0,0,1,2,0,0,4,2,1,2,0,0,0,7,0,0,0,1))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(c+2*rnorm(100,0,6))
  z<-cbind(z,x)
  fit1<-lm(x~.,z)
  xhat<-fit1$fitted.values
  fit2<-lm(y~xhat)
  pp<-summary(fit2)
  pp1<-pp$coefficients[2,4]
  tt1<-deviance(fit2)
  mse1<-mean(pp$residuals^2)
  list1[[length(list1)+1]]<-tt1
  list3[[length(list3)+1]]<-pp1
  
  z<-as.data.frame(mvrnorm(100,rep(2,40),s))
  b<-(c(3,1.5,0,0,2,0,0,0,2,0,0,1,2,5,0,0,0,0,1,1,0,0,0,0,1,2,0,0,4,2,1,2,0,0,0,7,0,0,0,1))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(c+2*rnorm(100,0,6))
  z<-cbind(z,x)
  xhat<-as.numeric(predict(fit1,as.data.frame(z[,1:40])))
  ypre<-as.numeric(predict(fit2,as.data.frame(xhat)))
  list5[[length(list5)+1]]<-mean((y-ypre)^2)
}
ap3<-unlist(list3)

list8<-list()
list9<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(100,rep(2,40),s)
  b<-(c(3,1.5,0,0,2,0,0,0,2,0,0,1,2,5,0,0,0,0,1,1,0,0,0,0,1,2,0,0,4,2,1,2,0,0,0,7,0,0,0,1))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(c+2*rnorm(100,0,6))
  fit5 <-  cv.glmnet(y=x,x=z,keep=TRUE)
  xhat = as.numeric(fit5$fit.preval[,which(fit5$lambda==fit5$lambda.min)[1]])
  fit6 = lm(y~xhat)
  pp4 = 1
  if(nrow(coef(summary(fit6)))==2) pp4 = coef(summary(fit6))[2,4]
  list8[[length(list8)+1]]<-pp4
  
  z<-mvrnorm(100,rep(2,40),s)
  b<-(c(3,1.5,0,0,2,0,0,0,2,0,0,1,2,5,0,0,0,0,1,1,0,0,0,0,1,2,0,0,4,2,1,2,0,0,0,7,0,0,0,1))
  c<-rnorm(100,0,2)
  x<-as.vector(b%*%t(z)+3*c+10*rnorm(100,0,10))
  y<-as.vector(c+2*rnorm(100,0,6))
  
  xhat<-as.numeric(predict(fit5 ,z, s=fit5$lambda.min))
  ppre<-as.numeric(predict(fit6,data.frame(xhat)))
  list9[[length(list9)+1]]<-mean((y-ppre)^2)
}
ap5<-unlist(list8)

tnn<-list4>0.05
tn<-length(tnn[tnn=="TRUE"])
fp<-100-tn
tnn1<-list3>0.05
tn1<-length(tnn1[tnn1=="TRUE"])
fp1<-100-tn1
tnn2<-list8>0.05
tn2<-length(tnn2[tnn2=="TRUE"])
fp2<-100-tn2

mcc(tp,tn,fp,fn)
mcc(tp1,tn1,fp1,fn1)
mcc(tp2,tn2,fp2,fn2)

#ROC plot
#ROC plot
app<-c(ap,ap1)
class = rep(c(1,0),c(length(ap1),length(ap)))
ord = order(app)
plot(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),type = 'l',main = 'ROC curve',xlab = 'FPR',ylab = 'TPR',col="darkblue",lwd=4)
abline(0,1)
p = c(ap2,ap3)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="skyblue",type = 'l',lwd=4)
p = c(ap4,ap5)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="darkgreen",type = 'l',lwd=4)
p = c(ap8,ap9)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="palegreen",type = 'l',lwd=4)
legend('bottomright',c('TSL','TSLS','TSPL','LIML'),col = c("darkblue","skyblue","darkgreen","palegreen"),lty = 1,cex=1,lwd=3)


#extra confounders sim1-showing lasso better
dta1<-data.frame(x=1:8)
for (j in 1:8){
  list1<-list()
  for (i in 1:8){
    tl<-0.5^abs(i-j)
    list1[[length(list1)+1]]<-tl
  }
  dta1<-cbind(dta1,unlist(list1))
}
s<-matrix(c(dta1[,2],dta1[,3],dta1[,4],dta1[,5],dta1[,6],dta1[,7],dta1[,8],dta1[,9]),ncol=8,nrow=8)

list1<-list()
for(i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+10*c+2*rnorm(20))
  y<-as.vector(x+5*c+2*rnorm(20,0,2))
  
  p<-ivmodel(Y=y,D=x,Z=z)
  ttp<-LIML(p)
  limlp<-ttp$p.value
  list1[[length(list1)+1]]<-limlp
}
ap8<-unlist(list1)
fnn3<-list1>0.05
fn3<-length(fnn3[fnn3=="TRUE"])
tp3<-100-fn3

list2<-list()
for(i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+10*c+2*rnorm(20))
  y<-as.vector(5*c+2*rnorm(20,0,2))
  
  p<-ivmodel(Y=y,D=x,Z=z)
  ttp<-LIML(p)
  limlp<-ttp$p.value
  list2[[length(list2)+1]]<-limlp
}
ap9<-unlist(list2)

tnn3<-list2>0.05
tn3<-length(tnn3[tnn3=="TRUE"])
fp3<-100-tn3


list2<-list()
list4<-list()
list6<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-b%*%t(z)+10*c+2*rnorm(20)
  y<-as.vector(x+5*c+2*rnorm(20,0,2))
  crossval <-  cv.glmnet(y=x,x=z)
  penalty <- crossval$lambda.min
  fit3 <-glmnet( y=x,x=z, alpha = 1, lambda = penalty )
  xhat<-predict(crossval, newx = z, s = penalty)
  xhat<-as.vector(xhat[,1])
  fit4<-lm(y~xhat)
  pp2<-summary(fit4)
  pp3<-1
  if(nrow(coef(summary(fit4)))==2) pp3 = coef(summary(fit4))[2,4]
  tt2<-deviance(fit4)
  list2[[length(list2)+1]]<-tt2
  list4[[length(list4)+1]]<-pp3
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-b%*%t(z)+10*c+2*rnorm(20)
  y<-as.vector(x+5*c+2*rnorm(20,0,2))
  
  xhat<-as.numeric(predict(fit3 ,z, s=fit3$lambda.min))
  ypre<-as.numeric(predict(fit4,as.data.frame(as.numeric(xhat))))
  list6[[length(list6)+1]]<-mean((y-ypre)^2)
}

ap<-unlist(list4)

list1<-list()
list3<-list()
list5<-list()
for (i in 1:100){
  set.seed(i)
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+10*c+2*rnorm(20))
  y<-as.vector(x+5*c+2*rnorm(20,0,2))
  z<-cbind(z,x)
  fit1<-lm(x~.,z)
  xhat<-fit1$fitted.values
  fit2<-lm(y~xhat)
  pp<-summary(fit2)
  pp1<-pp$coefficients[2,4]
  tt1<-deviance(fit2)
  mse1<-mean(pp$residuals^2)
  list1[[length(list1)+1]]<-tt1
  list3[[length(list3)+1]]<-pp1
  
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-b%*%t(z)+10*c+2*rnorm(20)
  y<-as.vector(x+5*c+2*rnorm(20,0,2))
  z<-cbind(z,x)
  xhat<-as.numeric(predict(fit1,as.data.frame(z[,1:8])))
  ypre<-as.numeric(predict(fit2,as.data.frame(xhat)))
  list5[[length(list5)+1]]<-mean((y-ypre)^2)
}
ap2<-unlist(list3)

list8<-list()
list9<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20)
  x<-b%*%t(z)+10*c+2*rnorm(20)
  y<-as.vector(x+5*c+2*rnorm(20,0,2))
  fit5 <-  cv.glmnet(y=x,x=z,keep=TRUE)
  xhat = as.numeric(fit5$fit.preval[,which(fit5$lambda==fit5$lambda.min)[1]])
  fit6 = lm(y~xhat)
  pp4 = 1
  if(nrow(coef(summary(fit6)))==2) pp4 = coef(summary(fit6))[2,4]
  list8[[length(list8)+1]]<-pp4
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-b%*%t(z)+10*c+2*rnorm(20)
  y<-as.vector(x+5*c+2*rnorm(20,0,2))
  
  xhat<-as.numeric(predict(fit5 ,z, s=fit5$lambda.min))
  ppre<-as.numeric(predict(fit6,data.frame(xhat)))
  list9[[length(list9)+1]]<-mean((y-ppre)^2)
}
ap4<-unlist(list8)

fnn<-list4>0.05
fn<-length(fnn[fnn=="TRUE"])
tp<-100-fn
fnn1<-list3>0.05
fn1<-length(fnn1[fnn1=="TRUE"])
tp1<-100-fn1
fnn2<-list8>0.05
fn2<-length(fnn2[fnn2=="TRUE"])
tp2<-100-fn2

##no association
set.seed(5)
dta1<-data.frame(x=1:8)
for (j in 1:8){
  list1<-list()
  for (i in 1:8){
    tl<-0.5^abs(i-j)
    list1[[length(list1)+1]]<-tl
  }
  dta1<-cbind(dta1,unlist(list1))
}
s<-matrix(c(dta1[,2],dta1[,3],dta1[,4],dta1[,5],dta1[,6],dta1[,7],dta1[,8],dta1[,9]),ncol=8,nrow=8)


list2<-list()
list4<-list()
list6<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20)
  x<-b%*%t(z)+10*c+2*rnorm(20)
  y<-as.vector(5*c+2*rnorm(20,0,2))
  crossval <-  cv.glmnet(y=x,x=z)
  penalty <- crossval$lambda.min
  fit3 <-glmnet( y=x,x=z, alpha = 1, lambda = penalty )
  xhat<-predict(crossval, newx = z, s = penalty)
  xhat<-as.vector(xhat[,1])
  fit4<-lm(y~xhat)
  pp2<-summary(fit4)
  pp3<-1
  if(nrow(coef(summary(fit4)))==2) pp3 = coef(summary(fit4))[2,4]
  tt2<-deviance(fit4)
  list2[[length(list2)+1]]<-tt2
  list4[[length(list4)+1]]<-pp3
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-b%*%t(z)+10*c+2*rnorm(20)
  y<-as.vector(5*c+2*rnorm(20,0,2))
  
  xhat<-as.numeric(predict(fit3 ,z, s=fit3$lambda.min))
  ypre<-as.numeric(predict(fit4,as.data.frame(as.numeric(xhat))))
  list6[[length(list6)+1]]<-mean((y-ypre)^2)
}
ap1<-unlist(list4)

list1<-list()
list3<-list()
list5<-list()
for (i in 1:100){
  set.seed(i)
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-as.vector(b%*%t(z)+10*c+2*rnorm(20))
  y<-as.vector(5*c+2*rnorm(20,0,2))
  z<-cbind(z,x)
  fit1<-lm(x~.,z)
  xhat<-fit1$fitted.values
  fit2<-lm(y~xhat)
  pp<-summary(fit2)
  pp1<-pp$coefficients[2,4]
  tt1<-deviance(fit2)
  mse1<-mean(pp$residuals^2)
  list1[[length(list1)+1]]<-tt1
  list3[[length(list3)+1]]<-pp1
  
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-b%*%t(z)+10*c+2*rnorm(20)
  y<-as.vector(5*c+2*rnorm(20,0,2))
  z<-cbind(z,x)
  xhat<-as.numeric(predict(fit1,as.data.frame(z[,1:8])))
  ypre<-as.numeric(predict(fit2,as.data.frame(xhat)))
  list5[[length(list5)+1]]<-mean((y-ypre)^2)
}
ap3<-unlist(list3)

list8<-list()
list9<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20)
  x<-b%*%t(z)+10*c+2*rnorm(20)
  y<-as.vector(5*c+2*rnorm(20,0,2))
  fit5 <-  cv.glmnet(y=x,x=z,keep=TRUE)
  xhat = as.numeric(fit5$fit.preval[,which(fit5$lambda==fit5$lambda.min)[1]])
  fit6 = lm(y~xhat)
  pp4 = 1
  if(nrow(coef(summary(fit6)))==2) pp4 = coef(summary(fit6))[2,4]
  list8[[length(list8)+1]]<-pp4
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,2)
  x<-b%*%t(z)+10*c+2*rnorm(20)
  y<-as.vector(5*c+2*rnorm(20,0,2))
  
  xhat<-as.numeric(predict(fit5 ,z, s=fit5$lambda.min))
  ppre<-as.numeric(predict(fit6,data.frame(xhat)))
  list9[[length(list9)+1]]<-mean((y-ppre)^2)
}
ap5<-unlist(list8)

tnn<-list4>0.05
tn<-length(tnn[tnn=="TRUE"])
fp<-100-tn
tnn1<-list3>0.05
tn1<-length(tnn1[tnn1=="TRUE"])
fp1<-100-tn1
tnn2<-list8>0.05
tn2<-length(tnn2[tnn2=="TRUE"])
fp2<-100-tn2

mcc(tp,tn,fp,fn)
mcc(tp1,tn1,fp1,fn1)
mcc(tp2,tn2,fp2,fn2)

app<-c(ap,ap1)
class = rep(c(1,0),c(length(ap1),length(ap)))
ord = order(app)
plot(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),type = 'l',main = 'ROC curve',xlab = 'FPR',ylab = 'TPR',col="yellow")
abline(0,1)
p = c(ap2,ap3)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="blue",type = 'l')
p = c(ap4,ap5)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="red",type = 'l')
app5<-c(ap8,ap9)
ord = order(app5)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),type = 'l',main = 'ROC curve',xlab = 'FPR',ylab = 'TPR',col="green")
abline(0,1)
legend('bottomright',c('TSL','TSLS','TSPL','LIML'),col = c("yellow","blue","red","green"),lty = 1,cex=0.5)

#ROC plot
app<-c(ap,ap1)
class = rep(c(1,0),c(length(ap1),length(ap)))
ord = order(app)
plot(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),type = 'l',main = 'ROC curve',xlab = 'FPR',ylab = 'TPR',col="darkblue",lwd=4)
abline(0,1)
p = c(ap2,ap3)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="skyblue",type = 'l',lwd=4)
p = c(ap4,ap5)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="darkgreen",type = 'l',lwd=4)
p = c(ap8,ap9)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="palegreen",type = 'l',lwd=4)
legend('bottomright',c('TSL','TSLS','TSPL','LIML'),col = c("darkblue","skyblue","darkgreen","palegreen"),lty = 1,cex=1,lwd=3)




#neg association x and y - preval better
dta1<-data.frame(x=1:8)
for (j in 1:8){
  list1<-list()
  for (i in 1:8){
    tl<-0.5^abs(i-j)
    list1[[length(list1)+1]]<-tl
  }
  dta1<-cbind(dta1,unlist(list1))
}
s<-matrix(c(dta1[,2],dta1[,3],dta1[,4],dta1[,5],dta1[,6],dta1[,7],dta1[,8],dta1[,9]),ncol=8,nrow=8)

list1<-list()
for(i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(1,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(-0.25*x+c+3*rnorm(20))
  
  p<-ivmodel(Y=y,D=x,Z=z)
  ttp<-LIML(p)
  limlp<-ttp$p.value
  list1[[length(list1)+1]]<-limlp
}
ap8<-unlist(list1)
fnn3<-list1>0.05
fn3<-length(fnn3[fnn3=="TRUE"])
tp3<-100-fn3

list2<-list()
for(i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(1,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(c+3*rnorm(20))
  
  p<-ivmodel(Y=y,D=x,Z=z)
  ttp<-LIML(p)
  limlp<-ttp$p.value
  list2[[length(list2)+1]]<-limlp
}
ap9<-unlist(list2)

tnn3<-list2>0.05
tn3<-length(tnn3[tnn3=="TRUE"])
fp3<-100-tn3


list2<-list()
list4<-list()
list6<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(1,1)
  x<-g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20)
  y<-as.vector(-0.25*x+c+3*rnorm(20))
  crossval <-  cv.glmnet(y=x,x=z)
  penalty <- crossval$lambda.min
  fit3 <-glmnet( y=x,x=z, alpha = 1, lambda = penalty )
  xhat<-predict(crossval, newx = z, s = penalty)
  xhat<-as.vector(xhat[,1])
  fit4<-lm(y~xhat)
  pp2<-summary(fit4)
  pp3<-1
  if(nrow(coef(summary(fit4)))==2) pp3 = coef(summary(fit4))[2,4]
  tt2<-deviance(fit4)
  list2[[length(list2)+1]]<-tt2
  list4[[length(list4)+1]]<-pp3
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(1,1)
  x<-g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20)
  y<-as.vector(-0.25*x+c+3*rnorm(20))
  
  xhat<-as.numeric(predict(fit3 ,z, s=fit3$lambda.min))
  ypre<-as.numeric(predict(fit4,as.data.frame(as.numeric(xhat))))
  list6[[length(list6)+1]]<-mean((y-ypre)^2)
}
ap<-unlist(list4)

list1<-list()
list3<-list()
list5<-list()
for (i in 1:100){
  set.seed(i)
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(1,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(-0.25*x+c+3*rnorm(20))
  z<-cbind(z,x)
  fit1<-lm(x~.,z)
  xhat<-fit1$fitted.values
  fit2<-lm(y~xhat)
  pp<-summary(fit2)
  pp1<-pp$coefficients[2,4]
  tt1<-deviance(fit2)
  mse1<-mean(pp$residuals^2)
  list1[[length(list1)+1]]<-tt1
  list3[[length(list3)+1]]<-pp1
  
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(1,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(-0.25*x+c+3*rnorm(20))
  z<-cbind(z,x)
  xhat<-as.numeric(predict(fit1,as.data.frame(z[,1:8])))
  ypre<-as.numeric(predict(fit2,as.data.frame(xhat)))
  list5[[length(list5)+1]]<-mean((y-ypre)^2)
}
ap2<-unlist(list3)

list8<-list()
list9<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(1,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(-0.25*x+c+3*rnorm(20))
  fit5 <-  cv.glmnet(y=x,x=z,keep=TRUE)
  xhat = as.numeric(fit5$fit.preval[,which(fit5$lambda==fit5$lambda.min)[1]])
  fit6 = lm(y~xhat)
  pp4 = 1
  if(nrow(coef(summary(fit6)))==2) pp4 = coef(summary(fit6))[2,4]
  list8[[length(list8)+1]]<-pp4
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(1,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(-0.25*x+c+3*rnorm(20))
  
  xhat<-as.numeric(predict(fit5 ,z, s=fit5$lambda.min))
  ppre<-as.numeric(predict(fit6,data.frame(xhat)))
  list9[[length(list9)+1]]<-mean((y-ppre)^2)
}

ap4<-unlist(list8)

fnn<-list4>0.05
fn<-length(fnn[fnn=="TRUE"])
tp<-100-fn
fnn1<-list3>0.05
fn1<-length(fnn1[fnn1=="TRUE"])
tp1<-100-fn1
fnn2<-list8>0.05
fn2<-length(fnn2[fnn2=="TRUE"])
tp2<-100-fn2

##no association
dta1<-data.frame(x=1:8)
for (j in 1:8){
  list1<-list()
  for (i in 1:8){
    tl<-0.5^abs(i-j)
    list1[[length(list1)+1]]<-tl
  }
  dta1<-cbind(dta1,unlist(list1))
}
s<-matrix(c(dta1[,2],dta1[,3],dta1[,4],dta1[,5],dta1[,6],dta1[,7],dta1[,8],dta1[,9]),ncol=8,nrow=8)


list2<-list()
list4<-list()
list6<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(1,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(c+3*rnorm(20))
  crossval <-  cv.glmnet(y=x,x=z)
  penalty <- crossval$lambda.min
  fit3 <-glmnet( y=x,x=z, alpha = 1, lambda = penalty )
  xhat<-predict(crossval, newx = z, s = penalty)
  xhat<-as.vector(xhat[,1])
  fit4<-lm(y~xhat)
  pp2<-summary(fit4)
  pp3<-1
  if(nrow(coef(summary(fit4)))==2) pp3 = coef(summary(fit4))[2,4]
  tt2<-deviance(fit4)
  list2[[length(list2)+1]]<-tt2
  list4[[length(list4)+1]]<-pp3
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(1,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(c+3*rnorm(20))
  
  xhat<-as.numeric(predict(fit3 ,z, s=fit3$lambda.min))
  ypre<-as.numeric(predict(fit4,as.data.frame(as.numeric(xhat))))
  list6[[length(list6)+1]]<-mean((y-ypre)^2)
}
ap1<-unlist(list4)

list1<-list()
list3<-list()
list5<-list()
for (i in 1:100){
  set.seed(i)
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(1,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(c+3*rnorm(20))
  z<-cbind(z,x)
  fit1<-lm(x~.,z)
  xhat<-fit1$fitted.values
  fit2<-lm(y~xhat)
  pp<-summary(fit2)
  pp1<-pp$coefficients[2,4]
  tt1<-deviance(fit2)
  mse1<-mean(pp$residuals^2)
  list1[[length(list1)+1]]<-tt1
  list3[[length(list3)+1]]<-pp1
  
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(1,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(c+3*rnorm(20))
  z<-cbind(z,x)
  xhat<-as.numeric(predict(fit1,as.data.frame(z[,1:8])))
  ypre<-as.numeric(predict(fit2,as.data.frame(xhat)))
  list5[[length(list5)+1]]<-mean((y-ypre)^2)
}

ap3<-unlist(list3)


list8<-list()
list9<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(1,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(c+3*rnorm(20))
  fit5 <-  cv.glmnet(y=x,x=z,keep=TRUE)
  xhat = as.numeric(fit5$fit.preval[,which(fit5$lambda==fit5$lambda.min)[1]])
  fit6 = lm(y~xhat)
  pp4 = 1
  if(nrow(coef(summary(fit6)))==2) pp4 = coef(summary(fit6))[2,4]
  list8[[length(list8)+1]]<-pp4
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(1,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(c+3*rnorm(20))
  
  xhat<-as.numeric(predict(fit5 ,z, s=fit5$lambda.min))
  ppre<-as.numeric(predict(fit6,data.frame(xhat)))
  list9[[length(list9)+1]]<-mean((y-ppre)^2)
}
ap5<-unlist(list8)


mcc(tp,tn,fp,fn)
mcc(tp1,tn1,fp1,fn1)
mcc(tp2,tn2,fp2,fn2)


#ROC plot
app<-c(ap,ap1)
class = rep(c(1,0),c(length(ap1),length(ap)))
ord = order(app)
plot(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),type = 'l',main = 'ROC curve',xlab = 'FPR',ylab = 'TPR',col="darkblue",lwd=4)
abline(0,1)
p = c(ap2,ap3)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="skyblue",type = 'l',lwd=4)
p = c(ap4,ap5)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="darkgreen",type = 'l',lwd=4)
p = c(ap8,ap9)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="palegreen",type = 'l',lwd=4)
legend('bottomright',c('TSL','TSLS','TSPL','LIML'),col = c("darkblue","skyblue","darkgreen","palegreen"),lty = 1,cex=1,lwd=3)


##two confounders

dta1<-data.frame(x=1:8)
for (j in 1:8){
  list1<-list()
  for (i in 1:8){
    tl<-0.5^abs(i-j)
    list1[[length(list1)+1]]<-tl
  }
  dta1<-cbind(dta1,unlist(list1))
}
s<-matrix(c(dta1[,2],dta1[,3],dta1[,4],dta1[,5],dta1[,6],dta1[,7],dta1[,8],dta1[,9]),ncol=8,nrow=8)

list1<-list()
for(i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  c1<-rnorm(20,0,3)
  g=c(1,1,2)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c,c1))+3*rnorm(20))
  y<-as.vector(x+-2*c+-5*c1+3*rnorm(20))
  
  p<-ivmodel(Y=y,D=x,Z=z)
  ttp<-LIML(p)
  limlp<-ttp$p.value
  list1[[length(list1)+1]]<-limlp
}
ap8<-unlist(list1)
fnn3<-list1>0.05
fn3<-length(fnn3[fnn3=="TRUE"])
tp3<-100-fn3

list2<-list()
for(i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  c1<-rnorm(20,0,3)
  g=c(1,1,2)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c,c1))+3*rnorm(20))
  y<-as.vector(x+-2*c+-5*c1+3*rnorm(20))
  
  p<-ivmodel(Y=y,D=x,Z=z)
  ttp<-LIML(p)
  limlp<-ttp$p.value
  list2[[length(list2)+1]]<-limlp
}
ap9<-unlist(list2)

tnn3<-list2>0.05
tn3<-length(tnn3[tnn3=="TRUE"])
fp3<-100-tn3


list2<-list()
list4<-list()
list6<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  c1<-rnorm(20,0,3)
  g=c(1,1,2)
  x<-g%*%t(cbind(as.numeric(b%*%t(z)),c,c1))+3*rnorm(20)
  y<-as.vector(x+-2*c+-5*c1+3*rnorm(20))
  crossval <-  cv.glmnet(y=x,x=z)
  penalty <- crossval$lambda.min
  fit3 <-glmnet( y=x,x=z, alpha = 1, lambda = penalty )
  xhat<-predict(crossval, newx = z, s = penalty)
  xhat<-as.vector(xhat[,1])
  fit4<-lm(y~xhat)
  pp2<-summary(fit4)
  pp3<-1
  if(nrow(coef(summary(fit4)))==2) pp3 = coef(summary(fit4))[2,4]
  tt2<-deviance(fit4)
  list2[[length(list2)+1]]<-tt2
  list4[[length(list4)+1]]<-pp3
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  c1<-rnorm(20,0,3)
  g=c(1,1,2)
  x<-g%*%t(cbind(as.numeric(b%*%t(z)),c,c1))+3*rnorm(20)
  y<-as.vector(x+-2*c+-5*c1+3*rnorm(20))
  
  xhat<-as.numeric(predict(fit3 ,z, s=fit3$lambda.min))
  ypre<-as.numeric(predict(fit4,as.data.frame(as.numeric(xhat))))
  list6[[length(list6)+1]]<-mean((y-ypre)^2)
}
ap<-unlist(list4)

list1<-list()
list3<-list()
list5<-list()
for (i in 1:100){
  set.seed(i)
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  c1<-rnorm(20,0,3)
  g=c(1,1,2)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c,c1))+3*rnorm(20))
  y<-as.vector(x+-2*c+-5*c1+3*rnorm(20))
  z<-cbind(z,x)
  fit1<-lm(x~.,z)
  xhat<-fit1$fitted.values
  fit2<-lm(y~xhat)
  pp<-summary(fit2)
  pp1<-pp$coefficients[2,4]
  tt1<-deviance(fit2)
  mse1<-mean(pp$residuals^2)
  list1[[length(list1)+1]]<-tt1
  list3[[length(list3)+1]]<-pp1
  
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  c1<-rnorm(20,0,3)
  g=c(1,1,2)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c,c1))+3*rnorm(20))
  y<-as.vector(x+-2*c+-5*c1+3*rnorm(20))
  z<-cbind(z,x)
  xhat<-as.numeric(predict(fit1,as.data.frame(z[,1:8])))
  ypre<-as.numeric(predict(fit2,as.data.frame(xhat)))
  list5[[length(list5)+1]]<-mean((y-ypre)^2)
}
ap2<-unlist(list3)

list8<-list()
list9<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  c1<-rnorm(20,0,3)
  g=c(1,1,2)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c,c1))+3*rnorm(20))
  y<-as.vector(x+-2*c+-5*c1+3*rnorm(20))
  fit5 <-  cv.glmnet(y=x,x=z,keep=TRUE)
  xhat = as.numeric(fit5$fit.preval[,which(fit5$lambda==fit5$lambda.min)[1]])
  fit6 = lm(y~xhat)
  pp4 = 1
  if(nrow(coef(summary(fit6)))==2) pp4 = coef(summary(fit6))[2,4]
  list8[[length(list8)+1]]<-pp4
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  c1<-rnorm(20,0,3)
  g=c(1,1,2)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c,c1))+3*rnorm(20))
  y<-as.vector(x+-2*c+-5*c1+3*rnorm(20))
  
  xhat<-as.numeric(predict(fit5 ,z, s=fit5$lambda.min))
  ppre<-as.numeric(predict(fit6,data.frame(xhat)))
  list9[[length(list9)+1]]<-mean((y-ppre)^2)
}

ap4<-unlist(list8)

fnn<-list4>0.05
fn<-length(fnn[fnn=="TRUE"])
tp<-100-fn
fnn1<-list3>0.05
fn1<-length(fnn1[fnn1=="TRUE"])
tp1<-100-fn1
fnn2<-list8>0.05
fn2<-length(fnn2[fnn2=="TRUE"])
tp2<-100-fn2

##no association
set.seed(5)
dta1<-data.frame(x=1:8)
for (j in 1:8){
  list1<-list()
  for (i in 1:8){
    tl<-0.5^abs(i-j)
    list1[[length(list1)+1]]<-tl
  }
  dta1<-cbind(dta1,unlist(list1))
}
s<-matrix(c(dta1[,2],dta1[,3],dta1[,4],dta1[,5],dta1[,6],dta1[,7],dta1[,8],dta1[,9]),ncol=8,nrow=8)


list2<-list()
list4<-list()
list6<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  c1<-rnorm(20,0,3)
  g=c(1,1,2)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c,c1))+3*rnorm(20))
  y<-as.vector(-2*c+-5*c1+3*rnorm(20))
  crossval <-  cv.glmnet(y=x,x=z)
  penalty <- crossval$lambda.min
  fit3 <-glmnet( y=x,x=z, alpha = 1, lambda = penalty )
  xhat<-predict(crossval, newx = z, s = penalty)
  xhat<-as.vector(xhat[,1])
  fit4<-lm(y~xhat)
  pp2<-summary(fit4)
  pp3<-1
  if(nrow(coef(summary(fit4)))==2) pp3 = coef(summary(fit4))[2,4]
  tt2<-deviance(fit4)
  list2[[length(list2)+1]]<-tt2
  list4[[length(list4)+1]]<-pp3
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  c1<-rnorm(20,0,3)
  g=c(1,1,2)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c,c1))+3*rnorm(20))
  y<-as.vector(-2*c+-5*c1+3*rnorm(20))
  
  xhat<-as.numeric(predict(fit3 ,z, s=fit3$lambda.min))
  ypre<-as.numeric(predict(fit4,as.data.frame(as.numeric(xhat))))
  list6[[length(list6)+1]]<-mean((y-ypre)^2)
}
ap1<-unlist(list4)

list1<-list()
list3<-list()
list5<-list()
for (i in 1:100){
  set.seed(i)
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  c1<-rnorm(20,0,3)
  g=c(1,1,2)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c,c1))+3*rnorm(20))
  y<-as.vector(-2*c+-5*c1+3*rnorm(20))
  z<-cbind(z,x)
  fit1<-lm(x~.,z)
  xhat<-fit1$fitted.values
  fit2<-lm(y~xhat)
  pp<-summary(fit2)
  pp1<-pp$coefficients[2,4]
  tt1<-deviance(fit2)
  mse1<-mean(pp$residuals^2)
  list1[[length(list1)+1]]<-tt1
  list3[[length(list3)+1]]<-pp1
  
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  c1<-rnorm(20,0,3)
  g=c(1,1,2)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c,c1))+3*rnorm(20))
  y<-as.vector(-2*c+-5*c1+3*rnorm(20))
  z<-cbind(z,x)
  xhat<-as.numeric(predict(fit1,as.data.frame(z[,1:8])))
  ypre<-as.numeric(predict(fit2,as.data.frame(xhat)))
  list5[[length(list5)+1]]<-mean((y-ypre)^2)
}

ap3<-unlist(list3)


list8<-list()
list9<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  c1<-rnorm(20,0,3)
  g=c(1,1,2)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c,c1))+3*rnorm(20))
  y<-as.vector(-2*c+-5*c1+3*rnorm(20))
  fit5 <-  cv.glmnet(y=x,x=z,keep=TRUE)
  xhat = as.numeric(fit5$fit.preval[,which(fit5$lambda==fit5$lambda.min)[1]])
  fit6 = lm(y~xhat)
  pp4 = 1
  if(nrow(coef(summary(fit6)))==2) pp4 = coef(summary(fit6))[2,4]
  list8[[length(list8)+1]]<-pp4
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  c1<-rnorm(20,0,3)
  g=c(1,1,2)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c,c1))+3*rnorm(20))
  y<-as.vector(-2*c+-5*c1+3*rnorm(20))
  
  xhat<-as.numeric(predict(fit5 ,z, s=fit5$lambda.min))
  ppre<-as.numeric(predict(fit6,data.frame(xhat)))
  list9[[length(list9)+1]]<-mean((y-ppre)^2)
}
ap5<-unlist(list8)

tnn<-list4>0.05
tn<-length(tnn[tnn=="TRUE"])
fp<-100-tn
tnn1<-list3>0.05
tn1<-length(tnn1[tnn1=="TRUE"])
fp1<-100-tn1
tnn2<-list8>0.05
tn2<-length(tnn2[tnn2=="TRUE"])
fp2<-100-tn2

mcc(tp,tn,fp,fn)
mcc(tp1,tn1,fp1,fn1)
mcc(tp2,tn2,fp2,fn2)

#ROC plot
app<-c(ap,ap1)
class = rep(c(1,0),c(length(ap1),length(ap)))
ord = order(app)
plot(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),type = 'l',main = 'ROC curve',xlab = 'FPR',ylab = 'TPR',col="darkblue",lwd=4)
abline(0,1)
p = c(ap2,ap3)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="skyblue",type = 'l',lwd=4)
p = c(ap4,ap5)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="darkgreen",type = 'l',lwd=4)
p = c(ap8,ap9)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="palegreen",type = 'l',lwd=4)
legend('bottomright',c('TSL','TSLS','TSPL','LIML'),col = c("darkblue","skyblue","darkgreen","palegreen"),lty = 1,cex=1,lwd=3)



##no informative IVs
dta1<-data.frame(x=1:8)
for (j in 1:8){
  list1<-list()
  for (i in 1:8){
    tl<-0.5^abs(i-j)
    list1[[length(list1)+1]]<-tl
  }
  dta1<-cbind(dta1,unlist(list1))
}
s<-matrix(c(dta1[,2],dta1[,3],dta1[,4],dta1[,5],dta1[,6],dta1[,7],dta1[,8],dta1[,9]),ncol=8,nrow=8)

list1<-list()
for(i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(.5,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(x+c+3*rnorm(20))
  p<-ivmodel(Y=y,D=x,Z=z)
  ttp<-LIML(p)
  limlp<-ttp$p.value
  list1[[length(list1)+1]]<-limlp
}
ap8<-unlist(list1)
fnn3<-list1>0.05
fn3<-length(fnn3[fnn3=="TRUE"])
tp3<-100-fn3

list2<-list()
for(i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(0,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(x+c+3*rnorm(20))
  
  p<-ivmodel(Y=y,D=x,Z=z)
  ttp<-LIML(p)
  limlp<-ttp$p.value
  list2[[length(list2)+1]]<-limlp
}
ap9<-unlist(list2)

tnn3<-list2>0.05
tn3<-length(tnn3[tnn3=="TRUE"])
fp3<-100-tn3


list2<-list()
list4<-list()
list6<-list()
list12<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(.5,1)
  x<-g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20)
  y<-as.vector(x+c+3*rnorm(20))
  crossval <-  cv.glmnet(y=x,x=z)
  penalty <- crossval$lambda.min
  fit3 <-glmnet( y=x,x=z, alpha = 1, lambda = penalty )
  xhat<-predict(crossval, newx = z, s = penalty)
  xhat<-as.vector(xhat[,1])
  fit4<-lm(y~xhat)
  rl<-fit4$coefficients[2]
  list12[[length(list12)+1]]<-rl
  pp2<-summary(fit4)
  pp3<-1
  if(nrow(coef(summary(fit4)))==2) pp3 = coef(summary(fit4))[2,4]
  tt2<-deviance(fit4)
  list2[[length(list2)+1]]<-tt2
  list4[[length(list4)+1]]<-pp3
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(.5,1)
  x<-g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20)
  y<-as.vector(x+c+3*rnorm(20))
  
  xhat<-as.numeric(predict(fit3 ,z, s=fit3$lambda.min))
  ypre<-as.numeric(predict(fit4,as.data.frame(as.numeric(xhat))))
  list6[[length(list6)+1]]<-mean((y-ypre)^2)
}
ap<-unlist(list4)

list1<-list()
list3<-list()
list5<-list()
list11<-list()
for (i in 1:100){
  set.seed(i)
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(.5,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(x+c+3*rnorm(20))
  z<-cbind(z,x)
  fit1<-lm(x~.,z)
  xhat<-fit1$fitted.values
  fit2<-lm(y~xhat)
  ro<-fit2$coefficients[2]
  list11[[length(list11)+1]]<-ro
  pp<-summary(fit2)
  pp1<-pp$coefficients[2,4]
  tt1<-deviance(fit2)
  mse1<-mean(pp$residuals^2)
  list1[[length(list1)+1]]<-tt1
  list3[[length(list3)+1]]<-pp1
  
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(.5,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(x+c+3*rnorm(20))
  z<-cbind(z,x)
  xhat<-as.numeric(predict(fit1,as.data.frame(z[,1:8])))
  ypre<-as.numeric(predict(fit2,as.data.frame(xhat)))
  list5[[length(list5)+1]]<-mean((y-ypre)^2)
}
ap2<-unlist(list3)

list8<-list()
list9<-list()
list10<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(.5,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(x+c+3*rnorm(20))
  fit5 <-  cv.glmnet(y=x,x=z,keep=TRUE)
  xhat = as.numeric(fit5$fit.preval[,which(fit5$lambda==fit5$lambda.min)[1]])
  fit6 = lm(y~xhat)
  rp<-fit6$coefficients[2]
  list10[[length(list10)+1]]<-rp
  pp4 = 1
  if(nrow(coef(summary(fit6)))==2) pp4 = coef(summary(fit6))[2,4]
  list8[[length(list8)+1]]<-pp4
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(.5,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(x+c+3*rnorm(20))
  
  xhat<-as.numeric(predict(fit5 ,z, s=fit5$lambda.min))
  ppre<-as.numeric(predict(fit6,data.frame(xhat)))
  list9[[length(list9)+1]]<-mean((y-ppre)^2)
}

ap4<-unlist(list8)
cd1<-unlist(list12)
cd2<-unlist(list11)
cd3<-unlist(list10)

fnn<-list4>0.05
fn<-length(fnn[fnn=="TRUE"])
tp<-100-fn
fnn1<-list3>0.05
fn1<-length(fnn1[fnn1=="TRUE"])
tp1<-100-fn1
fnn2<-list8>0.05
fn2<-length(fnn2[fnn2=="TRUE"])
tp2<-100-fn2

##no association
set.seed(5)
dta1<-data.frame(x=1:8)
for (j in 1:8){
  list1<-list()
  for (i in 1:8){
    tl<-0.5^abs(i-j)
    list1[[length(list1)+1]]<-tl
  }
  dta1<-cbind(dta1,unlist(list1))
}
s<-matrix(c(dta1[,2],dta1[,3],dta1[,4],dta1[,5],dta1[,6],dta1[,7],dta1[,8],dta1[,9]),ncol=8,nrow=8)


list2<-list()
list4<-list()
list6<-list()
list12<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(0,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(x+c+3*rnorm(20))
  crossval <-  cv.glmnet(y=x,x=z)
  penalty <- crossval$lambda.min
  fit3 <-glmnet( y=x,x=z, alpha = 1, lambda = penalty )
  xhat<-predict(crossval, newx = z, s = penalty)
  xhat<-as.vector(xhat[,1])
  fit4<-lm(y~xhat)
  rl<-fit4$coefficients[2]
  list12[[length(list12)+1]]<-rl
  pp2<-summary(fit4)
  pp3<-1
  if(nrow(coef(summary(fit4)))==2) pp3 = coef(summary(fit4))[2,4]
  tt2<-deviance(fit4)
  list2[[length(list2)+1]]<-tt2
  list4[[length(list4)+1]]<-pp3
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(0,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(x+c+3*rnorm(20))
  
  xhat<-as.numeric(predict(fit3 ,z, s=fit3$lambda.min))
  ypre<-as.numeric(predict(fit4,as.data.frame(as.numeric(xhat))))
  list6[[length(list6)+1]]<-mean((y-ypre)^2)
}
ap1<-unlist(list4)

list1<-list()
list3<-list()
list5<-list()
list11<-list()
for (i in 1:100){
  set.seed(i)
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(0,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(x+c+3*rnorm(20))
  z<-cbind(z,x)
  fit1<-lm(x~.,z)
  xhat<-fit1$fitted.values
  fit2<-lm(y~xhat)
  ro<-fit2$coefficients[2]
  list11[[length(list11)+1]]<-ro
  pp<-summary(fit2)
  pp1<-pp$coefficients[2,4]
  tt1<-deviance(fit2)
  mse1<-mean(pp$residuals^2)
  list1[[length(list1)+1]]<-tt1
  list3[[length(list3)+1]]<-pp1
  
  z<-as.data.frame(mvrnorm(20,rep(2,8),s))
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(0,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(x+c+3*rnorm(20))
  z<-cbind(z,x)
  xhat<-as.numeric(predict(fit1,as.data.frame(z[,1:8])))
  ypre<-as.numeric(predict(fit2,as.data.frame(xhat)))
  list5[[length(list5)+1]]<-mean((y-ypre)^2)
}

ap3<-unlist(list3)


list8<-list()
list9<-list()
list10<-list()
for (i in 1:100){
  set.seed(i)
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(0,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(x+c+3*rnorm(20))
  fit5 <-  cv.glmnet(y=x,x=z,keep=TRUE)
  xhat = as.numeric(fit5$fit.preval[,which(fit5$lambda==fit5$lambda.min)[1]])
  fit6 = lm(y~xhat)
  rp<-fit6$coefficients[2]
  list10[[length(list10)+1]]<-rp
  pp4 = 1
  if(nrow(coef(summary(fit6)))==2) pp4 = coef(summary(fit6))[2,4]
  list8[[length(list8)+1]]<-pp4
  
  z<-mvrnorm(20,rep(2,8),s)
  b<-(c(3,1.5,0,0,2,0,0,0))
  c<-rnorm(20,0,3)
  g=c(0,1)
  x<-as.vector(g%*%t(cbind(as.numeric(b%*%t(z)),c))+3*rnorm(20))
  y<-as.vector(x+c+3*rnorm(20))
  
  xhat<-as.numeric(predict(fit5 ,z, s=fit5$lambda.min))
  ppre<-as.numeric(predict(fit6,data.frame(xhat)))
  list9[[length(list9)+1]]<-mean((y-ppre)^2)
}
ap5<-unlist(list8)

cd4<-unlist(list12)
cd5<-unlist(list11)
cd6<-unlist(list10)

tnn<-list4>0.05
tn<-length(tnn[tnn=="TRUE"])
fp<-100-tn
tnn1<-list3>0.05
tn1<-length(tnn1[tnn1=="TRUE"])
fp1<-100-tn1
tnn2<-list8>0.05
tn2<-length(tnn2[tnn2=="TRUE"])
fp2<-100-tn2



#make MCC function
mcc<-function(a,b,c,d){
  (a*b-c*d)/sqrt((a+c)*(b+d)*(c+b)*(a+d))
}
#a=tp b=tn c=fp d=fn

mcc(tp,tn,fp,fn)
mcc(tp1,tn1,fp1,fn1)
mcc(tp2,tn2,fp2,fn2)

#ROC plot
app<-c(ap,ap1)
class = rep(c(1,0),c(length(ap1),length(ap)))
ord = order(app)
plot(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),type = 'l',main = 'ROC curve',xlab = 'FPR',ylab = 'TPR',col="darkblue",lwd=4)
abline(0,1)
p = c(ap2,ap3)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="skyblue",type = 'l',lwd=4)
p = c(ap4,ap5)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="darkgreen",type = 'l',lwd=4)
p = c(ap8,ap9)
ord = order(p)
points(cumsum(1-class[ord])/sum(1-class),cumsum(class[ord])/sum(class),col="palegreen",type = 'l',lwd=4)
legend('bottomright',c('TSL','TSLS','TSPL','LIML'),col = c("darkblue","skyblue","darkgreen","palegreen"),lty = 1,cex=1,lwd=3)
