### Workshop 3!!! Regression analysis!!!
##also dont forget that regression assumes a lot of noise 
## that might be useful for your paper 
## this was mentioned in the reading, beginning of linear regression

#set working directory
setwd("~/Uni/4th year/Bioinfo/Machine Learning")

#install and load packages
#install.packages(c('corrplot', 'brnn', 'nnet'))
library(corrplot)
library(brnn)
library(nnet)

#load and normalise data
X = read.csv('phi.coefficient.csv')
minimum = min(X$Phi)
maximum = max(X$Phi)
X$Phi = (X$Phi-minimum)/(maximum-minimum)

#preliminary correlation analysis
par(bg='grey90')
corrplot(cor(X))

#very cool^
#anyway, now lets do Linear Regression! 
model = lm(Phi~.,X)
print(model)

stats = summary(model)
coef.mat = data.frame(coef(stats))

r.sq.lm = stats$r.squared
F.lm = stats$fstatistic
p.lm = pf(F.lm[1], F.lm[2], F.lm[3], lower.tail = F)

#nice, got that sorted out, now lets get our CI bands and plot it 
bands = data.frame(predict(model,int='c'))
or=order(bands$fit)
X=X[or,]
bands=bands[or,]
scale=range(bands)
plot(X$Phi, ylim=scale, ylab='Phi', col='gray',pch=19)
lines(bands$fit)
lines(bands$upr,lty=2)
lines(bands$lwr,lty=2)
legend('topleft', pch=c(19,NA,NA),lty = c(NA,1,2), cex=.4,
       legend = c('Observations', 'Regression Model', 'Confidence Bands'))
legend('bottomright', legend = c(paste('p=', round(p.lm,digits=10)),
        paste('F statistic', round(F.lm[1],digits=3)),
        paste('Rsquare=', round(r.sq.lm, digits=3))), cex=.4)


########time to do some
##########NonLinear Neural Network Regression Analysis Using nnet 
model=nnet(Phi~.,data=X,size=6)
y.nnet=predict(model)
r.sq.nnet=1-var(y.nnet-X$Phi)/var(X$Phi)
or=order(y.nnet)
plot(X$Phi[or],ylab = 'Phi', col= 'grey',pch=19)
lines(y.nnet[or],lwd=2)
legend('topleft',paste('R squared=',round(r.sq.nnet,digits=6)), cex=.4)

par(cex=0.6)
model=nnet(Phi~.,data=X,size=6)
y.nnet=predict(model)
r.sq.nnet=1-var(y.nnet-X$Phi)/var(X$Phi)
or=order(y.nnet)
plot(X$Phi[or],ylab='Phi',col='grey',pch=19)
lines(y.nnet[or],lwd=2)
legend('topleft',paste('R.square=',round(r.sq.nnet,digits=6)))
