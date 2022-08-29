###Workshop 5 Classification Analysis


#set working directory
setwd("~/Uni/4th year/Bioinfo/Machine Learning")

#install and load packages
#install.packages(c('corrplot', 'brnn', 'nnet'))
library(corrplot)
library(brnn)
library(nnet)
library(VennDiagram)
library(MASS)

train=read.csv('seeds.train.csv')
test=read.csv('seeds.test.csv')
train.inp=train[,which(colnames(train)!='class')]
test.inp=test[,which(colnames(train)!='class')]

train$class=train$class-1
test$class=test$class-1

########data visualisation using PCA
pca=prcomp(train.inp)
eigen2=round(100*sum(pca$sdev[1:2]^2)/sum(pca$sdev^2),digits=2)
par(mfrow=c(1,2))
barplot(pca$sdev^2,beside=TRUE,las=2)
legend('topright',paste('Var(2)=',eigen2),cex=.4)
plot(pca$x[,1:2],pch=c(1,3)[train$class],xlab='PC1',ylab='PC2')

#######confusion matrix time!!! write a function and run it

confusion=function(y,y.hat)
{
  fd=table(y,y.hat)
  fd=cbind(fd,c(fd[1,1]/sum(fd[1,]),fd[2,2]/sum(fd[2,])))
  fd=rbind(fd,c(fd[1,1]/sum(fd[,1]),fd[2,2]/sum(fd[,2]),sum(diag(fd))/sum(fd)))
  return(fd)
}

######apply LDA 
model=lda(class~.,train)
y.lda=predict(model,newdata=test.inp)$class
print(confusion(test$class,y.lda))

######apply brnn
model=brnn(class~.,train,neurons=7)
y.brnn=predict(model,newdata=test.inp)
print(confusion(test$class,y.brnn>0.5))

#####apply nnet 
model=nnet(class~.,train,size=7)
y.nnet=predict(model,newdata=test.inp)
print(confusion(test$class,y.nnet>0.5))


######consensus analysis
sen.ven=list(LDA=which(test$class==1 & y.lda==1),
             BRNN=which(test$class==1 & y.brnn>0.5),
             NNET=which(test$class==1 & y.nnet>0.5))
venn.diagram(sen.ven,file='sensitivity.venn.tiff',fill=c('red','blue','green'),
             col=c('red','blue','green'),cex=2,cex.cat=2,mar=0.2)
spe.ven=list(LDA=which(test$class==0 & y.lda==0),
             BRNN=which(test$class==0 & y.brnn<=0.5),
             NNET=which(test$class==0 & y.nnet<=0.5))
venn.diagram(spe.ven,file='specificity.venn2.tiff',fill=c('red','blue','green'),
             col=c('red','blue','green'),cex=2,cex.cat=2,mar=0.2, ext.text=T)









