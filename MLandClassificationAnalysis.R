#### Machine learning and classification analysis

#set working directory
setwd("~/Uni/4th year/Bioinfo/Machine Learning")

#install and load packages
#install.packages(c('MASS', 'brnn', 'nnet', 'VennDiagram'))
library(MASS)
library(brnn)
library(nnet)
library(VennDiagram)

#load and name data variables
train=read.csv('seeds.train.csv', header = T)
test=read.csv('seeds.test.csv', header = T)
train.inp=train[,which(colnames(train)!='class')]
test.inp=test[,which(colnames(train)!='class')]

train$class <- train$class-1
test$class <- test$class-1 
table(1,2)
#Now lets do PCA analysis!
pca=prcomp(train.inp)
eigen.sqr=round(100*sum(pca$sdev[1:2]^2)/sum(pca$sdev^2),digits=2)
par(mfrow=c(1,2))
barplot(pca$sdev^2,beside=TRUE,las=2)
legend(x='topright',paste('Var(2)=',eigen.sqr,'%'), bty= 'n')
plot(pca$x[,1:2],pch=c(1,4)[(train$class+1)],xlab='PC 1',
     ylab='PC 2', col=c('red', 'blue')[train$class+1])

#Reusable function to generate visual confusion matrix
conf.matrix =function(y,y.hat) {
  
  fd=table(y,y.hat)
  TN <- fd[1,1]
  TP <- fd[2,2]
  FN <- fd[2,1]
  FP <- fd[1,2]
  
  spe= TN/sum(TN,FP) # specificity, TN/(TN+FP)
  sen= TP/sum(TP,FN) # sensitivity, TP/(TP+FN)
  tot=sum(TN,TP)/sum(TN,TP,FN,FP) # total accuracy, (TN+TP)/((TN+TP+FP+FN)
  npp= TN/sum(TN,FN) # negative prediction power: TN/(TN+FN)
  ppp= TP/sum(TP,FP) # positive prediction power: TP/(TP+FP)
  
  #create output with easily indexed performance metrics
  confusion.matrix <- matrix(c(TN,FN,npp,FP,TP,ppp,spe,sen,tot),nrow=3)
  
  return(confusion.matrix)
}
#Reusable function to generate easily called performance metrics
metrics =function(conf.matrix) {
  TN <- conf.matrix[1,1]
  FP <- conf.matrix[1,2]
  spe <- conf.matrix[1,3]
  FN <- conf.matrix[2,1]
  TP <- conf.matrix[2,2]
  sen <- conf.matrix[2,3]
  npp <- conf.matrix[3,1]
  ppp <- conf.matrix[3,2]
  tot <- conf.matrix[3,3]
  
  #create output that can be used to call each of 9 total measurements by name
  metrics <- list("TN"=TN,"TP"=TP,"FN"=FN,"FP"=FP,"sensitivity"=sen,
                  "specificity"=spe,"NPP"=npp,"PPP"=ppp,"total.accuracy" = tot)
  return(metrics)
}


#apply LDA, store confusion matrix and print performance metrics
model.lda=lda(class~.,train)
y.lda=predict(model.lda,newdata=test.inp)$class
conf.matrix.lda <- conf.matrix(test$class, y.lda)
metrics.lda <- metrics(conf.matrix.lda)
metrics.lda

#apply brnn ANN function, store confusion matrix and print performance metrics
model.brnn=brnn(class~.,train,neurons=5)
y.brnn=predict(model.brnn,newdata=test.inp)
conf.matrix.brnn <- conf.matrix(test$class, y.brnn>0.5)
metrics.brnn <- metrics(conf.matrix.brnn)
metrics.brnn

#apply nnet ANN function, store confusion matrix and print performance metrics
model.nnet=nnet(class~.,train,size=5)
y.nnet=predict(model.nnet,newdata=test.inp)
conf.matrix.nnet <- conf.matrix(test$class, y.nnet>0.5)
metrics.nnet <- metrics(conf.matrix.nnet)
metrics.nnet

#evaluate total accuracy between ANN functions and select model displaying higher accuracy
if (metrics.brnn$total.accuracy > metrics.nnet$total.accuracy) {
  conf.matrix.ann <- conf.matrix.brnn
  metrics.ann <- metrics.brnn
  y.ann <- y.brnn
  print("brnn ANN function showing higher accuracy than nnet")
} else {
  conf.matrix.ann <- conf.matrix.nnet
  metrics.ann <- metrics.nnet
  y.ann <- y.nnet
  print("nnet ANN function showing higher accuracy than brnn")
}

#consensus analysis between linear LDA and non-linear ann models

#convert y.ann to binary values
y.ann[y.ann>0.5]<-1
y.ann[y.ann<0.5]<-0

tot.acc.venn = list(LDA=which(test$class == y.lda),
                    ANN=which((test$class== y.ann)))
venn.diagram(tot.acc.venn, file= 'tot.acc.tiff', cat.col = c("red","blue"),
             col= "black", fill=(c("red","blue")), lty= c(2,1))



