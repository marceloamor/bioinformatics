###Bioinformatics Workshop 2
##DATA VISUALIATION FUNCTIONS
#PCA, Sammon, and SOM networks!

#set working directory
setwd("~/Uni/4th year/Bioinfo/Machine Learning")

#upload necessary libraries
#install.packages("MASS")
#install.packages("kohonen")
library(MASS)
library(kohonen)

##upload data from wd and perform necessary transformations
#remove definitely essential genes
#log transformation
#fiddle with symbols list
X1 = read.csv('essential.genes.csv', row.names = 1)
X1 = X1[which(X1$sites>0),]
dim(X)

buf=unlist(strsplit(rownames(X1),split = '#'))
symbol = buf[seq(2,length(buf),2)]

X1=log(1+X1)

cut=2
label=1*(X1$mf<cut)+1

print(c('Essentialgene:',sum(label==2)),quote=F)
print(c('Non-essentialgene:',sum(label==1)),quote=F)

#create spot essential and nonessential gene
hemC.pos=which(symbol=='hemC')
kdpD.pos=which(symbol=='kdpD')
print(c(hemC.pos,label[hemC.pos]))
print(c(kdpD.pos,label[kdpD.pos]))

plot(density(X$mf),lwd=2,main = '', xlab= 'log(Mf)')
text(X$mf[c(hemC.pos,kdpD.pos)],0.02,symbol[c(hemC.pos,kdpD.pos)],srt=45)


##PCA mapping
pca=prcomp(X1)
ratio=round(100*sum(pca$sdev[1:2]^2)/sum(pca$sdev^2),digits = 2)
print(ratio)

names(pca$sdev)=paste('PC',1:3,sep='')
par(mfrow=c(1,2))
barplot(pca$sdev^2,las =2)
legend('topright', paste('Var(2)=',ratio),cex=.5)

plot(pca$x[,1:2],pch=c(1,3)[label],cex=0.7,col='gray')
points(pca$x[hemC.pos,1],pca$x[hemC.pos,2],pch=2,cex=2)
points(pca$x[kdpD.pos,1],pca$x[kdpD.pos,2],pch=3,cex=2)

plot(pca$x[,1:2],pch=c(1,3)[label],cex=0.7,col='gray')
text(pca$x[hemC.pos,1],pca$x[hemC.pos,2],'hemC')
text(pca$x[kdpD.pos,1],pca$x[kdpD.pos,2],'kdpD')

#SOM mapping
label[hemC.pos]=4
label[kdpD.pos]=3
som.col=c('gray','gray','black','black')[label]
som.pch=c(2,3,2,3)[label]
som.cex=c(0.7,0.7,2,2)[label]
som.lwd=c(1,1,2,2)[label]
som.model=som(as.matrix(X),grid=somgrid(20,20),rlen=1000)
plot(som.model,type='mapping',pch=som.pch,cex=som.cex,col=som.col,lwd=som.lwd,main='')

plot(som.model,type='codes',main='')
plot(som.model,type='counts',main='')
plot(som.model,type='dist.neighbours',main='')
plot(som.model,type='changes',main='')








