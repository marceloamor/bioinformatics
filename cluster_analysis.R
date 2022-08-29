###workshop 3! 
##Cluster Analysis!

##alright were gonna start as usual with bringing in the essential
##gene information, loading libraries, and fiddling w the data#
library(VennDiagram)
library(kohonen)
#install.packages('mclust')
library(mclust)

setwd("~/Uni/4th year/Bioinfo/Machine Learning")

Z = readLines('DEG.genes.txt')
X = read.csv('essential.genes.csv', row.names = 1)
X = X[which(X$sites>0),]
X=log(1+X)

buf=unlist(strsplit(rownames(X),split = '#'))
symbol = buf[seq(2,length(buf),2)]

#####we start with Hierarchical Clustering! One feature at a time
mf.HC=hclust(dist(X$mf))
count.HC=hclust(dist(X$counts))
site.HC=hclust(dist(X$sites))

#then we plot all 3 side by side 
par(bg='grey90',mfrow=c(1,3))
plot(mf.HC,main='MF',xlab=''); grid(20,20,col='white')
plot(count.HC,main='count',xlab=''); grid(20,20,col='white')
plot(site.HC,main='site',xlab=''); grid(20,20,col='white')

###then we cut our trees to identify essential genes 
##using MF
mf.cluster=cutree(mf.HC,k=2)
cls1=which(mf.cluster==1)
cls2=which(mf.cluster==2)
mu1=mean(X$mf[cls1])
mu2=mean(X$mf[cls2])
ess.mf.cluster=ifelse(mu1 < mu2,1,2)
#get predicted genes and verified genes from this
predicted.mf=symbol[which(mf.cluster==ess.mf.cluster)]
verified.mf=predicted.mf[which(predicted.mf %in% Z)]
head(verified.mf)
rate.mf=length(verified.mf)/length(predicted.mf)
print(rate.mf)

#using counts
count.cluster=cutree(count.HC,k=2)
mu1=mean(X$counts[which(count.cluster==1)])
mu2=mean(X$counts[which(count.cluster==2)])
ess.count.cluster=ifelse(mu1 < mu2,1,2)

predicted.counts=symbol[which(count.cluster==ess.count.cluster)]
verified.counts=predicted.counts[which(predicted.counts %in% Z)]
head(verified.counts)
rate.counts=length(verified.counts)/length(predicted.counts)
print(rate.counts)
# using sites
site.cluster=cutree(site.HC,k=2)
mu1=mean(X$sites[which(site.cluster==1)])
mu2=mean(X$sites[which(site.cluster==2)])
ess.site.cluster=ifelse(mu1 < mu2,1,2)

predicted.sites=symbol[which(site.cluster==ess.count.cluster)]
verified.sites=predicted.sites[which(predicted.sites %in% Z)]
head(verified.sites)
rate.sites=length(verified.sites)/length(predicted.sites)
print(rate.sites)

##consensus analysis 
venn.list = list(MF=verified.mf, counts = verified.counts, 
                 sites = verified.sites)
venn.diagram(venn.list, filename = 'hclust.venn.tiff', cex=2,
             cat.cex=2, fill= c('red', 'blue', 'green'),
             col = c('red', 'blue', 'green'))


####HC analysis using all 3 features
hc.model = hclust(dist(X))
#get predicted genes
hc.cluster=cutree(hc.model,k=2)
mu1=mean(X$mf[which(hc.cluster==1)])
mu2=mean(X$mf[which(hc.cluster==2)])
ess.hc.cluster=ifelse(mu1 < mu2,1,2)
predicted.hc=symbol[which(hc.cluster==ess.hc.cluster)]
#get verified genes and rate
verified.hc= predicted.hc[which(predicted.hc %in% Z)]
rate.hc=length(verified.hc) / length(predicted.hc)
print(rate.hc)

####K means analysis using all 3 features
km.model = kmeans(X, centers = 2)
#get predicted genes
mu1=mean(X$mf[which(km.model$cluster==1)])
mu2=mean(X$mf[which(km.model$cluster==2)])
ess.km.cluster=ifelse(mu1 < mu2,1,2)
predicted.km=symbol[which(km.model$cluster==ess.km.cluster)]
# get verified essential genes for Kmeans model
verified.km=predicted.km[which(predicted.km %in% Z)]
rate.km=length(verified.km) / length(predicted.km)
head(verified.km)
print(rate.km)


#### Mixed Model for all 3 features
mm.model = Mclust(X, G=2)
mu1=mean(X$mf[which(mm.model$classification==1)])
mu2=mean(X$mf[which(mm.model$classification==2)])
ess.mm.cluster=ifelse(mu1 < mu2,1,2)
predicted.mm=symbol[which(mm.model$classification==ess.mm.cluster)]
# get verified essential genes for mixture model
verified.mm=predicted.mm[which(predicted.mm %in% Z)]
rate.mm=length(verified.mm) / length(predicted.mm)
head(verified.mm)
print(rate.mm)

##time to do some heckin consensus analysis
for.venn= list(Hclust=verified.hc, Kmeans=verified.km, Mclust= verified.mm)
venn.diagram(for.venn, filename = 'algorithms.venn.tiff',
             cex=2, cat.cex=2, fill = c('red','blue','green'),
             col= c('red','blue','green'))
model = som(as.matrix(X), grid= somgrid(30,30))
plot(model, type = 'mapping', labels= c('E', 'N')[mm.model$classification],
     col= c('red', 'blue')[mm.model$classification], cex=0.7, main='')



