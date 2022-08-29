###Bioinformatics Workshop 1
##DENSITY FUNCTIONS

#set working directory
setwd("~/Uni/4th year/Bioinfo/Machine Learning")


#uploadin my data
X0 = read.csv('essential.genes.csv', row.names=1)
Z = readLines('DEG.genes.txt')
dim(X0)
head(X0)

###Definite genes. Find and remove.
definite = which(X0$sites==0)
length(definite)

X=X0[-definite,]
dim(X)

#extracting gene symbols from X
head(rownames(X))
buf = unlist(strsplit(rownames(X), split = '#'))
head(buf)
symbol = buf[seq(2,length(buf),2)]
head(symbol)

#log transforming of MF
MF = log(1+X$mf)

par(mfrow=c(1,3))

#hist(X$mf, nclass=30, main= '', xlab= 'raw MF')
hist(MF, nclass=30, main= '', xlab= 'log MF')

##begin identifying essential genes using MF
plot(density(MF), lwd=2, main = '', xlab='log(MF)')

predicted.MF = which(MF<2)
predicted.MF.symbol = symbol[predicted.MF]
verified.MF.id = which(predicted.MF.symbol %in% Z)
verified.MF = predicted.MF.symbol[verified.MF.id]
head(verified.MF)
rate.MF=length(verified.MF)/length(predicted.MF.symbol)
print(rate.MF)


## identify essential genes based on counts
counts=log(1+X$counts)
par(mfrow=c(1,1))
hist(counts,nclass=30,xlab='log(counts)',main='',prob=TRUE)
lines(density(counts),lwd=2)
predicted.counts=which(counts<5)
predicted.counts.symbol=symbol[predicted.counts]
verified.counts.id=which(predicted.counts.symbol %in% Z)
verified.counts=predicted.counts.symbol[verified.counts.id]
rate.counts=length(verified.counts)/length(predicted.counts.symbol
)
print(rate.counts)


# identify essential genes based on sites
sites=log(1+X$sites)
par(mfrow=c(1,1))
hist(counts,nclass=30,xlab='log(sites)',main='',prob=TRUE)
lines(density(counts),lwd=2)
predicted.sites=which(sites<3)
predicted.sites.symbol=symbol[predicted.sites]
verified.sites.id=which(predicted.sites.symbol %in% Z)
verified.sites=predicted.sites.symbol[verified.sites.id]
rate.sites=length(verified.sites)/length(predicted.sites.symbol)
print(rate.sites)
print(c(rate.MF,rate.counts,rate.sites))

#consensus analysis
install.packages('VennDiagram')
library(VennDiagram)
inp =list(MF=verified.MF, counts=verified.counts, sites=verified.sites)
venn.diagram(inp,filename='essential.venn.tiff', cex=2, cat.cex=2)















