library(ggplot2)
library(gplots)
library(gridExtra)
library(car)
library(vegan)
library(reshape)
library(reshape2)
library(broom)
library(pvclust)
library(dendextend)

rm(list=ls())

#######   Read in data from master folder holding all sample folders - - - - - - - - - - - - -

col<-c(0,1,2,3,4)
stage<-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34)
metric<-c("Quant", "Sigma", "Connect", "Tau1", "Tau2", "NPQ", "qP", "ABQ")
metrix<-c(1,2,3,4,5,6,7,8)
h<-1
metaphys<-list()
for( I in 1:5)
{
	for (TT in 1:34)
	{
		for (FF in 1:8)
		{
			uni<-paste(metric[FF], stage[TT], col[I], sep=".")
			metaphys[[h]]<-data.frame("uni"=uni, "metric"=metric[FF], "metrix"=metrix[FF], "stage"=stage[TT], "col"=(col[I]+1))
			h<-h+1
		}
	}
}
metafluor<-data.frame(do.call(rbind, metaphys))


#######
#######
#######
#######
#######
#######

netdat<-read.csv("all60Cellular.csv", header=TRUE)
dim(netdat)
netdat2<-read.csv("reanalysisData.csv", header=TRUE)
rownames(netdat2)<-netdat2[,1]
netdat2<-t(netdat2[,-1])
dim(netdat2)

netdat3<-cbind(netdat,netdat2)[,-1]
dim(netdat3)

FinDa<-netdat3
dim(FinDa)

stupid<-list()
for(I in 1:ncol(FinDa))
{
	stupid[[I]]<-data.frame(as.numeric(FinDa[,I]))
}
FinDa2<-data.frame(do.call(cbind, stupid))
dim(FinDa2)
colnames(FinDa2)<-colnames(FinDa)
rownames(FinDa2)<-rownames(FinDa)

metadata2<-read.csv("metadataOhio.csv", fileEncoding = 'UTF-8-BOM')
metadata1<-metadata2[-73:-75,]
metadata<-subset(metadata1, use=="y")
metadata$Gen<-as.factor(metadata$Gen)
metadata$Symb<-as.factor(metadata$Symb)
metadata$combi<-paste(metadata$Species, metadata$Tank, sep="-")

metafeta<-metadata[match(rownames(FinDa2), metadata$uniq),]

FinDad<-aggregate(FinDa2, by=list(metafeta$combi), FUN=mean, na.rm=TRUE)
FinDa3<-FinDad[,-1]
rownames(FinDa3)<-FinDad[,1]

#FinDa3<-FinDa2
FinDa3<-FinDa3[,-12:-15]
FinDa3<-FinDa3[,-9]
FinDa3<-FinDa3[,-1:-3]
FinDa3<-FinDa3[,-7]
dim(FinDa3)
#Finner<-FinDa3[,4:8]

colnames(FinDa3[,1:15])

library(igraph)

cor_mat <- as.matrix(as.dist(cor(FinDa3, method = c("pearson"), use="complete.obs")))
cor_g <- graph_from_adjacency_matrix(cor_mat, mode='undirected', weighted = TRUE)
cor_edge_list <- as_data_frame(cor_g, 'edges')
dim(cor_edge_list)
summary(cor_edge_list)

possibles<-colnames(netdat)[2:12]
finset<-cor_edge_list
dim(finset)
finer<-list()
gg<-1
for (t in 1:dim(finset)[1])
{	if(finset[t,1] %in% possibles)
	{	if(!(finset[t,2] %in% possibles))
		{
			finer[[gg]]<-cbind(finset[t,], "lm"=1, "phys"=2)
			gg<-gg+1
		}
	}
	if(finset[t,2] %in% possibles)
	{	if(!(finset[t,1] %in% possibles))
		{
			finer[[gg]]<-cbind(finset[t,], "phys"=2, "lm"=1)
			gg<-gg+1
		}
	}
}
resulter1<-data.frame(do.call(rbind, finer))
dim(resulter1)
resulter<-subset(resulter1, abs(weight)>0.6)
dim(resulter)

dim(subset(resulter, from=="CN"))
dim(subset(resulter, from=="CP"))
dim(subset(resulter, from=="NP"))
dim(subset(resulter, from=="ChlA"))
dim(subset(resulter, from=="SSC"))
dim(subset(resulter, from=="Size"))

res1<-data.frame("vert"=resulter$from, "type"=resulter$phys, "weight"=resulter$weight)
res2<-data.frame("vert"=resulter$to, "type"=resulter$lm, "weight"=resulter$weight)
res3<-rbind(res1,res2)
res4<-res3[!duplicated(res3$vert),]

wavlen<-data.frame("uni"=(possibles), "metric"="cellular", "metrix"=10,"stage"=35,"col"=6) 
colnames(wavlen)[1]<-"uni"
aller<-rbind(metafluor, wavlen)
alls<-aller[match(res4$vert, aller$uni),]
res5<-data.frame(res4, alls)

new_g<-graph_from_data_frame(resulter[,1:3],F)
g6<-simplify(new_g, remove.multiple = T, remove.loops = T)

V(g6)$light<-res5$stage
V(g6)$type<-res5$type
V(g6)$clust<-res5$col
V(g6)$shape<-res5$metrix
E(g6)$color<-c(1,2)[as.factor(sign(E(g6)$weight))]
E(g6)$weight<-((((abs(E(g6)$weight)-range(abs(E(g6)$weight))[1])/(range(abs(E(g6)$weight))[2]-range(abs(E(g6)$weight))[1])))+0.75)^3


ll <- layout_with_dh(g6)

samp<-paste("pearson0.676-Ave-metric", ".pdf", sep="")
pdf(file = samp, width = 8, height = 8, bg="transparent")
par(mfrow=c(1,1), oma = c(4, 5.9, 1.1, 0.5), mar = c(1, 0.1, 0.1, 0.1))
plot(g6, layout=ll, vertex.shape=c("circle", "circle", "circle", "circle", "circle", "circle", "circle", "circle", "circle", "square")[V(g6)$shape], vertex.size=c(7, 11)[V(g6)$type], vertex.color=c("#01665e","#f5f5f5","#35978f","#80cdc1","#c7eae5","#f6e8c3","#dfc27d","#bf812d","#8c510a","grey")[V(g6)$shape], vertex.frame.color="black", vertex.frame.cex=1, vertex.label.color="black", vertex.label.cex=0.00001, vertex.label.dist=0.01, edge.width=E(g6)$weight, edge.curved=0.2, edge.color=c("orange","black")[E(g6)$color])
dev.off()

samp<-paste("pearson0.676-Ave-col", ".pdf", sep="")
pdf(file = samp, width = 6, height = 6, bg="transparent")
par(mfrow=c(1,1), oma = c(4, 5.9, 1.1, 0.5), mar = c(1, 0.1, 0.1, 0.1))
plot(g6, layout=ll, vertex.shape=c("circle", "circle", "circle", "circle", "circle", "circle", "circle", "circle", "circle", "square")[V(g6)$shape], vertex.size=c(7, 11)[V(g6)$type], vertex.color=c("#810f7c","#253494","#1d91c0","#7fcdbb","#41ab5d","grey")[V(g6)$clust], vertex.frame.color="black", vertex.frame.cex=1, vertex.label.color="black", vertex.label.cex=0.00001, vertex.label.dist=0.01, edge.width=E(g6)$weight, edge.curved=0.2, edge.color=c("orange","black")[E(g6)$color])
dev.off()

samp<-paste("pearson0.676-Ave-light", ".pdf", sep="")
pdf(file = samp, width = 6, height = 6, bg="transparent")
par(mfrow=c(1,1), oma = c(4, 5.9, 1.1, 0.5), mar = c(1, 0.1, 0.1, 0.1))
plot(g6, layout=ll, vertex.shape=c("circle", "circle", "circle", "circle", "circle", "circle", "circle", "circle", "circle", "square")[V(g6)$shape], vertex.size=c(7, 11)[V(g6)$type], vertex.color=c("black","black","#fee391","#fee391","#fee391","#fee391","#fee391","#fee391","#fee391","#fee391","#ffffe5","#ffffe5","#ffffe5","#ffffe5","#ffffe5","#ffffe5","#ffffe5","#ffffe5", "#fe9929", "#fe9929", "#fe9929", "#fe9929", "#fe9929", "#fe9929", "#fe9929", "#fe9929", "#636363", "#636363", "#636363", "#636363", "#636363", "#636363", "#636363", "#636363", "grey")[V(g6)$light], vertex.frame.color="black", vertex.frame.cex=1, vertex.label.color="black", vertex.label.cex=0.00001, vertex.label.dist=0.01, edge.width=E(g6)$weight, edge.curved=0.2, edge.color=c("orange","black")[E(g6)$color])
dev.off()

############
############
############
############
############

samp<-paste("F6-NP-Quant", ".pdf", sep="")
pdf(file = samp, width = 6, height = 4, bg="transparent")
par(mfrow=c(1,1), oma = c(2, 6, 0.5, 0.1), mar = c(0.1, 0.1, 0.1, 0.1))

#Sym<-metadata$Symb[match(rownames(FinDa3),metadata$uniq)]
Sym<-c("C3","C3","C3","C21","C21","C15","C1","C26","C15","C15","C15","C15","C15","D1","D1","C1","D1","D1","D1","D1")
highlow<-c(1,2,1,1,1,1,1,1,2,1,2,1,1,1,1,1,1,2,1,1)
y<-FinDa3$Quant.16.1
x<-FinDa3$NP
cor.test(x, y, method=c("pearson"))
plot(x, y, pch=c(17,19)[highlow], cex=2.5, col=c("#67001f", "#f46d43", "#fdae61", "#fee090", "#d73027", "#abd9e9")[as.factor(Sym)], axes=FALSE)
axis(1, at=, col.axis="black", las=1, cex.axis=2)
axis(2, at=, col.axis="black", las=1, cex.axis=2)
box(col="black")
abline(lm(y ~ x), col="black", lwd=3)
summary(lm(y ~ x))
dev.off()

############
############

samp<-paste("F6-NP-Tau2", ".pdf", sep="")
pdf(file = samp, width = 6, height = 4, bg="transparent")
par(mfrow=c(1,1), oma = c(2, 6, 0.5, 0.1), mar = c(0.1, 0.1, 0.1, 0.1))

Sym<-c("C3","C3","C3","C21","C21","C15","C1","C26","C15","C15","C15","C15","C15","D1","D1","C1","D1","D1","D1","D1")
highlow<-c(1,2,1,1,1,1,1,1,2,1,2,1,1,1,1,1,1,2,1,1)
y<-FinDa3$Tau2.6.2
x<-FinDa3$NP
cor.test(x, y, method=c("pearson"))
plot(x, y, pch=c(17,19)[highlow], cex=2.5, col=c("#67001f", "#f46d43", "#fdae61", "#fee090", "#d73027", "#abd9e9")[as.factor(Sym)], ylim=c(0.5e+06,6.5e+06), axes=FALSE)
axis(1, at=, col.axis="black", las=1, cex.axis=2)
axis(2, at=, col.axis="black", las=1, cex.axis=2)
box(col="black")
abline(lm(y ~ x), col="black", lwd=3)
summary(lm(y ~ x))
dev.off()

############
############

samp<-paste("F6-CP-NPQ", ".pdf", sep="")
pdf(file = samp, width = 6, height = 4, bg="transparent")
par(mfrow=c(1,1), oma = c(2, 6, 0.5, 0.1), mar = c(0.1, 0.1, 0.1, 0.1))

Sym<-c("C3","C3","C3","C21","C21","C15","C1","C26","C15","C15","C15","C15","C15","D1","D1","C1","D1","D1","D1","D1")
highlow<-c(1,2,1,1,1,1,1,1,2,1,2,1,1,1,1,1,1,2,1,1)
y<-FinDa3$qP.20.4
x<-FinDa3$CP
cor.test(x, y, method=c("pearson"))
plot(x, y, pch=c(17,19)[highlow], cex=2.5, col=c("#67001f", "#f46d43", "#fdae61", "#fee090", "#d73027", "#abd9e9")[as.factor(Sym)], axes=FALSE, ylim=c(0.4,1))
axis(1, at=, col.axis="black", las=1, cex.axis=2)
axis(2, at=, col.axis="black", las=1, cex.axis=2)
box(col="black")
abline(lm(y ~ x), col="black", lwd=3)
summary(lm(y ~ x))
dev.off()

############
############

samp<-paste("F6-CP-Connect", ".pdf", sep="")
pdf(file = samp, width = 6, height = 4, bg="transparent")
par(mfrow=c(1,1), oma = c(2, 6, 0.5, 0.1), mar = c(0.1, 0.1, 0.1, 0.1))

Sym<-c("C3","C3","C3","C21","C21","C15","C1","C26","C15","C15","C15","C15","C15","D1","D1","C1","D1","D1","D1","D1")
highlow<-c(1,2,1,1,1,1,1,1,2,1,2,1,1,1,1,1,1,2,1,1)
y<-FinDa3$Connect.15.3
x<-FinDa3$CP
cor.test(x, y, method=c("pearson"))
plot(x, y, pch=c(17,19)[highlow], cex=2.5, col=c("#67001f", "#f46d43", "#fdae61", "#fee090", "#d73027", "#abd9e9")[as.factor(Sym)], axes=FALSE)
axis(1, at=, col.axis="black", las=1, cex.axis=2)
axis(2, at=, col.axis="black", las=1, cex.axis=2)
box(col="black")
abline(lm(y ~ x), col="black", lwd=3)
summary(lm(y ~ x))
dev.off()

############
############

samp<-paste("F6-SSC-Quant", ".pdf", sep="")
pdf(file = samp, width = 6, height = 4, bg="transparent")
par(mfrow=c(1,1), oma = c(2, 6, 0.5, 0.1), mar = c(0.1, 0.1, 0.1, 0.1))

Sym<-c("C3","C3","C3","C21","C21","C15","C1","C26","C15","C15","C15","C15","C15","D1","D1","C1","D1","D1","D1","D1")
highlow<-c(1,2,1,1,1,1,1,1,2,1,2,1,1,1,1,1,1,2,1,1)
y<-FinDa3$Quant.28.0
x<-FinDa3$SSC
cor.test(x, y, method=c("pearson"))
plot(x, y, pch=c(17,19)[highlow], cex=2.5, col=c("#67001f", "#f46d43", "#fdae61", "#fee090", "#d73027", "#abd9e9")[as.factor(Sym)], ylim=c(0.15,0.45), axes=FALSE)
axis(1, at=, col.axis="black", las=1, cex.axis=2)
axis(2, at=, col.axis="black", las=1, cex.axis=2)
box(col="black")
abline(lm(y ~ x), col="black", lwd=3)
summary(lm(y ~ x))
dev.off()

############
############

samp<-paste("F6-Size-ABQ", ".pdf", sep="")
pdf(file = samp, width = 6, height = 4, bg="transparent")
par(mfrow=c(1,1), oma = c(2, 6, 0.5, 0.1), mar = c(0.1, 0.1, 0.1, 0.1))

Sym<-c("C3","C3","C3","C21","C21","C15","C1","C26","C15","C15","C15","C15","C15","D1","D1","C1","D1","D1","D1","D1")
highlow<-c(1,2,1,1,1,1,1,1,2,1,2,1,1,1,1,1,1,2,1,1)
y<-FinDa3$ABQ.5.0
x<-FinDa3$Size
cor.test(x, y, method=c("pearson"))
plot(x, y, pch=c(17,19)[highlow], cex=2.5, col=c("#67001f", "#f46d43", "#fdae61", "#fee090", "#d73027", "#abd9e9")[as.factor(Sym)], axes=FALSE)
axis(1, at=, col.axis="black", las=1, cex.axis=2)
axis(2, at=, col.axis="black", las=1, cex.axis=2)
box(col="black")
abline(lm(y ~ x), col="black", lwd=3)
summary(lm(y ~ x))
dev.off()

############
############

samp<-paste("F6-CN-Connect", ".pdf", sep="")
pdf(file = samp, width = 6, height = 4, bg="transparent")
par(mfrow=c(1,1), oma = c(2, 6, 0.5, 0.1), mar = c(0.1, 0.1, 0.1, 0.1))

Sym<-c("C3","C3","C3","C21","C21","C15","C1","C26","C15","C15","C15","C15","C15","D1","D1","C1","D1","D1","D1","D1")
highlow<-c(1,2,1,1,1,1,1,1,2,1,2,1,1,1,1,1,1,2,1,1)
y<-FinDa3$Connect.16.3
x<-FinDa3$CN
cor.test(x, y, method=c("pearson"))
plot(x, y, pch=c(17,19)[highlow], cex=2.5, col=c("#67001f", "#f46d43", "#fdae61", "#fee090", "#d73027", "#abd9e9")[as.factor(Sym)], axes=FALSE)
axis(1, at=, col.axis="black", las=1, cex.axis=2)
axis(2, at=, col.axis="black", las=1, cex.axis=2)
box(col="black")
abline(lm(y ~ x), col="black", lwd=3)
summary(lm(y ~ x))
dev.off()

############
############
############
############
############
############
############