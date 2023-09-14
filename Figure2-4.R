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
library(colorspace)
library(lme4)
library(lmerTest)
library(multcomp)

rm(list=ls())

#######   Read in data from master folder holding all sample folders - - - - - - - - - - - - - - - -

metadata2<-read.csv("metadataOhio.csv", fileEncoding = 'UTF-8-BOM')
metadata1<-metadata2[-73:-75,]
metadata<-subset(metadata1, use=="y")
metadata$Gen<-as.factor(metadata$Gen)
metadata$Symb<-as.factor(metadata$Symb)
#stat1<-read.csv("OhioTransposedPhotoOnly.csv", fileEncoding = 'UTF-8-BOM')
#stat<-read.csv("OhioTransposedPhotoOnly2.csv", fileEncoding = 'UTF-8-BOM')
stat<-read.csv("CF-annot-npq.csv", header=TRUE)
#stat<-read.csv("CF-annot-qm.csv", header=TRUE)
#stat<-read.csv("CF-annot-1-q.csv", header=TRUE)
#stat<-read.csv("CF-annot-npq.csv", header=TRUE)

#stat<-read.csv("TransOhioDataNAs.csv", header=TRUE)
rownames(stat)<-stat[,1]
rownames(stat)
colnames(stat)
stat<-stat[,-1]
dim(stat)

stat[] <- lapply(stat, function(x) as.numeric(as.character(x)))
#stat
sapply(stat, class)
head(stat)

stater<-t(stat)
stat1<-stat[match(metadata$Sample, rownames(stat)),]
rownames(stat1)
rownames(stat1)<-metadata$uniq[match(rownames(stat1), metadata$Sample)]
dim(stat1)
rownames(stat1)
stat<-t(stat1)
dim(stat)

#########Identify metrics of interest########
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
dim(metafluor)

chooserFluo<-subset(metafluor, metric=="Quant" | metric=="Sigma" | metric=="Connect" | metric=="Tau1" | metric=="Tau2" | metric=="NPQ" | metric=="qP" | metric=="ABQ")
dim(chooserFluo)
dim(stat)
stat<-stat[match(chooserFluo$uni, rownames(stat)),]
dim(stat)

write.csv(stat, "reanalysisData.csv")


rs<-rownames(stat)
cls<-colnames(stat)
Boo<-dim(stat)[[1]]
datum <- list()
dim(stat)

replicates<-metadata$Symb[match(colnames(stat), metadata$uniq)]

pvale<-0.05 #change based on bonferroni correction (divide by factor levels or just factor levels driving phenotype)
for (i in 1:Boo)
{	H<-rownames(stat)[i]	
	if(is.na(summary(aov(stat[i,]~replicates))[[1]][,5][1])==FALSE)
	{	if (as.matrix(shapiro.test(resid(aov(stat[i,]~replicates)))[[2]])>pvale)
		{	if(summary(aov(stat[i,]~replicates))[[1]][,5][1]<pvale)
			{	#print(summary(aov(t(stat[i,])~replicates))[[1]][,5][1])
				datum[[H]]<-((as.numeric(stat[i,]) - mean(as.numeric(stat[i,]), na.rm=TRUE))/sd(as.numeric(stat[i,]), na.rm=TRUE))
			}
		}
		else
		{	if(kruskal.test(stat[i,]~replicates)[[3]]<pvale)
			{	#print(kruskal.test(t(stat[i,])~replicates)[[3]])
				datum[[H]]<-((as.numeric(stat[i,]) - mean(as.numeric(stat[i,]), na.rm=TRUE))/sd(as.numeric(stat[i,]), na.rm=TRUE))
			}
		}
	}
}
allerie2<-data.matrix(do.call(rbind, datum))
colnames(allerie2)<-cls
allerie2<-na.omit(allerie2)
#summary(allerie2)
dim(allerie2)
#head(allerie2) #has dimensions

seqdata<-read.csv("SeqData727.csv", header=TRUE)
summary(seqdata)
data<-seqdata[,8:13]
rownames(data)<-paste(seqdata$SpecVarFrag)
summary(data)

# ss2<-t(allerie2) #transposes samples to rows
# ss3<-ss2[match(rownames(data), rownames(ss2)),] #tries to match qpcr or sequencing data to samples in allerie output, making unmatched ones NA
# allerie2<-t(na.omit(ss3)) #omits the NAs
# dim(allerie2)

#"euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
#"ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).


cb<-pvclust(allerie2, method.dist="canberra", method.hclust="ward.D", nboot=10) #change back to 1000 to run for real
par(mfrow=c(1,1), oma = c(3, 1, 1, 1), mar = c(1, .5, 1, 1))
plot(cb)

ca<-hclust(dist(allerie2, method="canberra"), method="ward.D")
caclusters<-cutree(ca, k=4)     ######### Chose how many sections you want k=....
cbclusters<-cutree(cb, k=4)     ######### Chose how many sections you want k=....

cb$hclust$order

dd<-flip_leaves(as.dendrogram(cb$hclust), c(43, 59, 58, 60, 45, 42, 44, 38, 40), c(46, 47, 48, 37, 39, 31, 41))
plot(dd)
order.dendrogram(dd)
dr<-flip_leaves(dd, c(37, 39), c(31, 41))
plot(dr)
order.dendrogram(dr)
de<-flip_leaves(dr, c(21, 19, 20), c(34, 9, 35, 7, 36))
plot(de)
order.dendrogram(de)
ds<-flip_leaves(de, c(22), c(49))
plot(ds)
order.dendrogram(ds)
da<-flip_leaves(ds, c(15), c(49, 22))
plot(da)
order.dendrogram(da)
daa<-flip_leaves(da, c(1,2,3), c(24, 50, 51, 49, 22, 15, 23, 28))
plot(daa)
order.dendrogram(daa)
das<-flip_leaves(daa, c(5, 6, 25, 26, 32, 4, 33, 8, 27), c(24, 50, 51, 49, 22, 15, 23, 28, 1, 2, 3))
plot(das)
order.dendrogram(das)

samp<-paste("DENDY811", ".pdf", sep="")
pdf(file = samp, width = 14, height = 8, bg="transparent")

gd<-as.dendrogram(das)
ca <- color_branches(ca, k = 4, col = c("#000000","#000000","#000000","#000000","#000000","#000000","#000000")) # add color to the lines
gd <- color_branches(gd, k = 4, col = c("black", "#ed1c24", "#21409a", "#00a14b")) # add color to the lines

par(mfrow=c(4,1), oma = c(3, 1, 1, 1), mar = c(1, .5, 1, 1))
plot(gd, print.pv=TRUE, hang=0.1, cex=0.5)
colsplot<-labels(gd)
colsplot2<-data[match(labels(gd), rownames(data)),]
data1<-colsplot2
da2<-colnames(data1)
order(da2)
barplot(t(data1), col=c("#67001f", "#f46d43", "#fdae61", "#fee090", "#d73027", "#abd9e9")[as.factor(da2)], border="white", space=0.02, axes=F, las=2, cex.names=0.005, legend=FALSE, args.legend = list(x = "topright", inset = c(.1, .1)))
axis(2, at=, col.axis="black", las=.5, cex.axis=1)

colsplot<-labels(gd)
colsplot3<-metadata$Gen[match(labels(gd), metadata$uniq)]
da1<-colsplot3

das1<-rowSums(data1)
barplot(das1, col=c("#543005", "#bf812d", "#f6e8c3", "#c7eae5", "#80cdc1", "#01665e", "#003c30")[as.factor(da1)], 
        border="white", space=0.02, axes=F, las=2, cex.names=0.5, legend.text=FALSE)
axis(2, at=, col.axis="black", las=0.5,cex.axis=1)

colsplot<-labels(gd)
colsplot4<-metadata$Tank[match(labels(gd), metadata$uniq)]
da3<-colsplot4

das1<-rowSums(data1)
barplot(das1, col=c("#636363","#bdbdbd")[as.factor(da3)], 
        border="white", space=0.02, axes=F, las=2, cex.names=0.5, legend.text=FALSE)
axis(2, at=, col.axis="black", las=0.5,cex.axis=1)


dev.off()

#factor(caclusters)

samp<-paste("HeatMap811", ".pdf", sep="")
pdf(file = samp, width = 7, height = 4, bg="transparent")

col_breaks = c(seq(-3.0,-0.5,length.out=115), seq(-0.4,0.4,length.out=200), seq(0.5,3.0,length.out=115))
my_palette <- colorRampPalette(c("#2166ac","white","#d73027"))(length(col_breaks)-1)

par(oma=c(2,0.5,0.5,0.5))
heatmap.2(allerie2, col=my_palette, breaks=col_breaks, margins = c(2,3), density.info="none", trace="none", dendrogram = "both", symm=FALSE, symkey=FALSE, symbreaks=FALSE, sepwidth=F, na.color = "pink", Colv=as.dendrogram(gd), Rowv=as.dendrogram(ca), labRow=ca$labels, RowSideColors=c("#d9d9d9", "#000000", "#bdbdbd", "#969696","#737373","#525252","#252525")[caclusters], scale="none", key=T, key.xlab="Z-score", key.title=NA, cexRow=0.00007, cexCol=0.001, srtCol=90, lmat=rbind(c(5,6,4), c(3,1,2)), lhei=c(2, 8), lwid=c(0.25,0.1, 3))
dev.off()


colsplot<-labels(gd)
colsplot2<-data[match(labels(gd), rownames(data)),]
data1<-colsplot2

samp<-paste("Symbionts811", ".pdf", sep="")
pdf(file = samp, width = 14, height = 4, bg="transparent")
da2<-colnames(data1)
order(da2)
par(mfrow=c(1,1), oma = c(4, 2, 5, 5), mar = c(1, .5, 3, 2))
# Get the stacked barplot
barplot(t(data1), col=c("#67001f", "#f46d43", "#fdae61", "#fee090", "#d73027", "#abd9e9")[as.factor(da2)], border="white", space=0.02, axes=F, las=2, cex.names=0.005, legend=FALSE, args.legend = list(x = "topright", inset = c(.1, .1)))
axis(4, at=, col.axis="black", las=.5, cex.axis=1)

dev.off()




colsplot<-labels(gd)
colsplot3<-metadata$Gen[match(labels(gd), metadata$uniq)]
da1<-colsplot3

samp<-paste("coralBar811", ".pdf", sep="")
pdf(file = samp, width = 14, height = 3, bg="transparent")
par(mar = c(1, .5, 3, 2))
das1<-rowSums(data1)
barplot(das1, col=c("#543005", "#bf812d", "#f6e8c3", "#c7eae5", "#80cdc1", "#01665e", "#003c30")[as.factor(da1)], 
        border="white", space=0.02, axes=F, las=2, cex.names=0.5, legend.text=FALSE)
axis(4, at=, col.axis="black", las=0.5,cex.axis=1)

dev.off()





colsplot<-labels(gd)
colsplot4<-metadata$Tank[match(labels(gd), metadata$uniq)]
da3<-colsplot4

samp<-paste("BarInvOut811", ".pdf", sep="")
pdf(file = samp, width = 14, height = 3, bg="transparent")
par(mar = c(1, .5, 3, 2))
das1<-rowSums(data1)
barplot(das1, col=c("#636363","#bdbdbd")[as.factor(da3)], 
        border="white", space=0.02, axes=F, las=2, cex.names=0.5, legend.text=FALSE)
axis(4, at=, col.axis="black", las=0.5,cex.axis=1)

dev.off()

##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################

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

caclust<-data.frame(caclusters)

ca1<-subset(caclust, caclusters=="1")
ca2<-subset(caclust, caclusters=="2")
ca3<-subset(caclust, caclusters=="3")
ca4<-subset(caclust, caclusters=="4")



dim(ca1) #147, Quant 72, qP 23, NPQ 5
dim(ca2) #167, 44 Sigma, 56 Tau1
dim(ca3) #blue	 100 Tau2
dim(ca4) #167, 44 NPQ, 33 qP, 20 ABQ, 3 connect




caa<-ca1
ca1Rich<-metafluor[match(rownames(caa), metafluor$uni),]
ca1Rich1<-data.frame(aggregate(ca1Rich[,2], by=list(ca1Rich$metric), FUN=NROW))
ca1m<-data.frame(ca1Rich1, "percent"=ca1Rich1$x/dim(ca1Rich)[1])
ca1Rich2<-data.frame(aggregate(ca1Rich[,5], by=list(ca1Rich$col), FUN=NROW))
ca1c<-data.frame(ca1Rich2, "percent"=ca1Rich2$x/dim(ca1Rich)[1])
ca1Rich3<-data.frame(aggregate(ca1Rich[,4], by=list(ca1Rich$stage), FUN=NROW))
ca1t<-data.frame(ca1Rich3, "percent"=ca1Rich3$x/dim(ca1Rich)[1])
ca1m
ca1c
ca1t

###NPQ,ABQ,Tau1,Quant,Tau2,
##########################
##########################
##########################
##########################
cbclust<-data.frame(cbclusters)

cb1<-subset(cbclust, cbclusters=="1") #grey
cb2<-subset(cbclust, cbclusters=="2") #red
cb3<-subset(cbclust, cbclusters=="3") #blue
cb4<-subset(cbclust, cbclusters=="4") #green

allerie4<-t(stat)
dim(allerie4)

phen1<-allerie4[match(rownames(cb1),rownames(allerie4)),]
dim(phen1)
rownames(phen3)

phen2<-allerie4[match(rownames(cb2),rownames(allerie4)),]
dim(phen2)

phen3<-allerie4[match(rownames(cb3),rownames(allerie4)),]
dim(phen3)

phen4<-allerie4[match(rownames(cb4),rownames(allerie4)),]
dim(phen4)

##########Phenotrace###########
##########Phenotrace###########
##########Phenotrace###########
##########Phenotrace###########
##########Phenotrace###########
##########Phenotrace###########
rownames(phen4)
stupid<-list()
for(I in 1:ncol(phen1))
{
	stupid[[I]]<-data.frame("phys"=colnames(phen1)[I], "mean"=mean(as.numeric(phen1[,I]), na.rm=TRUE), "sd"=sd(as.numeric(phen1[,I]), na.rm=TRUE)/sqrt(nrow(phen1)))
}
alleriephen1<-data.frame(do.call(rbind, stupid))

stupid<-list()
for(I in 1:ncol(phen2))
{
	stupid[[I]]<-data.frame("phys"=colnames(phen2)[I], "mean"=mean(as.numeric(phen2[,I]), na.rm=TRUE), "sd"=sd(as.numeric(phen2[,I]), na.rm=TRUE)/sqrt(nrow(phen2)))
}
alleriephen2<-data.frame(do.call(rbind, stupid))

stupid<-list()
for(I in 1:ncol(phen3))
{
	stupid[[I]]<-data.frame("phys"=colnames(phen3)[I], "mean"=mean(as.numeric(phen3[,I]), na.rm=TRUE), "sd"=sd(as.numeric(phen3[,I]), na.rm=TRUE)/sqrt(nrow(phen3)))

}
alleriephen3<-data.frame(do.call(rbind, stupid))

stupid<-list()
for(I in 1:ncol(phen4))
{
	stupid[[I]]<-data.frame("phys"=colnames(phen4)[I], "mean"=mean(as.numeric(phen4[,I]), na.rm=TRUE), "sd"=sd(as.numeric(phen4[,I]), na.rm=TRUE)/sqrt(nrow(phen4)))
}
alleriephen4<-data.frame(do.call(rbind, stupid))

P1<-data.frame(metafluor, alleriephen1[,2:3])
P2<-data.frame(metafluor, alleriephen2[,2:3])
P3<-data.frame(metafluor, alleriephen3[,2:3])
P4<-data.frame(metafluor, alleriephen4[,2:3])

metric<-c("Quant", "Sigma", "Connect", "Tau1", "Tau2", "NPQ", "qP", "ABQ")
sets<-c(0.5,7.5,1,4000,7000000,1.25,1,0.5)
osets<-c(0,0,0,0,0,-0.25,0,0)
TT<-c(1,2,6,7,4,5)


samp<-paste("Profiles-Ohio811-2", ".pdf", sep="")
pdf(file = samp, width = 8.5, height = 2, bg="transparent")

par(mfrow=c(1,4), oma = c(0.1, 5.9, 1.1, 0.5), mar = c(1, 0.1, 0.1, 0.1))

#for (G in 1:6)
#{  
	G<-3
	X<-TT[G]
	print(metric[X])
	ggg<-list(P3,P4,P2,P1)
	for (I in 1:4)
	{	col1<-subset(ggg[[I]], metric==metric[X] & col=="1")
		col2<-subset(ggg[[I]], metric==metric[X] & col=="2")
		col3<-subset(ggg[[I]], metric==metric[X] & col=="3")
		col4<-subset(ggg[[I]], metric==metric[X] & col=="4")
		col5<-subset(ggg[[I]], metric==metric[X] & col=="5")
		
		plot(col1$stage, col1$mean, type="l", lwd=1.25, col="purple", ylim=c(osets[X],sets[X]), axes=F)
		par(new=TRUE)
		plot(col1$stage, col1$mean+col1$sd, type="l", lwd=0.5, col="purple", ylim=c(osets[X],sets[X]), axes=F)
		par(new=TRUE)
		plot(col1$stage, col1$mean-col1$sd, type="l", lwd=0.5, col="purple", ylim=c(osets[X],sets[X]), axes=F)
		par(new=TRUE)
		plot(col2$stage, col2$mean, type="l", lwd=1.25, col="blue", ylim=c(osets[X],sets[X]), axes=F)
		par(new=TRUE)
		plot(col2$stage, col2$mean+col2$sd, type="l", lwd=0.5, col="blue", ylim=c(osets[X],sets[X]), axes=F)
		par(new=TRUE)
		plot(col2$stage, col2$mean-col2$sd, type="l", lwd=0.5, col="blue", ylim=c(osets[X],sets[X]), axes=F)
		par(new=TRUE)
		plot(col3$stage, col3$mean, type="l", lwd=1.25, col="light blue", ylim=c(osets[X],sets[X]), axes=F)
		par(new=TRUE)
		plot(col3$stage, col3$mean+col3$sd, type="l", lwd=0.5, col="light blue", ylim=c(osets[X],sets[X]), axes=F)
		par(new=TRUE)
		plot(col3$stage, col3$mean-col3$sd, type="l", lwd=0.5, col="light blue", ylim=c(osets[X],sets[X]), axes=F)
		par(new=TRUE)
		plot(col4$stage, col4$mean, type="l", lwd=1.25, col="cyan", ylim=c(osets[X],sets[X]), axes=F)
		par(new=TRUE)
		plot(col4$stage, col4$mean+col4$sd, type="l", lwd=0.5, col="cyan", ylim=c(osets[X],sets[X]), axes=F)
		par(new=TRUE)
		plot(col4$stage, col4$mean-col4$sd, type="l", lwd=0.5, col="cyan", ylim=c(osets[X],sets[X]), axes=F)
		par(new=TRUE)
		plot(col5$stage, col5$mean, type="l", lwd=1.25, col="green", ylim=c(osets[X],sets[X]), axes=F)
		par(new=TRUE)
		plot(col5$stage, col5$mean+col5$sd, type="l", lwd=0.5, col="green", ylim=c(osets[X],sets[X]), axes=F)
		par(new=TRUE)
		plot(col5$stage, col5$mean-col5$sd, type="l", lwd=0.5, col="green", ylim=c(osets[X],sets[X]), axes=F)
		if(I==1)
		{
		if(G==1)
		{ axis(2, at=c(0,0.25,0.5), col.axis="black", las=2, cex.axis=1.5)}
		else if(G==2)
		{ axis(2, at=c(0,2.5,7.5), col.axis="black", las=2, cex.axis=1.5)}
		else if(G==3)
		{ axis(2, at=c(0,0.6,1.2), col.axis="black", las=2, cex.axis=1.5)}
		else if(G==4)
		{ axis(2, at=c(0,0.5,1), col.axis="black", las=2, cex.axis=1.5)}
		else if(G==5)
		{ axis(2, at=c(0,2000,4000), col.axis="black", las=2, cex.axis=1.5)}
		else if(G==6)
		{ axis(2, at=c(0,3500000,7000000), col.axis="black", las=2, cex.axis=1.5)}
		}
		box(col="black")
		grid(col="gray")
		if(G==6)
		{ axis(1, at=c(0,10,30), col.axis="black", las=1, cex.axis=1.5)}
	}
	dev.off()

#}
#dev.off()

##########Stats###########
##########Stats###########
##########Stats###########
##########Stats###########
##########Stats###########
##########Stats###########


bbbb<-data.frame(metafluor, t(phen3), t(phen4), t(phen2), t(phen1))
dim(bbbb)

colnames(bbbb)

PP1<-data.frame("nn"=phen3[,1], "bin"=1)
PP2<-data.frame("nn"=phen4[,1], "bin"=2)
PP3<-data.frame("nn"=phen2[,1], "bin"=3)
PP4<-data.frame("nn"=phen1[,1], "bin"=4)

Factorss<-rbind(PP1,PP2,PP3,PP4)

factor<-Factorss[,2]
fact<-as.factor(Factorss[,2])
#isolating quant for each excitation wavelength
Bquant1<-subset(bbbb, metric=="Quant" & col=="1")
Bquant2<-subset(bbbb, metric=="Quant" & col=="2")
Bquant3<-subset(bbbb, metric=="Quant" & col=="3")
Bquant4<-subset(bbbb, metric=="Quant" & col=="4")
Bquant5<-subset(bbbb, metric=="Quant" & col=="5")

quant1<-Bquant5[,-1:-3]
quant<-quant1[,-2]
quan<-t(quant[,-1])
colnames(quan)<-quant[,1]
qua<-data.frame("fact"=fact, "sample"=rownames(quan), quan)
active7<-na.omit(melt(qua, id.vars=c("fact", "sample"), variable.name="cat", value.name="change"))
fit1 <- lmer(as.numeric(change) ~ fact + (1|sample) + (1/cat), data=active7, REML=FALSE)
anova(fit1) #for the quantum yield or whatever metric we are using, we are looking at excitation color X. Are there sig diff across 4 profiles
summary(glht(fit1, linfct=mcp(fact="Tukey")), test=adjusted("bonferroni"))#since yes, which differ from one another?


Bsig1<-subset(bbbb, metric=="Sigma" & col=="1")
Bsig2<-subset(bbbb, metric=="Sigma" & col=="2")
Bsig3<-subset(bbbb, metric=="Sigma" & col=="3")
Bsig4<-subset(bbbb, metric=="Sigma" & col=="4")
Bsig5<-subset(bbbb, metric=="Sigma" & col=="5")

quant1<-Bsig5[,-1:-3]
quant<-quant1[,-2]
quan<-t(quant[,-1])
colnames(quan)<-quant[,1]
qua<-data.frame("fact"=fact, "sample"=rownames(quan),quan)
active7<-na.omit(melt(qua, id.vars=c("fact", "sample"), variable.name="cat", value.name="change"))
fit1 <- lmer(as.numeric(change) ~ fact + (1|sample) + (1/cat), data=active7, REML=FALSE)
anova(fit1)
summary(glht(fit1, linfct=mcp(fact="Tukey")), test=adjusted("bonferroni"))





Bqp1<-subset(bbbb, metric=="qP" & col=="1")
Bqp2<-subset(bbbb, metric=="qP" & col=="2")
Bqp3<-subset(bbbb, metric=="qP" & col=="3")
Bqp4<-subset(bbbb, metric=="qP" & col=="4")
Bqp5<-subset(bbbb, metric=="qP" & col=="5")

quant1<-Bqp5[,-1:-3]
quant<-quant1[,-2]
quan<-t(quant[,-1])
colnames(quan)<-quant[,1]
qua<-data.frame("fact"=fact, "sample"=rownames(quan),quan)
active7<-na.omit(melt(qua, id.vars=c("fact", "sample"), variable.name="cat", value.name="change"))
fit1 <- lmer(as.numeric(change) ~ fact + (1|sample) + (1/cat), data=active7, REML=FALSE)
anova(fit1)
summary(glht(fit1, linfct=mcp(fact="Tukey")), test=adjusted("bonferroni"))





Bnpq1<-subset(bbbb, metric=="NPQ" & col=="1")
Bnpq2<-subset(bbbb, metric=="NPQ" & col=="2")
Bnpq3<-subset(bbbb, metric=="NPQ" & col=="3")
Bnpq4<-subset(bbbb, metric=="NPQ" & col=="4")
Bnpq5<-subset(bbbb, metric=="NPQ" & col=="5")

quant1<-Bnpq5[,-1:-3]
quant<-quant1[,-2]
quan<-t(quant[,-1])
colnames(quan)<-quant[,1]
qua<-data.frame("fact"=fact, "sample"=rownames(quan),quan)
active7<-na.omit(melt(qua, id.vars=c("fact", "sample"), variable.name="cat", value.name="change"))
fit1 <- lmer(as.numeric(change) ~ fact + (1|sample) + (1/cat), data=active7, REML=FALSE)
anova(fit1)
summary(glht(fit1, linfct=mcp(fact="Tukey")), test=adjusted("bonferroni"))




B1tau1<-subset(bbbb, metric=="Tau1" & col=="1")
B2tau1<-subset(bbbb, metric=="Tau1" & col=="2")
B3tau1<-subset(bbbb, metric=="Tau1" & col=="3")
B4tau1<-subset(bbbb, metric=="Tau1" & col=="4")
B5tau1<-subset(bbbb, metric=="Tau1" & col=="5")

quant1<-B5tau1[,-1:-3]
quant<-quant1[,-2]
quan<-t(quant[,-1])
colnames(quan)<-quant[,1]
qua<-data.frame("fact"=fact, "sample"=rownames(quan),quan)
active7<-na.omit(melt(qua, id.vars=c("fact", "sample"), variable.name="cat", value.name="change"))
fit1 <- lmer(as.numeric(change) ~ fact + (1|sample) + (1/cat), data=active7, REML=FALSE)
anova(fit1)
summary(glht(fit1, linfct=mcp(fact="Tukey")), test=adjusted("bonferroni"))



B1tau2<-subset(bbbb, metric=="Tau2" & col=="1")
B2tau2<-subset(bbbb, metric=="Tau2" & col=="2")
B3tau2<-subset(bbbb, metric=="Tau2" & col=="3")
B4tau2<-subset(bbbb, metric=="Tau2" & col=="4")
B5tau2<-subset(bbbb, metric=="Tau2" & col=="5")

quant1<-B5tau2[,-1:-3]
quant<-quant1[,-2]
quan<-t(quant[,-1])
colnames(quan)<-quant[,1]
qua<-data.frame("fact"=fact, "sample"=rownames(quan),quan)
active7<-na.omit(melt(qua, id.vars=c("fact", "sample"), variable.name="cat", value.name="change"))
fit1 <- lmer(as.numeric(change) ~ fact + (1|sample) + (1/cat), data=active7, REML=FALSE)
anova(fit1)
summary(glht(fit1, linfct=mcp(fact="Tukey")), test=adjusted("bonferroni"))






##########Statsacrosslightcolor###########
##########Statsacrosslightcolor###########
##########Statsacrosslightcolor###########
##########Statsacrosslightcolor###########
##########Statsacrosslightcolor###########
##########Statsacrosslightcolor###########
#we want to run ANOVA for each metric, across light colors within each phenotype

#isolating metrics for each excitation wavelength
B1<-data.frame(metafluor, t(phen3))
B2<-data.frame(metafluor, t(phen4))
B3<-data.frame(metafluor, t(phen2))
B4<-data.frame(metafluor, t(phen1))

#subset based on a single metric--will need to change each time I run to something different
Bquant1<-subset(B4, metric=="Tau2")

#remove the first few columns from the subset--will need to change subset used to run all 5
quant1<-Bquant1[,-1:-3] 
quantmelt<-melt(quant1, id.vars=c("stage", "col")) #melts data into long format
quantmelt$col<-as.factor(quantmelt$col)
quantmelt$value<-as.numeric(quantmelt$value)
fit1 <- lmer(value ~ col + (1|variable), data=quantmelt, REML=FALSE) #makes excitation wavelength and time important, makes coral ID random
anova(fit1) 
summary(fit1) #doesn't show info for each individual light color so maybe I did this wrong
summary(glht(fit1, linfct=mcp(col="Tukey")), test=adjusted("bonferroni")) #does not work

################
################
################
################

metaAve<-read.csv("all60Cellular.csv", fileEncoding = 'UTF-8-BOM')
head(metaAve)

phen1c<-data.frame(metaAve[match(rownames(phen3),metaAve$X),], phen="1")
dim(phen1c)

phen2c<-data.frame(metaAve[match(rownames(phen4),metaAve$X),], phen="2")
dim(phen2c)

phen3c<-data.frame(metaAve[match(rownames(phen2),metaAve$X),], phen="3")
dim(phen3c)

phen4c<-data.frame(metaAve[match(rownames(phen1),metaAve$X),], phen="4")
dim(phen4c)

CellPheno<-rbind(phen1c,phen2c,phen3c,phen4c)
CellPheno<-CellPheno[,-2:-4]
CellPheno<-CellPheno[,-7]
CellPheno<-CellPheno[,-8:-12]
colnames(CellPheno)

meanCP<-aggregate(CellPheno[,2:7], by=list(CellPheno$phen), FUN=mean, na.rm=TRUE)
sdCP<-aggregate(CellPheno[,2:7], by=list(CellPheno$phen), FUN=sd, na.rm=TRUE)



pdf(file = "phenocellBar.pdf", width = 3, height = 6, bg="transparent")

par(mfrow=c(3,2), oma = c(2, 4, 2, 4), mar = c(1, 0.1, 0.5, 0.7))

XX<-6 #FSC
bar3<-barplot(meanCP[,XX], axes=FALSE, col="gray", ylim=c(0,7)) 
arrows(bar3,(meanCP[,XX]+(sdCP[,XX]/sqrt(1))), bar3,(meanCP[,XX]-(sdCP[,XX]/sqrt(1))), length=0.01, angle=90, code=3,  col="black", lwd=1)
axis(2, at=, col.axis="black", las=2)
box(col = "black")

shapiro.test(CellPheno[,XX]) #SSC
kruskal.test(CellPheno[,XX] ~ CellPheno$phen)
pairwise.wilcox.test(CellPheno[,XX], CellPheno$phen, paired=FALSE, p.adjust.method="bonferroni")

XX<-7 #SSC
bar3<-barplot(meanCP[,XX], axes=FALSE, col="gray", ylim=c(0,5)) 
arrows(bar3,(meanCP[,XX]+(sdCP[,XX]/sqrt(1))), bar3,(meanCP[,XX]-(sdCP[,XX]/sqrt(1))), length=0.01, angle=90, code=3,  col="black", lwd=1)
axis(4, at=, col.axis="black", las=2)
box(col = "black")

shapiro.test(CellPheno[,XX]) #SSC
aovs<-aov(CellPheno[,XX] ~ CellPheno$phen)
summary(aovs)
TukeyHSD(aovs)

XX<-5 #chla
bar1<-barplot(meanCP[,XX], axes=FALSE, col="gray", ylim=c(0,40)) 
arrows(bar1,(meanCP[,XX]+(sdCP[,XX]/sqrt(1))), bar1,(meanCP[,XX]-(sdCP[,XX]/sqrt(1))), length=0.01, angle=90, code=3,  col="black", lwd=1)
axis(2, at=, col.axis="black", las=2)
box(col = "black")

shapiro.test(CellPheno[,XX]) #Chla
kruskal.test(CellPheno[,XX] ~ CellPheno$phen)
pairwise.wilcox.test(CellPheno[,XX], CellPheno$phen, paired=FALSE, p.adjust.method="bonferroni")

XX<-3 #NP
bar2<-barplot(meanCP[,XX], axes=FALSE, col="gray", ylim=c(0,45)) 
arrows(bar2,(meanCP[,XX]+(sdCP[,XX]/sqrt(1))), bar2,(meanCP[,XX]-(sdCP[,XX]/sqrt(1))), length=0.01, angle=90, code=3,  col="black", lwd=1)
axis(4, at=, col.axis="black", las=2)
box(col = "black")

shapiro.test(CellPheno[,XX]) #CN
aovs<-aov(CellPheno[,XX] ~ CellPheno$phen)
summary(aovs)
TukeyHSD(aovs)

XX<-4 #CP
bar3<-barplot(meanCP[,XX], axes=FALSE, col="gray", ylim=c(0,350)) 
arrows(bar3,(meanCP[,XX]+(sdCP[,XX]/sqrt(1))), bar3,(meanCP[,XX]-(sdCP[,XX]/sqrt(1))), length=0.01, angle=90, code=3,  col="black", lwd=1)
axis(2, at=, col.axis="black", las=2)
box(col = "black")

shapiro.test(CellPheno[,XX]) #CN
aovs<-aov(CellPheno[,XX] ~ CellPheno$phen)
summary(aovs)
TukeyHSD(aovs)

XX<-2 #CN
bar3<-barplot(meanCP[,XX], axes=FALSE, col="gray", ylim=c(0,12)) 
arrows(bar3,(meanCP[,XX]+(sdCP[,XX]/sqrt(1))), bar3,(meanCP[,XX]-(sdCP[,XX]/sqrt(1))), length=0.01, angle=90, code=3,  col="black", lwd=1)
axis(4, at=, col.axis="black", las=2)
box(col = "black")

shapiro.test(CellPheno[,XX]) #NP
aovs<-aov(CellPheno[,XX] ~ CellPheno$phen)
summary(aovs)
TukeyHSD(aovs)

dev.off()























