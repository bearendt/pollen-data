#Combining pollen data from multiple csv files, transforming and then running a CA
#edited for South Pavilion BA 11.13.2018

setwd('P:/Sediment Analyses/Pollen/Data/SouthPavilion')

# read the first three rows and transpose them into a data frame
labels <- read.table('montSouthPavilion2018.csv', header=F, skip=2, sep=",", nrows = 2, 
                     stringsAsFactors=F, colClasses='character')  
labels1 <- data.frame(t(labels),  stringsAsFactors=F)
labels2 <- labels1[2:nrow(labels1),] 
names(labels2) <- as.character(labels1[1,])
labels2$Elevation <- as.numeric(labels2$Elevation)

#pull the counts and species in
kdata<- read.table('montSouthPavilion2018.csv', header=F, skip= 6, sep=",", check.names=FALSE, na.strings="")
kdata1 <- data.frame(t(kdata),stringsAsFactors = FALSE)
colnames(kdata1)<-kdata1[1,]
kdata1 <- kdata1[-c(1),]
kdata3<-data.frame(sapply(kdata1[,1:58],as.numeric))

#join labels and data together
SPav<-data.frame(labels2,kdata3)


#delete Elevation and sample number collumns
SPav$Elevation<- NULL


library(reshape2)
#melt the data into three columns
SPavMelted <- melt(SPav, id= 'Context')


#######################################################################
# put the  data frames together

#remove Indeterminate
SPavMelted1 <- subset(SPavMelted, ! SPavMelted$variable %in% c('Total.Pollen.Sum', 'Indeterminate', 
                                        'Lycopodium.Tracers', 'Concentration.Value'))

contexts <- as.data.frame(table(SPavMelted1$Context), stringsAsFactors =F)
contexts <- contexts[order(contexts$Var1),]
# print out the context names and look for contexts that need to be deleted
contexts

# get rid of PZ and outlier contexts
SPavMelted2 <- subset(SPavMelted1, ! SPavMelted1$Context %in% c('2588F'))

#SPavMelted1 <- subset(SPavMelted, ! SPavMelted$variable %in%
  #                  c('Indeterminate', 'Ailanthus'))
# transpose into a matrix
SPDataTransp <- acast(SPavMelted2, Context ~ variable,
      sum,
      value.var='value')
#remove first row
SPDataTransp = SPDataTransp[-1, ]
#replace NA with 0
SPDataTransp[is.na(SPDataTransp)] <- 0

# check the sample sizes
rowSums(SPDataTransp)

# check the number of conrtext in which taxa occur
nContextsWithTaxa <- colSums(SPDataTransp >0) 

taxaToKeep <- nContextsWithTaxa > 1

SPDataTransp1 <-SPDataTransp[,taxaToKeep]

countMat <- SPDataTransp1


#################################################################
require(ca)


caResult <-ca(countMat)


broken.stick <- function(p)
  # Compute the expected values of the broken-stick distribution for 'p' pieces.
  # Example: broken.stick.out.20 = broken.stick(20)
  #             Pierre Legendre, April 2007
{
  result = matrix(0,p,2)
  colnames(result) = c("j","E(j)")
  for(j in 1:p) {
    E = 0
    for(x in j:p) E = E+(1/x)
    result[j,1] = j
    result[j,2] = E/p
  }
  return(result)
}

ProportionVariance<- prop.table(caResult$sv^2)
bs<-broken.stick(length(caResult$sv))
plot(1:length(caResult$sv), ProportionVariance, type="b", xlab="Dimension Order", 
     ylab= "Proportion of Inertia", cex=2, cex.axis=1.5, cex.lab=1.5)
lines(bs[,1],bs[,2], col="sky blue", lwd=2)


par(mar=c(5,8,4,2)) #increase left margins

par(mfrow=c(1,1))
barplot(prop.table(caResult$sv^2), axis.lty=1, names.arg=(1:length(caResult$sv)),
        xlab='Dimension', 
        ylab= 'Proportion of Inertia',
        cex.lab=1.5,
        cex.axis=1.5,
        cex.names=1.5)





summary(caResult)

op<-par(mfrow=c(1,1))




#######################################################
# scatter plots of dim scores
library (maptools)
plot(caResult$rowcoord[,1],caResult$rowcoord[,2],pch=21, cex=3, bg= adjustcolor("red", alpha.f=.5),
     xlab="Dimension 1", 
     ylab="Dimension 2", 
   #  xlim= c(-3, 3),
     cex.lab=1.5, cex.axis=1.5)
abline(h=0,v=0, col='grey')    
pointLabel(caResult$rowcoord[,1],caResult$rowcoord[,2], caResult$rownames)

plot(caResult$colcoord[,1],caResult$colcoord[,2], pch=21, cex=3, bg= adjustcolor("red", alpha.f=.5),
     xlab="Dimension 1", 
     ylab="Dimension 2",
     #xlim=c(-4.5,4.5),
     cex.lab=1.5, cex.axis=1.5)
abline(h=0,v=0, col= 'grey')
pointLabel(caResult$colcoord[,1],caResult$colcoord[,2], caResult$colnames, 
     pos=1,
     cex =.75)


# scatter plots of dim scores
library (maptools)
plot(caResult$rowcoord[,1],caResult$rowcoord[,3],pch=21, cex=3, bg= adjustcolor("red", alpha.f=.5),
     xlab="Dimension 1", 
     ylab="Dimension 3", 
     #  xlim= c(-3, 3),
     cex.lab=1.5, cex.axis=1.5)
abline(h=0,v=0, col='grey')    
pointLabel(caResult$rowcoord[,1],caResult$rowcoord[,3], caResult$rownames)

plot(caResult$colcoord[,1],caResult$colcoord[,3], pch=21, cex=3, bg= adjustcolor("red", alpha.f=.5),
     xlab="Dimension 1", 
     ylab="Dimension 3",
     #xlim=c(-4.5,4.5),
     cex.lab=1.5, cex.axis=1.5)
abline(h=0,v=0, col= 'grey')
pointLabel(caResult$colcoord[,1],caResult$colcoord[,3], caResult$colnames, 
           pos=1,
           cex =.75)
##########################################################
#New ggplot

# scatter plots of dim scores
rowscores <- data.frame(caResult$rowcoord[,1], caResult$rowcoord[,2])
colnames(rowscores) <- c("Dim1", "Dim2")

colscores <- data.frame(caResult$colcoord[,1], caResult$colcoord[,2])
colnames(colscores) <- c('Dim1', 'Dim2')

require(ggplot2)
library(ggrepel)
Plot1 <- ggplot(rowscores, aes(x=rowscores$Dim1,y=rowscores$Dim2))+
  geom_point(shape=21, size=5, colour="black", fill="cornflower blue")+
  #geom_text(aes(label=CA_MCD_Phase1$unit),vjust=-.6, cex=5)+
  geom_text_repel(aes(label=rownames(rowscores)), cex=6) +
  theme_classic()+
  labs(x="Dimension 1", y="Dimension 2")+
  theme(plot.title=element_text(size=rel(2.25), hjust=0.5),axis.title=element_text(size=rel(1.75)),
        axis.text=element_text(size=rel(1.5)))

Plot2 <- ggplot(colscores, aes(x=colscores$Dim1,y=colscores$Dim2))+
  geom_point(shape=21, size=5, colour="black", fill="cornflower blue")+
  #geom_text(aes(label=CA_MCD_Phase1$unit),vjust=-.6, cex=5)+
  geom_text_repel(aes(label=rownames(colscores)), cex=6) +
  theme_classic()+
  labs(x="Dimension 1", y="Dimension 2")+
  theme(plot.title=element_text(size=rel(2.25), hjust=0.5),axis.title=element_text(size=rel(1.75)),
        axis.text=element_text(size=rel(1.5)))

library(ggpubr)
ggarrange(Plot1, Plot2)

