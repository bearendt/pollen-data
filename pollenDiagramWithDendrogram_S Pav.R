
# code for a pollen diagram with  dendrogram
# may need to delete rare taxa
# also bstick produces odd results iof chi distance.
# 4.22.2016 FDN
#edited 11.12.18 BA



setwd('P:/Sediment Analyses/Pollen/Data/SouthPavilion')


#South Pavilion Kitchen


# read the first three rows and transpose them into a data frame
labels <- read.table('montSouthPavilion2018.csv', header=F,  sep=",", skip = 2, nrows = 2, 
                     stringsAsFactors=F, colClasses='character')  
labels1 <- data.frame(t(labels),  stringsAsFactors=F)
labels2 <- labels1[2:nrow(labels1),] 
names(labels2) <- as.character(labels1[1,])
labels2$Elevation <- as.numeric(labels2$Elevation)

#pull the counts and species in
kdata<- read.table('montSouthPavilion2018.csv', header=F, 
                   skip= 6, sep=",", check.names=FALSE, na.strings="")
# transpose
kdata1 <- data.frame(t(kdata),stringsAsFactors = FALSE)
# make the first row into col names
colnames(kdata1)<-kdata1[1,]
# delete the first row of names
kdata1 <- kdata1[-c(1),]
# convert the counts from character to numeric 
kdata3<-data.frame(sapply(kdata1[,1:59],as.numeric))

#join labels and data together
Spav<-data.frame(labels2,kdata3)
# get rid of the taxon columns that are not id'd taxon counts
Spav <- Spav[ ,!names(Spav) %in% c('Indeterminate'  ,'Total.Pollen.Sum', 'Lycopodium.Tracers',
                                   'Concentration.Value', 'X.')]
#fill NAs with 0
Spav[is.na(Spav)]<-0
# compute depth
Spav$depth <- round((max(Spav$Elevation) - Spav$Elevation), 2)  
# get rid of elevation
Spav$Elevation <- NULL

par(mfrow=c(1,1))

library(rioja)
# the matrix of counts: you need to change the index for each dataset!!!! 
taxonMatrix<- as.matrix(Spav[,2:55])
# make contexts the rownames

# drop rows with samples < 100
samplesThatAreBigEnough <- rowSums(taxonMatrix) >= 100
taxonMatrix <- taxonMatrix[samplesThatAreBigEnough,]

# a column vector of depths
depth<- Spav$depth[samplesThatAreBigEnough]
Spav$depth <-NULL


# convert taxonmatrix of counts to percents

taxonMatrixP <- prop.table(taxonMatrix,1)*100 


library(analogue) # needed for chi.distances
# compute chi-squared diatnces and do constrained cluster analysis
diss <- distance(taxonMatrixP, method = 'chi.distance', dist=T)
clust <- chclust(diss, method="coniss")



bstick(clust)


x <- strat.plot(taxonMatrixP, 
                xRight = 1, # amount of space for dendrogram 
                yvar = depth, 
                y.rev=TRUE,
                scale.percent=TRUE, 
                title="South Pavilion", 
                ylabel="Depth (ft)",
                cex.xlabel =.75, #scaling for taxon labels
                srt.xlabel = 45,  #rotation angle for variable names.
                cex.axis =.5, 
                plot.poly=F,
                plot.line=F,
                plot.bar=T,
                lwd.bar =10,    
                col.bar ='blue',
                clust=clust)
addClustZone(x, clust, n=3, col= 'grey', lwd=3)

