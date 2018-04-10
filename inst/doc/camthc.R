## ----setup, include = FALSE------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
knitr::knit_hooks$set(webgl = rgl::hook_webgl)

## ---- eval=FALSE-----------------------------------------------------------
#  rCAM <- CAM(data, K = 2:5, thres.low = 0.30, thres.high = 0.95)

## ---- echo=TRUE------------------------------------------------------------
library(CAMTHC)

data(ratMix3) 
#ratMix3$X: X matrix cotaining mixtue expression profiless to be analyzed
#ratMix3$A: ground truth A matrix cotaining proportions
#ratMix3$S: ground truth S matrix containing subpopulation-specific expression profiles

data <- ratMix3$X
#10000 genes * 21 tissues
#meet the input data requirements

## ---- echo=TRUE, message=FALSE---------------------------------------------
rCAM <- CAM(data, K = 2:5, thres.low = 0.30, thres.high = 0.95, seed = 111)
#CAM return three sub results: 
#rCAM$PrepResult contains details corresponding to data preprocessing.
#rCAM$MGResult contains details corresponding to marker gene clusters detection.
#rCAM$ASestResult contains details corresponding to A and S matrix estimation.

## ---- echo=TRUE------------------------------------------------------------
Aest <- rCAM$ASestResult$'3'$Aest
Sest <- rCAM$ASestResult$'3'$Sest

## ---- echo=TRUE------------------------------------------------------------
m <- which(names(rCAM$ASestResult) == '3')
Aest <- rCAM$ASestResult[[m]]$Aest
Sest <- rCAM$ASestResult[[m]]$Sest

## ---- echo=TRUE------------------------------------------------------------
MGlist <- MGsforA(rCAM, K = 3) #for three subpopulations

## ---- echo=TRUE------------------------------------------------------------
MGstat <- MGstatistic(data, Aest, boot.alpha = 0.05, nboot = 1000)
MGlist.FC <- lapply(1:3, function(x) 
    rownames(MGstat)[MGstat$idx == x & MGstat$OVE.FC > 10])
MGlist.FCboot <- lapply(1:3, function(x) 
    rownames(MGstat)[MGstat$idx == x & MGstat$OVE.FC.alpha > 10])

## ---- fig.height=6, fig.width=6--------------------------------------------
simplexplot(data, Aest, MGlist)
simplexplot(data, Aest, MGlist.FC)
simplexplot(data, Aest, MGlist.FCboot)

## ---- fig.height=6, fig.width=6--------------------------------------------
simplexplot(data, Aest, MGlist.FCboot, corner.order = c(2,1,3), 
            data.col = "blue", corner.col = c("red","orange","green"))

## ---- webgl=TRUE-----------------------------------------------------------
library(rgl)
Xp <- data %*% t(rCAM$PrepResult$W)
plot3d(Xp[, 1:3], col='gray', size=3, 
       xlim=range(Xp[,1]), ylim=range(Xp[,2:3]), zlim=range(Xp[,2:3]))
abclines3d(0,0,0, a=diag(3), col="black")
for(i in seq_along(MGlist)){
    points3d(Xp[MGlist[[i]], 1:3], col= rainbow(3)[i], size = 8)
}

## ---- webgl=TRUE, message=FALSE--------------------------------------------
library(rgl)
clear3d()
Xproj <- XWProj(data, rCAM$PrepResult$W)
Xp <- Xproj[,-1]
plot3d(Xp[, 1:3], col='gray', size=3, 
       xlim=range(Xp[,1:3]), ylim=range(Xp[,1:3]), zlim=range(Xp[,1:3]))
abclines3d(0,0,0, a=diag(3), col="black")
for(i in seq_along(MGlist)){
    points3d(Xp[MGlist[[i]], 1:3], col= rainbow(3)[i], size = 8)
}

## ---- echo=TRUE------------------------------------------------------------
MDL(rCAM)$mdls

## ---- fig.height=6, fig.width=6--------------------------------------------
plot(MDL(rCAM), data.term = TRUE)

## ---- echo=TRUE------------------------------------------------------------
cor(Aest, ratMix3$A)
cor(Sest, ratMix3$S)

## ---- echo=TRUE------------------------------------------------------------
unlist(lapply(1:3, function(k) {
    k.match <- which.max(cor(Aest[,k], ratMix3$A));
    mgk <- MGlist.FCboot[[k]];
    corr <- cor(Sest[mgk, k], ratMix3$S[mgk, k.match]);
    names(corr) <- colnames(ratMix3$A)[k.match];
    corr}))


## ---- echo=TRUE------------------------------------------------------------
#Data preprocession
rPrep <- CAMPrep(data, thres.low = 0.30, thres.high = 0.95, seed = 111)

#Marker gene cluster detection with a fixed K
rMGC <- CAMMGCluster(3, rPrep)

#A and S matrix estimation
rASest <- CAMASest(rMGC, rPrep, data)

#Obtain A and S matrix
Aest <- rASest$Aest
Sest <- rASest$Sest

#obtain marker gene list detected by CAM and used for A estimation
MGlist <- MGsforA(PrepResult = rPrep, MGResult = rMGC)

#obtain a full list of marker genes 
MGstat <- MGstatistic(data, Aest, boot.alpha = 0.05, nboot = 1000)
MGlist.FC <- lapply(1:3, function(x) 
    rownames(MGstat)[MGstat$idx == x & MGstat$OVE.FC > 10])
MGlist.FCboot <- lapply(1:3, function(x) 
    rownames(MGstat)[MGstat$idx == x & MGstat$OVE.FC.alpha > 10])

## ---- echo=TRUE, message=FALSE---------------------------------------------
#clustering
library(apcluster)
apres <- apclusterK(negDistMat(r=2), t(data),  K = 10) 
#You can also use apcluster(), but need to make sure the number of clusters is large than potentional subpopulation number.

data.clusterMean <- lapply(slot(apres,"clusters"), 
                           function(x) rowMeans(data[, x, drop = FALSE]))
data.clusterMean <- do.call(cbind, data.clusterMean)

rCAM <- CAM(data.clusterMean, K = 2:5, thres.low = 0.30, thres.high = 0.95, 
            cores = 30, seed = 111)
# or rPrep <- CAMPrep(data.clusterMean, thres.low = 0.30, thres.high = 0.95, seed = 111)

## ---- echo=TRUE------------------------------------------------------------
Sest <- rCAM$ASestResult$'3'$Sest
MGlist <- MGsforA(rCAM, K = 3)
Aest <- AfromMarkers(data, MGlist)

## ---- echo=TRUE, message=FALSE---------------------------------------------
#download data and phenotypes
library(GEOquery)
gsm<- getGEO('GSE11058')
pheno <- pData(phenoData(gsm[[1]]))$characteristics_ch1
mat <- exprs(gsm[[1]])
mat <- mat[-grep("^AFFX", rownames(mat)),]
mat.aggre <- sapply(unique(pheno), function(x) rowMeans(mat[,pheno == x]))
data <- mat.aggre[,5:8]

#running CAM
rCAM <- CAM(data, K = 4, thres.low = 0.70, thres.high = 0.95, 
            cores = 30, seed = 111)
Aest <- rCAM$ASestResult$'4'$Aest
Aest

#Use ground truth A to validate CAM-estimated A matrix
Atrue <- matrix(c(2.50, 0.50, 0.10, 0.02,
              1.25, 3.17, 4.95, 3.33,
              2.50, 4.75, 1.65, 3.33,
              3.75, 1.58, 3.30, 3.33), 4,4,
              dimnames = list(c("MixA", "MixB", "MixC","MixD"),
                               c("Jurkat", "IM-9", "Raji", "THP-1")))
Atrue <- Atrue / rowSums(Atrue)
Atrue
cor(Aest, Atrue)

## ---- echo=TRUE, eval=FALSE------------------------------------------------
#  Aest <- AfromMarkers(data, MGlist)
#  #MGlist is a list of vectors, each of which contains known markers for one subpopulation

## ---- echo=TRUE------------------------------------------------------------
data <- ratMix3$X
S <- ratMix3$S #known S matrix

pMGstat <- MGstatistic(S, c("Liver","Brain","Lung"))
pMGlist.FC <- lapply(c("Liver","Brain","Lung"), function(x) 
    rownames(pMGstat)[pMGstat$idx == x & pMGstat$OVE.FC > 10])

Aest <- AfromMarkers(data, pMGlist.FC)

## ---- echo=TRUE------------------------------------------------------------
data <- ratMix3$X
A <- ratMix3$A #known A matrix

Sest <- t(NMF::.fcnnls(A, t(data))$coef)
MGstat <- MGstatistic(data, A)

## ---- echo=TRUE, eval=FALSE------------------------------------------------
#  Aest <- AfromMarkers(data, MGlist)
#  #MGlist is a list of vectors, each of which contains known markers and/or CAM-detected markers for one subpopulation

