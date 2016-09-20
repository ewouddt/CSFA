# Project: Connectivity
# 
# Author: lucp8394
###############################################################################


## Function to check overlapping biclusters between ref and query compounds
# Input: fabia result + thresholdL + thresholdZ + index of ref compounds


fabia.overlap <- function(resFAB,ind.ref,type="loadings",thresZ=0.5,thresL=NULL){
	
	bicResult <- extractBic(resFAB,thresZ=thresZ,thresL=thresL)
	ind.bc.load <- bicResult$numn[,1]
	ind.bc.factor <- bicResult$numn[,2]
	
	overlap.BC <- c(1:length(ind.bc.load))
	n.overlap.ref <- c(1:length(ind.bc.load))
	n.overlap.query <- c(1:length(ind.bc.load))
	
	if(type=="loadings"){
		ind.total <- c(1:dim(resFAB@L)[1])
		ind.query <- ind.total[-ind.ref]
		
		for(i.bc in 1:length(ind.bc.load)){
			ind.bc.load.i <- ind.bc.load[[i.bc]]
			nref.temp <- length(intersect(ind.ref,ind.bc.load.i))
			nquery.temp <- length(intersect(ind.query,ind.bc.load.i))
			
			if((length(intersect(ind.ref,ind.bc.load.i))>=1)  & (length(intersect(ind.query,ind.bc.load.i))>=1)){
				overlap.BC[i.bc] <- i.bc
				n.overlap.ref[i.bc] <- nref.temp
				n.overlap.query[i.bc] <- nquery.temp
			}
			
#			# OR
#			if((sum(ind.ref %in% ind.bc.load.i)>=1 ) & (sum(ind.query %in% ind.bc.load.i)>=1 )  ){
#				overlap.BC <- c(overlap.BC,i.bc)
#			}
					
		}
	}
	
	if(type=="factors"){
		ind.total <- c(1:dim(resFAB@Z)[2])
		ind.query <- ind.total[-ind.ref]
		
		for(i.bc in 1:length(ind.bc.factor)){
			ind.bc.factor.i <- ind.bc.factor[[i.bc]]
			nref.temp <- length(intersect(ind.ref,ind.bc.factor.i))
			nquery.temp <- length(intersect(ind.query,ind.bc.factor.i))
			
			if((length(intersect(ind.ref,ind.bc.factor.i))>=1)  & (length(intersect(ind.query,ind.bc.factor.i))>=1)){
				overlap.BC[i.bc] <- i.bc
				n.overlap.ref[i.bc] <- nref.temp
				n.overlap.query[i.bc] <- nquery.temp
			}
			
		}
	}
	
	#bicResult$numn[1,1]  # row = bicluster ; col = 1 (compounds), col = 2 (genes)
	
	out <- data.frame(Bicluster=overlap.BC,Number.Ref=n.overlap.ref,Number.Query=n.overlap.query)
	
	return(out)
}


## Function to compute the weighted data (as done in MFA)   (author: Nolen)
#library(pls)
getWeightedDat <- function(Mat1,Mat2,scale.unit=TRUE){
	if(scale.unit){
		Mat1 <- stdize(Mat1,center=TRUE,scale=TRUE)
		Mat2 <- stdize(Mat2,center=TRUE,scale=TRUE) 
		resPCA1 <- PCA(Mat1,graph=FALSE)
		resPCA2 <- PCA(Mat2,graph=FALSE)
	}
	else{
		resPCA1 <- PCA(Mat1,scale.unit=FALSE,graph=FALSE)
		resPCA2 <- PCA(Mat2,scale.unit=FALSE,graph=FALSE)
	}
	getWeight1 <- resPCA1$eig[1,1]
	
	getWeight2 <- resPCA2$eig[1,1]
	Mat1.w <- Mat1*(1/sqrt(getWeight1))
	Mat2.w <- Mat2*(1/sqrt(getWeight2))
	dat <- list()
	dat[[1]] <- data.matrix(cbind(Mat1.w, Mat2.w)) 
	dat[[2]] <- 1/sqrt(getWeight1)
	dat[[3]] <- 1/sqrt(getWeight2)
	dat[[4]] <- resPCA1
	dat[[5]] <- resPCA2
	names(dat) <- c("weightedData", "weight1","weight2", "resPCA1", "resPCA2")
	return(dat)
}


