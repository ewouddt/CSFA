# Project: CSFA
# 
# Author: lucp8394
###############################################################################


# MAKE A FUNCTION TO DISPLAY GENES PROFILES BASED ON TOP COMPOUNDS AND GENES (SEE MCF7)
CSprofiles <- function(data,ref_index,gene.select,cmpd.select,profile.type=c("gene","cmpd"),cmpd.loadings,gene.scores,component.plot,gene.thresP,gene.thresN,basefilename,plot.type,thresP.col,thresN.col,main.base){
	
	## Plot-in and -out functions
	plot.in <- function(plot.type,name){
		if(plot.type=="pdf"){pdf(name)}
		if(plot.type=="device"){dev.new()}
		if(plot.type=="sweave"){}
	}
	plot.out <- function(plot.type){if(plot.type=="pdf"){dev.off()}}
	
	gene.index <- gene.select
	cmpd.index <- cmpd.select
	
	if(profile.type=="gene"){
		# Gene Profiles
		
		new.order <- order(abs(cmpd.loadings[,component.plot]),decreasing=TRUE) # Other compounds are ordered by their absolute loadings
		order.others <- setdiff(new.order,c(ref_index,cmpd.index)) 
		
		data.new <- data[gene.index,c(ref_index,cmpd.index,order.others),drop=FALSE]
		
		gene.colors <- distinctColorPalette(length(gene.index))
		
		plot.in(plot.type,paste0(basefilename,"_GeneProfiles.pdf"))
		plot(0,type="n",axes=FALSE,xlim=c(1,dim(data.new)[2]),ylim=c(min(data.new),max(data.new)),ylab="Gene Expression",xlab="Compounds",main=paste0(main.base," - Gene Profiles for ",paste0(rownames(data.new),collapse=", ")))
		for(i.gene in c(1:length(gene.index))){
			lines(c(1:dim(data.new)[2]),data.new[i.gene,],col=gene.colors[i.gene])
		}
		axis(side=2)
		axis(side=1,labels=FALSE)
		mtext(colnames(data.new),cex=0.7,side=1,las=3,line=1,at=c(1:dim(data.new)[2]),col=c(rep("blue",length(ref_index)),rep("red",length(cmpd.index)),rep("black",dim(data.new)[2]-length(ref_index)-length(cmpd.index))))
		legend("topright",title="Compound Labels",legend=c("Reference Compounds","Selected Compounds"),col=c("blue","red"),lty=1,bty="n")
		legend("topleft",title="Gene Labels",legend=rownames(data.new),lty=1,col=c(1:length(gene.index)),bty="n")
		plot.out(plot.type)
		
	}
	
	
	if(profile.type=="cmpd"){
		# Cmpd Profiles
		scores <- gene.scores
		# Re-order the genes first if thresholds have been set for the gene score.
		order.genes <- c(1:dim(data)[1])
		order.temp <- c()
		col.labels <- rep("black",dim(data)[1])
		
		if(!is.null(gene.thresP)){
			boolean_P <- (scores[,component.plot]>=gene.thresP)
			order.temp <- c(order.temp,order.genes[boolean_P])
			col.labels[boolean_P] <- thresP.col
		}
		else{
			boolean_P <- rep(FALSE,dim(data)[1])
		}
		if(!is.null(gene.thresN)){
			boolean_N <- (scores[,component.plot]<=gene.thresN)
			order.temp <- c(order.temp,order.genes[boolean_N])
			col.labels[boolean_N] <- thresN.col
		}
		else{
			boolean_N <- rep(FALSE,dim(data)[1])
		}
		
#		order.temp <- c(order.temp,order.genes[!boolean_P & !boolean_N])
		order.temp <- c(order.temp)
		col.labels <- col.labels[order.temp]
		
		data.new <- data[order.temp,c(ref_index,cmpd.index),drop=FALSE]
		lty.temp <- c(rep(1,length(ref_index)),rep(2,length(cmpd.index)))
		plot.in(plot.type,paste0(basefilename,"_CmpdProfiles.pdf"))
		plot(-1,0,type="n",main=paste0(main.base," - Compound Profiles for ",paste0(colnames(data.new),collapse=", ")),xaxt='n',xlab="Gene Index",ylab="Gene Expression",xlim=c(1,dim(data.new)[1]),ylim=c(min(data.new),max(data.new)))
		axis(1,at=c(1:dim(data.new)[1]),labels=FALSE)
		mtext(rownames(data.new),col=col.labels,at=c(1:dim(data.new)[1]),las=3,line=1,side=1)
		for(i.col in c(1:dim(data.new)[2])){
			lines(c(1:dim(data.new)[1]),data.new[,i.col],col=i.col,lty=lty.temp[i.col])
		}
		leg1.names <- colnames(data.new)
		leg1.names[c(1:length(ref_index))] <- paste0(leg1.names[c(1:length(ref_index))]," (Ref)")
		legend("topleft",leg1.names,col=c(1:dim(data.new)[2]),title="Compound Labels",lty=lty.temp,bty="n")
		legend("topright",c(paste0("Gene Score > ",gene.thresP),paste0("Gene Score < ",gene.thresN)),col=c(thresP.col,thresN.col),title="Gene Labels",lty=1,bty="n")
		
		plot.out(plot.type)
		
		
	}
	
	
	
}


plot_contributions <- function(CSresult,factor.plot,color.columns,legend.names=NULL,legend.cols=NULL,col.names,plot.type="device",basefilename="result"){
	
	if(class(CSresult)!="CSresult"){stop("Object is not of class CSresult",call.=FALSE)}
	if(CSresult@type!="CSmfa"){stop("Only possible for MFA results.",call.=FALSE)}
	
	if(!is.null(color.columns)){groupCol <- color.columns} else { groupCol <- "black"}
	
	## Plot-in and -out functions
	plot.in <- function(plot.type,name){
		if(plot.type=="pdf"){pdf(name)}
		if(plot.type=="device"){dev.new()}
		if(plot.type=="sweave"){}
	}
	plot.out <- function(plot.type){if(plot.type=="pdf"){dev.off()}}
	
	# note when taking contributions they are in the output in %
	ctr_loadings <- CSresult@object$quanti.var$contr[,factor.plot]/100
	ctr_scores <- CSresult@object$ind$contr[,factor.plot]/100
	
	# only loading contributions for now
	plot.in(plot.type,paste0(basefilename,"_MFALoadings_Ctr.pdf"))
	par(mfrow=c(1,1))
	plot(c(1:length(ctr_loadings)),ctr_loadings,  type="p",
			xlab="Compound Index", 
			ylab="Compound Contributions",
			pch=21,
			bg="grey",
			col=groupCol,
			cex=1,main=paste0("MFA ",factor.plot," - Compound Contributions")
	)
	text(c(1:length(ctr_loadings)),ctr_loadings, col.names,	pos=1,	cex=0.5,	col=groupCol)
#		if(length(legend.names)>0){legend((length(loadings[,factor.plot])*(1-0.45)),max(loadings),legend.names,pch=21,col=legend.cols,pt.bg="white",bty="n")}
	if(length(legend.names)>0){legend("topright",legend.names,pch=21,col=legend.cols,pt.bg="white",bty="n")}
	
	plot.out(plot.type)
}

#plot_contributions(out_MFA,factor.plot=1,color.columns=color.columns,legend.names=c("Ref","SP"),legend.cols=c("blue","red"),col.names=colnames(cbind(refMat,querMat)))

## some tests with simulated data (used in vignette)
#
#
#### VARIABLES ###
#
#head(out_MFA@object$quanti.var$coord)
##           Dim.1       Dim.2       Dim.3        Dim.4       Dim.5
## ref 1 0.8966614 -0.04849309 -0.09049056 -0.029259559 -0.01983110
## ref 2 0.8471277 -0.22926927  0.10776559 -0.004045711  0.01020555
## ref 3 0.8472827 -0.22255523  0.07060746  0.139170673  0.04289461
## ref 4 0.8450028 -0.21641341  0.14384994 -0.001313402 -0.02571653
## ref 5 0.8377648 -0.25332689  0.07001035 -0.041190796  0.02748658
## ref 6 0.8424596 -0.21580273  0.07969462  0.070208129  0.04355740
#head(out_MFA@object$quanti.var$contr)
##          Dim.1     Dim.2     Dim.3        Dim.4      Dim.5
## ref 1 9.250346 0.1721504 1.2646309 0.1390589296 0.06449031
## ref 2 8.256553 3.8480435 1.7935665 0.0026586000 0.01707945
## ref 3 8.259575 3.6259673 0.7699427 3.1460004259 0.30172172
## ref 4 8.215185 3.4285982 3.1957781 0.0002801941 0.10844909
## ref 5 8.075050 4.6979768 0.7569755 0.2755901509 0.12389176
## ref 6 8.165809 3.4092758 0.9808785 0.8006415435 0.31111783
#
#a = out_MFA@object$quanti.var$contr[,1] / out_MFA@object$quanti.var$coord[,1]^2
#unique(a)
## [1] 11.505381 11.505381  3.004126  3.004126  3.004126  3.004126
#a_w <- unique(a)[c(1,3)]
#a_w
## [1] 11.505381  3.004126
#
#alpha = unique(out_MFA@object$call$col.w)
#alpha
## [1] 0.21375952 0.05581393
#
##squared alpha is used for weighting if you want to do PCA on weighted matrix, otherwise alpha is used to weighting in generalized PCA
#sqrt(alpha)
## [1] 0.4623413 0.2362497
#
#a_w/alpha
## [1] 53.82395 53.82395
#
#alpha/a_w
## [1] 0.01857909 0.01857909
#
#out_MFA@object$eig[1,1]
## [1] 1.857909
#
## 'own' contributions
#
#ctr <- out_MFA@object$call$col.w * out_MFA@object$quanti.var$coord[,1]^2
#head(ctr)
##     ref 1     ref 2     ref 3     ref 4     ref 5     ref 6 
## 0.1718630 0.1533992 0.1534554 0.1526307 0.1500271 0.1517133 
#sum(ctr)
## [1] 1.857909
#
#
#### INDIVIDUALS ###
#
#
#head(out_MFA@object$ind$coord)
##            Dim.1     Dim.2       Dim.3        Dim.4       Dim.5
## gene-1 -2.166813 0.3094216  0.14659530  0.442442517 -0.01132318
## gene-2 -1.989373 0.3027026  0.61379498 -0.191495549  0.13904135
## gene-3 -2.200029 0.3307298  0.15092926  0.089350588 -0.37025398
## gene-4 -2.393666 0.3426540 -0.34154920 -0.005151102 -0.37498029
## gene-5 -2.282853 0.2280999  0.07620256 -0.350156120  0.09307208
## gene-6 -2.318218 0.1296621  0.23230141 -0.095681107 -0.56679636
#
#head(out_MFA@object$ind$contr)
##            Dim.1       Dim.2       Dim.3        Dim.4        Dim.5
## gene-1 0.2527076 0.032788689 0.015526441 1.487481e-01 9.835863e-05
## gene-2 0.2130139 0.031380150 0.272193947 2.786475e-02 1.483077e-02
## gene-3 0.2605148 0.037460147 0.016458061 6.066423e-03 1.051659e-01
## gene-4 0.3083918 0.040210023 0.084282682 2.016221e-05 1.078680e-01
## gene-5 0.2804992 0.017818572 0.004195377 9.316683e-02 6.645298e-03
## gene-6 0.2892571 0.005757701 0.038988407 6.956491e-03 2.464505e-01
#
#m = out_MFA@object$ind$contr[,1] / out_MFA@object$ind$coord[,1]^2
#unique(m)
## [1] 0.05382395 0.05382395 0.05382395 0.05382395 0.05382395
#m_w <- unique(m)[1]
#
#
#M = unique(out_MFA@object$call$row.w)
#M
## [1] 0.001
#
#M/m_w
## [1] 0.01857909
#

compare_CS_CSRankScores <- function(CSresult,color.columns=NULL,plot.type="device",basefilename=""){
	
	
	## Plot-in and -out functions
	plot.in <- function(plot.type,name){
		if(plot.type=="pdf"){pdf(name)}
		if(plot.type=="device"){dev.new()}
		if(plot.type=="sweave"){}
	}
	plot.out <- function(plot.type){if(plot.type=="pdf"){dev.off()}}	
	
	if(is.null(color.columns)){
		color.columns <- rep("black",dim(CSresult@CS$CS.query)[1])
	}
	else{
		color.columns <- color.columns[-c(1:dim(CSresult@CS$CS.ref)[1])]
		if(length(color.columns)!=dim(CSresult@CS$CS.query)[1]){stop("color.columns wrong lengths")}
	}
	
	if(CSresult@type=="CSfabia"){
		
		for(i.plot in 1:length(CSresult@CSRankScores)){
			loadings <- CSresult@CS$CS.query[,i.plot]
			scores <- CSresult@CSRankScores[[i.plot]]$CSRankScores
			
			plot.in(plot.type,basefilename)
			plot(loadings,scores,col=color.columns,bg="grey",pch=21,xlab="Loadings",ylab="CSRankScores",main=paste0(names(CSresult@CSRankScores)[i.plot],": Loadings VS CSRankScores"))
#			text(loadings,scores, rownames(CSresult@CS$CS.query),pos=1,cex=0.5,col=color.columns)
			plot.out(plot.type)
		}
		
	}else{
		loadings <- CSresult@CS$CS.query[,1]
		scores <- CSresult@CSRankScores[[1]]$CSRankScores
		
		plot.in(plot.type,basefilename)
		plot(loadings,scores,col=color.columns,bg="grey",pch=21,xlab="Loadings",ylab="CSRankScores",main=paste0(names(CSresult@CSRankScores)[1],": Loadings VS CSRankScores"))
#		text(loadings,scores, rownames(CSresult@CS$CS.query),pos=1,cex=0.5,col=color.columns)
		plot.out(plot.type)
		
	}
	
}


CSgrouploadings <- function(loadings,grouploadings.labels,grouploadings.cutoff,ref.index,method.name,component.name,basefilename,plot.type,plot=FALSE){
	
	## Plot-in and -out functions
	plot.in <- function(plot.type,name){
		if(plot.type=="pdf"){pdf(name,width=14)}
#		if(plot.type=="device"){dev.new()}
		if(plot.type=="sweave"){}
	}
	plot.out <- function(plot.type){if(plot.type=="pdf"){dev.off()}}	
	
		
	## Make Factor Labels 
	if(method.name=="FABIA"){component.short <- "BC"}else{component.short <- "F"}
	
	if(is.null(grouploadings.cutoff)){grouploadings.cutoff <- round(0.1*nrow(loadings))}
	
	
	labelmatrix <- matrix(FALSE,nrow=nrow(loadings),ncol=ncol(loadings))
	
	for(i in 1:ncol(labelmatrix)){
		labelmatrix[order(abs(loadings[,i]),decreasing=TRUE)[1:grouploadings.cutoff],i] <- TRUE
	}
	
	factor.labels <- apply(labelmatrix,MARGIN=1,FUN=function(x){
				if(sum(x)==0){return("None")}else{return(paste0(component.short,which(x),collapse="-"))}
			})
	
	# If other labels not provided; based on cutoff
		
	if(is.null(grouploadings.labels)){
		grouploadings.labels <- factor.labels
		unique.labels <- unique(factor.labels)
		unique.labels <- unique.labels[-which(unique.labels=="None")]
		unique.labels <- c(unique.labels,"None")
		
		# TO DO: still need some special ordering so None is as last
		# TO DO: RETURN LABELS TOO IN EXTRA! -> WANT TO USE THIS TO COMPARE WITH ACTUAL GROUPING -> MAKING FUNCTION FOR TO COMPARE THIS WITH ACTUAL LABELS -> FROM UNSUPERVISE TO SUPERVISE (1 plot for each label, then color with true labels)
	}else{
		unique.labels <- unique(grouploadings.labels)
	}
	
	
	## Make plots
	if(plot){
		
		inset.temp <- ifelse(plot.type=="pdf",0.1,0.25)
		
		# Correct colors
		colorpalette <- distinctColorPalette(length(unique.labels))
		names(colorpalette) <- unique.labels
		if(names(colorpalette)[length(colorpalette)]=="None"){colorpalette[length(colorpalette)] <- "grey"}
		
		col.temp <- rep("black",nrow(loadings))
		order.index <- rep(0,nrow(loadings))
		current.at <- 1
		
		for(i.label in unique.labels){
			label.index <- which(grouploadings.labels==i.label)
			order.index[current.at:(current.at+length(label.index)-1)] <- label.index
			current.at <- current.at+length(label.index)
			
			col.temp[label.index] <- colorpalette[i.label]
		}
		bg.temp <- col.temp
		col.temp[ref.index] <- "blue"

		# Sort the data + colors to get separate groups (in order of unique.labels)
		
		plot.in(plot.type,name=paste0(basefilename,"_grouploadings.pdf"))
		
		
		for(i in 1:ncol(loadings)){
			if(plot.type=="device"){dev.new()}
			par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
			
			plot(loadings[order.index,i],pch=21,col=col.temp[order.index],bg=bg.temp[order.index],xlab="Index",ylab="Loadings",main=paste0(method.name," - ",component.name," ",i," - Grouped Loadings"))
#			abline(h=0,lty=2)
			lines(c(0,nrow(loadings)+1),c(0,0),lty=2)
			
			legend("topright",c("Ref.",names(colorpalette)),col=c("blue",colorpalette),pt.bg=c("white",colorpalette),inset=c(-inset.temp,0),bty="n",pch=21)
		}
		
		
		plot.out(plot.type)
		
	}
	
	return(factor.labels)
}


get.loadings <- function(CSresult){
	type <- CSresult@type
	
	if(type=="CSfabia"){
		loadings <- CSresult@extra$object@L
	}
	else if(type=="CSmfa"){
		loadings <- CSresult@extra$object$quanti.var$cor
	}
	else if(type =="CSpca"){
		loadings <- CSresult@extra$object$var$cor
	}
	else if(type =="CSsmfa"){
		loadings <- CSresult@extra$object$loadings
	}
	else{
		stop("Result type not recognised")
	}
	return(loadings)	
	
}	




#' Compare Automatic Factor Labels with Manual Provided Labels.
#' 
#' With this function you can compare the automatic created labels based of the absolute loadings in \emph{CSanalysis} (\code{which=9}) with your own provided labels to investigate if there is relation between them.\cr
#' See the \code{type} parameter which two plots can be created.\cr
#' Note that the automatic created factor labels in \emph{CSanalysis} denote which factors this loading has a high/low value and these can be regenerated (with different a different cut-off) by simply running \emph{CSanalysis} again. Providing \code{resultavailable} will skip the analysis computation step and only regenerate the labels.
#'  
#' @export
#' @param CSresult Object of CSresult S4 Class.
#' @param labels Provide a vector with labels. (Length should be the number of references and queries together)
#' @param type 
#' \itemize{
#' \item \code{type="factorlabels"}:\cr
#' A K number of plots will be created (K = number of components in the analysis). Each plot will have the loadings on the y-axis and the original automatic generated factor labels on the x-axis.
#' The loadings are plotted for these factor labels (with jitter) and are colored according to the manual provided labels (\code{labels}) which is shown in the legend. The coloring also shows which loadings were in the reference set.
#'
#' \item \code{type="factors"}:\cr
#' A K number of plots will be created (K = number of components in the analysis). Each plot will have the loadings on the y-axis and factor labels on the x-axis.
#' These factor labels are not exactly the generated labels, but simply "Factor 1", "Factor 2",..., "None" or "BC 1", "BC 2",..., "None". This means that should a loading be high/low in multiple factors, it will appear multiple times on this plot, namely for each corresponding factor.
#' The loadings are plotted for these factor labels (with jitter) and are colored according to the manual provided labels (\code{labels}) which is shown in the legend. The coloring also shows which loadings were in the reference set.
#' }
#' Note that if none of the loadings is high/low in multiple factors, the types of plots should be identical.
#' @param basefilename Base of the filename when saving the graph as a pdf (\code{plot.type="pdf"})
#' @param plot.type How should the plots be outputted? \code{"pdf"} to save them in pdf files, \code{device} to draw them in a graphics device (default), \code{sweave} to use them in a sweave or knitr file.
#' @examples
#' \dontrun{
#' data("dataSIM",package="CSFA")
#' Mat1 <- dataSIM[,c(1:6)]
#' Mat2 <- dataSIM[,-c(1:6)]
#' 
#' MFA_out <- CSanalysis(Mat1,Mat2,"CSmfa",component.plot=1,which=c())
#' 
#' labels <- rep("Noise",ncol(dataSIM))
#' labels[c(1:31,332:341)] <- "Signal"
#' 
#' CSlabelscompare(CSresult=MFA_out,labels=labels,type="factors")
#' CSlabelscompare(CSresult=MFA_out,labels=labels,type="factorlabels")
#' }
CSlabelscompare <- function(CSresult,labels,type="factors",basefilename="CSanalysis",plot.type="device"){
	
	if(class(CSresult)!="CSresult"){stop("CSresult is not of the correct class type",call.=FALSE)}
	if(CSresult@type=="CSzhang"){stop("CSlabelscompare is not available for CSzhang results",call.=FALSE)}
	
	if(CSresult@type=="CSfabia"){
		component.name <- label.name <- "BC"
	}else{
		label.name <- "F"
		component.name <- "Factor"
	}
	
	factorlabels <- CSresult@extra$samplefactorlabels
	if(length(factorlabels)!=length(labels)){stop("Provided labels parameter does not have the correct length.",call.=FALSE)}
	
	loadings <- get.loadings(CSresult)
	
	ref.index <- 1:CSresult@call$dimensions$col[1]
			
			# jitter scatter plots (special function?)
	
	
	# make 2 plots? Based on factors and based on factorlabels?
	# Coloring always doen with true labels?
	# factors should also contain "other"
	
	## Plot-in and -out functions
	plot.in <- function(plot.type,name){
		if(plot.type=="pdf"){pdf(name,width=14)}
#		if(plot.type=="device"){dev.new()}
		if(plot.type=="sweave"){}
	}
	plot.out <- function(plot.type){if(plot.type=="pdf"){dev.off()}}	
	
	
	# General Preparation
	inset.temp <- ifelse(plot.type=="pdf",0.1,0.25)
	
	unique.truenames <- unique(labels)		
	colors.truenames <- distinctColorPalette(length(unique.truenames))
	bg.true <- colors.true <- sapply(labels,FUN=function(x){return(colors.truenames[which(x==unique.truenames)])})
	colors.true[ref.index] <- "blue"
	
	
	if(type=="factors"){
		
		
		# Get indices for new labels
		factor.list <- vector("list",ncol(loadings))
		for(i in 1:length(factor.list)){
			factor.list[[i]] <- which(sapply(factorlabels,FUN=function(x){return(grepl(paste0(label.name,i),x))}))
		}
		# Add None label
		factor.list[[length(factor.list)+1]] <- (1:nrow(loadings))[-unlist(factor.list[1:ncol(loadings)])]
		names(factor.list) <- c(paste0(label.name,1:ncol(loadings)),"None")
		
		xlim.temp <- c(1,length(factor.list))
		ylim.temp <- c(min(loadings),max(loadings))
		
		
		# there will be some duplicate stuff... should be taken into account for true labels
		plot.in(plot.type,paste0(basefilename,"_labelscompare_factors.pdf"))
		for(i in 1:ncol(loadings)){
			if(plot.type=="device"){dev.new()}
			par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
			
			plot(0,type="n",axes=FALSE,xlim=xlim.temp,ylim=ylim.temp,xlab="Factor Labels",ylab="Loadings",main=paste0(CSresult@type," - ",component.name," ",i))
		
			for(j in 1:length(factor.list)){
				data <- loadings[factor.list[[j]],i]
			
				if(length(data)==1){x.temp <- j}else{x.temp <- jitter(rep(j,length(data)),amount=0.2)}
			
				points(x.temp,data,pch=21,col=colors.true[factor.list[[j]]],bg=bg.true[factor.list[[j]]])
				lines(c(0.8,length(factor.list)+0.2),c(0,0),lty=2)
			}
			axis(1,at=1:length(factor.list),labels=names(factor.list))
			axis(2,at=seq(-1,1,0.2))
			legend("topright",c("Ref.",unique.truenames),col=c("blue",colors.truenames),pt.bg=c("white",colors.truenames),inset=c(-inset.temp,0),bty="n",pch=21)
		
		}
		plot.out(plot.type)
	}
	
	if(type=="factorlabels"){
	
		unique.factornames <- unique(factorlabels)
		if("None"%in%unique.factornames){
			unique.factornames <- unique.factornames[-which(unique.factornames=="None")]
			unique.factornames <- c(unique.factornames,"None")
		}
		
		xlim.temp <- c(1,length(unique.factornames))
		ylim.temp <- c(min(loadings),max(loadings))
		
		plot.in(plot.type,paste0(basefilename,"_labelscompare_factorlabels.pdf"))
		for(i in 1:ncol(loadings)){
			if(plot.type=="device"){dev.new()}
			par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
			
			plot(0,type="n",axes=FALSE,xlim=xlim.temp,ylim=ylim.temp,xlab="Factor Labels",ylab="Loadings",main=paste0(CSresult@type," - ",component.name," ",i))
			
			for(j in 1:length(unique.factornames)){
				current.select <- which(factorlabels==unique.factornames[j])
				data <- loadings[current.select,i]
				
				if(length(data)==1){x.temp <- j}else{x.temp <- jitter(rep(j,length(data)),amount=0.2)}
				
				points(x.temp,data,pch=21,col=colors.true[current.select],bg=bg.true[current.select])
				lines(c(0.8,length(unique.factornames)+0.2),c(0,0),lty=2)
			}
			axis(1,at=1:length(unique.factornames),labels=unique.factornames)
			axis(2,at=seq(-1,1,0.2))
			legend("topright",c("Ref.",unique.truenames),col=c("blue",colors.truenames),pt.bg=c("white",colors.truenames),inset=c(-inset.temp,0),bty="n",pch=21)
			
		}
		plot.out(plot.type)
	}
}


