# Project: CSFA
# 
# Author: lucp8394
###############################################################################

# NOTE: TO DO: Check if cor.rank also works when no plots are drawn


#' Compare CS Results.
#' 
#' After applying different CSanalysis on the same data, you can compare 2 different results of connectivity loadings, connectivity ranking scores and gene scores
#' Unless the result came from a Zhang and Gant analysis, you choose from which component (factor, PC, bicluster) the scores should be derived.
#' Further, for Zhang and Gant analysis, the "CRanking Scores" and "CLoadings" will be the same as the ZG Score as well as the p-values.
#' 
#' @export
#' @param CSresult1 First result.
#' @param CSresult2 Second result.
#' @param component1.plot If you are using a non-Zhang&Gant result, specify the bicluster, factor or principal component which should be used to derive connectivity scores from for the \emph{first} result.
#' @param component2.plot If you are using a non-Zhang&Gant result, specify the bicluster, factor or principal component which should be used to derive connectivity scores from for the \emph{second} result.
#' @param which Choose one or both plots which should be created.
#' \enumerate{
#' \item CS Comparison Plot
#' \item GS Comparison Plot
#' \item CSRankScores (Normal CS for CSzhang) Comparison Plot
#' \item CS p-values comparison plot (Raw & Adjusted).
#' \item CRankScores p-values comparison plot (Raw & Adjusted).
#' }
#' @param color.columns Vector of colors for the reference and query columns (compounds). If \code{NULL}, blue will be used for reference and black for query. Use this option to highlight reference columns and query columns of interest.
#' @param gene.thresP Vector of length 2 containing the positive gene thresholds for \code{CSresult1} and \code{CSresult2}. Genes above the threshold will be colored. (e.g. \code{c(1,2)})
#' @param gene.thresN Vector of length 2 containing the negative gene thresholds for \code{CSresult1} and \code{CSresult2}. Genes below the threshold will be colored. (e.g. \code{c(-1,-2)})
#' @param thresP.col Vector of length 2 containing the colors for the high gene scores for \code{CSresult1} and \code{CSresult2} (e.g. \code{c("blue","light blue")}).
#' @param thresN.col Vector of length 2 containing the colors for the low gene scores for \code{CSresult1} and \code{CSresult2} (e.g. \code{c("red","pink")}).
#' @param legend.names Option to draw a legend (about the highlights in \code{color.columns}) in the CS plot. If \code{NULL}, only references are in the legend.
#' @param legend.cols Colors to be used for the \code{legend.names}.
#' @param legend.pos The location of the legend: \code{"bottomright"}, \code{"bottom"}, \code{"bottomleft"}, \code{"left"}, \code{"topleft"}, \code{"top"}, \code{"topright"}, \code{"right"} and \code{"center"}.
#' @param plot.type How should the plots be outputted? \code{"pdf"} to save them in pdf files, \code{device} to draw them in a graphics device (default), \code{sweave} to use them in a sweave or knitr file.
#' @param basefilename Basename of the graphs if saved in pdf files
#' @param threshold.pvalues If both CSresult1 and CSresult contain pvalues (and adjusted pvalues), this threshold will be used to compare the number of overlapping significant results. 
#' @return A list object with 2 slotes. In the first slot, Pearson and Spearman correlation between the results (CLoadings, Gene Scores, CRanking Scores, (adjusted) p-values) can be found. The second slot, if permutation was applied, contaisn a small comparison between the significant results based on \code{threshold.pvalues}.
#' @examples
#' \dontrun{
#' data("dataSIM",package="CSFA")
#' Mat1 <- dataSIM[,c(1:6)]
#' Mat2 <- dataSIM[,-c(1:6)]
#' 
#' MFA_analysis <- CSanalysis(Mat1,Mat2,"CSmfa")
#' ZHANG_analysis <- CSanalysis(Mat1,Mat2,"CSzhang")
#' 
#' CScompare(MFA_analysis,ZHANG_analysis,1)
#' }
CScompare <- function(CSresult1,CSresult2,component1.plot,component2.plot,threshold.pvalues=0.05,which=c(1,2,3),color.columns=NULL,gene.thresP=NULL,gene.thresN=NULL,thresP.col=c("blue","light blue"),thresN.col=c("red","pink"),legend.names=NULL,legend.cols=NULL,legend.pos="topright",plot.type="device",basefilename="CScompare"){
	

	if(class(CSresult1)!="CSresult"){stop("CSresult1 is not of the correct class type")}
	if(class(CSresult2)!="CSresult"){stop("CSresult2 is not of the correct class type")}
	
	
	refdim1 <- CSresult1@call$dimensions$col[1]
	refdim2 <- CSresult2@call$dimensions$col[1]
	querdim1 <- CSresult1@call$dimensions$col[2]
	querdim2 <- CSresult2@call$dimensions$col[2]
	
	
	if(refdim1!=refdim2){stop("Using 2 different reference matrices",call.=FALSE)}
	if(querdim1!=querdim2){stop("Using 2 different query matrices",call.=FALSE)}
	if(CSresult1@call$dimensions$row != CSresult2@call$dimensions$row){warning("Different amount of genes between 2 results. No correlation computation or scatter plot will be done for the Gene Scores.",call.=FALSE)}
	if(!all(rownames(CSresult1@CS[[1]]$CS.query)==rownames(CSresult2@CS[[1]]$CS.query))){stop("Different rownames for 2 results",call.=FALSE)}
	
	CSGS1 <- get.CS.GS(CSresult1,component1.plot,refdim1)
	CSGS2 <- get.CS.GS(CSresult2,component2.plot,refdim2)
	
	loadings1 <- CSGS1$CS
	loadings2 <- CSGS2$CS
	scores1 <- CSGS1$GS
	scores2 <- CSGS2$GS
	name1 <- CSGS1$name
	name2 <- CSGS2$name
	axename1 <- CSGS1$axename
	axename2 <- CSGS2$axename
	rankscores1 <- CSGS1$CSRank
	rankscores2 <- CSGS2$CSRank
	names(rankscores1) <- names(rankscores2) <- names(loadings1)[-c(1:refdim1)]
	
	
	
	if(CSresult1@type!="CSzhang" & CSresult2@type!="CSzhang"){
		if(is.null(legend.names) & is.null(legend.cols) &is.null(color.columns)){
			
			color.columns <- c(rep("blue",refdim1),rep("black",querdim1))
			legend.names <- c("References")
			legend.cols <- c("blue")
			
			legend.names.rank <- c()
			legend.cols.rank <- c()
		}else{
			if(is.null(legend.names)){legend.names <- c()}
			if(is.null(legend.cols)){legend.cols <- c()}
			if(is.null(color.columns)){color.columns <- rep("black",querdim1+refdim1)}
						
			legend.names.rank <- legend.names
			legend.cols.rank <- legend.cols
		}
	
	}else{
		if(is.null(legend.names)){legend.names <- c()}
		if(is.null(legend.cols)){legend.cols <- c()}
		if(is.null(color.columns)){color.columns <- rep("black",querdim1+refdim1)}
		
		legend.names.rank <- legend.names
		legend.cols.rank <- legend.cols
		
	}
	
	color.columns.loadings <- color.columns
	
	# add a check here if one of the results is Zhang, if so (unless 2 zhang), delete ref.index + also use correct default colors (no refblue)
	if(length(loadings1)!=length(loadings2)){
		if(length(loadings1)<length(loadings2)){
			loadings2 <- loadings2[-c(1:refdim2)]
		}else{
			loadings1 <- loadings1[-c(1:refdim1)]
		}
		
		if(!is.null(color.columns)){color.columns.loadings <- color.columns[-c(1:refdim1)]}
	}
	
	if(length(legend.names)!=length(legend.cols)){stop("Error in legend parameters, different length of legend names and colors.",call.=FALSE)}
	
	
	
		# below needs to be done correctly (e.g. correct color columns length)
	
	scores_correlation <- matrix(NA,nrow=2,ncol=3,dimnames=list(c("Correlation_Pearson","Correlation_Spearman"),c("CLoadings","CRankScores","GeneScores")))
	
	
	if(1%in%which){
		compare.CS.plot(loadings1=loadings1,loadings2=loadings2,name1=name1,name2=name2,axename1=axename1,axename2=axename2,nref=refdim1,color.columns=color.columns.loadings,legend.names=legend.names,legend.cols=legend.cols,legend.pos=legend.pos,plot.type=plot.type,basefilename=basefilename)
	}
	scores_correlation[1,1] <- cor(loadings1,loadings2,use="complete.obs")
	scores_correlation[2,1] <- cor(loadings1,loadings2,use="complete.obs",method="spearman")
	
	
	if(!(is.null(scores1)|is.null(scores2)|!(2 %in% which))){
		if(CSresult1@call$dimensions$row == CSresult2@call$dimensions$row){
				compare.GS.plot(scores1=scores1,scores2=scores2,name1=name1,name2=name2,axename1=axename1,axename2=axename2,nref=refdim1,gene.thresP=gene.thresP,gene.thresN=gene.thresN,thresP.col=thresP.col,thresN.col=thresN.col,plot.type=plot.type,basefilename=basefilename)
		}
	}
	if(!(is.null(scores1)|is.null(scores2))){
		if(CSresult1@call$dimensions$row == CSresult2@call$dimensions$row){
			scores_correlation[1,3] <- cor(scores1,scores2,use="complete.obs")
			scores_correlation[2,3] <- cor(scores1,scores2,use="complete.obs",method="spearman")
		}else{
			scores_correlation[1,3] <- NA
			scores_correlation[2,3] <- NA
		}
	}
	
	
	
	if(3%in%which){
		if(length(color.columns)==(refdim1+querdim1)){color.columns.rank <- color.columns[-c(1:refdim1)]}else{color.columns.rank <-color.columns}
		
		compare.CSRank.plot(rankscores1=rankscores1,rankscores2=rankscores2,name1=name1,name2=name2,axename1=axename1,axename2=axename2,color.columns=color.columns.rank,legend.names=legend.names.rank,legend.cols=legend.cols.rank,legend.pos=legend.pos,plot.type=plot.type,basefilename=basefilename)
	}
	scores_correlation[1,2] <- cor(rankscores1,rankscores2,use="complete.obs")
	scores_correlation[2,2] <- cor(rankscores1,rankscores2,use="complete.obs",method="spearman")
	
	out_pval_compare <- NULL
	cor.pvalues <- vector("list",2)
		
	if(!is.null(CSresult1@permutation.object) & !is.null(CSresult2@permutation.object)){
		
		# NEED TO ADD A CHECK THAT CHOSEN FACTOR IS SAME AS ONE IN PERMUTATION OJECT
		correct.pvalues <- TRUE
		if(CSresult1@type!="CSzhang" & correct.pvalues){
			correct.pvalues <- ifelse(CSresult1@permutation.object$extra.parameter$mfa.factor==component1.plot,TRUE,FALSE)
		}
		if(CSresult2@type!="CSzhang" & correct.pvalues){
			correct.pvalues <- ifelse(CSresult2@permutation.object$extra.parameter$mfa.factor==component2.plot,TRUE,FALSE)
		}
	
		if(correct.pvalues){
			if(CSresult1@type=="CSzhang"){CSresult1@permutation.object$CSRank.pval.dataframe <- CSresult1@permutation.object$CS.pval.dataframe}
			if(CSresult2@type=="CSzhang"){CSresult2@permutation.object$CSRank.pval.dataframe <- CSresult2@permutation.object$CS.pval.dataframe}
						
			
			list.pval.dataframe <- list(CSresult1@permutation.object$CS.pval.dataframe,CSresult2@permutation.object$CS.pval.dataframe)
			list.pval.dataframe.rank <- list(CSresult1@permutation.object$CSRank.pval.dataframe,CSresult2@permutation.object$CSRank.pval.dataframe)
			
			
			
			#p-values compare plot (2 plots) ### ADD A WAY TO GIVE c(TRUE TRUE) to plot to give either CS or CSRank....  + check if adjusted inside plot=TRUE
			if((4 %in% which) | (5 %in% which)){			
				plot <- c((4 %in% which) , (5 %in% which))
				cor.pvalues <- compare.pvalues.plot(list.pval.dataframe,list.pval.dataframe.rank,nref=refdim1,name1=name1,name2=name2,axename1=axename1,axename2=axename2,color.columns=color.columns[-c(1:refdim1)],legend.names=legend.names,legend.cols=legend.cols,legend.pos=legend.pos,plot.type=plot.type,basefilename=paste0(basefilename,"_CS"),plot=plot)
			}
			else{
				cor.pvalues <- compare.pvalues.plot(list.pval.dataframe,list.pval.dataframe.rank,nref=refdim1,name1=name1,name2=name2,axename1=axename1,axename2=axename2,color.columns=color.columns[-c(1:refdim1)],legend.names=legend.names,legend.cols=legend.cols,legend.pos=legend.pos,plot.type=plot.type,basefilename=paste0(basefilename,"_CS"),plot=c(FALSE,FALSE))
			}
			
			out_pval_compare <- pvalue2_compare(list.pval.dataframe,threshold=threshold.pvalues)
			
		}else{
			warning("P-values are available, but not for the chosen component.plot.")
		}
		

	}

	out.correlation <- list(scores=scores_correlation,pvalues=cor.pvalues[[1]],adj.pvalues=cor.pvalues[[2]])
	
	return(list(correlation=out.correlation,comparison=out_pval_compare))
}



get.CS.GS <- function(CSresult,component.plot,refdim){
	type <- CSresult@type
	
	if(type=="CSfabia"){
		loadings <- CSresult@extra$object@L[,component.plot]
		scores <- t(CSresult@extra$object@Z)[,component.plot]
		rankscores <- CSrank2(CSresult@extra$object@L,1:refdim,plot=FALSE,component.plot=component.plot)[,1]
		
		
		return(list(CS=loadings,GS=scores,name="FABIA",axename=paste0("Fabia BC ",component.plot),CSRank=rankscores))
			
	}
	else if(type=="CSmfa"){
		loadings <- CSresult@extra$object$quanti.var$cor[,component.plot]		
		scores <- CSresult@extra$object$ind$coord[,component.plot]	
		
		rankscores <- CSrank2(CSresult@extra$object$quanti.var$cor,1:refdim,plot=FALSE,component.plot=component.plot)[,1]
		
		return(list(CS=loadings,GS=scores,name="MFA",axename=paste0("MFA Factor ",component.plot),CSRank=rankscores))
	}
	else if(type =="CSpca"){
		loadings <- CSresult@extra$object$var$cor[,component.plot]
		scores <- CSresult@extra$object$ind$coord[,component.plot]
		rankscores <- CSrank2(CSresult@extra$object$var$cor,1:refdim,plot=FALSE,component.plot=component.plot)[,1]
		
		
		return(list(CS=loadings,GS=scores,name="PCA",axename=paste0("PCA PC ",component.plot),CSRank=rankscores))
	}
	else if(type =="CSsmfa"){
		loadings <- CSresult@extra$object$loadings[,component.plot]
		scores <- CSresult@extra$object$scores[,component.plot]
		rankscores <- CSrank2(CSresult@extra$object$loadings,1:refdim,plot=FALSE,component.plot=component.plot)[,1]
		
		
		return(list(CS=loadings,GS=scores,name="sMFA",axename=paste0("sMFA Factor ",component.plot),CSRank=rankscores))
	}
	else if(type == "CSzhang"){
		loadings <- CSresult@CS$CS.query[,1]
		names(loadings) <- rownames(CSresult@CS$CS.query)
		
		
		return(list(CS=loadings,GS=NULL,name="ZHANG",axename="Zhang CS",CSRank=loadings))
	}
	else{
		stop("Result type not recognised")
	}
}
	


compare.CS.plot <- function(loadings1,loadings2,name1,name2,axename1,axename2,nref,color.columns,legend.names,legend.cols,legend.pos,plot.type,basefilename){
	
	## Plot-in and -out functions
	plot.in <- function(plot.type,name){
		if(plot.type=="pdf"){pdf(name)}
		if(plot.type=="device"){dev.new()}
		if(plot.type=="sweave"){}
	}
	plot.out <- function(plot.type){if(plot.type=="pdf"){dev.off()}}	
	
	##
	if(!is.null(color.columns)){groupCol <- color.columns} else { groupCol <- "black"}
	##	
	
	
	minX <- min(-1,loadings1,na.rm=TRUE)
	maxX <- max(1,loadings1,na.rm=TRUE)
	minY <- min(-1,loadings2,na.rm=TRUE)
	maxY <- max(1,loadings2,na.rm=TRUE)
	
	plot.in(plot.type,paste0(basefilename,"_CS_",name1,"_VS_",name2,".pdf"))
	par(mfrow=c(1,1))
	plot(loadings1,loadings2,xlim=c(minX,maxX),ylim=c(minY,maxY),col=groupCol,bg="grey",main=paste0(name1," VS ",name2," CScores - ",nref," Ref Compound"),xlab=paste0(axename1," Connectivity Scores"),ylab=paste0(axename2," Connectivity Scores"),pch=21)
	text(loadings1,loadings2, names(loadings1),	pos=1,	cex=0.5,	col=groupCol)
	if(length(legend.names)>0){legend(legend.pos,legend.names,pch=21,col=legend.cols,bty="n")}
	plot.out(plot.type)

}



compare.GS.plot <- function(scores1,scores2,name1,name2,axename1,axename2,nref,gene.thresP,gene.thresN,thresP.col,thresN.col,plot.type,basefilename){
	
	## Plot-in and -out functions
	plot.in <- function(plot.type,name){
		if(plot.type=="pdf"){pdf(name)}
		if(plot.type=="device"){dev.new()}
		if(plot.type=="sweave"){}
	}
	plot.out <- function(plot.type){if(plot.type=="pdf"){dev.off()}}	
	
	## Gene Scores Coloring function
	.give.gene.color <- function(data,P1,N1,P2,N2,col.P1,col.N1,col.P2,col.N2){
		data1 <- data[1]
		data2 <- data[2]
		
		if(data2>=P2){
			if(data1<=N1){return(c(col.P2,col.N1))}
			else if(data1>=P1){return(c(col.P2,col.P1))}
			else if(data1>N1 & data1<P1){return(c(col.P2,col.P2))}
			else {return(c("grey","grey"))}
		}
		else if(data2<=N2){
			if(data1<=N1){return(c(col.N2,col.N1))}
			else if(data1>=P1){return(c(col.N2,col.P1))}
			else if(data1>N1 & data1<P1){return(c(col.N2,col.N2))}
			else {return(c("grey","grey"))}
		}
		else if(data1<=N1){return(c(col.N1,col.N1))}
		else if(data1>=P1){return(c(col.P1,col.P1))}
		else {return(c("grey","grey"))}
	}
		
	minX <- min(-1,scores1,na.rm=TRUE)
	maxX <- max(1,scores1,na.rm=TRUE)
	minY <- min(-1,scores2,na.rm=TRUE)
	maxY <- max(1,scores2,na.rm=TRUE)
	
	# gene colors
	if(!is.null(gene.thresP) | !is.null(gene.thresN)){
		if(is.null(gene.thresP)){gene.thresP <- c(99999,99999)}
		if(is.null(gene.thresN)){gene.thresN <- c(-99999,-99999)}
		
		scores <- rbind(scores1,scores2)
		list.scores <- as.list(as.data.frame(scores))
		list.colors <- lapply(X=list.scores,FUN=.give.gene.color,P1=gene.thresP[1],P2=gene.thresP[2],N1=gene.thresN[1],N2=gene.thresN[2],col.P1=thresP.col[1],col.P2=thresP.col[2],col.N1=thresN.col[1],col.N2=thresN.col[2])
		list.colors <- t(as.data.frame(list.colors))
		colnames(list.colors) <- c("bg","col")
	}
	else{
		list.colors <- matrix("grey",ncol=2,nrow=length(scores2)) 
	}
	
	plot.in(plot.type,paste0(basefilename,"_GS_",name1,"_VS_",name2,".pdf"))
	par(mfrow=c(1,1))
	plot(scores1,scores2,xlim=c(minX,maxX),ylim=c(minY,maxY),col=list.colors[,2],bg=list.colors[,1],main=paste0(name1," VS ",name2," GScores - ",nref," Ref Compound"),xlab=paste0(axename1," Gene Scores"),ylab=paste0(axename2," Gene Scores"),pch=21)
	text(scores1,scores2, names(scores1),	pos=1,	cex=0.5,	col=list.colors[,2])
	if(!is.null(gene.thresP)){
		abline(v=gene.thresP[1],lty=3)
		abline(h=gene.thresP[2],lty=3)
	}
	if(!is.null(gene.thresN)){
		abline(v=gene.thresN[1],lty=3)
		abline(h=gene.thresN[2],lty=3)
	}
	plot.out(plot.type)

}

compare.CSRank.plot <- function(rankscores1,rankscores2,name1,name2,axename1,axename2,color.columns,legend.names,legend.cols,legend.pos,plot.type,basefilename){
	
	## Plot-in and -out functions
	plot.in <- function(plot.type,name){
		if(plot.type=="pdf"){pdf(name)}
		if(plot.type=="device"){dev.new()}
		if(plot.type=="sweave"){}
	}
	plot.out <- function(plot.type){if(plot.type=="pdf"){dev.off()}}	
	
	##
	if(!is.null(color.columns)){groupCol <- color.columns} else { groupCol <- "black"}
	##	
	
	minX <- min(-1,rankscores1,na.rm=TRUE)
	maxX <- max(1,rankscores1,na.rm=TRUE)
	minY <- min(-1,rankscores2,na.rm=TRUE)
	maxY <- max(1,rankscores2,na.rm=TRUE)
	
	
	axename1.temp <- ifelse(axename1=="Zhang CS","Zhang & Gant Score",paste(axename1," Connectivity Rank Score"))
	axename2.temp <- ifelse(axename2=="Zhang CS","Zhang & Gant Score",paste(axename2," Connectivity Rank Score"))
	
		
	plot.in(plot.type,paste0(basefilename,"_CSRank_",name1,"_VS_",name2,".pdf"))
	par(mfrow=c(1,1))
	plot(rankscores1,rankscores2,xlim=c(minX,maxX),ylim=c(minY,maxY),col=groupCol,bg="grey",main=paste0(name1," VS ",name2," CRankScores"),xlab=axename1.temp,ylab=axename2.temp,pch=21)
	text(rankscores1,rankscores2, names(rankscores1),	pos=1,	cex=0.5,	col=groupCol)
	if(length(legend.names)>0){legend(legend.pos,legend.names,pch=21,col=legend.cols,bty="n")}
	plot.out(plot.type)
	
	
}


compare.pvalues.plot <- function(list.pval.dataframe,list.pval.dataframe.rank,nref,name1,name2,axename1,axename2,color.columns,legend.names,legend.cols,legend.pos,plot.type="device",basefilename="",plot=TRUE){
	## Plot-in and -out functions
	plot.in <- function(plot.type,name){
		if(plot.type=="pdf"){pdf(name)}
		if(plot.type=="device"){dev.new()}
		if(plot.type=="sweave"){}
	}
	plot.out <- function(plot.type){if(plot.type=="pdf"){dev.off()}}	
	
	pvalues_correlation <- matrix(NA,nrow=2,ncol=2,dimnames=list(c("Correlation_Pearson","Correlation_Spearman"),c("CLoadings","CRankScores")))
	adjustedpvalues_correlation <- matrix(NA,nrow=2,ncol=2,dimnames=list(c("Correlation_Pearson","Correlation_Spearman"),c("CLoadings","CRankScores")))
	
	list.temp <- list(list.pval.dataframe,list.pval.dataframe.rank)
	
	rank <- FALSE
	for(i in 1:length(list.temp)){
		list.pval.dataframe <- list.temp[[i]]
		if(!is.null(list.pval.dataframe[[1]]) & !is.null(list.pval.dataframe[[2]])){
			if(i==2){rank <- TRUE}
			
			adjusted.available <- unlist(lapply(list.pval.dataframe,FUN=function(x){"pvalues.adjusted"%in%colnames(x)}))
			
			if(sum(adjusted.available)==2){
				use.adjust <- TRUE
			}else{
				use.adjust <- FALSE
			}
			
			
			pvalues1 <- -log(list.pval.dataframe[[1]]$pvalues)
			pvalues2 <- -log(list.pval.dataframe[[2]]$pvalues)
			pvalues.name <- "p-values"
			pvalues_correlation[1,i] <- cor(pvalues1,pvalues2,use="complete.obs")
			pvalues_correlation[2,i] <- cor(pvalues1,pvalues2,use="complete.obs",method="spearman")
			
			if(plot[i]){
				
				rank.name <- ifelse(rank,"CRank","")
				
				axename1.temp <- ifelse(axename1=="Zhang CS","Zhang & Gant Score",paste(axename1," Connectivity ",rank.name," Score"))
				axename2.temp <- ifelse(axename2=="Zhang CS","Zhang & Gant Score",paste(axename2," Connectivity ",rank.name," Score"))
				
				plot.in(plot.type,paste0(basefilename,"_",paste0(rank.name,"CScore",name1,"VS",name2,"_pvalues.pdf")))
				plot(pvalues1,pvalues2,col=color.columns,bg="grey",main=paste0(name1," VS ",name2," -log(",pvalues.name,") of CS",rank.name),xlab=axename1.temp,ylab=axename2.temp,pch=21)
				text(pvalues1,pvalues2,as.character(list.pval.dataframe[[1]][,1]),col=color.columns,pos=1,cex=0.5)
				if(length(legend.names)>0){legend(legend.pos,legend.names,pch=21,col=legend.cols,bty="n")}
				plot.out(plot.type)
				
#				plot.in(plot.type,paste0(basefilename,"_",paste0(name1,"VS",name2,"_rankpvalues.pdf")))
#				plot(rank(pvalues1),rank(pvalues2),col=color.columns,bg="grey",main=paste0(name1," VS ",name2," ranks of -log(",pvalues.name,") of CS",rank.name),xlab=axename1.temp,ylab=axename2.temp,pch=21)
#				text(rank(pvalues1),rank(pvalues2),as.character(list.pval.dataframe[[1]][,1]),col=color.columns,pos=1,cex=0.5)
#				if(length(legend.names)>0){legend(legend.pos,legend.names,pch=21,col=legend.cols,bty="n")}
#				plot.out(plot.type)
			
			}
			
			if(use.adjust){
				pvalues1 <- -log(list.pval.dataframe[[1]]$pvalues.adjusted)
				pvalues2 <- -log(list.pval.dataframe[[2]]$pvalues.adjusted)			
				adjustedpvalues_correlation[1,i] <- cor(pvalues1,pvalues2,use="complete.obs")
				adjustedpvalues_correlation[2,i] <- cor(pvalues1,pvalues2,use="complete.obs",method="spearman")
				
				pvalues.name <- "adjusted p-values"
						
			
				if(plot[i]){
					rank.name <- ifelse(rank,"CRank","")
				
					axename1.temp <- ifelse(axename1=="Zhang CS","Zhang & Gant Score",paste(axename1," Connectivity ",rank.name," Score"))
					axename2.temp <- ifelse(axename2=="Zhang CS","Zhang & Gant Score",paste(axename2," Connectivity ",rank.name," Score"))
				
					plot.in(plot.type,paste0(basefilename,"_",paste0(rank.name,"CScore",name1,"VS",name2,"_pvalues.pdf")))
					plot(pvalues1,pvalues2,col=color.columns,bg="grey",main=paste0(name1," VS ",name2," -log(",pvalues.name,") of CS",rank.name),xlab=axename1.temp,ylab=axename2.temp,pch=21)
					text(pvalues1,pvalues2,as.character(list.pval.dataframe[[1]][,1]),col=color.columns,pos=1,cex=0.5)
					if(length(legend.names)>0){legend(legend.pos,legend.names,pch=21,col=legend.cols,bty="n")}
					plot.out(plot.type)
				}
			}
		}
	}
	return(list(pval=pvalues_correlation,adjpval=adjustedpvalues_correlation))
}




