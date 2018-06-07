# Project: Connectivity
# 
# Author: lucp8394
###############################################################################

## EXAMPLE DATA ##
#' Simulated Microarray Data
#'
#' A matrix containing some simulated example microarray data. 
#' The first 6 columns of this matrix make up the reference matrix part.
#'
#' @format A matrix with 1000 rows and 341 columns.
#' @name dataSIM
NULL


## IMPORTS ##

#' @import methods stats graphics grDevices
#' @importFrom fabia fabia extractBic showSelected
#' @importFrom FactoMineR PCA MFA
#' @importFrom pls stdize
#' @importFrom elasticnet spca arrayspc
#' @importFrom randomcoloR distinctColorPalette
#' @importFrom parallel detectCores splitIndices clusterCall clusterExport
#' @importFrom snowFT clusterApplyFT makeClusterFT stopClusterFT clusterSetupRNG.FT

## CLASSES ##
setClass("CSfabia",slots=c(call="call"))
setClass("CSmfa",slots=c(call="call"))
setClass("CSpca",slots=c(call="call"))
setClass("CSsmfa",slots=c(call="call"))
setClass("CSzhang",slots=c(call="call"))

#' An S4 class in which the results of the Connectivity Scores by Factor Analysis are stored.
#' 
#' @export
#' @slot type A character string containing the analysis type.
#' @slot CS List of any number of lists (depending on how many components were selected) which contain the connectivity loadings and ranking scores for the query (and reference loadings). If permutation was applied, will also contain p-values.
#' @slot GS Dataframe containing the gene scores.
#' @slot extra List which contains \code{CSRank_Full} (contains all intermediate values while calculating the CS Ranking Score), \code{Object} (contains the complete original FA or Zhang result) and \code{samplefactorlabels} (contains thresholded labels based on the factor loadings, see plot \code{which=8}).
#' @slot permutation.object Contains CS for permuted data (matrix) and a dataframe with the p-values (only for MFA and Zhang).
#' @slot call List object containing the original call of \code{CSanalysis} as well as the parameters for the chosen method.
setClass("CSresult",slots=list(type="character",CS="list",GS="data.frame",extra="list",permutation.object="ANY",call="ANY"))


# CHECK IF THESE STILL NECESSARY
#setClass("CSzhangCompare",slots=c(CSresult1="CSresult",CSresult2="CSresult"))
#setClass("CSfabiaCompare",slots=c(CSresult1="CSresult",CSresult2="CSresult"))


## METHODS  - CSANALYSIS ##
#' Connectivity Score Analysis.
#' 
#' Doing a CS analysis, interactively generating graphs. See specific type for additional parameteres.\cr
#' Types:\cr
#' \itemize{
#' \item \code{\link[=CSanalysis,matrix,matrix,CSzhang-method]{Zhang and Gant}}
#' \item \code{\link[=CSanalysis,matrix,matrix,CSmfa-method]{MFA}}
#' \item \code{\link[=CSanalysis,matrix,matrix,CSpca-method]{PCA}}
#' \item \code{\link[=CSanalysis,matrix,matrix,CSsmfa-method]{Sparse MFA}}
#' \item \code{\link[=CSanalysis,matrix,matrix,CSfabia-method]{FABIA}}
#' }
#' 
#' @export
#' @param refMat Reference matrix (Rows = genes and columns = compounds)
#' @param querMat Query matrix
#' @param type Type of Factor Analysis or Zhang & Gant ( \code{"CSfabia"}, \code{"CSmfa"}, \code{"CSpca"}, \code{"CSsmfa"} or \code{"CSzhang"})
#' @param ... Additional parameters for analysis
#' @return An object of the S4 Class \code{\link{CSresult-class}}.
#' @examples
#' \dontshow{
#' data("dataSIM",package="CSFA")
#' Mat1 <- dataSIM[,c(1:6)]
#' Mat2 <- dataSIM[,-c(1:6)]
#' ZHANG_analysis <- CSanalysis(Mat1,Mat2,"CSzhang")
#' }
#'  
#' \donttest{
#' data("dataSIM",package="CSFA")
#' Mat1 <- dataSIM[,c(1:6)]
#' Mat2 <- dataSIM[,-c(1:6)]
#' 
#' MFA_analysis <- CSanalysis(Mat1,Mat2,"CSmfa")
#' FABIA_analysis <- CSanalysis(Mat1,Mat2,"CSfabia")
#' ZHANG_analysis <- CSanalysis(Mat1,Mat2,"CSzhang")
#' }
setGeneric('CSanalysis', function(refMat,querMat,type, ...){standardGeneric('CSanalysis')})


#' Connectivity Score Analysis.
#' 
#' Doing a CS analysis, interactively generating graphs. See specific type for additional parameteres.\cr
#' Types:\cr
#' \itemize{
#' \item \code{\link[=CSanalysis,matrix,matrix,CSzhang-method]{Zhang and Gant}}
#' \item \code{\link[=CSanalysis,matrix,matrix,CSmfa-method]{MFA}}
#' \item \code{\link[=CSanalysis,matrix,matrix,CSpca-method]{PCA}}
#' \item \code{\link[=CSanalysis,matrix,matrix,CSsmfa-method]{Sparse MFA}}
#' \item \code{\link[=CSanalysis,matrix,matrix,CSfabia-method]{FABIA}}
#' }
#' 
#' @export
#' @param refMat Reference matrix (Rows = genes and columns = compounds)
#' @param querMat Query matrix
#' @param type Type of Factor Analysis or Zhang & Gant ( \code{"CSfabia"}, \code{"CSmfa"}, \code{"CSpca"}, \code{"CSsmfa"} or \code{"CSzhang"})
#' @param ... Additional parameters for analysis
#' @return An object of the S4 Class \code{\link{CSresult-class}}.
#' @examples
#' \dontshow{
#' data("dataSIM",package="CSFA")
#' Mat1 <- dataSIM[,c(1:6)]
#' Mat2 <- dataSIM[,-c(1:6)]
#' ZHANG_analysis <- CSanalysis(Mat1,Mat2,"CSzhang")
#' }
#' \donttest{
#' data("dataSIM",package="CSFA")
#' Mat1 <- dataSIM[,c(1:6)]
#' Mat2 <- dataSIM[,-c(1:6)]
#' 
#' MFA_analysis <- CSanalysis(Mat1,Mat2,"CSmfa")
#' FABIA_analysis <- CSanalysis(Mat1,Mat2,"CSfabia")
#' ZHANG_analysis <- CSanalysis(Mat1,Mat2,"CSzhang")
#' }
setMethod('CSanalysis', c('matrix','matrix','character'),
		function(refMat,querMat,type, ...) {
			if(type %in% c("CSfabia","CSmfa","CSpca","CSsmfa","CSzhang")){
				type <- new(type,call=match.call())
				CSanalysis(refMat,querMat,type,...)
			}
			else{
				stop("This method type is not available.")
			}
		})


# FABIA

#' "CSfabia"
#' 
#' Doing interactive CS analysis with FABIA (Factor Analysis for Bicluster Acquisition). One or multiple reference compounds are possible in this analysis.
#' 
#' @export 
#' @param refMat Reference matrix (Rows = genes and columns = compounds)
#' @param querMat Query matrix
#' @param type \code{"CSfabia"}
#' @param p \emph{Fabia Parameter:} number of hidden factors = number of biclusters; default = 13
#' @param alpha \emph{Fabia Parameter:} sparseness loadings (0 - 1.0); default = 0.01
#' @param cyc \emph{Fabia Parameter:} number of iterations; default = 500 
#' @param spl \emph{Fabia Parameter:} sparseness prior loadings (0 - 2.0); default = 0 (Laplace)
#' @param spz \emph{Fabia Parameter:} sparseness factors (0.5 - 2.0); default = 0.5 (Laplace)
#' @param non_negative \emph{Fabia Parameter:} Non-negative factors and loadings if non_negative > 0; default = 0
#' @param random \emph{Fabia Parameter:} <=0: by SVD, >0: random initialization of loadings in [-random,random]; default = 1.0
#' @param center \emph{Fabia Parameter:} data centering: 1 (mean), 2 (median), > 2 (mode), 0 (no); default = 2
#' @param norm \emph{Fabia Parameter:} data normalization: 1 (0.75-0.25 quantile), >1 (var=1), 0 (no); default = 1
#' @param scale \emph{Fabia Parameter:} loading vectors are scaled in each iteration to the given variance. 0.0 indicates non scaling; default = 0.0
#' @param lap \emph{Fabia Parameter:} minimal value of the variational parameter; default = 1.0
#' @param nL \emph{Fabia Parameter:} maximal number of biclusters at which a row element can participate; default = 0 (no limit)
#' @param lL \emph{Fabia Parameter:} maximal number of row elements per bicluster; default = 0 (no limit)
#' @param bL \emph{Fabia Parameter:} cycle at which the nL or lL maximum starts; default = 0 (start at the beginning)
#' @param which Choose one or more plots to draw: 
#' \enumerate{
#' \item Information Content for Bicluster (Only available for "CSfabia")
#' \item Loadings for reference compounds
#' \item Loadings for Component (Factor/Bicluster) \code{component.plot}
#' \item Gene Scores for Component (Factor/Bicluster) \code{component.Plot}
#' \item Connectivity Ranking Scores for Component \code{component.plot}
#' \item Component \code{component.plot} VS Other Component : Loadings & Genes 
#' \item Profile plot (see \code{profile.type})
#' \item Group Loadings Plots for all components (see \code{grouploadings.labels}).
#' }
#' @param component.plot Which components (Factor/Bicluster) should be investigated? Can be a vector of multiple (e.g. \code{c(1,3,5)}). If \code{NULL}, you can choose components of interest interactively from reference loadings plot.
#' @param CSrank.refplot Logical value deciding if the CS Rank Scores (\code{which=5}) should also be plotted per reference (instead of only the weighted mean).
#' @param column.interest Numeric vector of indices of query columns which should be in the profiles plots (\code{which=7}). If \code{NULL}, you can interactively select genes on the Compound Loadings plot (\code{which=3}).
#' @param row.interest Numeric vector of gene indices to be plotted in gene profiles plot (\code{which=7}, \code{profile.type="gene"}). If \code{NULL}, you can interactively select them in the gene scores plot (\code{which=4}).
#' @param profile.type Type of \code{which=7} plot:
#' \itemize{
#' \item \code{"gene"}: Gene profiles plot of selected genes in \code{row.interest} with the reference compounds and those selected in \code{column.interest} ordered first on the x axis. The other compounds are ordered in decreasing CScore. 
#' \item \code{"cmpd"}: Compound profiles plot of reference and selected compounds (\code{column.interest}) and only those genes on the x-axis which beat the thresholds (\code{gene.thresP}, \code{gene.thresN})
#' }
#' @param color.columns Vector of colors for the reference and query columns (compounds). If \code{NULL}, blue will be used for reference and black for query. Use this option to highlight reference columns and query columns of interest.
#' @param gene.highlight Single numeric vector or list of maximum 5 numeric vectors. This highlights gene of interest in gene scores plot (\code{which=4}) up to 5 different colors. (e.g. You can use this to highlight genes you know to be differentially expressed)
#' @param gene.thresP Threshold for genes with a high score (\code{which=4}).
#' @param gene.thresN Threshold for genes with a low score (\code{which=4}).
#' @param thresP.col Color of genes above \code{gene.thresP}.
#' @param thresN.col Color of genes below \code{gene.thresN}.
#' @param grouploadings.labels This parameter used for the Group Loadings Plots (\code{which=8}). In general this plot will contain the loadings of all factors, grouped and colored by the labels given in this parameter.
#' \itemize{
#' \item If \code{grouploadings.labels!=NULL}:\cr
#' Provide a vector for all samples (ref + query) containing labels on which the plot will be based on.
#' 
#' \item If \code{grouploadings.labels=NULL}: \cr
#' If no labels are provided when choosing \code{which=8}, automatic labels ("Top Samples of Component 1, 2....") will be created. These labels are given to the top \code{grouploadings.cutoff}  number of samples based on the absolute values of the loadings. 
#' }
#' Plot \code{which=8} can be used to check 2 different situations. Either to check if your provided labels coincide with the discovered structure in the analysis. The other aim is to find new interesting structures (of samples) which strongly appear in one or multiple components. A subsequent step could be to take some strong samples/compounds of these compounds and use them as a new reference set in a new CS analysis to check its validity or to find newly connected compounds.
#' 
#' Please note that even when \code{group.loadings.labels!=NULL}, that the labels based on the absolute loadings of all the factors (the top \code{grouploadings.cutoff}) will always be generated and saved in \code{samplefactorlabels} in the \code{extra} slot of the \code{CSresult} object. 
#' This can then later be used for the \code{CSlabelscompare} function to compare them with your true labels.
#' @param grouploadings.cutoff Parameter used in plot \code{which=8}. See \code{grouploadings.labels=NULL} for more information. If this parameter is not provided, it will be automatically set to 10\% of the total number of loadings.
#' @param legend.names Option to draw a legend of for example colored columns in Compound Loadings plot (\code{which=3}). If \code{NULL}, only "References" will be in the legend.
#' @param legend.cols Colors to be used in legends. If \code{NULL}, only blue for "References is used".
#' @param legend.pos Position of the legend in all requested plots, can be \code{"topright"}, \code{"topleft"}, \code{"bottomleft"}, \code{"bottomright"}, \code{"bottom"}, \code{"top"}, \code{"left"}, \code{"right"}, \code{"center"}.
#' @param result.available You can a previously returned object by \code{CSanalysis} in order to only draw graphs, not recompute the scores.
#' @param result.available.update Logical value. If \code{TRUE}, the CS and GS will be overwritten depending on the new \code{component.plot} choice. This would also delete the p-values if \code{permutation.object} was available.
#' @param plot.type How should the plots be outputted? \code{"pdf"} to save them in pdf files, \code{device} to draw them in a graphics device (default), \code{sweave} to use them in a sweave or knitr file.
#' @param basefilename Directory including filename of the graphs if saved in pdf files
#' @return An object of the S4 Class \code{\link{CSresult-class}}.
setMethod('CSanalysis',c('matrix','matrix','CSfabia'),
		function(refMat,querMat,type="CSfabia",p=13,alpha=0.01,cyc=500,spl=0,spz=0.5,non_negative=0,random=1.0,center=2,norm=1,scale=0.0,lap=1.0,nL=0,lL=0,bL=0
				,which=c(2,3,4,5)
				,component.plot=NULL,CSrank.refplot=FALSE
				,column.interest=NULL,row.interest=NULL,profile.type="gene"
				,color.columns=NULL,gene.highlight=NULL
				,gene.thresP=1,gene.thresN=-1,thresP.col="blue",thresN.col="red"
				,grouploadings.labels=NULL,grouploadings.cutoff=NULL
				,legend.names=NULL,legend.cols=NULL,legend.pos="topright"
				,result.available=NULL,result.available.update=FALSE
				,plot.type="device",basefilename=NULL
		) {
		  
		  check_filename(plot.type,basefilename)

			data <- cbind(as.matrix(refMat),as.matrix(querMat))


			colour.columns <- color.columns
			if(!(legend.pos %in% c("topright", "topleft", "bottomleft", "bottomright", "bottom", "top", "left", "right", "center"))){stop(paste0("legend.pos can not be \"",legend.pos,"\""),call.=FALSE)}
			
			analysis.pm.temp <- list(p=p,alpha=alpha,cyc=cyc,spl=spl,spz=spz,non_negative=non_negative,random=random,center=center,norm=norm,scale=scale,lap=lap,nL=nL,lL=lL,bL=bL)
			
			if(!is.null(result.available)){
				if(class(result.available) != "CSresult"){
					stop("result.available is not of the correct class") 
				}
				else if(result.available@type != class(type)){
					stop("result.available is from a different analysis")
				}
				else if(!identical(result.available@call$analysis.pm,analysis.pm.temp)){  # TO AVOID TO REWRITE THE PARAMETERS...
					stop("result.available used other analysis parameters")
				}
			}
		
			
			if(!is.null(column.interest)){column.interest <- column.interest + dim(refMat)[2]}

			out <- analyse_FA(matchcall=type@call,modeltype=class(type),ref.index=c(1:dim(refMat)[2]),
					data=data,p=p,alpha=alpha,cyc=cyc,spl=spl,spz=spz,non_negative=non_negative,random=random,center=center,norm=norm,scale=scale,lap=lap,nL=nL,lL=lL,bL=bL,
					basefilename=basefilename,
#					weighted.data=TRUE,
					component.plot=component.plot,
					column.interest=column.interest,row.interest=row.interest,gene.thresP=gene.thresP,gene.thresN=gene.thresN,
					colour.columns=colour.columns,legend.pos=legend.pos,legend.names=legend.names,legend.cols=legend.cols,thresP.col=thresP.col,thresN.col=thresN.col,
					result.available=result.available,plot.type=plot.type,CSrank.refplot=CSrank.refplot,
					which=which,
					gene.highlight=gene.highlight,profile.type=profile.type,
					grouploadings.labels=grouploadings.labels,grouploadings.cutoff=grouploadings.cutoff,
					result.available.update=result.available.update)


			return(out)
			
		})




# MFA
#' "CSmfa"
#' 
#' @description Doing interactive CS analysis with MFA (Multiple Factor Analysis). Should use multiple references for this analysis.
#' Uses the \code{\link[FactoMineR]{MFA}} function.
#' 
#' @export 
#' @param refMat Reference matrix (Rows = genes and columns = compounds)
#' @param querMat Query matrix
#' @param type \code{"CSmfa"}
#' @param ncp \emph{MFA Parameter:} Number of dimensions kept in the results (by default 5).
#' @param weight.col.mfa \emph{MFA Parameter:} Vector of weights, useful for HMFA method (by default, \code{NULL} and an MFA is performed).
#' @param row.w \emph{MFA Parameter:} An optional row weights (by default, a vector of 1 for uniform row weights).
#' @param mfa.type \emph{MFA Parameter:} The type of column variables (compounds) in both the Query and Reference matrix. "c" or "s" (= default) for quantitative variables (the difference is that for "s" variables are scaled to unit variance), "n" for categorical variables and "f" for frequencies (from a contingency tables)
#' @param which Choose one or more plots to draw: 
#' \enumerate{
#' \item Information Content for Bicluster (Only available for "CSfabia")
#' \item Loadings for reference compounds
#' \item Loadings for Component (Factor/Bicluster) \code{component.plot}
#' \item Gene Scores for Component (Factor/Bicluster) \code{component.Plot}
#' \item Connectivity Ranking Scores for Component \code{component.plot}
#' \item Component \code{component.plot} VS Other Component : Loadings & Genes 
#' \item Profile plot (see \code{profile.type})
#' \item Group Loadings Plots for all components (see \code{grouploadings.labels}).
#' }
#' @param component.plot Which components (Factor/Bicluster) should be investigated? Can be a vector of multiple (e.g. \code{c(1,3,5)}). If \code{NULL}, you can choose components of interest interactively from reference loadings plot.
#' @param CSrank.refplot Logical value deciding if the CS Rank Scores (\code{which=5}) should also be plotted per reference (instead of only the weighted mean).
#' @param column.interest Numeric vector of indices of query columns which should be in the profiles plots (\code{which=7}). If \code{NULL}, you can interactively select genes on the Compound Loadings plot (\code{which=3}).
#' @param row.interest Numeric vector of gene indices to be plotted in gene profiles plot (\code{which=7}, \code{profile.type="gene"}). If \code{NULL}, you can interactively select them in the gene scores plot (\code{which=4}).
#' @param profile.type Type of \code{which=7} plot:
#' \itemize{
#' \item \code{"gene"}: Gene profiles plot of selected genes in \code{row.interest} with the reference compounds and those selected in \code{column.interest} ordered first on the x axis. The other compounds are ordered in decreasing CScore. 
#' \item \code{"cmpd"}: Compound profiles plot of reference and selected compounds (\code{column.interest}) and only those genes on the x-axis which beat the thresholds (\code{gene.thresP}, \code{gene.thresN})
#' }
#' @param color.columns Vector of colors for the reference and query columns (compounds). If \code{NULL}, blue will be used for reference and black for query. Use this option to highlight reference columns and query columns of interest.
#' @param gene.highlight Single numeric vector or list of maximum 5 numeric vectors. This highlights gene of interest in gene scores plot (\code{which=4}) up to 5 different colors. (e.g. You can use this to highlight genes you know to be differentially expressed)
#' @param gene.thresP Threshold for genes with a high score (\code{which=4}).
#' @param gene.thresN Threshold for genes with a low score (\code{which=4}).
#' @param thresP.col Color of genes above \code{gene.thresP}.
#' @param thresN.col Color of genes below \code{gene.thresN}.
#' @param grouploadings.labels This parameter used for the Group Loadings Plots (\code{which=8}). In general this plot will contain the loadings of all factors, grouped and colored by the labels given in this parameter.
#' \itemize{
#' \item If \code{grouploadings.labels!=NULL}:\cr
#' Provide a vector for all samples (ref + query) containing labels on which the plot will be based on.
#' 
#' \item If \code{grouploadings.labels=NULL}: \cr
#' If no labels are provided when choosing \code{which=8}, automatic labels ("Top Samples of Component 1, 2....") will be created. These labels are given to the top \code{grouploadings.cutoff}  number of samples based on the absolute values of the loadings. 
#' }
#' Plot \code{which=8} can be used to check 2 different situations. Either to check if your provided labels coincide with the discovered structure in the analysis. The other aim is to find new interesting structures (of samples) which strongly appear in one or multiple components. A subsequent step could be to take some strong samples/compounds of these compounds and use them as a new reference set in a new CS analysis to check its validity or to find newly connected compounds.
#' 
#' Please note that even when \code{group.loadings.labels!=NULL}, that the labels based on the absolute loadings of all the factors (the top \code{grouploadings.cutoff}) will always be generated and saved in \code{samplefactorlabels} in the \code{extra} slot of the \code{CSresult} object. 
#' This can then later be used for the \code{CSlabelscompare} function to compare them with your true labels.
#' @param grouploadings.cutoff Parameter used in plot \code{which=8}. See \code{grouploadings.labels=NULL} for more information. If this parameter is not provided, it will be automatically set to 10\% of the total number of loadings.
#' @param legend.names Option to draw a legend of for example colored columns in Compound Loadings plot (\code{which=3}). If \code{NULL}, only "References" will be in the legend.
#' @param legend.cols Colors to be used in legends. If \code{NULL}, only blue for "References is used".
#' @param legend.pos Position of the legend in all requested plots, can be \code{"topright"}, \code{"topleft"}, \code{"bottomleft"}, \code{"bottomright"}, \code{"bottom"}, \code{"top"}, \code{"left"}, \code{"right"}, \code{"center"}.
#' @param result.available You can a previously returned object by \code{CSanalysis} in order to only draw graphs, not recompute the scores.
#' @param result.available.update Logical value. If \code{TRUE}, the CS and GS will be overwritten depending on the new \code{component.plot} choice. This would also delete the p-values if \code{permutation.object} was available.
#' @param plot.type How should the plots be outputted? \code{"pdf"} to save them in pdf files, \code{device} to draw them in a graphics device (default), \code{sweave} to use them in a sweave or knitr file.
#' @param basefilename Directory including filename of the graphs if saved in pdf files
#' @return An object of the S4 Class \code{\link{CSresult-class}}.
setMethod("CSanalysis",c("matrix","matrix","CSmfa"),function(
				refMat,querMat,type="CSmfa",ncp=5,weight.col.mfa=NULL,row.w=NULL,
				mfa.type="s",
				which=c(2,3,4,5)
				,component.plot=NULL,CSrank.refplot=FALSE
				,column.interest=NULL,row.interest=NULL,profile.type="gene",
				color.columns=NULL,gene.highlight=NULL,				
				gene.thresP=1,gene.thresN=-1,thresP.col="blue",thresN.col="red",
				grouploadings.labels=NULL,grouploadings.cutoff=NULL,
				legend.names=NULL,legend.cols=NULL,legend.pos="topright",				
				result.available=NULL,result.available.update=FALSE,
				plot.type="device",basefilename=NULL
				){
  
      if(length(mfa.type)!=1){stop("mfa.type needs to be of length 1")}
      if(class(mfa.type)!="character"){stop("mfa.type needs to be a character")}
      if(!(mfa.type %in% c("c","s","n","f"))){stop("mfa.type needs to be \"c\", \"s\", \"n\" or \"f\"")}
  
      check_filename(plot.type,basefilename)
				
			if(!(dim(refMat)[2]>1)){
				stop("Reference matrix should have more than 1 reference")
			}
			if(dim(querMat)[2]==1){
				stop("Query matrix only has 1 query.")
			}
			
			data <- cbind(as.matrix(refMat),as.matrix(querMat))
			
			
			colour.columns <- color.columns
			if(!(legend.pos %in% c("topright", "topleft", "bottomleft", "bottomright", "bottom", "top", "left", "right", "center"))){stop(paste0("legend.pos can not be \"",legend.pos,"\""),call.=FALSE)}
			
			
			analysis.pm.temp <- list(ncp=ncp,weight.col.mfa=weight.col.mfa,row.w=row.w)
			
			if(!is.null(result.available)){
				if(class(result.available) != "CSresult"){
					stop("result.available is not of the correct class") 
				}
				else if(result.available@type != class(type)){
					stop("result.available is from a different analysis")
				}
				else if(!identical(result.available@call$analysis.pm,analysis.pm.temp)){
					stop("result.available used other analysis parameters")
				}
			}

			
			if(!is.null(column.interest)){column.interest <- column.interest + dim(refMat)[2]}
			
			out <- analyse_FA(matchcall=type@call,modeltype=class(type),ref.index=c(1:dim(refMat)[2]),
					data=data,type.mfa=rep(mfa.type,2),ind.sup=NULL,ncp=ncp,name.group=c("Reference","Query"),num.group.sup=NULL,graph=FALSE,weight.col.mfa=weight.col.mfa,row.w=row.w,axes=c(1,2),tab.comp=NULL,
					basefilename=basefilename,
					component.plot=component.plot,column.interest=column.interest,row.interest=row.interest,gene.thresP=gene.thresP,gene.thresN=gene.thresN,
					colour.columns=colour.columns,legend.pos=legend.pos,legend.names=legend.names,legend.cols=legend.cols,thresP.col=thresP.col,thresN.col=thresN.col,
					result.available=result.available,plot.type=plot.type,CSrank.refplot=CSrank.refplot,
					which=which,
					grouploadings.labels=grouploadings.labels,grouploadings.cutoff=grouploadings.cutoff,
					gene.highlight=gene.highlight,profile.type=profile.type,result.available.update=result.available.update)
			
			out <- trim_object(out)
			
		return(out)
							
		})



# PCA
#' "CSpca"
#' 
#' @description Doing interactive CS analysis with PCA (Principal Component Analysis). This analysis is meant for 1 reference signature.
#' Uses the \code{\link[FactoMineR]{PCA}} function.
#' 
#' @export 
#' @param refMat Reference matrix (Rows = genes and columns = compounds)
#' @param querMat Query matrix
#' @param type \code{"CSpca"}
#' @param ncp \emph{PCA Parameter:} Number of dimensions kept in the results (by default 5).
#' @param scale.unit \emph{PCA Parameter:} A boolean, if TRUE (value set by default) then data are scaled to unit variance.
#' @param row.w \emph{PCA Parameter:} An optional row weights (by default, a vector of 1 for uniform row weights).
#' @param col.w \emph{PCA Parameter:} An optional column weights (by default, uniform column weights).
#' @param which Choose one or more plots to draw: 
#' \enumerate{
#' \item Information Content for Bicluster (Only available for "CSfabia")
#' \item Loadings for reference compounds
#' \item Loadings for Component (Factor/Bicluster) \code{component.plot}
#' \item Gene Scores for Component (Factor/Bicluster) \code{component.Plot}
#' \item Connectivity Ranking Scores for Component \code{component.plot}
#' \item Component \code{component.plot} VS Other Component : Loadings & Genes 
#' \item Profile plot (see \code{profile.type})
#' \item Group Loadings Plots for all components (see \code{grouploadings.labels}).
#' }
#' @param component.plot Which components (Factor/Bicluster) should be investigated? Can be a vector of multiple (e.g. \code{c(1,3,5)}). If \code{NULL}, you can choose components of interest interactively from reference loadings plot.
#' @param CSrank.refplot Logical value deciding if the CS Rank Scores (\code{which=5}) should also be plotted per reference (instead of only the weighted mean).
#' @param column.interest Numeric vector of indices of query columns which should be in the profiles plots (\code{which=7}). If \code{NULL}, you can interactively select genes on the Compound Loadings plot (\code{which=3}).
#' @param row.interest Numeric vector of gene indices to be plotted in gene profiles plot (\code{which=7}, \code{profile.type="gene"}). If \code{NULL}, you can interactively select them in the gene scores plot (\code{which=4}).
#' @param profile.type Type of \code{which=7} plot:
#' \itemize{
#' \item \code{"gene"}: Gene profiles plot of selected genes in \code{row.interest} with the reference compounds and those selected in \code{column.interest} ordered first on the x axis. The other compounds are ordered in decreasing CScore. 
#' \item \code{"cmpd"}: Compound profiles plot of reference and selected compounds (\code{column.interest}) and only those genes on the x-axis which beat the thresholds (\code{gene.thresP}, \code{gene.thresN})
#' }
#' @param color.columns Vector of colors for the reference and query columns (compounds). If \code{NULL}, blue will be used for reference and black for query. Use this option to highlight reference columns and query columns of interest.
#' @param gene.highlight Single numeric vector or list of maximum 5 numeric vectors. This highlights gene of interest in gene scores plot (\code{which=4}) up to 5 different colors. (e.g. You can use this to highlight genes you know to be differentially expressed)
#' @param gene.thresP Threshold for genes with a high score (\code{which=4}).
#' @param gene.thresN Threshold for genes with a low score (\code{which=4}).
#' @param thresP.col Color of genes above \code{gene.thresP}.
#' @param thresN.col Color of genes below \code{gene.thresN}.
#' @param grouploadings.labels This parameter used for the Group Loadings Plots (\code{which=8}). In general this plot will contain the loadings of all factors, grouped and colored by the labels given in this parameter.
#' \itemize{
#' \item If \code{grouploadings.labels!=NULL}:\cr
#' Provide a vector for all samples (ref + query) containing labels on which the plot will be based on.
#' 
#' \item If \code{grouploadings.labels=NULL}:\cr 
#' If no labels are provided when choosing \code{which=8}, automatic labels ("Top Samples of Component 1, 2....") will be created. These labels are given to the top \code{grouploadings.cutoff}  number of samples based on the absolute values of the loadings. 
#' }
#' Plot \code{which=8} can be used to check 2 different situations. Either to check if your provided labels coincide with the discovered structure in the analysis. The other aim is to find new interesting structures (of samples) which strongly appear in one or multiple components. A subsequent step could be to take some strong samples/compounds of these compounds and use them as a new reference set in a new CS analysis to check its validity or to find newly connected compounds.
#' 
#' Please note that even when \code{group.loadings.labels!=NULL}, that the labels based on the absolute loadings of all the factors (the top \code{grouploadings.cutoff}) will always be generated and saved in \code{samplefactorlabels} in the \code{extra} slot of the \code{CSresult} object. 
#' This can then later be used for the \code{CSlabelscompare} function to compare them with your true labels.
#' @param grouploadings.cutoff Parameter used in plot \code{which=8}. See \code{grouploadings.labels=NULL} for more information. If this parameter is not provided, it will be automatically set to 10\% of the total number of loadings.
#' @param legend.names Option to draw a legend of for example colored columns in Compound Loadings plot (\code{which=3}). If \code{NULL}, only "References" will be in the legend.
#' @param legend.cols Colors to be used in legends. If \code{NULL}, only blue for "References is used".
#' @param legend.pos Position of the legend in all requested plots, can be \code{"topright"}, \code{"topleft"}, \code{"bottomleft"}, \code{"bottomright"}, \code{"bottom"}, \code{"top"}, \code{"left"}, \code{"right"}, \code{"center"}.
#' @param result.available You can a previously returned object by \code{CSanalysis} in order to only draw graphs, not recompute the scores.
#' @param result.available.update Logical value. If \code{TRUE}, the CS and GS will be overwritten depending on the new \code{component.plot} choice. This would also delete the p-values if \code{permutation.object} was available.
#' @param plot.type How should the plots be outputted? \code{"pdf"} to save them in pdf files, \code{device} to draw them in a graphics device (default), \code{sweave} to use them in a sweave or knitr file.
#' @param basefilename Directory including filename of the graphs if saved in pdf files.
#' @return An object of the S4 Class \code{\link{CSresult-class}}.
setMethod("CSanalysis",c("matrix","matrix","CSpca"),function(
				refMat,querMat,type="CSpca",ncp=5,scale.unit=TRUE,row.w=NULL,col.w=NULL,
				which=c(2,3,4,5),
				component.plot=NULL,CSrank.refplot=FALSE,
				column.interest=NULL,row.interest=NULL,profile.type="gene",
				color.columns=NULL,gene.highlight=NULL,				
				gene.thresP=1,gene.thresN=-1,thresP.col="blue",thresN.col="red",
				grouploadings.labels=NULL,grouploadings.cutoff=NULL,
				legend.names=NULL,legend.cols=NULL,legend.pos="topright",
				result.available=NULL,result.available.update=FALSE,
				plot.type="device",basefilename=NULL
				){
  
      check_filename(plot.type,basefilename)
  
			if((dim(refMat)[2]!=1)){
				stop("Reference matrix should have only 1 reference")
			}
			if(dim(querMat)[2]==1){
				stop("Query matrix only has 1 query.")
			}
					
			data <- cbind(as.matrix(refMat),as.matrix(querMat))
									
			colour.columns <- color.columns
			if(!(legend.pos %in% c("topright", "topleft", "bottomleft", "bottomright", "bottom", "top", "left", "right", "center"))){stop(paste0("legend.pos can not be \"",legend.pos,"\""),call.=FALSE)}
			
					
			analysis.pm.temp <- list(ncp=ncp,scale.unit=scale.unit,row.w=row.w,col.w=col.w)
			
			if(!is.null(result.available)){
				if(class(result.available) != "CSresult"){
					stop("result.available is not of the correct class") 
				}
				else if(result.available@type != class(type)){
					stop("result.available is from a different analysis")
				}
				else if(!identical(result.available@call$analysis.pm,analysis.pm.temp)){
					stop("result.available used other analysis parameters")
				}
			}
	
			
			if(!is.null(column.interest)){column.interest <- column.interest + dim(refMat)[2]}
				
			out <- analyse_FA(matchcall=type@call,modeltype=class(type),ref.index=c(1:dim(refMat)[2]),
					data=data, 	scale.unit = scale.unit, ncp = ncp, ind.sup = NULL,
					quanti.sup = NULL, quali.sup = NULL, row.w = row.w,
					col.w = col.w, graph = FALSE, axes = c(1,2),
					basefilename=basefilename,
					component.plot=component.plot,column.interest=column.interest,gene.thresP=gene.thresP,gene.thresN=gene.thresN,
					colour.columns=colour.columns,legend.pos=legend.pos,legend.names=legend.names,legend.cols=legend.cols,thresP.col=thresP.col,thresN.col=thresN.col,
					result.available=result.available,plot.type=plot.type,which=which,CSrank.refplot=CSrank.refplot,
					profile.type=profile.type,gene.highlight=gene.highlight,
					grouploadings.labels=grouploadings.labels,grouploadings.cutoff=grouploadings.cutoff,
					result.available.update=result.available.update)
			
			out <- trim_object(out)
			
			return(out)
			
					
		})

		
# sMFA
#' "CSsmfa"
#' 
#' @description
#' Doing interactive CS analysis with sMFA (Sparse Multiple Factor Analysis). Should use multiple references for this analysis.
#' Either \code{\link[elasticnet]{spca}} or \code{\link[elasticnet]{arrayspc}} is used.
#' 
#' @export 
#' @param refMat Reference matrix (Rows = genes and columns = compounds)
#' @param querMat Query matrix
#' @param type \code{"CSsmfa"}
#' @param K \emph{sMFA Parameters:} Number of components.
#' @param para \emph{sMFA Parameters:} A vector of length K. All elements should be positive. If \code{sparse="varnum"}, the elements integers. 
#' @param lambda \emph{sMFA Parameters:} Quadratic penalty parameter. Default value is 1e-6. If the target dimension of the sparsness is higher than the other dimension (p > n), it is advised to put \code{lambda} to \code{Inf} which uses the \code{arrayspc} algorithm optimized for this case. For the other case, p < n, a zero or positive \code{lambda} is sufficient and will utilize the normal \code{spca} algorithm.
#' @param sparse.dim \emph{sMFA Parameters:} Which dimension should be sparse? 1: Rows, 2: Columns (default) (Note: For Connectivity Scores it is advised to apply sparsity on the compounds/columns)
#' @param sparse \emph{sMFA Parameters (\code{lambda < Inf} only):} If \code{sparse="penalty"}, \code{para} is a vector of 1-norm penalty parameters. If \code{sparse="varnum"}, \code{para} defines the number of sparse loadings to be obtained.
#' @param max.iter \emph{sMFA Parameters:} Maximum number of iterations.
#' @param eps.conv \emph{sMFA Parameters:} Convergence criterion.
#' @param which Choose one or more plots to draw: 
#' \enumerate{
#' \item Information Content for Bicluster (Only available for "CSfabia")
#' \item Loadings for reference compounds
#' \item Loadings for Component (Factor/Bicluster) \code{component.plot}
#' \item Gene Scores for Component (Factor/Bicluster) \code{component.Plot}
#' \item Connectivity Ranking Scores for Component \code{component.plot}
#' \item Component \code{component.plot} VS Other Component : Loadings & Genes 
#' \item Profile plot (see \code{profile.type})
#' \item Group Loadings Plots for all components (see \code{grouploadings.labels}).
#' }
#' @param component.plot Which components (Factor/Bicluster) should be investigated? Can be a vector of multiple (e.g. \code{c(1,3,5)}). If \code{NULL}, you can choose components of interest interactively from reference loadings plot.
#' @param CSrank.refplot Logical value deciding if the CS Rank Scores (\code{which=5}) should also be plotted per reference (instead of only the weighted mean).
#' @param column.interest Numeric vector of indices of query columns which should be in the profiles plots (\code{which=7}). If \code{NULL}, you can interactively select genes on the Compound Loadings plot (\code{which=3}).
#' @param row.interest Numeric vector of gene indices to be plotted in gene profiles plot (\code{which=7}, \code{profile.type="gene"}). If \code{NULL}, you can interactively select them in the gene scores plot (\code{which=4}).
#' @param profile.type Type of \code{which=7} plot:
#' \itemize{
#' \item \code{"gene"}: Gene profiles plot of selected genes in \code{row.interest} with the reference compounds and those selected in \code{column.interest} ordered first on the x axis. The other compounds are ordered in decreasing CScore. 
#' \item \code{"cmpd"}: Compound profiles plot of reference and selected compounds (\code{column.interest}) and only those genes on the x-axis which beat the thresholds (\code{gene.thresP}, \code{gene.thresN})
#' }
#' @param color.columns Vector of colors for the reference and query columns (compounds). If \code{NULL}, blue will be used for reference and black for query. Use this option to highlight reference columns and query columns of interest.
#' @param gene.highlight Single numeric vector or list of maximum 5 numeric vectors. This highlights gene of interest in gene scores plot (\code{which=4}) up to 5 different colors. (e.g. You can use this to highlight genes you know to be differentially expressed)
#' @param gene.thresP Threshold for genes with a high score (\code{which=4}).
#' @param gene.thresN Threshold for genes with a low score (\code{which=4}).
#' @param thresP.col Color of genes above \code{gene.thresP}.
#' @param thresN.col Color of genes below \code{gene.thresN}.
#' @param grouploadings.labels This parameter used for the Group Loadings Plots (\code{which=8}). In general this plot will contain the loadings of all factors, grouped and colored by the labels given in this parameter.
#' \itemize{
#' \item If \code{grouploadings.labels!=NULL}:\cr
#' Provide a vector for all samples (ref + query) containing labels on which the plot will be based on.
#' 
#' \item If \code{grouploadings.labels=NULL}: \cr
#' If no labels are provided when choosing \code{which=8}, automatic labels ("Top Samples of Component 1, 2....") will be created. These labels are given to the top \code{grouploadings.cutoff}  number of samples based on the absolute values of the loadings. 
#' }
#' Plot \code{which=8} can be used to check 2 different situations. Either to check if your provided labels coincide with the discovered structure in the analysis. The other aim is to find new interesting structures (of samples) which strongly appear in one or multiple components. A subsequent step could be to take some strong samples/compounds of these compounds and use them as a new reference set in a new CS analysis to check its validity or to find newly connected compounds.
#' 
#' Please note that even when \code{group.loadings.labels!=NULL}, that the labels based on the absolute loadings of all the factors (the top \code{grouploadings.cutoff}) will always be generated and saved in \code{samplefactorlabels} in the \code{extra} slot of the \code{CSresult} object. 
#' This can then later be used for the \code{CSlabelscompare} function to compare them with your true labels.
#' @param grouploadings.cutoff Parameter used in plot \code{which=8}. See \code{grouploadings.labels=NULL} for more information. If this parameter is not provided, it will be automatically set to 10\% of the total number of loadings.
#' @param legend.names Option to draw a legend of for example colored columns in Compound Loadings plot (\code{which=3}). If \code{NULL}, only "References" will be in the legend.
#' @param legend.cols Colors to be used in legends. If \code{NULL}, only blue for "References is used".
#' @param legend.pos Position of the legend in all requested plots, can be \code{"topright"}, \code{"topleft"}, \code{"bottomleft"}, \code{"bottomright"}, \code{"bottom"}, \code{"top"}, \code{"left"}, \code{"right"}, \code{"center"}.
#' @param result.available You can a previously returned object by \code{CSanalysis} in order to only draw graphs, not recompute the scores.
#' @param result.available.update Logical value. If \code{TRUE}, the CS and GS will be overwritten depending on the new \code{component.plot} choice. This would also delete the p-values if \code{permutation.object} was available.
#' @param plot.type How should the plots be outputted? \code{"pdf"} to save them in pdf files, \code{device} to draw them in a graphics device (default), \code{sweave} to use them in a sweave or knitr file.
#' @param basefilename Directory including filename of the graphs if saved in pdf files#' @return An object of the S4 Class \code{CSresult}.
#' @return An object of the S4 Class \code{\link{CSresult-class}}.
setMethod("CSanalysis",c("matrix","matrix","CSsmfa"),function(
				refMat,querMat,type="Csmfa",K=15,para,lambda=1e-6,sparse.dim=2,sparse="penalty",max.iter=200,eps.conv=1e-3,
				which=c(2,3,4,5),
				component.plot=NULL,CSrank.refplot=FALSE,
				column.interest=NULL,row.interest=NULL,profile.type="gene",
				color.columns=NULL,gene.highlight=NULL,
				gene.thresP=1,gene.thresN=-1,thresP.col="blue",thresN.col="red",
				grouploadings.labels=NULL,grouploadings.cutoff=NULL,
				legend.names=NULL,legend.cols=NULL,legend.pos="topright",
				result.available=NULL,result.available.update=FALSE,
				plot.type="device",basefilename=NULL
				){
			
      check_filename(plot.type,basefilename)
			
      if(!(dim(refMat)[2]>1)){
				stop("Reference matrix should have more than 1 reference")
			}
					
			data <- cbind(as.matrix(refMat),as.matrix(querMat))
					
			colour.columns <- color.columns
			if(!(legend.pos %in% c("topright", "topleft", "bottomleft", "bottomright", "bottom", "top", "left", "right", "center"))){stop(paste0("legend.pos can not be \"",legend.pos,"\""),call.=FALSE)}
			
					
			analysis.pm.temp <- list(K=K,lambda=lambda,sparse=sparse,max.iter=max.iter,eps.conv=eps.conv,para=para,sparse.dim=sparse.dim)
			
			if(!is.null(result.available)){
				if(class(result.available) != "CSresult"){
					stop("result.available is not of the correct class") 
				}
				else if(result.available@type != class(type)){
					stop("result.available is from a different analysis")
				}
				else if(!identical(result.available@call$analysis.pm,analysis.pm.temp)){
					stop("result.available used other analysis parameters")
				}
			}

			
			if(!is.null(column.interest)){column.interest <- column.interest + dim(refMat)[2]}
			
			out <- analyse_FA(matchcall=type@call,modeltype=class(type),ref.index=c(1:dim(refMat)[2]),
					data=data,K=K,para=para,type.smfa="predictor",sparse=sparse,use.corr=FALSE,lambda=lambda,max.iter=max.iter,trace=FALSE,eps.conv=eps.conv,sparse.dim=sparse.dim,
					basefilename=basefilename,
					component.plot=component.plot,column.interest=column.interest,gene.thresP=gene.thresP,gene.thresN=gene.thresN,
					colour.columns=colour.columns,legend.pos=legend.pos,legend.names=legend.names,legend.cols=legend.cols,thresP.col=thresP.col,thresN.col=thresN.col,
					result.available=result.available,plot.type=plot.type,CSrank.refplot=CSrank.refplot,
					which=which,profile.type=profile.type,gene.highlight=gene.highlight,
					grouploadings.labels=grouploadings.labels,grouploadings.cutoff=grouploadings.cutoff,
					row.interest=row.interest,result.available.update=result.available.update)
			
			return(out)
						
		})


# Zhang and Gant
#' "CSzhang"
#' 
#' Compute the Connectivity Scores by Zhang and Gant (2008). One or multiple reference compounds are possible in this analysis.
#' 
#' @export 
#' @param refMat Reference matrix (Rows = genes and columns = compounds)
#' @param querMat Query matrix
#' @param type \code{"CSzhang"}
#' @param nref \emph{Zhang Parameter:} Number of top up- and downregulated genes in reference signature. If \code{NULL}, all rows (genes) are used.
#' @param nquery \emph{Zhang Parameter:} Number of top up- and downregulated genes in query signature. If \code{NULL}, all rows (genes) are used. (Note that \eqn{nref >= nquery})
#' @param ord.query \emph{Zhang Parameter:} Logical value. Should the query signature be treated as ordered?
#' @param which Choose plot to draw.
#' \enumerate{
#' \item Zhang and Gant Scores Plot
#' }
# @param permute \emph{Zhang Parameter:} Logical value. Should p-values be computed through permutation?
# @param B \emph{Zhang Parameter:} Number of permutations for p-value calculation.
# @param ntop.pvalues \emph{Zhang Parameter:} Number of top p-values to be reported first. 
#' @param ntop.scores \emph{Zhang Parameter:} Number of top positive and negative CS to be reported first.
#' @param color.query Vector of colors for the query columns. You can use this option to highlight columns(compounds) of interest in the CS plot. (This does not include the reference columns since they are not included in the CS plot.)
#' @param legend.names Option to draw a legend (about the highlights in \code{color.query}) in the CS plot. If \code{NULL}, no legend will be drawn.
#' @param legend.cols Colors to be used for the \code{legend.names}.
#' @param legend.pos Position of the legend in all requested plots, can be \code{"topright"}, \code{"topleft"}, \code{"bottomleft"}, \code{"bottomright"}, \code{"bottom"}, \code{"top"}, \code{"left"}, \code{"right"}, \code{"center"}.
#' @param result.available You can a previously returned object by \code{CSanalysis} in order to only draw graphs, not recompute the scores. If this object also contains the permutation object, in the score plot the values with a (adjusted) pvalue smaller than 0.05 will be colored purple.
#' @param result.available.update Logical value. If \code{TRUE}, the CS and GS will be overwritten depending on the new \code{component.plot} choice. This would also delete the p-values if \code{permutation.object} was available.
#' @param plot.type How should the plots be outputted? \code{"pdf"} to save them in pdf files, \code{device} to draw them in a graphics device (default), \code{sweave} to use them in a sweave or knitr file.
#' @param basefilename Directory including filename of the graphs if saved in pdf files
#' @return An object of the S4 Class \code{\link{CSresult-class}}. The CS slot will also contain the top positive and negative scores as well as the top p-values. The GS slot will be empty for Zhang and Gant.
setMethod("CSanalysis",c("matrix","matrix","CSzhang"),function(refMat,querMat,type="CSzhang",
				nref=NULL,nquery=NULL,ord.query=TRUE,ntop.scores=20,
				which=c(1),
#				B=100000,ntop.pvalues=20,permute=FALSE,
				color.query=NULL,
				legend.names=NULL,legend.cols=NULL,legend.pos="topright",
				result.available=NULL,result.available.update=FALSE,
				plot.type="device",basefilename=NULL
				){
			
      check_filename(plot.type,basefilename)
		
    	colour.query <- color.query
			if(is.null(legend.cols)){legend.cols <- "black"}
			
			if(!(legend.pos %in% c("topright", "topleft", "bottomleft", "bottomright", "bottom", "top", "left", "right", "center"))){stop(paste0("legend.pos can not be \"",legend.pos,"\""),call.=FALSE)}
			
			
			if(!is.null(result.available)){
				if(class(result.available) != "CSresult"){
					stop("result.available is not of the correct class") 
				}
				else if(result.available@type != class(type)){
					stop("result.available is from a different analysis")
				}
			}

			
			out <- analyse_zhang(dataref=refMat,dataquery=querMat,nref=nref,nquery=nquery,ord.query=ord.query,ntop.pvalues=20,ntop.scores=ntop.scores,
#					permute=permute,B=B,ntop.pvalues=ntop.values,
					basefilename=basefilename,
					colour.query=colour.query,legend.names=legend.names,legend.cols=legend.cols,legend.pos=legend.pos,
					result.available=result.available,plot.type=plot.type,
					which=which)
			
			
			CS <- list(CS.query=data.frame(
							ZGscore=out$All[,1],
							ZGrank=as.integer(rank(-out$All[,1])),
							ZGabsrank=as.integer(rank(-abs(out$All[,1]))),
							row.names=rownames(out$All)
							))
			
			
			call.object <- list(match.call=type@call,analysis.pm=list(nref=nref,nquery=nquery,ord.query=ord.query,ntop.scores=ntop.scores))
			call.object$dimensions <- list(row=dim(refMat)[1],col=c(ref=dim(refMat)[2],query=dim(querMat)[2]))
			
			
			if(!is.null(result.available)){ 
				if(result.available.update){ # Chance to update object
					if(!is.null(result.available@permutation.object) ){	warning("P-values were deleted in the CS and permutation.object slot ")}
					return(new("CSresult",type=class(type),CS=CS,GS=data.frame(),extra=list(CSRank.Full=NULL,object=out),permutation.object=NULL,call=call.object))
				}else{
					return(new("CSresult",type=class(type),CS=result.available@CS,GS=data.frame(),extra=result.available@extra,permutation.object=result.available@permutation.object,call=result.available@call))
				}
			}else{
				return(new("CSresult",type=class(type),CS=CS,GS=data.frame(),extra=list(CSRank.Full=NULL,object=out),permutation.object=NULL,call=call.object))
			}
			
		})

