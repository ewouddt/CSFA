# Project: CSFA
# 
# Author: lucp8394
###############################################################################

# @importFrom parallel detectCores splitIndices clusterCall clusterExport
# @importFrom snowFT clusterApplyFT makeClusterFT stopClusterFT
# @importfrom snow


#' Permute CS results
#' 
#' Apply permutation on MFA or Zhang results to obtain p-values of 1 of the components. 
#' The function asks for a CSresult object which is returned by CSanalysis. The CSpermute function will return the same CSresult object with added information such as p-values.
#' If asked, the CSpermute function will also draw a volcanoplot and/or histograms of the p-values. If you simply want to redraw these plots, simply use the returned CSresult object by CSpermute again in the CSpermute function.
#' If the number of permutations was not changed, this will prevent the entire permutation analysis from being redone.
#' 
#' \bold{IMPORTANT!} For MFA, \code{CSpermute} should \emph{only} be used to compute the p-values of the Component in which the structure (loadings) of the references is the strongest.
#' This because in each permutation the factor with the highest average reference loadings will be chosen. 
#' The ability to compute p-values of other factors (in which the reference set also increased loadings) will be added in a later release. 
#' 
#' @export
#' @param refMat Reference matrix (Rows = genes and columns = compounds).
#' @param querMat Query matrix
#' @param CSresult A CSresult class object.
#' @param B Number of permutations.
#' @param mfa.factor If permuting a CSmfa result, mfa.factor will decide of which factor the p-values should be computed. If \code{NULL}, the factor chosen in CSanalysis will be chosen (the factor chosen in the CS slot of the CSresult). NOTE: If the mfa.factor is different from the factor in the CS slot, the CS slot will be overwritten with this new factor.
#' @param method.adjust Correction method of multiplicity adjusted p-values: "none", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY" or "fdr". (Raw p-values are also always provided)
#' @param verbose If \code{TRUE}, progression dots	 of the permutation analysis will be printed.
# @param querMat.Perm Possible to provide user-created permuted querMat data. Should be a list object of B times permuted querMat matrices.
# @param save.querMat.Perm If \code{TRUE}, the list of permuted querMat matrices will be saved in the permutation.object slot of the CSresult.
#' @param which Choose which plot to draw:
#' \enumerate{
#' \item A volcano plot of the -log(p-values) versus the observed connection scores.
#' \item A histogram of the permuted connection scores under the null hypothesis for a specific compound. A vertical line(s) is added for the observed CS and its p-value. The \code{cmpd.hist} parameter determines which compounds are drawn like this.
#' \item Analog to \code{which=1}, but for CSRankScores.
#' \item Analog to \code{which=2}, but for CSRankScores.
#' }
#' @param cmpd.hist Query index vector which decides which query compounds are plotted for the histogram distribution under null hypothesis (\code{which=2}). If \code{NULL}, you can select which compounds you want interactively on the volcano plot.
#' @param plot.type How should the plots be outputted? \code{"pdf"} to save them in pdf files, \code{device} to draw them in a graphics device (default), \code{sweave} to use them in a sweave or knitr file.
#' @param color.columns Option to color the compounds on the volcano plot (\code{which=1}). Should be a vector of colors with the length of number of queries.
#' @param basefilename Basename of the graphs if saved in pdf files
#' @param MultiCores Logical value parallelisation should be used for permutation. \code{FALSE} by default. (This option uses \code{\link[snowFT]{clusterApplyFT}} in order to provide load balancing and reproducible results with \code{\link{MultiCores.seed}})
#' @param MultiCores.number Number of cores to be used for \code{MultiCores=TRUE}. By default total number of physical cores.
#' @param MultiCores.seed Seed to be used for \code{MultiCores=TRUE} using (\code{\link[snowFT]{see clusterSetupRNG.FT}})
#' @return Returns the same CSresult object with added p-values to the CS slot and added information to the permutation.object slot. This CSresult can be reused in CSpermute to redraw the plots without calculation.
#' @examples
#' \dontrun{
#' data("dataSIM",package="CSFA")
#' Mat1 <- dataSIM[,c(1:6)]
#' Mat2 <- dataSIM[,-c(1:6)]
#' 
#' MFA_analysis <- CSanalysis(Mat1,Mat2,"CSmfa")
#' MFA_analysis <- CSpermute(Mat1,Mat2,MFA_analysis,B=200)
#' }
CSpermute <- function(refMat,querMat,CSresult,B=500,mfa.factor=NULL,method.adjust="none",verbose=TRUE,
#		querMat.Perm=NULL,save.querMat.Perm=FALSE,
		which=c(1,3),cmpd.hist=NULL,color.columns=NULL,plot.type="device",basefilename="CSpermute",
    MultiCores=FALSE,MultiCores.number=detectCores(logical=FALSE),MultiCores.seed=NULL
  ){
	
	if(class(CSresult)!="CSresult"){stop("CSresult is not a class object of CSresult")}
	if(!(method.adjust %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"))){stop("Incorrect method.adjust. Should we one of the following 'none', 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr'")}
	
	type <- CSresult@type
	ref.index <- c(1:dim(refMat)[2])
	
	if(!(type %in% c("CSmfa","CSzhang"))){stop("Permutation is only available for MFA and Zhang results")}
	

	if(MultiCores){
	  if(MultiCores.number==1){
	    warning("Only 1 core was chosen, no parallelisation will be used.")
	    MultiCores <- FALSE
	  }
	}
	if(MultiCores & is.null(MultiCores.seed)){MultiCores.seed <- sample(1:9999999,1)}
	
	
	
	# redo permutation and analysis if object not available of different number of permutations is asked
	if(is.null(CSresult@permutation.object) | (length(CSresult@permutation.object[[1]])!=B)){ 
		
		permutation.object <- list()
		
		
		###### APPLYING ANALYSIS TO DATA ########
		
		
		## MFA Analysis ##
		if(type=="CSmfa"){
			
			## Some MFA pre-amble ##
			
			#mfa.factor cannot be NULL (correction, choose the one with highest ref loadings, but give warning) + check if a new addition to the CS object will have to be made or not
			if(is.null(mfa.factor)){
				ref.loadings <- apply(CSresult@extra$object$quanti.var$cor,MARGIN=2,FUN=function(x){mean(x[ref.index])})
				mfa.factor <- which(abs(ref.loadings)==max(abs(ref.loadings)))
				# cat(paste0("Since mfa.factor was not given, Factor ",mfa.factor," has been chosen based on the CSresult object.\n\n"))
				warning(paste0("Since mfa.factor was not given, Factor ",mfa.factor," has been chosen based on the CSresult object."))
			}
			
			if(!(mfa.factor %in% CSresult@call$component.select)){
				CS.extrafactor <- TRUE
			}else{
				CS.extrafactor <- FALSE
			}
			
			
			##

			if(verbose){cat("Analysing Permuted Data with MFA\n(Factor is chosen based on highest average reference loadings):\n")}
			
			CS.Perm <- vector("list",B)
			MFA.factor.Perm <- c(1:B)
			
			
			# TO DO: Add parallellisation option (B/nCPU models per core)
			
			
			if(MultiCores){
			  message("Parallelisation: ",MultiCores.number," cores used")
			  current_environment <- environment()
			  
			  worker_divide <- splitIndices(B,MultiCores.number)
			  gentype <- "RNGstream"
			  
			  
			  
			  cl <- do.call("makeClusterFT", c(list(min(MultiCores.number, length(worker_divide)), 
			                                          "SOCK", ft_verbose = FALSE), NULL))


			  
			  clusterCall(cl, eval, substitute(library(FactoMineR)), env = .GlobalEnv)
			  clusterExport(cl, c("worker_divide","querMat","refMat","CSresult","ref.index"), envir=current_environment)
			  
			  clusterSetupRNG.FT(cl, type = gentype, streamper = "replicate", 
			                       seed = MultiCores.seed, n = length(worker_divide), prngkind = "default")
			  

			  
			  
			  res <- clusterApplyFT(cl, x=1:length(worker_divide),
			                        fun=function(x){
			                          
			                          CS.Perm <- vector("list",length(worker_divide[[x]]))
			                          MFA.factor.Perm <- c(1:length(worker_divide[[x]])) 
			                          
			                          for(i in 1:length(worker_divide[[x]])){
			                            ## Permuting Query Matrix
			                            querMat.Perm <- matrix(sample(as.vector(querMat),(dim(querMat)[1]*dim(querMat)[2]),replace=FALSE),nrow=dim(querMat)[1],ncol=dim(querMat)[2])
			                            
			                            rownames(querMat.Perm) <- rownames(refMat)
			                            colnames(querMat.Perm) <- colnames(querMat)
			                            
			                            #########################
			                            
			                            data_comb <- cbind(refMat,querMat.Perm)
			                            rownames(data_comb) <- rownames(refMat)
			                            colnames(data_comb) <- c(colnames(refMat),colnames(querMat))
			                            
			                            out_MFA <- MFA(base=data_comb,group=c(ncol(refMat),ncol(querMat)),type=c("s","s"),ind.sup=NULL,row.w=CSresult@call$analysis.pm$row.w,weight.col.mfa=CSresult@call$analysis.pm$weight.col.mfa,ncp=CSresult@call$analysis.pm$ncp,name.group=c("Reference","Query"),graph=FALSE)
			                            
			                            ## FACTOR IS CHOSEN BASED ON HIGHEST AVERAGE REFERENCE LOADINGS  (put other options here later)
			                            ref.loadings <- apply(out_MFA$quanti.var$cor,MARGIN=2,FUN=function(x){mean(x[ref.index])})
			                            factor.choose <- which(abs(ref.loadings)==max(abs(ref.loadings)))
			                            
			                            
			                            MFA.factor.Perm[i] <- factor.choose
			                            CS.Perm[[i]] <- out_MFA$quanti.var$cor[,factor.choose]
			                            
			                          }
			                          return(list(CS.Perm=CS.Perm,MFA.factor.Perm=MFA.factor.Perm))
			                          
			                          
			                        }
			                        ,gentype = gentype, seed = MultiCores.seed, prngkind = "default",mngtfiles=c("","",""))[[1]]
			  
			  stopClusterFT(cl)
			  
			  # res <- performParallel(MultiCores.number,1:length(worker_divide),
			  #                        seed=MultiCores.seed,cltype="SOCK",
			  #                        initexpr=library(FactoMineR),
			  #                        export=c("worker_divide","querMat","refMat","CSresult","ref.index"),
			  #                        
			  #                        
			  #                        fun=function(x){
			  #   
			  #  CS.Perm <- vector("list",length(worker_divide[[x]]))
			  #  MFA.factor.Perm <- c(1:length(worker_divide[[x]])) 
			  #  
			  #  for(i in 1:length(worker_divide[[x]])){
			  #    ## Permuting Query Matrix
			  #    querMat.Perm <- matrix(sample(as.vector(querMat),(dim(querMat)[1]*dim(querMat)[2]),replace=FALSE),nrow=dim(querMat)[1],ncol=dim(querMat)[2])
			  #    
			  #    rownames(querMat.Perm) <- rownames(refMat)
			  #    colnames(querMat.Perm) <- colnames(querMat)
			  #    
			  #    #########################
			  #    
			  #    data_comb <- cbind(refMat,querMat.Perm)
			  #    rownames(data_comb) <- rownames(refMat)
			  #    colnames(data_comb) <- c(colnames(refMat),colnames(querMat))
			  #    
			  #    out_MFA <- MFA(base=data_comb,group=c(ncol(refMat),ncol(querMat)),type=c("s","s"),ind.sup=NULL,row.w=CSresult@call$analysis.pm$row.w,weight.col.mfa=CSresult@call$analysis.pm$weight.col.mfa,ncp=CSresult@call$analysis.pm$ncp,name.group=c("Reference","Query"),graph=FALSE)
			  #    
			  #    ## FACTOR IS CHOSEN BASED ON HIGHEST AVERAGE REFERENCE LOADINGS  (put other options here later)
			  #    ref.loadings <- apply(out_MFA$quanti.var$cor,MARGIN=2,FUN=function(x){mean(x[ref.index])})
			  #    factor.choose <- which(abs(ref.loadings)==max(abs(ref.loadings)))
			  #    
			  #    
			  #    MFA.factor.Perm[i] <- factor.choose
			  #    CS.Perm[[i]] <- out_MFA$quanti.var$cor[,factor.choose]
			  #    
			  #  }
			  #  return(list(CS.Perm=CS.Perm,MFA.factor.Perm=MFA.factor.Perm))
			  #  
			  #   
			  # })
			  
			  MFA.factor.Perm <- unlist(lapply(res,FUN=function(x){x$MFA.factor.Perm}))
			  CS.Perm <- do.call(c,lapply(res,FUN=function(x){x$CS.Perm}))
			  
			}else{
			  for(i in 1:B){
			    
			    ## Permuting Query Matrix
			    querMat.Perm <- matrix(sample(as.vector(querMat),(dim(querMat)[1]*dim(querMat)[2]),replace=FALSE),nrow=dim(querMat)[1],ncol=dim(querMat)[2])
			    
			    rownames(querMat.Perm) <- rownames(refMat)
			    colnames(querMat.Perm) <- colnames(querMat)
			    
			    #########################
			    
			    data_comb <- cbind(refMat,querMat.Perm)
			    rownames(data_comb) <- rownames(refMat)
			    colnames(data_comb) <- c(colnames(refMat),colnames(querMat))
			    
			    out_MFA <- MFA(base=data_comb,group=c(ncol(refMat),ncol(querMat)),type=c("s","s"),ind.sup=NULL,row.w=CSresult@call$analysis.pm$row.w,weight.col.mfa=CSresult@call$analysis.pm$weight.col.mfa,ncp=CSresult@call$analysis.pm$ncp,name.group=c("Reference","Query"),graph=FALSE)
			    
			    ## FACTOR IS CHOSEN BASED ON HIGHEST AVERAGE REFERENCE LOADINGS  (put other options here later)
			    ref.loadings <- apply(out_MFA$quanti.var$cor,MARGIN=2,FUN=function(x){mean(x[ref.index])})
			    factor.choose <- which(abs(ref.loadings)==max(abs(ref.loadings)))
			    
			    
			    MFA.factor.Perm[i] <- factor.choose
			    CS.Perm[[i]] <- out_MFA$quanti.var$cor[,factor.choose]
			    
			    if(verbose){
			      cat(".")
			      if(i%%100==0){cat(" ",i,"\n")}
			    }
			  }
			  if(verbose){cat("\nDONE\n")}
			}
			
			

			
			permutation.object$CLoadings.Perm <- CS.Perm
			permutation.object$ChosenFactor.Perm <- MFA.factor.Perm
			
			CS.Perm.rank <- lapply(CS.Perm,FUN=function(x){
						return(CSrank2(as.data.frame(x),ref.index=ref.index,plot=FALSE,component.plot=1)$CSRankScores)
					})
			permutation.object$CRankScores.Perm <- CS.Perm.rank
		}
		
		## Zhang Analysis ##
		if(type=="CSzhang"){
			if(verbose){cat("Analysing Permuted Data with Zhang and Gant:\n")}
			
			CS.Perm <- vector("list",B)
			
			if(MultiCores){
			  message("Parallelisation: ",MultiCores.number," cores used")
			  
			  worker_divide <- splitIndices(B,MultiCores.number)
			  
			  
			  
			  current_environment <- environment()
			  

			  gentype <- "RNGstream"
			  
			  cl <- do.call("makeClusterFT", c(list(min(MultiCores.number, length(worker_divide)), 
			                                        "SOCK", ft_verbose = FALSE), NULL))
			  
			  clusterCall(cl, eval, substitute(library(CSFA)), env = .GlobalEnv)

			  clusterExport(cl, c("worker_divide","querMat","refMat","CSresult","ref.index"), envir=current_environment)
			  clusterSetupRNG.FT(cl, type = gentype, streamper = "replicate", 
			                     seed = MultiCores.seed, n = length(worker_divide), prngkind = "default")
			  
			  
			  res <- clusterApplyFT(cl, x=1:length(worker_divide),
			                        fun=function(x){
			                          
			                          CS.Perm <- vector("list",length(worker_divide[[x]]))

			                          for(i in 1:length(worker_divide[[x]])){
			                            ## Permuting Query Matrix
			                            querMat.Perm <- matrix(sample(as.vector(querMat),(dim(querMat)[1]*dim(querMat)[2]),replace=FALSE),nrow=dim(querMat)[1],ncol=dim(querMat)[2])
			                            ##########################
			                            
			                            rownames(querMat.Perm) <- rownames(refMat)
			                            colnames(querMat.Perm) <- colnames(querMat)
			                            
			                            out_zhang <- CSFA:::analyse_zhang(refMat,querMat.Perm,
			                                                              nref=CSresult@call$analysis.pm$nref,nquery=CSresult@call$analysis.pm$nquery,ord.query=CSresult@call$analysis.pm$ord.query,ntop.scores=CSresult@call$analysis.pm$ntop.scores,
			                                                              basefilename="analyseZhang",which=c(),plot.type="device",print.top=FALSE)
			                            
			                            CS.Perm[[i]] <- out_zhang$All[,1]
			                            names(CS.Perm[[i]]) <- rownames(out_zhang$All)
			                            
			                          }
			                          return(CS.Perm)
			                          
			                          
			                        }
			                        ,gentype = gentype, seed = MultiCores.seed, prngkind = "default",mngtfiles=c("","",""))[[1]]
			  
			  stopClusterFT(cl)
			  
			  CS.Perm <- do.call(c,res)
			  
			}else{
			  
			  for(i in 1:B){
			    
			    ## Permuting Query Matrix
			    querMat.Perm <- matrix(sample(as.vector(querMat),(dim(querMat)[1]*dim(querMat)[2]),replace=FALSE),nrow=dim(querMat)[1],ncol=dim(querMat)[2])
			    ##########################
			    
			    rownames(querMat.Perm) <- rownames(refMat)
			    colnames(querMat.Perm) <- colnames(querMat)
			    
			    out_zhang <- analyse_zhang(refMat,querMat.Perm,
			                               nref=CSresult@call$analysis.pm$nref,nquery=CSresult@call$analysis.pm$nquery,ord.query=CSresult@call$analysis.pm$ord.query,ntop.scores=CSresult@call$analysis.pm$ntop.scores,
			                               basefilename="analyseZhang",which=c(),plot.type="device",print.top=FALSE)
			    
			    CS.Perm[[i]] <- out_zhang$All[,1]
			    names(CS.Perm[[i]]) <- rownames(out_zhang$All)
			    
			    if(verbose){
			      cat(".")
			      if(i%%100==0){cat(" ",i,"\n")}
			    }
			  }
			  
			  if(verbose){cat("\nDONE\n")}
			  
			}	
			

			
			permutation.object$ZGscore <- CS.Perm
			CS.Perm.rank <- NULL
			
			mfa.factor <- NULL
			CS.extrafactor <- FALSE
		}
	}
	else{
		permutation.object <- CSresult@permutation.object # If the previous part was not necessary to do, extract the existing permutation.objects
	
		if(type=="CSzhang"){
			CS.Perm <- CSresult@permutation.object$ZGscore
			CS.Perm.rank <- NULL
		}else{
			CS.Perm <- CSresult@permutation.object$CLoadings.Perm
			CS.Perm.rank <- CSresult@permutation.object$CRankScores.Perm
			mfa.factor <- CSresult@permutation.object$extra.parameters$mfa.factor
		
			CS.extrafactor <- FALSE
			
#			if(!(mfa.factor %in% CSresult@call$component.select)){
#				CS.extrafactor <- TRUE
#			}else{
#				CS.extrafactor <- FALSE
#			}
		}
		
		
	}
	
	
	##### COMPUTING THE P-VALUES  + FINISHING PERMUTATION.OBJECT #####
	
	# Note: pvalues will always be re-computed to allow for different mfa.factor
	
#	# Getting Pvalues from CS and CSrank
	pval.dataframe.temp <- pvalue_compute(obs.result=CSresult,ref.index=ref.index,list.h0.result=CS.Perm,list.h0.result.rank=CS.Perm.rank,mfa.factor=mfa.factor,method.adjust=method.adjust)
	pval.dataframe <- pval.dataframe.temp[[1]]
	pval.dataframe.rank <- pval.dataframe.temp[[2]]
	
	permutation.object$CS.pval.dataframe <- pval.dataframe
	permutation.object$CSRank.pval.dataframe <- pval.dataframe.rank
	
	
	CSresult@permutation.object <- permutation.object # Saving the updated permutation.object in CSresult
	
	
	#### UPDATE THE CSresult object which will be returned in the end
	if(type=="CSzhang"){
		
		## Adding p-values to CS slot
		
		CS <- CSresult@CS
		
		CS$CS.query$pvalues <- pval.dataframe$pvalues
		col.order <- c("ZGscore","pvalues","ZGrank","ZGabsrank")
		if(method.adjust!="none"){
			CS$CS.query$pvalues.adjusted <- pval.dataframe$pvalues.adjusted
			col.order <- c("ZGscore","pvalues","pvalues.adjusted","ZGrank","ZGabsrank")
		}
		CS$CS.query <- CS$CS.query[,col.order]
		
		CSresult@CS <- CS
				
	}
	if(type=="CSmfa"){
		
		CS <- CSresult@CS
		
		if(CS.extrafactor){
			warning("Due to choice of mfa.factor, another factor is added to the CS and GS slot.")
			
			loadings <- CSresult@extra$object$quanti.var$coord[,mfa.factor]
			scores <- CSresult@extra$object$ind$coord[,mfa.factor]	
			CSRank <- pval.dataframe.rank$observed
			
			if(method.adjust!="none"){
				CS.query.temp <- data.frame(
						CLoadings=loadings[-c(1:length(ref.index))],
						CLpvalues=pval.dataframe$pvalues,
						CLpvalues.adjusted=pval.dataframe$pvalues.adjusted,
						CRankScores=CSRank,
						CRpvalues=pval.dataframe.rank$pvalues,
						CRpvalues.adjusted=pval.dataframe.rank$pvalues.adjusted,
						CLrank=as.integer(rank(-loadings[-c(1:length(ref.index))])),
						CLabsrank=as.integer(rank(-abs(loadings[-c(1:length(ref.index))]))),
						CRrank=as.integer(rank(-CSRank)),
						CRabsrank=as.integer(rank(-abs(CSRank))),
						row.names=rownames(loadings)[-c(1:length(ref.index))]
				)
				
			}else{
				CS.query.temp <- data.frame(
						CLoadings=loadings[-c(1:length(ref.index))],
						CLpvalues=pval.dataframe$pvalues,
						CRankScores=CSRank,
						CRpvalues=pval.dataframe.rank$pvalues,
						CLrank=as.integer(rank(-loadings[-c(1:length(ref.index))])),
						CLabsrank=as.integer(rank(-abs(loadings[-c(1:length(ref.index))]))),
						CRrank=as.integer(rank(-CSRank)),
						CRabsrank=as.integer(rank(-abs(CSRank))),
						row.names=rownames(loadings)[-c(1:length(ref.index))]
				)
				
			}
			
			CS[[length(CS)+1]] <- list(
					CS.ref=data.frame(CLoadings=loadings[c(1:length(ref.index))],row.names=rownames(loadings)[c(1:length(ref.index))]),
					CS.query=CS.query.temp
				)
			
			names(CS)[length(CS)] <- paste0("Factor",mfa.factor)
			
			GS$extrafactor <- scores
			colnames(GS)[dim(GS)[2]] <- paste0("Factor",mfa.factor)
			
			CSresult@CS <- CS
			CSresult@GS <- GS
			
			
		}else{
			factor.index <- which(paste0("Factor",mfa.factor)==names(CS))
			CS[[factor.index]]$CS.query$CLpvalues <- pval.dataframe$pvalues
			CS[[factor.index]]$CS.query$CRpvalues <- pval.dataframe.rank$pvalues
			col.order <- c("CLoadings","CLpvalues","CRankScores","CRpvalues","CLrank","CLabsrank","CRrank","CRabsrank")
			
			if(method.adjust!="none"){
				CS[[factor.index]]$CS.query$CLpvalues.adjusted <- pval.dataframe$pvalues.adjusted
				CS[[factor.index]]$CS.query$CRpvalues.adjusted <- pval.dataframe.rank$pvalues.adjusted
				col.order <- c("CLoadings","CLpvalues","CLpvalues.adjusted","CRankScores","CRpvalues","CRpvalues.adjusted","CLrank","CLabsrank","CRrank","CRabsrank")
			}
			CS[[factor.index]]$CS.query <- CS[[factor.index]]$CS.query[,col.order]
			
			CSresult@CS <- CS
		}
		
		
	}
	
	
	##### POSSIBLE PLOTS #####
	
	# NOTE: Still need to check all cases with 2 plot.types 
	
	hist.drawn <- FALSE
	
	# Volcano plot
	if(1 %in% which){
		
		if((2 %in% which) & (is.null(cmpd.hist))){
			# Volcano plot with selecting compounds
			pvalue_volc(pval.dataframe=permutation.object$CS.pval.dataframe,type=type,color.columns=color.columns,list.h0.result=CS.Perm,list.h0.result.rank=CS.Perm.rank,make.hist=TRUE,plot.type=plot.type,plot.type.hist=plot.type,basefilename=paste0(basefilename,"_volcanoplot"))
			
			hist.drawn <- TRUE
		}else{
			# Volcano plots without selecting compounds
			pvalue_volc(pval.dataframe=permutation.object$CS.pval.dataframe,type=type,color.columns=color.columns,list.h0.result=NULL,list.h0.result.rank=NULL,make.hist=FALSE,plot.type=plot.type,plot.type.hist=plot.type,basefilename=paste0(basefilename,"_volcanoplot"))
		}
		
	}
	
	
	# Histogram plot
	if((2 %in% which) & !hist.drawn){
		if(is.null(cmpd.hist)){
			#volc with plot device
			pvalue_volc(pval.dataframe=permutation.object$CS.pval.dataframe,type=type,color.columns=color.columns,list.h0.result=CS.Perm,list.h0.result.rank=CS.Perm.rank,make.hist=TRUE,plot.type="device",plot.type.hist=plot.type,basefilename=paste0(basefilename))
			
		}else{
			# just hist
			pvalue_hist(pval.dataframe=permutation.object$CS.pval.dataframe[cmpd.hist,],list.h0.result=CS.Perm,list.h0.result.rank=CS.Perm.rank,type=type,plot.type=plot.type,basefilename=paste0(basefilename,"_histogram"))
			
		}
		
		hist.drawn <- TRUE
	}
	
	## SAME PLOTS FOR CSRANKSCORES
	
	hist.drawn.rank <- FALSE
	
	if(!is.null(permutation.object$CSRank.pval.dataframe)){ #zhang does not have this
		
		# Volcano plot
		if(3 %in% which){
		
			if((4 %in% which) & (is.null(cmpd.hist))){
				# Volcano plot with selecting compounds
				pvalue_volc(pval.dataframe=permutation.object$CSRank.pval.dataframe,CSRank=TRUE,type=type,color.columns=color.columns,list.h0.result=CS.Perm,list.h0.result.rank=CS.Perm.rank,make.hist=TRUE,plot.type=plot.type,plot.type.hist=plot.type,basefilename=paste0(basefilename,"_volcanoplot_CSRankScores"))
			
				hist.drawn.rank <- TRUE
			}else{
				# Volcano plots without selecting compounds
				pvalue_volc(pval.dataframe=permutation.object$CSRank.pval.dataframe,CSRank=TRUE,type=type,color.columns=color.columns,list.h0.result=NULL,list.h0.result.rank=NULL,make.hist=FALSE,plot.type=plot.type,plot.type.hist=plot.type,basefilename=paste0(basefilename,"_volcanoplot_CSRankScores"))
			
			}
		
		}
	
	
		# Histogram plot
		if((4 %in% which) & !hist.drawn.rank){
			if(is.null(cmpd.hist)){
				#volc with plot device
				pvalue_volc(pval.dataframe=permutation.object$CSRank.pval.dataframe,CSRank=TRUE,type=type,color.columns=color.columns,list.h0.result=CS.Perm,list.h0.result.rank=CS.Perm.rank,make.hist=TRUE,plot.type="device",plot.type.hist=plot.type,basefilename=paste0(basefilename,"_CSRankScores"))
			
			}else{
				# just hist
				pvalue_hist(pval.dataframe=permutation.object$CSRank.pval.dataframe[cmpd.hist,],CSRank=TRUE,list.h0.result=CS.Perm,list.h0.result.rank=CS.Perm.rank,type=type,plot.type=plot.type,basefilename=paste0(basefilename,"_histogram_CSRankScores"))
			
			}
		
			hist.drawn.rank <- TRUE
		}
	
	}
	
	##### RETURN OBJECT ######
	

	
	# add p-value information + for which factor it was done
	CSresult@permutation.object$extra.parameters <- list(mfa.factor=mfa.factor,method.adjust=method.adjust) #change to component.select
	
	return(CSresult)
	# Contains:
	# - updated CS with pvalues 
	# - permutation.object slot with permuted CS, pval.dataframe and if asked the permuted data.
	
	
}



# add adjusted pvalues
pvalue_compute <- function(obs.result,list.h0.result,list.h0.result.rank,ref.index=1,mfa.factor=1,method.adjust=c("none","holm","hochberg","hommel","bonferroni","BH","BY","fdr")){
	if(class(obs.result)!="CSresult"){stop("obs.result is not a \"CSresult\" class object")}
	
	type <- obs.result@type
	
	if(!(type %in% c("CSzhang","CSmfa"))){stop("Only Zhang and MFA results can be used")}
	
	if(type=="CSzhang"){
		
		obs.scores <- obs.result@CS$CS.query$ZGscore
		
		pval.dataframe <- data.frame(Cmpd=rownames(obs.result@CS$CS.query))
		pval.dataframe$Cmpd <- as.character(pval.dataframe$Cmpd)
		temp.pval <- c(1:length(obs.scores))
		
		for(i.cmpd in c(1:length(obs.scores))){
			h0.data <- unlist(lapply(list.h0.result,FUN=function(x){return(x[i.cmpd])}))
			obs <- obs.scores[i.cmpd]
			pvalue <- (1+sum(abs(h0.data)>=abs(obs)))/(length(h0.data)+1)
			
			temp.pval[i.cmpd] <- pvalue
		}
		
		pval.dataframe$observed <- obs.scores
		pval.dataframe$pvalues <- temp.pval
		if(method.adjust!="none"){pval.dataframe$pvalues.adjusted <- p.adjust(temp.pval,method=method.adjust)}
		
		
		return(list(pval.dataframe,NULL))
	}
	
	if(type=="CSmfa"){
		
		
		# Prep for CS pvalues
		obs.scores <- obs.result@extra$object$quanti.var$coord[-ref.index,mfa.factor]
		pval.dataframe <- data.frame(Cmpd=names(obs.scores))
		pval.dataframe$Cmpd <- as.character(pval.dataframe$Cmpd)
		temp.pval <- c(1:length(obs.scores))
		
		# Prep for CSrank pvalues
		obs.scores.rank <- CSrank2(obs.result@extra$object$quanti.var$coord,ref.index=ref.index,plot=FALSE,component.plot=mfa.factor)$CSRankScores
		pval.dataframe.rank<- data.frame(Cmpd=names(obs.scores))
		temp.pval.rank <- c(1:length(obs.scores))
		
		for(i.cmpd in c(1:length(obs.scores))){
			
			## P-values of CS
			h0.data <- unlist(lapply(list.h0.result,FUN=function(x){return(x[i.cmpd+length(ref.index)])}))
			obs <- obs.scores[i.cmpd]
			pvalue <- (1+sum(abs(h0.data)>=abs(obs)))/(length(h0.data)+1)
			temp.pval[i.cmpd] <- pvalue
			
			# P-values of CSRank
			h0.data.rank <- unlist(lapply(list.h0.result.rank,FUN=function(x){return(x[i.cmpd])}))
			obs.rank <- obs.scores.rank[i.cmpd]
			pvalue.rank <- (1+sum(abs(h0.data.rank)>=abs(obs.rank)))/(length(h0.data.rank)+1)
			temp.pval.rank[i.cmpd] <- pvalue.rank
						
		}
		pval.dataframe$pvalues <- temp.pval
		if(method.adjust!="none"){pval.dataframe$pvalues.adjusted <- p.adjust(temp.pval,method=method.adjust)}
		pval.dataframe$observed <- obs.scores
		
		pval.dataframe.rank$pvalues <- temp.pval.rank
		if(method.adjust!="none"){pval.dataframe.rank$pvalues.adjusted <- p.adjust(temp.pval.rank,method=method.adjust)}
		pval.dataframe.rank$observed <- obs.scores.rank
		
		return(list(pval.dataframe,pval.dataframe.rank))
		
	}
}

# function for all pvalues and h0 given
pvalue_hist <- function(pval.dataframe,list.h0.result,list.h0.result.rank,type=c("CSmfa","CSzhang"),plot.type="device",basefilename="pvalue_hist",CSRank=FALSE){
	
	plot.in <- function(plot.type,name){
		if(plot.type=="pdf"){pdf(paste0(name,".pdf"))}
		if(plot.type=="device"){dev.new()}
		if(plot.type=="sweave"){}
	}
	plot.out <- function(plot.type){if(plot.type=="pdf"){dev.off()}}
	
	
	cmpd.names <- pval.dataframe$Cmpd
	obs.scores <- pval.dataframe$observed
	
	
	if(type=="CSzhang"){
		
		for(i in 1:dim(pval.dataframe)[1]){
			
			h0.data <- unlist(lapply(list.h0.result,FUN=function(x){return(x[ which(names(list.h0.result[[1]])==cmpd.names[i])   ])}))
			
			
			obs <- obs.scores[i]
			pvalue <- pval.dataframe$pvalues[i]
			max.y <- max(hist(h0.data,nclass=15,plot=FALSE)$counts)
			
			plot.in(plot.type,paste0(basefilename,"_cmpd_",cmpd.names[i]))
			hist(h0.data,xlim=c(min(c(h0.data,obs,-obs)),max(c(h0.data,obs,-obs))),nclass=15,col='lightblue',main=paste0('Distribution of Cmpd ',cmpd.names[i]," under H0"),xlab="Zhang Score")
			abline(v=obs,lw=2,col="red")
			abline(v=-obs,lw=2,lty=2,col="red")
			if(obs>=0){text(obs,0.9*max.y, paste0("Observed Value (p-value=",round(pvalue,digits=4),")"),pos=2,col="red",offset=1)} else{text(obs,0.9*max.y, paste0("Observed Value (p-value=",round(pvalue,digits=4),")"),pos=4,col="red",offset=1)}
			
			if("pvalues.adjusted"%in%names(pval.dataframe)){
				if(obs>=0){text(obs,0.75*max.y, paste0("(adjusted p-value=",round(pval.dataframe$pvalues.adjusted[i],digits=4),")"),pos=2,col="red",offset=1)} else{text(obs,0.75*max.y, paste0("(adjusted p-value=",round(pval.dataframe$pvalues.adjusted[i],digits=4),")"),pos=4,col="red",offset=1)}
				
			}
			
			plot.out(plot.type)
		}
		
	}
	
	if(type=="CSmfa"){
		
		
		for(i in 1:dim(pval.dataframe)[1]){
			obs <- obs.scores[i]
			
			if(CSRank){
				len.ref <- length(list.h0.result[[1]]) - length(list.h0.result.rank[[1]])
				
				i.temp <- which(names(list.h0.result[[1]])==cmpd.names[i]) - len.ref
				
				h0.data.rank <- unlist(lapply(list.h0.result.rank,FUN=function(x){return(x[i.temp])}))
				
				
				pvalue <- pval.dataframe$pvalues[i]
				max.y <- max(hist(h0.data.rank,nclass=15,plot=FALSE)$counts)
				
				
				plot.in(plot.type,paste0(basefilename,"_cmpd_",cmpd.names[i]))
				hist(h0.data.rank,xlim=c(min(c(h0.data.rank,obs,-obs)),max(c(h0.data.rank,obs,-obs))),nclass=15,col='lightblue',main=paste0('Distribution of Cmpd ',cmpd.names[i]," under H0"),xlab="MFA CSRankscore")
				abline(v=obs,lw=2,col="red")
				abline(v=-obs,lw=2,lty=2,col="red")
				if(obs>=0){text(obs,0.9*max.y, paste0("Observed Value"),pos=2,col="red",offset=1)} else{text(obs,0.9*max.y, paste0("Observed Value"),pos=4,col="red",offset=1)}
				text(obs,0.75*max.y,paste0("P-Value = ",round(pvalue,digits=4)),col="red",pos=ifelse(obs>=0,2,4))
				if("pvalues.adjusted"%in%names(pval.dataframe)){
					if(obs>=0){text(obs,0.6*max.y, paste0("Adjusted P-value = ",round(pval.dataframe$pvalues.adjusted[i],digits=4)),pos=2,col="red",offset=1)} else{text(obs,0.6*max.y, paste0("Adjusted P-value = ",round(pval.dataframe$pvalues.adjusted[i],digits=4)),pos=4,col="red",offset=1)}
					
				}
				plot.out(plot.type)
			}
			else{
				
				h0.data <- unlist(lapply(list.h0.result,FUN=function(x){return(x[cmpd.names[i]])}))
				

				pvalue <- pval.dataframe$pvalues[i]
				max.y <- max(hist(h0.data,nclass=15,plot=FALSE)$counts)
				
				
				plot.in(plot.type,paste0(basefilename,"_cmpd_",cmpd.names[i]))
				hist(h0.data,xlim=c(min(c(h0.data,obs,-obs)),max(c(h0.data,obs,-obs))),nclass=15,col='lightblue',main=paste0('Distribution of Cmpd ',cmpd.names[i]," under H0"),xlab="MFA CScore")
				abline(v=obs,lw=2,col="red")
				abline(v=-obs,lw=2,col="red",lty=2)
				
				if(obs>=0){text(obs,0.9*max.y, paste0("Observed Value (p-value=",round(pvalue,digits=4),")"),pos=2,col="red",offset=1)} else{text(obs,0.9*max.y, paste0("Observed Value (p-value=",round(pvalue,digits=4),")"),pos=4,col="red",offset=1)}
				
				if("pvalues.adjusted"%in%names(pval.dataframe)){
					if(obs>=0){text(obs,0.75*max.y, paste0("Adjusted P-value = ",round(pval.dataframe$pvalues.adjusted[i],digits=4)),pos=2,col="red",offset=1)} else{text(obs,0.75*max.y, paste0("Adjusted P-value = ",round(pval.dataframe$pvalues.adjusted[i],digits=4)),pos=4,col="red",offset=1)}
					
				}
				plot.out(plot.type)
			}
			
		}
		
	}
	
	
}


# Use compute pvalue plot
pvalue_volc <- function(pval.dataframe,type=c("CSmfa","CSzhang"),CSRank=FALSE,color.columns=NULL,list.h0.result=NULL,list.h0.result.rank=NULL,make.hist=FALSE,plot.type="device",plot.type.hist="device",basefilename="pvalue_volc"){
	
	plot.in <- function(plot.type,name){
		if(plot.type=="pdf"){pdf(paste0(name,".pdf"))}
		if(plot.type=="device"){dev.new()}
		if(plot.type=="sweave"){}
	}
	plot.out <- function(plot.type){if(plot.type=="pdf"){dev.off()}}
	
	if("pvalues.adjusted" %in% names(pval.dataframe)){use.adjusted <- TRUE}else{use.adjusted <- FALSE}
	
	if(use.adjusted){
		pvalues <- pval.dataframe$pvalues.adjusted
	}else{
		pvalues <- pval.dataframe$pvalues
	}
		
	obs.scores <- pval.dataframe$observed
	temp.names <- pval.dataframe$Cmpd
	
	plot.pvalues <- -log(pvalues)

	
	if(is.null(color.columns)){color.columns <- "black"}
	ylab.temp <- ifelse(use.adjusted,"-log(adjusted pvalues)","-log(pvalues)")
	main.temp <- ifelse(CSRank,"CSRankScores","CScores")
	xlab.temp <- ifelse(CSRank,"Observed CSRank for Cmpds","Observed CS for Cmpds")
	
	plot.in(plot.type,basefilename)
	plot(obs.scores,plot.pvalues,col=color.columns,bg="grey",pch=21,xlab=xlab.temp,ylab=ylab.temp,main=paste0("Volcano Plot for ",type," result - ",main.temp))
	text(obs.scores,plot.pvalues,temp.names,pos=1,col=color.columns)
	abline(h=-log(0.05),col="red",lty=3)
	text(min(obs.scores),-log(0.05)+0.1,paste0(ifelse(use.adjusted,"adjusted p-value","p-value")," = 0.05"),pos=4,col="red")
	plot.out(plot.type)
	
	
	if(make.hist){
		if(!is.null(list.h0.result)){
			
			if(plot.type!="device"){
				dev.new()
				plot(obs.scores,plot.pvalues,col=color.columns,bg="grey",pch=21,xlab="Observed CS for cmpds",ylab="-log(pvalues)",main=paste0("Volcano Plot for ",type," result"))
				text(obs.scores,plot.pvalues,temp.names,pos=1,col=color.columns)
			}
			
			cat("Please choose one or more compounds with left-click.\n To end the selection, press right-click.")
			choose.cmpd <- identify(obs.scores,plot.pvalues,n=9999,labels=temp.names,col="slateblue3") 
			if(!CSRank){choose.cmpd <- choose.cmpd }

			pvalue_hist(pval.dataframe[choose.cmpd,],list.h0.result,list.h0.result.rank,type,plot.type=plot.type.hist,basefilename=paste(basefilename,"_histogram"),CSRank=CSRank)
			
		}
	}
}

# add this to the CScompare with a pvalcompare option
# special case for only 2 pval.dataframe and when transforming already happened before
pvalue2_compare <- function(list.pval.dataframe,threshold=0.05){
	if(length(list.pval.dataframe)!=2){stop("Need exactly 2 pval.dataframe results")}
	
	pvalues <- lapply(list.pval.dataframe,FUN=function(x){x$pvalues})
	m.pvalues <- matrix(unlist(pvalues),ncol=length(pvalues))
	compare.val <- apply(m.pvalues,MARGIN=1,FUN=function(x){paste0("S",c(1:length(pvalues))[x<=threshold],collapse="")})
	compare.val <- sapply(compare.val,FUN=function(x){if(x=="S"){return("NS")}else{return(x)}})
	
	table.compare.val <- table(compare.val)
	
	
	m11 <- ifelse(!is.na(table.compare.val["S1S2"]),table.compare.val["S1S2"],0)
	m21 <- ifelse(!is.na(table.compare.val["S1"]),table.compare.val["S1"],0)
	m12 <- ifelse(!is.na(table.compare.val["S2"]),table.compare.val["S2"],0)
	m22 <- ifelse(!is.na(table.compare.val["NS"]),table.compare.val["NS"],0)
	
	pval.matrix <- matrix(c(m11,m21,m12,m22),nrow=2)
	colnames(pval.matrix) <- c("Result1.Sign","Result1.NotSign")
	rownames(pval.matrix) <- c("Result2.Sign","Result2.NotSign")
		
	# adjusted pvalues
	adjusted.available <- unlist(lapply(list.pval.dataframe,FUN=function(x){"pvalues.adjusted"%in%colnames(x)}))
	if(sum(adjusted.available)==1){warning("Only 1 of the results has adjusted p-values.",call.=FALSE)}
	if(sum(adjusted.available)==2){
		
		pvalues.adj <- lapply(list.pval.dataframe,FUN=function(x){x$pvalues.adjusted})
		m.pvalues.adj <- matrix(unlist(pvalues.adj),ncol=length(pvalues.adj))
		compare.val.adj <- apply(m.pvalues.adj,MARGIN=1,FUN=function(x){paste0("S",c(1:length(pvalues.adj))[x<=threshold],collapse="")})
		compare.val.adj <- sapply(compare.val.adj,FUN=function(x){if(x=="S"){return("NS")}else{return(x)}})
		
		table.compare.val.adj <- table(compare.val.adj)
		
		
		m11 <- ifelse(!is.na(table.compare.val.adj["S1S2"]),table.compare.val.adj["S1S2"],0)
		m21 <- ifelse(!is.na(table.compare.val.adj["S1"]),table.compare.val.adj["S1"],0)
		m12 <- ifelse(!is.na(table.compare.val.adj["S2"]),table.compare.val.adj["S2"],0)
		m22 <- ifelse(!is.na(table.compare.val.adj["NS"]),table.compare.val.adj["NS"],0)
		
		pval.adj.matrix <- matrix(c(m11,m21,m12,m22),nrow=2)
		colnames(pval.adj.matrix) <- c("Result1.Sign","Result1.NotSign")
		rownames(pval.adj.matrix) <- c("Result2.Sign","Result2.NotSign")
		
	}
	else{pval.adj.matrix <- NULL}
		
		
	# return 1 (or 2) matrices and the new list.pval.dataframes
	return(list(pvalues=pval.matrix,adj.pvalues=pval.adj.matrix))
	
	
}


