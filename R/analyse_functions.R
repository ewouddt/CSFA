# Project: CSFA
# 
# Author: lucp8394
###############################################################################


#### GENERAL FUNCTION TO DO FA MODEL ####



analyse_FA <- function(
		# GENERAL
		data,modeltype="CSmfa",matchcall,
		ref.index=c(1:5),
		
		# MFA SPECIFIC
#		group=c(5,100),
		type.mfa=rep("s",2),ind.sup=NULL,ncp=5,name.group=NULL,
		num.group.sup=NULL,graph=FALSE,weight.col.mfa=NULL,row.w=NULL,
		axes=c(1,2),tab.comp=NULL,
		
		# PCA SPECIFIC
#		ncp = 5,ind.sup = NULL,row.w = NULL,axes = c(1,2),graph = FALSE, 
		scale.unit = TRUE,  
		quanti.sup = NULL, quali.sup = NULL, 
		col.w = NULL, 
		
		
		# SMFA SPECIFC	
		K=15,para=NULL,type.smfa=c("predictor","Gram"),sparse=c("penalty","varnum"),
		use.corr=FALSE,lambda=1e-6,max.iter=200,trace=FALSE,eps.conv=1e-3,
		sparse.dim=2,
		
		# FABIA SPECIFIC
		p=13,alpha=0.01,cyc=500,spl=0,spz=0.5,non_negative=0,random=1.0,
		center=2,norm=1,scale=0.0,lap=1.0,nL=0,lL=0,bL=0,
		
		
		# COMPONENT SELECTION AND PLOTS
		component.plot=1,
		which=c(1,2,3,4,5),
				
		column.interest=NULL,row.interest=NULL,gene.highlight=NULL,
		
		colour.columns=NULL,labels=TRUE,
		legend.pos="topright",legend.names=NULL,legend.cols=unique(colour.columns),
		thresP.col="blue",thresN.col="red",gene.thresP=NULL,gene.thresN=NULL,
			
		CSrank.refplot=FALSE,profile.type="gene",grouploadings.labels=NULL,grouploadings.cutoff=NULL,
		basefilename="analyseMFA",result.available=NULL,result.available.update=FALSE,plot.type="pdf"
		
		){
			#####################
			### DATA ANALYSIS ###
			#####################
			
			
			### MFA ###
			
			if(modeltype=="CSmfa"){
				if(is.null(result.available)){
					
					group <- c(length(ref.index),(dim(data)[2]-length(ref.index)))
					# Analysis
					result <- MFA(data,	group=group, type=type.mfa, ind.sup=ind.sup,ncp=ncp,name.group=name.group,num.group.sup=num.group.sup,graph=graph,weight.col.mfa=weight.col.mfa,row.w=row.w,axes=axes,tab.comp=tab.comp)
				}
				else{
					result <- result.available@extra$object
				}
				
				loadings <- result$quanti.var$cor	# For compounds	
				scores <- result$ind$coord			# For genes
				
				call.object <- list(match.call=matchcall,analysis.pm=list(ncp=ncp,weight.col.mfa=weight.col.mfa,row.w=row.w))
				
			}
			
						
			### PCA ###
			if(modeltype=="CSpca"){
				# Checking reference index is correct
				if(length(ref.index)>1){stop("There can only be 1 reference compound.",call.=FALSE)}
				# Analysis
				if(is.null(result.available)){
					result <- PCA(X=data,scale.unit=scale.unit,ncp=ncp,ind.sup=ind.sup,quanti.sup=quanti.sup,quali.sup=quali.sup,row.w=row.w,col.w=col.w,graph=graph,axes=axes)
				}
				else{
					result <- result.available@extra$object
				}
				
				loadings <- result$var$cor
				scores <- result$ind$coord
				
				call.object <- list(match.call=matchcall,analysis.pm=list(ncp=ncp,scale.unit=scale.unit,row.w=row.w,col.w=col.w))
				
			}
			
			
			### SMFA ###
			if(modeltype=="CSsmfa"){
				
				## Weighting the data
				total.index <- c(1:dim(data)[2])
				Mat1 <- data[,ref.index,drop=FALSE]
				Mat2 <- data[,total.index[-ref.index]]
				
				if(length(ref.index)>1){
					data.new <- getWeightedDat(Mat1,Mat2,scale.unit=TRUE)[[1]]
				}else{
					data.new <- cbind(Mat1,Mat2)
				}
			
				
				## Doing sMFA Analysis
				if(is.null(result.available)){
					
					if(!(sparse.dim %in% c(1,2))){stop("Please use a correct sparse.dim")}
					
					if(sparse.dim==1){data.spca <- t(data.new)}else{data.spca <- data.new}
					
					if(lambda==Inf){
						if(dim(data.spca)[2]<dim(data.spca)[1]){warning("The to-be-reduced-with-sparsness dimension  is larger than the other dimension. Consider putting lambda at 0.")}
												
						result <- arrayspc(x=data.spca,K=K,para=para,use.corr=use.corr,max.iter=max.iter,trace=trace,eps=eps.conv)
					}
					else{
						if(dim(data.spca)[2]>dim(data.spca)[1]){warning("The to-be-reduced-with-sparsness dimension  is larger than the other dimension. Consider putting lambda at Inf.")}
						
						result <- spca(x=data.spca,K=K,para=para,type=type.smfa,sparse=sparse,use.corr=use.corr,lambda=lambda,max.iter=max.iter,trace=trace,eps.conv=eps.conv)
					}
				}
				else{
					result <- result.available@extra$object
				}
				
				if(sparse.dim==2){
					loadings <- result$loadings	# For compounds	
					scores <- data.new %*% loadings	# For genes
					
				}
				else if(sparse.dim==1){
					scores <- result$loadings	# For genes	
					loadings <- data.spca %*% scores	# For compounds
				}
				
				call.object <- list(match.call=matchcall,analysis.pm=list(K=K,lambda=lambda,sparse=sparse,max.iter=max.iter,eps.conv=eps.conv,para=para,sparse.dim=sparse.dim))
				
			}
			
			
			### FABIA ###
			if(modeltype=="CSfabia"){
				total.index <- c(1:dim(data)[2])
				Mat1 <- data[,ref.index,drop=FALSE]
				Mat2 <- data[,total.index[-ref.index]]
				
				data.new <- getWeightedDat(Mat1,Mat2,scale.unit=TRUE)[[1]]
				
				if(is.null(result.available)){
					result <- fabia(X=t(data.new),p=p,alpha=alpha,cyc=cyc,spl=spl,spz=spz,non_negative=non_negative,random=random,center=center,norm=norm,scale=scale,lap=lap,nL=nL,lL=lL,bL=bL)
				}
				else{
					result <- result.available@extra$object
				}
								
				loadings <- result@L
				scores <- t(result@Z)
			
				call.object <- list(match.call=matchcall,analysis.pm=list(p=p,alpha=alpha,cyc=cyc,spl=spl,spz=spz,non_negative=non_negative,random=random,center=center,norm=norm,scale=scale,lap=lap,nL=nL,lL=lL,bL=bL))
				
			}
			
			######################################
			### APPLY PLOTS AND COMPUTE CSRANK ###
			######################################
			
			
			out <- analyse_FA2(data=data,result=result,loadings=loadings,scores=scores,
					ref.index=ref.index,modeltype=modeltype,
					component.plot=component.plot,
					which=which,labels=labels,
					column.interest=column.interest,row.interest=row.interest,gene.highlight=gene.highlight,
					colour.columns=colour.columns,
					legend.pos=legend.pos,legend.names=legend.names,legend.cols=legend.cols,
					thresP.col=thresP.col,thresN.col=thresN.col,gene.thresP=gene.thresP,gene.thresN=gene.thresN,
					CSrank.refplot=CSrank.refplot,profile.type=profile.type,grouploadings.labels=grouploadings.labels,
					grouploadings.cutoff=grouploadings.cutoff,
					basefilename=basefilename,result.available=result.available,plot.type=plot.type
			)
			
			
			component.select <- out$component.select
			
			loadings <- out$loadings
			scores <- out$scores
			CSRank <- out$CSRank
			sample.factorlabels <- out$factorlabels
			
			CS <- vector("list",length(component.select))
			
			for(i.component in 1:length(component.select)){
				CS[[i.component]] <- list(
						CS.ref=data.frame(CLoadings=loadings[c(1:length(ref.index)),component.select[i.component]],row.names=rownames(loadings)[c(1:length(ref.index))]),
						CS.query=data.frame(
								CLoadings=loadings[-c(1:length(ref.index)),component.select[i.component]],
								CRankScores=CSRank[[i.component]][,1],
								CLrank=as.integer(rank(-loadings[-c(1:length(ref.index)),component.select[i.component]])),
								CLabsrank=as.integer(rank(-abs(loadings[-c(1:length(ref.index)),component.select[i.component]]))),
								CRrank=as.integer(rank(-CSRank[[i.component]][,1])),
								CRabsrank=as.integer(rank(-abs(CSRank[[i.component]][,1]))),
								row.names=rownames(loadings)[-c(1:length(ref.index))]
						)
				)
			}
			names(CS) <- paste0(out$component.name,component.select)
			
			
			GS <- data.frame(scores[,component.select])
			rownames(GS) <- rownames(scores)
			colnames(GS) <- paste0(out$component.name,component.select)
			
			call.object$component.select=component.select
			
			call.object$dimensions <- list(row=dim(data)[1],col=c(ref=length(ref.index),query=(dim(data)[2]-length(ref.index))))
			
			if(!is.null(result.available)){ 
				if(result.available.update){ # Chance to update object
					
					if(!is.null(result.available@permutation.object) ){ # Case for when CSpermute was used
						if(!identical(unname(sort(component.select)),unname(sort(c(result.available@call$component.select,result.available@permutation.object$extra.parameters$component.select))))){
							warning("CS and GS slot will be overwritten due to a different Component choice. \n Any p-values were also deleted together with the permutation.object slot. ")
							return(new("CSresult",type=modeltype,CS=CS,GS=GS,extra=list(CSRank.Full=CSRank,object=result,samplefactorlabels=sample.factorlabels),permutation.object=NULL,call=call.object))
						}else{
							return(new("CSresult",type=modeltype,CS=result.available@CS,GS=result.available@GS,extra=result.available@extra,permutation.object=result.available@permutation.object,call=result.available@call))
						}
					}else{ # Case for when no permutation was used
						if(!identical(unname(sort(component.select)),unname(sort(result.available@call$component.select)))){warning("CS and GS slot will be overwritten due to a different Component choice.")}
						return(new("CSresult",type=modeltype,CS=CS,GS=GS,extra=list(CSRank.Full=CSRank,object=result,samplefactorlabels=sample.factorlabels),permutation.object=NULL,call=call.object))
					}
					
				}else{
					return(new("CSresult",type=modeltype,CS=result.available@CS,GS=result.available@GS,extra=result.available@extra,permutation.object=result.available@permutation.object,call=result.available@call))
				}
				
			}else{
				return(new("CSresult",type=modeltype,CS=CS,GS=GS,extra=list(CSRank.Full=CSRank,object=result,samplefactorlabels=sample.factorlabels),permutation.object=NULL,call=call.object))
			}
				
}





#### GENERAL FUNCTION TO MAKE PLOTS AS WELL AS CSRANKSCORES ####

## which:	-1: Percentage Variance Explained by factors
##			-2: Loadings for reference compounds
##			-3: Loadings for 'component.plot'
##			-4: Genes for 'component.plot'
##			-5: CS Rank Scores for 'factor.plot'
##			-6: Factor 'factor.plot' VS Factor '?' : Loadings & Genes
##			-7: Compound Profiles (Select if necessary)

analyse_FA2 <- function(data,result,loadings,scores,ref.index,modeltype,
		component.plot,which,column.interest,row.interest,gene.highlight,
		colour.columns,legend.pos,legend.names,legend.cols,thresP.col,
		thresN.col,gene.thresP,gene.thresN,CSrank.refplot,labels,
		profile.type,grouploadings.labels,grouploadings.cutoff,basefilename,result.available,plot.type){
	
	
	## Plot-in and -out functions ##
	plot.in <- function(plot.type,name){
		if(plot.type=="pdf"){pdf(name)}
		if(plot.type=="device"){dev.new()}
		if(plot.type=="sweave"){}
	}
	plot.out <- function(plot.type){if(plot.type=="pdf"){dev.off()}}
		
	
	legend.names.loadings <- NULL
	legend.cols.loadings <- NULL
	if(is.null(legend.names) & is.null(legend.cols) &is.null(colour.columns)){
		
		colour.columns <- c(rep("blue",length(ref.index)),rep("black",(dim(data)[2]-length(ref.index))))
		legend.names.loadings <- c("References")
		legend.cols.loadings <- c("blue")
	}
	if(is.null(legend.names)){legend.names <- c()}
	if(is.null(legend.cols)){legend.cols <- c()}
	if(is.null(colour.columns)){colour.columns <- rep("black",dim(data)[2])}
	
	if(is.null(legend.names.loadings)){legend.names.loadings <- legend.names}
	if(is.null(legend.cols.loadings)){legend.cols.loadings <- legend.cols}
	
	
	## IS PVALUE DATA AVAILABLE? ##
	if(!is.null(result.available)){
		if(!is.null(result.available@permutation.object)){
			add.pvalue.color <- TRUE
		}else{
			add.pvalue.color <- FALSE
		}
	}else{
		add.pvalue.color <- FALSE
	}	

	
	## Colors ##
	if(!is.null(colour.columns)){groupCol <- colour.columns} else { groupCol <- "black"}
	
	
	## Specific Naming ##
	if(modeltype=="CSfabia"){
		component.name <- "BC"
		method.name <- "FABIA"
	}
	if(modeltype=="CSpca"){
		component.name <- "PC"
		method.name <- "PCA"	
	}
	if(modeltype=="CSmfa"){
		component.name <- "Factor"
		method.name <- "MFA"	
	}
	if(modeltype=="CSsmfa"){
		component.name <- "Factor"
		method.name <- "sMFA"	
	}
	
	
		
	## PERCENTAGE VARIANCE EXPLAINED / INFORMATION CONTENT ##
	if(1 %in% which){
		if(modeltype=="CSfabia"){
			plot.in(plot.type,paste0(basefilename,"_IC.pdf"))
			showSelected(result,which=c(1))
			plot.out(plot.type)
		}
	  
	#   else{
	# 		if(modeltype=="CSsmfa"){
	# 			perc.var <- result$pev			
	# 		}else{
	# 			perc.var <- result$eig[,2]
	# 		}
	# 		
	# 		plot.in(plot.type,paste0(basefilename,"_percvar.pdf"))
	# 		plot(perc.var,main="Percentage of Variance Explained",xlab="Number of Components",ylab="Perc. Var. Explained")
	# 	}
		# plot.out(plot.type)
	}
	
	## LOADINGS FOR REFERENCE COMPOUNDS ##
	if(2 %in% which){
		plot.in(plot.type,paste0(basefilename,"_RefLoadings.pdf"))
		plot(0,0,type="n",main=paste0(method.name," Loadings for Ref ",paste(ref.index,collapse=",")),xlab=paste0(component.name," Index"),ylab="Loadings",ylim=c(min(loadings[ref.index,]),max(loadings[ref.index,])),xlim=c(1,dim(loadings)[2]))
		for(i.ref in ref.index){
			points(c(1:dim(loadings)[2]),loadings[i.ref,],col=i.ref)
		}
		abline(0,0,lty=3)
		legend(legend.pos,colnames(data)[ref.index],col=ref.index,bty="n",pch=21)
		plot.out(plot.type)
	}
	
	
	# Select Components from 'loadings for reference compounds'	(redraw if necessary) #
	if(is.null(component.plot) 
#			& 
#			((3 %in% which) | (4 %in% which)| (5 %in% which) | (6 %in% which)  | 		
#			(is.null(column.interest)&(7 %in% which)) | ((7 %in% which) & ( (!is.null(gene.thresP)) |(!is.null(gene.thresN)))) )
		){
		
		if(plot.type=="pdf" | !(2 %in% which)){
			dev.new()
			plot(0,0,type="n",main=paste0(method.name," Loadings for Ref ",paste(ref.index,collapse=",")),xlab=paste0(component.name," Index"),ylab="Loadings",ylim=c(min(loadings[ref.index,]),max(loadings[ref.index,])),xlim=c(1,dim(loadings)[2]))
			for(i.ref in ref.index){
				points(c(1:dim(loadings)[2]),loadings[i.ref,],col=i.ref)
			}
			abline(0,0,lty=3)
			legend(legend.pos,colnames(data)[ref.index],col=ref.index,bty="n",pch=21)
		}
		if(modeltype=="CSfabia"){
			out.overlap <- fabia.overlap(result,ref.index)
			cat("Bicluster Suggestions (BC's which overlap with Ref & Query with default thresholds):\n")
			cat("------------------------------------------------------------------------------------\n\n")
			print(out.overlap)
		}

		cat("\n",paste0("Please select with left mousebutton which ",component.name,"s should be investigated.") ,"\n If multiple reference were used, click in the middle of the group of points.\nRight-click to end selection procedure.\n")
		y.mean <- apply(loadings[ref.index,,drop=FALSE],MARGIN=2,FUN=mean)
		component.plot <- identify(x=c(1:dim(loadings)[2]),y=y.mean,plot=TRUE,n=dim(loadings)[2],tolerance=0.5)
	}
	
	if(length(component.plot)==0){stop("Select or input 1 or more components!",call.=FALSE)}
	
	## LOADINGS FOR COMPOUNDS ##
	column.interest.list <- vector("list",max(component.plot))
	
	for(i.component in c(1:length(component.plot))){
		
		if(3 %in% which){
			
			signCol <- "grey"		
			if(add.pvalue.color){ # MFA.FACTOR NEEDS TO BE CHANGED!!!
				if(result.available@permutation.object$extra.parameters$mfa.factor==i.component){
					pvaltype <- ifelse(result.available@permutation.object$extra.parameters$method.adjust=="none","pvalues","pvalues.adjusted")
					signCol <- ifelse(result.available@permutation.object$CS.pval.dataframe[,pvaltype]<=0.05,"purple","grey")
					signCol <- c(rep("grey",length(ref.index)),signCol)
				}
			}
			
			plot.in(plot.type,paste0(basefilename,"_component",component.plot[i.component],"Loadings.pdf"))
			par(mfrow=c(1,1))
			plot(c(1:length(loadings[,component.plot[i.component]])),loadings[,component.plot[i.component]],  type="p",
					xlab="Compound Index", 
					ylab="Compound Loadings",
					pch=21,
					bg=signCol,
					col=groupCol,
					cex=1,main=paste0(method.name," - ",component.name," ",component.plot[i.component]," - Compound Loadings")
			)
			if(labels){
			  text(c(1:length(loadings[,component.plot[i.component]])),loadings[,component.plot[i.component]], colnames(data),	pos=1,	cex=0.5,	col=groupCol)
			 }
			
			legend.names.loadings2 <- legend.names.loadings
			legend.cols.loadings2 <- legend.cols.loadings
			
			legend.bg2 <- "white"
			if(add.pvalue.color){
				if(result.available@permutation.object$extra.parameters$mfa.factor==i.component){
					legend.names.loadings2 <- c(legend.names.loadings2,paste0(result.available@permutation.object$extra.parameters$method.adjust," adj. p-value <= 0.05"))
					legend.bg2 <- c(rep("white",length(legend.names.loadings2)-1),"purple")
					legend.cols.loadings2 <- c(legend.cols.loadings2,"white")	
				}
			}
			
			if(length(legend.names.loadings)>0){legend(legend.pos,legend.names.loadings2,pch=21,col=legend.cols.loadings2,pt.bg=legend.bg2,bty="n")}
			plot.out(plot.type)
		}
		
		# Selecting compounds of interest (replot in device if necessary)
		if(7 %in% which){
			if(is.null(column.interest)){
				if(plot.type=="pdf" | !(3 %in% which)){
					dev.new()
					par(mfrow=c(1,1))
					plot(c(1:length(loadings[,component.plot[i.component]])),loadings[,component.plot[i.component]],  type="p",
							xlab="Compound Index", 
							ylab="Compound Loadings",
							pch=21,
							bg="grey",
							col=groupCol,
							cex=1,main=paste0(method.name," - ",component.name," ",component.plot[i.component]," - Compound Loadings")
					)
					text(c(1:length(loadings[,component.plot[i.component]])),loadings[,component.plot[i.component]], colnames(data),	pos=1,	cex=0.5,	col=groupCol)
					if(length(legend.names.loadings)>0){legend(legend.pos,max(loadings),legend.names.loadings,pch=21,col=legend.cols.loadings,pt.bg="white",bty="n")}
					
				}
				
				cat("Select as many compounds as desired with left mouse button. Right-click to end selection procedure.\n\n")
				id.temp <- identify(c(1:length(loadings[,component.plot[i.component]])),loadings[,component.plot[i.component]],n=999,labels=colnames(data))
				if(!(length(id.temp)==0)){
					column.interest.list[[component.plot[i.component]]] <- id.temp
				}
#				else if(profile.type=="gene" & (7 %in% which)){
#					
#				}
				
				
			}
			else{
				column.interest.list[[component.plot[i.component]]] <- column.interest
				
			}
		}
		
	}
		
	## SCORES FOR GENES ##
	row.interest.list <- vector("list",max(component.plot))
	
	for(i.component in component.plot){
		if(4 %in% which){
			
			bg.temp <-  col.temp <- rep("grey",length(scores[,i.component]))
			if(!is.null(gene.thresP)){
				temp.boolean <- (scores[,i.component]>=gene.thresP)
				bg.temp[temp.boolean] <- thresP.col
				col.temp[temp.boolean] <- thresP.col
			}
			if(!is.null(gene.thresN)){					
				temp.boolean <- (scores[,i.component]<=gene.thresN)
				bg.temp[temp.boolean] <- thresN.col
				col.temp[temp.boolean] <- thresN.col
			}
			
			# highlighting genes
			if(!is.null(gene.highlight)){
				if(class(gene.highlight)=="numeric" | class(gene.highlight)=="integer"){
					col.temp[gene.highlight] <- "green"
				}
				if(class(gene.highlight)=="list"){
					if(length(gene.highlight)>5){stop("Too many different gene.highlight",call.=FALSE)}
					col.pool <- c("green","deeppink","darkorchid4","gold3","tan1")
					for(i.list in 1:length(gene.highlight)){
						col.temp[gene.highlight[[i.list]]] <- col.pool[i.list]
					}
				}
			}
			
			
			## Plot
			plot.in(plot.type,paste0(basefilename,"_component",i.component,"Genescores.pdf"))
			plot(c(1:length(scores[,i.component])),scores[,i.component],  type="p",
					xlab="Gene Index", 
					ylab="Gene Scores",
					pch=21,
					bg=bg.temp,
					col=col.temp,
					cex=1,main=paste0(method.name," - ",component.name," ",i.component," - Gene Factor Scores")
			)
			if(labels){
			  text(c(1:length(scores[,i.component])),scores[,i.component], rownames(data),	pos=1,	cex=0.5,	col=col.temp)
			}
			if(!is.null(gene.thresP)){abline(gene.thresP,0,lty=3)}
			if(!is.null(gene.thresN)){abline(gene.thresN,0,lty=3)}
			plot.out(plot.type)
			
		}
		
		# SELECT GENES FOR PROFILES PLOT
		if(7 %in% which & profile.type=="gene"){
			if(is.null(row.interest)){
				if(plot.type=="pdf" | !(4 %in% which)){
					bg.temp <-  col.temp <- rep("grey",length(scores[,i.component]))
					if(!is.null(gene.thresP)){
						temp.boolean <- (scores[,i.component]>=gene.thresP)
						bg.temp[temp.boolean] <- thresP.col
						col.temp[temp.boolean] <- thresP.col
					}
					if(!is.null(gene.thresN)){					
						temp.boolean <- (scores[,i.component]<=gene.thresN)
						bg.temp[temp.boolean] <- thresN.col
						col.temp[temp.boolean] <- thresN.col
					}
					
					# highlighting genes
					if(!is.null(gene.highlight)){
						if(class(gene.highlight)=="numeric" | class(gene.highlight)=="integer"){
							col.temp[gene.highlight] <- "green"
						}
						if(class(gene.highlight)=="list"){
							if(length(gene.highlight)>5){stop("Too many different gene.highlight",call.=FALSE)}
							col.pool <- c("green","deeppink","darkorchid4","gold3","tan1")
							for(i.list in 1:length(gene.highlight)){
								col.temp[gene.highlight[[i.list]]] <- col.pool[i.list]
							}
						}
					}
					dev.new()
					plot(c(1:length(scores[,i.component])),scores[,i.component],  type="p",
							xlab="Gene Index", 
							ylab="Gene Scores",
							pch=21,
							bg=bg.temp,
							col=col.temp,
							cex=1,main=paste0(method.name," - ",component.name," ",i.component," - Gene Factor Scores")
					)
					text(c(1:length(scores[,i.component])),scores[,i.component], rownames(data),	pos=1,	cex=0.5,	col=col.temp)
					if(!is.null(gene.thresP)){abline(gene.thresP,0,lty=3)}
					if(!is.null(gene.thresN)){abline(gene.thresN,0,lty=3)}
				}
				cat("Select as many genes as desired with left mouse button. Right-click to end selection procedure.\n\n")
				id.temp <- identify(c(1:length(scores[,i.component])),scores[,i.component],n=999,labels=rownames(data))
				if(!(length(id.temp)==0)){
					row.interest.list[[i.component]] <- id.temp
				}
				
			}
			else{
				row.interest.list[[i.component]] <- row.interest
			}
		}
		
	}
	
	## CSRANK SCORES ##
	
	if(5%in%which){
		out_CS_rank <- vector("list",length(component.plot))
		
		legend.bg <- "white"
		legend.names.csrank <- legend.names
		legend.cols.csrank <- legend.cols
		
		signCol <- "grey"		
		if(add.pvalue.color){ # MFA.FACTOR NEEDS TO BE CHANGED!!!
			if(result.available@permutation.object$extra.parameters$mfa.factor==i.component){
				pvaltype <- ifelse(result.available@permutation.object$extra.parameters$method.adjust=="none","pvalues","pvalues.adjusted")
				signCol <- ifelse(result.available@permutation.object$CSRank.pval.dataframe[,pvaltype]<=0.05,"purple","grey")
			
				legend.names.csrank <- c(legend.names.csrank,paste0(result.available@permutation.object$extra.parameters$method.adjust," adj. p-value <= 0.05"))
				legend.bg <- c(rep("white",length(legend.names.csrank)-1),"purple")
				legend.cols.csrank <- c(legend.cols.csrank,"white")
			}
		}
		
		
		for(i.component in c(1:length(component.plot))){
			out_CS_rank[[i.component]] <- CSrank2(loadings,ref.index,color.columns=groupCol,signCol=signCol,ref.plot=CSrank.refplot,legend.pos=legend.pos,loadings_names=colnames(data),component.plot=component.plot[i.component],type.component=component.name,plot=TRUE,plot.type=plot.type,basefilename=basefilename,legend.bg=legend.bg,legend.names=legend.names.csrank,legend.cols=legend.cols.csrank,labels=labels)
		}
		names(out_CS_rank) <- paste0(component.name,component.plot)
	}
	else{
		out_CS_rank <- replicate(length(component.plot),list)
		
		for(i.component in c(1:length(component.plot))){
			out_CS_rank[[i.component]] <- CSrank2(loadings,ref.index,color.columns=groupCol,ref.plot=CSrank.refplot,loadings_names=colnames(data),component.plot=component.plot[i.component],type.component=component.name,plot=FALSE)
		}
		names(out_CS_rank) <- paste0(component.name,component.plot)
	}
	
	
	## COMPONENT VS COMPONENT (LOADINGS & SCORES) ##

	
	if(6 %in% which){
		
		if(length(component.plot)>1){
			component1.plot <- component.plot[1]
			component2.plot <- component.plot[2]
		}
		else{
			component1.plot <- component.plot
			component2.plot <- ifelse(component1.plot==1,2,1)
		}
		
		##  Loadings
		plot.in(plot.type,paste0(basefilename,"_",component.name,component1.plot,component.name,component2.plot,"_loadings.pdf"))
		par(mfrow=c(1,1))
		plot(loadings[,component1.plot],loadings[,component2.plot],  type="p",xlim=c(min(loadings[,component1.plot]),max(loadings[,component1.plot])),ylim=c(min(loadings[,component2.plot]),max(loadings[,component2.plot])),
				xlab=paste0(component.name," ",component1.plot," - Compound Loadings"), 
				ylab=paste0(component.name," ",component2.plot," - Compound Loadings"),
				pch=21,
				bg="grey",
				col=groupCol,
				cex=1,main=paste0(method.name," - ",component.name," ",component1.plot," vs ",component2.plot)
		)
		if(labels){
		  text(loadings[,component1.plot],loadings[,component2.plot], colnames(data),	pos=1,	cex=0.5,	col=groupCol)
		}
		if(length(legend.names)>0){legend(legend.pos,legend.names,pch=21,col=legend.cols,pt.bg="white",bty="n")}
		plot.out(plot.type)
		
		
		##  Factors
		plot.in(plot.type,paste0(basefilename,"_",component.name,component1.plot,component.name,component2.plot,"_factors.pdf"))
		plot(scores[,component1.plot],scores[,component2.plot],  type="p",
				xlab=paste0(component.name," ",component1.plot,"- Gene Factor Scores"), 
				ylab=paste0(component.name," ",component2.plot," - Gene Factor Scores"),
				pch=21,
				bg="grey",
				col="grey",
				cex=1,main=paste0(method.name," - ",component.name," ",component1.plot," vs ",component2.plot)
		)
		if(labels){
		  text(scores[,component1.plot],scores[,component2.plot], rownames(data),	pos=1,	cex=0.5,	col="grey")
		}
		plot.out(plot.type)
	}
	
	## PROFILE PLOTS ##
	if( (7 %in% which) & (length(column.interest.list)>0) ){# Excluding the case in which you don't select anything + no pm
		
		if(!(profile.type=="gene" & !length(row.interest.list)>0)){ # case of gene profiles but no selection not allowed
			
			for(i.component in component.plot){
				cmpds_interest <- column.interest.list[[i.component]]
				
				if(length(row.interest.list)==0){
					genes_interest <- vector("list",length(cmpds_interest))
				}
				else{
					genes_interest <- row.interest.list[[i.component]]
					
				}
				
				if(!is.null(cmpds_interest)){
					if(!(profile.type=="gene" & is.null(genes_interest))){
												
						base.name <- paste0(basefilename,"_",component.name,i.component)
						main.base <- paste0(method.name," ",component.name,i.component)
						
						CSprofiles(data=data,ref_index=ref.index,gene.select=genes_interest,cmpd.select=cmpds_interest,profile.type=profile.type,cmpd.loadings=loadings,gene.scores=scores,component.plot=i.component,gene.thresP=gene.thresP,gene.thresN=gene.thresN,basefilename=base.name,plot.type=plot.type,thresP.col=thresP.col,thresN.col=thresN.col,main.base=main.base)
						
					}
				}
			}
		}
	}
	
	## TREND PLOTS ##
#	if((8 %in% which & ....)){
#		
#	}

	## GROUPED LOADINGS PLOTS ##
	if(8 %in% which){		
		sample.factorlabels <- CSgrouploadings(loadings=loadings,grouploadings.labels=grouploadings.labels,grouploadings.cutoff=grouploadings.cutoff,ref.index=ref.index,method.name=method.name,component.name=component.name,basefilename=basefilename,plot.type=plot.type,plot=TRUE)
	}else{
		sample.factorlabels <- CSgrouploadings(loadings=loadings,grouploadings.labels=grouploadings.labels,grouploadings.cutoff=grouploadings.cutoff,ref.index=ref.index,method.name=method.name,component.name=component.name,basefilename=basefilename,plot.type=plot.type,plot=FALSE)
	}
	
	## RETURN OBJECT ##
	out <- list(component.select=component.plot,result=result,loadings=loadings,scores=scores,CSRank=out_CS_rank,component.name=component.name,factorlabels=sample.factorlabels)	
	return(out)
	
	


}



#' @title Internal Function: Zhang and Gant Analysis
#' @description \strong{Internal function not meant for end-users.} Internal function involved in the computation of the Zhang and Gant score when using multiple cores with \code{parallel}.
#' @export
#' @keywords internal
#' @param dataref NA
#' @param dataquery NA
#' @param nref NA
#' @param nquery NA
#' @param ord.query NA
#' @param permute NA
#' @param B NA
#' @param ntop.pvalues NA
#' @param ntop.scores NA
#' @param basefilename NA
#' @param colour.query NA
#' @param legend.names NA
#' @param legend.cols NA
#' @param legend.pos NA
#' @param result.available NA
#' @param plot.type NA
#' @param print.top NA
#' @param which NA
#' @return NA
analyse_zhang <- function(dataref,dataquery,nref=NULL,nquery=NULL,ord.query=TRUE,permute=FALSE,B=100000,ntop.pvalues=20,ntop.scores=20,
		basefilename="analyseZhang",
		#column.interest=NULL,
		colour.query=NULL,legend.names=NULL,legend.cols=unique(colour.query),legend.pos="topright",
		result.available=NULL,plot.type="pdf",print.top=TRUE,labels=TRUE,
		which=c(1)){
	
	
	## Plot-in and -out functions			
	plot.in <- function(plot.type,name){
		if(plot.type=="pdf"){pdf(name)}
		if(plot.type=="device"){dev.new()}
		if(plot.type=="sweave"){}
	}
	plot.out <- function(plot.type){if(plot.type=="pdf"){dev.off()}}
	
	
	## Doing Zhang Analysis
	if(is.null(result.available)){
		zhang_result <- zhangscore(dataref=dataref,dataquery=dataquery,nref=nref,nquery=nquery,ord.query=ord.query,permute=permute,B=B,ntop.pvalues=ntop.pvalues,ntop.scores=ntop.scores,pvalue.method="BH")
		add.pvalue.color <- FALSE
		
	}
	else{
		zhang_result <- result.available@extra$object
		
		if(!is.null(result.available@permutation.object)){
			add.pvalue.color <- TRUE	
		}else{
			add.pvalue.color <- FALSE
		}
		
	}
	
	## Printing Top Connections Scores + Pvalues
#	if(print.top){print(zhang_result$Top)}
#	if(permute==TRUE){print(zhang_result$Toppvalues)}
	
	##
	if(!is.null(colour.query)){groupCol <- colour.query} else { groupCol <- "black"}
	if(add.pvalue.color){
		pvaltype <- ifelse(result.available@permutation.object$extra.parameters$method.adjust=="none","pvalues","pvalues.adjusted")
		signCol <- ifelse(result.available@permutation.object$CS.pval.dataframe[,pvaltype]<=0.05,"purple","grey")
	}
	else{
		signCol <- "grey"
	}
	
	
	ncol_q <- ncol(dataquery)
	##
	
	if(1 %in% which){
		## PLOT: Zhang Score Plot
		plot.in(plot.type,paste0(basefilename,"_zhangscore.pdf"))
		par(mfrow=c(1,1))
		plot(zhang_result$All[,1],xlim=c(1,ncol_q),ylim=c(-1,1),pch=21,bg=signCol,col=groupCol,main=paste0("Zhang Score"),xlab="Compound Index",ylab="Connection Score")
		if(labels){
		  text(c(1:length(zhang_result$All[,1])),zhang_result$All[,1],rownames(zhang_result$All),	pos=1,	cex=0.5,	col=groupCol)
		}
		abline(0,0,lty=3)
		legend.bg <- c()
		if(add.pvalue.color){
			legend.names <- c(legend.names,paste0(result.available@permutation.object$extra.parameters$method.adjust," adj. p-value <= 0.05"))
			legend.bg <- c(rep("white",length(legend.names)-1),"purple")
			legend.cols <- c(legend.cols,"white")
		}
		else{
			legend.bg <- "white"
		}
#		if(is.null(legend.x)){legend.x <- (ncol_q-0.45*ncol_q)}
#		if(is.null(legend.y)){legend.y <- 1}
		if(length(legend.names)>0){legend(legend.pos,legend.names,pch=21,col=legend.cols,pt.bg=legend.bg,bty="n")}
		plot.out(plot.type)
	}	
	
	return(zhang_result)
}

CSrank2 <- function(loadings,ref.index,color.columns=NULL,ref.plot=FALSE,loadings_names=NULL,component.plot,type.component="Factor",plot=TRUE,plot.type="pdf",basefilename="base",signCol="grey",legend.pos="topright",legend.bg="grey",legend.names="",legend.cols="black",labels=TRUE){
	
	term1.w <- 0.5
	term2.w <- 0.5
	
	## Plot-in and -out functions
	plot.in <- function(plot.type,name){
		if(plot.type=="pdf"){pdf(name)}
		if(plot.type=="device"){dev.new()}
		if(plot.type=="sweave"){}
	}
	plot.out <- function(plot.type){if(plot.type=="pdf"){dev.off()}}
	
	names <- loadings_names
	loadings <- loadings[,component.plot]
	
	RefCS.matrix <- matrix(0,nrow=length(loadings[-ref.index]),ncol=length(ref.index))
	colnames(RefCS.matrix) <- rep("",length(ref.index))
	
	query.loadings <- loadings[-ref.index]
	
	
	# |q|/|R| for each ref  
	for(i.ref in 1:length(ref.index)){
		ref.loading <- loadings[ref.index[i.ref]]	
		
		
		temp <- sapply(query.loadings,FUN=function(x){
					
					if(abs(ref.loading)>=abs(x)){
						return(abs(x)/abs(ref.loading))
					}
					else{
						return(abs(ref.loading)/abs(x))
					}
				})
		
		RefCS.matrix[,i.ref] <- temp
		colnames(RefCS.matrix)[i.ref] <- paste0(names(loadings)[ref.index[i.ref]]," (Ref ",i.ref,")")
	}
	
	weights <- abs(loadings[ref.index])
	term1 <- apply(RefCS.matrix,MARGIN=1,FUN=weighted.mean,w=weights)
	
	
	# Extra term for each q
	query.mean <- mean(abs(query.loadings))
	term2.temp <- abs(query.loadings)-query.mean
	
	term2 <- term2.temp/max(abs(term2.temp))
	
	# Combine both terms with (weighted) mean  +  change to 0 if result is negative
	combinetest1 <- sapply(1:length(term1),FUN=function(x){weighted.mean(c(term1[x],term2[x]),c(term1.w,term2.w))})
	combinetest2 <- sapply(combinetest1,FUN=function(x){if(x<0){return(0)}else{return(x)}})
	
	# Add Sign to final Score
	sign.temp <- sapply(query.loadings,FUN=function(x){
				ifelse(sign(mean(loadings[ref.index]))==sign(x),1,-1)
			})
	
	# Final
	combinetest3 <- combinetest2*sign.temp
	
	out <- as.data.frame(cbind(combinetest3,term1,RefCS.matrix,term2,combinetest1,combinetest2))
	colnames(out) <- c("CSRankScores","Term1.Weighted",paste0("Term1.Ref",c(1:length(ref.index))),"Term2","Combine1.wmean","Combine2.zero")
	rownames(out) <- names[-ref.index]
	
	
	## PLOTTING
	
	if(plot){
		
		if(ref.plot){
			
			nref <- dim(RefCS.matrix)[2]
			nplots <- dim(out)[2]-1
			
			col.plots <- ifelse(nplots<4,nplots,4)
			row.plots <- ifelse(nplots%%4==0,nplots%/%4,nplots%/%4 +1)
			par(mfrow=c(row.plots,col.plots))
			
			plot.in(plot.type,paste0(basefilename,"_CSRankExtra.pdf"))
			for(i.plot in 2:dim(out)[2]){
				plot(out[,i.plot],xlab="Query Compound Index",ylab="Score",col=color.columns[-ref.index],bg="grey",pch=21,main=paste0(type.component," ",component.plot," - ",colnames(out)[i.plot]))
				text(out[,i.plot],names[-ref.index],col=color.columns[-ref.index],pos=2)
			}
			plot.out(plot.type)
			
		}
		
		par(mfrow=c(1,1))
		
		## Plot Final
		main.temp <- ifelse(ref.plot,"CS Rank Score (Final)","CS Rank Score")
		plot.in(plot.type,paste0(basefilename,"_CSRank.pdf"))
		plot(combinetest3,col=color.columns[-ref.index],xlab="Query Compound Index",ylab="CS Rankscore",bg=signCol,pch=21,main=paste0(type.component," ",component.plot," - ",main.temp))
		
		if(label){
		  text(combinetest3,labels=names[-ref.index],col=color.columns[-ref.index],pos=2)
		}
		if(length(legend.names)>0){legend(legend.pos,legend.names,pch=21,col=legend.cols,pt.bg=legend.bg,bty="n")}
		
		plot.out(plot.type)
	}
	
	
	
	return(out)
	
}