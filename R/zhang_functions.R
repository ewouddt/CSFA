# Project: Connectivity
# 
# Author: lucp8394
###############################################################################

#load("data/data1ref_ref.RData")
#load("data/data1ref_query.RData")
#load("data/data6ref_ref.RData")
#load("data/data6ref_query.RData")

## ZHANG SCORE 

## Data format: genes in rows, compounds in columns (this is the data on which the genes are ranked)
## 			-> Column and row names will be used!!

# dataref: 1 or more compounds as a reference 
# dataquery: matrix containing 1 or more compounds (score will be computed for each compound wrt the dataref

# nref: Number of top genes to pick from the reference profiles (same for 1 or all of set)
# nquery: Number of top genes to pick from query compounds. This should be a vector with a threshold for each compound.
## IF THE NUMBER IS NOT FILLED IN, ALL GENES WILL BE USED

#ord.query: if set to FALSE, it is assumed the queries are not ordered. "ranks" will be -1 or 1


#zhangscore.support_OLD <- function(dataref,dataquery,nref=NULL,nquery=NULL,ord.query=TRUE,ntop.scores=10){
#	
#	if(is.null(rownames(dataref))){stop("Provide rownames for dataref",call.=FALSE)}
#	if(is.null(rownames(dataquery))){stop("Provide rownames for dataquery",call.=FALSE)}
#	if(ord.query==FALSE & !is.null(nquery)){stop("Cannot take top of non-ordered query",call.=FALSE)}
#
#	# Make ranks for dataref
#	refSign <- sign(dataref)
#	dataref_Rank <- apply(abs(dataref), 2, rank) * refSign
#	
#	
#	# Make ranks for query
#	querySign <- sign(dataquery)
#	if(ord.query==TRUE){dataquery_Rank <- apply(abs(dataquery),2,rank)*querySign} else {dataquery_Rank <- querySign}
#	
#	if(is.null(nref)){nref <- nrow(dataref)}
#	if(is.null(nquery)){nquery <- nrow(dataquery)} 
#	if(sum(nref<nquery)>0){stop("Number of chosen genes larger in query than in reference. ",call.=FALSE)}
#	
#	scores_query <- NULL
#	
#	# Go into a for-loop for each column in the query
#	for(i.query in 1:ncol(dataquery)){
#				
#		# Go into for-loop for when more than 1 compound is in the reference set. (Mean will be taken in the end)
#		refset_scores <- NULL
#		for(i.ref in 1:ncol(dataref)){
#					
#			# Sort and take threshold
#			sort.ref <- dataref_Rank[order(abs(dataref_Rank[,i.ref]),decreasing=TRUE),i.ref]
#			sort.query <- dataquery_Rank[order(abs(dataquery_Rank[,i.query]),decreasing=TRUE),i.query]
#			
#			thres.ref <- sort.ref[1:nref]
#			thres.query <- sort.query[1:nquery]
#			
#			# Multiply + Sum common genes
#			
#			if(length(intersect(names(thres.ref),names(thres.query)))>0){
#				if(ord.query==TRUE){
#					maxTheoreticalScore <- sum(abs(thres.ref)[1:length(thres.query)]*abs(thres.query))
#				}else{
#					maxTheoreticalScore <- sum(abs(thres.ref)[1:length(thres.query)])
#				}
#				
#				connscore <- sum(thres.query * thres.ref[names(thres.query)],na.rm=TRUE)/maxTheoreticalScore
#				refset_scores[i.ref] <- connscore 
#			}
#			else{refset_scores[i.ref] <- NA} # If Top N genes of ref do not share genes with top m genes of query
#		}
#		scores_query[i.query] <- mean(refset_scores,na.rm=TRUE)
#		
#	}
#	
#	if(!is.null(colnames(dataquery))){names(scores_query)=colnames(dataquery)}
#	output = data.frame(ZhangScore=scores_query)
#	
#	if(ncol(dataquery)<ntop.scores){ntop <- ncol(dataquery)} else { ntop <- ntop.scores}
#	scores_query_sort <- sort(scores_query,decreasing=TRUE)
#	toppos <- scores_query_sort[1:ntop]
#	topneg <- scores_query_sort[length(scores_query_sort):(length(scores_query_sort)-(ntop-1))]
#	
#	if(!is.null(colnames(dataquery))){
#	topoutput <- data.frame(posname=names(toppos),posscore=toppos,negname=names(topneg),negscore=topneg)
#	}
#	else{
#		topoutput <- data.frame(posscore=toppos,negscore=topneg)
#		
#	}
#	
#	rownames(topoutput) <- NULL
#	
#	listout <- list(Top=topoutput,All=output)
#	return(listout)
#	
#}

zhangscore.support <- function(dataref,dataquery,nref=NULL,nquery=NULL,ord.query=TRUE,ntop.scores=10){
	
	if(is.null(rownames(dataref))){stop("Provide rownames for dataref",call.=FALSE)}
	if(is.null(rownames(dataquery))){stop("Provide rownames for dataquery",call.=FALSE)}
	if(ord.query==FALSE & !is.null(nquery)){stop("Cannot take top of non-ordered query",call.=FALSE)}
	
			
	if(is.null(nref)){nref <- nrow(dataref)}
	if(is.null(nquery)){nquery <- nrow(dataquery)} 
	if(sum(nref<nquery)>0){stop("Number of chosen genes larger in query than in reference. ",call.=FALSE)}
	
	scores_query <- NULL
	
	# Go into a for-loop for each column in the query
	for(i.query in 1:ncol(dataquery)){
		
		# Go into for-loop for when more than 1 compound is in the reference set. (Mean will be taken in the end)
		refset_scores <- NULL
		for(i.ref in 1:ncol(dataref)){
			
			# Sort and take threshold
			sort.ref <- dataref[order(abs(dataref[,i.ref]),decreasing=TRUE),i.ref]
			sort.query <- dataquery[order(abs(dataquery[,i.query]),decreasing=TRUE),i.query]
			
			thres.ref <- sort.ref[1:nref]
			thres.query <- sort.query[1:nquery]
			
			# Make ranks for dataref
			refSign <- sign(thres.ref)
			thres.ref.rank <- rank(abs(thres.ref)) * refSign
			
			# Make ranks for query
			querySign <- sign(thres.query)
			if(ord.query==TRUE){thres.query.rank <- rank(abs(thres.query))*querySign} else {thres.query.rank <- querySign}
			
			
			# Multiply + Sum common genes
			
			if(length(intersect(names(thres.ref.rank),names(thres.query.rank)))>0){
				if(ord.query==TRUE){
							maxTheoreticalScore <- sum(abs(thres.ref.rank)[1:length(thres.query.rank)]*abs(thres.query.rank))
						}else{
							maxTheoreticalScore <- sum(abs(thres.ref.rank)[1:length(thres.query.rank)])
						}
						
				connscore <- sum(thres.query.rank * thres.ref.rank[names(thres.query.rank)],na.rm=TRUE)/maxTheoreticalScore
				refset_scores[i.ref] <- connscore 
			}
			else{refset_scores[i.ref] <- NA} # If Top N genes of ref do not share genes with top m genes of query
		}
		scores_query[i.query] <- mean(refset_scores,na.rm=TRUE)
		
	}
	
	if(!is.null(colnames(dataquery))){names(scores_query)=colnames(dataquery)}
	output = data.frame(ZhangScore=scores_query)
	
	if(ncol(dataquery)<ntop.scores){ntop <- ncol(dataquery)} else { ntop <- ntop.scores}
	scores_query_sort <- sort(scores_query,decreasing=TRUE)
	toppos <- scores_query_sort[1:ntop]
	topneg <- scores_query_sort[length(scores_query_sort):(length(scores_query_sort)-(ntop-1))]
	
	if(!is.null(colnames(dataquery))){
		topoutput <- data.frame(posname=names(toppos),posscore=toppos,negname=names(topneg),negscore=topneg)
	}
	else{
		topoutput <- data.frame(posscore=toppos,negscore=topneg)
		
	}
	
	rownames(topoutput) <- NULL
	
	listout <- list(Top=topoutput,All=output)
	return(listout)
	
}









zhangscore.support2 <- function(dataref,dataquery,nref=NULL,nquery=NULL,ord.query=TRUE,ntop.scores=10){
  
  if(is.null(rownames(dataref))){stop("Provide rownames for dataref",call.=FALSE)}
  if(is.null(rownames(dataquery))){stop("Provide rownames for dataquery",call.=FALSE)}
  if(ord.query==FALSE & !is.null(nquery)){stop("Cannot take top of non-ordered query",call.=FALSE)}
  
  
  if(is.null(nref)){nref <- nrow(dataref)}
  if(is.null(nquery)){nquery <- nrow(dataquery)} 
  if(sum(nref<nquery)>0){stop("Number of chosen genes larger in query than in reference. ",call.=FALSE)}
  
  scores_query <- NULL
  
  # Go into a for-loop for each column in the query
  for(i.ref in 1:ncol(dataref)){
    
    # Go into for-loop for when more than 1 compound is in the reference set. (Mean will be taken in the end)
    querset_scores <- NULL
    for(i.query in 1:ncol(dataquery)){
      
      # Sort and take threshold
      sort.ref <- dataref[order(abs(dataref[,i.ref]),decreasing=TRUE),i.ref]
      sort.query <- dataquery[order(abs(dataquery[,i.query]),decreasing=TRUE),i.query]
      
      thres.ref <- sort.ref[1:nref]
      thres.query <- sort.query[1:nquery]
      
      # Make ranks for dataref
      refSign <- sign(thres.ref)
      thres.ref.rank <- rank(abs(thres.ref)) * refSign
      
      # Make ranks for query
      querySign <- sign(thres.query)
      if(ord.query==TRUE){thres.query.rank <- rank(abs(thres.query))*querySign} else {thres.query.rank <- querySign}
      
      
      # Multiply + Sum common genes
      
      if(length(intersect(names(thres.ref.rank),names(thres.query.rank)))>0){
        if(ord.query==TRUE){
          maxTheoreticalScore <- sum(abs(thres.ref.rank)[1:length(thres.query.rank)]*abs(thres.query.rank))
        }else{
          maxTheoreticalScore <- sum(abs(thres.ref.rank)[1:length(thres.query.rank)])
        }
        
        connscore <- sum(thres.query.rank * thres.ref.rank[names(thres.query.rank)],na.rm=TRUE)/maxTheoreticalScore
        querset_scores[i.query] <- connscore 
      }
      else{querset_scores[i.query] <- NA} # If Top N genes of ref do not share genes with top m genes of query
    }
    scores_query[i.ref] <- mean(querset_scores,na.rm=TRUE)
    
  }
  
  if(!is.null(colnames(dataquery))){names(scores_query)=colnames(dataquery)}
  output = data.frame(ZhangScore=scores_query)
  
  if(ncol(dataquery)<ntop.scores){ntop <- ncol(dataquery)} else { ntop <- ntop.scores}
  scores_query_sort <- sort(scores_query,decreasing=TRUE)
  toppos <- scores_query_sort[1:ntop]
  topneg <- scores_query_sort[length(scores_query_sort):(length(scores_query_sort)-(ntop-1))]
  
  if(!is.null(colnames(dataquery))){
    topoutput <- data.frame(posname=names(toppos),posscore=toppos,negname=names(topneg),negscore=topneg)
  }
  else{
    topoutput <- data.frame(posscore=toppos,negscore=topneg)
    
  }
  
  rownames(topoutput) <- NULL
  
  listout <- list(Top=topoutput,All=output)
  return(listout)
  
}













# This function takes a vector with connection scores (function provided by Ziv from... ? )
threshSamp <- function(connScore, thresholdSD)
{
	sdev=1
	f=NULL
	sdev.all=NULL
	sdev.length=NULL
	thresh = median(connScore) + thresholdSD*(sd(connScore))
	while( (median(connScore) + sdev*(sd(connScore))) <= 1)
	{
		z = (median(connScore) + sdev*(sd(connScore)))
		f = c(f, z)
		sdev.length = c(sdev.length, length(which(connScore > z)))
		sdev.all = c(sdev.all, sdev)
		sdev=sdev+1
	}
	finalTable = cbind(sdev.all, f, sdev.length)
	colnames(finalTable)= c("S.D.", "Connection Score", "# Samples")
	finalTable = finalTable[which(finalTable[,3] > 0),]
	if(class(finalTable) == "numeric")
	{
		finalTable = t(as.matrix(finalTable))
	}
	return(finalTable)
}

### SOME TESTING
#zhangscore.support(data1ref_ref,data1ref_query,nref=100,nquery=50)
#zhangscore.support(data6ref_ref,data6ref_query)
#
#test1 <- as.matrix(data1ref_query[,1])
#rownames(test1) <- rownames(data1ref_query)
#zhangscore.support(data1ref_ref,test1) 
#
#
#
#score = zhangscore.support(data1ref_ref,data1ref_query,ntop.scores=15)
#threshSamp(score$All[,1],1)
###

## Full function which can also compute the p-values
zhangscore <- function(dataref,dataquery,nref=NULL,nquery=NULL,ord.query=TRUE,permute=FALSE,B=100000,ntop.pvalues=10,ntop.scores=10,pvalue.method="BH"){
	
	if(is.null(rownames(dataref))){stop("Provide rownames for dataref",call.=FALSE)}
	if(is.null(rownames(dataquery))){stop("Provide rownames for dataquery",call.=FALSE)}
	if(ord.query==FALSE & !is.null(nquery)){stop("Cannot take top of non-ordered query",call.=FALSE)}
	
	
	
	if(permute==FALSE){
		return(zhangscore.support(dataref,dataquery,nref,nquery,ord.query,ntop.scores))
	}
	else{
		
		observed.scores.full <- zhangscore.support(dataref,dataquery,nref,nquery,ord.query,ntop.scores)
		observed.scores <- observed.scores.full$All
		
		if(is.null(nref)){nref <- nrow(dataref)}
		if(is.null(nquery)){nquery <- nrow(dataquery)} 
		if(sum(nref<nquery)>0){stop("Number of chosen genes larger in query than in reference. ",call.=FALSE)}
		
		# Pool from which will be sampled (if more than 1 ref, random choice from which sample will be taken)
		
		sample.pool.names <- matrix(0,nrow=nref,ncol=ncol(dataref))
		
		for(i.pool in 1:ncol(dataref)){
			temp.sort <- sort(abs(dataref[,i.pool]),decreasing=TRUE)
			
			sample.pool.names[,i.pool] <- names(temp.sort)[1:nref]
		}
		
		boot.values <- NULL
		for(i.B in 1:B){
			which.pool <- sample(1:ncol(dataref),1)
			index.bootsample <- sample(1:nref,nquery)
	
		
			bootsample <- as.matrix((nquery:1)*sample(c(-1,1),nquery,replace=T))
			rownames(bootsample) <- sample.pool.names[index.bootsample,i.pool]
			
			
			out.boot.temp <- zhangscore.support(dataref,bootsample,nref,nquery,ord.query)$All[,1]
			boot.values <- c(boot.values,out.boot.temp)
			
		}
		
		connpval <- function(score,values){
			B <- length(boot.values)
			return(sum(abs(score)<=abs(values))/B)
		}
		all.p.values <- sapply(observed.scores[,1],FUN=connpval,values=boot.values)
		
		if(pvalue.method %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr")){
			all.p.values.adjusted <- p.adjust(all.p.values,method=pvalue.method)
			p.adj <- TRUE
		}
		else{p.adj <- FALSE}
		
		# All output
		if(p.adj){
			output <- cbind(observed.scores,all.p.values,all.p.values.adjusted)
		}
		else{
			output <- cbind(observed.scores,all.p.values)
		}
		
		# Top p values
		if(p.adj){
			Toppvalues <- output[order(output[,3])[1:ntop.pvalues],]
		}
		else{
			Toppvalues <- output[order(output[,2])[1:ntop.pvalues],]
		}
		
		oldtop <- observed.scores.full$Top
		# Adding the pvalues to the top connectivity ones
		posname <- as.character(oldtop[,1])
		negname <- as.character(oldtop[,3])
		
		index.pos <- sapply(posname,FUN=function(x,output){which(x==rownames(output))},output=output)
		index.neg <- sapply(negname,FUN=function(x,output){which(x==rownames(output))},output=output)
		p.pos <- output[index.pos,2]
		p.neg <- output[index.neg,2]
		if(p.adj){
			p.pos.adj <- output[index.pos,3]
			p.neg.adj <- output[index.neg,3]
		}
		
		if(p.adj){
			newtop <- data.frame(posname=oldtop[,1],posscore=oldtop[,2],pospvalue=p.pos,pospvalue.adj=p.pos.adj,negname=oldtop[,3],negscore=oldtop[,4],negpvalue=p.neg,negvalue.adj=p.neg.adj)
			
		}
		else{
			newtop <- data.frame(posname=oldtop[,1],posscore=oldtop[,2],pospvalue=p.pos,negname=oldtop[,3],negscore=oldtop[,4],negpvalue=p.neg)
			
		}
		
		
		listout <- list(Toppvalues=Toppvalues ,Top=newtop ,All=output)
		return(listout)
	}
}

