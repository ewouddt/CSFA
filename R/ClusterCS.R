
## EXAMPLE
# load("C:/Users/lucp8394/Documents/StatFiles/CenStat/Connectivity Score/Research/ClusterCS/MCF7Data.RData")
# load("C:/Users/lucp8394/Documents/StatFiles/CenStat/Connectivity Score/Research/ClusterCS/cluster_target.RData")
# data <- geneMat
# clusterlabels <- cutree(cluster_target$Clust,k=7)
# out <- CScluster(data,clusterlabels,WithinSave=FALSE)
# out <- CScluster(data,clusterlabels,WithinSave=TRUE)
# out <- CScluster(data,clusterlabels,WithinSave=TRUE,Within=c(1,2,5))
# out <- CScluster(data,clusterlabels,WithinSave=TRUE,Between=c(1,2,5))
# out <- CScluster(data,clusterlabels,WithinSave=TRUE,Within=c(1,2,5),Between=c(1,2,5))


# out <- CScluster(data,clusterlabels,WithinSave=TRUE,type="CSfabia",cyc=10)
# out <- CScluster(data,clusterlabels,WithinSave=TRUE,type="CSzhang")

# To do: if time left: allow for parallelisation (over the i loop, although it heavily depends on cluster size... so might not be efficient)



#' @importFrom utils capture.output
#' @title CScluster
#' @description Apply the Connectivity Scores to a \emph{K} clustering result. More information can be found in the Details section below.
#' @details After applying cluster analysis on the additional data matrix, \emph{K} clusters are obtained.
#' Each cluster will be seen as a potential query set (for \code{\link{CSanalysis}}) for which 2 connectivity score metrics can be computed, the \emph{Within-Cluster CS} and the \emph{Between-Cluster CS}.
#' 
#' \emph{Within-Cluster CS}\cr
#' This metric will answer the question if the \emph{k}th cluster is connected on a gene expression
#' level (in addition to the samples being similar based on the other data source). The
#' Within-Cluster CS for a cluster is computed as following:
#' \enumerate{
#' \item Repeatedly for the \emph{i}th sample in the \emph{k}th cluster, apply CSMFA with:
#' \itemize{
#' \item \emph{Query Set}: All cluster samples \strong{excluding} the \emph{i}th sample.
#' \item \emph{Reference}: All samples including the \emph{i}th sample of the \emph{k}th cluster.
#' \item Retrieve the CS of the \emph{i}th sample in the cluster.
#' }
#' \item The Within-Cluster CS for cluster \emph{k} is now defined as the average of all retrieved CS.
#' }
#' The concept of this metric is to investigate the connectivity for each compound with
#' the cluster. The average of the 'leave-one-out' connectivity scores, the Within-Cluster
#' CS, gives an indication of the gene expression connectivity of this cluster. A high
#' Within-Cluster CS implies that the cluster is both similar on the external data source
#' and on the gene expression level. A low score indicates that the cluster does not share a
#' similar latent gene profile structure.
#'
#' \emph{Between-Cluster CS}\cr
#' In this stage of the analysis, we focus on the \emph{l}th cluster and use all compounds in this
#' cluster as the query set. A CSMFA is performed in which all other clusters are the
#' reference set. Next, the connectivity scores are calculated for all reference compounds
#' and averaged over the clusters (=the between connectivity score).
#' A high Between-Cluster CS between the \emph{l}th and \emph{j}th clusters implies that, while the two
#' clusters are not similar based on the other data source, they do share a latent structure
#' when considering the gene expression data.
#'     
#' @author Ewoud De Troyer
#' @param data A gene expression matrix with the compounds in the columns.
#' @param clusterlabels A vector of integers that represents the cluster grouping of the columns (compounds) in \code{data}. The labels should be integers starting from 1 to the total number of clusters. (e.g. the output of \code{\link[stats]{cutree}})
#' @param type Type of CS anaylsis (default=\code{"CSmfa"}): 
#' \itemize{
#' \item \code{"CSmfa"} (MFA or PCA)
#' \item \code{"CSsmfa"} (Sparse MFA or Sparse PCA)
#' \item \code{"CSfabia"} (Fabia)
#' \item \code{"CSzhang"} (Zhang and Gant)
#' }
#' In the first two options, either MFA or PCA is used depending on the cluster size. If the query set only contains a single compound, the latter is used. 
#' Also note that if a cluster only contains a single compound, no \emph{Within-CS} can be computed.
#' 
#' @param WithinABS Boolean value to take the mean of the absolute values in the final step of the \emph{Within-Cluster CS} (default=\code{TRUE}).
#' @param BetweenABS Boolean value to take the mean of the absolute values in the final step of the \emph{Between-Cluster CS} (default=\code{TRUE}).
#' @param FactorABS Boolean value to take the absolute value of the query loadings when determining the best factor (= factor with highest query loadings) in a \code{\link{CSanalysis}} application (default=\code{FALSE}). This option might be helpful if the `best factor` contains large positive and negative query loading which would average to zero.
#' @param verbose Boolean value to output warnings and information about which factor is chosen in a CS analysis (if applicable).
#' @param Within A vector for which cluster numbers the \emph{Within-Cluster CS} should be computed. By default (=\code{NULL}) all within-cluster scores are computed, but this might not be feasible for larger data in which a single \code{\link{CSanalysis}} run might already take a sufficient amount of computation time.
#' @param Between A vector fir which cluster numbers the \emph{Beween-Cluster CS} (with the cluster as a query set) should be computed. By default (=\code{NULL}) all between-cluster scores are computed, but this might not be feasible for larger data in which a single \code{\link{CSanalysis}} run might already take a sufficient amount of computation time.
#' @param WithinSave Boolean value to save the \code{Within} object in the \code{Save} slot of the returned list (default=\code{FALSE}).
#' @param BetweenSave Boolean value to save the \code{Between} object in the \code{Save} slot of the returned list (default=\code{TRUE}).
#' @param ... Additional parameters given to \code{\link{CSanalysis}} specific to a certain \code{type} of CS analysis.
#' @return A list object with components:
#' \itemize{
#' \item{\code{CSmatrix}: }{A \emph{K}\eqn{\times}\emph{K} matrix containing the Within scores on the diagonal and the Between scores elsewhere with the rows being the query set clusters (e.g. \eqn{m_{13}=} Between CS between cluster 1 (as query set) and cluster 3).}
#' \item{\code{CSRankmatrix}: }{The same as \code{CSmatrix}, but with connectivity ranking scores (if applicable).}
#' \item{\code{clusterlabels}: }{The provided \code{clusterlabels}}
#' \item{\code{Save}: }{A list with components:}
#' \itemize{
#' \item{\code{Within}: }{A list with a component for each cluster \emph{k} that contains:}
#' \itemize{
#' \item{\code{LeaveOneOutCS}: }{Each leave-one-out connectivity score for cluster \emph{k}.}
#' \item{\code{LeaveOneOutCSRank}: }{Each leave-one-out connectivity ranking score for cluster \emph{k} (if applicable).}
#' \item{\code{factorselect}: }{A vector containing which factors/BCs were selected in each leave-one-out CS analysis (if applicable).}
#' \item{\code{CS}: }{A (columns (compounds) \eqn{\times} size of cluster \emph{k}) matrix that contains all the connectivity scores in a leave-one-out CS analysis for each left out compound.}
#' \item{\code{CSRank}: }{The same as \code{CS}, but with connectivity ranking scores (if applicable).}
#' }
#' \item{\code{Between}: }{List:}
#' \itemize{
#' \item{\code{DataBetweenCS}: }{A (columns (compounds) \eqn{\times} clusters) matrix containing all compound connectivity scores for each query cluster set.}
#' \item{\code{DataBetweenCSRank}: }{The same as \code{DataBetweenCS}, but with connectivity ranking scores (if applicable).}
#' \item{\code{queryindex}: }{The column indices for each query set in all CS analyses.}
#' \item{\code{factorselect}: }{A vector containing which factors/BCs were selected in each CS analysis (if applicable).}
#' }
#' }
#' }
#' @examples 
#' \dontshow{
#'   # Example Data Set
#'   data("dataSIM",package="CSFA")
#'   # Remove some no-connectivity compounds
#'   nosignal <- sapply(colnames(dataSIM),FUN=function(x){grepl("c-",x)})
#'   data <- dataSIM[,-which(nosignal)[1:250]]
#'   
#'   # Toy example with random cluster assignment:
#'   clusterlabels <- sample(1:10,size=ncol(data),replace=TRUE)
#'   
#'   result2 <- CScluster(data,clusterlabels,type="CSzhang",Within=c(),Between=1)
#'   
#' } 
#' \donttest{
#'   # Example Data Set
#'   data("dataSIM",package="CSFA")
#'   # Remove some no-connectivity compounds
#'   nosignal <- sapply(colnames(dataSIM),FUN=function(x){grepl("c-",x)})
#'   data <- dataSIM[,-which(nosignal)[1:250]]
#'   
#'   # Toy example with random cluster assignment:
#'   # Note: clusterlabels can be acquired through cutree(hclust(...))
#'   clusterlabels <- sample(1:10,size=ncol(data),replace=TRUE)
#'   
#'   result1 <- CScluster(data,clusterlabels,type="CSmfa")
#'   result2 <- CScluster(data,clusterlabels,type="CSzhang")
#'   
#'   result1$CSmatrix
#'   result1$CSRankmatrix
#'   
#'   result2$CSmatrix
#' }
#' @export
CScluster <- function(data,clusterlabels,type="CSmfa",WithinABS=TRUE,BetweenABS=TRUE,FactorABS=FALSE,verbose=FALSE,
                      Within=NULL,Between=NULL,WithinSave=FALSE,BetweenSave=TRUE,...){
  

  
  if(!(type %in% c("CSmfa","CSsmfa","CSfabia","CSzhang"))){stop("type not available. Please choose one of: \"CSmfa\", \"CSsmfa\", \"CSfabia\" or \"CSzhang\".")}
  
  if(length(clusterlabels)!=ncol(data)){stop("Length of clusterlabels is different than the number of columns in data.")}
  if(is.null(rownames(data))){rownames(data) <- paste0("Row",1:nrow(data))}
  if(is.null(colnames(data))){colnames(data) <- paste0("Col",1:ncol(data))}
  
  nclusters <- length(unique(clusterlabels))
  
  if(!identical(1:nclusters,as.integer(sort(unique(clusterlabels))))){stop(paste0("The clusterlabels are not a sequence of integers from 1 to ",nclusters))}
  
  pm <- match.call()
  if(is.null(Within) & deparse(pm$Within)!="c()"){
    Within <- 1:nclusters
  }
  if(is.null(Between) & deparse(pm$Between)!="c()"){
    Between <- 1:nclusters
  }
  
  if(!all(Within<=nclusters & Within>=1)){stop(paste0("Within contains a cluster number not in the range of [1,",nclusters,"]."))}
  if(!all(Between<=nclusters & Between>=1)){stop(paste0("Between contains a cluster number not in the range of [1,",nclusters,"]."))}
  
  
  CSmatrix <- CSRankmatrix <- matrix(NA,nrow=nclusters,ncol=nclusters,dimnames=list(paste0("C",1:nclusters),paste0("C",1:nclusters)))

  if(WithinSave){
    Save.Within <- vector("list",nclusters)
    names(Save.Within) <- paste0("C",1:nclusters)
  }else{
    Save.Within <- NULL
  }
  if(BetweenSave){
    Save.Between <- NULL
    DataBetweenCS <- DataBetweenCSRank <- matrix(NA,nrow=ncol(data),ncol=nclusters,dimnames=list(colnames(data),paste0("C",1:nclusters)))
    refindex_between <- vector("list",nclusters)
    names(refindex_between) <- colnames(DataBetweenCS)
    factorselect_between <- rep(NA,nclusters)
  }else{
    Save.Between <- NULL
  }
  
  
  

  for(i in 1:nclusters){
    
    
    
    i.cluster <- i
    other.cluster <- (1:nclusters)[-i]
    
    cluster.index <- which(clusterlabels==i.cluster)
    

    
    
    ## COMPUTE WITHIN SCORES ##
    
    if(i %in% Within){
      
      # Withing-CS = NA when only 1 compound
      if(length(cluster.index)==1){
        
        if(verbose){
          warning(paste0("Cluster ",i," only has a single compound. No Within-CS will be computed"))
        }
        
        
      }else{
        
        
        extradata1.temp <- matrix(0,nrow=dim(data)[2],ncol=length(cluster.index),dimnames=list(colnames(data),colnames(data)[cluster.index]))
        extradata2.temp <- matrix(NA,nrow=dim(data)[2],ncol=length(cluster.index),dimnames=list(colnames(data),colnames(data)[cluster.index]))
        
        CompoundCS.FactorSelect <-  CSlist1.temp <- CSlist2.temp  <- rep(NA,length(cluster.index))
        
        
        for(j in 1:length(cluster.index)){
          
          original.colindex <- 1:ncol(data)
          original.colindex <- c(original.colindex[cluster.index[-j]],original.colindex[cluster.index[j]],original.colindex[-cluster.index])
          
          querdata <- data[,cluster.index[-j],drop=FALSE]
          refdata <- cbind(data[,cluster.index[j],drop=FALSE],data[,-cluster.index])
          
          
          if(type=="CSmfa"){
            method.type <- ifelse(dim(querdata)[2]==1,"CSpca","CSmfa")
          }else if(type=="CSsmfa"){
            method.type <- ifelse(dim(querdata)[2]==1,"CSspca","CSsmfa")
          }else{
            method.type <- type
          }
          
          
          # ZHANG AND GANT
          if(method.type=="CSzhang"){
            out_ZG <- CSanalysis(querMat=querdata,refMat=refdata,method.type,which=c(),...)
            CSlist1.temp[j] <- out_ZG@CS$CS.ref[1,1]
            
          # FACTOR ANALYSIS METHODS    
          }else{
            temp_cap <- capture.output({out_MFA <- CSanalysis(querMat=querdata,refMat=refdata,method.type,which=c(),component.plot=1,...)})
            
            
            # Check if Factor 1 is best
            factor.select <- best.factor(out_MFA,FactorABS=FactorABS)
            CompoundCS.FactorSelect[j] <- factor.select
            
            if(factor.select!=1){
              
              # This redo currently does not work in CSFA. Alternatively just grab the values from $object
              temp_cap <- capture.output({out_MFA <- CSanalysis(querMat=querdata,refMat=refdata,method.type,which=c(),component.plot=factor.select,result.available=out_MFA,result.available.update=TRUE,...)})
              
              if(verbose){
                warning(paste0("MFA for Compound nr ",colnames(data)[cluster.index[j]]," in Cluster ",i,": Factor ",factor.select," was choosen."),call.=FALSE)
                
              }
            }
            
            
            CSlist1.temp[j] <- out_MFA@CS[[1]]$CS.ref[1,1]
            CSlist2.temp[j] <- out_MFA@CS[[1]]$CS.ref[1,2]
            
          }# end of factor analysis specific part
          
          
         
          if(WithinSave){
            
            if(method.type=="CSzhang"){
              extradata1.temp[,j] <- c(rep(NA,ncol(querdata)),out_ZG@CS$CS.ref$ZGscore)[order(original.colindex)]
              extradata2.temp <- NULL
            }else{
              extradata1.temp[,j] <- get.loadings(out_MFA)[order(original.colindex),factor.select]
              extradata2.temp[,j] <- c(rep(NA,ncol(querdata)),out_MFA@CS[[1]]$CS.ref$CRankScores)[order(original.colindex)]
            }
            
          }else{
            extradata1.temp <- NULL
            extradata2.temp <- NULL
          }
          
        } # end of j loop
        names(CSlist1.temp) <- names(CSlist2.temp) <- colnames(data)[cluster.index]
        
        
        if(WithinABS){
          CSmatrix[i,i] <- mean(abs(CSlist1.temp))  # Average CS is based on absolute values!
          if(method.type!="CSzhang"){CSRankmatrix[i,i] <- mean(abs(CSlist2.temp))}
        }else{
          CSmatrix[i,i] <- mean((CSlist1.temp))  # Average CS is NOT based on absolute values!
          if(method.type!="CSzhang"){CSRankmatrix[i,i] <- mean((CSlist2.temp))}
        }
        
        
        if(WithinSave){
          Save.Within[[i]] <- list(
            LeaveOneOutCS=CSlist1.temp,
            LeaveOneOutCSRank=CSlist2.temp,
            factorselect=CompoundCS.FactorSelect,
            CS=extradata1.temp,
            CSRank=extradata2.temp
          )
        }
      }
      
      
      
    }

    


    ## COMPUTE BETWEEN SCORES ##
    
    if(!(i %in% Between)){next}
    
    if(type=="CSmfa"){
      method.type <- ifelse(length(cluster.index)==1,"CSpca","CSmfa")
    }else if(type=="CSsmfa"){
      method.type <- ifelse(length(cluster.index)==1,"CSspca","CSsmfa")
    }else{
      method.type <- type
    }
    
    
    original.colindex <- 1:ncol(data)
    original.colindex <- c(original.colindex[cluster.index],original.colindex[-cluster.index])
    
    
    
    # ZHANG AND GANT
    if(method.type=="CSzhang"){
      
      out_ZG <- CSanalysis(querMat=data[,cluster.index,drop=FALSE],refMat=data[,-cluster.index],method.type,which=c(),...)
      
      temp_CS <- c(rep(NA,length(cluster.index)),out_ZG@CS$CS.ref$ZGscore)[order(original.colindex)]
      
      
    # FACTOR ANALYSIS METHODS  
    }else{
      
      temp_cap <- capture.output({out_MFA <- CSanalysis(querMat=data[,cluster.index,drop=FALSE],refMat=data[,-cluster.index],method.type,which=c(),component.plot=1,...)})
      
      # Check if Factor 1 is best
      factor.select <- best.factor(out_MFA,FactorABS=FactorABS)
      factorselect_between[i] <- factor.select
      
      if(factor.select!=1){
        temp_cap <- capture.output({out_MFA <- CSanalysis(querMat=data[,cluster.index,drop=FALSE],refMat=data[,-cluster.index],method.type,which=c(),component.plot=factor.select,result.available=out_MFA,result.available.update=TRUE,...)})
        
        if(verbose){
          warning(paste0("MFA for Cluster ",i,": Factor ",factor.select," was choosen."),call.=FALSE)
        }
      }
      
      temp_CS <- get.loadings(out_MFA)[order(original.colindex),factor.select]
      temp_CSRank <- c(rep(NA,length(cluster.index)),out_MFA@CS[[1]]$CS.ref$CRankScores)[order(original.colindex)]
    }
    
    

    if(BetweenSave){
      DataBetweenCS[,i] <- temp_CS
      if(method.type!="CSzhang"){DataBetweenCSRank[,i] <- temp_CSRank}
     
      refindex_between[[i]] <- cluster.index
      
    }
    
    


    if(BetweenABS){
      for(ii in 1:length(other.cluster)){
        ii.cluster <- other.cluster[[ii]]
        CSmatrix[i,ii.cluster] <- mean(abs(temp_CS[which(clusterlabels==ii.cluster)]))
        if(method.type!="CSzhang"){CSRankmatrix[i,ii.cluster] <- mean(abs(temp_CSRank[which(clusterlabels==ii.cluster)]))}
      }
      
    }else{
      for(ii in 1:length(other.cluster)){
        ii.cluster <- other.cluster[[ii]]
        CSmatrix[i,ii.cluster] <- mean((temp_CS[which(clusterlabels==ii.cluster)]))
        if(method.type!="CSzhang"){CSRankmatrix[i,ii.cluster] <- mean((temp_CSRank[which(clusterlabels==ii.cluster)]))}
      }
    }
  }
  

  if(BetweenSave){ 
    Save.Between <- list(
      DataBetweenCS=DataBetweenCS,
      DataBetweenCSRank=DataBetweenCSRank,
      queryindex=refindex_between,
      factorselect=factorselect_between
    )
  }
  
  

  # New out
  out <- list(
    CSmatrix=CSmatrix,
    CSRankmatrix=CSRankmatrix,
    clusterlabels=clusterlabels,
    Save=list(
      Within=Save.Within,
      Beween=Save.Between
    )
  )
  
  return(out)
}



# Function to check if Factor 1 was ok
best.factor <- function(CSresult,FactorABS=FALSE){
  
  
  loadings_temp <- get.loadings(CSresult)
  
  
  if(FactorABS){
    mean.quer.loadings <- abs(colMeans(abs(loadings_temp[c(1:nrow(CSresult@CS[[1]]$CS.query)),])))
    
  }else{
    mean.quer.loadings <- abs(colMeans(loadings_temp[c(1:nrow(CSresult@CS[[1]]$CS.query)),]))
  }
  
  return(which(mean.quer.loadings==max(mean.quer.loadings)))
}


# 
# 
# # Plot function which takes:
# # - ClusterCS
# # - which cluster
# # - which compound (NA when compare with all)
# # - ordered or not
# # - CS or CSRank
# 
# # To do: CLusterCS needs to save more information for some of these compound != NA
# plot.ClusterCS <- function(ClusterCS,cluster,compound=NA,type="CS",colors=NULL){
#   
#   
#   
#   
#   if(!(type %in% c("CS","CSRank"))){stop("type parameter incorrect",call.=FALSE)}
#   
#   if(type=="CS"){
#     plotdata <- ClusterCS$CS
#   }
#   if(type=="CSRank"){
#     plotdata <- ClusterCS$CSRank
#   }
#   
#   nr.clusters <- length(plotdata)-1
#   if(!(cluster %in% c(1:nr.clusters))){stop("cluster parameter incorrect",call.=FALSE)}
#   
#   # General Plot (Cluster vs other Clusters)
#   if(is.na(compound)){
#     
#     #		devtools::install_git("https://github.com/ronammar/randomcoloR")
#     if(is.null(colors)){
#       require(randomcoloR)
#       set.seed(1)
#       col.palette <- distinctColorPalette(nr.clusters)
#     }else{
#       col.palette <- colors
#     }
#     
#     
#     plotdata_cluster <- plotdata[[cluster]]
#     
#     if(type=="CS"){
#       d1 <- plotdata_cluster$CStop[order(plotdata_cluster$CStop$Cluster),]
#       d2 <- plotdata_cluster$RefLoadings
#       
#       col1 <- c(rep("blue",length(d2)),rep("black",dim(d1)[1]))
#       colbg <- c(rep("blue",length(d2)),col.palette[d1$Cluster])
#       
#       d3 <- c(d2,d1$Score)
#       names(d3) <- c(names(d2),rownames(d1))
#       ylab.temp <- "CS"
#       
#     }
#     if(type=="CSRank"){
#       d1 <- plotdata_cluster$CStop[order(plotdata_cluster$CStop$Cluster),]
#       col1 <- c(rep("black",dim(d1)[1]))
#       colbg <- c(col.palette[d1$Cluster])
#       d3 <- d1$Score
#       names(d3) <- rownames(d1)
#       
#       ylab.temp <- "CS RankScores"
#     }
#     
#     plot(d3,pch=21,col=col1,bg=colbg,xlab="Compounds",ylab=ylab.temp,main=paste0("Cluster ",names(plotdata)[cluster]," vs Other Clusters"))
#     text(d3,names(d3),col=colbg,pos=3)
#     abline(h=0,lty=2)
#     
#     if(type=="CS"){legend("topright",c(paste0("Query Cluster ",names(plotdata)[cluster]),"Other Clusters"),pch=21,col=c("blue","black"),pt.bg=c("blue","white"),bty="n")}
#     if(type=="CSRank"){legend("topright",c("Other Clusters"),pch=21,col=c("black"),pt.bg=c("white"),bty="n")}
#     
#   }else{# Inside-Cluster plot (all-1 vs 1)
#     
#     
#     
#     plotdata_cluster <- plotdata[[cluster]]
#     if(!(compound %in% c(1:length(plotdata_cluster$CompoundCS)))){stop("compound parameter incorrect",call.=FALSE)}
#     
#     loadings <- plotdata_cluster$CompoundCS.extradata[,compound]
#     
#     ref.index <- sapply(names(plotdata_cluster$CompoundCS[-compound]),FUN=function(x){which(x==rownames(ClusterCS$CS[[1]]$CompoundCS.extradata))}) # take first cluster as example, because same for all
#     quer1.index <- which(names(plotdata_cluster$CompoundCS[compound])==rownames(ClusterCS$CS[[1]]$CompoundCS.extradata))
#     # check if rownames same in all clusters -> okay
#     
#     # Put Color
#     col.loadings <- rep("black",length(loadings))
#     col2.loadings <- rep("grey",length(loadings))
#     col.loadings[ref.index] <- col2.loadings[ref.index] <- "blue"
#     col.loadings[quer1.index] <- col2.loadings[quer1.index] <- "red"
#     
#     # Make Data
#     if(type=="CS"){
#       d <- c(loadings[ref.index],loadings[-ref.index])
#       col1 <- c(col.loadings[ref.index],col.loadings[-ref.index])
#       col2 <- c(col2.loadings[ref.index],col2.loadings[-ref.index])
#       ylab.temp <- "CS"
#     }
#     if(type=="CSRank"){
#       d <- c(loadings[-ref.index])
#       col1 <- c(col.loadings[-ref.index])
#       col2 <- c(col2.loadings[-ref.index])
#       ylab.temp <- "CS RankScores"
#     }
#     
#     
#     
#     
#     plot(d,pch=21,bg=col2,col=col1,ylab=ylab.temp,xlab="Compounds",main=paste0("Cluster ",names(plotdata)[cluster]," - Compound \"",names(plotdata_cluster$CompoundCS)[compound],"\""))
#     text(d,labels=names(d),col=col1,pos=3)	
#     abline(h=0,lty=2)
#     
#     
#     if(type=="CS"){legend("topright",c("Query","Leave-one-out Cmpd"),pch=21,pt.bg=c("blue","red"),col=c("blue","red"),bty="n")}
#     if(type=="CSRank"){legend("topright",c("Leave-one-out Cmpd"),pch=21,pt.bg="red",col="red",bty="n")}
#     
#   }
#   
#   
#   
# }
