
## EXAMPLE
# load("C:/Users/lucp8394/Documents/StatFiles/CenStat/Connectivity Score/Research/ClusterCS/MCF7Data.RData")
# load("C:/Users/lucp8394/Documents/StatFiles/CenStat/Connectivity Score/Research/ClusterCS/cluster_target.RData")
# data <- geneMat
# clusterlabels <- cutree(cluster_target$Clust,k=7)










## STILL NEED TO INCLUDE PCA FOR WHEN 1 COMPOUND
## Currently: computes all withing and all between scores -> maybe problematic for larger data

# TO DO: change appending + what if cluster is only 1 compound -> NA withing CS


# Disclaimers:
# - factor selection is done through biggest mean. So it could be issue when query set contains both high and negative together -> add abs option (WITHIN-CS), analysis in chapter is with abs=FALSE



ClusterCS <- function(data,clusterlabels,WithinABS=FALSE,verbose=TRUE){
  
  
  # Add some stuff before such as: add row/column names if not in data
  
  
  nclusters <- length(unique(clusterlabels))
  
  
  CSmatrix1 <- matrix(0,nrow=nclusters,ncol=nclusters,dimnames=list(paste0("C",unique(clusterlabels)),paste0("C",unique(clusterlabels))))
  CSlist1 <- rep( list(list()), nclusters) 
  
  CSmatrix2 <- matrix(0,nrow=nclusters,ncol=nclusters,dimnames=list(paste0("C",unique(clusterlabels)),paste0("C",unique(clusterlabels))))
  CSlist2 <- rep( list(list()), nclusters ) 
  
  for(i in 1:length(unique(clusterlabels))){
    
    
    
    i.cluster <- unique(clusterlabels)[i]
    other.cluster <- unique(clusterlabels)[-i]
    
    cluster.index <- which(clusterlabels==i.cluster)
    
    CSlist1.temp <- CSlist2.temp  <- vector("list",length(cluster.index))

    
    extradata1.temp <- matrix(0,nrow=dim(data)[2],ncol=length(cluster.index),dimnames=list(colnames(data),colnames(data)[cluster.index]))
    extradata2.temp <- matrix(NA,nrow=dim(data)[2],ncol=length(cluster.index),dimnames=list(colnames(data),colnames(data)[cluster.index]))
    
    CompoundCS.FactorSelect <- rep(NA,length(cluster.index))
    
    ## COMPUTE WITHIN SCORES
    
    # TO DO !!!!
    
    # # Withing-CS = NA when only 1 compound
    # if(length(cluster.index)==1){
    #   
    # }else{
    #   
    # }
    
    for(j in 1:length(cluster.index)){
      
      querdata <- data[,cluster.index[-j],drop=FALSE]
      refdata <- cbind(data[,cluster.index[j],drop=FALSE],data[,-cluster.index])
      
      method.type <- ifelse(dim(querdata)[2]==1,"CSpca","CSmfa")
      
      out_MFA <- CSanalysis(querMat=querdata,refMat=refdata,method.type,which=c(),component.plot=1)
      
      # Check if Factor 1 is best
      factor.select <- best.factor(out_MFA)
      CompoundCS.FactorSelect[j] <- factor.select
      
      if(factor.select!=1){
        
        # This redo currently does not work in CSFA. Alternatively just grab the values from $object
        out_MFA <- CSanalysis(querMat=querdata,refMat=refdata,method.type,which=c(),component.plot=factor.select,result.available=out_MFA,result.available.update=TRUE)
        
        if(verbose){
          warning(paste0("MFA for Compound nr ",j," in Cluster nr ",i,": Factor ",factor.select," was choosen."),call.=FALSE)
          
        }
      }
      
      # I AM HERE
      
      CSlist1.temp <- c(CSlist1.temp,out_MFA@CS$CS.query[1,1])
      CSlist2.temp <- c(CSlist2.temp,out_MFA@CSRankScores[[1]][1,1])
      
      extradata1.temp[,j] <- rbind(out_MFA@CS$CS.ref[,1,drop=FALSE],out_MFA@CS$CS.query[,1,drop=FALSE])[rownames(extradata1.temp),] 
      
      extra.temp <- out_MFA@CSRankScores[[1]][,1,drop=FALSE]
      extradata2.temp[intersect(rownames(extradata2.temp),rownames(extra.temp)),j] <- extra.temp[intersect(rownames(extradata2.temp),rownames(extra.temp)),]
      
    }
    names(CSlist1.temp) <- names(CSlist2.temp) <- colnames(data)[cluster.index]
    
    CSlist1[[i]]$CompoundCS <- CSlist1.temp
    CSlist2[[i]]$CompoundCS <- CSlist2.temp
    
    CSlist1[[i]]$CompoundCS.FactorSelect <- CompoundCS.FactorSelect
    CSlist2[[i]]$CompoundCS.FactorSelect <- CompoundCS.FactorSelect
    
    # To do:  (also something similar down below)
    CSlist1[[i]]$CompoundCS.extradata <- extradata1.temp
    CSlist2[[i]]$CompoundCS.extradata <- extradata2.temp
    
    CSmatrix1[i,i] <- mean(abs(CSlist1.temp))  # Average CS is based on absolute values!
    CSmatrix2[i,i] <- mean(abs(CSlist2.temp))
    
    
    method.type <- ifelse(length(cluster.index)==1,"CSpca","CSmfa")
    
    # COMPUTE BETWEEN SCORES
    out_MFA <- CSanalysis(querMat=data[,cluster.index,drop=FALSE],refMat=data[,-cluster.index],method.type,which=c(),component.plot=1)
    
    # Check if Factor 1 is best
    factor.select <- best.factor(out_MFA)
    
    if(factor.select!=1){
      out_MFA <- CSanalysis(querMat=data[,cluster.index,drop=FALSE],refMat=data[,-cluster.index],method.type,which=c(),component.plot=factor.select,result.available=out_MFA)
      
      warning(paste0("MFA for Cluster nr ",i,": Factor ",factor.select," was choosen."),call.=FALSE)
    }
    
    CStop1.temp <- data.frame(Score=out_MFA@CS$CS.query[,1],Cluster=clusterlabels[-cluster.index])
    CStop2.temp <- data.frame(Score=out_MFA@CSRankScores[[1]][,1],Cluster=clusterlabels[-cluster.index])
    rownames(CStop1.temp) <- rownames(CStop2.temp) <- colnames(data)[-cluster.index]
    
    CStop1.temp <- CStop1.temp[order(abs(CStop1.temp$Score),decreasing=TRUE),]
    CStop2.temp <- CStop2.temp[order(abs(CStop2.temp$Score),decreasing=TRUE),]
    
    CSlist1[[i]]$CStop <- CStop1.temp
    CSlist2[[i]]$CStop <- CStop2.temp
    
    CSlist1[[i]]$RefLoadings <- out_MFA@CS$CS.ref[,1]
    CSlist2[[i]]$RefLoadings <- out_MFA@CS$CS.ref[,1]
    
    CSlist1[[i]]$Factor.Select <- factor.select
    CSlist2[[i]]$Factor.Select <- factor.select
    
    names(CSlist1[[i]]$RefLoadings) <- names(CSlist2[[i]]$RefLoadings) <- colnames(data)[cluster.index]
    
    for(ii in 1:length(other.cluster)){
      ii.cluster <- other.cluster[ii]
      
      CSmatrix1[i,which(unique(clusterlabels)==ii.cluster)] <- mean(abs(CStop1.temp[CStop1.temp$Cluster==ii.cluster,1])) # Average CS is based on absolute values!
      CSmatrix2[i,which(unique(clusterlabels)==ii.cluster)] <- mean(abs(CStop2.temp[CStop2.temp$Cluster==ii.cluster,1]))
      
      
    }
    
    
    
  }
  
  # don't forget to name CSlist1 en 2
  names(CSlist1) <- paste0("C",unique(clusterlabels))
  names(CSlist2) <- paste0("C",unique(clusterlabels))
  
  
  CSlist1$CSmatrix <- CSmatrix1
  CSlist2$CSmatrix <- CSmatrix2
  
  out <- list(CS=CSlist1,CSRank=CSlist2)
  return(out)
  
}



# Function to check if Factor 1 was ok
best.factor <- function(CSresult){
  mean.quer.loadings <- abs(colMeans(CSresult@extra$object$quanti.var$cor[c(1:nrow(CSresult@CS[[1]]$CS.query)),]))
  return(which(mean.quer.loadings==max(mean.quer.loadings)))
}




# Plot function which takes:
# - ClusterCS
# - which cluster
# - which compound (NA when compare with all)
# - ordered or not
# - CS or CSRank

# To do: CLusterCS needs to save more information for some of these compound != NA
plot.ClusterCS <- function(ClusterCS,cluster,compound=NA,type="CS",colors=NULL){
  
  if(!(type %in% c("CS","CSRank"))){stop("type parameter incorrect",call.=FALSE)}
  
  if(type=="CS"){
    plotdata <- ClusterCS$CS
  }
  if(type=="CSRank"){
    plotdata <- ClusterCS$CSRank
  }
  
  nr.clusters <- length(plotdata)-1
  if(!(cluster %in% c(1:nr.clusters))){stop("cluster parameter incorrect",call.=FALSE)}
  
  # General Plot (Cluster vs other Clusters)
  if(is.na(compound)){
    
    #		devtools::install_git("https://github.com/ronammar/randomcoloR")
    if(is.null(colors)){
      require(randomcoloR)
      set.seed(1)
      col.palette <- distinctColorPalette(nr.clusters)
    }else{
      col.palette <- colors
    }
    
    
    plotdata_cluster <- plotdata[[cluster]]
    
    if(type=="CS"){
      d1 <- plotdata_cluster$CStop[order(plotdata_cluster$CStop$Cluster),]
      d2 <- plotdata_cluster$RefLoadings
      
      col1 <- c(rep("blue",length(d2)),rep("black",dim(d1)[1]))
      colbg <- c(rep("blue",length(d2)),col.palette[d1$Cluster])
      
      d3 <- c(d2,d1$Score)
      names(d3) <- c(names(d2),rownames(d1))
      ylab.temp <- "CS"
      
    }
    if(type=="CSRank"){
      d1 <- plotdata_cluster$CStop[order(plotdata_cluster$CStop$Cluster),]
      col1 <- c(rep("black",dim(d1)[1]))
      colbg <- c(col.palette[d1$Cluster])
      d3 <- d1$Score
      names(d3) <- rownames(d1)
      
      ylab.temp <- "CS RankScores"
    }
    
    plot(d3,pch=21,col=col1,bg=colbg,xlab="Compounds",ylab=ylab.temp,main=paste0("Cluster ",names(plotdata)[cluster]," vs Other Clusters"))
    text(d3,names(d3),col=colbg,pos=3)
    abline(h=0,lty=2)
    
    if(type=="CS"){legend("topright",c(paste0("Query Cluster ",names(plotdata)[cluster]),"Other Clusters"),pch=21,col=c("blue","black"),pt.bg=c("blue","white"),bty="n")}
    if(type=="CSRank"){legend("topright",c("Other Clusters"),pch=21,col=c("black"),pt.bg=c("white"),bty="n")}
    
  }else{# Inside-Cluster plot (all-1 vs 1)
    
    
    
    plotdata_cluster <- plotdata[[cluster]]
    if(!(compound %in% c(1:length(plotdata_cluster$CompoundCS)))){stop("compound parameter incorrect",call.=FALSE)}
    
    loadings <- plotdata_cluster$CompoundCS.extradata[,compound]
    
    ref.index <- sapply(names(plotdata_cluster$CompoundCS[-compound]),FUN=function(x){which(x==rownames(ClusterCS$CS[[1]]$CompoundCS.extradata))}) # take first cluster as example, because same for all
    quer1.index <- which(names(plotdata_cluster$CompoundCS[compound])==rownames(ClusterCS$CS[[1]]$CompoundCS.extradata))
    # check if rownames same in all clusters -> okay
    
    # Put Color
    col.loadings <- rep("black",length(loadings))
    col2.loadings <- rep("grey",length(loadings))
    col.loadings[ref.index] <- col2.loadings[ref.index] <- "blue"
    col.loadings[quer1.index] <- col2.loadings[quer1.index] <- "red"
    
    # Make Data
    if(type=="CS"){
      d <- c(loadings[ref.index],loadings[-ref.index])
      col1 <- c(col.loadings[ref.index],col.loadings[-ref.index])
      col2 <- c(col2.loadings[ref.index],col2.loadings[-ref.index])
      ylab.temp <- "CS"
    }
    if(type=="CSRank"){
      d <- c(loadings[-ref.index])
      col1 <- c(col.loadings[-ref.index])
      col2 <- c(col2.loadings[-ref.index])
      ylab.temp <- "CS RankScores"
    }
    
    
    
    
    plot(d,pch=21,bg=col2,col=col1,ylab=ylab.temp,xlab="Compounds",main=paste0("Cluster ",names(plotdata)[cluster]," - Compound \"",names(plotdata_cluster$CompoundCS)[compound],"\""))
    text(d,labels=names(d),col=col1,pos=3)	
    abline(h=0,lty=2)
    
    
    if(type=="CS"){legend("topright",c("Query","Leave-one-out Cmpd"),pch=21,pt.bg=c("blue","red"),col=c("blue","red"),bty="n")}
    if(type=="CSRank"){legend("topright",c("Leave-one-out Cmpd"),pch=21,pt.bg="red",col="red",bty="n")}
    
  }
  
  
  
}
