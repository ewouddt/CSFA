
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



ClusterCS <- function(data,clusterlabels,WithinABS=TRUE,BetweenABS=TRUE,FactorABS=FALSE,verbose=TRUE,WithinSave=FALSE,BetweenSave=TRUE,
                      Within=NULL,Between=NULL){
  
  # TO DOOOO!!!!!
  # Add some stuff before such as: add row/column names if not in data
  # Check if Within and Between do not contain integers higher or lower than in clusterlabels + still include these options
  # Clusterlabels need to be integers
  # Make sure unique(clusterlabels) is a sequence
  # identical(1:nclusters,sort(unique(clusterlabels)))
  
  
  
  
  nclusters <- length(unique(clusterlabels))
  
  if(is.null(Within)){
    Within <- 1:nclusters
  }
  if(is.null(Between)){
    Between <- 1:nclusters
  }
  
  
  CSmatrix <- CSRankmatrix <- matrix(NA,nrow=nclusters,ncol=nclusters,dimnames=list(paste0("C",sort(unique(clusterlabels))),paste0("C",unique(clusterlabels))))

  if(WithinSave){
    Save.Within <- vector("list",nclusters)
  }else{
    Save.Within <- NULL
  }
  if(BetweenSave){
    Save.Between <- NULL
    DataBetweenCS <- DataBetweenCSRank <- matrix(NA,nrow=nrow(data),ncol=nclusters,dimnames=list(rownames(data),paste0("C",sort(unique(clusterlabels)))))
    refindex_between <- vector("list",nclusters)
    names(refindex_between) <- colnames(DataBetweenCS)
    factorselect_between <- rep(NA,nclusters)
  }else{
    Save.Between <- NULL
  }
  
  
  

  for(i in 1:length(unique(clusterlabels))){
    
    
    
    i.cluster <- unique(clusterlabels)[i]
    other.cluster <- unique(clusterlabels)[-i]
    
    cluster.index <- which(clusterlabels==i.cluster)
    

    
    
    ## COMPUTE WITHIN SCORES ##
    

    # Withing-CS = NA when only 1 compound
    if(length(cluster.index)==1){

      if(verbose){
        warning(paste0("Cluster ",unique(clusterlabels)[i]," only has a single compound. No Within-CS will be computed"))
      }
      
      
    }else{
      
      
      extradata1.temp <- matrix(0,nrow=dim(data)[2],ncol=length(cluster.index),dimnames=list(colnames(data),colnames(data)[cluster.index]))
      extradata2.temp <- matrix(NA,nrow=dim(data)[2],ncol=length(cluster.index),dimnames=list(colnames(data),colnames(data)[cluster.index]))
      
      CompoundCS.FactorSelect <-  CSlist1.temp <- CSlist2.temp  <- rep(NA,length(cluster.index))
      
      
      for(j in 1:length(cluster.index)){
        
        if(!(j %in% Within)){break}
        
        original.colindex <- 1:ncol(data)
        original.colindex <- c(original.colindex[cluster.index[-j]],original.colindex[cluster.index[j]],original.colindex[-cluster.index])
        
        querdata <- data[,cluster.index[-j],drop=FALSE]
        refdata <- cbind(data[,cluster.index[j],drop=FALSE],data[,-cluster.index])
        
        method.type <- ifelse(dim(querdata)[2]==1,"CSpca","CSmfa")
        
        out_MFA <- CSanalysis(querMat=querdata,refMat=refdata,method.type,which=c(),component.plot=1)
        
        # Check if Factor 1 is best
        factor.select <- best.factor(out_MFA,FactorABS=FactorABS)
        CompoundCS.FactorSelect[j] <- factor.select
        
        if(factor.select!=1){
          
          # This redo currently does not work in CSFA. Alternatively just grab the values from $object
          out_MFA <- CSanalysis(querMat=querdata,refMat=refdata,method.type,which=c(),component.plot=factor.select,result.available=out_MFA,result.available.update=TRUE)
          
          if(verbose){
            warning(paste0("MFA for Compound nr ",colnames(data)[cluster.index[j]]," in Cluster nr ",unique(clusterlabels)[i],": Factor ",factor.select," was choosen."),call.=FALSE)
            
          }
        }
        
        
        CSlist1.temp[j] <- out_MFA@CS[[1]]$CS.ref[1,1]
        CSlist2.temp[j] <- out_MFA@CS[[1]]$CS.ref[1,2]
        
        
        if(WithinSave){
          extradata1.temp[,j] <- out_MFA@extra$object$quanti.var$cor[order(original.colindex),factor.select]
          extradata2.temp[,j] <- c(rep(NA,ncol(querdata)),out_MFA@CS[[1]]$CS.ref$CRankScores)[order(original.colindex)]
        }else{
          extradata1.temp <- NULL
          extradata2.temp <- NULL
        }
        
      }
      names(CSlist1.temp) <- names(CSlist2.temp) <- colnames(data)[cluster.index]
      
      
      if(WithinABS){
        CSmatrix[i,i] <- mean(abs(CSlist1.temp))  # Average CS is based on absolute values!
        CSRankmatrix[i,i] <- mean(abs(CSlist2.temp))
      }else{
        CSmatrix[i,i] <- mean((CSlist1.temp))  # Average CS is NOT based on absolute values!
        CSRankmatrix[i,i] <- mean((CSlist2.temp))
      }
      
      
      if(WithinSave){
        Save.Within[[i]] <- list(
          LeaveOneOutCS=CSlist1.temp,
          LeaveOneOutCSRank=CSlist2.temp,
          FactorSelect=CompoundCS.FactorSelect,
          CSMFA=extradata1.temp,
          CSRankMFA=extradata2.temp
        )
      }
    }
    


    ## COMPUTE BETWEEN SCORES ##
    
    if(!(i %in% Between)){break}
    
    method.type <- ifelse(length(cluster.index)==1,"CSpca","CSmfa")
    original.colindex <- 1:ncol(data)
    original.colindex <- c(original.colindex[cluster.index],original.colindex[-cluster.index])
    
    
    
    out_MFA <- CSanalysis(querMat=data[,cluster.index,drop=FALSE],refMat=data[,-cluster.index],method.type,which=c(),component.plot=1)
    
    # Check if Factor 1 is best
    factor.select <- best.factor(out_MFA,FactorABS=FactorABS)
    
    if(factor.select!=1){
      out_MFA <- CSanalysis(querMat=data[,cluster.index,drop=FALSE],refMat=data[,-cluster.index],method.type,which=c(),component.plot=factor.select,result.available=out_MFA,result.available.update=TRUE)
      
      if(verbose){
        warning(paste0("MFA for Cluster nr ",unique(clusterlabels)[i],": Factor ",factor.select," was choosen."),call.=FALSE)
      }
    }
    
    
    # I AM HERE
    
    temp_CS <- out_MFA@extra$object$quanti.var$cor[order(original.colindex),factor.select]
    temp_CSRank <- c(rep(NA,length(cluster.index)),out_MFA@CS[[1]]$CS.ref$CRankScores)[order(original.colindex)] 
    
    if(BetweenSave){
      DataBetweenCS[,i] <- temp_CS
      DataBetweenCSRank[,i] <- temp_CSRank
     
      refindex_between[[i]] <- cluster.index
      factorselect_between[i] <- factor.select
    }
    
    


    if(BetweenABS){
      for(ii in 1:length(other.cluster)){
        ii.cluster <- other.cluster[[ii]]
        CSmatrix[i,ii.cluster] <- mean(abs(temp_CS[which(clusterlabels==ii.cluster)]))
        CSRankmatrix[i,ii.cluster] <- mean(abs(temp_CSRank[which(clusterlabels==ii.cluster)]))
      }
      
    }else{
      for(ii in 1:length(other.cluster)){
        ii.cluster <- other.cluster[[ii]]
        CSmatrix[i,ii.cluster] <- mean((temp_CS[which(clusterlabels==ii.cluster)]))
        CSRankmatrix[i,ii.cluster] <- mean((temp_CSRank[which(clusterlabels==ii.cluster)]))
      }
    }
  }
  

  if(BetweenSave){ 
    Save.Between <- list(
      DataBetweenCS=DataBetweenCS,
      DataBetweenCSRank=DataBetweenCSRank,
      refindex=refindex_between,
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
  
  if(FactorABS){
    mean.quer.loadings <- abs(colMeans(abs(CSresult@extra$object$quanti.var$cor[c(1:nrow(CSresult@CS[[1]]$CS.query)),])))
    
  }else{
    mean.quer.loadings <- abs(colMeans(CSresult@extra$object$quanti.var$cor[c(1:nrow(CSresult@CS[[1]]$CS.query)),]))
  }
  
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
