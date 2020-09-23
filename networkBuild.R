
library(ggplot2)
library(mice)
library(ggthemes)
library(gplots)

data.gene_expression <- as.matrix(read.table("F:/project/2.13/end_touse/data_gene_expression_NM.txt"))
#data.gene_expression<-apply(data.gene_expression, 1, scale)
colnames <- read.table("F:/project/newkang_brainNet/kang/new_code_document/soft/column_names.txt",header=F)
colnames(data.gene_expression) <- colnames$V1
#data.NM<-read.table("F:/project/newkang_brainNet/kang/new_code_document/soft/NM_name.txt")

data.FG_NM<- read.table("F:/project/2.13/end_touse/FG_NM.txt",header=F,stringsAsFactors = F)
data.NM<-data.FG_NM$V2
factor.data.NM <- factor(data.NM)
data.mean.gene_expression <- apply(data.gene_expression,2,function(x){tapply(x,factor.data.NM,mean)})

data.column <- read.table("F:/project/newkang_brainNet/kang/column.txt",sep="\t",header=T,row.names=1)
data.ncx.left.column<-data.column[which(data.column$Sample_label>=11&data.column$Sample_position=='cortical'&data.column$Sample_characteristics_hemisphere=='L'),]
# data.ncx.left.column<-data.ncx.left.column1[which(data.ncx.left.column1$Sample_characteristics_brain_code!="HSB133"&data.ncx.left.column1$Sample_characteristics_brain_code!="HSB125"),]
# 

order.data.ncx.left.column <- data.ncx.left.column[order(data.ncx.left.column$Sample_label,data.ncx.left.column$Sample_characteristics_region,data.ncx.left.column$Sample_characteristics_hemisphere,data.ncx.left.column$Sample_characteristics_brain_code,decreasing=F),]

data.gene_expression.left.leave<-data.mean.gene_expression[,rownames(order.data.ncx.left.column)]
#colnames(data.gene_expression.left.leave)<-rownames(data.ncx.left.column)
#rownames(data.gene_expression.left.leave)<-rownames(data.mean.gene_expression)

#data.gene_expression.left.leave<-apply(data.gene_expression.left.leave, 1, scale)
#脑网络与基因网络进行构建，构建网络有两种方法#
##gene network build #
geneNetbuild<-function(expressionMatrix,sampleinfo){
  #computer the correlation
  expressionMatrix<-sampledata
  sampleinfo<-sampleinfo
  data.my_cor <- function(x){
    #x<-data.gene.mat
    length <- ncol(x)
    data.cor <- as.vector(matrix(NA,(length*(length-1))/2,1))
    m=1
    for(i in c(1:(length-1))){
    #  i=1
      for(j in c((i+1):length)){
       # j=2
        data <- x[,c(i,j)]
        data.is.na <- apply(data,1,function(x){any(is.na(x))})
        data.na.delete <- data[which(data.is.na==FALSE),]
        if(nrow(data.na.delete)>2){
          data.cor[m] <- cor(data.na.delete[,1],data.na.delete[,2])
        }else{
          data.cor[m] <- cor(data[,1],data[,2])
        }
        m=m+1
      }
    }
    data.cor<-atanh(data.cor)##fisher change, keep the normal distribution
    return(data.cor)
  }
  geneCornetwork<-function(expressionMatrix,sampleinfo){
    #expressionMatrix<-sampledata
    colnames(sampleinfo)<-c('subject','sample')
    #get the subject number info
    subjectNum <- factor(sampleinfo$subject)
    levels.subjectNum.L <- levels(subjectNum)
    num.label.Subject <- length(levels.subjectNum.L)
    
    factor.data.column.label.region.L <- factor(sampleinfo$sample)
    levels.sample.label <- levels(factor.data.column.label.region.L)
    num.label.sample <- length(levels.sample.label)
    
    data.normal.mat.list <- matrix(NA,nrow(expressionMatrix),ncol(expressionMatrix))
    colname<-vector()
    k=1
    for (i in levels.data.column.label.HSB.L){
      ##normaliza the data
      data_num1<-sampleinfo[which(sampleinfo$subject==i),]
      data_num_geneExpression <- expressionMatrix[,rownames(data_num1)]
      #####process the normallization，first step normalize the transcript in 11 brain region
      data_Expression_noramalization1<-t(apply(data_num_geneExpression,2, scale))
      #data_Expression_noramalization22<-scale(data_num_geneExpression)
      data_Expression_noramalization<-t(apply(data_Expression_noramalization1,1, scale))
      data.normal.mat.list[,k:(k+nrow(data_num1)-1)]<-data_Expression_noramalization
      k=k+nrow(data_num1)
    }
    rownames(data.normal.mat.list)<-rownames(expressionMatrix)
    colnames(data.normal.mat.list)<-colnames(expressionMatrix)
    data.label.L<-data.normal.mat.list
    
    data.mean.mat.list <- list()
    data.mean.mat.list <- list()
    
    for(k in c(1:num.label.Subject)){
      # k=1
      data.mean.mat <- matrix(0,nrow(expression),num.label.sample)
      colnames(data.mean.mat) <- levels.sample.label
      rownames(data.mean.mat) <- row.names(expression)
      
      data.column <- sampleinfo[which(sampleinfo$subject==levels.subjectNum.L[k]),]
      data <- as.matrix(expression[,row.names(data.column)])
      
      colnames(data) <- data.column$sample
      
      data.mean.mat.list <- c(data.mean.mat.list, list(data))
    }
    
    label.roi <- levels.data.column.label.region.L
    len.label.roi <- length(label.roi)
    data.cor.mat.gene <- matrix(NA,nrow(data.label.L),(len.label.roi*(len.label.roi-1)/2))
    rownames(data.cor.mat.gene) <- row.names(data.label.L)
    names <- c(1:len.label.roi)
    k=1
    for(i in c(1:(len.label.roi-1))){
      for(j in c((i+1):len.label.roi)){
        names[k] <- paste(label.roi[i],' & ',label.roi[j],sep='');
        k=k+1
      }
    }
    
    label.roi <- levels.sample.label
    len.label.roi <- length(label.roi)
    data.cor.mat.gene <- matrix(NA,nrow(data.label.L),(len.label.roi*(len.label.roi-1)/2))
    rownames(data.cor.mat.gene) <- row.names(data.label.L)
    colnames(data.cor.mat.gene) <- names
    
    for(i in c(1:nrow(data.cor.mat.gene))){
      #i=1
      data.gene.mat <- matrix(NA,num.label.Subject,len.label.roi)
      colnames(data.gene.mat) <- label.roi
      rownames(data.gene.mat) <- levels.subjectNum.L
      
      genename <- row.names(data.cor.mat.gene)[i]
      
      for(j in c(1:num.label.Subject)){
        sub.name <- levels.subjectNum.L[j]
        roi.name <- colnames(data.mean.mat.list[[j]])
        data.gene.mat[sub.name,roi.name] <- data.mean.mat.list[[j]][genename,]
      }
      data.gene.sub.cor <- data.my_cor(data.gene.mat)
      data.cor.mat.gene[genename,] <- data.gene.sub.cor
    }
    return(data.cor.mat.gene)  
  }
  return(geneCornetwork(sampledata,sampleinfo))
}
####build the individual network
geneNetbuilNet<-geneNetbuild(sampledata,sampleinfo)
#build the functional connectome,对每个人建模#

##分组构建脑网络 creat the data
functionaM<-data.rest1.bold[1:110,]
regioniD<-rownames(functionaM)
subject<-c()
for(i in c(1:10)){
 # i=1
  temp<-rep(paste('sub',i,sep='-'),11)
  subject<-c(subject,temp)
}
roidata<-c(rep(label.roi,10))
group<-c(rep("group1",55),rep('group2',55))
dataframInfo<-cbind(regioniD,roidata,subject,group)%>%as.matrix()
#functionaM<-data.sampel.bold[1:100,]
sampleRegioninfo<-dataframInfo
##
data.fmri.networkBuild <- function(functionaM,sampleRegioninfo,choseRegion=c(),mean=TRUE){
  # functionaM  -- the functional connectome of one people(include the l eft hemisphere and right hemisphere), row is regional sample. colnale is variables   
  # sampleRegion -- selecting regions, the rowdata
  # chosing region -- a vector for the brain region which you will use
  #choseRegion<-c('V1C','DFC')
  
  if(length(choseRegion)>0){
    tempdata<-matrix(0,1,ncol(sampleRegioninfo) )
    for (i in choseRegion) {
    #  i=choseRegion[1]
      tempdata<-rbind(tempdata,sampleRegioninfo[which(sampleRegioninfo[,2]==i),])
    }
    tempdata<-tempdata[-1,]
    sampleRegioninfo<-tempdata
  }
  roiNum<-sampleRegioninfo[,2]%>%unique()%>%length
  subjectNum<-sampleRegioninfo[,3]%>%unique()%>%length
  data.sample.bold.network <- matrix(NA,subjectNum,(roiNum*(roiNum-1))/2)
  ##对每个人构建脑功能连接网络
  for(i in c(1:subjectNum)){
    #=1
    ##select the responsed data
    subjectData<-sampleRegioninfo[,3]%>%unique()
    choseSubject<-sampleRegioninfo[which(sampleRegioninfo[,3]==subjectData[i]),]# select the subjct data
    choseData<-functionaM[choseSubject[,1],]
    roi_name<-choseSubject[,2]%>%unique()
    data.sample.bold.L.11 <- matrix(NA,roiNum,ncol(functionaM))
    ## if one region have some sample, using the mean value 
    for(j in c(1:roiNum)){
      #j=1
      data.roi.label <- choseSubject[which(choseSubject[,2]==roi_name[j]),]
      #data.roi.label <- choseSubject[which(choseSubject[,4]=="group1"),]
      if(!is.null(nrow(data.roi.label))){
        labeldata<-choseData[data.roi.label[,1],]
        label.index <- apply(labeldata,2,mean)
        data.sample.bold.L.11[j,]<-label.index
      }else{
        labeldata<-choseData[data.roi.label[1],]
        data.sample.bold.L.11[j,]<-as.matrix(labeldata)
      }
    }
    data.cor <- cor(t(data.sample.bold.L.11),method="pearson")
    #
    data.sample.bold.network[i,] <- data.cor[lower.tri(data.cor)]
  }
  if(mean){
    data.sample.bold.network<-apply(data.sample.bold.network,2,mean)
    data.sample.bold.networkEnd<-atanh(data.sample.bold.network)
  }
  data.sample.bold.networkEnd<-atanh(t(data.sample.bold.network))#fisher z transfrom
  return(data.sample.bold.network)
}
brainNet<-data.fmri.networkBuild(functionaM,sampleRegioninfo)
comCorrelatio<-function(brainNet,geneNetbuilNet,stablegene=c(),corMethod='pearson',adjmeh='fdr',cutTd=0.05){
  # brainNet -- the input of brain network[full connect]
  # geneNet -- the input of gene network.
  # stablegene --  gene vector which used to select the gene
  # corMethod --the method for computer the weight of edge
  # adjmeh-- the correct method for the p value of gene CCA
  # cutTd the threshold for select end gene
  if(length(stablegene)>0){
    networkPvalue<- apply(geneNetbuilNet[stablegene,],
                          1,function(y){cor.test(y,brainNet,
                          method = corMethod)$'p.value'})
  }else{
    networkPvalue<- apply(geneNetbuilNet,
                          1,function(y){cor.test(y,brainNet,
                          method = corMethod)$'p.value'})
  } 
  ##adjust the value
  if(length(adjmeh)>0){
    p_adjvalue<-p.adjust(networkPvalue, method =adjmeh)
  }
  ##end result
  endRes<-p_adjvalue[which(p_adjvalue<=cutTd)]
  return(endRes)
}

res<-comCorrelatio(brainNet,geneNetbuilNet,cutTd=0.5)






