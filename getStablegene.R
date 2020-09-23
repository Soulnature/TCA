

data.gene_expression <- read.table("F:/project/2.13/end_touse/data_gene_expression_NM.txt")
data.column <- read.table("F:/project/newkang_brainNet/kang/column.txt",sep="\t",header=T,row.names=1)
data.ncx.left.column<-data.column[which(data.column$Sample_label>=11&data.column$Sample_position=='cortical'&data.column$Sample_characteristics_hemisphere=='L'),]
colnames(data.gene_expression) <-rownames(data.ncx.left.column)
#data.NM<-read.table("F:/project/newkang_brainNet/kang/new_code_document/soft/NM_name.txt")
data.FG_NM<- read.table("F:/project/2.13/end_touse/FG_NM.txt",header=F)
data.NM<-data.FG_NM$V2
NM_info <- read.table('F:/project/newkang_brainNet/kang/new_code_document/Gene_Annotation/NM_info.txt')
NM_info_1<-NM_info[data.NM,]
factor.data.NM <- factor(NM_info_1$V2)
data.mean.gene_expression<-apply(data.gene_expression,2,function(x){tapply(x,factor.data.NM,mean)})

library(mice)
sampleinfo<-data.column[which(data.column$Sample_label=='11'|data.column$Sample_label=='12'),]
sampleinfo<-sampleinfo[which(sampleinfo$Sample_position=='cortical'&sampleinfo$Sample_characteristics_hemisphere=='L'),]
sampledata<-data.mean.gene_expression[1:300,rownames(sampleinfo)]
write.table(sampledata,'Stabaleinputsampledata.txt',quote = FALSE)
sampleinfo1<-sampleinfo[,2:3]
#1. find age stable gene
getStablegene<-function(expression,sampleinfo,cutoff=0.1,fixed=FALSE,Fixmethod="pmm"){
  expression<-sampledata
  sampleinfo<-sampleinfo1
  library(dplyr)
  library(mice)
  # expression --- expression data,row are gene ,and col are sample
  # sampleinfo --- sample info, include the mapping relations between the sample and subject
  # compute -- the correlation between the sample in different subjects
  # cutoff -- the  threshold  for selecting  the candidate gene, range from 0 to 1,default is 0.1
  # fixed --  whether data needs to be completed
  # Fixmethod -- fixed method based on mice('pmm', 'logreg', 'polyreg' ,'polr'),default is pmm 
  data.my_cor <- function(x,fixMethod,Fixmethod){
   # x<-data.gene.mat
    if(fixMethod){
      data<-mice(x,m=5,method = Fixmethod,maxit =20,seed=2)
      data <- complete(data,action=3)
      }
    data.cor.vecs <- c()
    len <- nrow(x)
    for(i in c(1:(len-1))){
      for(j in c((i+1):len)){
        data <- rbind(x[i,],x[j,])
        data.is.na <- apply(data,2,function(x){any(is.na(x))})
        if(any(data.is.na)){
          data.l <- data[,-which(data.is.na)]
        }else{
          data.l <- data
        }
        if(ncol(data.l)>2){
          data.cor.vecs <- c(data.cor.vecs,cor(data.l[1,],data.l[2,],method="pearson"))
        }
      }
    }
    data.cor.mean <- mean(data.cor.vecs)
    return(data.cor.mean)
  }
  #####
  data.cor.label.L <- function(expression,sampleinfo){
    
   
    colnames(sampleinfo)<-c('subject','sample')
    #get the subject number info
    subjectNum <- factor(sampleinfo$subject)
    levels.subjectNum.L <- levels(subjectNum)
    num.label.Subject <- length(levels.subjectNum.L)

    factor.data.column.label.region.L <- factor(sampleinfo$sample)
    levels.sample.label <- levels(factor.data.column.label.region.L)
    num.label.sample <- length(levels.sample.label)
    
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
    names(data.mean.mat.list)=paste("subject_",levels.subjectNum.L[1:num.label.Subject],sep='')
    
    label.roi <- levels.sample.label
    len.label.roi <- length(label.roi)
    data.cor.mat.gene <- matrix(NA,nrow(expression),1)
    rownames(data.cor.mat.gene) <- row.names(expression)
    
    for(i in c(1:nrow(data.cor.mat.gene))){
     # i=1
      data.gene.mat <- matrix(NA,num.label.Subject,len.label.roi)
      colnames(data.gene.mat) <- label.roi
      rownames(data.gene.mat) <- levels.subjectNum.L
      
      genename <- row.names(data.cor.mat.gene)[i]
      
      for(j in c(1:num.label.Subject)){
        sub.name <- levels.subjectNum.L[j]
        roi.name <- colnames(data.mean.mat.list[[j]])
        data.gene.mat[sub.name,roi.name] <- data.mean.mat.list[[j]][genename,]
      }
      data.gene.sub.cor <- data.my_cor(data.gene.mat,fixed,Fixmethod)
      data.cor.mat.gene[genename,] <- data.gene.sub.cor
    }
    return(data.cor.mat.gene)  
  }
  
 outputdata<-data.cor.label.L(c,sampleinfo) 
 ##sort the stable gene value 
 seqData<-outputdata[order(outputdata),]
 ##Screening the candidate genes
 cutoff=0.1
 cutoffvalue<-seqData[length(seqData)*cutoff]
 selectValue<-seqData[which(seqData>cutoffvalue)]%>%as.matrix()
 return(selectValue)
}


#data.mean.gene_expression <- apply(data.gene_expression,2,function(x){tapply(x,factor.data.NM,mean)})

data.column <- read.table("F:/project/newkang_brainNet/kang/column.txt",sep="\t",header=T,row.names=1)
data.ncx.left.column <- data.column[which(data.column$Sample_label>=11&data.column$Sample_position=='cortical'&data.column$Sample_characteristics_hemisphere=='L'),]
order.data.ncx.left.column <- data.ncx.left.column[order(data.ncx.left.column$Sample_label,data.ncx.left.column$Sample_characteristics_region,data.ncx.left.column$Sample_characteristics_hemisphere,data.ncx.left.column$Sample_characteristics_brain_code,decreasing=F),]

data.gene_expression.left.leave <- data.mean.gene_expression[,row.names(order.data.ncx.left.column)]
#data.gene_expression.left.leave.raw <- data.mean.gene_expression[,row.names(order.data.ncx.left.column)]
# data.gene_expression.left.leave.mean <- as.matrix(apply(data.gene_expression.left.leave.raw,1,mean))
# data.gene_expression.left.leave <- data.gene_expression.left.leave.raw[data.gene_expression.left.leave.mean>6,]

data.my_cor <- function(x){
  #x<-data.gene.mat
  # data<-mice(x,m=5,method = "pmm",maxit = 100,seed=1)
  data.cor.vecs <- c()
  len <- nrow(x)
  for(i in c(1:(len-1))){
    for(j in c((i+1):len)){
      data <- rbind(x[i,],x[j,])
      data.is.na <- apply(data,2,function(x){any(is.na(x))})
      if(any(data.is.na)){
        data.l <- data[,-which(data.is.na)]
      }else{
        data.l <- data
      }
      if(ncol(data.l)>2){
        data.cor.vecs <- c(data.cor.vecs,cor(data.l[1,],data.l[2,],method="pearson"))
      }
    }
  }
  data.cor.mean <- mean(data.cor.vecs)
  return(data.cor.mean)
}


data.cor.label.L <- function(index.label1,index.label2){
  #index.label1<-13
  #index.label2<-13
  data.column.label.L <- order.data.ncx.left.column[which(order.data.ncx.left.column$Sample_label>=index.label1&order.data.ncx.left.column$Sample_label<index.label2),]
  data.label.L <- data.gene_expression.left.leave[,row.names(data.column.label.L)]
  
  factor.data.column.label.HSB.L <- factor(data.column.label.L$Sample_characteristics_brain_code)
  levels.data.column.label.HSB.L <- levels(factor.data.column.label.HSB.L)
  num.label.HSB <- length(levels.data.column.label.HSB.L)
  
  factor.data.column.label.region.L <- factor(data.column.label.L$Sample_characteristics_region)
  levels.data.column.label.region.L <- levels(factor.data.column.label.region.L)
  num.label.region <- length(levels.data.column.label.region.L)
  
  data.mean.mat.list <- list()
  
  for(k in c(1:num.label.HSB)){
    data.mean.mat <- matrix(0,nrow(data.label.L),num.label.region)
    colnames(data.mean.mat) <- levels.data.column.label.region.L
    rownames(data.mean.mat) <- row.names(data.label.L)
    
    data.column <- data.column.label.L[which(data.column.label.L$Sample_characteristics_brain_code==levels.data.column.label.HSB.L[k]),]
    data <- as.matrix(data.label.L[,row.names(data.column)])
    
    colnames(data) <- data.column$Sample_characteristics_region
    
    data.mean.mat.list <- c(data.mean.mat.list, list(data))
  }
  names(data.mean.mat.list)=paste("subject_",levels.data.column.label.HSB.L[1:num.label.HSB],sep='')
  
  label.roi <- levels.data.column.label.region.L
  len.label.roi <- length(label.roi)
  data.cor.mat.gene <- matrix(NA,nrow(data.label.L),1)
  rownames(data.cor.mat.gene) <- row.names(data.label.L)
  
  for(i in c(1:nrow(data.cor.mat.gene))){
    data.gene.mat <- matrix(NA,num.label.HSB,len.label.roi)
    colnames(data.gene.mat) <- label.roi
    rownames(data.gene.mat) <- levels.data.column.label.HSB.L
    
    genename <- row.names(data.cor.mat.gene)[i]
    
    for(j in c(1:num.label.HSB)){
      sub.name <- levels.data.column.label.HSB.L[j]
      roi.name <- colnames(data.mean.mat.list[[j]])
      data.gene.mat[sub.name,roi.name] <- data.mean.mat.list[[j]][genename,]
    }
    
    data.gene.sub.cor <- data.my_cor(data.gene.mat)
    data.cor.mat.gene[genename,] <- data.gene.sub.cor
    
  }
  return(data.cor.mat.gene)  
}

#index.label1=11;index.label2=11
data.gene.cor.sample.5ageGroup.p1.L <- data.cor.label.L(11,13)
data.gene.cor.sample.5ageGroup.p2.L <- data.cor.label.L(13,14)
data.gene.cor.sample.5ageGroup.p3.L <- data.cor.label.L(14,16)
#data.gene.cor.sample.5ageGroup.p5.L <- data.cor.label.L(15,15)
write.table(data.gene.cor.sample.5ageGroup.p1.L,'F:/project/3.15/age_stableGene/age1_stable.txt')
write.table(data.gene.cor.sample.5ageGroup.p2.L,'F:/project/3.15/age_stableGene/age2_stable.txt')
write.table(data.gene.cor.sample.5ageGroup.p3.L,'F:/project/3.15/age_stableGene/age3_stable.txt')

#write.table(data.gene.cor.sample.5ageGroup.p5.L,'./lifespan/report20170120/ageDSgene/pearson/ageSampleCor/data.gene.cor.sample.5ageGroup.p5.L.txt')
data.gene.cor.sample.5ageGroup.p1.L<-as.matrix(read.table("F:/project/3.15/age_stableGene/age1_stable.txt"))
data.gene.cor.sample.5ageGroup.p2.L<-as.matrix(read.table("F:/project/3.15/age_stableGene/age1_stable.txt"))
data.gene.cor.sample.5ageGroup.p3.L<-as.matrix(read.table("F:/project/3.15/age_stableGene/age1_stable.txt"))
#data.gene.cor.sample.5ageGroup.p1.L<-data.gene.cor.sample.5ageGroup.p1.L$V1
#data.gene.cor.sample.5ageGroup.p2.L<-data.gene.cor.sample.5ageGroup.p2.L$V1
#data.gene.cor.sample.5ageGroup.p3.L<-data.gene.cor.sample.5ageGroup.p3.L$V1
#data.gene.cor.sample.5ageGroup.p4.L<-data.gene.cor.sample.5ageGroup.p4.L$V1
data.gene.cor.sample.5ageGroup.p1.L.sort <- sort(data.gene.cor.sample.5ageGroup.p1.L,decreasing=T)
data.gene.cor.sample.5ageGroup.p2.L.sort <- sort(data.gene.cor.sample.5ageGroup.p2.L,decreasing=T)
data.gene.cor.sample.5ageGroup.p3.L.sort <- sort(data.gene.cor.sample.5ageGroup.p3.L,decreasing=T)
data.gene.cor.sample.5ageGroup.p4.L.sort <- sort(data.gene.cor.sample.5ageGroup.p4.L,decreasing=T)
# data.gene.cor.sample.5ageGroup.p1.L.sort <- sort(abs(data.gene.cor.sample.5ageGroup.p1.L),decreasing=T)
# data.gene.cor.sample.5ageGroup.p2.L.sort <- sort(abs(data.gene.cor.sample.5ageGroup.p2.L),decreasing=T)
# data.gene.cor.sample.5ageGroup.p3.L.sort <- sort(abs(data.gene.cor.sample.5ageGroup.p3.L),decreasing=T)
# data.gene.cor.sample.5ageGroup.p4.L.sort <- sort(abs(data.gene.cor.sample.5ageGroup.p4.L),decreasing=T)
# data.gene.cor.sample.5ageGroup.p5.L.sort <- sort(abs(data.gene.cor.sample.5ageGroup.p5.L),decreasing=T)

#index.label1=11;index.label2=11;data.DS=data.gene.cor.sample.5ageGroup.p1.L;age.label="p1: 6-12Y"
data.pca <- function(index.label1,index.label2,data.DS,age.label){
  
  data.column.label.L <- order.data.ncx.left.column[which(order.data.ncx.left.column$Sample_label>=index.label1&order.data.ncx.left.column$Sample_label<=index.label2),]
  data.label.L <- data.gene_expression.left.leave[,row.names(data.column.label.L)]
  
  factor.data.column.label.HSB.L <- factor(data.column.label.L$Sample_characteristics_brain_code)
  levels.data.column.label.HSB.L <- levels(factor.data.column.label.HSB.L)
  num.label.HSB <- length(levels.data.column.label.HSB.L)
  
  factor.data.column.label.region.L <- factor(data.column.label.L$Sample_characteristics_region)
  levels.data.column.label.region.L <- levels(factor.data.column.label.region.L)
  num.label.region <- length(levels.data.column.label.region.L)
  
  label.roi <- levels.data.column.label.region.L
  
  data.l.mat <- matrix()
  for(k in c(1:num.label.HSB)){
    data.mean.mat <- matrix(0,nrow(data.label.L),num.label.region)
    colnames(data.mean.mat) <- levels.data.column.label.region.L
    rownames(data.mean.mat) <- row.names(data.label.L)
    
    data.column <- data.column.label.L[which(data.column.label.L$Sample_characteristics_brain_code==levels.data.column.label.HSB.L[k]),]
    data <- as.matrix(data.label.L[,row.names(data.column)])
    
    colnames(data) <- data.column$Sample_characteristics_region
    
    data.l <- data[rownames(as.matrix(data.DS)),]
    
    if(k==1){
      data.l.mat <- t(data.l)
    }else{
      data.l.mat <- rbind(data.l.mat,t(data.l))
    }
  }
  
  data.label.l <- factor(rownames(data.l.mat))
  data.label <- apply(as.matrix(data.label.l),1,function(x){return(which(x==label.roi))})
  
  library("vegan")
  #cols <- rainbow(11)
  #cols <- c("white","black","red","orange","yellow","green","blue","gray","purple","cyan","pink")
  cols <- c("white","black","red","orange","yellow","green","blue","gray","purple","cyan","pink")
  data.col <- apply(as.matrix(data.label),1,function(x){cols[x]})
  
  sol <-metaMDS(data.l.mat,distance="bray")
  plot(sol, type = "p",display = c("sites"),main=age.label)
  points(sol, "sites", pch=21,  bg=data.col)
  #text(sol)
  #legend("topright",label.roi,col=cols,pch=21,box.col="black",cex=0.4,text.width=0.4)
}

par(mfrow=c(2,3));
data.pca(11,11,data.gene.cor.sample.5ageGroup.p1.L,"p1: 6-12Y")
data.pca(12,12,data.gene.cor.sample.5ageGroup.p2.L,"p2: 12-20Y")
data.pca(13,13,data.gene.cor.sample.5ageGroup.p3.L,"p3: 20-40Y")
data.pca(14,14,data.gene.cor.sample.5ageGroup.p4.L,"p4: 40-60Y")
data.pca(15,15,data.gene.cor.sample.5ageGroup.p5.L,"p5: 60-Y")

percentage <- c(0.1,0.3,0.5,0.75)
#percentage <- c(0.9,0.1,0.15,0.25,0.5,0.75)
#percentage <- c(0.9,0.75,0.5,0.25)
for(k in c(1:4)){
  #window();
  k=1
  threshold.p1 <- data.gene.cor.sample.5ageGroup.p1.L.sort[length(which(data.gene.cor.sample.5ageGroup.p1.L.sort>0))*percentage[k]]
  threshold.p2 <- data.gene.cor.sample.5ageGroup.p2.L.sort[length(which(data.gene.cor.sample.5ageGroup.p2.L.sort>0))*percentage[k]]
  threshold.p3 <- data.gene.cor.sample.5ageGroup.p3.L.sort[length(which(data.gene.cor.sample.5ageGroup.p3.L.sort>0))*percentage[k]]
  #threshold.p4 <- data.gene.cor.sample.5ageGroup.p3.L.sort[length(which(data.gene.cor.sample.5ageGroup.p4.L.sort>0))*percentage[k]]
  #threshold.p4 <- data.gene.cor.sample.5ageGroup.p4.L.sort[length(which(data.gene.cor.sample.5ageGroup.p4.L.sort>0))*percentage[k]]
  #threshold.p5 <- data.gene.cor.sample.5ageGroup.p5.L.sort[length(which(data.gene.cor.sample.5ageGroup.p5.L.sort>0))*percentage[k]]
  
  data.DS.5ageGroup.p1.L <- data.gene.cor.sample.5ageGroup.p1.L[which(data.gene.cor.sample.5ageGroup.p1.L>=threshold.p1),]
  data.DS.5ageGroup.p2.L <- data.gene.cor.sample.5ageGroup.p2.L[which(data.gene.cor.sample.5ageGroup.p2.L>=threshold.p2),]
  data.DS.5ageGroup.p3.L <- data.gene.cor.sample.5ageGroup.p3.L[which(data.gene.cor.sample.5ageGroup.p3.L>=threshold.p3),]
  #data.DS.5ageGroup.p4.L <- data.gene.cor.sample.5ageGroup.p3.L[which(data.gene.cor.sample.5ageGroup.p4.L>=threshold.p4),]
  #data.DS.5ageGroup.p4.L <- data.gene.cor.sample.5ageGroup.p4.L[which(data.gene.cor.sample.5ageGroup.p4.L>=threshold.p4),]
  #data.DS.5ageGroup.p5.L <- data.gene.cor.sample.5ageGroup.p5.L[which(data.gene.cor.sample.5ageGroup.p5.L>=threshold.p5),]
  write.table(data.DS.5ageGroup.p1.L,paste('F:/project/3.15/age_stableGene/',percentage[k],'/data.DS.5ageGroup.p11.L.txt',sep=""))
  write.table(data.DS.5ageGroup.p2.L,paste('F:/project/3.15/age_stableGene/',percentage[k],'/data.DS.5ageGroup.p22.L.txt',sep=""))
  write.table(data.DS.5ageGroup.p3.L,paste('F:/project/3.15/age_stableGene/',percentage[k],'/data.DS.5ageGroup.p33.L.txt',sep=""))
  #write.table(data.DS.5ageGroup.p4.L,paste('F:/project/2.13/ageStable_gene/four_LR/',percentage[k],'/data.DS.5ageGroup.p44.L.txt',sep=""))
  ## write.table(data.DS.5ageGroup.p5.L,paste('F:/lifespan/report20170120/ageDSgene/pearson/ageDSgeneCor',percentage[k],'/data.DS.5ageGroup.p5.L.txt',sep=""))
  
  # data.DS.5ageGroup.p1.L <- as.matrix(read.table(paste('F:/lifespan/report20170120/ageDSgene/pearson/ageDSgeneCor',percentage[k],'/data.DS.5ageGroup.p1.L.txt',sep="")))
  # data.DS.5ageGroup.p2.L <- as.matrix(read.table(paste('F:/lifespan/report20170120/ageDSgene/pearson/ageDSgeneCor',percentage[k],'/data.DS.5ageGroup.p2.L.txt',sep="")))
  # data.DS.5ageGroup.p3.L <- as.matrix(read.table(paste('F:/lifespan/report20170120/ageDSgene/pearson/ageDSgeneCor',percentage[k],'/data.DS.5ageGroup.p3.L.txt',sep="")))
  # data.DS.5ageGroup.p4.L <- as.matrix(read.table(paste('F:/lifespan/report20170120/ageDSgene/pearson/ageDSgeneCor',percentage[k],'/data.DS.5ageGroup.p4.L.txt',sep="")))
  # data.DS.5ageGroup.p5.L <- as.matrix(read.table(paste('F:/lifespan/report20170120/ageDSgene/pearson/ageDSgeneCor',percentage[k],'/data.DS.5ageGroup.p5.L.txt',sep="")))
  # 
  # par(mfrow=c(2,3));
  # data.pca(11,11,data.DS.5ageGroup.p1.L,"p1: 6-12Y")
  # data.pca(12,12,data.DS.5ageGroup.p2.L,"p2: 12-20Y")
  # data.pca(13,13,data.DS.5ageGroup.p3.L,"p3: 20-40Y")
  # data.pca(14,14,data.DS.5ageGroup.p4.L,"p4: 40-60Y")
  # data.pca(15,15,data.DS.5ageGroup.p5.L,"p5: 60-Y")
}


