fullAge1<-read.table("Full_age1.txt")
fullAge2<-read.table("Full_age2.txt")
fullAge3<-read.table("Full_age3.txt")
#unNormalize build the gene network
data.my_cor <- function(x){
  length <- ncol(x)
  data.cor <- as.vector(matrix(NA,55,1))
  k=1
  for(i in c(1:(length-1))){
    for(j in c((i+1):length)){
      data <- x[,c(i,j)]
      data.cor[k] <- cor(data[,1],data[,2])
      k=k+1
    }
  }
  #data.cor<-atanh(data.cor)
  return(data.cor)
}
geneNet_Build<-function(data){
  #data<-fullAge1
  len<-nrow(data)
  num<-len/11
  
  gene_net<-apply(data, 2, function(x){
   # x<-data[,1]
    cor_matrix<-matrix(NA,num,11)
    k=1
    for(i in c(1:num)){
      cor_matrix[i,]<-x[k:(k+10)]
      k=k+11
    }
   data_cor<-data.my_cor(cor_matrix)
   }
  )
data_cor1<-t(gene_net)
}
age1_net<-geneNet_Build(fullAge1)
age2_net<-geneNet_Build(fullAge2)
age3_net<-geneNet_Build(fullAge3)
#把填充的脑网络写入到文件中去
write.table(age1_net,"/home1/zhaoxz/genenet/full_agenet1.txt")
write.table(age2_net,"/home1/zhaoxz/genenet/full_agenet2.txt")
write.table(age3_net,"/home1/zhaoxz/genenet/full_agenet3.txt")

####

data.gene.cor.label.5ageGroup.p1.L<-read.table("/home1/zhaoxz/genenet/full_agenet1.txt")
data.gene.cor.label.5ageGroup.p2.L<-read.table("/home1/zhaoxz/genenet/full_agenet2.txt")
data.gene.cor.label.5ageGroup.p3.L<-read.table("/home1/zhaoxz/genenet/full_agenet3.txt")
########build the image net

data.fmri.sample.cor.network <- function(file_path,column_index){
  
  
  data.rest1.bold.LR <- read.table(paste("./report20160227/",file_path,"/fmri/REST1/COLUMN_",column_index,"_rfMRI_REST1_LR_Atlas_hp2000_clean.dtseries.txt",sep=""),header=F)
  data.rest1.bold.RL <- read.table(paste("./report20160227/",file_path,"/fmri/REST1/COLUMN_",column_index,"_rfMRI_REST1_RL_Atlas_hp2000_clean.dtseries.txt",sep=""),header=F)
  #data.rest1.bold <- cbind(data.rest1.bold.RL)
  data.rest1.bold <- cbind(data.rest1.bold.LR)
  #data.rest1.bold <- cbind(data.rest1.bold.RL,data.rest1.bold.LR)
  data.column.L <- read.table(paste("./report20160227/",file_path,"/fmri/COLUMN_",column_index,"_L.txt",sep=""))
  data.column.R <- read.table(paste("./report20160227/",file_path,"/fmri/COLUMN_",column_index,"_R.txt",sep=""))
  each <- nrow(data.column.L) + nrow(data.column.R) + 19
  
  data.label.L <- read.table(paste("./report20160227/",file_path,"_L.txt",sep=""))
  data.label_to_11.L <- read.table(paste("./report20160227/",file_path,"_L_1.txt",sep=""))
  data.sample.bold.network <- matrix(NA,27,55)
  for(i in c(1:27)){
    #i=1
    data.sampel.bold <- data.rest1.bold[(each*(i-1)+1):(each*i),]
    data.sample.bold.L.11 <- matrix(NA,11,ncol(data.sampel.bold))
    for(j in c(1:11)){
      data.roi.label <- as.matrix(data.label_to_11.L[which(data.label_to_11.L[,2]==j),1])
      label.index <- apply(data.roi.label,1,function(x){
        return(which(data.label.L[,3]==x))
      })
      data.sample.bold.L.11[j,] <- apply(t(data.sampel.bold[label.index,]),1,mean)
    }
    
    #data.sample.bold.L.11<-apply(data.sample.bold.L.11,1,scale)
    data.cor <- cor(t(data.sample.bold.L.11),method="pearson")
    data.sample.bold.network[i,] <- data.cor[lower.tri(data.cor)]
  }
  
  data.sample.bold.network<-atanh(data.sample.bold.network)
  return(data.sample.bold.network)
}



#1. brain networks
## For parcellations_Vgd11b COLUMN_2 Label
data.sample.ba <- data.fmri.sample.cor.network("parcellations_VGD11b",2)
#data.sample.ba <- data.fmri.sample.cor.network("RSN-networks",2)
#For five age group:|11|12|13|14|15|
data.label.fmri.cor.5ageGroup.p1 <- apply(data.sample.ba[1:12,],2,mean)
#data.label.fmri.cor.5ageGroup.p2 <- apply(data.sample.ba[7:12,],2,mean)
data.label.fmri.cor.5ageGroup.p2 <- apply(data.sample.ba[13:17,],2,mean)
data.label.fmri.cor.5ageGroup.p3 <- apply(data.sample.ba[18:27,],2,mean)
#data.label.fmri.cor.5ageGroup.p5 <- apply(data.sample.ba[23:27,],2,mean)

age.DS.percentage <- 0.1
data.DS.5ageGroup.p1.L <- as.matrix(read.table(paste('./',age.DS.percentage,'/data.DS.5ageGroup.p11.L.txt',sep="")))
data.DS.5ageGroup.p2.L <- as.matrix(read.table(paste('./',age.DS.percentage,'/data.DS.5ageGroup.p22.L.txt',sep="")))
data.DS.5ageGroup.p3.L <- as.matrix(read.table(paste('./',age.DS.percentage,'/data.DS.5ageGroup.p33.L.txt',sep="")))
data.gene.cor.label.5ageGroup.p1.L<-atanh(data.gene.cor.label.5ageGroup.p1.L)
data.gene.cor.label.5ageGroup.p2.L<-atanh(data.gene.cor.label.5ageGroup.p2.L)
data.gene.cor.label.5ageGroup.p3.L<-atanh(data.gene.cor.label.5ageGroup.p3.L)

data.label.brain.gene.cor.5ageGroup.p1.pvalue <- apply(data.gene.cor.label.5ageGroup.p1.L[rownames(data.DS.5ageGroup.p1.L),],1,function(y){cor.test(y,data.label.fmri.cor.5ageGroup.p1,method = "pearson")$'p.value'})
data.label.brain.gene.cor.5ageGroup.p2.pvalue <- apply(data.gene.cor.label.5ageGroup.p2.L[rownames(data.DS.5ageGroup.p2.L),],1,function(y){cor.test(y,data.label.fmri.cor.5ageGroup.p2,method = "pearson")$'p.value'}) 
data.label.brain.gene.cor.5ageGroup.p3.pvalue <- apply(data.gene.cor.label.5ageGroup.p3.L[rownames(data.DS.5ageGroup.p3.L),],1,function(y){cor.test(y,data.label.fmri.cor.5ageGroup.p3,method = "pearson")$'p.value'}) 





# data.label.brain.gene.cor.5ageGroup.p1.pvalue <- apply(data.gene.cor.label.5ageGroup.p1.L,1,function(y){cor.test(y,data.label.fmri.cor.5ageGroup.p1,method = "pearson")$'p.value'})
# data.label.brain.gene.cor.5ageGroup.p2.pvalue <- apply(data.gene.cor.label.5ageGroup.p2.L,1,function(y){cor.test(y,data.label.fmri.cor.5ageGroup.p2,method = "pearson")$'p.value'})
# data.label.brain.gene.cor.5ageGroup.p3.pvalue <- apply(data.gene.cor.label.5ageGroup.p3.L,1,function(y){cor.test(y,data.label.fmri.cor.5ageGroup.p3,method = "pearson")$'p.value'})
NM_info <- read.table('./NM_info.txt',header =T,row.names=1,sep='\t')
data.cor.gene.5ageGroup.p1 <- data.label.brain.gene.cor.5ageGroup.p1.pvalue[which(data.label.brain.gene.cor.5ageGroup.p1.pvalue<=0.05)]
data.NM_info.label.cor.5ageGroup.p1 <- as.vector(na.omit(unique(NM_info[rownames(as.matrix(data.cor.gene.5ageGroup.p1)),])))
data.cor.gene.5ageGroup.p2 <- data.label.brain.gene.cor.5ageGroup.p2.pvalue[which(data.label.brain.gene.cor.5ageGroup.p2.pvalue<=0.05)]
data.NM_info.label.cor.5ageGroup.p2 <- as.vector(na.omit(unique(NM_info[rownames(as.matrix(data.cor.gene.5ageGroup.p2)),])))
data.cor.gene.5ageGroup.p3 <- data.label.brain.gene.cor.5ageGroup.p3.pvalue[which(data.label.brain.gene.cor.5ageGroup.p3.pvalue<=0.05)]
data.NM_info.label.cor.5ageGroup.p3 <- as.vector(na.omit(unique(NM_info[rownames(as.matrix(data.cor.gene.5ageGroup.p3)),])))
length(data.NM_info.label.cor.5ageGroup.p1)
length(data.NM_info.label.cor.5ageGroup.p2)
length(data.NM_info.label.cor.5ageGroup.p3)
####all fisher content
write.table(data.NM_info.label.cor.5ageGroup.p1,'./endAllFisherage1.txt',row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(data.NM_info.label.cor.5ageGroup.p2,'./endAllFisherage2.txt',row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(data.NM_info.label.cor.5ageGroup.p3,'./endAllFisherage3.txt',row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(as.matrix(data.cor.gene.5ageGroup.p1),'./endAllFisherTransage1.txt')
write.table(as.matrix(data.cor.gene.5ageGroup.p2),'./endAllFisherTransage2.txt')
write.table(as.matrix(data.cor.gene.5ageGroup.p3),'./endAllFisherTransage3.txt')
intersect(age2Cor_geneSymbol,data.NM_info.label.cor.5ageGroup.p2)
intersect(age3Cor_geneSymbol,data.NM_info.label.cor.5ageGroup.p3)
age1Cor_geneSymbol<-data.NM_info.label.cor.5ageGroup.p1
age2Cor_geneSymbol<-data.NM_info.label.cor.5ageGroup.p2
age3Cor_geneSymbol<-data.NM_info.label.cor.5ageGroup.p3
####
write.table(data.NM_info.label.cor.5ageGroup.p1,"./end_touse/ztrans/ageCoregene11.txt")
write.table(data.NM_info.label.cor.5ageGroup.p2,"./end_touse/ztrans/ageCoregene21.txt")
write.table(data.NM_info.label.cor.5ageGroup.p3,"./end_touse/ztrans/ageCoregene31.txt")
write.table(data.cor.gene.5ageGroup.p1,"./end_touse/ztrans/ageCoretrans11.txt")
write.table(data.cor.gene.5ageGroup.p2,"./end_touse/ztrans/ageCoretrans21.txt")
write.table(data.cor.gene.5ageGroup.p3,"./end_touse/ztrans/ageCoretrans31.txt")



