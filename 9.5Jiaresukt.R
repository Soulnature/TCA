library(dplyr)
library(export)
library(ggplot2)
##read the find gene
age1Cor_geneSymbol1<-read.table("F://project/11.15/endFisher/endAllFisherage1.txt",stringsAsFactors = F)
age2Cor_geneSymbol1<-read.table("F://project/11.15/endFisher/endAllFisherage2.txt",stringsAsFactors = F)
age3Cor_geneSymbol1<-read.table("F://project/11.15/endFisher/endAllFisherage3.txt",stringsAsFactors = F)
age1Cor_geneSymbol<-age1Cor_geneSymbol1$V1
age2Cor_geneSymbol<-age2Cor_geneSymbol1$V1
age3Cor_geneSymbol<-age3Cor_geneSymbol1$V1
allgene<-union(age1Cor_geneSymbol,union(age2Cor_geneSymbol,age3Cor_geneSymbol))%>%unique()#369个gene
###针对这315个gene进行差异表达基因的检验
data.gene_expression <- as.matrix(read.table("F://project/2.13/end_touse/data_gene_expression_NM.txt"))
colnames <- read.table("F://project/newkang_brainNet/kang/new_code_document/soft/column_names.txt",header=F)
colnames(data.gene_expression) <- colnames$V1
data.FG_NM<- read.table("F://project/2.13/end_touse/FG_NM.txt",header=F)
data.NM<-data.FG_NM$V2
factor.data.NM <- factor(data.NM)
data.mean.gene_expression <- apply(data.gene_expression,2,function(x){tapply(x,factor.data.NM,mean)})
data.column <- read.table("F://project/newkang_brainNet/kang/column.txt",sep="\t",header=T,row.names=1)
data.ncx.left.column <- data.column[which(data.column$Sample_label>=11&data.column$Sample_position=='cortical'&data.column$Sample_characteristics_hemisphere=='L'),]
order.data.ncx.left.column <- data.ncx.left.column[order(data.ncx.left.column$Sample_label,data.ncx.left.column$Sample_characteristics_region,data.ncx.left.column$Sample_characteristics_hemisphere,data.ncx.left.column$Sample_characteristics_brain_code,decreasing=F),]
data.gene_expression.left.leave <- data.mean.gene_expression[,row.names(order.data.ncx.left.column)]

#index.label1=11;index.label2=11
data.cor.label.L <- function(index.label1,index.label2){
  #index.label1=11
  #index.label2=12
  data.column.label.L <- order.data.ncx.left.column[which(order.data.ncx.left.column$Sample_label>=index.label1&order.data.ncx.left.column$Sample_label<=index.label2),]
  data.label.L <- data.gene_expression.left.leave[,row.names(data.column.label.L)]
  
  factor.data.column.label.HSB.L <- factor(data.column.label.L$Sample_characteristics_brain_code)
  levels.data.column.label.HSB.L <- levels(factor.data.column.label.HSB.L)
  num.label.HSB <- length(levels.data.column.label.HSB.L)
  
  factor.data.column.label.region.L <- factor(data.column.label.L$Sample_characteristics_region)
  levels.data.column.label.region.L <- levels(factor.data.column.label.region.L)
  num.label.region <- length(levels.data.column.label.region.L)
  
  data.mean.mat.list <- list()
  ####
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
  
  data.giExp.raw = do.call(cbind,data.mean.mat.list)
  
  roi.label <- c("A1C","DFC","IPC","ITC","M1C","MFC","OFC","S1C","STC","V1C","VFC")
  data.reg_label = sapply(colnames(data.giExp.raw),function(x){return(which(roi.label==x))})
  colnames(data.giExp.raw) <- c(1:ncol(data.giExp.raw))
  
  data.giExp <- data.frame(label_age=index.label1,label_region=data.reg_label,giExp=t(data.giExp.raw))#这个地方将第一个年龄段设置成11，第二个设置成12，第三个设置成13
  colnames(data.giExp) <- c("ages","regions",rownames(data.giExp.raw))
  
  return(data.giExp)
}

#-----single factor anova-----Age-----------------------------
#index.label1=11;index.label2=11
data3.gene.cor.sample.5ageGroup.p1.L <- data.cor.label.L(11,12)
data3.gene.cor.sample.5ageGroup.p2.L <- data.cor.label.L(13,13)
data3.gene.cor.sample.5ageGroup.p3.L <- data.cor.label.L(14,15)
#data.gene.cor.sample.5ageGroup.p4.L <- data.cor.label.L(14,14)
#data.gene.cor.sample.5ageGroup.p5.L <- data.cor.label.L(15,15)
#x<-data.gene_aov_mean[1,]

#data.gene.cor.sample.5ageGroup.p1.L <- data3.gene.cor.sample.5ageGroup.p1.L[,rownames(data3.DS.5ageGroup.p1.L)
data.gene.cor.sample.5ageGroup.L <- rbind(data3.gene.cor.sample.5ageGroup.p1.L,data3.gene.cor.sample.5ageGroup.p2.L,data3.gene.cor.sample.5ageGroup.p3.L,deparse.level = 1)
#将所有年龄段取均值，得出3个年龄段的均值，为后面得出了anova分析的结果后，证明基因是该年龄段的差异表达基因（该年龄段的基因表达值>其他四个年龄段的表达的均值，可认为这个基因是差异
#表达的）
anovaResult<-function(age){
# age<-11
  data.age.label <- data.gene.cor.sample.5ageGroup.L[,1]
  data.giExp.len <- ncol(data.gene.cor.sample.5ageGroup.L)
  #####
  data.gene.roi.mean.5ageGroup_NM.p1<-apply(data.gene.cor.sample.5ageGroup.L[,3:data.giExp.len],1,function(x){tapply(x,factor(NM_info.all),mean)}) 
  NM_info.all<-NM_info[rownames(data.gene_expression.left.leave),]
  findgene<-data.gene.roi.mean.5ageGroup_NM.p1[allgene,]
  print(findgene)
  if(age==11){
    data.age.label[which(data.age.label==13)]<-14
  }else if(age==13){
    data.age.label[which(data.age.label==11)]<-14
  }else if(age==14){
    data.age.label[which(data.age.label==13)]<-11
  }
  #sample()
  anovaData<-rbind(data.age.label,findgene)
  #DS expression
  data.age.exp.single.aov.pval <- as.matrix(apply(anovaData[2:nrow(anovaData),],1,function(data.age.exp){
    # data.age.exp<-anovaData[2,]
    data.sgi.exp <- data.frame(data.age.exp,data.age.label=factor(data.age.label))
    fit <- aov(data.age.exp~data.age.label,data=data.sgi.exp)
    result_1<- summary(fit)
    result_1.pr <- sapply(result_1, function(v) return(v[5]))$`Pr(>F)`
    result_1.pval <- result_1.pr[1]
    # result<-p.adjust(result_1.pval,method = "bonferroni",length(result_1.pval))
    # if(result<0.05)
    #{
    #  return(result)
    # }
    # else
    #  return(0)
    return(result_1.pval)
    ###################
  }))
  #data_age_exp_adjust_pval<-p.adjust(data.age.exp.single.aov.pval,method ="fdr",n=length(data.age.exp.single.aov.pval))%>%as.matrix()
  #rownames(data_age_exp_adjust_pval)<-rownames(data.age.exp.single.aov.pval)
  ###取出anova分析中P-value小于0.01的作为差异表达基因
  #data.age.single.aov.gi <- data.age.exp.single.aov.pval[which(data_age_exp_adjust_pval<0.01),]
 # data_age_aov_Diff_gene<-as.matrix(data.age.single.aov.gi)##228 gene
  return(data.age.exp.single.aov.pval)
}
AovResultage1<-anovaResult(11)
AovResultage2<-anovaResult(13)
AovResultage3<-anovaResult(14)
#permutation result### no gene can remain after permutaiton filite
PermutationanovaResult<-function(age){
  #age<-11
  data.age.label <- data.gene.cor.sample.5ageGroup.L[,1]
  data.giExp.len <- ncol(data.gene.cor.sample.5ageGroup.L)
  #####
  data.gene.roi.mean.5ageGroup_NM.p1<-apply(data.gene.cor.sample.5ageGroup.L[,3:data.giExp.len],1,function(x){tapply(x,factor(NM_info.all),mean)}) 
  NM_info.all<-NM_info[rownames(data.gene_expression.left.leave),]
  findgene<-data.gene.roi.mean.5ageGroup_NM.p1[allgene,]
  print(findgene)
  if(age==11){
    data.age.label[which(data.age.label==13)]<-14
  }else if(age==13){
    data.age.label[which(data.age.label==11)]<-14
  }else if(age==14){
    data.age.label[which(data.age.label==13)]<-11
  }
  savematrix<-matrix(NA,1000,length(allgene)) 
  for (i in c(1:1000)) {
    example<-sample(232,232)
    findgene<-findgene[,example]
    anovaData<-rbind(data.age.label,findgene)
    #DS expression
    data.age.exp.single.aov.pval <- as.matrix(apply(anovaData[2:nrow(anovaData),],1,function(data.age.exp){
      # data.age.exp<-anovaData[2,]
      data.sgi.exp <- data.frame(data.age.exp,data.age.label=factor(data.age.label))
      fit <- aov(data.age.exp~data.age.label,data=data.sgi.exp)
      result_1<- summary(fit)
      result_1.pr <- sapply(result_1, function(v) return(v[5]))$`Pr(>F)`
      result_1.pval <- result_1.pr[1]
      # result<-p.adjust(result_1.pval,method = "bonferroni",length(result_1.pval))
      # if(result<0.05)
      #{
      #  return(result)
      # }
      # else
      #  return(0)
      return(result_1.pval)
      ###################
    }))
   # i=1
    savematrix[i,]<-data.age.exp.single.aov.pval
  }
  return(savematrix)
}
 permutation1<-PermutationanovaResult(11)
 permutation2<-PermutationanovaResult(13)
 permutation3<-PermutationanovaResult(14)


######permutation 过滤一下gene
test_permutation<-function(real,permutationRes){
 # real<-AovResultage2
 # permutationRes<-permutation2
  colnames(permutationRes)<-allgene
  TrueGene<-vector()
  saveMatrix<-matrix(NA,length(real),1)
  rownames(saveMatrix)<-allgene
  #permutationRes<-permutationRes[,rownames(real)]
  for (i in c(1:length(real))) {
    #i<-2
    i
    tempdta<-permutationRes[,i]
    sob<-real[i]
    tempdtanum<-tempdta[which(tempdta<sob)]%>%length()
    tempValue<-tempdtanum/1000
    saveMatrix[i]<-tempValue
  }
  return(saveMatrix)
}
####全过了permutation了
permutationResage1<-test_permutation(AovResultage1,permutation1)
permutationResage2<-test_permutation(AovResultage2,permutation2)
permutationResage3<-test_permutation(AovResultage3,permutation3)

#####
AovResultage1<-permutationResage1[which(permutationResage1<0.01),]##201
AovResultage2<-permutationResage2[which(permutationResage2<0.01),]##99
AovResultage3<-permutationResage3[which(permutationResage3<0.01),]##180
####将最后的permutation的结果放在project/paperResult/permutationAOV下面
dir.create("./permutationAov")
write.table(AovResultage1,"AovResultage1.txt")
write.table(AovResultage2,"AovResultage2.txt")
write.table(AovResultage3,"AovResultage3.txt")
#########oadsadzrignal DS gene
# AovResultage1<-age1anovaSymbol
# AovResultage2<-age2anovaSymbol
# AovResultage3<-age3anovaSymbol


#####
AovResultage1<-read.table("AovResultage1.txt")
AovResultage2<-read.table("AovResultage2.txt")
AovResultage3<-read.table("AovResultage3.txt")
# ###看下有无重复的，去掉重复的
# intersect(inter1,inter2)
# intersect(inter1,00inter3)
####
write.table(rownames(AovResultage1),'AoVAge1.txt',row.names = FALSE, col.names = FALSE,quote = FALSE)





####

# inter1<-intersect(AovResultage1,age1Cor_geneSymbol)#57
# inter2<-intersect(AovResultage1,age2Cor_geneSymbol)#46
# inter3<-intersect(AovResultage3,age3Cor_geneSymbol)#73
length(inter1)/length(age1Cor_geneSymbol)
length(inter2)/length(age1Cor_geneSymbol)
length(inter3)/length(age1Cor_geneSymbol)
##
inter1<-intersect(rownames(AovResultage1),age1Cor_geneSymbol)#80
inter2<-intersect(rownames(AovResultage2),age2Cor_geneSymbol)#35
inter3<-intersect(rownames(AovResultage3),age3Cor_geneSymbol)#74
#####the gene which related to connectome but not in DS gene set
connectomeAge1<-setdiff(age1Cor_geneSymbol,inter1)
connectomeAge2<-setdiff(age2Cor_geneSymbol,inter2)
connectomeAge3<-setdiff(age3Cor_geneSymbol,inter3)
####将这两组基因写入到文件中去
dir.create('Twogroups')
write.table(inter1,'DEgenesAgeAndCon1.txt',row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(inter2,'DEgenesAgeAndCon2.txt',row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(inter3,'DEgenesAgeAndCon3.txt',row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(connectomeAge1,'OthersAge1.txt',row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(connectomeAge2,'OthersAge2.txt',row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(connectomeAge3,'OthersAge3.txt',row.names = FALSE,col.names = FALSE,quote = FALSE)
######



### data visualization of the two gene group
ageGroup<-c(rep("8-20",2),rep("20-40",2),rep("40-",2))
levels(ageGroup)<-c("8-20","20-40","40-")
condition<-rep(c("DEgene","connectome"),3)
value<-c(length(inter1),length(connectomeAge1),length(inter2),length(connectomeAge2),length(inter3),length(connectomeAge3))

plotData<-data.frame(ageGroup,condition,value)
#levels(plotData$ageGroup)<-factor(c("8-20","20-40","40-"))
ggplot(plotData,aes(x=ageGroup,y=value,fill=condition))+
  geom_bar(stat="identity")+
  geom_text(aes(label=value),position=position_stack(vjust=0.5))#3 vjust参数用来调整标签的为重，vjust=0.5将标签放到对应部位的中部
##
graph2ppt(file="twoGroup.pptx", width=10, height=7)
###test draw the histogram
# library(ggplot2)
# with(mpg,table(class,year))
# p <- ggplot(data=mpg,aes(x=class,fill=factor(year)))
# p + geom_bar(position='dodge')
# p + geom_bar(position='stack')
# p + geom_bar(position='fill')
# p + geom_bar(position='identity',alpha=0.3)
###
#######correlation with structureal feature
gene_cor<-function(name,age,x)
{  
  #x=1
 #age=inter1
 #name<-"thickness"
  if(x==1){
    expressionData<-data.gene.roi.mean.5ageGroup_NM.p1[age,]
  }else if(x==2){
    expressionData<-data.gene.roi.mean.5ageGroup_NM.p2[age,]
  }else if(x==3){
    expressionData<-data.gene.roi.mean.5ageGroup_NM.p3[age,]
  }
  #$#age<-data.gene.5ageGroup.p1
  #x=1
  age_area1<-as.matrix(read.table(paste('F://project/Smri_feature/',name,'.txt',sep = ""),header = T))
  age1_roi<-age_area1[x,]
  NM_gene_data<-t(expressionData)
  pca_gene<-prcomp(NM_gene_data)
  a<-cor.test(pca_gene$x[,1],age1_roi,method="pearson")
  #return(scc_value)
  print(a)
  }


######list r-z transform
gene_corZvalue<-function(name,age,x)
{  
 #  x=1
  # age=inter1
    name<-"thickness"
  if(x==1){
    expressionData<-data.gene.roi.mean.5ageGroup_NM.p1[age,]
  }else if(x==2){
    expressionData<-data.gene.roi.mean.5ageGroup_NM.p2[age,]
  }else if(x==3){
    expressionData<-data.gene.roi.mean.5ageGroup_NM.p3[age,]
  }
  #$#age<-data.gene.5ageGroup.p1
  #x=1
  age_area1<-as.matrix(read.table(paste('F://project/Smri_feature/',name,'.txt',sep = ""),header = T))
  age1_roi<-age_area1[x,]
  rvalue<-apply(expressionData, 1, function(x){
    cor(x,age1_roi)
  })
  ###tranform the r to z
  #rvalueZ<-0.5 * log( (1+rvalue)/(1-rvalue) )
  return(rvalue)
}




###两组数据没有显著的差异性
z1<-gene_corZvalue("thickness",inter1,1)#t-test result 
z2<-gene_corZvalue("thickness",connectomeAge1,1)#r=0.745 
t.test(z1,z2,paired = FALSE) #p=0.04
z1<-gene_corZvalue("thickness",inter2,2)#t-test result 
z2<-gene_corZvalue("thickness",connectomeAge2,2)#r=0.745 
t.test(z1,z2,paired = FALSE) #p=6.861e-05
 z1<-gene_corZvalue("thickness",inter3,2)#t-test result 
 z2<-gene_corZvalue("thickness",connectomeAge3,2)#r=0.745 p=0.011
 t.test(z1,z2,paired = FALSE)
boxplot(z1,z2)
##



##只在第二个年龄段上是两组数据是显著不同的

####
z1<-gene_corZvalue("SmoothedMyelinMap",inter1,1)#t-test result 
z2<-gene_corZvalue("SmoothedMyelinMap",connectomeAge1,1)#r=0.745 p=0.011
t.test(z1,z2,paired = FALSE)
z1<-gene_corZvalue("SmoothedMyelinMap",inter3,3)#t-test result 
z2<-gene_corZvalue("SmoothedMyelinMap",connectomeAge3,3)#r=0.745 p=0.011
t.test(z1,z2,paired = FALSE)
###########



###
####test the result in structure
gene_cor("thickness",inter1,1)#r=0.7909091 p=.006
gene_cor("thickness",connectomeAge1,1)#r=0.745 p=0.011  


gene_cor("thickness",inter2,2)#r=0.709 p=0.01
gene_cor("thickness",connectomeAge2,2)#r=0.5727 p=0.07
gene_cor("thickness",inter3,3)#r=0.718 p=0.017
gene_cor("thickness",connectomeAge3,3)#r=0.7 p=0.02
##直接看两组r值进行t。test是不显著的
gene_cor("SmoothedMyelinMap",inter1,1)#-0.945  0
gene_cor("SmoothedMyelinMap",connectomeAge1,1)#-0.92727 0
gene_cor("SmoothedMyelinMap",inter2,2)#-0.9182 0
gene_cor("SmoothedMyelinMap",connectomeAge2,2)#-0.891 0.0004
gene_cor("SmoothedMyelinMap",inter3,3)#-0.9364 0
gene_cor("SmoothedMyelinMap",connectomeAge3,3)#-0.982 0

#######
grep("LPNH3",connectomeAge3,fixed = TRUE)
intersect("LPHN3",connectomeAge2)####而且LPHN3是在connectome相关的基因里面的
###########Disease test   
diseaseTest<-function(disease_name,number){
  #disease_name<-"AUTgene"
  disease_info <- read.table(paste("F://project/validation/diease_gene/",disease_name,"/",tolower(disease_name),".txt",sep=""),header=F,sep="\t",stringsAsFactors = F)
  staM<-list()
  if(number==1){
    numberT1<-intersect(inter1,disease_info$V1)
    numberT2<- intersect(connectomeAge1,disease_info$V1)
    staM[[1]]<-numberT1
    staM[[2]]<-numberT2
  }else if(number==2){
    numberT1<-intersect(inter2,disease_info$V1)
    numberT2<-intersect(connectomeAge2,disease_info$V1)
    staM[[1]]<-numberT1
    staM[[2]]<-numberT2
  }else if(number==3){
    numberT1<- intersect(inter3,disease_info$V1)
    numberT2<-intersect(connectomeAge3,disease_info$V1)
    staM[[1]]<-numberT1
    staM[[2]]<-numberT2
  }
  return(staM)
}
##conclusion：
Pd1<-diseaseTest("PDgene",2)##4 4
Pd2<-diseaseTest("PDgene",3) #7 9 ##PD在第二个年龄段是与connectome gene相关的多，在第三个年龄段两组基因相交的数量是一样的
Pd3<-diseaseTest("PDgene",1)
###7   15
diseaseTest("AUTgene",1)#4 8 #ASD第一个年龄段是与connectome相关的多，第三个年龄段是与差异表达的基因相关的多
diseaseTest("AUTgene",3)#9 5
diseaseTest("AUTgene",2)#9 5
##13   13
diseaseTest("Depression",2)# 0 4 #只与connectome相关有交集
diseaseTest("ALzgene",3)# 6 4  #与特异的相交的多
###6 4
diseaseTest("SZgene",2)#4 10 #与connectome相交的多
diseaseTest("SZgene",1)#7 5
diseaseTest("SZgene",3)#9 7
####19  22
diseaseTest("SZgene",3)#9 7
diseaseTest("ADHD",2)# 0 2 #只有与connectome的有相交的多
##0 2
ageGroup<-c(rep("PDgene",2),rep("AUTgene",2),rep("Depression",2),rep("ALzgene",2),rep("SZgene",2),rep("ADHD",2))
condition<-rep(c("DSgene","connectome"),3)
number<-c(7,15,12,14,2,2,6,4,9,7,0,2)
plotData<-data.frame(ageGroup,condition,number)
ggplot(plotData,aes(x=ageGroup,y=number,fill=condition))+
  geom_bar(stat="identity")+
  geom_text(aes(label=number),position=position_stack(vjust=0.5))#3 vjust参数用来调整标签的为重，vjust=0.5将标签放到对应部位的中部
##

###connectome的相关性
data.fmri.sample.cor.network <- function(file_path,column_index){
  
  #file_path<-"parcellations_VGD11b"
  #column_index<-2
  
  data.rest1.bold.LR <- read.table(paste("F://project/newkang_brainNet/report20160227/",file_path,"/fmri/REST1/COLUMN_",column_index,"_rfMRI_REST1_LR_Atlas_hp2000_clean.dtseries.txt",sep=""),header=F)
  data.rest1.bold.RL <- read.table(paste("F://project/newkang_brainNet/report20160227/",file_path,"/fmri/REST1/COLUMN_",column_index,"_rfMRI_REST1_RL_Atlas_hp2000_clean.dtseries.txt",sep=""),header=F)
  #data.rest1.bold <- cbind(data.rest1.bold.LR,data.rest1.bold.RL)
  data.rest1.bold <- cbind(data.rest1.bold.RL)
  data.column.L <- read.table(paste("F://project/newkang_brainNet/report20160227/",file_path,"/fmri/COLUMN_",column_index,"_L.txt",sep=""))
  data.column.R <- read.table(paste("F://project/newkang_brainNet/report20160227/",file_path,"/fmri/COLUMN_",column_index,"_R.txt",sep=""))
  each <- nrow(data.column.L) + nrow(data.column.R) + 19
  
  data.label.L <- read.table(paste("F://project/newkang_brainNet/report20160227/",file_path,"_L.txt",sep=""))
  data.label_to_11.L <- read.table(paste("F://project/newkang_brainNet/report20160227/",file_path,"_L_1.txt",sep=""))
  data.sample.bold.network <- matrix(NA,27,55)
  for(i in c(1:27)){
    data.sampel.bold <- data.rest1.bold[(each*(i-1)+1):(each*i),]
    data.sample.bold.L.11 <- matrix(NA,11,ncol(data.sampel.bold))
    for(j in c(1:11)){
      data.roi.label <- as.matrix(data.label_to_11.L[which(data.label_to_11.L[,2]==j),1])
      label.index <- apply(data.roi.label,1,function(x){
        return(which(data.label.L[,3]==x))
      })
      data.sample.bold.L.11[j,] <- apply(t(data.sampel.bold[label.index,]),1,mean)
    }
    
    data.cor <- cor(t(data.sample.bold.L.11),method="pearson")
    #
    data.sample.bold.network[i,] <- data.cor[lower.tri(data.cor)]
  }
  
 # data.sample.bold.network<-atanh(data.sample.bold.network)
  
  return(data.sample.bold.network)
  
}
#1. brain networks
## For parcellations_Vgd11b COLUMN_2 Label
data.sample.ba <- data.fmri.sample.cor.network("parcellations_VGD11b",2)
#data.sample.ba <- data.fmri.sample.cor.network("RSN-networks",2)
#For five age group:|11|12|13|14|15|
data.label.fmri.cor.5ageGroup.p1 <- apply(data.sample.ba[1:12,],2,mean)
data.label.fmri.cor.5ageGroup.p2 <- apply(data.sample.ba[13:17,],2,mean)
data.label.fmri.cor.5ageGroup.p3 <- apply(data.sample.ba[18:27,],2,mean)
####画一下这个关联矩阵的
write.table(data.label.fmri.cor.5ageGroup.p1,"age1.txt")
write.table(data.label.fmri.cor.5ageGroup.p2,"age2.txt")
write.table(data.label.fmri.cor.5ageGroup.p3,"age3.txt")
###
drawThecor<-function(geneExp){
  data1<-matrix(NA,11,11)
  up1<-upper.tri(data1)
  data1[upper.tri(data1)]<-geneExp
  data1[lower.tri(data1)]<-geneExp
  diag(data1)<-1
  colnames(data1)<-colnames(data.gene.roi.mean.5ageGroup_NM.p1)
  rownames(data1)<-colnames(data.gene.roi.mean.5ageGroup_NM.p1)
  nameOrder<-c("ITC","A1C","STC","MFC","V1C","M1C","S1C","DFC","IPC","OFC","VFC")
  data1<-data1[nameOrder,nameOrder]
 # heatmap(data1)
  heatmap.2(data1,scale ="none",
            dendrogram = "col",
            Colv =FALSE,
            Rowv = FALSE,
            trace ="none",
            col=rev(colorRampPalette(brewer.pal(10, "RdBu"))(10)),
            breaks = seq(0,1,0.1),
            density.info=c("none"),
            cexRow =1.2, cexCol = 1.2   # 行/列名字体大小
  )
}
drawThecor(data.label.fmri.cor.5ageGroup.p1)
drawThecor(data.label.fmri.cor.5ageGroup.p2)
drawThecor(data.label.fmri.cor.5ageGroup.p3)

NM_info <- read.table('F://project/newkang_brainNet/kang/new_code_document/Gene_Annotation/NM_info.txt',header =T,row.names=1,sep='\t')
NM_info<-as.matrix(NM_info)
######build the co-expression data gene在 11 个脑区上的表达值
age1inter<-data.gene.roi.mean.5ageGroup_NM.p1[inter1,]
age1connect<-data.gene.roi.mean.5ageGroup_NM.p1[connectomeAge1,]
age2RoidataMean<-data.gene.roi.mean.5ageGroup_NM.p2[age2Cor_geneSymbol,]
age3RoidataMean<-data.gene.roi.mean.5ageGroup_NM.p3[age3Cor_geneSymbol,]
######co-expression
co_expression<-cor(t(age1inter))
up1<-upper.tri(co_expression)
up_matrix<-co_expression[up1]








##
library(stringi)
library(stringr)
findTheageRelatedtrans<-function(allgene){
  #allgene<-inter1
  #number<-list()
  allgeneVector<-list()
  for(i in c(1:length(allgene))){
   # i=2
    listNum<-vector()
    temp<-allgene[i]
    for(j in c(1:nrow(NM_info))){
      if(as.character(NM_info[j,])==temp){
        #print(j)
        listNum<-c(listNum,j)
      }
    }
    
    nameTrans<-rownames(NM_info)[listNum]
    allgeneVector[[i]]<-as.matrix(nameTrans)
  }
  allgeneVector<-unlist(allgeneVector)%>%unique()
  return(allgeneVector)
}
#####to the theory
age1Network1<-findTheageRelatedtrans(inter1)
age1Network2<-findTheageRelatedtrans(connectomeAge1)
age2Network1<-findTheageRelatedtrans(inter2)
age2Network2<-findTheageRelatedtrans(connectomeAge2)
age3Network1<-findTheageRelatedtrans(inter3)
age3Network2<-findTheageRelatedtrans(connectomeAge3)
#####extract network
inter1Network<-age1transNet[age1Network1,]%>%as.matrix()
connect1Network<-age1transNet[age1Network2,]%>%as.matrix()
inter2Network<-age2transNet[age2Network1,]%>%as.matrix()
connect2Network<-age2transNet[age2Network2,]%>%as.matrix()
inter3Network<-age3transNet[age3Network1,]%>%as.matrix()
connect3Network<-age3transNet[age3Network2,]%>%as.matrix()

########test the correlation
corTest1<-apply(inter3Network,1,function(x){cor(x,data.label.fmri.cor.5ageGroup.p3,method ="pearson")})
corTest2<-apply(connect3Network,1,function(x){cor(x,data.label.fmri.cor.5ageGroup.p3,method ="pearson")})
#tranform r to z
rtoZvalue<-function(value){
  rvalueZ<-0.5 * log( (1+value)/(1-value) )
  return(rvalueZ)
}
corTest1<-rtoZvalue(corTest1)
corTest2<-rtoZvalue(corTest2)
t.test(corTest1,corTest2,paired = FALSE)
########second pca analysis####
pcaAnalysis<-function(interData,ConnectData){
  pca1<-prcomp(t(interData),scale. = FALSE)
  pca2<-prcomp(t(ConnectData),scale. = FALSE)
  #percentVar<-pca1$sdev^2/sum(pca1$sdev^2)
  data1<-pca1$x[,1:13]
  data2<-pca2$x[,1:13]
  corTest11<-apply(data1,2,function(x){cor(x,data.label.fmri.cor.5ageGroup.p3,method ="pearson")})
  corTest22<-apply(data2,2,function(x){cor(x,data.label.fmri.cor.5ageGroup.p3,method ="pearson")})
  print(t.test(corTest11,corTest1))
}

###enrichment###
KEGG_GO <- function(annotation,pvalueCutoff,genes){
  
  library("org.Hs.eg.db")
  library("GSEABase")
  library("GOstats")
  #genes <- c("AREG", "FKBP5", "CXCL13", "KLF9", "ZC3H12A", "P4HA1", "TLE1", "CREB3L2", "TXNIP", "PBX1", "GJA1", "ITGB8", "CCL3", "CCND2", "KCNJ15", "CFLAR", "CXCL10", "CYSLTR1", "IGFBP7", "RHOB", "MAP3K5", "CAV2", "CAPN2", "AKAP13", "RND3", "IL6ST", "RGS1", "IRF4", "G3BP1", "SEL1L", "VEGFA", "SMAD1", "CCND1", "CLEC3B", "NEB", "AMD1", "PDCD4", "SCD", "TM2D3", "BACH2", "LDLR", "BMPR1B", "RFXAP", "ASPH", "PTK2B", "SLC1A5", "ENO2", "TRPM8", "SATB1", "MIER1", "SRSF1", "ATF3", "CCL5", "MCM6", "GCH1", "CAV1", "SLC20A1")
  
  ##GO BP enrichment analysis
  goAnn <- get("org.Hs.egGO")
  universe <- Lkeys(goAnn)
  
  genes <- na.omit(genes)
  entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
  entrezIDs <- as.character(entrezIDs)
  entrezIDs <- na.omit(entrezIDs)
  #index.entrez <- duplicated(entrezIDs)
  #entrezIDs.fliter <- as.vector(entrezIDs[which(index.entrez==FALSE)])
  
  params <- new("GOHyperGParams",
                geneIds=entrezIDs,
                universeGeneIds=universe,
                annotation="org.Hs.eg.db",
                ontology=annotation,
                pvalueCutoff=pvalueCutoff,
                conditional=FALSE,
                testDirection="over")
  over <- hyperGTest(params)
  
  library(Category)
  glist <- geneIdsByCategory(over)
  glist <- sapply(glist, function(.ids) {
    .sym <- mget(.ids, envir=org.Hs.egSYMBOL, ifnotfound=NA)
    .sym[is.na(.sym)] <- .ids[is.na(.sym)]
    paste(.sym, collapse=";")
  })
  
  bp <- summary(over)
  
  if(annotation=="BP"){
    bp$Symbols <- glist[as.character(bp$GOBPID)]
  }else if(annotation=="MF"){
    bp$Symbols <- glist[as.character(bp$GOMFID)]
  }else if(annotation=="CC"){
    bp$Symbols <- glist[as.character(bp$GOCCID)]
  }
  return(bp)
}

function.GO.analysis <- function(gene_symbol){
 # gene_symbol<-inter1
  index.gene <- duplicated(gene_symbol)
  genes <- as.vector(gene_symbol[which(index.gene==FALSE)])
  
  #KEGG_GO <- function(annotation,pvalueCutoff,genes)
  
  KEGG_GO.BP <- KEGG_GO("BP",1,genes)
  #KEGG_GO.BP2<-KEGG_GO.BP[KEGG_GO.BP$Size>200,]
  #KEGG_GO.BP2$Pvalue<-p.adjust(KEGG_GO.BP2$Pvalue,method = "fdr",n=length(KEGG_GO.BP2$Pvalue))
  # KEGG_GO.BP1<-KEGG_GO.BP[which(KEGG_GO.BP$Size>200),]
 # fdr_keeg<-KEGG_GO.BP2[which(KEGG_GO.BP2$Pvalue<=0.05),]
  #resultConnect1$Pvalue<-p.adjust(resultConnect1$Pvalue,method = "fdr",n=length(resultConnect1$Pvalue))
  # return(KEGG_GO.BP)
  return(KEGG_GO.BP)
  #KEGG_GO.MF <- KEGG_GO("MF",0.05,genes)
  #write.table(KEGG_GO.MF,paste('./kang/new-BA11/Gene_GO/',parent_path,'/specific_gene_0',i,'/',name.roi1.roi2.gene,'.KEGG_GO.MF.diff.txt',sep=''),row.names = F,col.names = T,quote = F,sep='\t')
  
  #KEGG_GO.CC <- KEGG_GO("CC",0.05,genes)
  #write.table(KEGG_GO.CC,paste('./kang/new-BA11/Gene_GO/',parent_path,'/specific_gene_0',i,'/',name.roi1.roi2.gene,'.KEGG_GO.CC.diff.txt',sep=''),row.names = F,col.names = T,quote = F,sep='\t')
}
#########gene enrich analysis
resultInter1<- function.GO.analysis(inter1)
resultInter2<- function.GO.analysis(inter2)
resultInter3<- function.GO.analysis(inter3)
resultConnect1<- function.GO.analysis(connectomeAge1)
resultConnect2<- function.GO.analysis(connectomeAge2)
resultConnect3<- function.GO.analysis(connectomeAge3)
resultInter1Filet<-resultInter1[resultInter1$Size>200,]
resultInter2Filet<-resultInter2[resultInter2$Size>200,]
resultInter3Filet<-resultInter3[resultInter3$Size>200,]
resultConnect1Filet<-resultConnect1[resultConnect1$Size>200,]
resultConnect2Filet<-resultConnect2[resultConnect2$Size>200,]
resultConnect3Filet<-resultConnect3[resultConnect3$Size>200,]
###correct the result 
resultInter1Filet$Pvalue<-p.adjust(resultInter1Filet$Pvalue,method = "fdr",n=nrow(resultInter1Filet))
resultInter2Filet$Pvalue<-p.adjust(resultInter2Filet$Pvalue,method = "fdr",n=nrow(resultInter2Filet))
resultInter3Filet$Pvalue<-p.adjust(resultInter3Filet$Pvalue,method = "fdr",n=nrow(resultInter3Filet))
resultConnect1Filet$Pvalue<-p.adjust(resultConnect1Filet$Pvalue,method = "fdr",n=nrow(resultConnect1Filet))
resultConnect2Filet$Pvalue<-p.adjust(resultConnect2Filet$Pvalue,method = "fdr",n=nrow(resultConnect2Filet))
resultConnect3Filet$Pvalue<-p.adjust(resultConnect3Filet$Pvalue,method = "fdr",n=nrow(resultConnect3Filet))
#####exact the end result
resultInter1End<-resultInter1Filet[which(resultInter1Filet$Pvalue<=0.2),]
resultInter2End<-resultInter2Filet[which(resultInter2Filet$Pvalue<=0.2),]
resultInter3End<-resultInter3Filet[which(resultInter3Filet$Pvalue<=0.2),]
resultConnect1End<-resultConnect1Filet[which(resultConnect1Filet$Pvalue<=0.2),]
resultConnect2End<-resultConnect2Filet[which(resultConnect2Filet$Pvalue<=0.2),]
resultConnect3End<-resultConnect3Filet[which(resultConnect3Filet$Pvalue<=0.2),]
dir.create("CorrectEnrich")
write.csv(resultInter1End,"inter1CorEnrich.csv")
write.csv(resultInter2End,"inter2CorEnrich.csv")
write.csv(resultInter3End,"inter3CorEnrich.csv")

#####keeg result ##


###########把所有的gene放在一起

plotEnrichres<-function(geneEnrich){
  library(ggplot2)
  #rank
  resultInter1<-resultInter1[order(resultInter1$Pvalue,decreasing =TRUE),]
  term<-resultInter1$Term
  resultInter1$Term<-factor(resultInter1$Term,levels =resultInter1$Term )
  p = ggplot(resultInter1,aes(OddsRatio,Term))
  p=p + geom_point()
  # change the size of point
  p=p + geom_point(aes(size=Count))
  # show
  pbubble = p + geom_point(aes(size=Count,color=Pvalue))
  # plot the pathway result
  pr = pbubble + scale_colour_gradient(low="red",high="green") + 
    labs(color=expression(Pvalue),size="Gene number",
         x="GeneRatio",y="Pathway name",title="Top20 of pathway enrichment")
  print(pr)
}

# 改变图片的样式（主题）
age1enrich<-read.csv("plotage1.csv",header = TRUE,stringsAsFactors = FALSE)
age2enrich<-read.csv("plotenrich2.csv",header = TRUE,stringsAsFactors = FALSE)
age3enrich<-read.csv("plotage3.csv",header = TRUE,stringsAsFactors = FALSE)
plotEnrichres(age2enrich)
plotEnrichres(age3enrich)


#####look the TFs in two group
interTfSet1<-fisherTrs_res(inter1)
interTfSet11<-interTfSet1[which(interTfSet1<=0.05),]
ConnectTfSet1<-fisherTrs_res(connectomeAge1)
ConnectTfSet11<-ConnectTfSet1[which(ConnectTfSet1<=0.05),]
interTfSet2<-fisherTrs_res(inter2)
interTfSet2<-fisherTrs_res(inter2)

ConnectTfSet1<-fisherTrs_res(connectomeAge1)
ConnectTfSet1<-fisherTrs_res(connectomeAge1)
tese1<-KEEG_pathway(names(interTfSet11))
tese2<-KEEG_pathway(names(ConnectTfSet11))
barplot(tese1,showCategory = 30)  
barplot(tese2,showCategory = 30)  
###gene enrichment
ss<-c('MEF2C','HNRNPH3','FOXP4','DLX1', 'FGF2','LHX6','LEPTIN','SMAD6','DBX2','NAV1','MPZL2')
brainDe<-read.table('brainDevelopment')
intersect(brainDe$V1,Tf1$x)
intersect(brainDe$V1,Tf2$x)
intersect(brainDe$V1,Tf3$x)
intersect(brainDe$V1,age1Cor_geneSymbol)
intersect(brainDe$V1,age2Cor_geneSymbol)
intersect(brainDe$V1,age3Cor_geneSymbol)
#
intersect(ss,as.character(age1Cor_geneSymbol))
intersect(ss,as.character(age2Cor_geneSymbol))
intersect(ss,as.character(age3Cor_geneSymbol))


#
Tf1<-read.table('Age1AllTFAfterCorrect.txt')
Tf2<-read.table('Age2AllTFAfterCorrect.txt')
Tf3<-read.table('Age3AllTFAfterCorrect.txt')
