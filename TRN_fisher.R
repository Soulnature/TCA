###my geneset in TRN network
library(dplyr)
TrnModel<-read.csv("TRN.csv",header = TRUE,stringsAsFactors = FALSE)
TFs<-TrnModel$X...node1%>%unique()
targetGene<-TrnModel$node2%>%unique()
#####all gene
allGeneSymbol<-Reduce(union,list(age1Cor_geneSymbol,age2Cor_geneSymbol,age3Cor_geneSymbol))
##
interTFs<-intersect(allGeneSymbol,TFs)###17个转录因子
#[1] "NFIB"    "TSHZ2"   "MAFK"    "NPAS2"   "RFX3"    "STAT4"   "RXRG"    "PKNOX2"  "SCRT2"   "BHLHE40" "FOSL2"   "POU3F2" 
#[13] "NEUROD6" "ESRRG"   "CREBL2"  "ARNT2"   "KLF16"  
##统计下边的数目
findTrs<-TrnModel[which(TrnModel$X...node1==interTFs),]
#############
intersect(findTrs$node2%>%unique(),allGeneSymbol)
degree<-read.table("outdegree")
allEdgeW<-intersect(targetGene,allGeneSymbol)######转录调控网络里面有259个node是重叠的
########take the first age group
TFs_age1<-intersect(TFs,age1Cor_geneSymbol)
TFs_age2<-intersect(TFs,age2Cor_geneSymbol)
TFs_age3<-intersect(TFs,age3Cor_geneSymbol)
########extract the subclass,and analysis the class
setAgeOne<-TrnModel[which(TrnModel$X...node1==TFs_age1),]
AgeOnegroupTrs<-union(unique(setAgeOne$node2),TFs_age1)
setAgetwo<-TrnModel[which(TrnModel$X...node1==TFs_age2),]
AgetwogroupTrs<-union(unique(setAgeOne$node2),TFs_age2)
setAgeThree<-TrnModel[which(TrnModel$X...node1==TFs_age3),]
AgeThreegroupTrs<-union(unique(setAgeThree$node2),TFs_age3)
###fisher exact test for the our gene set and Trs model set
age.fisher.test <- function(netage_gene,anova){
  # netage_gene=age1Cor_geneSymbol
  # anova=tempTar
  data.NM_info.label.cor.3ageGroup.p1<-read.table("/Users/zhaoxingzhong/Desktop/project/newkang_brainNet/geneNet/Full_trans_Net/age1.txt")
  NM_info1 <- unique(NM_info[rownames(data.NM_info.label.cor.3ageGroup.p1),])
  intersect.anova<- intersect(anova,NM_info1)
  ##----disease-----------------------------------
  all.NM <- unique(NM_info1);
  all.diseaseg <- unique(anova)
  intersect.NM.is.diseaseg <- intersect(all.NM,all.diseaseg)
  intersect.NM.not.diseaseg <- setdiff(all.NM,all.diseaseg)
  #Ag-Ag' Ag-not Ag'
  intersect.diseaseg.is.5ageGroup.p1 <- intersect(intersect.NM.is.diseaseg,netage_gene)
  #print(intersect.diseaseg.is.5ageGroup.p1)
  intersect.diseaseg.not.5ageGroup.p1 <- setdiff(intersect.NM.is.diseaseg,netage_gene)
  #notAg-notAg' notAg-notnotAg'
  intersect.ndiseaseg.is.5ageGroup.p1 <- intersect(intersect.NM.not.diseaseg,netage_gene)
  intersect.ndiseaseg.not.5ageGroup.p1 <- setdiff(intersect.NM.not.diseaseg,netage_gene)
  
  result.disease.NM_info.label.cor.5ageGroup.p1 <- rbind(c(length(intersect.diseaseg.is.5ageGroup.p1),length(intersect.diseaseg.not.5ageGroup.p1)),c(length(intersect.ndiseaseg.is.5ageGroup.p1),length(intersect.ndiseaseg.not.5ageGroup.p1)))
 # print(result.disease.NM_info.label.cor.5ageGroup.p1)
  res1<-fisher.test(result.disease.NM_info.label.cor.5ageGroup.p1,alternative ="greater")
  return(res1$p.value)
}
fisherTrs_res<-function(ageCorgeneT){
  ageSavematrix<-matrix(NA,length(TFs),1)
  rownames(ageSavematrix)<-TFs
  for(i in c(1:length(TFs))){
    temp<-TFs[i]
    tempTar<-TrnModel[TrnModel$X...node1==temp,]$node2%>%unique()
    ageSavematrix[i]<-age.fisher.test(ageCorgeneT,tempTar)
  }
  return(ageSavematrix)
}
age1_Fisher_Res<-fisherTrs_res(age1Cor_geneSymbol)
age2_Fisher_Res<-fisherTrs_res(age2Cor_geneSymbol)
age3_Fisher_Res<-fisherTrs_res(age3Cor_geneSymbol)
####
age1_adjust<-p.adjust(age1_Fisher_Res,method = "BH")
age2_adjust<-p.adjust(age2_Fisher_Res,method = "BH")
age3_adjust<-p.adjust(age3_Fisher_Res,method = "BH")
#age1_Fisher_Res<-age1_adjust
names(age1_adjust)<-rownames(age1_Fisher_Res)
names(age2_adjust)<-rownames(age2_Fisher_Res)
names(age3_adjust)<-rownames(age3_Fisher_Res)
Age1SignificantCor<-age1_adjust[which(age1_adjust<=0.05)]
Age2SignificantCor<-age2_adjust[which(age2_adjust<=0.05)]
Age3SignificantCor<-age3_adjust[which(age3_adjust<=0.05)]
dir.create('Fisher_TFs')
write.table(as.matrix(Age1SignificantCor),'age1FisherTfs.txt',quote = FALSE)
write.table(as.matrix(Age2SignificantCor),'age2FisherTfs.txt',quote = FALSE)
write.table(as.matrix(Age3SignificantCor),'age3FisherTfs.txt',quote = FALSE)
###### Without correct ####
Age1Significant<-age1_Fisher_Res[age1_Fisher_Res<=0.05,]#52
Age2Significant<-age2_Fisher_Res[age2_Fisher_Res<=0.05,]#49
Age3Significant<-age3_Fisher_Res[age3_Fisher_Res<=0.05,]#71
write.table(Age1Significant,"Age1SignificantWithoutCor.txt")
write.table(Age2Significant,"Age2SignificantWithoutCor.txt")
write.table(Age3Significant,"Age3SignificantWithoutCor.txt")
interallSe1<-intersect(names(Age1Significant),unique(c(names(Age2Significant)),names(Age3Significant)))
######enrichment analysis####
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
enrich<-function(gene){
  data(geneList, package="DOSE")
  #gene<-names(Age2Significant)  
  test1 = bitr(gene, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")
  ego_ALL <- enrichGO(gene = test1$ENTREZID, 
                      universe = names(geneList), #背景基因集
                      OrgDb = org.Hs.eg.db, #没有organism="human"，改为OrgDb=org.Hs.eg.db
                      #keytype = 'ENSEMBL',
                      ont = "BP", #也可以是 CC  BP  MF中的一种
                      pAdjustMethod = "fdr", #矫正方式 holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”中的一种
                      pvalueCutoff = 0.05, #P值会过滤掉很多，可以全部输出
                      qvalueCutoff = 1,
                      readable = TRUE)
  return(ego_ALL)
}
age1TfsEnrich<-enrich(names(Age1Significant))
age2TfsEnrich<-enrich(names(Age2Significant))
age3TfsEnrich<-enrich(names(Age3Significant))
#TFsAgeOne<-function.GO.analysis()plot the result
barplot(age1TfsEnrich,showCategory = 30,title="EnrichmentGO_BP_top30")
barplot(age2TfsEnrich,showCategory = 30,title="EnrichmentGO_BP_top30")
barplot(age3TfsEnrich,showCategory = 30,title="EnrichmentGO_BP_top30")
####
re11<-enrich(age1Cor_geneSymbol)
##

#####对这三个set进行一下富集分析，看下结果
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
data(geneList, package="DOSE") #富集分析的背景基因集
transfrom<-bitr(connectomeAge1,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

########对两组基因做KEEG
TrnEnrichage1<-KEEG_pathway(names(Age1Significant))
TrnEnrichage2<-KEEG_pathway(names(Age2Significant))
TrnEnrichage3<-KEEG_pathway(names(Age3Significant))
GrnEnrichage1<-KEEG_pathway(names(Age1SignificantGrn))
GrnEnrichage2<-KEEG_pathway(names(Age2SignificantGrn))
GrnEnrichage3<-KEEG_pathway(names(Age3SignificantGrn))
barplot(GrnEnrichage1,showCategory = 30)
barplot(GrnEnrichage2,showCategory = 30)
barplot(GrnEnrichage3,showCategory = 30)
########kEEG intersection
interAge1Tfs<-intersect(names(Age1Significant),names(Age1SignificantGrn))
interAge2Tfs<-intersect(names(Age2Significant),names(Age2SignificantGrn))
interAge3Tfs<-intersect(names(Age3Significant),names(Age3SignificantGrn))
##########union set test by correct gene set###
ageOneCorunion<-union(names(Age1SignificantCor),names(Age1SignificantCorGrn))#18
ageTwoCorunion<-union(names(Age2SignificantCor),names(Age2SignificantCorGrn))#81
ageThreeCorunion<-union(names(Age3SignificantCor),names(Age3SignificantCorGrn))#52
###将找到的TFs写入到文件去
write.table(ageOneCorunion,'Age1AllTFAfterCorrect.txt')
write.table(ageTwoCorunion,'Age2AllTFAfterCorrect.txt')
write.table(ageThreeCorunion,'Age3AllTFAfterCorrect.txt')
###
#######keeg enrichment###
unionAge1<-KEEG_pathway(ageOneCorunion)
unionAge2<-KEEG_pathway(ageTwoCorunion)
unionAge3<-KEEG_pathway(ageThreeCorunion)
barplot(unionAge1,showCategory = 30,title="EnrichmentGO_BP_top30")
barplot(unionAge2,showCategory = 30,title="EnrichmentGO_BP_top30")
barplot(unionAge3,showCategory = 30,title="EnrichmentGO_BP_top30")
###GO enrichment
unionAge1Go<-enrich(ageOneCorunion)
unionAge2Go<-enrich(ageTwoCorunion)
unionAge3Go<-enrich(ageThreeCorunion)
barplot(unionAge1Go,showCategory = 30,title="Ageonethe top 30 BP term")
barplot(unionAge2Go,showCategory = 30,title="AgeTwothe top 30 BP term")
barplot(unionAge3Go,showCategory = 30,title="AgeThreethe top 30 BP term")

######disease analysis#####
disease_TFs.overlap.test <- function(disease_name,target_set){
  #disease_name<-"PDgene"
  disease_info <- read.table(paste("/Users/zhaoxingzhong/Desktop/project/validation/diease_gene/",disease_name,"/",tolower(disease_name),".txt",sep=""),header=F,sep="\t",stringsAsFactors = F)
  print(nrow(disease_info))
  intersect(disease_info$V1,target_set)
}

disease_TFs.overlap.test("AUTgene",keyGene)
disease_TFs.overlap.test("Depression",keyGene)
disease_TFs.overlap.test("ALzgene",keyGene)
disease_TFs.overlap.test("SZgene",keyGene)
disease_TFs.overlap.test("PDgene",keyGene)
disease_TFs.overlap.test("ADHD",keyGene)
####算下fisher test 
