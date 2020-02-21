####Gene Ratio  

##the formula of computing connections ratio Ratio=(C_in/C_all)/(N_out/N_all)|
#|anotation:{C_in 代表的是子模块中的连接数.C_out代表的是模块中的所有的节点的数量.N_out代表的是模块中与外部的连接的数量
#| N_out代表的是外部的所有的节点数量}
###
statisticGeneratio<-function(ChosegeneList,TrnModel){
  #ChosegeneList<-allgene
  #allGeneSymbol<-allgene
  #ChosegeneList<-age1Cor_geneSymbol
  #allGeneSymbol<-age1Cor_geneSymbol
  allnetworkData<-vector()
  allTargetgene<-TrnModel$node2%>%unique()
  otherGene<-setdiff(allTargetgene,ChosegeneList)
  for (gene in ChosegeneList) {
    #tem=ChosegeneList[1]
    needData<-TrnModel[which(TrnModel$node2==gene),]
    allnetworkData<-rbind(allnetworkData,needData)
  }
  ##compute the ratio,extract the subclass node
  theGeneinTrn<-allnetworkData$node2%>%unique()
  geneFind<-apply(allnetworkData, 1, function(x){
   # x<-allnetworkData[1,]
    temJ<-intersect(theGeneinTrn,x[1])
    if(length(temJ)>0){
      return(1)
    }else{
      return(0)
    }
  })
  inEdge<-geneFind[which(geneFind==1)]%>%length()
  outEdge<-geneFind[which(geneFind==0)]%>%length()
  ratio<-(inEdge/length(theGeneinTrn))/(outEdge/length(otherGene))###1.49215内部连接是大于外部连接的
  return(ratio)
}


TrnModel<-read.csv("/Users/zhaoxingzhong/Desktop/project/TRS/TRN.csv",header = TRUE,stringsAsFactors = FALSE)
age1ratio<-statisticGeneratio(age1Cor_geneSymbol,TrnModel)
age2ratio<-statisticGeneratio(age2Cor_geneSymbol,TrnModel)
age3ratio<-statisticGeneratio(age3Cor_geneSymbol,TrnModel)
#GrnModel<-read.csv("/Users/zhaoxingzhong/Desktop/project/TRS/INT-11_ElasticNet_Filtered_Cutoff_0.1_GRN_1 (1).csv",header = TRUE,stringsAsFactors = FALSE)
#GRN  0.6319719
###1.062036 mean
###new 1.373375
###trnmodel
trnmodel<-TrnModel[,c(1,2,4)]
write.table(trnmodel,'trnmodel.txt',sep = '\t',row.names = FALSE,quote = FALSE)
TFs<-TrnModel$X...node1
 intersect(TFs,'PKNOX1') 
 