library(Rmisc)
library(export)
getDataFromregion<-function(set1,set2)
{
  saveMatrix<-matrix(NA,11,2)
  #set1<-age1RoidataMean
  #set2<-data.gene.roi.mean.5ageGroup_NM.p1
  set1First<-apply(set1, 2, mean)
  set2Second<-apply(set2,2,mean)
  print(shapiro.test(set1First))
  print(shapiro.test(set2Second))
  #print(wilcox.test(set1,set2))
  #print(t.test(set1First,set2Second,paired = TRUE))
  print(mean(set1First)/mean(set2Second))
  groupallRegion<-cbind(set1First,set2Second)
  saveMatrix[,1]<-set1First
  saveMatrix[,2]<-set2Second
 # saveMatrix[,1]<-apply(groupallRegion, 2, mean)
  #saveMatrix[,2]<-apply(groupallRegion, 2, sd)
  return(saveMatrix)
}

age1<-getDataFromregion(age1RoidataMean,data.gene.roi.mean.5ageGroup_NM.p1)
age2<-getDataFromregion(age2RoidataMean,data.gene.roi.mean.5ageGroup_NM.p2) 
age3<-getDataFromregion(age3RoidataMean,data.gene.roi.mean.5ageGroup_NM.p3)  
            ###柱状图
# allage<-rbind(age1,age2,age3) 
# xname<-c(rep("8-20",2),rep("20-40",2),rep("40-",2))
# group_split<-rep(c('ourGeneset','allGeneset'),3)
# plotDframall<- data.frame(x=xname,y=allage[,1],sd=allage[,2],supDivide=group_split) 
# plotDframall$x<-factor(plotDframall$x)
#  ggplot(plotDframall, aes(x=x, y=y, fill=supDivide)) + 
#     geom_bar(position=position_dodge(), stat="identity",
#              colour="black", # Use black outlines,
#              size=.3) +      # Thinner lines
#     geom_errorbar(aes(ymin=y-sd, ymax=y+sd),
#                   size=.3,    # Thinner lines
#                   width=.2,
#                   position=position_dodge(.9)) +
#     xlab("brainRegion") +
#     ylab("Gene Expression") +
#     scale_fill_hue(name="Supplement type", # Legend label, use darker colors
#                    breaks=c("ourGeneset", "allGeneset"),
#                    labels=c("ourGeneset", "allGeneset")) +
#     scale_y_continuous(breaks=0:20*4) +
#     theme_bw()
###箱形图
realY<-c(as.vector(age1),as.vector(age2),as.vector(age3))
realx<-rep(c(rep("ourGeneset",11),rep("allGeneset",11)),3)
groupDe<-c(rep("adolescences",22),rep("adult",22),rep("latter adult ",22))
Data<-data.frame(Genegroup=realx,Gene_expresison=realY,Agegroup=groupDe)
p<-ggplot(data=Data, aes(x=Agegroup,y=Gene_expresison))+geom_boxplot(aes(fill=Genegroup))
p 
##
graph2ppt(file='expressionDiff1.pptx')




