#*******************************************************************************#
#R version 4.0.2 (2020-06-22) -- "Taking Off Again"
#Copyright (C) 2020 The R Foundation for Statistical Computing
#Platform: x86_64-apple-darwin17.0 (64-bit)


#*******************************************************************************#
#Load the R Packages
library("matrixStats")
#‘0.58.0’

library("data.table")
#‘1.14.3’


#*******************************************************************************#
#Setup the work path, Your data path
setwd(YourDir)
rm(list=ls())


#*******************************************************************************#
#Initialize variable
#Tissue name
TissTitle=c("eWAT","sWAT","Liver","Muscle")

#Comparison list
comg=c("Chow_FastvsChow_AdLibitum","Chow_RefeedvsChow_Fast",
       "HFD_FastvsHFD_AdLibitum","HFD_RefeedvsHFD_Fast")


#*******************************************************************************#
#Input the expression dynamics and statistical significance information
for(tissi in TissTitle){
  for(tcomgi in comg){
    tFCFDR=fread(paste0(tissi,"/",tcomgi,".log2FC.FDR.txt"),sep = "\t",data.table = F)
    rownames(tFCFDR)=tFCFDR$Geneid
    assign(paste0(tissi,".",tcomgi,".FCFDR"),tFCFDR)
  }
}


#*******************************************************************************#
#Calculate the Gene Elasticity Score (GElaS)
for(tissi in TissTitle){
  #Calculate the GElaS in CD and HFD mice
  for(tagei in c("Chow","HFD")){
    #Variable name for Comparison
    tFN.name=paste0(tissi,".",tagei,"_Fast","vs",tagei,"_AdLibitum.FCFDR")
    tRF.name=paste0(tissi,".",tagei,"_Refeed","vs",tagei,"_Fast.FCFDR")
    
    #Create the data.frame to store the GElaS
    tScoreMat=get(tFN.name)[,c("Geneid", "Symbol")]
    #Gene List
    tgid=tScoreMat$Geneid
    
    #Expression dynamics in FastvsAdLibitum and RefeedvsFast
    tFN.FC = get(tFN.name)[tgid,"log2FC"]
    tRF.FC = get(tRF.name)[tgid,"log2FC"]
    
    #Sign of expression dynamics in the AdLibitum-Fast-Refeed cycle
    tMultilog2FCSign = -1*sign(tFN.FC * tRF.FC)
    tMultilog2FCSign[tMultilog2FCSign <= 0] = 0
    
    #Absolute expression dynamics matrix
    tAbslog2FC.Mat=data.frame(Abslog2FC.FN=abs(get(tFN.name)[tgid,"log2FC"]),
                              Abslog2FC.RF=abs(get(tRF.name)[tgid,"log2FC"]),
                              stringsAsFactors = F)
    
    #Calculate restoration extent of gene expression in the AdLibitum-Fast-Refeed cycle
    tlog2FC.Ratio=rowMins(as.matrix(tAbslog2FC.Mat))/rowMaxs(as.matrix(tAbslog2FC.Mat))
    tlog2FC.Ratio[is.na(tlog2FC.Ratio)]=0
    
    #Calculate the weight for FDR (statistical significance)
    tFDR.FN=-log10(get(tFN.name)[tgid,"adj.P.Val"])
    tFDR.RF=-log10(get(tRF.name)[tgid,"adj.P.Val"])
    tFDR.FN[tFDR.FN <= -log10(0.05)]=tFDR.FN[tFDR.FN <= -log10(0.05)]/(-log10(0.05))
    tFDR.FN[tFDR.FN > -log10(0.05)]=1
    tFDR.RF[tFDR.RF <= -log10(0.05)]=tFDR.RF[tFDR.RF <= -log10(0.05)]/(-log10(0.05))
    tFDR.RF[tFDR.RF > -log10(0.05)]=1
    
    #Integrate the expression dynamics, statistical significance, and restoration extent to 
    #calculate the GElaS
    tAbslog2FC.FDRSum=(tAbslog2FC.Mat$Abslog2FC.FN*tFDR.FN+tAbslog2FC.Mat$Abslog2FC.RF*tFDR.RF)
    tGElaS= tMultilog2FCSign*tlog2FC.Ratio*tAbslog2FC.FDRSum
    tScoreMat[,paste0(tagei,"_GElaS")]=data.frame(tGElaS,stringsAsFactors = F)
    
    #Output the GElaS
    fwrite(tScoreMat,file = paste0(tissi,"/",tagei,".GElaS.txt"),sep = "\t")
  }
}

#*******************************************************************************#


#gene Type
'gtype=read.csv(file = "~/Postdoc/pub/GeneCode/vM26/gencode.vM26.annotation.gtype.txt",
               sep="\t",header = T,stringsAsFactors = F)
LNCgtype=read.csv(file = "~/Postdoc/pub/GeneCode/vM26/gencode.vM26.long_noncoding_RNAs.gtype.txt",
                  sep="\t",header = T,stringsAsFactors = F)
PCgtype=subset(gtype,Type=="protein_coding")'

rm(list = ls())
#diff plasticity score
lncTypeNum=data.frame(matrix(nrow = 0,ncol=4))
colnames(lncTypeNum)=c("ToptalGene","TotalUpDw","HFDvsCD.GElaSUp",
                       "HFDvsCD.GElaS")
NoDiffTypeNum=data.frame(matrix(nrow = 0,ncol=2))
colnames(NoDiffTypeNum)=c("ToptalGene","TotalNoDiff")

#We used PSFCCut=2 and PSDFCut=0.5
ScoFCCut=2 #1.5
ScoDFCut=0.5 #0.45
ScoPlus=0.1 #0.2
ScoCut=0.2 #
ScoFCNo=0.25 #for the persistent gene
ScoDFNo=0.1 #for the persistent gene 

for(tissi in c("eWAT","sWAT","Liver","Muscle")){
  tTissYGElaS=fread(file = paste0(tissi,".CD.GElaS.txt"),data.table = F)
  rownames(tTissYGElaS)=tTissYGElaS$Geneid
  tTissAGElaS=fread(file = paste0(tissi,".Aged.GElaS.txt"),data.table = F)
  rownames(tTissAGElaS)=tTissAGElaS$Geneid
  tTissGElaS=cbind(tTissYGElaS,subset(tTissAGElaS[tTissYGElaS$Geneid,],select = -c(1:2)))
  
  #tTissGElaS$AgedvsYoung.Plus0.1log2FC.PS=log2((tTissGElaS$Aged.log2FC.PlaSco+0.1)/(tTissGElaS$Young.log2FC.PlaSco+0.1))
  tFCScoName=paste0("AgedvsYoung.log2FC.GElaS")
  tDiffScoName=paste0("AgedvsYoung.Diff.GElaS")
  
  tTissGElaS[,tFCScoName] = log2((tTissGElaS$Aged_GElaS+ScoPlus)/(tTissGElaS$Young_GElaS+ScoPlus))
  tTissGElaS[,tDiffScoName] = tTissGElaS$Aged_GElaS - tTissGElaS$Young_GElaS
  tTissGElaSFilt = subset(tTissGElaS,Aged_GElaS > ScoCut | Young_GElaS > ScoCut)
  
  'write.table(tTissGElaS,file = paste0(tissi,".AgedvsYoung.TPScore.xls"),sep = "\t",quote = F,row.names = F)
  write.table(tTissGElaSFilt,file = paste0(tissi,".AgedvsYoung.TPS.Filt",ScoCut,".xls"),
              sep = "\t",quote = F,row.names = F)
  write.table(subset(tTissGElaSFilt,Type=="protein_coding"),file = paste0(tissi,".AgedvsYoung.TPS.Filt",PSCut,".PC.xls"),
              sep = "\t",quote = F,row.names = F)'
  
  #log2FC
  tTissGElaSFiltFCUpDw = tTissGElaSFilt[abs(tTissGElaSFilt[[tFCScoName]]) > ScoFCCut,]#1.5 or 2
  tTissGElaSFiltFCUp = tTissGElaSFilt[tTissGElaSFilt[[tFCScoName]] > ScoFCCut,]
  tTissGElaSFiltFCDw = tTissGElaSFilt[tTissGElaSFilt[[tFCScoName]] < -ScoFCCut,]#
  
  tTissGElaSFiltFCNo = subset(tTissGElaSFilt[abs(tTissGElaSFilt[[tFCScoName]]) < ScoFCNo,],Young_GElaS > ScoCut & Aged_GElaS > ScoCut)#No difference
  print(paste(tissi," log2FC,AvsY.Up:",nrow(tTissGElaSFiltFCUp),"AvsY.Down:",nrow(tTissGElaSFiltFCDw),
              "No Diff:",nrow(tTissGElaSFiltFCNo)))

  lncTypeNum[paste0(tissi,".log2FC",ScoDFCut),]=c(nrow(tTissGElaSFilt),nrow(tTissGElaSFiltFCUpDw),
                                                nrow(tTissGElaSFiltFCUp),nrow(tTissGElaSFiltFCDw))
  NoDiffTypeNum[paste0(tissi,".log2FC",ScoFCNo),]=c(nrow(tTissGElaSFilt),nrow(tTissGElaSFiltFCNo))
  
  'write.table(subset(tTissGElaSFilt,Type=="protein_coding")[,c("Symbol",tFCScoName)],
              file = paste0(tissi,".AgedvsYoung.TPS.Filt",ScoCut,".Plus",ScoPlus,"log2FC.Pro.rnk"),
              sep = "\t",quote = F,row.names = F)
  write.table(subset(tTissGElaSFiltFCUpDw),file = paste0(tissi,".AgedvsYoung.TPS.Filt",ScoCut,".Plus",ScoPlus,"log2FC",ScoFCCut,"UpDw.xls"),
              sep = "\t",quote = F,row.names = F)
  write.table(subset(tTissGElaSFiltFCUp,Type=="protein_coding"),file = paste0(tissi,".AgedvsYoung.TPS.Filt",ScoCut,".Plus",ScoPlus,"log2FC",ScoFCCut,"Up.PC.xls"),
              sep = "\t",quote = F,row.names = F)
  write.table(subset(tTissGElaSFiltFCDw,Type=="protein_coding"),file = paste0(tissi,".AgedvsYoung.TPS.Filt",ScoCut,".Plus",ScoPlus,"log2FC",ScoFCCut,"Dw.PC.xls"),
              sep = "\t",quote = F,row.names = F)
  write.table(tTissGElaSFiltFCUpLnc,file = paste0(tissi,".AgedvsYoung.TPS.Filt",ScoCut,".Plus",ScoPlus,"log2FC",ScoFCCut,"Up.LNC.xls"),
              sep = "\t",quote = F,row.names = F)
  write.table(tTissGElaSFiltFCDwLnc,file = paste0(tissi,".AgedvsYoung.TPS.Filt",ScoCut,".Plus",ScoPlus,"log2FC",ScoFCCut,"Dw.LNC.xls"),
              sep = "\t",quote = F,row.names = F)'
  
  #Top 500 and Bottom 500
  'write.table(tTissGElaSFiltFCUpTop500PC,file = paste0(tissi,".AgedvsYoung.TPS.Filt",ScoCut,".Plus",ScoPlus,"log2FC",ScoFCCut,"UpTop500.PC.xls"),
              sep = "\t",quote = F,row.names = F)
  write.table(tTissGElaSFiltFCDwTop500PC,file = paste0(tissi,".AgedvsYoung.TPS.Filt",ScoCut,".Plus",ScoPlus,"log2FC",ScoFCCut,"DwTop500.PC.xls"),
              sep = "\t",quote = F,row.names = F)
  write.table(tTissGElaSFiltFCNCBot500PC,file = paste0(tissi,".AgedvsYoung.TPS.Filt",ScoCut,".Plus",ScoPlus,"log2FC",ScoFCCut,"NCBot500.PC.xls"),
              sep = "\t",quote = F,row.names = F)'
  
  #No difference genes
  'write.table(tTissGElaSFiltFCNo,file = paste0(tissi,".AgedvsYoung.TPS.Filt",ScoCut,".Plus",ScoPlus,"log2FC",ScoFCNo,"NoDiff.xls"),
              sep = "\t",quote = F,row.names = F)
  write.table(subset(tTissGElaSFiltFCNo,Type=="protein_coding"),file = paste0(tissi,".AgedvsYoung.TPS.Filt",ScoCut,".Plus",ScoPlus,"log2FC",ScoFCNo,"NoDiff.PC.xls"),
              sep = "\t",quote = F,row.names = F)
  write.table(tTissGElaSFiltFCNoLnc,file = paste0(tissi,".AgedvsYoung.TPS.Filt",ScoCut,".Plus",ScoPlus,"log2FC",ScoFCNo,"NoDiff.LNC.xls"),
              sep = "\t",quote = F,row.names = F)'
  
  
  
  #Diff
  tTissGElaSFiltDiffUpDw = tTissGElaSFilt[abs(tTissGElaSFilt[[tDiffScoName]]) > ScoDFCut,]#
  tTissGElaSFiltDiffUp = tTissGElaSFilt[tTissGElaSFilt[[tDiffScoName]] > ScoDFCut,]
  tTissGElaSFiltDiffDw = tTissGElaSFilt[tTissGElaSFilt[[tDiffScoName]] < -ScoDFCut,]#
  
  tTissGElaSFiltDiffNo = subset(tTissGElaSFilt[abs(tTissGElaSFilt[[tDiffScoName]]) < ScoDFNo,],Young_GElaS > ScoCut & Aged_GElaS > ScoCut)#no differenc
  
  print(paste(tissi," Diff,AvsY.Up:",nrow(tTissGElaSFiltDiffUp),"AvsY.Down:",nrow(tTissGElaSFiltDiffDw),
              "No Diff:",nrow(tTissGElaSFiltDiffNo)))

  lncTypeNum[paste0(tissi,".Diff",ScoDFCut),]=c(nrow(tTissGElaSFilt),nrow(tTissGElaSFiltDiffUpDw),
                                         nrow(tTissGElaSFiltDiffUp),nrow(tTissGElaSFiltDiffDw))
  NoDiffTypeNum[paste0(tissi,".Diff",ScoDFNo),]=c(nrow(tTissGElaSFilt),nrow(tTissGElaSFiltDiffNo))
  
  'write.table(subset(tTissGElaSFilt,Type=="protein_coding")[,c("Symbol",tDiffScoName)],
              file = paste0(tissi,".AgedvsYoung.TPS.Filt",ScoCut,".Diff.Pro.rnk"),
              sep = "\t",quote = F,row.names = F)
  write.table(subset(tTissGElaSFiltDiffUpDw),file = paste0(tissi,".AgedvsYoung.TPS.Filt",ScoCut,".Diff",ScoDFCut,"UpDw.xls"),
              sep = "\t",quote = F,row.names = F)
  write.table(subset(tTissGElaSFiltDiffUp,Type=="protein_coding"),file = paste0(tissi,".AgedvsYoung.TPS.Filt",ScoCut,".Diff",ScoDFCut,"Up.PC.xls"),
              sep = "\t",quote = F,row.names = F)
  write.table(subset(tTissGElaSFiltDiffDw,Type=="protein_coding"),file = paste0(tissi,".AgedvsYoung.TPS.Filt",ScoCut,".Diff",ScoDFCut,"Dw.PC.xls"),
              sep = "\t",quote = F,row.names = F)
  write.table(tTissGElaSFiltDiffUpLnc,file = paste0(tissi,".AgedvsYoung.TPS.Filt",ScoCut,".Diff",ScoDFCut,"Up.LNC.xls"),
              sep = "\t",quote = F,row.names = F)
  write.table(tTissGElaSFiltDiffDwLnc,file = paste0(tissi,".AgedvsYoung.TPS.Filt",ScoCut,".Diff",ScoDFCut,"Dw.LNC.xls"),
              sep = "\t",quote = F,row.names = F)
  
  #Top 500 and Bottom 500
  write.table(tTissGElaSFiltDiffUpTop500PC,file = paste0(tissi,".AgedvsYoung.TPS.Filt",ScoCut,".Diff",ScoDFCut,"UpTop500.PC.xls"),
               sep = "\t",quote = F,row.names = F)
  write.table(tTissGElaSFiltDiffDwTop500PC,file = paste0(tissi,".AgedvsYoung.TPS.Filt",ScoCut,".Diff",ScoDFCut,"DwTop500.PC.xls"),
              sep = "\t",quote = F,row.names = F)
  write.table(tTissGElaSFiltDiffNCBot500PC,file = paste0(tissi,".AgedvsYoung.TPS.Filt",ScoCut,".Diff",ScoDFCut,"NCBot500.PC.xls"),
              sep = "\t",quote = F,row.names = F)'
  
  'write.table(tTissGElaSFiltDiffNo,file = paste0(tissi,".AgedvsYoung.TPS.Filt",ScoCut,".Diff",ScoDFNo,"NoDiff.xls"),
              sep = "\t",quote = F,row.names = F)
  write.table(subset(tTissGElaSFiltDiffNo,Type=="protein_coding"),file = paste0(tissi,".AgedvsYoung.TPS.Filt",ScoCut,".Diff",ScoDFNo,"NoDiff.PC.xls"),
              sep = "\t",quote = F,row.names = F)
  write.table(tTissGElaSFiltDiffNoLnc,file = paste0(tissi,".AgedvsYoung.TPS.Filt",ScoCut,".Diff",ScoDFNo,"NoDiff.LNC.xls"),
              sep = "\t",quote = F,row.names = F)'
  
  #plot
}

'write.table(lncTypeNum,file = paste0("Filt",ScoCut,".Plus",ScoPlus,"log2FC",ScoFCCut,".Diff",ScoDFCut,".DTScoNumber.txt"),
            sep = "\t",quote = F,col.names = NA)
write.table(NoDiffTypeNum,file = paste0("Filt",ScoCut,".Plus",ScoPlus,"log2FC",ScoFCNo,".Diff",ScoDFNo,".NoDiffNumber.txt"),
            sep = "\t",quote = F,col.names = NA)'

ls()
