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

