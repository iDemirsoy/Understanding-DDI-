#### Uses OFFSIDES to get meddra similarity 

setwd("/Users/idemirsoy/1DriveMyFsu/OneDrive - Florida State Students/925/MyoKardiaDownloads/PvSummerCodes/SimilarityMatrices")


#  TwoSide=read.csv("TWOSIDES.csv.gz",header=T,sep=",")
#  OneSide=read.csv("OFFSIDES.csv.gz",header=T,sep=",")
#  save.image("TwoSide.RData")

#load("dall.RData")
#dru1Inquote="DB01242"

load("TwoSide.RData")

OneS1=dplyr::select(OneSide,drug_rxnorn_id,condition_meddra_id)%>%data.frame()
names(OneS1)=c("identifier","MedDRA")



exT=dall$external_identifiers_drug
exTrec=dplyr::filter(exT,resource=="RxCUI")%>%data.frame()


OneSideDB=inner_join(OneS1,exTrec,by="identifier")
OneSiDB=dplyr::select(OneSideDB,MedDRA,parent_key)

table(OneSiDB$MedDRA)

DclasssDF=data.frame(table(OneSiDB$MedDRA))
DclasssDF$IDF=log(sum(DclasssDF$Freq)/DclasssDF$Freq)
names(DclasssDF)[names(DclasssDF)=="Var1"]=c("MedDRA")
DrugClass1=left_join(OneSiDB,DclasssDF,by="MedDRA")
DrugClass12=dplyr::select(DrugClass1,c(MedDRA,parent_key,IDF))

DrugClass12=data.frame(DrugClass12)   ## problem was theos


pathSim=function(DrugNameInQoute,data1){
  M1=matrix(0,length(unique(data1$parent_key)),8)
  M1[,1]=DrugNameInQoute
  M1[,2]=dim(data1[(data1$parent_key)==DrugNameInQoute,])[1]
  M1[,7]=sum(data1[(data1$parent_key)==DrugNameInQoute,3])                  ## sum of terms in chosen 
  for (i in 1:length(unique(data1$parent_key))){
    #M1[i,1]=DrugNameInQoute
    #M1[i,2]=dim(data1[(data1$parent_key)==DrugNameInQoute,])[1]
    M1[i,3]=unique(data1$parent_key)[i]
    M1[i,4]=dim(data1[data1$parent_key==unique(data1$parent_key)[i],])[1]
    M1[i,5]=sum(data1[(data1$parent_key)==DrugNameInQoute,1] %in% data1[data1$parent_key==unique(data1$parent_key)[i],1])  # number of intersecting terms 
    M1[i,6]=sum(data1[(data1$parent_key)==DrugNameInQoute,][((data1[(data1$parent_key)==DrugNameInQoute,1] %in% data1[data1$parent_key==unique(data1$parent_key)[i],1])=="TRUE"),3]) #intersecting terms total
    #M1[i,7]=sum(data1[(data1$parent_key)==DrugNameInQoute,3])                  ## sum of terms in chosen 
    M1[i,8]= sum(data1[data1$parent_key==unique(data1$parent_key)[i],3])  ## sum of second drug
    print(paste0("i = ", i))
  }
  return(M1)
}


MedD1=pathSim(dru1Inquote,DrugClass12)  ## 
MedD1=as.data.frame(MedD1)
MedD1$V2=as.numeric(MedD1$V2)
MedD1$V4=as.numeric(MedD1$V4)
MedD1$V5=as.numeric(MedD1$V5)
MedD1$V6=as.numeric(MedD1$V6)
MedD1$V7=as.numeric(MedD1$V7)
MedD1$V8=as.numeric(MedD1$V8)


MedD1$CosineIDF=MedD1$V6/sqrt(MedD1$V7*MedD1$V8)   ### Cosine similarity 
MedD1$tcIDF=MedD1$V6/(MedD1$V7+MedD1$V8-MedD1$V6)    ### Tanimoto similarity 

MedD1[order(MedD1$CosineIDF,decreasing = T),]


MedD2=dplyr::select(MedD1,c(V1,V3,tcIDF))
names(MedD2)=c("drug_1","drug_2","MedSim")

# fdat=read.csv("FinalData2.csv",header=T,sep=",")
# fdat1=dplyr::filter(fdat, !is.na(ATCCosineIDF))
# fdat2=dplyr::select(fdat1,-c(drug_1.y))

# EE12=inner_join(fdat2,E12,by="drug_2")

# write.csv(E12,"MedDRAsim.csv")


Comb3=inner_join(MedD2,Comb2,by="drug_2")%>%dplyr::select(.,-drug_1.y)








