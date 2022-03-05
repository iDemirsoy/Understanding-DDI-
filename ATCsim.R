setwd("/Users/idemirsoy/1DriveMyFsu/OneDrive - Florida State Students/925/MyoKardiaDownloads/PvSummerCodes/SimilarityMatrices")

library(ggplot2)
library(dplyr)

#load(file = "dbpars.RData")
#dru1Inquote="DB00176"
load("dall.RData")

Datc=dall$atc_codes_drug

Datc[Datc$`drugbank-id`==dru1Inquote,]
Datc[Datc$code_4=="B",]

length(unique(Datc$code_1))

df1=as.data.frame(table(Datc$code_4))
df1$Var1=as.character(df1$Var1)
ggplot(df1,aes(x=Var1,y=Freq ,fill=Var1))+
  geom_bar(stat = "identity")+  geom_text(aes(label = Freq), vjust = -0.3)+
  xlab("ATC groups") +
  ylab("Frequency")


data.frame(Datc[Datc$code_2=="B01A",])[,c(1,3,5,7,9,10)]


DatcDF1=data.frame(table(Datc$code_1))
DatcDF1$IDFc1=log(sum(DatcDF1$Freq)/DatcDF1$Freq)
names(DatcDF1)[names(DatcDF1)=="Var1"]=c("code_1")
names(DatcDF1)[names(DatcDF1)=="Freq"]=c("FreqC1")


DatcDF2=data.frame(table(Datc$code_2))
DatcDF2$IDFc2=log(sum(DatcDF2$Freq)/DatcDF2$Freq)
names(DatcDF2)[names(DatcDF2)=="Var1"]=c("code_2")
names(DatcDF2)[names(DatcDF2)=="Freq"]=c("FreqC2")


DatcDF3=data.frame(table(Datc$code_3))
DatcDF3$IDFc3=log(sum(DatcDF3$Freq)/DatcDF3$Freq)
names(DatcDF3)[names(DatcDF3)=="Var1"]=c("code_3")
names(DatcDF3)[names(DatcDF3)=="Freq"]=c("FreqC3")


DatcDF4=data.frame(table(Datc$code_4))
DatcDF4$IDFc4=log(sum(DatcDF4$Freq)/DatcDF4$Freq)
names(DatcDF4)[names(DatcDF4)=="Var1"]=c("code_4")
names(DatcDF4)[names(DatcDF4)=="Freq"]=c("FreqC4")


DrugA=dplyr::select(Datc,-c(level_1,level_2,level_3,level_4))
DatcDF1$code_1=as.character(DatcDF1$code_1)
DatcDF2$code_2=as.character(DatcDF2$code_2)
DatcDF3$code_3=as.character(DatcDF3$code_3)
DatcDF4$code_4=as.character(DatcDF4$code_4)


DrugATC1=left_join(DrugA,DatcDF1,by="code_1")
DrugATC2=left_join(DrugATC1,DatcDF2,by="code_2")
DrugATC3=left_join(DrugATC2,DatcDF3,by="code_3")
DrugATC4=left_join(DrugATC3,DatcDF4,by="code_4")



#DrugAcomb=dplyr::select(DrugATC4,c(codecode_2,`drugbank-id`,IDFc1,IDFc2,IDFc3,IDFc4))
DrugAcomb=data.frame(DrugATC4)   ## problem was theos

names(DrugAcomb)[names(DrugAcomb)=="drugbank.id"]=c("parent_key")
splitted <- t(sapply(DrugAcomb$atc_code, function(x) substring(x, first=c(1,1,1,1), last=c(1,3,4,5))))
DBB=dplyr::select(DrugAcomb,c(code_4,code_3,code_2,code_1,IDFc1,IDFc2,IDFc3,IDFc4,parent_key))
DBB=as.data.frame(DBB)


CA1=combn(unique(DBB$parent_key),2)
Ata=c(which(CA1[2,]==dru1Inquote),which(CA1[1,]==dru1Inquote))
CA2=CA1[,Ata]


A=matrix(0,dim(MolSim1)[1],1)
A1=matrix(0,dim(MolSim1)[1],1)
for(i in 1:dim(MolSim1)[1]){
  if (length(which(CA2==as.character(MolSim1$drug_2)[i]))!=0)
  { 
    if(length(which(CA2[1,]==as.character(MolSim1$drug_2)[i]))!=0)
    {
      A[i]=which(CA2[1,]==as.character(MolSim1$drug_2)[i])
      A1[i]=as.character(MolSim1$drug_2)[i]
    } else {
      A[i]=which(CA2[2,]==as.character(MolSim1$drug_2)[i])
      A1[i]=as.character(MolSim1$drug_2)[i]
    }
  } else {
    A[i]=0
    A1[i]=as.character(MolSim1$drug_2)[i]
  }
}

Abb=cbind(A,A1)


L=NULL
for(k in 1:length(A)){
  if (A[k]==0) {
    L[[k]]=0
  } else {
id1 = dplyr::filter(DBB,parent_key==CA2[,A[k]][1])
id2 = dplyr::filter(DBB,parent_key==CA2[,A[k]][2])
mat = matrix(0L, nrow = dim(id1)[1], ncol = dim(id2)[1])
for(i in 1:dim(id1)[1]){
  for(j in 1:dim(id2)[1]){
    if (id1[i,4]==id2[j,4]) {
      mat[i, j] = 1
    } else if (id1[i,3]==id2[j,3]) {
    c=sum(id1[i,c(6:8)])   
    a=sum(id1[i,c(5:8)])                  ## sum of terms in chosen 
    b=sum(id2[j,c(5:8)])   
    mat[i, j] = c/sqrt(a*b)
    } else if (id1[i,2]==id2[j,2]) {
      c=sum(id1[i,c(7:8)])   
      a=sum(id1[i,c(5:8)])                  ## sum of terms in chosen 
      b=sum(id2[j,c(5:8)])   
      mat[i, j] = c/sqrt(a*b)
    }  else if (id1[i,1]==id2[j,1]) {
      c=sum(id1[i,c(8:8)])   
      a=sum(id1[i,c(5:8)])                  ## sum of terms in chosen 
      b=sum(id2[j,c(5:8)])   
      mat[i, j] = c/sqrt(a*b)
    } else {
      mat[i, j] = 0
    }
  }
}
L[[k]]=mat
print(paste0("k = ",k))
}
}

ATCss = sapply(L, max)

ATCAB=cbind(A,A1,ATCss)%>%data.frame()
names(ATCAB)=c("X1","drug_2","ATCsim")

Comb2=left_join(MolSim1,ATCAB,by="drug_2")%>%filter(X1!=0)%>%dplyr::select(.,-c(X1))




################## Ignore   ###############


# ATCsim[order(ATCsim$ATCCosineIDF,decreasing =T),]
# dim(ATCsim%>%filter(ATCCosineIDF>0))

# FinalData2=left_join(FinalData,ATCsim,by="drug_2")

# FinalData2[FinalData2$Y==0,]%>%filter(is.na(ATCCosineIDF))
## found them online to update 
# FinalData2[FinalData2$drug_2=="DB01583",]$ATCCosineIDF=0
# FinalData2[FinalData2$drug_2=="DB00974",]$ATCCosineIDF=0
# FinalData2[FinalData2$drug_2=="DB01013",]$ATCCosineIDF=0
# FinalData2[FinalData2$drug_2=="DB03894",]$ATCCosineIDF=0
# FinalData2[FinalData2$drug_2=="DB08804",]$ATCCosineIDF=0
# FinalData2[FinalData2$drug_2=="DB02379",]$ATCCosineIDF=0.2784326 ## got from DB00145 they have same code_1







