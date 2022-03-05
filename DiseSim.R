setwd("/Users/idemirsoy/1DriveMyFsu/OneDrive - Florida State Students/925/MyoKardiaDownloads/PvSummerCodes/SimilarityMatrices")

library(dbparser)
library(dplyr)

#load("dall.RData")
#dru1Inquote="DB01242"

categ=dall$categories_drug
categ1=dplyr::filter(categ,`mesh-id`!="")
dbid = dplyr::select(categ1,-category)
names(dbid)=c("IDs","parent_key")
chem = read.csv('CTD_chemicals.csv.gz', comment.char = "#",header=F)
dise = read.csv('CTD_chemicals_diseases.csv.gz', comment.char = "#",header=F)
omim = read.table('MimMiner_Exp_AC_T_TXCS_basedonACMESH_filt_RW.mat.gz', 
                  sep = '\t', row.names = 1, header = FALSE)
colnames(omim) = row.names(omim)

chem1 = chem[, c(2, 5)]
names(chem1)=c("ChemicalID", "IDs")

chem2=chem1 %>%
  mutate(IDs = gsub('MESH:', '', IDs)) %>%
  tidyr::separate_rows(IDs, sep = '\\|')%>%
  mutate(ChemicalID = gsub('MESH:', '', ChemicalID)) %>%data.frame()



#chem11 = chem1[which(chem1[, 2] != ''), ]
#chem11[, 1] = gsub('MESH:', '', chem11[, 1])
#names(chem11)=c("ChemicalID", "IDs")

# expand duplicated DrugBank IDs
#chemtmp = vector('list', nrow(chem11))
#for (i in 1:nrow(chem11)) chemtmp[[i]] = as.character(chem11[i, 2])
#names(chemtmp) = as.character(chem11[, 1])
#chemtmp = sapply(chemtmp, strsplit, split = '\\|')
#tmp = unlist(chemtmp)
#chem2 = data.frame(IDs = tmp, ChemicalID = names(tmp))
#chem2[, 1] = gsub('MESH:', '', chem2[, 1])

dise101 = dise[, c(2, 9)]
dise1 = dise101[which(dise101[, 2] != ''), ]
names(dise1)=c("ChemicalID", "DisesOMIM")


dise2=dise1 %>%
  tidyr::separate_rows(DisesOMIM, sep = '\\|') %>%data.frame()


# expand duplicated DrugBank IDs
#disetmp = vector('list', nrow(dise1))
#for (i in 1:nrow(dise1)) disetmp[[i]] = as.character(dise1[i, 2])
#names(disetmp) = as.character(dise1[, 1])
#disetmp = sapply(disetmp, strsplit, split = '\\|')
#dtmp = unlist(disetmp)
#dise2 = data.frame(DisesOMIM = dtmp, ChemicalID = names(dtmp))


String1=function(data){
  B1=data
  B2=NULL
  for (i in 1:length(B1)) {if(length(B1[[i]]==7))
  {
    B2[[i]]=paste0(B1[[i]][1:7],collapse = "")
  } else{
    B2[[i]]=paste0(B1[[i]][1:7],collapse = "")
  }
  }
  B2=unlist(B2)
  return(B2)
}

dise2$ChemicalID=as.character(dise2$ChemicalID)

DisChe=inner_join(chem2,dise2,by="ChemicalID")

DCdb1=inner_join(dbid,DisChe,by="IDs")

DiseBank=dplyr::select(DCdb1,parent_key,DisesOMIM)


DiseMat=as.data.frame(table(DiseBank$DisesOMIM))
DiseMat$IDF=log(sum(DiseMat$Freq)/DiseMat$Freq)
names(DiseMat)[names(DiseMat)=="Var1"]=c("DisesOMIM")

DiseMat$DisesOMIM=as.character(DiseMat$DisesOMIM)

DiseP1=left_join(DiseBank,DiseMat,by="DisesOMIM")
DiseP2=dplyr::select(DiseP1,c(DisesOMIM,parent_key,IDF))

DiseP2=as.data.frame(DiseP2)


DRUG2=c(Comb3$drug_2,dru1Inquote)   ### to decrease the time of the code, we subset the drugs 

DDi=NULL
for(i in 1:length(DRUG2)){ DDi[[i]]=which(DiseP2$parent_key==DRUG2[i])}
names(DDi) = as.character(DRUG2)
tmpDise = unlist(DDi)
DiseDF = data.frame(bi1 = names(tmpDise),bi2 = tmpDise)
DiseDF$bi1=as.character(DiseDF$bi1)
DiseDF$bi2=as.character(DiseDF$bi2)

DB1=strsplit(DiseDF$bi1,split = '')
DiseDF$bi1=String1(DB1)

New12=DiseP2[DiseDF$bi2,]
NewMat=distinct(New12,DisesOMIM, parent_key,IDF)


pathSim=function(DrugNameInQoute,data1){
  M1=matrix(0,length(unique(data1$parent_key)),8)
  M1[,1]=DrugNameInQoute
  M1[,2]=dim(data1[(data1$parent_key)==DrugNameInQoute,])[1]
  M1[,7]=sum(data1[(data1$parent_key)==DrugNameInQoute,3])                       ## sum of terms in chosen 
  for(i in 1:length(unique(data1$parent_key))){
    # M1[i,1]=DrugNameInQoute
    # M1[i,2]=dim(data1[(data1$parent_key)==DrugNameInQoute,])[1]
    M1[i,3]=unique(data1$parent_key)[i]
    M1[i,4]=dim(data1[data1$parent_key==unique(data1$parent_key)[i],])[1]
    M1[i,5]=sum(data1[(data1$parent_key)==DrugNameInQoute,1] %in% data1[data1$parent_key==unique(data1$parent_key)[i],1])  # number of intersecting terms 
    M1[i,6]=sum(data1[(data1$parent_key)==DrugNameInQoute,][((data1[(data1$parent_key)==DrugNameInQoute,1] %in% data1[data1$parent_key==unique(data1$parent_key)[i],1])=="TRUE"),3]) #intersecting terms total
    # M1[i,7]=sum(data1[(data1$parent_key)==DrugNameInQoute,3])                  ## sum of terms in chosen 
    M1[i,8]= sum(data1[data1$parent_key==unique(data1$parent_key)[i],3])         ## sum of second drug
    print(paste0("i = ", i))
  }
  return(M1)
}

Disease=pathSim(dru1Inquote,NewMat)  ## 
Disease=as.data.frame(Disease)
Disease$V2=as.numeric(Disease$V2)
Disease$V4=as.numeric(Disease$V4)
Disease$V5=as.numeric(Disease$V5)
Disease$V6=as.numeric(Disease$V6)
Disease$V7=as.numeric(Disease$V7)
Disease$V8=as.numeric(Disease$V8)


Disease$DiseIDF=Disease$V6/(Disease$V7+Disease$V8-Disease$V6)                    ### tanimoto similarity 
DiseaseDF=dplyr::select(Disease,V1,V3,DiseIDF)
names(DiseaseDF)=c("drug_1","drug_2","DiseaseSim")

Comb4=inner_join(Comb3,DiseaseDF,by="drug_2")%>%dplyr::select(.,-drug_1.x)

#write.csv(Comb4,"1242Comb4.csv")

