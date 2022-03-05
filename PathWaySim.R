###################333
#
# categories of drugs comes from DrugBank
# CTD_chemicals and CTD_chem_pathways_enriched are coming from http://ctdbase.org/downloads/#cg
# Combined chemical and pathways and drugbank for IDs and ChemicalIDs 
# Applied IDF for weight of pathways 
# Finally applied tacomoto similarity
# the code takes     hours to run 3pm started -- 

library(dplyr)
setwd("/Users/idemirsoy/1DriveMyFsu/OneDrive - Florida State Students/925/MyoKardiaDownloads/PvSummerCodes/SimilarityMatrices")
#dru1Inquote="DB00321"

#Comb6=read.csv("DB00321/Comb6321.csv",header=T,sep=",")%>%dplyr::select(.,-X)

#load("dall.RData")
categ=dall$categories_drug
categ1=dplyr::filter(categ,`mesh-id`!="")

dbid = dplyr::select(categ1,-category)
names(dbid)=c("IDs","parent_key")
chem = read.csv('CTD_chemicals.csv.gz',comment.char = "#",header = F)
path = read.csv('CTD_chem_pathways_enriched.csv.gz', comment.char = "#",header=F)

chem1 = chem[, c(2, 5)]

path1 = path[, c(2, 5)]

names(path1)=c("ChemicalID", "path")
names(chem1)=c("ChemicalID", "IDs")

chem2=chem1 %>%
  mutate(IDs = gsub('MESH:', '', IDs)) %>% 
  tidyr::separate_rows(IDs, sep = '\\|')

chem3=chem2 %>%
  mutate(ChemicalID = gsub('MESH:', '', ChemicalID)) 


LL42=inner_join(chem3,path1,by="ChemicalID")
dbch1=inner_join(dbid,LL42,by="IDs")


dim(dbch1)
###

dbankT=dplyr::select(dbch1,c(path,parent_key))


TheOnes=Comb6$drug_2
TheOnes1=c(TheOnes,dru1Inquote)


f1=NULL
for(i in 1:length(TheOnes1))
{
  f1[[i]]=dplyr::filter(dbankT,parent_key==TheOnes1[i])
  print(paste0("i = ", i))
}

f33=do.call(rbind, f1)
f44=distinct(f33,path,parent_key)   ### drop the dublicated inside drugs 

PathMed=as.data.frame(table(f44$path))
PathMed$IDF=log(sum(PathMed$Freq)/PathMed$Freq)
names(PathMed)[names(PathMed)=="Var1"]=c("path")

Ppath1=left_join(f44,PathMed,by="path")
Ppath2=dplyr::select(Ppath1,-Freq)

Ppath2=as.data.frame(Ppath2)

pathSim=function(DrugNameInQoute,data1){
  M1=matrix(0,length(unique(data1$parent_key)),8)
  M1[,1]=DrugNameInQoute
  M1[,2]=dim(data1[(data1$parent_key)==DrugNameInQoute,])[1]
  M1[,7]=sum(data1[(data1$parent_key)==DrugNameInQoute,3])                  ## sum of terms in chosen 
  for(i in 1:length(unique(data1$parent_key))){
    # M1[i,1]=DrugNameInQoute
    # M1[i,2]=dim(data1[(data1$parent_key)==DrugNameInQoute,])[1]
    M1[i,3]=unique(data1$parent_key)[i]
    M1[i,4]=dim(data1[data1$parent_key==unique(data1$parent_key)[i],])[1]
    M1[i,5]=sum(data1[(data1$parent_key)==DrugNameInQoute,1] %in% data1[data1$parent_key==unique(data1$parent_key)[i],1])  # number of intersecting terms 
    M1[i,6]=sum(data1[(data1$parent_key)==DrugNameInQoute,][((data1[(data1$parent_key)==DrugNameInQoute,1] %in% data1[data1$parent_key==unique(data1$parent_key)[i],1])=="TRUE"),3]) #intersecting terms total
    # M1[i,7]=sum(data1[(data1$parent_key)==DrugNameInQoute,3])                  ## sum of terms in chosen 
    M1[i,8]= sum(data1[data1$parent_key==unique(data1$parent_key)[i],3])  ## sum of second drug
    print(paste0("i = ", i))
  }
  return(M1)
}


Path11=pathSim(dru1Inquote,Ppath2)  ## 
Path11=as.data.frame(Path11)
Path11$V2=as.numeric(Path11$V2)
Path11$V4=as.numeric(Path11$V4)
Path11$V5=as.numeric(Path11$V5)
Path11$V6=as.numeric(Path11$V6)
Path11$V7=as.numeric(Path11$V7)
Path11$V8=as.numeric(Path11$V8)

Path11$PathtcIDF=Path11$V6/(Path11$V7+Path11$V8-Path11$V6)  

names(Path11)[names(Path11)=="V1"]=c("drug_1")
names(Path11)[names(Path11)=="V3"]=c("drug_2")

PathSel=dplyr::select(Path11,c(drug_1,drug_2,PathtcIDF))


Comb7=inner_join(PathSel,Comb6,by="drug_2")%>%dplyr::select(.,-drug_1.y)

#write.csv(Comb7,"DB00321/Amitriptyline7.csv")
#Gene196=read.csv("DB196/geneSimDB196.csv",header=T,sep=",")%>%dplyr::select(.,-X)
#DB00196=inner_join(Gene196,Comb7,by="drug_2")%>%select(.,-c(drug_1.y,MolSim.y,ATCSim))
#write.csv(DB00196,"DB196/fluconazole.csv")

