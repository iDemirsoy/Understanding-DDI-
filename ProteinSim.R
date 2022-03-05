###############################################################################
## Dall is dragbank dataset 
## uniprot_sprot.fasta.gz is downloaded from uniprot.org - Reviewed(Swiss-Prot)
## 


#BiocManager::install("Biostrings")
#install.packages("seqinr")
setwd("/Users/idemirsoy/1DriveMyFsu/OneDrive - Florida State Students/925/MyoKardiaDownloads/PvSummerCodes/SimilarityMatrices")
#Comb4=read.csv("1242Comb4.csv",header=T,sep=",")%>%select(-X)
memory.limit(98732423532)
memory.size(98732423532)
#dru1Inquote="DB01242"
#load("dall.RData")

library(seqinr)
library(dplyr)
library(Biostrings)



#Comb3=read.csv("Comb4.csv",sep=",",header=T)
#Comb4=dplyr::select(Comb3,-X)

## this file has 2 rows 
# 1 row is name of drug_1
# 2 row is name of drug_2

## fatsa file 
fasta2 = read.fasta('uniprot_sprot.fasta.gz',seqtype = "AA")
fasta22=names(fasta2)
B1= strsplit(fasta22,split='')

## target gene information
PK3=dall$targ_drug[,c("id","parent_key")]
PI3=dall$polypeptide_target_drug[,c("id","parent_id")]
names(PI3)=c("GeneID","id")
PIK=left_join(PI3,PK3,by="id")

data(BLOSUM62)


String2=function(data,from,to)
{
  B2=NULL
  for(i in 1:length(data)){B2[[i]]=paste0(data[[i]][from:to],collapse = "")}
  B2=unlist(B2)
  return(B2)
}
GeneID=String2(B1,4,9)

fastaString=function(data){
  B22=NULL
  for(i in 1:length(data)){
    d=length(data[[i]]);
    B22[[i]]=paste0(data[[i]][1:d],collapse = "")
  }
  B22=unlist(B22)
  return(B22)
}

FasSeq=fastaString(fasta2)

FasSeqID=cbind(GeneID,FasSeq)
# dall$polypeptide_target_drug[757,]$amino_acid_sequence #  same "P30556"

# FasALL=left_join(PIK,FasSeqID,by="GeneID",copy=TRUE)
Gen=matrix(0L,length(PIK$GeneID),1)

for(i in 1:length(PIK$GeneID)){
  if (length(which(GeneID==PIK$GeneID[i]))==0)
  {
    Gen[i]=0 
  }
  else {Gen[i]=which(GeneID==PIK$GeneID[i]) }
  print(paste0("i = ",i))}


### takes a bit time to run 
PIK1=as.data.frame(PIK)
PIK1$Amino=0
length(PIK1$GeneID)
for(i in 1:length(PIK1$GeneID)){
  if (length(FasSeq[Gen[i]])!=0)
  {
    PIK1[i,]$Amino=FasSeq[Gen[i]]
  }
  else {    PIK1[i,]$Amino=0}
  print(paste0("i = ",i))}
PIK1[PIK1$parent_key=="DB00625",]$Amino="PISPIETVPVKLKPGMDGPKVKQWPLTEEKIKALVEICTEMEKEGKISKIGPENPYNTPVFAIKKKDSTKWRKLVDFRELNKRTQDFWEVQLGIPHPAGLKKRKSVTVLDVGDAYFSVPLDEDFRKYTAFTIPSINNETPGIRYQYNVLPQGWKGSPAIFQSSMTKILEPFRKQNPDIVIYQYMDDLYVGSDLEIGQHRTKIEELRQHLLRWGLTTPDKKHQKEPPFLWMGYELHPDKWTVQPIVLPEKDSWTVNDIQKLVGKLNWASQIYPGIRVRQLCKLLRGTKALTEVIPLTEEAELELAENREILKEPVHGVYYDPSKDLIAEIQKQGQGQWTYQIYQEPFKNLRTGKYARMRGAHTNDVKQLTEAVQKITTESIVIWGKTPKFKLPIQKETWETWWTEYWQATWIPEWEFVNTPPLVKLWYQLEKEPIVGAETFYVDGAANRETKLGKAGYVTNRGRQKVVTLTDTTNQKTELQAIYLALQDSGLEVNIVTDSQYALGIIQAQPDQSESELVNQIIEQLIKKEKVYLAWVPAHKGIGGNEQVDKLVSAGIRKVLFLDGID"

PIK4=dplyr::filter(PIK1,Amino!=0)    #### drop the ones w/o amino code
PIK2=dplyr::filter(PIK4,GeneID!="P49908" & GeneID!="P36969" & GeneID!="P59796" & GeneID!="P07203" & GeneID!="Q16881", GeneID!="P18283", GeneID!="P22352"  )   ### this amino-acid has U in and it gives an error, therefore we drop it.


CP1=combn(unique(PIK2$parent_key),2)
B=c(which(CP1[2,]==dru1Inquote),which(CP1[1,]==dru1Inquote))
CP2=CP1[,B]

OrdMat=function(CombMat,DrugNameMat){
  CP2=CombMat
  fdat2=DrugNameMat
  AAA=matrix(0,dim(fdat2)[1],1)
  AAA1=matrix(0,dim(fdat2)[1],1)
  for(i in 1:dim(fdat2)[1]){
    if (length(which(CP2==as.character(fdat2$drug_2)[i]))!=0)
    { 
      if(length(which(CP2[1,]==as.character(fdat2$drug_2)[i]))!=0)
      {
        AAA[i]=which(CP2[1,]==as.character(fdat2$drug_2)[i])
        AAA1[i]=as.character(fdat2$drug_2)[i]
      } else {
        AAA[i]=which(CP2[2,]==as.character(fdat2$drug_2)[i])
        AAA1[i]=as.character(fdat2$drug_2)[i]
      }
    } else {
      AAA[i]=0
      AAA1[i]=as.character(fdat2$drug_2)[i]
    }
  }
  list(Numbers=AAA,Names=AAA1)
}


AAA=OrdMat(CP2,Comb4)

  EnzSimFun=NULL
  PIK44= PIK2
  BBB=AAA$Numbers

  for(k in 1:length(BBB)){  ## 729 is length of AAA --- second column in FinalData2
    if (BBB[k]==0) {
      EnzSimFun[[k]]=0
    } else {
      filter1 = dplyr::filter(PIK44,parent_key==CP2[,BBB[k]][1])  ### filter to chosen drug and other drug names 
      filter2 = dplyr::filter(PIK44,parent_key==CP2[,BBB[k]][2])  ### filter to chosen drug and other drug names 
      mat = matrix(0L, nrow = dim(filter1)[1], ncol = dim(filter2)[1])
      for(i in 1:dim(filter1)[1]){
        for(j in 1:dim(filter2)[1]){
          
          seq1=AAString(filter1[,4][i])
          seq2=AAString(filter2[,4][j])
          seq12 = pairwiseAlignment(seq1, seq2, type = "local", 
                                    substitutionMatrix = BLOSUM62, scoreOnly = TRUE)
          seq11 = pairwiseAlignment(seq1, seq1, type = "local", 
                                    substitutionMatrix = BLOSUM62, scoreOnly = TRUE)
          seq22 = pairwiseAlignment(seq2, seq2, type = "local", 
                                    substitutionMatrix = BLOSUM62, scoreOnly = TRUE)
          
          if ( is.numeric(seq12) == FALSE | is.numeric(seq11) == FALSE | is.numeric(seq22) == FALSE ) {
            mat[i, j] = 0.0
          } else if ( abs(seq11) < .Machine$double.eps | abs(seq22) < .Machine$double.eps ) {
            mat[i, j] = 0.0
          } else {
            mat[i, j] = seq12/sqrt(seq11 * seq22)      ### cosine similarity
          }
        }
      }
      EnzSimFun[[k]]=mat
    }
    print(paste0("k = ",k))}


ProSMax = sapply(EnzSimFun, max)    ### get max cosine similarity out of others 


ProteinA=cbind(AAA$Numbers,AAA$Names,ProSMax)%>%data.frame()
names(ProteinA)=c("X1","drug_2","ProteinSim")

Comb5=left_join(Comb4,ProteinA,by="drug_2")%>%filter(X1!=0)%>%dplyr::select(.,-c(X1))





