###############################################################################
## Dall is dragbank dataset 
## uniprot_sprot.fasta.gz is downloaded from uniprot.org - Reviewed(Swiss-Prot)
## 


#BiocManager::install("Biostrings")
#install.packages("seqinr")
library(seqinr)
library(dplyr)
library(Biostrings)

setwd("/Users/idemirsoy/1DriveMyFsu/OneDrive - Florida State Students/925/MyoKardiaDownloads/PvSummerCodes/SimilarityMatrices")

#dru1Inquote="DB01242"
#load("dall.RData")

## this file has 2 rows 
# 1 row is name of drug_1
# 2 row is name of drug_2


## fatsa file 
#fasta2 = read.fasta('uniprot_sprot.fasta.gz',seqtype = "AA")
#fasta22=names(fasta2)
#B1= strsplit(fasta22,split='')
## target gene information

PK4=dall$enzymes_drug[,c("id","parent_key")]
PI4=dall$polypeptides_enzyme_drug[,c("id","parent_id")]
names(PI4)=c("GeneID","id")
PIK4=left_join(PI4,PK4,by="id")

data(BLOSUM62)

#String2=function(data,from,to){
#  B2=NULL
#  for(i in 1:length(data)){B2[[i]]=paste0(data[[i]][from:to],collapse = "")}
#  B2=unlist(B2)
#  return(B2)
#}
# GeneID=String2(B1,4,9)

# fastaString=function(data){
#  B22=NULL
#  for(i in 1:length(data)){
#    d=length(data[[i]]);
#    B22[[i]]=paste0(data[[i]][1:d],collapse = "")
#  }
#  B22=unlist(B22)
#  return(B22)
#}

#FasSeq=fastaString(fasta2)

#FasSeqID=cbind(GeneID,FasSeq)
# dall$polypeptide_target_drug[757,]$amino_acid_sequence #  same "P30556"

# FasALL=left_join(PIK,FasSeqID,by="GeneID",copy=TRUE)
Enz=matrix(0L,length(PIK4$GeneID),1)

for(i in 1:length(PIK4$GeneID)){
  if (length(which(GeneID==PIK4$GeneID[i]))==0)
  {
    Enz[i]=0 
  }
  else {Enz[i]=which(GeneID==PIK4$GeneID[i]) }
  print(paste0("i =",i))
}


PIK41=as.data.frame(PIK4)
PIK41$Amino=0
for(i in 1:length(PIK41$GeneID)){
  if (length(FasSeq[Enz[i]])!=0)
  {
    PIK41[i,]$Amino=FasSeq[Enz[i]]
  }
  else {    PIK41[i,]$Amino=0}
}


PIK433=dplyr::filter(PIK41,Amino!=0)    #### drop the ones w/o amino code
PIK44=dplyr::filter(PIK433,GeneID!="P05164" & GeneID!="P07203" & GeneID!="Q16881" & GeneID!="Q9NZV6")   ### this amino-acid has U in and it gives an error, therefore we drop it.


CE1=combn(unique(PIK44$parent_key),2)
E1=c(which(CE1[2,]==dru1Inquote),which(CE1[1,]==dru1Inquote))
CE2=CE1[,E1]


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

BBB=OrdMat(CE2,Comb5)

### this function computes cosine similarity for enzyme gene similarity between 
BBn=BBB$Numbers
EnzSimFunc=list()
for(k in 1:length(BBn)){  ## 729 is length of AAA --- second column in FinalData2
  if (BBn[k]==0) {
    EnzSimFunc[[k]]=0
  } else {
    filter1 = dplyr::filter(PIK44,parent_key==CE2[,BBn[k]][1])  ### filter to chosen drug and other drug names 
    filter2 = dplyr::filter(PIK44,parent_key==CE2[,BBn[k]][2])  ### filter to chosen drug and other drug names 
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
    EnzSimFunc[[k]]=mat}
    print(paste0("k = ",k))
  }



#EnzSimFunMax = sapply(EnzSimFun, max)    ### get max cosine similarity out of others 
EnzSimFunAvg = sapply(EnzSimFunc, mean)
#EnzSimFunmin = sapply(EnzSimFun, min)
#EnzSimFunsum = sapply(EnzSimFun, sum)
#EnzSimFunMax = sapply(EnzSimFun, max)


EnzymeA=cbind(BBB$Numbers,BBB$Names,EnzSimFunAvg)%>%data.frame()
names(EnzymeA)=c("X1","drug_2","EnzymeSim")

Comb6=left_join(Comb5,EnzymeA,by="drug_2")%>%filter(X1!=0)%>%dplyr::select(.,-c(X1))

#write.csv(Comb6,"DB00321/Comb6321.csv")

#Comb6=inner_join(Comb5,EnzSimilarity,by="drug_2")
