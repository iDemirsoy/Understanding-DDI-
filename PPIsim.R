## this file used to named as GeneSimSon.R 
## Tje mamed changed to PPIsim.R


##### this R file computes Gene Similarity between 2 drugs ##### 
## give second drugs name fdat2$drug_2 

## head(fdat2[,c("drug_1.x","drug_2")])
## drug_1.x  drug_2
## 1  DB00208 DB00035
## 2  DB00208 DB00091
## 3  DB00208 DB00104
## 4  DB00208 DB01303
## 5  DB00208 DB00130
#
#


setwd("/Users/idemirsoy/1DriveMyFsu/OneDrive - Florida State Students/925/MyoKardiaDownloads/PvSummerCodes/SimilarityMatrices")
library(dplyr)
library(igraph)

#dru1Inquote="DB00208"
load("dall.RData")

#chem = read.csv('CTD_chemicals.csv.gz', comment.char = "#",header=F)

gene = read.csv('CTD_genes.csv.gz', comment.char = "#",header=F)

gene1=dplyr::filter(gene,!is.na(V3))



gene2= gene1[, c(3, 8)]

gene3=gene2[gene2$V8!="",]

names(gene3)=c( "geneID","UniprotID")



geneL=gene3 %>%tidyr::separate_rows(UniprotID, sep = '\\|')%>%data.frame()





PK1=dall$carriers_drug[,c("id","parent_key")]

PK2=dall$transporters_drug[,c("id","parent_key")]

PK3=dall$targ_drug[,c("id","parent_key")]

PK4=dall$enzymes_drug[,c("id","parent_key")]



PKcmob=rbind(PK1,PK2,PK3,PK4)

PKnoDup=distinct(PKcmob,id,parent_key)





PI1=dall$polypeptides_carrier_drug[,c("id","parent_id")]

PI2=dall$polypeptides_transporter_drug[,c("id","parent_id")]

PI3=dall$polypeptide_target_drug[,c("id","parent_id")]

PI4=dall$polypeptides_enzyme_drug[,c("id","parent_id")]





PIcmob=rbind(PI1,PI2,PI3,PI4)

PInoDup=distinct(PIcmob,id,parent_id)



names(PInoDup)=c("UniprotID","id")



PiPk=left_join(PInoDup,PKnoDup,by="id")



GeneComb1=inner_join(geneL,PiPk,by="UniprotID")                   ## inner_join will do the job







## bioGrid gene interactions

##

biogrid = read.table('B1O2.txt',
                     sep = '\t', header = T, fill = TRUE,
                     stringsAsFactors = FALSE, quote = '')



biogrid = biogrid[,c("SWISS.PROT.Accessions.Interactor.A","SWISS.PROT.Accessions.Interactor.B")]

biogrid1=biogrid[biogrid$SWISS.PROT.Accessions.Interactor.A!="",]

biogrid.1=biogrid1[biogrid1$SWISS.PROT.Accessions.Interactor.B!="-",]



names(biogrid.1)=c("bi1","bi2")





BB=biogrid.1 %>%
  
  tidyr::separate_rows(bi1, sep = '\\|')%>%  tidyr::separate_rows(bi2, sep = '\\|')%>%data.frame()





# expand duplicated DrugBank IDs

#biotmp1 = vector('list', nrow(biogrid.1))

#for (i in 1:nrow(biogrid.1)) biotmp1[[i]] = as.character(biogrid.1[i, 1])

#names(biotmp1) = as.character(biogrid.1[, 2])

#biotmp1 = sapply(biotmp1, strsplit, split = '\\|')

#tmpBio1 = unlist(biotmp1)

#Bio1L = data.frame(bi1 = tmpBio1,bi2 = names(tmpBio1))



#biotmp2 = vector('list', nrow(Bio1L))

#for (i in 1:nrow(Bio1L)) biotmp2[[i]] = as.character(Bio1L[i, 2])

#names(biotmp2) = as.character(Bio1L[, 1])

#biotmp2 = sapply(biotmp2, strsplit, split = '\\|')

#tmpBio2 = unlist(biotmp2)

#Bio1L2 = data.frame(bi1 = names(tmpBio2),bi2 = tmpBio2)





#B1=strsplit(Bio1L2[,1],split = '')

#B3=strsplit(Bio1L2[,2],split = '')

String1=function(data){
  
  B1=data
  
  B2=NULL
  
  for (i in 1:length(B1)) {if(length(B1[[i]]==6))
    
  {
    
    B2[[i]]=paste0(B1[[i]][1],B1[[i]][2],B1[[i]][3],B1[[i]][4],B1[[i]][5],B1[[i]][6])
    
  } else{
    
    B2[[i]]=paste0(B1[[i]][1],B1[[i]][2],B1[[i]][3],B1[[i]][4],B1[[i]][5],B1[[i]][6])
    
  }
    
  }
  
  B2=unlist(B2)
  
  return(B2)
  
}



#BB=cbind(String1(B1),String1(B3))



BB=dplyr::filter(BB,bi1!="-")







g = graph.data.frame(BB, directed = FALSE)

nodes = unique(as.vector(as.matrix(BB)))



CP1=combn(unique(GeneComb1$parent_key),2)

B=c(which(CP1[2,]==dru1Inquote),which(CP1[1,]==dru1Inquote))

CP2=CP1[,B]



#Comb7=read.csv("Comb7.csv",header=T,sep=",")





A=matrix(0,dim(Comb7)[1],1)

A1=matrix(0,dim(Comb7)[1],1)

for(i in 1:dim(Comb7)[1]){
  
  if (length(which(CP2==as.character(Comb7$drug_2)[i]))!=0)
    
  {
    
    if(length(which(CP2[1,]==as.character(Comb7$drug_2)[i]))!=0)
      
    {
      
      A[i]=which(CP2[1,]==as.character(Comb7$drug_2)[i])
      
      A1[i]=as.character(Comb7$drug_2)[i]
      
    } else {
      
      A[i]=which(CP2[2,]==as.character(Comb7$drug_2)[i])
      
      A1[i]=as.character(Comb7$drug_2)[i]
      
    }
    
  } else {
    
    A[i]=0
    
    A1[i]=as.character(Comb7$drug_2)[i]
    
  }
  
}





AA = 0.9 * exp(1)



Begin=Sys.time()



L=NULL

for (k in 1:length(A)){  ## 729 is length of A --- second column in FinalData2
  
  if (A[k]==0) {
    
    L[[k]]=0
    
  } else {
    
    id1 = dplyr::filter(GeneComb1,parent_key==CP2[,A[k]][1])
    
    id2 = dplyr::filter(GeneComb1,parent_key==CP2[,A[k]][2])
    
    mat = matrix(0L, nrow = dim(id1)[1], ncol = dim(id2)[1])
    
    for(i in 1:dim(id1)[1]){
      
      for(j in 1:dim(id2)[1]){
        
        gid1 = as.character(id1[,2][i])
        
        gid2 = as.character(id2[,2][j])
        
        if (gid1 == gid2) {
          
          mat[i, j] = 1
          
        } else if ( (gid1 %in% nodes) & (gid2 %in% nodes) ) {
          
          spath = length((get.shortest.paths(g, from = gid1, to = gid2, output = 'epath')$epath)[[1]])
          
          mat[i, j] = AA * ( exp(1)^(-spath) )
          
        } else {
          
          mat[i, j] = 0
          
        }
        
      }
      
    }
    
    L[[k]]=mat
    
    print(paste0("k = ",k))
    
  }
  
}

End=Sys.time()

End-Begin



GeneSim = sapply(L, mean)

GenAB=cbind(A,A1,GeneSim)%>%data.frame()

names(GenAB)=c("X1","drug_2","GeneSim")



GenAB$drug_2=as.character(GenAB$drug_2) 

Comb7$drug_2=as.character(Comb7$drug_2) 



Comb8=left_join(Comb7,GenAB,by="drug_2")%>%filter(X1!=0)%>%dplyr::select(.,-c(X1))

names(Comb8)[names(Comb8)=="drug_1.x"]=c("drug_1")





#write.csv(Comb8,"Esomeprazole.csv")