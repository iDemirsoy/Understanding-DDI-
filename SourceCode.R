#############################################################################
###                                                                         #
###      This is a source code to run all R codes in order                  #
###      Use structure.py to get two csv files for the dru1Inquote          #
###      dru1Inquote is the drug's drugbank code number                     #
###      And use them in Mol1 and Mol2                                      #
###      The memory of the laptop is important - it might not run last two  #
###                                                                         #
###                                                               iDemirsoy #
#############################################################################


setwd("/Users/idemirsoy/1DriveMyFsu/OneDrive - Florida State Students/925/MyoKardiaDownloads/PvSummerCodes/SimilarityMatrices/")
library(dplyr)

Mol1=read.csv("DB00335A.csv",sep=",",header = T)
Mol2=read.csv("DB00335B.csv",sep=",",header = T)
names(Mol2)=c("X","compound1","compound0","similarity")
Mol2=dplyr::select(Mol2,c(X,compound0,compound1,similarity))
Molec=rbind(Mol1,Mol2)
names(Molec)[names(Molec)=="compound0"]=c("drug_1")
names(Molec)[names(Molec)=="compound1"]=c("drug_2")
names(Molec)[names(Molec)=="similarity"]=c("MolSim")
MolSim1=dplyr::select(Molec,drug_1,drug_2,MolSim)

dru1Inquote="DB00335"

source("ATCsim.R")      ## Comb2
source("OneSiDB.R")     ## Comb3
source("DiseSim.R")     ## Comb4
source("ProteinSim.R")  ## Comb5
source("EnzymeSim.R")   ## Comb6
source("PathWaySim.R")  ## Comb7
source("PPIsim.R")      ## Comb8 


#write.csv(Comb8,"Atenolol.csv")

## add interactions 
druginter=dall$interactions_drug

drug_1=druginter[druginter$parent_key==dru1Inquote,]$parent_key
drug_2=druginter[druginter$parent_key==dru1Inquote,]$`drugbank-id`

DD=cbind(drug_1,drug_2)
DD=as.data.frame(DD)
DD$Y=1

DD1=left_join(Comb8,DD,by="drug_2")%>%dplyr::select(-drug_1.y)
#names(DD1)[names(DD1) == "drug_1.x"] <- "drug_1" # -- in case we need 
DD1$Y=ifelse(is.na(DD1$Y),0,1)

write.csv(DD1,"DB00335/Atenolol.csv")






