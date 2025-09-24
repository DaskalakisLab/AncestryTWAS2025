library(dplyr)
library(metaRNASeq)

addGroup <- function(bb)
{bb$group <- 0
bb$group[grepl("chr[0-9][p,q]",bb$Name)|
           grepl("chr[0-9][0-9][p,q]",bb$Name)] <- "cytoband"
bb$group[grepl("REACTOME",bb$Name)] <- "REACTOME"
bb$group[grepl("GOCC",bb$Name)] <- "GOCC"
bb$group[grepl("GOMF",bb$Name)] <- "GOMF"
bb$group[grepl("GOBP",bb$Name)] <- "GOBP"
bb$group[grepl("GO:",bb$Name)] <- "JustGO"
bb$group[grepl("KEGG",bb$Name)] <- "KEGG"
bb
}
# for numbers to not convert in scientific notation 
options(scipen=10)
setwd('~/mountMicc/ajajoo/')
# where the results are 
folderPath <- '~/mountMicc/ajajoo/PEC/PECJPEGMIX/ResultsAllcisWithPathway202306/'
# where JEPEGMIX results will be divded by models 
folderPathForMetalFormat <- '~/mountMicc/ajajoo/PEC/PECJPEGMIX/202404MetalFormat/Pathways/'
# where meta results will be stored
folderPathForMetalFiles <- '~/mountMicc/ajajoo/PEC/PECMetaOfPathways/202506LatestGWAS//'

system(paste0("mkdir ",folderPathForMetalFiles))
system(paste0("mkdir ",folderPathForMetalFiles,"/res"))
system(paste0("mkdir ",folderPathForMetalFiles,"/script"))

files <- basename(system(paste0('ls ',folderPath,'/*agn.txt'),intern = TRUE))
system(paste0('mkdir ',folderPathForMetalFormat))
filesWithMetalFormat <- basename(system(paste0('ls ',folderPathForMetalFormat,'/*.txt'),intern = TRUE))

# this needs more work it is not working 
files <- setdiff(files,filesWithMetalFormat)


source('~/mountEris/jajoo/PEC/PECAncestryPaper/scriptExplicityForPaper/BaseScript.R')

files[gsub(".agn.txt","",files) %in% GWASn[,1]]

for (ifile in files[gsub(".agn.txt","",files) %in% GWASn[,1]])
{
  GTAs <- read.table(paste0(folderPath,ifile),header=TRUE,stringsAsFactors = FALSE)
  GWASname <- gsub(".agn.txt","",ifile)
  GTAs$N <- GWASn[,2][match(GWASname,GWASn[,1])]
  GTAs$A1 <- 0
  GTAs$A2 <- 0 
  modelNames <- unique(GTAs$Tissue)
  for (imodel in modelNames)
  {
    modelGTAs <- GTAs %>% filter(Tissue == imodel)
    print('uncomment if running for the first time')
    write.table(modelGTAs,file=paste0(folderPathForMetalFormat,GWASname,"_",imodel,".txt"),quote =  FALSE,row.names = FALSE)
  }
}

## Run the following just to get modelNames 
for (ifile in files[1])
{
  GTAs <- read.table(paste0(folderPath,ifile),header=TRUE,stringsAsFactors = FALSE)
  GWASname <- gsub(".agn.txt","",ifile)
  GTAs$N <- GWASn[,2][match(GWASname,GWASn[,1])]
  GTAs$A1 <- 0
  GTAs$A2 <- 0 
  modelNames <- unique(GTAs$Tissue)
}
modelNames <- modelNames[!grepl("5",modelNames)]


fisherMeta <- function(x)
{
  rawpval <- list("pval1"=x$Pval.x,"pval2"=x$Pval.y)
  cbind(data.frame(fishercomb(rawpval, BHth = 0.05)[-1]),x)
  
}
# EU GWAS and EU enet model first 
################################
###-----READ THIS-----------####
# Check order of GWAS, it assumes first GWAS is EA and second AA
################################

ModelAppend <- "Control6"

for (Disease in c("BIP","PTSD","SCZ","MDD"))
{
  PairGWAS <- GWASn[GWASn[,4]==Disease,1]
  
  PairModel <- grep(ModelAppend,modelNames,value=TRUE)
  #PairModel <- c("EURAll9cisSNP","AAAll5cisSNP")
  #             EA_Enet    AA_Enet
  #   GWAS_EA    3           1
  #   GWAS_AA    4           2
  
  GWAS_Model_Combo<- paste0(folderPathForMetalFormat,
                            as.vector(outer(PairGWAS, PairModel, paste, sep="_")),
                            ".txt")
  comboInput <- cbind(c(1,1,3,3),c(2,4,2,4),paste0(ModelAppend,"_",
                                                   c("AA","mis","anc","EA")))
  
  # read template file 
  metalTemplate <- readLines('METAL/ASHG_METAL/METAL_SSB/script/Template.txt')
  
  for (icombo in 1:4)
  {
    outputfile <- paste0(folderPathForMetalFiles,"/",Disease,"_",comboInput[icombo,3],".csv")
    TWAS1 <- read.table(GWAS_Model_Combo[as.numeric(comboInput[icombo,1])],sep=" ",header = TRUE )
    TWAS2 <- read.table(GWAS_Model_Combo[as.numeric(comboInput[icombo,2])],sep=" ",header = TRUE )
    PrepareForMeta <- inner_join(TWAS1,TWAS2,by=c("Name"))
    PrepareForMeta <- addGroup(PrepareForMeta)
    fisherP <- PrepareForMeta %>% group_by(group) %>% group_map(~ fisherMeta(.x),keep=TRUE)
    fisherP <- as.data.frame(do.call(rbind,lapply(fisherP,as.data.frame)))
    write.csv(fisherP,file=outputfile)
  }
}


################################
###-----READ THIS-----------####
# Check order of GWAS, it assumes first GWAS is EA and second AA
################################

ModelAppend <- "ControlMax"

for (Disease in c("BIP","PTSD","SCZ","MDD"))
{
  PairGWAS <- GWASn[GWASn[,4]==Disease,1]
  PairModel <- c("AAControl6cisSNP","EURControl9cisSNP")
  #PairModel <- c("EURAll9cisSNP","AAAll5cisSNP")
  #             EA_Enet    AA_Enet
  #   GWAS_EA    3           1
  #   GWAS_AA    4           2
  
  GWAS_Model_Combo<- paste0(folderPathForMetalFormat,
                            as.vector(outer(PairGWAS, PairModel, paste, sep="_")),
                            ".txt")
  comboInput <- cbind(c(1,1,3,3),c(2,4,2,4),paste0(ModelAppend,"_",
                                                   c("AA","mis","anc","EA")))
  
  # read template file 
  metalTemplate <- readLines('METAL/ASHG_METAL/METAL_SSB/script/Template.txt')
  
  for (icombo in 1:4)
  {
    outputfile <- paste0(folderPathForMetalFiles,Disease,"_",comboInput[icombo,3],".csv")
    TWAS1 <- read.table(GWAS_Model_Combo[as.numeric(comboInput[icombo,1])],sep=" ",header = TRUE )
    TWAS2 <- read.table(GWAS_Model_Combo[as.numeric(comboInput[icombo,2])],sep=" ",header = TRUE )
    PrepareForMeta <- inner_join(TWAS1,TWAS2,by=c("Name"))
    PrepareForMeta <- addGroup(PrepareForMeta)
    fisherP <- PrepareForMeta %>% group_by(group) %>% group_map(~ fisherMeta(.x),keep=TRUE)
    fisherP <- as.data.frame(do.call(rbind,lapply(fisherP,as.data.frame)))
    write.csv(fisherP,file=outputfile)
    
    
  }
}


################################
###-----READ THIS-----------####
# Check order of GWAS, it assumes first GWAS is EA and second AA
################################

ModelAppend <- "AllMax"

for (Disease in c("BIP","PTSD","SCZ","MDD"))
{
  PairGWAS <- GWASn[GWASn[,4]==Disease,1]
  PairModel <- c("AAAll6cisSNP","EURAll9cisSNP")
  #PairModel <- c("EURAll9cisSNP","AAAll5cisSNP")
  #             EA_Enet    AA_Enet
  #   GWAS_EA    3           1
  #   GWAS_AA    4           2
  
  GWAS_Model_Combo<- paste0(folderPathForMetalFormat,
                            as.vector(outer(PairGWAS, PairModel, paste, sep="_")),
                            ".txt")
  comboInput <- cbind(c(1,1,3,3),c(2,4,2,4),paste0(ModelAppend,"_",
                                                   c("AA","mis","anc","EA")))
  
  # read template file 
  metalTemplate <- readLines('METAL/ASHG_METAL/METAL_SSB/script/Template.txt')
  
  for (icombo in 1:4)
  {
    outputfile <- paste0(folderPathForMetalFiles,Disease,"_",comboInput[icombo,3],".csv")
    TWAS1 <- read.table(GWAS_Model_Combo[as.numeric(comboInput[icombo,1])],sep=" ",header = TRUE )
    TWAS2 <- read.table(GWAS_Model_Combo[as.numeric(comboInput[icombo,2])],sep=" ",header = TRUE )
    PrepareForMeta <- inner_join(TWAS1,TWAS2,by=c("Name"))
    PrepareForMeta <- addGroup(PrepareForMeta)
    fisherP <- PrepareForMeta %>% group_by(group) %>% group_map(~ fisherMeta(.x),keep=TRUE)
    fisherP <- as.data.frame(do.call(rbind,lapply(fisherP,as.data.frame)))
    write.csv(fisherP,file=outputfile)
    
  }
}



################################
###-----READ THIS-----------####
# Check order of GWAS, it assumes first GWAS is EA and second AA
################################

# Run for each model on both GWAS 

for (Disease in c("BIP","PTSD","SCZ","MDD"))
{
  PairGWAS <- GWASn[GWASn[,4]==Disease,1]
  # for each model run on both GWAS 
  
  for (imodel in modelNames[!grepl("5",modelNames)])
  {  
    GWAS_Model_Combo<- paste0(folderPathForMetalFormat,
                              as.vector(outer(PairGWAS, imodel, paste, sep="_")),
                              ".txt")
    typeOfCombo <- ifelse(grepl("Mix",imodel),"Mix",ifelse(grepl("AA",imodel) ,"AA","EUR"))
    comboInput <- cbind(c(1),c(2),paste0(imodel,"_",typeOfCombo))
    
    # read template file 
    metalTemplate <- readLines('METAL/ASHG_METAL/METAL_SSB/script/Template.txt')
    
    for (icombo in 1)
    {
      outputfile <- paste0(folderPathForMetalFiles,Disease,"_",comboInput[icombo,3],".csv")
      TWAS1 <- read.table(GWAS_Model_Combo[as.numeric(comboInput[icombo,1])],sep=" ",header = TRUE )
      TWAS2 <- read.table(GWAS_Model_Combo[as.numeric(comboInput[icombo,2])],sep=" ",header = TRUE )
      PrepareForMeta <- inner_join(TWAS1,TWAS2,by=c("Name"))
      PrepareForMeta <- addGroup(PrepareForMeta)
      fisherP <- PrepareForMeta %>% group_by(group) %>% group_map(~ fisherMeta(.x),keep=TRUE)
      fisherP <- as.data.frame(do.call(rbind,lapply(fisherP,as.data.frame)))
      write.csv(fisherP,file=outputfile)
    } 
  }
}
