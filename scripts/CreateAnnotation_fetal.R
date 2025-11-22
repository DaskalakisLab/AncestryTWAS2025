library(dplyr)
weightsFolder <- "~/mountMicc/ajajoo/GandalLab/fetalFromSynapse/WEIGHTS/"
positionFile <- read.table("~/mountMicc/ajajoo/GandalLab/fetalFromSynapse/WEIGHTS//gene_all.pos", header = TRUE)

annotationTable <- data.frame()
weightsFiles <- list.files(weightsFolder, pattern = "\\.wgt\\.RDat$")
print(length(weightsFiles))


for (weightFile in weightsFiles) {
  
  geneWeights <- load(paste0(weightsFolder, weightFile))
  
  if (!is.na(cv.performance[1, 2]) && !is.na(cv.performance[2, 2])) {
    if (cv.performance[1, 2] > 0.01 & cv.performance[2, 2] < 0.05) {
      
      wgt.matrix <- data.frame(wgt.matrix)
      wgt.matrix <- wgt.matrix %>% filter(enet != 0) %>% select(enet)
      
      if (nrow(wgt.matrix > 0)){
        wgt.matrix$SNP <- row.names(wgt.matrix)
        SNPweights <- merge(wgt.matrix, snps, by.x = "SNP",by.y = "V2")
        for (j in c(1:nrow(SNPweights))){
          {
            if(j==1){
              if(SNPweights$enet[j]>0){
                ax=paste0("+",SNPweights$enet[j],"*",SNPweights$SNP[j],"_",SNPweights$V5[j],"_",SNPweights$V6[j])
              }
              if(SNPweights$enet[j]<0){
                ax=paste0(SNPweights$enet[j],"*",SNPweights$SNP[j],"_",SNPweights$V5[j],"_",SNPweights$V6[j]) 
              }
            }
            if(j!=1){
              if(SNPweights$enet[j]>0){
                ax=paste0(ax,"+",paste0(SNPweights$enet[j],"*",SNPweights$SNP[j],"_",SNPweights$V5[j],"_",SNPweights$V6[j]))
              }
              if(SNPweights$enet[j]<0){
                ax=paste0(ax,paste0(SNPweights$enet[j],"*",SNPweights$SNP[j],"_",SNPweights$V5[j],"_",SNPweights$V6[j]))
              }
            }  
          }
          
        }
        
        geneID <- weightFile
        positionFile_sub <- positionFile %>% filter(WGT == geneID)
        
        if (nrow(positionFile_sub) > 0) {
          annotationRow <- data.frame(
            Type = "Gene",
            Tissue = "fetal",
            gene_name = positionFile_sub$ID,
            chr = positionFile_sub$CHR,
            bp_start = positionFile_sub$P0,
            bp_end = positionFile_sub$P1,
            Q = 0.05,
            Formula = ax,
            gene_id = gsub(".wgt.RDat", "", weightFile)
          )
          
          annotationTable <- rbind(annotationTable, annotationRow)
        }
      }
    }
  }
}

write.table(annotationTable, "~/mountEris/jajoo/PEC/PECAncestryPaper/AncestryRevision/data/GandalLab/AnnotationFiles/fetalAnnotationFile.txt",
            quote = FALSE, row.names = FALSE)
