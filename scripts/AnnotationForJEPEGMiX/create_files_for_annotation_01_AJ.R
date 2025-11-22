create_files_for_annotation <- function(db_file, output_folder){
#prepare environment.
library(dplyr)
library(RSQLite)
library(data.table)
#library(tidyverse)

  my_y=c("SNP", "CHR", "BP" , "A1" , "A2","weight" )
  
  # Read in DB file
  driver <- dbDriver('SQLite')
  print(db_file)
  in_conn <- dbConnect(driver, db_file)
  print("connected")
  # Extract weights per SNP per Gene
  weights_PEC <- dbGetQuery(in_conn, 'select * from weights')
  # Extract stats per Gene
  extra_PEC <- dbGetQuery(in_conn, 'select * from extra')
  
  # Extract weights per gene per idp field
  idp <- data.frame(matrix(ncol=6,nrow=0))
  names(idp) <- c("gene","rsid","varID","ref_allele","eff_allele","weight")
  #idp<-dbGetQuery(in_conn, 'select * from bioW')
  dbDisconnect(in_conn)
  
  # Create list of genes from the weights/extra tables
  a1=extra_PEC$gene
  a2=unique(weights_PEC$gene)
  
  # Determine the genes which are not in both
  dd=setdiff(a1,a2)
  
  # Set Extra rownames to the genes
  row.names(extra_PEC)=a1
  # Extract only the common genes
  extra_PEC=extra_PEC[a2,]
  
  # For all the genes in the weights file
  for(i in 1:length(a2)){
    
    # Extract the weights per snp for this gene
    dtemp=subset(weights_PEC,weights_PEC$gene==a2[i])
    
    # Extract the idp fields that get weights for this gene and the weights
    dtemp2=subset(idp,idp$gene==a2[i])
    
    # If the gene has weights for an idp field
    if(nrow(dtemp2)!=0){
      print("Case1")
      
      # Loop through all weighted fields
      for(j in 1:nrow(dtemp2)){
        
        # Extract the number of characters before "_imputed" for the idp field
        pos=gregexpr(pattern ="_imputed",dtemp2$rsid[j])
        pos=unlist(pos)
        # Extract idp field name without "_imputed" could be lines 58-61 could be replaced with gsub
        nn=substr(dtemp2$rsid[j],1,pos[1]-1)
        
        # Read in the imputed weight files for the given idp field as a data frame
        #dx=fread(paste0("~/sshfs/micc/cdipietro/JEPEGMIX/PEC/Data/IDPs/",nn,".txt"),header=T)
        dx=fread(paste0("/data/cdipietro/JEPEGMIX/PEC/Data/IDPs/",nn,".txt"),header=T)
        dx=as.data.frame(dx)
        # Multiply all weights by the weights_PEC for this gene in this idp field
        dx$weight=dx$weight*dtemp2$weight[j]
        
        # Add to final data frame object
        if(j==1){dfin=dx}
        if(j!=1){dfin=rbind(dfin,dx)}
        
      }
      
    }
    
    # Extract the chromosome and BP for each SNP for the given gene
    bp=0;chr=0
    for(j in 1:nrow(dtemp)){
      
      pos=gregexpr(pattern ="_",dtemp$varID[j])
      pos=unlist(pos)
      chr[j]=substr(dtemp$varID[j],1,pos[1]-1)
      bp[j]=substr(dtemp$varID[j],pos[1]+1,pos[2]-1)
      
    }
    
    # Add chromosome and BP for each SNP
    dtemp$CHR=chr
    dtemp$bp=bp
    # Change ref_allele and eff_allele to A2 and A1
    colnames(dtemp)[c(4,5)]=c("A2","A1")
    # Remove gene and varID columns
    dtemp=dtemp[,-c(1,3)]
    # Reorder the columns to be "rsid", "CHR", "bp", "A1", "A2", "weight"
    dtemp=dtemp[,c(1,5,6,3,2,4)]
    
    # If the gene has weights for an idp field
    if(nrow(dtemp2)!=0){
      print("Case2")
      # Remove Feature column
      dfin=dfin[,-ncol(dfin)]
      # Rename columns
      colnames(dtemp)=colnames(dfin)
      # Combine data frames
      dg=rbind(dtemp,dfin)
      # Extract list of unique SNPs
      am=unique(dg$SNP)
      # For each SNP
      for(j in 1:length(am)){
        # Extract SNP row
        dtemp=subset(dg,dg$SNP==am[j])
        if(nrow(dtemp)!=0){
          # Sum together weights if more than 1
          dtemp$weight[1]=sum(dtemp$weight)
          dxm=dtemp[1,]
        } 
        if(nrow(dtemp)==0){
          dxm=dtemp
        }
        # Combine in data frame
        if(j==1){dfinal=dxm}
        if(j!=1){dfinal=rbind(dfinal,dxm)}
      }
    }
    
    # If the gene has no weights for any idp fields just sum weights from weights_PEC
    if(nrow(dtemp2)==0){
      colnames(dtemp)=my_y
      
      dg=dtemp
      
      am=unique(dg$SNP)
      
      for(j in 1:length(am)){
        
        dtemp=subset(dg,dg$SNP==am[j])
        
        if(nrow(dtemp)!=1){
          
          dtemp$weight[1]=sum(dtemp$weight)
          dxm=dtemp[1,]
        } 
        
        if(nrow(dtemp)==1){
          dxm=dtemp
        } 
        
        if(j==1){dfinal=dxm}
        if(j!=1){dfinal=rbind(dfinal,dxm)}
      }
    }
    
    # Reorder columns "CHR", "SNP", "BP", "A1", "A2", "weight"
    # Convert to numeric where applicable
    dfinal=dfinal[,c(2,1,3,4,5,6)]
    dfinal$CHR=as.numeric(dfinal$CHR)
    dfinal$BP=as.numeric(dfinal$BP)
    dfinal$weight=as.numeric(dfinal$weight)
    
    # Order by chromosome and BP
    dfinal=dfinal[order(dfinal$CHR,dfinal$BP),]
    dir.create(paste0(output_folder, "/for_plink/"), recursive = T)
    dir.create(paste0(output_folder, "/RDATA/"), recursive = T)
    
    
    ff=paste0(output_folder, "/for_plink/",extra_PEC$gene[i],"_or_",extra_PEC$genename[i],".txt")
    fff=paste0(output_folder, "/RDATA/",extra_PEC$gene[i],"_or_",extra_PEC$genename[i],".RData")
    
    # Write dfinal to txt file
    write.table(dfinal,ff,sep="\t",quote=F ,row.names=F,col.names=F,append=F)
    # Add gene and Symbol info and save to RData object
    dfinal$gene1=extra_PEC$gene[i]
    dfinal$gene2=extra_PEC$genename[i]
    dfinal$gene_type=extra_PEC$gene_type[i]
    save(dfinal,file=fff)
    
  }
}


