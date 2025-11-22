
create_annotation <- function(pwd11,modelName){
library(data.table) 
 out_folder <- pwd11
  pwd11 <- paste0(pwd11,"/RDATA/")
  
  atble11=list.files(pwd11)
  formula=0
  gene=0
  Type=0
  bp_start=0
  bp_end=0
  chr=0
  
  
  for(i in 1:length(atble11)){
    load(paste0(pwd11,atble11[i]))
    ax=0
    for(j in 1:nrow(dfinal)){
      if(j==1){
        ax=paste0(dfinal$weight[j],"*",dfinal$SNP[j],"_",dfinal$A1[j],"_",dfinal$A2[j])
      }
      if(j!=1){
        if(dfinal$weight[j]>0){
          ax=paste0(ax,"+",paste0(dfinal$weight[j],"*",dfinal$SNP[j],"_",dfinal$A1[j],"_",dfinal$A2[j]))
        }
        if(dfinal$weight[j]<0){
          ax=paste0(ax,paste0(dfinal$weight[j],"*",dfinal$SNP[j],"_",dfinal$A1[j],"_",dfinal$A2[j]))
        }
      }  }
    
    formula[i]=ax
    gene[i]=dfinal$gene2[1]
    Type[i]=dfinal$gene_type[1]
    chr[i]=dfinal$CHR[1]
    bp_start[i]=min(as.numeric(dfinal$BP))
    bp_end[i]=max(as.numeric(dfinal$BP))
   
  }
  
  
  
  
  
  # Type Tissue Name chr bp_start bp_end Q Formula
  # protein_coding DLPFC_PEC
  df=as.data.frame(cbind(gene,formula,Type,bp_start,bp_end,chr))
  #df$Type="protein_coding"
  df$Tissue=modelName
  df$Q=0.05
  
  df=df[,c("Type","Tissue","gene","chr","bp_start", "bp_end", "Q", "formula")]
  
  filename=paste0(out_folder, "/annotation_cis.txt")
  write.table(df,filename,sep=" ",quote=F ,row.names=F,col.names=T,append=F)  
  
}


