
finPaths <- list.dirs(path="/data/ajajoo/PEC/PECEnet/Results/",recursive=TRUE,full.names=TRUE)
finPaths<- paste0(finPaths[grep("fin$",finPaths)],"/")
finPaths<- finPaths[grep("cisSNP",finPaths)]
finPaths<- finPaths[-grep("EUR1",finPaths)]
finPaths<- finPaths[!grepl("toDelete",finPaths)]
finPaths<- finPaths[!grepl("MixAllcisSNP",finPaths)]
finPaths<- finPaths[!grepl("MixControlcisSNP",finPaths)]
finPaths<- finPaths[!grepl("geneOnly",finPaths)]


source('/data/ajajoo/PEC/PECJPEGMIX/script/create_files_for_annotation_01_AJ.R')
source('/data/ajajoo/PEC/PECJPEGMIX/script/create_files_for_annotation_cis_AJGeneName.R')
#resultFolder <- ReplaceResultFolder
resultFolder <- "/data/ajajoo/PEC/PECJPEGMIX/June2023/"
if (file.exists(resultFolder)) break 

system(paste0("mkdir ",resultFolder))

patternToIdentifyModel <- "cisSNP"
for (pathToDb in finPaths)
{
temp <- unlist(strsplit(pathToDb,split="/"))  
modelName <- temp[grep(patternToIdentifyModel,temp)]
out_dir <- paste0(resultFolder,modelName)
system(paste0("mkdir ",out_dir))
create_files_for_annotation(paste0(pathToDb,"/dbs/bulk_filtered.PF10.db"), out_dir)
create_annotation(out_dir,modelName)
}

# take header 
system(paste0("cat ",resultFolder,"/*/annotation*.txt | head -n 1 > ",resultFolder,"/AnnotationForAll.txt"))
# merge all annotation with modelName in tissue name to differentiate them 
system(paste0("ls ",resultFolder,"/*/annotation*.txt | xargs -I {} tail -n +2 {} >> ",resultFolder,"/AnnotationForAll.txt"))

# Add pathway annotation 
system(paste0("cat /data/cchatzinakos/ancestry/Pathways_updatated_v1.txt >> ",resultFolder,"/AnnotationForAll.txt"))
