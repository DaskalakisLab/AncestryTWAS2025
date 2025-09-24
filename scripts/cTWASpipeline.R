library(ctwas)
library(Rfast)
library(data.table)

# PTSD EUR GWAS Z-scores 
gwas <- read.table("/data/vbalakundi/GWAS/PTSD/Nievergelt_2024_NatGen/eur_ptsd_pcs_v4_aug3_2021_TWAS.txt", sep = " ", header =  T)
head(gwas)

z_snp <- read_gwas(gwas, id = 'SNP', A1 = 'A1', A2 = 'A2', z = 'Z')
head(z_snp)
gwas_n = 1222882

# AA Prediction model 
AACt6cisSNP_filepath <- "/data/ajajoo/PEC/PECEnet/Results/AAControl6cisSNP/fin/dbs/bulk_filtered.PF10.db"

# Reference data
genome_version <- get_predictdb_genome_build(AACt6cisSNP_filepath)

region_info <- readRDS("/data/vbalakundi/cTWAS/results/PTSD_AACt6cisSNP_UKBBLD_dbVariance_Union_10/region_info.Rds")

# Mapping reference SNPs and LD matrices to regions
# Making use of UKBB LD and SNP files 
# Downloaded from https://uchicago.app.box.com/s/jqocacd2fulskmhoqnasrknbt59x3xkn/folder/258913941342?page=1

LD_dir <- "/data/vbalakundi/cTWAS/data/LD_EUR/LD_UKBB_hg37"

LD_filestem <- sprintf("%s/ukb_%s_0.1_chr%s.R_snp.%s_%s", region_info$chrom,genome_version, region_info$chrom, region_info$start, region_info$stop)
region_metatable <- region_info
region_metatable$LD_file <- file.path(LD_dir, paste0(LD_filestem, ".RDS"))
region_metatable$SNP_file <- file.path(LD_dir, paste0(LD_filestem, ".Rvar"))

res <- create_snp_LD_map(region_metatable)
region_info <- res$region_info
snp_map <- res$snp_map
LD_map <- res$LD_map

ResultFolder <- "/data/vbalakundi/cTWAS/results/PTSD_AACt6cisSNP_UKBBLD_LDVariance_Union_10/"

saveRDS(z_snp, file = paste0(ResultFolder, "z_snp.Rds"))
saveRDS(snp_map, file = paste0(ResultFolder, "snp_map.Rds"))
saveRDS(LD_map, file = paste0(ResultFolder, "LD_map.Rds"))

z_snp <- readRDS(paste0(ResultFolder, "z_snp.Rds"))
snp_map <- paste0(ResultFolder, "snp_map.Rds")
LD_map <- readRDS(LD_map_file)

convert_to_ukb_varIDs <- function(varIDs) {
  varID_list <- strsplit(varIDs, split = "_|:")
  new_varIDs <- vapply(seq_along(varID_list), function(i) {
    x <- varID_list[[i]]
    if (length(x) >= 4) {
      sprintf("%s:%s_%s_%s", x[1], x[2], x[3], x[4])
    } else {
      varIDs[i] 
    }
  }, character(1))
  
  return(new_varIDs)
}


# Harmonizing GWAS z-scores and the reference data
z_snp <- preprocess_z_snp(z_snp, snp_map, 
                          drop_multiallelic = TRUE, 
                          drop_strand_ambig = TRUE,
                          varID_converter_fun = convert_to_ukb_varIDs)

# Harmonizing prediction models and the reference data
# default of options are used here
weights <- preprocess_weights( weight_file = weight_file,
                               region_info,
                               gwas_snp_ids = z_snp$id,
                               snp_map = snp_map,
                               LD_map = LD_map,
                               type = "count",
                               context = "DLPFC",
                               weight_name = "AACt6cisSNP",
                               weight_format = "PredictDB",
                               drop_strand_ambig = TRUE,
                               scale_predictdb_weights = TRUE,
                               load_predictdb_LD = load_predictdb_LD,
                               filter_protein_coding_genes = TRUE,
                               varID_converter_fun = convert_to_ukb_varIDs,
                               ncore = 6)

print("Processed Data Successfully")
# Before running the below command, make sure to delete the cor_matrix folder
# In the new update 05/01/2025, the "shared_all" is default option for group_prior_var_structure, previously it was "shared_type"
# thin variable default option is also changed from 0.1 to 1 (make use of all SNPs)
# min_group_size is 100 as default, I have reduced to 10 
ctwas_res <- ctwas_sumstats(z_snp, 
                            weights, 
                            region_info, 
                            LD_map, 
                            snp_map, 
                            thin = 1,
                            maxSNP = 20000,
                            min_group_size = 10, 
                            group_prior_var_structure = "shared_all", 
                            filter_L = TRUE,
                            filter_nonSNP_PIP = FALSE,
                            min_nonSNP_PIP = 0.5,
                            min_abs_corr = 0.1, 
                            ncore = 8, 
                            ncore_LD = 8,
                            save_cor = TRUE,
                            cor_dir =  paste0(ResultFolder,"/cor_matrix"),
                            force_compute_cor = FALSE)

z_gene <- ctwas_res$z_gene
param <- ctwas_res$param
finemap_res <- ctwas_res$finemap_res
susie_alpha_res <- ctwas_res$susie_alpha_res
boundary_genes <- ctwas_res$boundary_genes
region_data <- ctwas_res$region_data
screen_res <- ctwas_res$screen_res

print("Successful run of cTWAS")

library(EnsDb.Hsapiens.v75)
library(ggplot2)

convergence_plot <- make_convergence_plots(param, gwas_n)
ggsave(paste0(ResultFolder,"/convergence_plot.pdf"), plot = convergence_plot, width = 8, height = 6, dpi = 300)

ctwas_parameters <- summarize_param(param, gwas_n)
ctwas_parameters

# Inspecting and summarizing the cTWAS results
# Extract pval from zval
finemap_res$pval <- z2p(finemap_res$z)
head(finemap_res)

#subset(finemap_res, group != "SNP" & susie_pip > 0.8 & !is.na(cs))

# Adding gene annotations to cTWAS results
ens_db <- EnsDb.Hsapiens.v75
finemap_gene_res <- subset(finemap_res, group != "SNP")
gene_ids <- unique(finemap_gene_res$molecular_id)
gene_annot <- get_gene_annot_from_ens_db(ens_db, gene_ids)
colnames(gene_annot)[colnames(gene_annot) == "gene_id"] <- "molecular_id"
head(gene_annot)


finemap_res <- anno_finemap_res(finemap_res,
                                snp_map = snp_map,
                                mapping_table = gene_annot,
                                add_gene_annot = TRUE,
                                map_by = "molecular_id",
                                drop_unmapped = TRUE,
                                add_position = TRUE,
                                use_gene_pos = "mid")

print("Successful fine-mapping")
#subset(finemap_res, group != "SNP" & gene_type == "protein_coding" & susie_pip > 0.8 & !is.na(cs))
save(finemap_res,file = paste0(ResultFolder,"/FineMapResults.Rdata" ))

region_id_list <- c(unique(finemap_res$region_id))

for (i in 1:length(region_id_list)) {
  # Get the current region_id
  region_id <- region_id_list[i]
  
  # Generate the locusplot for the current region_id
  locusplot <- make_locusplot(finemap_res,
                              region_id = region_id,
                              ens_db = ens_db,
                              weights = weights,
                              highlight_pip = 0.8,
                              filter_protein_coding_genes = TRUE,
                              filter_cs = TRUE,
                              color_pval_by = "cs",
                              color_pip_by = "cs")
  
  ggsave(paste0(ResultFolder,"/locusplot_", region_id, ".pdf"), 
         plot = locusplot, width = 15, height = 6, dpi = 300)
}
