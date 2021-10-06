#!/usr/bin/env Rscript
# options(encoding = "UTF-8")
# work_dir = "~/Documents/work/ORCESTRA/NCI60/"
# setwd(work_dir)
library(PharmacoGx)
library(org.Hs.eg.db)
library(GenomicRanges)
library(SummarizedExperiment)

input_dir <- '/pfs/getNCI60/'
out_dir <- '/pfs/out/'
# input_dir <- ''
# out_dir <- ''

drug_all <- read.csv("/pfs/downAnnotations/drugs_with_ids.csv", na.strings=c("", " ", "NA"))

saveRDS(drug_all, '/pfs/out/test.rds')

# drug.obj <- readRDS(paste0(input_dir, 'drug_obj.rds'))
# cell.obj <- readRDS(paste0(input_dir, 'cell_obj.rds'))
# cell.obj.sen <- readRDS(paste0(input_dir, 'cell.obj.sen.rds'))
# info.sensitivity <- readRDS(paste0(input_dir, 'info.sensitivity.rds'))
# profile.sensitivity <- readRDS(paste0(input_dir, 'profile.sensitivity.rds'))
# raw.sensitivity <- readRDS(paste0(input_dir, 'raw.sensitivity.rds'))
# drug_with_ids <- readRDS(paste0(input_dir, 'drug_with_ids.rds'))

# RNA_SE <- readRDS(paste0(input_dir, 'RNA_SE.rds'))
# miRNA_SE <- readRDS(paste0(input_dir, 'miRNA_SE.rds'))
# RNA_seq_comp_SE <- readRDS(paste0(input_dir, 'RNA_seq_comp_SE.rds'))
# RNA_seq_iso_SE <- readRDS(paste0(input_dir, 'RNA_seq_iso_SE.rds'))


# ########################################### Curation ##########################################################
# ######### Drug #########
# cur.drug <- data.frame(unique.drugid = drug.obj$drugid,
#                        NCI60.drugid = drug.obj$NCI60.drugid,
#                        row.names = drug.obj$drugid)

# ######### Cell #########
# cur.cell <- data.frame(unique.cellid = cell.obj$cellid,
#                        NCI60.cellid = cell.obj$NCI60.cellid,
#                        row.names= cell.obj$cellid)

# ######### Tissue #########
# cur.tissue <- data.frame(unique.tissueid = cell.obj$tissueid, 
#                          NCI60.tissueid = cell.obj$NCI60.tissueid , 
#                          row.names = cell.obj$cellid)

# ########################### Sensitivity - continue ##################
# info.sensitivity <- merge(info.sensitivity , cell.obj.sen[, c("cellid" ,"tissueid", "CELL_NAME")], by="CELL_NAME" , all.x=T) # Cellid column is from cell.obj so contains unique.ids

# # Mapping drugids back to UNcontatenated NSC id to be able to map to senesitivity objects
# uniq.nsc.map <- drug_with_ids[ !is.na(drug_with_ids$NSC_id), c("unique.drugid" , "NSC_id", "cid" )]

# while (length(setdiff(unique(info.sensitivity$NSC) , uniq.nsc.map$NSC_id))>1){
  
#   for (i in seq(nrow(uniq.nsc.map))){
    
#     nsc = uniq.nsc.map$NSC_id[i]
    
#     if(grepl("///" , nsc)){
#       all.nsc = unlist(strsplit(nsc , "///"))
#       print(all.nsc)
#       for(j in 1:length(all.nsc)){
#         uniq.nsc.map[1+nrow(uniq.nsc.map) , ] = uniq.nsc.map[i , ]
#         uniq.nsc.map$NSC_id[nrow(uniq.nsc.map)] = all.nsc[j]
#       }
#       uniq.nsc.map<- uniq.nsc.map[-i,]
#     }
#   }
# }

# colnames(uniq.nsc.map) [colnames(uniq.nsc.map) == "NSC_id"] <- "NSC" # To merge with info.sensitivity
# colnames(uniq.nsc.map) [colnames(uniq.nsc.map) == "unique.drugid"] <- "drugid" # To merge with info.sensitivity
# any(duplicated(uniq.nsc.map$NSC))#FALSE

# # Mapping NSC# to unique drug ids from uniq.nsc.map object
# info.sensitivity <- merge(info.sensitivity, uniq.nsc.map , by="NSC" , all.x=T)

# rownames(info.sensitivity) <- info.sensitivity$EXP_details
# info.sensitivity <- info.sensitivity[rownames(raw.sensitivity),]


# ########################################### PSet Curation ##########################################################
# NCI60_PSet<- PharmacoGx::PharmacoSet("NCI60_PSet",
#                                      molecularProfiles = list( "rna" = RNA_SE, "mirna" = miRNA_SE,
#                                                                "rnaseq.comp"= RNA_seq_comp_SE,
#                                                                "rnaseq.iso" = RNA_seq_iso_SE),
                                     
                                     
#                                      cell = cell.obj,
#                                      drug = drug.obj,
#                                      sensitivityInfo = info.sensitivity,
#                                      sensitivityRaw = raw.sensitivity,
#                                      sensitivityProfiles <- profile.sensitivity,
#                                      curationDrug = cur.drug,
#                                      curationCell = cur.cell,
#                                      curationTissue = cur.tissue,
#                                      datasetType = "sensitivity",
#                                      verify = TRUE)


# NCI60_PSet@annotation$notes <- "This PSet includes drug-dose information from NCI60 dataset. Molecular profile layers are processed data fetched directly from \"https://discover.nci.nih.gov/cellminer/\". Dose values in sensitivity data are reported in micromolar."  

# ########################################### Making the PSet compatible with PHarmacoDB ###########################################
# # Converting SummarizedExperiment to RangedSummarizedExperiments
# SE_list <- molecularProfilesSlot(NCI60_PSet)
# rowDataL <- lapply(SE_list, FUN=rowData)
# gRangesL <- lapply(rowDataL, FUN=makeGRangesFromDataFrame, 
#                    keep.extra=TRUE)
# SE_list <- Map(f=`rowRanges<-`, x=SE_list, value=gRangesL)
# molecularProfilesSlot(NCI60_PSet) <- SE_list

# ####### Mapping genomic co-ordinates from mol profiles to ensemble ids #######
# ## =====================
# orgDB <- org.Hs.eg.db
# for (i in seq_along(molecularProfilesSlot(NCI60_PSet))) {
#   SE <- molecularProfilesSlot(NCI60_PSet)[[i]]
#   rRanges <- rowRanges(SE)
#   ## TODO:: Implement mapping for mirna if possible?
#   # Skip mirna
#   if (metadata(SE)$annotation == 'mirna') next
#   # -- 6.2 Try look-up with Symbols
#   symbol <- as.character(rRanges$Gene.name)
#   # Entrez multimaps to Ensembl gene and trascript, try taking the first result
#   cols <- c('GENENAME', 'ENSEMBL')
#   res <- as.data.table(AnnotationDbi::select(orgDB, keys=symbol, columns=cols, keytype='SYMBOL'))
#   bySymbolDT <- res[, lapply(.SD, first), by=SYMBOL]
#   # Join and make sure they match
#   rangeMColDT <- as.data.table(mcols(rRanges))
#   mergeMColDT <- merge.data.table(rangeMColDT, bySymbolDT,
#                                   by.x='Gene.name', by.y='SYMBOL', all.x=TRUE)
#   # Reassign to mcols
#   mcols(rRanges) <- as(mergeMColDT, 'DataFrame')
#   # -- 6.3 Retry look-up with entrez IDs
#   mcolsDT <- as.data.table(mcols(rRanges))
#   entrez <- as.character(mcolsDT[is.na(ENSEMBL), Entrez.gene.id])
#   res1 <- as.data.table(AnnotationDbi::select(orgDB, keys=entrez, columns=c('ENSEMBL'), 
#                                               keytype='ENTREZID'))
#   moreMappingsDT <- res1[!is.na(ENSEMBL), lapply(.SD, first), by=ENTREZID]
#   moreMappingsDT[, ENTREZID := as.numeric(ENTREZID)]
#   # Do a join with update by reference
#   mcolsDT[moreMappingsDT, ENSEMBL := i.ENSEMBL, on=c("Entrez.gene.id==ENTREZID")]
#   # Reassign to mcols
#   mcols(rRanges) <- as(mcolsDT, 'DataFrame')
#   # -- 6.4 Map to ENSEMBL transcripts
#   ensembl <- na.omit(mcols(rRanges)$ENSEMBL)
#   tx_ids <- as.data.table(AnnotationDbi::select(orgDB, keys=ensembl, columns=c('ENSEMBLTRANS'), 
#                                                 keytype='ENSEMBL'))
#   txByEnsembl <- tx_ids[!is.na(ENSEMBLTRANS), lapply(.SD, paste, collapse='|'), 
#                         by=ENSEMBL]
#   mcolsDT <- merge.data.table(mcolsDT, txByEnsembl, by='ENSEMBL', all.x=TRUE)
#   # Adding gene_id column to the mol-profile data
#   setnames(mcolsDT,
#            old=c('ENSEMBL', 'ENSEMBLTRANS', 'Entrez.gene.id', 'Cytoband', 'Gene.name',
#                  'Gene.name.url', 'Entrez.gene.id.url', 'Genomic.coordinate.url',
#                  'GENENAME'),
#            new=c('gene_id', 'ensembl_tid', 'entrez_gid', 'cytoband', 'hugo_symbol', 
#                  'gene_name_url', 'entrez_gid_url', 'genomic_coord_url', 
#                  'gene_description'),
#            skip_absent=TRUE)
#   mcols(rRanges) <- as(mcolsDT, 'DataFrame')
#   # -- 6.5 Assign back to the PSet
#   rowRanges(SE) <- rRanges
#   molecularProfilesSlot(NCI60_PSet)[[i]] <- SE
# }


# # Proportions of NAs in mod profile data
# for (i in seq_along(molecularProfilesSlot(NCI60_PSet))) {
#   SE <- molecularProfilesSlot(NCI60_PSet)[[i]]
#   print(metadata(SE)$annotation)
#   print(paste(length(which(is.na(mcols(rowRanges(SE))$gene_id))) / length(rowRanges(SE)) , "NA Ensemb_ids", sep =" "))
# }

# # output
# saveRDS(NCI60_PSet, paste0(out_dir,"NCI60_PSet.rds"))

#output ORCESTRA_ID and Pachyderm commit id
# write.table(dataset, file="/pfs/out/dataset.txt", row.names = F ,quote = F, sep = "\t", col.names = F)
# write.table(ORCESTRA_ID, file="/pfs/out/orcestra_id.txt", row.names = F ,quote = F, sep = "\t", col.names = F)				   
# pach_commit_id <- Sys.getenv("PACH_OUTPUT_COMMIT_ID")
# write.table(pach_commit_id, file="/pfs/out/commit_id.txt", row.names = F ,quote = F, sep = "\t", col.names = F)