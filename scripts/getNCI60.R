options(encoding = "UTF-8")
library(PharmacoGx)
library(readxl)
library(data.table)
library(stringr)
library(SummarizedExperiment)
library(GenomicRanges)
library(org.Hs.eg.db)

args <- commandArgs(trailingOnly = TRUE)
work_dir <- args[[1]]
filename <- args[[2]]

RNA_SE <- readRDS(paste0(work_dir, "se/RNA_SE.rds"))
miRNA_SE <- readRDS(paste0(work_dir, "se/miRNA_SE.rds"))
RNA_seq_comp_SE <- readRDS(paste0(work_dir, "se/RNA_seq_comp_SE.rds"))
RNA_seq_iso_SE <- readRDS(paste0(work_dir, "se/RNA_seq_iso_SE.rds"))
cell.obj <- readRDS(paste0(work_dir, "curation/cell_obj.rds"))
drug.obj <- readRDS(paste0(work_dir, "curation/drug_obj.rds"))
cur.drug <- readRDS(paste0(work_dir, "curation/curation_drug.rds"))
cur.cell <- readRDS(paste0(work_dir, "curation/curation_cell.rds"))
cur.tissue <- readRDS(paste0(work_dir, "curation/curation_tissue.rds"))
info.sensitivity <- readRDS(paste0(work_dir, "sensitivity/info_sensitivity.rds"))
raw.sensitivity <- readRDS(paste0(work_dir, "sensitivity/raw_sensitivity.rds"))
profile.sensitivity <- readRDS(paste0(work_dir, "sensitivity/profile_sensitivity.rds"))

########################################### PSet Curation ##########################################################
NCI60_PSet<- PharmacoGx::PharmacoSet("NCI60",
                                     molecularProfiles = list( "rna" = RNA_SE, "mirna" = miRNA_SE,
                                                               "rnaseq.comp"= RNA_seq_comp_SE,
                                                               "rnaseq.iso" = RNA_seq_iso_SE),
                                     
                                     
                                     cell = cell.obj,
                                     drug = drug.obj,
                                     sensitivityInfo = info.sensitivity,
                                     sensitivityRaw = raw.sensitivity,
                                     sensitivityProfiles = profile.sensitivity,
                                     curationDrug = cur.drug,
                                     curationCell = cur.cell,
                                     curationTissue = cur.tissue,
                                     datasetType = "sensitivity",
                                     verify = TRUE)


NCI60_PSet@annotation$notes <- "This PSet includes drug-dose information from NCI60 dataset. Molecular profile layers are processed data fetched directly from CellMiner website. Dose values in sensitivity objects are reported in micromolar. One Chinese hamster cell line (CHO) and two mouse cell lines (P388 and P388/ADR) are identified among the cell lines included in this PSet."
#qsave(NCI60_PSet, "~/nikta/out_nci60/NCI60_PSet.qs")
#NCI60_PSet <- qread("~/nikta/out_nci60/NCI60_PSet.qs")

########################################### Making the PSet compatible with PHarmacoDB ###########################################
# Converting SummarizedExperiment to RangedSummarizedExperiment
#NCI <- qread("~/nikta/out_nci60/NCI60_PSet.qs")
SE_list <- molecularProfilesSlot(NCI60_PSet)
rowDataL <- lapply(SE_list, FUN=rowData)
gRangesL <- lapply(rowDataL, FUN=makeGRangesFromDataFrame, 
                   keep.extra=TRUE)
SE_list <- Map(f=`rowRanges<-`, x=SE_list, value=gRangesL)
molecularProfilesSlot(NCI60_PSet) <- SE_list

####### Mapping genomic co-ordinates from mol profiles to ensemble ids #######
#### Chris's code :

## =====================
## ---- 0. Load the data
#NCI <- qread('NCI60_PSet_Ranges.qs')
## ====================
## ---- 6. Using OrgDB
# -- 6.1 Load the database package
orgDB <- org.Hs.eg.db
for (i in seq_along(molecularProfilesSlot(NCI60_PSet))) {
  SE <- molecularProfilesSlot(NCI60_PSet)[[i]]
  rRanges <- rowRanges(SE)
  ## TODO:: Implement mapping for mirna if possible?
  # Skip mirna
  if (metadata(SE)$annotation == 'mirna') next
  # -- 6.2 Try look-up with Symbols
  symbol <- as.character(rRanges$Gene.name)
  # Entrez multimaps to Ensembl gene and trascript, try taking the first result
  cols <- c('GENENAME', 'ENSEMBL')
  res <- as.data.table(AnnotationDbi::select(orgDB, keys=symbol, columns=cols, keytype='SYMBOL'))
  bySymbolDT <- res[, lapply(.SD, first), by=SYMBOL]
  # Join and make sure they match
  rangeMColDT <- as.data.table(mcols(rRanges))
  mergeMColDT <- merge.data.table(rangeMColDT, bySymbolDT,
                                  by.x='Gene.name', by.y='SYMBOL', all.x=TRUE, sort=FALSE)
  # Reassign to mcols
  mcols(rRanges) <- as(mergeMColDT, 'DataFrame')
  # -- 6.3 Retry look-up with entrez IDs
  mcolsDT <- as.data.table(mcols(rRanges))
  entrez <- as.character(mcolsDT[is.na(ENSEMBL), Entrez.gene.id])
  res1 <- as.data.table(AnnotationDbi::select(orgDB, keys=entrez, columns=c('ENSEMBL'), 
                                              keytype='ENTREZID'))
  moreMappingsDT <- res1[!is.na(ENSEMBL), lapply(.SD, first), by=ENTREZID]
  moreMappingsDT[, ENTREZID := as.numeric(ENTREZID)]
  # Do a join with update by reference
  mcolsDT[moreMappingsDT, ENSEMBL := i.ENSEMBL, on=c("Entrez.gene.id==ENTREZID")]
  # Reassign to mcols
  mcols(rRanges) <- as(mcolsDT, 'DataFrame')
  # -- 6.4 Map to ENSEMBL transcripts
  ensembl <- na.omit(mcols(rRanges)$ENSEMBL)
  tx_ids <- as.data.table(AnnotationDbi::select(orgDB, keys=ensembl, columns=c('ENSEMBLTRANS'), 
                                                keytype='ENSEMBL'))
  txByEnsembl <- tx_ids[!is.na(ENSEMBLTRANS), lapply(.SD, paste, collapse='|'), 
                        by=ENSEMBL]
  mcolsDT <- merge.data.table(mcolsDT, txByEnsembl, by='ENSEMBL', all.x=TRUE, sort=FALSE)
  # Adding gene_id column to the mol-profile data
  setnames(mcolsDT,
           old=c('ENSEMBL', 'ENSEMBLTRANS', 'Entrez.gene.id', 'Cytoband', 'Gene.name',
                 'Gene.name.url', 'Entrez.gene.id.url', 'Genomic.coordinate.url',
                 'GENENAME'),
           new=c('gene_id', 'ensembl_tid', 'entrez_gid', 'cytoband', 'hugo_symbol', 
                 'gene_name_url', 'entrez_gid_url', 'genomic_coord_url', 
                 'gene_description'),
           skip_absent=TRUE)
  mcols(rRanges) <- as(mcolsDT, 'DataFrame')
  # -- 6.5 Assign back to the PSet
  rowRanges(SE) <- rRanges
  molecularProfilesSlot(NCI60_PSet)[[i]] <- SE
}


# Getting proportions of NAs in mod profile data
for (i in seq_along(molecularProfilesSlot(NCI60_PSet))) {
  SE <- molecularProfilesSlot(NCI60_PSet)[[i]]
  print(metadata(SE)$annotation)
  print(paste(length(which(is.na(mcols(rowRanges(SE))$gene_id))) / length(rowRanges(SE)) , "NA Ensemb_ids", sep =" "))
}


saveRDS(NCI60_PSet, paste0(work_dir, filename))
