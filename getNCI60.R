options(encoding = "UTF-8")
# work_dir = "~/nikta/out_nci60/"
# setwd(work_dir)
library(PharmacoGx)
library(readxl)
library(data.table)
library(stringr)
library(SummarizedExperiment)
library(GenomicRanges)
library(qs)
library(org.Hs.eg.db)

annotation_dir <- "/pfs/downAnnotations/" # Stores annotation files used in the script. Kept up to date using Pachyderm's downAnnotations pipeline. 
sens_dir <- "/pfs/downloadNCI60SensData/" # Stores sensitivity data files
mol_dir <- "/pfs/downloadNCI60MolData/" # Stores molecular data files
cell_dir <- "/pfs/downloadNCI60CellData/" # Stores cell metadata
out_dir <- "/pfs/out/" # Stores the end product (NCI60 PSet).

# Lab cell names
lab.cell.names <- read.csv(paste0(annotation_dir, "cell_annotation_all.csv") , na.strings = "") # this file is "cell_annotation_all.csv" from pachy annotation 

# Removing "///" from nci60.cellid for merge
while(any(grepl("///" , lab.cell.names$NCI60.cellid))){
  for (i in seq(nrow(lab.cell.names))){
    name = lab.cell.names$NCI60.cellid[i]
    if(grepl("///" , name)){
      
      all.name = unlist(strsplit(name, "///"))
      print(all.name)
      for(j in 1:length(all.name)){
        lab.cell.names[1+nrow(lab.cell.names) , ] = lab.cell.names[i , ]
        lab.cell.names$NCI60.cellid[nrow(lab.cell.names)] = all.name[j]
      }
      lab.cell.names <- lab.cell.names[-i,]
    }
  }
}


### Cell-line information
cell<-read.delim(paste0(cell_dir, "NCI60_CELL_LINE_METADATA.txt") , skip = 7, check.names =F)[c(1:60),]
#colnames(cell) [colnames(cell) =="tissue of origin (a)"] <- "phen_tissue"

######################################################### pheno.data function ######################################################### 
phen.func <- function(assay , cell.published , lab.names) {
  
  phen <- data.frame(cell.line = colnames(assay)) #cell.line column is the published names
  phen $ modified.name <-  sub("." , ":" , phen $ cell.line   , fixed=T)
  phen $ modified.name <- gsub("." , "_" , phen $ modified.name, fixed = T)
  
  print(all(phen$cell.line %in% lab.cell.names$NCI60.cellid))
  
  # For merge to published cell info:
  phen $modified.name [phen $modified.name == "BR:HS_578T"] <- "BR:HS578T"
  phen $modified.name [phen $modified.name == "BR:T_47D"] <- "BR:T47D"
  phen $modified.name [phen $modified.name == "CO:COLO_205"] <- "CO:COLO205"
  phen $modified.name [phen $modified.name == "LC:A549_ATCC"] <- "LC:A549" 
  phen $modified.name [phen $modified.name == "ME:LOX_IMVI"] <- "ME:LOXIMVI"
  
  # Adding unique names from the lab names
  temp <- merge(lab.cell.names[ , c("unique.cellid" , "NCI60.cellid", "unique.tissueid")] , phen , by.x="NCI60.cellid" , by.y= "cell.line" , all.y=T )
  colnames(temp) <- c("NCI60.cellid","cellid", "tissueid","modified.name" )
  phen <- merge(temp , cell.published , by.x="modified.name" , by.y ="Cell Line Name" , all.x =T)[,-1] # removing first column which is "modified.name"
  rownames(phen) <- phen$cellid 
  rm(temp)
  
  phen$batchid <- NA
  return(phen)
}

######################################################### eSetToSE ######################################################
# A function converting ExpressionSet to SummarizedExperiment

eSetToSE <- function(eSet , annot_name) {
  
  BiocGenerics::annotation(eSet) <- annot_name
  stopifnot(all(rownames(fData(eSet)) == rownames(exprs(eSet))))
  stopifnot(all(rownames(pData(eSet)) == colnames(exprs(eSet))))
  
  # Build summarized experiment from eSet
  SE <- SummarizedExperiment::SummarizedExperiment(
    assays=SimpleList(as.list(Biobase::assayData(eSet))),
    # Switch rearrange columns so that IDs are first, probes second
    rowData=S4Vectors::DataFrame(Biobase::fData(eSet)),
    colData=S4Vectors::DataFrame(Biobase::pData(eSet)),
    metadata=list("experimentData" = eSet@experimentData, 
                  "annotation" = Biobase::annotation(eSet), 
                  "protocolData" = Biobase::protocolData(eSet)))
  # Extract names from expression set                  
  SummarizedExperiment::assayNames(SE) <- Biobase::assayDataElementNames(eSet)
  
  stopifnot(all(rownames(colData(SE)) == rownames(pData(eSet))))
  stopifnot(all(rownames(rowData(SE)) == rownames(fData(eSet))))
  stopifnot(all(rownames(colData(SE)) == colnames(assay(SE))))
  
  return(SE)
}

########################################### Molecular Profiles ##########################################################
############ RNA:Affy HG-U133 Plus 2.0 - RMA ############
# rna.file <- tempfile()
# download.file("https://discover.nci.nih.gov/cellminerdata/normalizedArchives/nci60_RNA__Affy_HG_U133_Plus_2.0_RMA.zip", rna.file)

rna <- read_excel(paste0(mol_dir, "RNA__Affy_HG_U133_Plus_2.0_RMA.xls") , skip =10)

# Assay data rna
assay.rna<- data.frame(rna[,c(8:ncol(rna))] , row.names= rna$"Identifier c") 

# Feature data rnaseq
feat.rna <- data.frame(rna[,c(1:7)] , row.names = rna$"Identifier c")
colnames(feat.rna) <- str_sub( colnames(feat.rna) , end=-3)
feat.rna$Gene.name.url <- "http://www.genenames.org/"
feat.rna$Entrez.gene.id.url <- "http://www.ncbi.nlm.nih.gov/gene"
feat.rna$Genomic.coordinate.url <- "https://www.affymetrix.com/analysis/netaffx/xmlquery.affx?netaffx=mapping"
feat.rna$Probe.info.url <- "http://www.affymetrix.com/estore/esearch/search.jsp?pd=131459&N=4294967292"
feat.rna$BEST <- NA
feat.rna$Symbol <- feat.rna$Gene.name
feat.rna$Entrez.gene.id <- ifelse(feat.rna$Entrez.gene.id==0, NA, feat.rna$Entrez.gene.id)

# Pheno data rna
phen.rna <- phen.func(assay = assay.rna , cell.published =cell, lab.names = lab.cell.names)  

# ESet rna
assay.rna<- assay.rna[, phen.rna$NCI60.cellid] # making sure the order is the same with pheno data
colnames(assay.rna) <- phen.rna$cellid
assay.rna <-assay.rna [rownames(feat.rna),rownames(phen.rna)] # rearrangement
rna.eSet<- ExpressionSet(assayData = as.matrix(assay.rna), phenoData = AnnotatedDataFrame(phen.rna), 
                         featureData = AnnotatedDataFrame(feat.rna)) 

############ RNA:Agilent Human microRNA(v2)  ############
# mirna.file <- tempfile()
# download.file("https://discover.nci.nih.gov/cellminerdata/normalizedArchives/nci60_RNA__Agilent_Human_microRNA_(V2)_GeneSpringGX.zip", mirna.file)

mirna <- read_excel(paste0(mol_dir, "RNA__Agilent_Human_microRNA_(V2)_GeneSpringGX.xls") ,skip = 10)

# Assay data mirna
assay.mirna <- data.frame(mirna[,c(10:ncol(mirna))] , row.names= mirna$`Identifier c`) 

# Feature data mirna
feat.mirna <- data.frame(mirna[,c(1:6)] , row.names = mirna$`Identifier c`)
colnames(feat.mirna) <- str_sub( colnames(feat.mirna) , end=-3)
feat.mirna$Gene.name.url <- "http://www.genenames.org/"
feat.mirna$Entrez.gene.id.url <- "http://www.ncbi.nlm.nih.gov/gene"
feat.mirna$Genomic.coordinate.url <- "https://www.ncbi.nlm.nih.gov/gene?Db=gene&Cmd=ShowDetailView"
feat.mirna$miRNA.accession.url <- "http://www.mirbase.org/help/nomenclature.shtml"
feat.mirna$MirBase.name.url <- "http://www.mirbase.org/help/nomenclature.shtml"
feat.mirna$gene_id <-rownames(feat.mirna)
feat.mirna$Entrez.gene.id <- ifelse(feat.mirna$Entrez.gene.id==0, NA, feat.mirna$Entrez.gene.id)

# Pheno data mirna
phen.mirna <- phen.func(assay = assay.mirna , cell.published =cell, lab.names = lab.cell.names)  

# ESet mirna
assay.mirna<- assay.mirna[, phen.mirna$NCI60.cellid] # making sure the order is the same with pheno data
colnames(assay.mirna) <- phen.mirna$cellid
assay.mirna <-assay.mirna [rownames(feat.mirna),rownames(phen.mirna)] # rearrangement
mirna.eSet<- ExpressionSet(assayData = as.matrix(assay.mirna), phenoData = AnnotatedDataFrame(phen.mirna), 
                           featureData = AnnotatedDataFrame(feat.mirna)) 

############ RNA:RNA_Seq_composite ############
# rnaseq.comp.file <- tempfile()
# download.file("https://discover.nci.nih.gov/cellminerdata/normalizedArchives/nci60_RNA__RNA_seq_composite_expression.zip", rnaseq.comp.file)

rnaseq.comp <- read_excel(paste0(mol_dir, "RNA__RNA_seq_composite_expression.xls") ,skip = 10)

# Assay data rnaseq.comp
assay.rnaseq.comp <- data.frame(rnaseq.comp[,c(7:ncol(rnaseq.comp))] , row.names= rnaseq.comp$`Gene name d`) 

# Feature data rnaseq.comp
feat.rnaseq.comp <- data.frame(rnaseq.comp[,c(1:6)] , row.names = rnaseq.comp$`Gene name d`)
colnames(feat.rnaseq.comp) <- str_sub( colnames(feat.rnaseq.comp) , end=-3)
feat.rnaseq.comp$Gene.name.url <- "http://www.genenames.org/"
feat.rnaseq.comp$Entrez.gene.id.url <- "http://www.ncbi.nlm.nih.gov/gene"
feat.rnaseq.comp$Genomic.coordinate.url <- "https://www.ncbi.nlm.nih.gov/gene"
feat.rnaseq.comp$Entrez.gene.id <- ifelse(feat.rnaseq.comp$Entrez.gene.id==0, NA, feat.rnaseq.comp$Entrez.gene.id)

# Pheno data rnaseq.comp
phen.rnaseq.comp <- phen.func(assay = assay.rnaseq.comp , cell.published =cell, lab.names = lab.cell.names)  

# ESet rnaseq.comp
assay.rnaseq.comp <- assay.rnaseq.comp[, phen.rnaseq.comp$NCI60.cellid] # making sure the order is the same with pheno data
colnames(assay.rnaseq.comp) <- phen.rnaseq.comp$cellid
assay.rnaseq.comp <-assay.rnaseq.comp [rownames(feat.rnaseq.comp),rownames(phen.rnaseq.comp)] # rearrangement
rnaseq.comp.eSet<- ExpressionSet(assayData = as.matrix(assay.rnaseq.comp), phenoData = AnnotatedDataFrame(phen.rnaseq.comp), 
                                 featureData = AnnotatedDataFrame(feat.rnaseq.comp)) 

############ RNA:RNA_Seq_isoform ############
# rnaseq.iso.file <- tempfile()
# download.file("https://discover.nci.nih.gov/cellminerdata/normalizedArchives/nci60_RNA__RNA_seq_isoforms.zip", rnaseq.iso.file)

rnaseq.iso <- read_excel(paste0(mol_dir, "RNA__RNA_seq_isoforms.xls") ,skip = 10)

# Assay data rnaseq.iso
assay.rnaseq.iso <- data.frame(rnaseq.iso[,c(9:ncol(rnaseq.iso))] , row.names= rnaseq.iso$`Identifier c`) 

# Feature data rnaseq.iso
feat.rnaseq.iso <- data.frame(rnaseq.iso[,c(1:8)] , row.names = rnaseq.iso$`Identifier c`)
colnames(feat.rnaseq.iso) <- str_sub( colnames(feat.rnaseq.iso) , end=-3)
colnames(feat.rnaseq.iso)[colnames(feat.rnaseq.iso) == "RefSeq..protein."] <- "RefSeq..protein"
feat.rnaseq.iso$Gene.name.url <- "http://www.genenames.org/"
feat.rnaseq.iso$Entrez.gene.id.url <- "http://www.ncbi.nlm.nih.gov/gene"
feat.rnaseq.iso$Genomic.coordinate.url <- "https://www.ncbi.nlm.nih.gov/gene"
feat.rnaseq.iso$Transcript.reference.url <- "http://www.ncbi.nlm.nih.gov/nucleotide"
feat.rnaseq.iso$Entrez.gene.id <- ifelse(feat.rnaseq.iso$Entrez.gene.id==0, NA, feat.rnaseq.iso$Entrez.gene.id)

# Pheno data rnaseq.iso
phen.rnaseq.iso <- phen.func(assay = assay.rnaseq.iso , cell.published =cell, lab.names = lab.cell.names)  

# ESet rnaseq.iso
assay.rnaseq.iso <- assay.rnaseq.iso[, phen.rnaseq.iso$NCI60.cellid] # making sure the order is the same with pheno data
colnames(assay.rnaseq.iso) <- phen.rnaseq.iso$cellid
assay.rnaseq.iso <-assay.rnaseq.iso [rownames(feat.rnaseq.iso),rownames(phen.rnaseq.iso)] # rearrangement
rnaseq.iso.eSet<- ExpressionSet(assayData = as.matrix(assay.rnaseq.iso), phenoData = AnnotatedDataFrame(phen.rnaseq.iso), 
                                featureData = AnnotatedDataFrame(feat.rnaseq.iso)) 

########################################### SE objects ##########################################################
#Checks are included in the eSetToSE function
RNA_SE <- eSetToSE(rna.eSet,annot_name="rna")
miRNA_SE <- eSetToSE(mirna.eSet,annot_name="mirna")
RNA_seq_comp_SE <- eSetToSE(rnaseq.comp.eSet,annot_name="rnaseq.comp")
RNA_seq_iso_SE <- eSetToSE(rnaseq.comp.eSet,annot_name="rnaseq.iso")

########################################### Cell object ##########################################################
dose.resp<-fread(paste0(sens_dir,"DOSERESP.csv"))

# cell-line information
cols <- c( "RELEASE_DATE","PREFIX","PANEL_NUMBER" , "CELL_NUMBER","PANEL_NAME","CELL_NAME","PANEL_CODE" )
dose.resp <- data.frame(dose.resp[!duplicated(CELL_NAME),.SD,.SDcols = cols])                             
colnames(dose.resp)[colnames(dose.resp) == "CELL_NAME"] <- "NCI60.cellid"

all(dose.resp$NCI60.cellid %in% lab.cell.names$NCI60.cellid)#TRUE

# Cell object based on dose.resp cell names:
#lab.cell.names$CELL_NAME <- lab.cell.names$NCI60.cellid
cell.obj.sen <- merge(lab.cell.names[,c("unique.cellid","unique.tissueid","NCI60.cellid", "NCI60.tissueid")], dose.resp, by= "NCI60.cellid", all.y=T)
colnames(cell.obj.sen)[colnames(cell.obj.sen) == "unique.tissueid"] <-"tissueid"
colnames(cell.obj.sen)[colnames(cell.obj.sen) == "unique.cellid"] <-"cellid"

# Adding NCI60.cellids from pheno data
all.phen <- as.data.frame(unique(rbindlist( list(phen.rna, phen.mirna, phen.rnaseq.comp, phen.rnaseq.iso))))
#colnames(all.phen)[colnames(all.phen) == "NCI60.tissueid"] <-"phen_tissue"

# Creating cell object
cell.obj <- merge(cell.obj.sen  , all.phen, by=c("cellid","tissueid","NCI60.cellid"), all=T)

for( cell in unique(cell.obj$cellid)){
  ind = which(cell.obj$cellid == cell)
  if (length(ind)>1){
    print(ind)
    for(c in colnames(cell.obj)){
      if (c == "cellid"){next}
      id = paste(na.omit(unique(cell.obj[ind, c])) , collapse="///")
      cell.obj[ind, c] = id
    }
  }
}

cell.obj[cell.obj == ""] <- NA
cell.obj <- unique(cell.obj)
cell.obj <- cell.obj[,colnames(cell.obj)[!colnames(cell.obj)=="tissue of origin (a)"]]

all(all.phen$cellid %in% cell.obj$cellid) #TRUE
all(cell.obj.sen$cellid %in% cell.obj$cellid) #TRUE
rownames(cell.obj) <- cell.obj$cellid

# removing (a), (b), ... from colnames
colnames(cell.obj) <- gsub("\\s*\\([^\\)]+\\)","", colnames(cell.obj))


########################################### Drug object ##########################################################
drug_with_ids <- read.csv(paste0(annotation_dir,"drugs_with_ids.csv"), stringsAsFactors = F , na.strings = "") # this file is "drug_with_ids.csv" from pachy-annotations
drug.obj <- drug_with_ids[ !is.na(drug_with_ids$NCI60.drugid), c("unique.drugid" , "NCI60.drugid", "NSC_number", "cid","smiles","inchikey")]
colnames(drug.obj)[colnames(drug.obj) =="unique.drugid"]<-"drugid"
rownames(drug.obj) <- drug.obj$drugid

########################################### Curation ##########################################################
######### Drug #########
cur.drug <- data.frame(unique.drugid = drug.obj$drugid,
                       NCI60.drugid = drug.obj$NCI60.drugid,
                       row.names = drug.obj$drugid)

######### Cell #########
cur.cell <- data.frame(unique.cellid = cell.obj$cellid,
                       NCI60.cellid = cell.obj$NCI60.cellid,
                       row.names= cell.obj$cellid)

######### Tissue #########
cur.tissue <- data.frame(unique.tissueid = cell.obj$tissueid, 
                         NCI60.tissueid = cell.obj$NCI60.tissueid , 
                         row.names = cell.obj$cellid)

########################################### Sensitivity ##########################################################
# ##### h4h job(1) #### 
# library(data.table)
# dose.resp<-fread("DOSERESP.csv")
# # length(unique(dose.resp $"NSC")) #55601
# 
# 
# # Defining unique exp_detail column
# # EXP_details = NSC + CELL_NAME + PANEL_NAME + EXPID 
# dose.resp $ PANEL_NAME  <- sub(" .*", "", dose.resp $ PANEL_NAME) # Removing "cancer" from the names
# dose.resp $ EXP_details <- paste(dose.resp$CELL_NAME, dose.resp$PANEL_NAME , sep="_")
# dose.resp $ EXP_details <- paste(dose.resp $NSC, dose.resp $EXP_details , sep ="_")
# dose.resp $ EXP_details <- paste(dose.resp $EXP_details, dose.resp $EXPID , sep ="_")
# 
# length(unique(dose.resp$ EXP_details))#4568480
# 
# ## Will keep only data with 4+ points
# sub.dose.resp <- sub.dose.resp[!is.na(AVERAGE_PTC)]
# sub.dose.resp <- merge(sub.dose.resp[,.N, EXP_details][N>3], sub.dose.resp)
# raw.sensitivity <- array(data = NA_real_,
#                          dim= c(length(unique(sub.dose.resp$EXP_details)),max(sub.dose.resp$N), 2))
# 
# colnames(raw.sensitivity) <- paste0('dose:', seq_len(NCOL(raw.sensitivity)))
# rownames(raw.sensitivity) <- unique(sub.dose.resp$EXP_details)
# dimnames(raw.sensitivity)[[3]] <- c("Dose", "Viability") 
# setkey(sub.dose.resp, EXP_details) # defining a"key' column for sorting
# setorder(sub.dose.resp, EXP_details, CONCENTRATION)
# 
# for(exp in unique(sub.dose.resp[,EXP_details])){
#   raw.sensitivity[exp, seq_len(sub.dose.resp[EXP_details==exp, unique(N)]),] <- data.matrix(sub.dose.resp[EXP_details==exp, .(CONCENTRATION, AVERAGE_PTC)])
# }
# 
# raw.sensitivity[,,"Dose"] <- 10^raw.sensitivity[,,"Dose"]
# saveRDS(raw.sensitivity , "raw.sensitivity_v3.rds") # raw.sensitivity_v3 is without any cell or drug names mapping
# 
#### data-slicing - conducted on h4h4 ####
# raw.sensitivity <- readRDS("/cluster/home/nfeizi/raw.sensitivity_v3.rds") 
# raw.sensitivity[,,"Dose"] <- raw.sensitivity[,,"Dose"]* 1000000 #Convert to micromolar
# dir.create("/cluster/projects/bhklab/Data/NCI60/outmicro/")
# dir.create("/cluster/projects/bhklab/Data/NCI60/outmicro/slices/")
# 
# setwd("/cluster/projects/bhklab/Data/NCI60/outmicro/")
# raw.sensitivity.x <- parallel::splitIndices(nrow(raw.sensitivity), floor(nrow(raw.sensitivity)/100000))
# for(i in seq_along(raw.sensitivity.x)){
#   slce <- raw.sensitivity[raw.sensitivity.x[[i]],,]
#   saveRDS(slce, file=paste0("slices/NCI60_raw_sens_", i, ".rds"))
# }
# 
# #### h4h job(2) ####
# library(PharmacoGx)
# args <- commandArgs(trailingOnly = TRUE)
# jobIndex <- as.numeric(args[[1]])
# myfn <- list.files("/cluster/projects/bhklab/Data/NCI60/outmicro/slices/", full.names=TRUE)[[jobIndex]]
# print(myfn)
# mybasenm <- basename(myfn)
# slice <- readRDS(myfn)
# res <- PharmacoGx:::.calculateFromRaw(slice)
# saveRDS(res, file=paste0("/cluster/projects/bhklab/Data/NCI60/outmicro/", gsub(mybasenm, pattern = ".rds", replacement="_recomp.rds", fixed=TRUE)))

#### output from job2 ####
# myfn <- list.files("~/nikta/out_nci60/outmicro", full.names=TRUE, pattern = ".rds")
# slices <- list()
# 
# for(fn in myfn){
#   temp <- readRDS(fn)
#   parTable <- do.call(rbind,temp[[3]])
#   # print(head(rownames(parTable)))
#   # print(str(temp[[3]]))
#   n <- cbind("aac_recomputed" = as.numeric(unlist(temp[[1]]))/100,
#              "ic50_recomputed" = as.numeric(unlist(temp[[2]])),
#              "HS" = as.numeric(unlist(parTable[,1])),
#              "E_inf" = as.numeric(unlist(parTable[,2])),
#              "EC50" = as.numeric(unlist(parTable[,3])))
#   print(head(rownames(n)))
#   rownames(n) <- names(temp[[3]])
#   slices[[fn]] <- n
# }
# 
# profile.sensitivity <- as.data.frame(do.call(rbind, slices))
# saveRDS(profile.sensitivity , "~/nikta/out_nci60/profile.sensitivity_v3.rds")

##################### Sensitivity - all #####################
raw.sensitivity <-readRDS(paste0(sens_dir,"raw.sensitivity_v3.rds"))
raw.sensitivity[,,"Dose"] <- raw.sensitivity[,,"Dose"]* 1000000 # convert to micromolar

profile.sensitivity <-readRDS(paste0(sens_dir,"profile.sensitivity_v3.rds"))
profile.sensitivity <- profile.sensitivity[rownames(raw.sensitivity), ] # Rearranging the rownames 


dose.resp<-fread(paste0(sens_dir,"DOSERESP.csv"))
# length(unique(dose.resp $"NSC")) #55601

# Defining unique exp_detail column
# EXP_details = NSC + CELL_NAME + PANEL_NAME + EXPID 
dose.resp $ PANEL_NAME  <- sub(" .*", "", dose.resp $ PANEL_NAME) # Removing "cancer" from the names
dose.resp $ EXP_details <- paste(dose.resp$CELL_NAME, dose.resp$PANEL_NAME , sep="_") 
dose.resp $ EXP_details <- paste(dose.resp $NSC, dose.resp $EXP_details , sep ="_")
dose.resp $ EXP_details <- paste(dose.resp $EXP_details, dose.resp $EXPID , sep ="_")

# Sensitivity info object
cols <- c("EXP_details" , "NSC" , "CELL_NAME" )
info.sensitivity <- data.frame(dose.resp[!duplicated(EXP_details),.SD,.SDcols = cols])                             
rm(dose.resp)

# Because we only kept >3 points
info.sensitivity <- info.sensitivity[info.sensitivity$EXP_details %in% rownames(raw.sensitivity) , ] 
colnames(info.sensitivity)[colnames(info.sensitivity) == "CELL_NAME"] <- "NCI60.cellid"
#saveRDS(info.sensitivity , "~/nikta/out_nci60/info.sensitivity_v3.rds")# drugs and cells from info.sensitivity will be used for creating drug and cell objects

# Mapping cellids to unique ids from cell.obj.sen
info.sensitivity <- merge(info.sensitivity , cell.obj.sen[, c("cellid" ,"tissueid", "NCI60.cellid")], by= "NCI60.cellid", all.x=T) # Cellid column is from cell.obj so contains unique.ids

#length(unique(info.sensitivity$cellid)) #162
#length(unique(info.sensitivity$tissueid)) #18 (one is NA)
#any(is.na(info.sensitivity$cellid)) #FALSE


# Mapping drugids back to UNcontatenated NSC id to be able to map to senesitivity objects
uniq.nsc.map <- drug_with_ids[ !is.na(drug_with_ids$NSC_number), c("unique.drugid" , "NSC_number", "cid" )]

while (length(setdiff(unique(info.sensitivity$NSC) , uniq.nsc.map$NSC_number))>1){
  
  for (i in seq(nrow(uniq.nsc.map))){
    
    nsc = uniq.nsc.map$NSC_number[i]
    
    if(grepl("///" , nsc)){
      all.nsc = unlist(strsplit(nsc , "///"))
      print(all.nsc)
      for(j in 1:length(all.nsc)){
        uniq.nsc.map[1+nrow(uniq.nsc.map) , ] = uniq.nsc.map[i , ]
        uniq.nsc.map$NSC_number[nrow(uniq.nsc.map)] = all.nsc[j]
      }
      uniq.nsc.map<- uniq.nsc.map[-i,]
    }
  }
}

colnames(uniq.nsc.map) [colnames(uniq.nsc.map) == "NSC_number"] <- "NSC" # To merge with info.sensitivity
colnames(uniq.nsc.map) [colnames(uniq.nsc.map) == "unique.drugid"] <- "drugid" # To merge with info.sensitivity
any(duplicated(uniq.nsc.map$NSC))#FALSE

# Mapping NSC# to unique drug ids from uniq.nsc.map object
info.sensitivity <- merge(info.sensitivity, uniq.nsc.map , by="NSC" , all.x=T)
#any(is.na(info.sensitivity$NSC)) #FALSE

rownames(info.sensitivity) <- info.sensitivity$EXP_details
info.sensitivity <- info.sensitivity[rownames(raw.sensitivity), c("EXP_details", "cellid" , "tissueid" , "drugid" , "NSC",  "NCI60.cellid")]

########################################### PSet Curation ##########################################################
NCI60_PSet<- PharmacoGx::PharmacoSet("NCI60",
                                     molecularProfiles = list( "rna" = RNA_SE, "mirna" = miRNA_SE,
                                                               "rnaseq.comp"= RNA_seq_comp_SE,
                                                               "rnaseq.iso" = RNA_seq_iso_SE),
                                     
                                     
                                     cell = cell.obj,
                                     drug = drug.obj,
                                     sensitivityInfo = info.sensitivity,
                                     sensitivityRaw = raw.sensitivity,
                                     sensitivityProfiles <- profile.sensitivity,
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
                                  by.x='Gene.name', by.y='SYMBOL', all.x=TRUE)
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
  mcolsDT <- merge.data.table(mcolsDT, txByEnsembl, by='ENSEMBL', all.x=TRUE)
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


saveRDS(NCI60_PSet, paste0(out_dir,"NCI60.rds"))