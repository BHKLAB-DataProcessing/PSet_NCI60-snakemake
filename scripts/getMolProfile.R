options(encoding = "UTF-8")
library(readxl)
library(data.table)
library(stringr)
library(SummarizedExperiment)
library(GenomicRanges)

args <- commandArgs(trailingOnly = TRUE)
work_dir <- args[1]

annotation_dir <- paste0(work_dir, "annotation") # Stores annotation files used in the script. Kept up to date using Pachyderm's downAnnotations pipeline.
mol_dir <- paste0(work_dir, "moldata") # Stores molecular data files
cell_dir <- paste0(work_dir, "celldata") # Stores cell metadata

# Lab cell names
lab.cell.names <- read.csv(file.path(annotation_dir, "cell_annotation_all.csv") , na.strings = "") # this file is "cell_annotation_all.csv" from pachy annotation 

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
cell<-read.delim(file.path(cell_dir, "NCI60_CELL_LINE_METADATA.txt") , skip = 7, check.names =F)[c(1:60),]
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

rna <- read_excel(file.path(mol_dir, "RNA__Affy_HG_U133_Plus_2.0_RMA.xls") , skip =10)

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
rna.eSet<- ExpressionSet(assayData = data.matrix(assay.rna), phenoData = AnnotatedDataFrame(phen.rna), 
                         featureData = AnnotatedDataFrame(feat.rna)) 

############ RNA:Agilent Human microRNA(v2)  ############
# mirna.file <- tempfile()
# download.file("https://discover.nci.nih.gov/cellminerdata/normalizedArchives/nci60_RNA__Agilent_Human_microRNA_(V2)_GeneSpringGX.zip", mirna.file)

mirna <- read_excel(file.path(mol_dir, "RNA__Agilent_Human_microRNA_V2_GeneSpringGX.xls") ,skip = 10)

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

rnaseq.comp <- read_excel(file.path(mol_dir, "RNA__RNA_seq_composite_expression.xls") ,skip = 10)

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

rnaseq.iso <- read_excel(file.path(mol_dir, "RNA__RNA_seq_isoforms.xls") ,skip = 10)

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
saveRDS(eSetToSE(rna.eSet,annot_name="rna"), paste0(work_dir, "se/RNA_SE.rds"))
saveRDS(eSetToSE(mirna.eSet,annot_name="rna"), paste0(work_dir, "se/miRNA_SE.rds"))
saveRDS(eSetToSE(rnaseq.comp.eSet,annot_name="rnaseq"), paste0(work_dir, "se/RNA_seq_comp_SE.rds"))
saveRDS(eSetToSE(rnaseq.iso.eSet,annot_name="rnaseq"), paste0(work_dir, "se/RNA_seq_iso_SE.rds"))

#other data to be used in downstream processing
saveRDS(lab.cell.names, paste0(work_dir, "common/lab_cell_names.rds"))
saveRDS(phen.rna, paste0(work_dir, "common/phen_rna.rds"))
saveRDS(phen.mirna, paste0(work_dir, "common/phen_mirna.rds"))
saveRDS(phen.rnaseq.comp, paste0(work_dir, "common/phen_rnaseq_comp.rds"))
saveRDS(phen.rnaseq.iso, paste0(work_dir, "common/phen_rnaseq_iso.rds"))
