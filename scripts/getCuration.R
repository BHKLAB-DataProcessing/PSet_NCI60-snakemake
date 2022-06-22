options(encoding = "UTF-8")
library(data.table)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
work_dir <- args[1]

annotation_dir <- paste0(work_dir, "annotation") # Stores annotation files used in the script. Kept up to date using Pachyderm's downAnnotations pipeline. 
sens_dir <- paste0(work_dir, "sensdata") # Stores sensitivity data files
common_dir <- paste0(work_dir, "common") # Stores depedent data

# load dependent data
lab.cell.names <- readRDS(file.path(common_dir, "lab_cell_names.rds"))
phen.rna <- readRDS(file.path(common_dir, "phen_rna.rds"))
phen.mirna <- readRDS(file.path(common_dir, "phen_mirna.rds"))
phen.rnaseq.comp <- readRDS(file.path(common_dir, "phen_rnaseq_comp.rds"))
phen.rnaseq.iso <- readRDS(file.path(common_dir, "phen_rnaseq_iso.rds"))

########################################### Cell object ##########################################################
dose.resp<-fread(file.path(sens_dir,"DOSERESP.csv"))

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
drug_with_ids <- read.csv(file.path(annotation_dir,"drugs_with_ids.csv"), stringsAsFactors = F , na.strings = "") # this file is "drug_with_ids.csv" from pachy-annotations
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

saveRDS(cur.drug, paste0(work_dir, "curation/curation_drug.rds"))
saveRDS(cur.cell, paste0(work_dir, "curation/curation_cell.rds"))
saveRDS(cur.tissue, paste0(work_dir, "curation/curation_tissue.rds"))
saveRDS(cell.obj, paste0(work_dir, "curation/cell_obj.rds"))
saveRDS(drug.obj, paste0(work_dir, "curation/drug_obj.rds"))
saveRDS(cell.obj.sen, paste0(work_dir, "curation/cell_obj_sen.rds"))