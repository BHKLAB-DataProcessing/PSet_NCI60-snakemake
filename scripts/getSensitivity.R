options(encoding = "UTF-8")
library(data.table)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
work_dir <- args[1]

annotation_dir <- paste0(work_dir, "annotation") # Stores annotation files used in the script. Kept up to date using Pachyderm's downAnnotations pipeline. 
sens_dir <- paste0(work_dir, "sensdata") # Stores sensitivity data files

drug_with_ids <- read.csv(file.path(annotation_dir, "drugs_with_ids.csv"), stringsAsFactors = F , na.strings = "") 
cell.obj.sen <- readRDS(paste0(work_dir, "curation/cell_obj_sen.rds"))

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
raw.sensitivity <-readRDS(file.path(sens_dir,"raw.sensitivity_v3.rds"))
raw.sensitivity[,,"Dose"] <- raw.sensitivity[,,"Dose"]* 1000000 # convert to micromolar

profile.sensitivity <-readRDS(file.path(sens_dir,"profile.sensitivity_v3.rds"))
profile.sensitivity <- profile.sensitivity[rownames(raw.sensitivity), ] # Rearranging the rownames 


dose.resp<-fread(file.path(sens_dir,"DOSERESP.csv"))
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

saveRDS(info.sensitivity, paste0(work_dir, "sensitivity/info_sensitivity.rds"))
saveRDS(profile.sensitivity, paste0(work_dir, "sensitivity/profile_sensitivity.rds"))
saveRDS(raw.sensitivity, paste0(work_dir, "sensitivity/raw_sensitivity.rds"))