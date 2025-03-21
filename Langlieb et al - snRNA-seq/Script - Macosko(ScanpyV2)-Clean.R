#Only use first build of environment
install.packages("renv")
renv::init()

#Pick the reticulate directory
renv::use_python()

#Install packages into renv
renv::install("anndata")
renv::install("reticulate")
renv::install("data.table")
renv::install("dplyr")
renv::install("xlsx")
renv::install("tidyr")

reticulate::py_install("scanpy")

library("anndata")
library("data.table")
library("dplyr")
library("tidyr")
library("xlsx")

renv::snapshot()

sc <- reticulate::import("scanpy")

#Changes R default programming to not use scientific notation, can ruin cell_labels
#Can reset it by setting it to 0
options(scipen=999)

adfile <- read_h5ad("Macosko_Mouse_Atlas_Single_Nuclei.Use_Backed.h5ad", backed = "r")

#Associate ensemble gene names with genes of interest
dopagenes <- data.frame(genenum = c("ENSMUSG00000021478", "ENSMUSG00000032259","ENSMUSG00000024397"),
                        genesym = c("DRD1","DRD2","IBA1"))

#Isolate for dopamine receptor genes and IBA1
adfiledopa <- sc$get$obs_df(adfile,
                            keys = dopagenes$genenum)

#Convert to data.table
adfiledopa_Table <- data.table(adfiledopa,keep.rownames=TRUE)

#Split off donor library
adfiledopa_Table <- adfiledopa_Table %>% separate(rn, c('rn', 'cell'))

#Get cell metadata
metadata <- fread("Library_Metadata.tsv")

#remove unnecessary columns
colnames(metadata)
metadata <- metadata[,!c("dissectate","donor_id","sex")]

#-------------------------------------------------------------------------------
#STRIATUM (No LSX) ANALYSIS
#-------------------------------------------------------------------------------

#Remove all cells not in the striatum
unique(metadata$brain_struct)
#Possible brain structs
#"Isocortex" "CTXsp"     "CB"        "HPF"       "STR"       "PAL"       "HY"        "P"         "MB"       
#"MY"        "OLF"       "TH"    

unique(metadata$region)
#Possible brain regions
#[1] "MOp"  "RSP"  "ACA"  "AMY"  "CB"   "CTX"  "AUD"  "HPF"  "STRd" "DCN"  "ENT"  "PALd" "HY"   "BS"   "MB"   "LSX"  "SUB" 
#[18] "PALm" "NTS"  "OLF"  "S1"   "PALv" "TH"   "VIS"  "VISP" "BNST"

unique(metadata$sub_region)
#Possible brain sub regions
#[1] "A"       "ACA"     "AMY"     "AN1"     "AN2"     "ANG"     "AUD"     "CA1"     "CA2"     "CA2CA3"  "CA3"    
#[12] "COP"     "CP"      "CUL"     "DCN"     "DEC"     "DG1"     "DG2"     "DG3"     "DLFC"    "ENT"     "F"      
#[23] "GP"      "HY"      "I"       "II"      "III"     "IV"      "IX"      "L1"      "L2"      "L3"      "L4"     
#[34] "L5"      "L6"      "L7"      "L8"      "LSNSH"   "M"       "M1"      "M2"      "M3"      "M4"      "M5"     
#[45] "M6"      "M7"      "M8"      "MSC"     "MTG"     "NTS"     "OLF"     "P"       "PF"      "PRM"     "S1"     
#[56] "SIM"     "SIMANDB" "SNrSNc"  "T1"      "T2"      "T3"      "T4"      "T5"      "T6"      "T7"      "T8"     
#[67] "TRSSF"   "V"       "VI"      "VII"     "VIII"    "VIS"     "VISP"    "X"       ""       

#Filtering for sub_region 'CP' achieves this in this data set
STRNOLSX_metadata_Libs <- metadata %>% filter(sub_region %in% "CP")

#Filter cells for metadata provided library
STRNOLSX_adfiledopa_Table <- adfiledopa_Table %>% filter(rn %in% STRNOLSX_metadata_Libs$library)

#Filter cells to remove non-microglia
# >0 can be used in place of != 0. Both work as there should be no cells with a negative number of genes
STRNOLSX_adfiledopa_Table <- STRNOLSX_adfiledopa_Table %>% filter(ENSMUSG00000024397 != 0)

#View and ensure you know which col is D1 and D2, then change col names
colSums(STRNOLSX_adfiledopa_Table != 0)
dopagenes

#Rename cols
colnames(STRNOLSX_adfiledopa_Table) <- c("STR_Cells","Cell_code","Drd1","Drd2","IBA1")
STRNOLSX_Counts <- colSums(STRNOLSX_adfiledopa_Table != 0)
STRNOLSX_Counts

#Send to excel
write.xlsx(STRNOLSX_Counts, file="Macosko-Output-Raw.xlsx", sheetName="Striatum_No_LSX", append=TRUE)

#-------------------------------------------------------------------------------
#FRONTAL CORTEX ANALYSIS
#-------------------------------------------------------------------------------

#Remove all cells not in the cortex

unique(metadata$brain_struct)
#Possible brain structs
#"Isocortex" "CTXsp"     "CB"        "HPF"       "STR"       "PAL"       "HY"        "P"         "MB"       
#"MY"        "OLF"       "TH"    

unique(metadata$region)
#Possible brain regions
#[1] "MOp"  "RSP"  "ACA"  "AMY"  "CB"   "CTX"  "AUD"  "HPF"  "STRd" "DCN"  "ENT"  "PALd" "HY"   "BS"   "MB"   "LSX"  "SUB" 
#[18] "PALm" "NTS"  "OLF"  "S1"   "PALv" "TH"   "VIS"  "VISP" "BNST"

#Use MOp

unique(metadata$sub_region)
#Possible brain sub regions
#[1] "A"       "ACA"     "AMY"     "AN1"     "AN2"     "ANG"     "AUD"     "CA1"     "CA2"     "CA2CA3"  "CA3"    
#[12] "COP"     "CP"      "CUL"     "DCN"     "DEC"     "DG1"     "DG2"     "DG3"     "DLFC"    "ENT"     "F"      
#[23] "GP"      "HY"      "I"       "II"      "III"     "IV"      "IX"      "L1"      "L2"      "L3"      "L4"     
#[34] "L5"      "L6"      "L7"      "L8"      "LSNSH"   "M"       "M1"      "M2"      "M3"      "M4"      "M5"     
#[45] "M6"      "M7"      "M8"      "MSC"     "MTG"     "NTS"     "OLF"     "P"       "PF"      "PRM"     "S1"     
#[56] "SIM"     "SIMANDB" "SNrSNc"  "T1"      "T2"      "T3"      "T4"      "T5"      "T6"      "T7"      "T8"     
#[67] "TRSSF"   "V"       "VI"      "VII"     "VIII"    "VIS"     "VISP"    "X"       ""       

#Further subdivide for MOp sub regions (A,M), excludes parts of prefrontal cortex

#Filter for Isocortex
ISO_metadata_Libs <- metadata %>% filter(brain_struct %in% "Isocortex")

#Filter for MOp
ISO_metadata_Libs <- ISO_metadata_Libs %>% filter(region %in% "MOp")

#Filter for sub regions 'A' and 'M'
ISO_metadata_Libs <- ISO_metadata_Libs %>% filter(sub_region %in% c("A","M"))

#Filter cells for metadata provided library
ISO_adfiledopa_Table <- adfiledopa_Table %>% filter(rn %in% ISO_metadata_Libs$library)

#Filter cells to remove non-microglia
# >0 can be used in place of != 0. Both work as there should be no cells with a negative number of genes
ISO_adfiledopa_Table <- ISO_adfiledopa_Table %>% filter(ENSMUSG00000024397 != 0)

#View and ensure you know which col is D1 and D2, then change col names
colSums(ISO_adfiledopa_Table != 0)
dopagenes

#Rename cols
colnames(ISO_adfiledopa_Table) <- c("FrtCTX_Cells","Cell_code","Drd1","Drd2","IBA1")
FrtCTX_Counts <- colSums(ISO_adfiledopa_Table != 0)
FrtCTX_Counts

#send to excel
write.xlsx(FrtCTX_Counts, file="Macosko-Output-Raw.xlsx", sheetName="FrtCTX", append=TRUE)
