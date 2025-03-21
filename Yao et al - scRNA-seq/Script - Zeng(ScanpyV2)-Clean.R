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


adfile <- read_h5ad("WMB-10Xv3-STR-log2.h5ad", backed = "r")

#Associate ensemble gene names with genes of interest
dopagenes <- data.frame(genenum = c("ENSMUSG00000021478", "ENSMUSG00000032259","ENSMUSG00000024397"),
                        genesym = c("DRD1","DRD2","IBA1"))

#Isolate for dopamine receptor genes and IBA1
adfiledopa <- sc$get$obs_df(adfile,
                            keys = dopagenes$genenum)

#Convert to data.table
adfiledopa_Table <- data.table(adfiledopa,keep.rownames=TRUE)

#Get cell metadata, coordinates, and parcellation index of post QC cells
metadata <- fread("cell_metadata_with_cluster_annotation.csv")
regions <- fread("region_of_interest_metadata.csv")

#Remove unnecessary cols
colnames(metadata)
metadata <- metadata[,-c(2:8,10:18,20:28)]

colnames(regions)
regions <- regions[,-c(1,4:5)]

#Merge meta data
metadata <- merge(metadata, regions, by.x = "region_of_interest_acronym", by.y = "acronym", all.x = TRUE)

#-------------------------------------------------------------------------------
#STRIATUM (DORSAL) ANALYSIS
#-------------------------------------------------------------------------------

#Remove all cells not in the dorsal striatum

unique(metadata$region_of_interest_acronym)
#Possible brain structs
#[1] "RHP"              "RSP"              "ACA"              "PL-ILA-ORB"       "AUD-TEa-PERI-ECT" "SS-GU-VISC"      
#[7] "MO-FRP"           "PAL"              "sAMY"             "CTXsp"            "HY"               "STRv"            
#[13] "OLF"              "LSX"              "AI"               "STRd"             "VIS-PTLp"         "VIS"             
#[19] "TH"               "MOp"              "ENT"              "HIP"              "P"                "MB"              
#[25] "MY"               "CB"               "AUD"              "SSp"              "TEa-PERI-ECT"       

unique(metadata$name)
#Possible brain regions
# [1] "Anterior cingulate area"                       "Agranular insular area"                       
#[3] "Auditory areas"                                "Auditory/temporal/perirhinal/ectorhinal areas"
#[5] "Cerebellum"                                    "Cortical subplate"                            
#[7] "Entorhinal area"                               "Hippocampal region"                           
#[9] "Hypothalamus"                                  "Lateral septal complex"                       
#[11] "Midbrain"                                      "Somatomotor - Frontal pole"                   
#[13] "Primary motor area"                            "Medulla"                                      
#[15] "Olfactory areas"                               "Pons"                                         
#[17] "Pallidum"                                      "Prelimbic/infralimbic/orbital areas"          
#[19] "Retrohippocampal region"                       "Retrosplenial area"                           
#[21] "Somatosensory/gustatory/visceral areas"        "Primary somatosensory area"                   
#[23] "Striatum dorsal region"                        "Striatum ventral region"                      
#[25] "Temporal/perirhinal/ectorhinal areas"          "Thalamus"                                     
#[27] "Visual areas"                                  "Visual/posterior parietal  areas"             
#[29] "Striatum-like amygdalar nuclei"

#Use "STRd"
STRDOR_metadata_Libs <- metadata %>% filter(region_of_interest_acronym %in% "STRd")

#Filter cells for metadata provided striatal cells
STRDOR_adfiledopa_Table <- adfiledopa_Table %>% filter(rn %in% STRDOR_metadata_Libs$cell_label)

#Filter cells to remove non-microglia
# >0 can be used in place of != 0. Both work as there should be no cells with a negative number of genes
STRDOR_adfiledopa_Table <- STRDOR_adfiledopa_Table %>% filter(ENSMUSG00000024397 != 0)

#View and ensure you know which col is D1 and D2, then change col names
colSums(STRDOR_adfiledopa_Table != 0)
dopagenes

#Rename cols
colnames(STRDOR_adfiledopa_Table) <- c("STR_Cells","Drd1","Drd2","IBA1")
StrDOR_Counts <- colSums(STRDOR_adfiledopa_Table != 0)
StrDOR_Counts

#Send to excel
write.xlsx(StrDOR_Counts, file="Zeng-Output-Raw.xlsx", sheetName="StriatumDOR", append=TRUE)

#-------------------------------------------------------------------------------
#STRIATUM (VENTRAL) ANALYSIS
#-------------------------------------------------------------------------------

#Remove all cells not in the ventral striatum

unique(metadata$region_of_interest_acronym)
#Possible brain structs
#[1] "RHP"              "RSP"              "ACA"              "PL-ILA-ORB"       "AUD-TEa-PERI-ECT" "SS-GU-VISC"      
#[7] "MO-FRP"           "PAL"              "sAMY"             "CTXsp"            "HY"               "STRv"            
#[13] "OLF"              "LSX"              "AI"               "STRd"             "VIS-PTLp"         "VIS"             
#[19] "TH"               "MOp"              "ENT"              "HIP"              "P"                "MB"              
#[25] "MY"               "CB"               "AUD"              "SSp"              "TEa-PERI-ECT"       

unique(metadata$name)
#Possible brain regions
# [1] "Anterior cingulate area"                       "Agranular insular area"                       
#[3] "Auditory areas"                                "Auditory/temporal/perirhinal/ectorhinal areas"
#[5] "Cerebellum"                                    "Cortical subplate"                            
#[7] "Entorhinal area"                               "Hippocampal region"                           
#[9] "Hypothalamus"                                  "Lateral septal complex"                       
#[11] "Midbrain"                                      "Somatomotor - Frontal pole"                   
#[13] "Primary motor area"                            "Medulla"                                      
#[15] "Olfactory areas"                               "Pons"                                         
#[17] "Pallidum"                                      "Prelimbic/infralimbic/orbital areas"          
#[19] "Retrohippocampal region"                       "Retrosplenial area"                           
#[21] "Somatosensory/gustatory/visceral areas"        "Primary somatosensory area"                   
#[23] "Striatum dorsal region"                        "Striatum ventral region"                      
#[25] "Temporal/perirhinal/ectorhinal areas"          "Thalamus"                                     
#[27] "Visual areas"                                  "Visual/posterior parietal  areas"             
#[29] "Striatum-like amygdalar nuclei"

#Use "STRv"
STRVEN_metadata_Libs <- metadata %>% filter(region_of_interest_acronym %in% "STRv")

#Filter cells for metadata provided striatal cells
STRVEN_adfiledopa_Table <- adfiledopa_Table %>% filter(rn %in% STRVEN_metadata_Libs$cell_label)

#Filter cells to remove non-microglia
# >0 can be used in place of != 0. Both work as there should be no cells with a negative number of genes
STRVEN_adfiledopa_Table <- STRVEN_adfiledopa_Table %>% filter(ENSMUSG00000024397 != 0)

#View and ensure you know which col is D1 and D2, then change col names
colSums(STRVEN_adfiledopa_Table != 0)
dopagenes

#Rename cols
colnames(STRVEN_adfiledopa_Table) <- c("STR_Cells","Drd1","Drd2","IBA1")
StrVEN_Counts <- colSums(STRVEN_adfiledopa_Table != 0)
StrVEN_Counts

#Send to excel
write.xlsx(StrVEN_Counts, file="Zeng-Output-Raw.xlsx", sheetName="StriatumVEN", append=TRUE)
