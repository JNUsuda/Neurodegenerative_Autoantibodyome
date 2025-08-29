# 1. Data processing -----------------------------------------------------------

## 1.1 Load packages and functions ---------------------------------------------

source("./codes/my_packages.R")
source("./codes/my_functions.R")

## 1.2 sample selection and metadata  ------------------------------------------
# Nagele et al. datasets 
# Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0


# list of datasets by Nagele et al.
# platform 5.0
Nagele_GSE_list50 = c("GSE62283", "GSE39087", "GSE74763")
# platform 5.1
Nagele_GSE_list51 = c("GSE95718", "GSE137422")

# use my_metadata_function() to retrieve metadata from GEO
sample_list50 = lapply(Nagele_GSE_list50, my_metadata_function)
sample_list51 = lapply(Nagele_GSE_list51, my_metadata_function)

names(sample_list50) = Nagele_GSE_list50
names(sample_list51) = Nagele_GSE_list51
df <- do.call("rbind", c(sample_list50, sample_list51))

# Merged metadata was filtered according to:
# "Duplicates, non-healthy controls, samples with unknown age or sex, and 
# samples with conflicting information between datasets".
df2 = df %>% 
  dplyr::select(-filename) %>% 
  # remove samples without age and sex information
  dplyr::filter(Age != "unk" & Age != "unknown" & Sex != "unk" &
  # remove non healthy control and non neurodegenerative samples
                Group != "cdr=0 plasma" & Group != "breast cancer") %>% 
  # remove duplicates
  dplyr::distinct(PA_ID, Age, Sex, Group, .keep_all = T) %>% 
  # remove samples with conflicting information between datasets
  dplyr::filter(PA_ID != "PA25475") %>% 
  # add major group column
  dplyr::mutate(grupao = 
                if_else(Group %in% c("alzheimer's disease",
                                     "cdr=0.5 plasma"), true = "AD",
                false =if_else(Group == "parkinson's disease", true = "PD",
                false =if_else(Group == "multiple sclerosis", true = "MS",
                                          false = "control"))))
# rename groups in Group column
df2$Group[grep(df2$Group_original, pattern = "alzheimer's disease")] = "AD"
df2$Group[grep(df2$Group_original, pattern = "^parkinson's disease$")] = "PD"
df2$Group[grep(df2$Group_original, pattern = "^early-stage parkinson's disease$")] = "earlyPD"
df2$Group[grep(df2$Group_original, pattern = "^multiple sclerosis$")] = "MS"
df2$Group[grep(df2$Group_original, pattern = "^secondary progressive multiple sclerosis$")] = "spMS"
df2$Group[grep(df2$Group_original, pattern = "^relapsing-remitting multiple sclerosis$")] = "rrMS"
df2$Group[grep(df2$Group_original, pattern = "cdr=0.5 plasma")] = "cdr05"
df2$Group[grep(df2$Group_original, pattern = "control")] = "control"
df2$Group[grep(df2$title, pattern = "LMCI")] = "LMCI"
df2$Group[grep(df2$title, pattern = "EMCI")] = "EMCI"

# add Group2 column
# "EMCI", "LMCI" and "cdr05" were joined as the "MCI" category.
df2 = df2 %>% 
  mutate(Group2 = gsub(Group, pattern = "cdr05|EMCI|LMCI", replacement = "MCI")) %>% 
# add dataset column
  mutate(dataset = if_else(GSE %in% Nagele_GSE_list50, true = "5.0", false = "5.1")) %>% 
# rename and reorder columns
  dplyr::mutate(FileName = gsub(supplementary_file, 
                                pattern = ".*/|\\.gz$", replacement = "")) %>% 
  dplyr::relocate(FileName, .before = title) %>% 
  dplyr::rename(ArrayID = title) %>% 
  dplyr::relocate(Group, Group_original, .after = ArrayID) %>% 
  dplyr::relocate(PA_ID, dataset, Group2, grupao, .after = GSM) %>% 
  dplyr::select(-supplementary_file) %>% 
  dplyr::arrange(Group, Sex, Age, PA_ID)

write_excel_csv(x = df2, file = "./data/merged_filtered_metadata.csv")

# Resulting table of selected samples are available in the Supplementary tables.

## 1.3 Prepare for processing --------------------------------------------------
# prepare folders with raw data and targets file (for PAA)

# all .gpr files from Nagele datasets are in:"./data/Nagele_datasets/merged_raw"
# gz_files = list.files(path = "./data/Nagele_datasets/merged_raw", full.names = T, recursive = T)

# list of datasets by Nagele et al.
Nagele_GSE_list = c("GSE62283", "GSE39087", "GSE74763", "GSE95718", "GSE137422")

for (item in Nagele_GSE_list) {
  my_GSE = item
  getwd()
  
  # set dataset folder as working directory
  my_folder = list.files(pattern = my_GSE, path = "./data", 
                         recursive = F, full.names = T)
  # if folder doesn't exist, create folder:
  if (is_empty(my_folder)) {
    my_path = paste0("./data/", my_GSE)
    dir.create(path = my_path)
    my_folder = my_path
  }
  
  setwd(my_folder)
  getwd() # check if correctly set
  
  ## get folder names
  raw_folder = paste0("./", my_GSE, "_RAW") # ex. "./GSE50866_RAW"
  raw_name = paste0(my_GSE, "_RAW")  # ex. "GSE50866_RAW"
  raw_tar_name = paste0(my_GSE, "_RAW.tar")  # ex. "GSE50866_RAW.tar"
  targets_file = paste0("./", my_GSE, "_targets.txt") # ex. "./GSE50866_targets.txt"
  
  # download data and make targets file for PAA
  
  # if "GSE137422_RAW.tar" is NOT in the files in the folder
  # download supplementary data from GEO, including raw data, and decompress
  if (!(raw_tar_name %in% list.files(path = "./"))) {
    # load data from GEO
    gset <- getGEO(my_GSE, GSEMatrix =TRUE, AnnotGPL=FALSE)
    gset <- gset[[1]]
    # load supplementary data from GEO
    suppl_files = getGEOSuppFiles(my_GSE, makeDirectory = FALSE, baseDir = getwd(),
                                  fetch_files = T, filter_regex = "RAW")
    
    # decompress raw data from .tar and .gz
    untar(raw_tar_name, exdir = raw_name)
    gz_files = list.files(path = raw_folder, full.names = T)
    for (x in gz_files) {gunzip(filename = x)}
  }
  
  # get metadata with my_metadata_function
  meta = read.csv(file = "../../data/merged_filtered_metadata.csv")
  
  ## targets file for PAA
  # get the .gpr files and their GSM
  FileName = list.files(path = raw_folder, pattern = "[.]gpr$", full.names = F)
  targets2 = tibble(FileName)
  targets2 = targets2 %>% dplyr::mutate(
    GSM = gsub(x = targets2$FileName, pattern = "_.*", replacement = ""))
  # join the tables by GSM
  #targets3 = full_join(targets2, meta, by = "GSM")
  targets3 = inner_join(meta, targets2, by = "GSM")
  # subset the columns for the filename, array ID and group, rename columns accordingly
  targets4 = targets3[,c("FileName.x", "ArrayID", "Group","Group_original","Age", 
                         "Sex", "GSE", "GSM", "PA_ID", "dataset", "Group2", "grupao")]
  colnames(targets4) = c("FileName", "ArrayID", "Group", "Group_original", "Age", 
                         "Sex", "GSE", "GSM", "PA_ID","dataset", "Group2", "grupao")
  # export targets .txt for PAA
  write_tsv(x = targets4, file = targets_file)
  
  setwd("../../")
}


# merge datasets of the same platform:

# list of datasets by Nagele et al.
# platform 5.0
Nagele_GSE_list50 = c("GSE62283", "GSE39087", "GSE74763")
# platform 5.1
Nagele_GSE_list51 = c("GSE95718", "GSE137422")

for (my_list in list(Nagele_GSE_list50, Nagele_GSE_list51)) {
  # select the platform
  if (identical(my_list, Nagele_GSE_list50)) {
    platform = "5.0"} else {platform = "5.1"}
  # select folder to save
  my_folder = list.dirs(path = "./data", recursive = F, full.names = T)
  # if folder doesn't exist, create folder:
  if (any(grepl(my_folder, pattern = "5.0", fixed = T)) == F) {
    my_path = paste0("./data/", "5.0")
    dir.create(path = my_path)}
  if (any(grepl(my_folder, pattern = "5.1", fixed = T)) == F) {
    my_path = paste0("./data/", "5.1")
    dir.create(path = my_path)}
  
  # get targets tables  
  targets_files = list()
  for (item in my_list) {
    targets_file = paste0("./data/", item, "/", item, "_targets.txt") # ex. "./GSE50866/GSE50866_targets.txt"
    targets = read.delim(targets_file)
    targets_files[[item]] = targets
  }
  
  # join targets tables
  merged_targets = do.call("rbind", targets_files)
  merged_targets_file = paste0("./data/", platform, "/merged_", platform, 
                               "_filtered_targets.txt")
  write_delim(file = merged_targets_file, x = merged_targets, delim = "\t", col_names = T)
  
  # move raw files to platform folder:
  merged_targets = merged_targets %>% 
    dplyr::mutate(paths = paste0("./data/", GSE, "/", GSE, "_RAW/", FileName))
  # select folder to save
  my_folder = list.dirs(path = paste0("./data/", platform), recursive = F, full.names = T)
  # if folder doesn't exist, create folder:
  if (any(grepl(my_folder, pattern = "merged_raw")) == F) {
    my_path = paste0("./data/", platform, "/merged_raw_", platform, "_filtered")
    dir.create(path = my_path)}
  file.copy(from = merged_targets$paths, 
            to = paste0("./data/", platform, "/merged_raw_", platform, "_filtered/"))
}

## 1.4 Data processing ----------------------------------------------------------

# Processing of raw (.gpr, fluorescence) data of protein arrays 
# (ProtoArray, HuProt) to 'counts' table.
# Protein array R packages: PAA and/or paweR.

# Use the code under "PAA" section for cyclic loess normalization, (recommended)
# or use the code under "paweR" section for rlm normalization.

# steps:
# 1. .gpr files loading
# 2. background correction
# 3. normalization and aggregation and gene symbol conversion 

# install packages:
## paweR:
#library(remotes)
#install_gitlab(host = "https://gl.cs.ut.ee/", repo = "biit/paweR")

## PAA:
#library("BiocManager")
#BiocManager::install("PAA", dependencies=TRUE)


# load libraries
library(PAA)
library(limma)
library(paweR)

# set palette
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
          "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))


# set directories - REPEAT FOR 5.1
setwd("./data/5.0")
getwd() # check if correctly set

raw_folder = "./merged_raw_5.0_filtered"
targets_file = "./merged_5.0_filtered_targets.txt"
metadata = read.delim(file = targets_file, header = T)

### PAA --------------

# Step 1. load .gpr files with PAA:

raw_folder # check folder
list.files(path = raw_folder) # check contents (.gpr files)
length(list.files(path = raw_folder)) # check number of files

# check .gpr files for "description" column 
## (information which spot is a feature, control or to be discarded)

# no aggregation performed in this step (only after normalization)

# Step 1. Load files
# load .gpr files (PAA):
my_elist = PAA::loadGPR(gpr.path = raw_folder,
                        array.type = "ProtoArray",
                        aggregation = "none",
                        targets.path = targets_file, 
                        array.columns = list(E = "F635 Median", Eb = "B635 Median"),
                        array.annotation = c("Block", "Column", "Row", 
                                             "Description", "Name", "ID"))

# visual inspection of arrays for artifacts (PAA):
# !!! saves images to working directory 
# background
plotArray(elist=my_elist, idx="all", data.type="bg", log=FALSE, normalized=FALSE,
          aggregation="none", colpal="topo.colors",
          graphics.device = "png", output.path = getwd())
# foreground
plotArray(elist=my_elist, idx="all", data.type="fg", log=FALSE, normalized=FALSE,
          aggregation="none", colpal="topo.colors",
          graphics.device = "png", output.path = getwd())
# rename folder
file.rename(from = "./array_plots" , to = "./array_plots_raw")


# Step 2. Background correction (limma)
# set desired offset
my_elist2 <- backgroundCorrect(my_elist, method="normexp",
                               normexp.method="saddle", offset = 50)


# check arrays after bg correction (PAA)
# !!! saves images to working directory 
plotArray(elist=my_elist2, idx="all", data.type="fg", log=FALSE, normalized=FALSE,
          aggregation="none", colpal="topo.colors",
          graphics.device = "png", output.path = getwd())
# rename folder to bg-corrected
file.rename(from = "./array_plots" , to = "./array_plots_bg-corrected")


# Plot boxplots with different normalization methods, to choose the method
# !!! saves images to working directory 
PAA::plotNormMethods(elist=my_elist2, include.rlm = T, controls = "internal",
                     output.path = getwd())
#plotMAPlots(elist=my_elist2, idx=1, include.rlm = T, controls = "internal")


# Step 3. Normalization (PAA)
# for rlm (robust linear modelling) normalization, 
# internal controls (HumanIgG and Anti-HumanIgG), as default in PAA.
# controls are removed
# no aggregation perfomed
my_elist3 = PAA::normalizeArrays(elist = my_elist2, 
                                 method = "cyclicloess"
                                 #, controls = "internal"
                                 #,output.path = getwd()
)

# check arrays after normalization (PAA)
# !!! saves images to working directory 
plotArray(elist=my_elist3, idx="all", data.type="fg", log=TRUE, normalized=TRUE,
          aggregation="none", colpal="topo.colors",
          graphics.device = "png", output.path = getwd())
# rename folder to add normalization
file.rename(from = "./array_plots" , to = "./array_plots_bg-corrected_normalized")


# Step 3. Convert IDs and aggregate duplicates (paweR)

# convert RefSeq IDs into ENSG IDs and Gene Names
my_elist4 = paweR::cleanNames(my_elist3)

# aggregate duplicates by median, based on ENSG
my_elist5 = aggregateElist(my_elist4, method = "median", by = "ENSG")

# skip to plots section

### PAWER --------------
# Step 1. Load files
# load .gpr files (paweR):
my_elist_paweR = paweR::readGPR(path = raw_folder,
                                backMethod = NULL,
                                source = "genepix.median",
                                columns = list(E = "F635 Median", Eb = "B635 Median"),
                                annotations = c("Block", "Column", "Row", 
                                                "Description", "Name", "ID"))
#my_elist_paweR$printer = my_elist$printer

# Step 2. Background correction (limma)
# set desired offset.
my_elist_paweR2 <- backgroundCorrect(my_elist_paweR, method="normexp",
                                     normexp.method="saddle", offset = 50)

# Step 3. Normalization 
# rlm (robust linear modelling) normalization
# internal controls (HumanIgG and Anti-HumanIgG), as default in PAA
# controls are removed
# aggregation by mean 

# control list:
control_list = c("HumanIgG1~N/A", 
                 "HumanIgG2~N/A", 
                 "HumanIgG3~N/A", 
                 "HumanIgG4~N/A",
                 "Anti-HumanIgG1~N/A",
                 "Anti-HumanIgG2~N/A",
                 "Anti-HumanIgG3~N/A",
                 "Anti-HumanIgG4~N/A")
control_list2 = as.list(control_list)
names(control_list2) = rep("Name", length(control_list2))

# Convert IDs and aggregate duplicates (paweR)

# convert RefSeq IDs into ENSG IDs and Gene Names
my_elist_paweR3 = paweR::cleanNames(my_elist_paweR2)

# normalize with rlm and aggregate duplicates by mean, based on ENSG
my_elist_paweR4 = paweR::normalise(elist = my_elist_paweR3, verbose = TRUE,
                                   aggMethod = "mean", aggBy = "ENSG", maxit = 10,
                                   allControls = F, controls = control_list2)

# !!! be careful if also running the section under 'PAA', will overwrite results:
my_elist5 = my_elist_paweR4 # for plots

### plots -----

# set directory to download results
newfolder = "./cyclicloess-normalized_5.0"
dir.create(path = newfolder)
setwd(newfolder)
getwd()

## boxplot - check distribution
# !!! saves images to working directory - set desired name
groups = c("Group", "GSE")
for (item in groups) {
  nome = paste0("boxplot_processed_", item, ".png")
  png(nome, width = 10000, height = 2000, pointsize = 16, res = 300)
  par(mar=c(1,3,1,1))
  boxplot(my_elist5$E, col = as.factor(metadata[,item]),
          boxwex=0.6, notch=T, outline=FALSE, las=2, xaxt='n', ann=FALSE)
  legend("topleft", levels(as.factor(metadata[,item])), fill=palette(), bty="n")
  dev.off()
} 

## PCAs
pca = prcomp(t(my_elist5$E), scale. = T)
# set desired variables from metadata
groups = c("Group", "Group_original", "Age", "Sex", "GSE")
tabela = cbind(metadata[,groups], pca$x)

# age (continuous variable) effect
tabela %>% 
  ggplot(aes(x = PC1, y=PC2, col= Age)) + 
  geom_point() + scale_colour_gradientn(colours=rainbow(5))

# !!! saves images to working directory - set desired name
ggsave(filename = "PCA_Age.png",  width = 4000, height = 3000,
       device = "png", units = "px", dpi = 500)

# categorical variables:
# set desired variables from metadata
groups = c("Group", "Group_original", "Sex", "GSE")
for (item in groups) {
  group <- factor(tabela[,item])
  fviz_pca_ind(pca,
               col.ind = group, # color by groups
               geom = c("point"),
               #palette = palette(),
               addEllipses = TRUE, # Concentration ellipses
               ellipse.type = "confidence",
               legend.title = "Groups",
               repel = TRUE
  ) + theme(legend.position = "bottom") + ggtitle("PCA") 
  
  nome = paste0("PCA_", item, ".png")
  ggsave(filename = nome,  width = 4000, height = 3000,
         device = "png", units = "px", dpi = 500)
}

## download tables

# download EList
tab_E = rownames_to_column(my_elist5$E, "rowname")
tab_targets = rownames_to_column(my_elist5$targets, "rowname")
tab_genes = rownames_to_column(my_elist5$genes, "rowname")

# check if columns (samples) are in the same order in the table and metadata
if (all(tab_targets$ArrayID == colnames(tab_E[-1]))) {
  print("columns match - safe to proceed")
} else {
  print("ATTENTION - SAMPLE ORDER DOESN'T MATCH")}  


write.table(x = tab_E, file = "./merged_normalized_E.tsv", 
            col.names = T, row.names = F, sep = "\t", quote = T)
write.table(x = tab_targets, file = "./merged_normalized_targets.tsv", 
            col.names = T, row.names = F, sep = "\t", quote = T)
write.table(x = my_elist5$genes, file = "./merged_normalized_genes.tsv", 
            col.names = T, row.names = F, sep = "\t", quote = T)


## arrange tables (for DEA)
# normalized values by patient 
norm_table = merge(x = my_elist5$genes[,c("ENSG", "GeneName")], 
                   y = my_elist5$E, 
                   by.x = "ENSG", by.y = "row.names") # by row names

any(duplicated(norm_table$ENSG)) # check for duplicates
any(duplicated(norm_table$GeneName)) # check for duplicates
anyNA(norm_table) # check for NA
#noNA_norm_table = na.omit(norm_table) # remove lines with no gene symbol/ensembl gene

if (all(tab_targets$ArrayID == colnames(norm_table)[-(1:2)])) {
  print("columns match - safe to proceed")
  colnames(norm_table) = c("ENSG", "GeneName", tab_targets$GSM)
} else {
  print("ATTENTION - SAMPLE ORDER DOESN'T MATCH")}  

write.table(x = norm_table, file = "./merged_normalized_table.tsv", 
            col.names = T, row.names = F, sep = "\t", quote = T)

setwd("../../../")
# (repeat for 5.1)

## 1.5 Batch effect correction -------------------------------------------------

# create directory to save files
dir.create("./data/GSEbatch")

# remove batch effect from normalized tables

# 'expression' tables
tabela1 <- read.table(file = "./data/5.0/cyclicloess-normalized_5.0/merged_normalized_E.tsv", header = T , 
                      sep = "\t", dec = ".", row.names = 1, 
                      comment.char = "", check.names = F)
tabela2 <- read.table(file = "./data/5.1/cyclicloess-normalized_5.1/merged_normalized_E.tsv", header = T , 
                      sep = "\t", dec = ".", row.names = 1, 
                      comment.char = "", check.names = F)

# merge tables
tabela3 = merge(tabela1, tabela2, by = "row.names")
tabela3 = column_to_rownames(tabela3, "Row.names")


# metadata tables
meta1 <- read.table(file = "./data/5.0/cyclicloess-normalized_5.0/merged_normalized_targets.tsv", header = T , 
                    sep = "\t", dec = ".", row.names = 1, comment.char = "")
meta2 <- read.table(file = "./data/5.1/cyclicloess-normalized_5.1/merged_normalized_targets.tsv", header = T , 
                    sep = "\t", dec = ".", row.names = 1, comment.char = "")

# merge metadata tables
meta1$dataset = "5.0"
meta2$dataset = "5.1"
meta3 = rbind(meta1, meta2)


if (all(meta3$ArrayID == colnames(tabela3))) {
  print("columns match - safe to proceed")
} else {
  print("ATTENTION - COLUMNS DON'T MATCH")} 
# update colnames to GSM
colnames(tabela3) = meta3$GSM

matriz <- as.matrix(tabela3)
## set batch (in this case, to GSE column)
lote = as.numeric(as.factor(meta3$GSE))
# batch effect correction with parametric adjustment
combat1 = ComBat(dat=matriz, batch=lote, mod=NULL, par.prior=TRUE, prior.plots=FALSE)

# PCA
# set groups for PCA
group = as.factor(meta3$Group2)

for (item in list(matriz, combat1)) {
  pca = prcomp(t(item), scale. = T)
  fviz_pca_ind(pca,
               col.ind = group, # color by groups
               geom = c("point"),
               #palette = palette(),
               addEllipses = TRUE, # Concentration ellipses
               ellipse.type = "confidence",
               legend.title = "Groups",
               repel = TRUE
  ) + theme(legend.position = "bottom") + ggtitle("PCA") 
  
  if (identical(item, matriz)) {
    nome = paste0("./data/GSEbatch/PCA_beforebatch.png")} else {
      nome = paste0("./data/GSEbatch/PCA_afterbatch.png")}
    ggsave(filename = nome,  width = 4000, height = 3000,
           device = "png", units = "px", dpi = 500)
  ## adjust names to add the grouping variable
}

# get gene tables (for conversion)
tabelagenes1 <- read.table(file = "./data/5.0/cyclicloess-normalized_5.0/merged_normalized_genes.tsv", header = T , 
                           sep = "\t", dec = ".", comment.char = "")
tabelagenes2 <- read.table(file = "./data/5.1/cyclicloess-normalized_5.1/merged_normalized_genes.tsv", header = T , 
                           sep = "\t", dec = ".", comment.char = "")

# merge tables
tabelagenes3 = rbind(tabelagenes1[,c("ENSG", "GeneName")],tabelagenes2[,c("ENSG", "GeneName")])
tabelagenes3 = unique(tabelagenes3)


# join with combat table to get genenames
combat2 = as.data.frame(combat1)
combat2 = rownames_to_column(combat2, var = "ENSG")
combat3 = right_join(tabelagenes3, combat2, by = join_by(ENSG))


# save tables
setwd("./data/GSEbatch/")
getwd()
write.table(x = combat3, file = "./merged_normalized_combat_table.tsv", sep = "\t", row.names = F, col.names = T)
write.table(x = meta3, file = "./merged_normalized_targets.tsv", sep = "\t", row.names = F, col.names = T)

setwd("../../")
