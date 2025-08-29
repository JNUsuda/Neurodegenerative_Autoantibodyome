# 2. Differential Reactivity Analysis ------------------------------------------

## 2.1 Load packages and functions ---------------------------------------------

source("./codes/my_packages.R")
source("./codes/my_functions.R")

## 2.2 Differential reactivity analysis ----------------------------------------

# with limma package
# compare levels of each subgroup with control, filter by adj. p-value and FC
# return complete tables and summarized lists of the differentially reactive 
# autoantibodies (DRAs)

setwd("./data/GSEbatch/")
getwd()

tabela0 = read.table(file = "./merged_normalized_combat_table.tsv", 
                     header = T, sep = "\t", check.names = F)
metadata = read.table(file = "./merged_normalized_targets.tsv", 
                      header = T, sep = "\t", check.names = F)

# set rownames and remove gene columns
tabela1 = column_to_rownames(tabela0, "ENSG")
tabela = tabela1[,-c(1)]

# check if columns (samples) are in the same order in the table and metadata
if (all(metadata$GSM == colnames(tabela))) {
  print("columns match - safe to proceed")
} else {
  print("ATTENTION - COLUMNS DON'T MATCH")}  

# set group factor
ordered_groups = levels(factor(metadata$Group2))
ordered_groups
# set name of the reference group
nome_control = "control"
groups = ordered_groups[! ordered_groups %in% c(nome_control)] 
Group = as.factor(metadata$Group2)
# set control group as reference
Group <- relevel(Group, ref=nome_control)
levels(Group)

design <- stats::model.matrix(~ Group + metadata$Sex + metadata$Age)
colnames(design)
colnames(design) = gsub(colnames(design), pattern = "Group|metadata\\$", replacement = "")
colnames(design) = gsub(colnames(design), pattern = ":", replacement = "_")
colnames(design)


fit <- limma::lmFit(tabela, design)
fit2 <- limma::eBayes(fit)
upper_res_list = list()
meta_list = list()
title = "GSEbatch"
salvar = T


i=2 # skip 1, the (intercept) coefficient
while (i <= length(colnames(design))) {
  
  full_results = topTable(fit2, number=Inf, coef = i, adjust.method = "BH", confint = T)
  full_results = tibble::rownames_to_column(full_results, "ENSG")
  # get the gene names - merge with the gene table
  full_results2 = merge(x = tabela0[,c(1:2)], y = full_results, 
                        by = "ENSG", all.y = T)
  meta_list[[(colnames(design)[i])]] = full_results2
  
  # filter results by adjusted p-value and FC
  filtered_results = full_results2 %>% dplyr::filter(adj.P.Val < 0.05)
  filtered_results_UP = filtered_results %>% dplyr::filter(logFC > 0)
  filtered_results_DOWN = filtered_results %>% dplyr::filter(logFC < -0)
  
  res_list = list()
  res_list[[paste0(colnames(design)[i],"_TOTAL")]] = c(filtered_results_UP$GeneName, 
                                                       filtered_results_DOWN$GeneName)
  res_list[[paste0(colnames(design)[i],"_UP")]] = filtered_results_UP$GeneName
  res_list[[paste0(colnames(design)[i],"_DOWN")]] = filtered_results_DOWN$GeneName
  # convert list to df
  df = l2df(res_list, byrow = F)
  
  # add results from the condition to the list with all conditions
  upper_res_list[[paste0(colnames(design)[i],"_TOTAL")]] = res_list[[paste0(colnames(design)[i],"_TOTAL")]]
  upper_res_list[[paste0(colnames(design)[i],"_UP")]] = res_list[[paste0(colnames(design)[i],"_UP")]]
  upper_res_list[[paste0(colnames(design)[i],"_DOWN")]] = res_list[[paste0(colnames(design)[i],"_DOWN")]]
  #upper_res_list[[condicao]] = res_list
  
  if (salvar == T) {
    # # export the full results table
    full_nome = paste0("./DRAresFull_", colnames(design)[i],"_", title, ".tsv")
    write.table(full_results2, file = full_nome, sep = "\t", 
                row.names = F, col.names = T, na = "")
    
    # export the filtered DRA results table by condition
    degs_nome = paste0("./DRAs_",colnames(design)[i],"_", title, ".tsv")
    degs_nome
    write.table(df, file = degs_nome, sep = "\t", 
                row.names = F, col.names = T, na = "")
  }
  
  i=i+1
}

# export the DRAs results table of all conditions together
df3 = l2df(upper_res_list, byrow = F)
degs_all_nome = paste0("./DRAs_ALL_", title, ".tsv")
write.table(df3, file = degs_all_nome, sep = "\t", 
            row.names = F, col.names = T, na = "")


# 
# 
# sequ = seq(from = 1, to = length(upper_res_list), by = 3)
# names(upper_res_list)[sequ]
# 
# sets = which(names(upper_res_list)[sequ] %in% c("AD_TOTAL", "MCI_TOTAL")) 
# sets = which(names(upper_res_list)[sequ] %in% c( "earlyPD_TOTAL", "PD_TOTAL")) 
# sets = which(names(upper_res_list)[sequ] %in% c("rrMS_TOTAL", "spMS_TOTAL", "MS_TOTAL")) 
# 
# names(upper_res_list)[sets*3-2]
# total_sets = sets*3-2
# 
# group_total = c()
# group_up = c()
# group_down = c()
# 
# for (i in total_sets) {
#   group_total = c(group_total, upper_res_list[[i]])
#   group_total = unique(group_total)
#   
#   group_up = c(group_up, upper_res_list[[i+1]])
#   group_up = unique(group_up)
#   
#   group_down = c(group_down, upper_res_list[[i+2]])
#   group_down = unique(group_down)
# }
# 
# #joined_list = list()
# 
# condicao = "MS-group"
# joined_list[[paste0(condicao,"_TOTAL")]] = group_total
# joined_list[[paste0(condicao,"_UP")]] = group_up
# joined_list[[paste0(condicao,"_DOWN")]] = group_down
# 
# joined_table = l2df(joined_list, byrow = F)
# # write.table(joined_table, row.names = F, col.names = T, sep = "\t",
# #             file = "./joined_GSEbatch.tsv", na = "")
# 

## 2.3 Diverging bars ----------------------------------------------------------

tab_net = read.delim(file = "./DRAs_ALL_GSEbatch.tsv", header = T, na.strings = c("", "NA"))
indices = grep(colnames(tab_net), pattern = "TOTAL|cdr00|Sexmale|Age", invert = T) 
tab_net = tab_net[indices]

list_bars = as.list(tab_net)

list_bars_len = lapply(X = list_bars, FUN = function(x){
  soma = sum(complete.cases(x))
  return(soma)
})
df_total = t(as.data.frame(list_bars_len))
colnames(df_total) = c("total")
df = as.data.frame(df_total)
df = rownames_to_column(df, var = "Groups")

df = df %>% dplyr::mutate(direction = gsub(Groups, pattern = ".*_", replacement = ""))
df = df %>% dplyr::mutate(disease = gsub(Groups, pattern = "_.*", replacement = ""))

df = df %>% 
  mutate(Classification = if_else(disease %in% c("MS","rrMS", "spMS", "metaMS"), true = "MS",
                   false =if_else(disease %in% c("PD","earlyPD", "metaPD"), true = "PD",
                   false = if_else(disease %in% c("AD","MCI", "metaAD"), true = "AD",
                   false = "others"))))

df2 = df
colnames(df2) = c("Groups","total","UPDOWN","Names", "Classification")
# adjust names and order
neworder <- rev(c("AD","MCI", "PD","earlyPD", "MS","rrMS","spMS")) #"Sexmale","Age"
df2 <- dplyr::arrange(transform(
  df2, Names=factor(Names,levels=neworder)), Names)
df2$Names = recode(df2$Names, 
                   PD = "Parkinson's Disease (PD)",
                   earlyPD = "Early-stage Parkinson's Disease (earlyPD)",
                   MS = "Multiple Sclerosis (MS)",
                   rrMS = "Relapsing-remitting Multiple Sclerosis (RRMS)",
                   spMS = "Secondary progressive Multiple Sclerosis (SPMS)",
                   AD = "Alzheimer's Disease (AD)",
                   MCI = "Mild Cognitive Impairment (MCI)"
                   #,Sexmale = "Sex", Age = "Age"
)

df2 <- df2 %>% mutate(Classification = 
                        recode(Classification, 
                               AD = "AD group", PD = "PD group", MS = "MS group"))
df2 <- dplyr::arrange(transform(
  df2, Classification = factor(
    Classification,levels=c("AD group", "PD group", "MS group"))), Classification) #, "others"

df2 <- df2 %>% dplyr::mutate(updown = tolower(UPDOWN))

tabela_final = df2

ggplot(data = tabela_final) +
  #downregulated total DEGs, and labels
  geom_bar(data = subset(tabela_final, (UPDOWN == "DOWN")), 
           position = position_dodge(), stat = "identity", 
           aes(fill = UPDOWN, y = -total, x = Names),
           alpha = 1, width = 0.7) +
  geom_text(data = subset(tabela_final, (UPDOWN == "DOWN")),
            aes(y = -total, x = Names, label = total, hjust = 1.2),
            size = 4.5) +
  #upregulated total DEGs, and labels
  geom_bar(data = subset(tabela_final, (UPDOWN == "UP")), 
           position = position_dodge(), stat = "identity", 
           aes(fill = UPDOWN, y = total, x = Names),
           alpha = 1, width = 0.7) +
  geom_text(data = subset(tabela_final, (UPDOWN == "UP")),
            aes(y = total, x = Names, label = total, hjust = - 0.2),
            size = 4.5) +
  theme(legend.position = "bottom", 
        legend.title =  element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_text(color = "black"),
        panel.background = element_rect( fill = "white", color = "#003f5c"),
        panel.grid = element_blank(),
        strip.text = element_text(vjust = 0.7, hjust = 0.5, 
                                  colour = "black", size = 14), 
        strip.background = element_rect(colour="#003f5c", fill="#fff"),
        text = element_text(size = 18, family = "Helvetica")) +
  labs(x = element_blank(), y = element_blank()) +
  coord_flip() +
  facet_grid(rows = vars(Classification), scale = "free_y", space = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) +
  scale_fill_manual(values=c("#2c4875", "#d74a49"),
                    labels = c("Down-regulated DRAs", "Up-regulated DRAs")) 

# save file (!!!)
ggsave(filename = "../../figures/diverging_bar_plots/bars_DRAs_gsebatch.svg", 
       width = 3500, height = 1500,  #width = 2800, height = 1400, 
       units = "px", dpi = 300, device = "svg", scale = 1)
