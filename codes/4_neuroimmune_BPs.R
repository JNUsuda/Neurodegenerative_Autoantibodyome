# 4. Neuroimmune BPs -----------------------------------------------------------

## 4.1 Load packages and functions ---------------------------------------------

source("./codes/my_packages.R")
source("./codes/my_functions.R")

# Load classification lists:
# neuroimmune classification - genes:
neuro = scan(file = "./data/proteinatlas_nervous.txt", what = "character")
immune = scan(file = "./data/proteinatlas_immune.txt", what = "character")

# neuroimmune classification - BPs:
immune_BPs = read.table(file = "./data/immune-system_UBERON_0002405_GO-BPs.txt", header = T)
immune_BPs = immune_BPs$id
neuro_BPs = read.table(file = "./data/nervous-system_UBERON_0001016_GO-BPs.txt", header = T)
neuro_BPs = neuro_BPs$id

# set working directory
setwd("./data/GSEbatch/")
getwd()

## 4.2 Appyter neuroimmune BPs -------------------------------------------------

# Load clustering table from Appyter
data.cluster = read.table(file = "../appyter_base.txt", header = T, sep = "\t", 
                          quote = "", comment.char = "")
# Create GO IDs column
data.cluster <- data.cluster %>%
  mutate(ID = str_extract(term, "(GO:[:digit:]{7})"))

tab_net = read.delim(file = "./DRAs_ALL_GSEbatch.tsv", header = T, 
                     na.strings = c("", "NA"))

# only get groups of interest and UP or DOWN cols, not the total 
indices = grep(colnames(tab_net), pattern = "TOTAL|cdr00|Sexmale|Age", invert = T) 
tab_net = tab_net[indices]
list_net = as.list(tab_net)
list_net = sapply(X = list_net, na.omit)

list_net_enrich1 = lapply(X = list_net, FUN = function(genes){
  
  ego <- enrichGO(gene          = genes,
                  OrgDb         = org.Hs.eg.db,
                  keyType = "SYMBOL", 
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 1,
                  qvalueCutoff  = 1,
                  readable      = TRUE)
  ego2 = ego@result
  return(ego2)
})

list_net_enrich2 = list()
i = 1
while (i <= length(list_net_enrich1)) {
  # add group column
  nome = names(list_net_enrich1[i])
  df = list_net_enrich1[[i]]
  df$Group_original = nome
  df = df %>% mutate(Groups = gsub(Group_original, pattern = "_.+", replacement = ""))
  df = df %>% mutate(Direction = gsub(Group_original, pattern = ".+_", replacement = ""))
  
  # merge with 'base' table from appyter
  l = full_join(data.cluster, df[, c("ID", "Description", "pvalue", "p.adjust", 
                                     "qvalue", "Groups","Direction" )], by = "ID")
  l = na.omit(l) #remove NAs
  l$cluster = factor(l$cluster) # set cluster as factors
  # set NAs to p value 1
  l$pvalue[is.na(l$pvalue)] <- 1
  l$p.adjust[is.na(l$p.adjust)] <- 1
  l$qvalue[is.na(l$qvalue)] <- 1
  l$qvalue[is.na(l$qvalue)] <- 1
  # significant p-values
  l$sig = ifelse(l$pvalue < 0.05, 1, 0)
  # get number os significant BPs by cluster
  l = l |>
    group_by(cluster) |>
    mutate(filtered = sum(sig))
  
  # select neuroimmune BPs 
  l = l %>% mutate(immune = if_else(ID %in% immune_BPs, T, F))
  l = l %>% mutate(neuro = if_else(ID %in% neuro_BPs, T, F))
  
  l = l %>% mutate(neuroORimmune = (immune|neuro))
  l = l %>% mutate(neuroANDimmune = (immune&neuro))
  l = l %>% mutate(System = if_else(neuroANDimmune, "immune\nand nervous", 
                            if_else(immune == T & neuro == F, "immune", 
                            if_else(immune == F & neuro == T, "nervous", 
                            false = "othersystem"))))
  list_net_enrich2[[i]] = l
  i = i+1
}

tab_cluster = do.call(what = rbind, list_net_enrich2)
tab_cluster$Groups = factor(tab_cluster$Groups, 
                            levels = c("AD", "MCI", "PD", "earlyPD", 
                                       "MS", "rrMS", "spMS"))


# plot 
ggplot(tab_cluster, aes(x, y)) +
  geom_point( 
    shape = 16,
    aes(color = System, size = -log10(pvalue))) +
  theme(panel.grid = element_blank(), 
        axis.ticks = element_line(linewidth = 1),
        axis.ticks.length = unit(0.2, "cm"),
        panel.background = element_rect(colour="black", fill="white"), 
        strip.text = element_text(vjust = 0.7, hjust = 0.5, colour = "black"), 
        strip.background = element_rect(colour="white", fill="white"),
        text = element_text(size = 36, family = "Helvetica")
  ) +
  labs(x = "UMAP1", y = "UMAP2") +
  # geom_label_repel(data = data_all,
  #                  mapping = aes(x = x_mean, y = y_mean,
  #                                label = paste(cluster)),
  #                  colour = 'black',
  #                  size = 6,
  #                  nudge_x = -1,  # Ajuste manual de deslocamento horizontal
  #                  nudge_y = 1.5)+
  #Ajuste manual de deslocamento vertical +
  facet_wrap(~ Groups, nrow = 1, dir = "h") +
  scale_color_manual(values = c("darkgreen","#80d353", "#d1ec57", "grey80")) #  "#80d353"


# save file (!!!)
ggsave(filename = "../../figures/appyter/appyter_NI.svg", width = 10400, 
       height = 1800, units = "px", dpi = 300, device = "svg", scale = 1)

## 4.3 Barplot neuroimmune BPs -------------------------------------------------

list_net_enrich3 = list()
i = 1
while (i <= length(list_net_enrich1)) {
  # add group column
  nome = names(list_net_enrich1[i])
  df = list_net_enrich1[[i]]
  df$Group_original = nome
  df = df %>% mutate(Groups = gsub(Group_original, pattern = "_.+", replacement = ""))
  df = df %>% mutate(Direction = gsub(Group_original, pattern = ".+_", replacement = ""))
  # filter for significant p-values
  l = df %>% subset(pvalue < 0.05)
  # nervous/immune classification
  l = l %>% mutate(immune = if_else(ID %in% immune_BPs, T, F))
  l = l %>% mutate(neuro = if_else(ID %in% neuro_BPs, T, F))
  l = l %>% mutate(neuroORimmune = (immune|neuro))
  l = l %>% mutate(neuroANDimmune = (immune&neuro))
  l = l %>% mutate(System = if_else(neuroANDimmune, "immune and\nnervous BP", 
                            if_else(immune == T & neuro == F, "immune BP", 
                            if_else(immune == F & neuro == T, "nervous BP", 
                            false = "othersystem"))))

  list_net_enrich3[[i]] = l
  i = i+1
}

tab_cluster2 = do.call(what = rbind, list_net_enrich3)
tab_cluster2$Groups = factor(tab_cluster2$Groups, 
                             levels = c("AD", "MCI", "PD", "earlyPD", 
                                        "MS", "rrMS", "spMS"))

bars = tab_cluster2 |> 
  subset(System != "othersystem") |>
  group_by(Groups, System) |> 
  summarise( n = n())

bars$Groups = factor(x = bars$Groups, levels = c("AD","MCI","PD","earlyPD",
                                                 "MS","rrMS","spMS"))
bars$System = factor(x = bars$System, levels = c("immune BP","nervous BP",
                                                 "immune and\nnervous BP"))

ggplot(bars, aes(fill=System, y=n, x=System)) + 
  geom_bar(position="dodge", stat="identity", width = 0.7) + 
  facet_wrap(~ Groups, ncol = 2, dir = "h", scales = "free_y", axes = "all") +
  #theme_classic() + 
  scale_fill_manual(values = c("darkgreen", "#d1ec57","#80d353")) +
  theme(panel.grid = element_blank(), 
        panel.background = element_blank(), #element_rect(colour="black", fill="white"), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line = element_line(color= "black", linewidth = 0.3),
        axis.title = element_blank(), 
        strip.text = element_text(size = 12),
        strip.background = element_rect(colour="white", fill="white"), 
        legend.position = "right"
  ) +   
  geom_text(aes(label = n), vjust = -0.5) + 
  scale_y_continuous(expand = expansion(mult = c(0, .3)))

# save file (!!!)
ggsave(filename = "../../../figures/barplot/barplot_NI_BPs.svg", width = 1200,
       height = 1600, units = "px", dpi = 300, device = "svg", scale = 1)

## 4.4 Dotplots top neuroimmune BPs --------------------------------------------

# select immune BPs
tab_cluster_N = tab_cluster2 %>% subset(ID %in% immune_BPs)

# repeat for nervous
#tab_cluster_N = tab_cluster2 %>% subset(ID %in% neuro_BPs)

# get the top unique BPs by disease (by p-value)
top_unique = tab_cluster_N %>%
  group_by_at(vars(ID)) %>%
  mutate(num_rows = sum(n())) %>%
  filter(num_rows == 1) %>% 
  ungroup() %>% group_by(Groups) %>% 
  arrange(pvalue) %>% 
  slice_head(n = 1) 
unique(top_unique$Description)

# get the top shared BPs by disease (by p-value)
top_common = tab_cluster_N %>%
  filter(System == "immune BP") %>% #filter unique to system
  group_by_at(vars(ID)) %>%
  mutate(num_rows = sum(n())) %>%
  filter(num_rows > 1) %>% 
  ungroup() %>% 
  arrange(pvalue) %>% 
  distinct(ID, .keep_all = T) %>% 
  slice_head(n = 13) 
unique(top_common$Description)

top_tab_cluster_N = tab_cluster2 %>% 
  subset(ID %in% c(top_unique$ID, top_common$ID)) %>% 
  arrange(Groups) %>% 
  dplyr::mutate(Description = factor(Description, 
                              levels=rev(c(unique(top_unique$Description), 
                                           unique(top_common$Description)))))

ggplot(top_tab_cluster_N) +
  geom_point(aes(x = Groups, y = Description, color = Groups, 
                 size = -log10(pvalue)), shape = 19) + 
  theme_bw() +
  theme(axis.text.y = element_text(vjust = 0.5, debug = F, size = 10), 
        axis.text.x = element_text(vjust = 1, hjust = 1, size = 10, angle = 45),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        panel.grid = element_blank(),
        panel.background = element_blank()) +
  scale_color_manual(values = c("#4a2e66","#8a508f","grey30","#787878",
                                "#ff8531","#ffd331","#ffa600"))

# save file (!!!)
ggsave(filename = "../../figures/dotplot/hm_I_BPs.png", width = 2000, 
       height = 1300, units = "px", dpi = 300, device = "png", scale = 1)

## 4.5 Network neuroimmune BPs -------------------------------------------------

tab_net = read.delim(file = "./DRAs_ALL_GSEbatch.tsv", header = T, 
                     na.strings = c("", "NA"))

# only get groups of interest and UP or DOWN cols, not the total 
indices = grep(colnames(tab_net), pattern = "TOTAL|cdr00|Sexmale|Age", invert = T) 
tab_net = tab_net[indices]
list_net = as.list(tab_net)
list_net = sapply(X = list_net, na.omit)

list_net_enrich = lapply(X = list_net, function(genes){
  ego <- enrichGO(gene          = genes,
                  OrgDb         = org.Hs.eg.db,
                  keyType = "SYMBOL", 
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  readable      = TRUE)
  ego2 = ego@result
  ego2 = subset(ego2, pvalue < 0.05) 
  return(ego2)
})

i = 1
while (i <= length(list_net_enrich)) {
  nome = names(list_net_enrich[i])
  list_net_enrich[[i]]$Groups = nome
  i = i+1
}

tab_net = do.call(what = rbind, list_net_enrich)
tab_net = tab_net %>% dplyr::filter(pvalue < 0.01)

tab_net = tab_net %>% 
  dplyr::mutate(direction = gsub(Groups, pattern = ".*_", replacement = "")) %>% 
  dplyr::mutate(disease = gsub(Groups, pattern = "_.*", replacement = ""))
tab_net = na.omit(tab_net)
#write.table(tab_net, file = "./clusterprofiler_enrichment-res.txt", row.names = F, sep = "\t")
tab_net$Group = tab_net$disease
tab_net$Term = tab_net$ID
tab_net$Genes = tab_net$geneID

# BP classification by system
tab_net = tab_net %>% 
  subset(Term %in% unique(c(neuro_BPs, immune_BPs))) %>% 
  mutate(Term_immune = if_else(Term %in% immune_BPs, T, F)) %>%
  mutate(Term_neuro = if_else(Term %in% neuro_BPs, T, F)) %>% 
  mutate(Term_System = 
           if_else(Term_immune == T & Term_neuro == F, "immune_BP", 
           if_else(Term_immune == F & Term_neuro == T, "neuro_BP", 
           if_else(Term_immune == T & Term_neuro == T, "neuroimmune_BP", 
           false = "other_BP"))))
dfnet = tab_net 

dfnet <- dfnet %>%
  separate_rows(Genes, sep = "/")

gene_anno = as.data.frame(dfnet$Genes)
colnames(gene_anno) = "Genes"
gene_anno = gene_anno %>% 
  mutate(immune = if_else(Genes %in% immune, T, F)) %>% 
  mutate(neuro = if_else(Genes %in% neuro, T, F)) %>% 
  mutate(neuroORimmune = (immune|neuro)) %>% 
  mutate(neuroANDimmune = (immune&neuro)) %>% 
  mutate(System = if_else(neuroANDimmune, "neuroimmune", 
                  if_else(immune == T & neuro == F, "immune", 
                  if_else(immune == F & neuro == T, "neuro", 
                  false = "othersystem"))))

gene_anno2 = unique(gene_anno[,c("Genes", "System")])
dfnet = left_join(x = dfnet, y = gene_anno2, by = "Genes")

#dfnet =  subset(dfnet, System != "othersystem")
dfnet <- dfnet %>% mutate(cTerm = paste(Term, Group, sep = " "))
dfnet <- dfnet %>% mutate(cGenes = paste(Genes, Group, sep = " "))
net_data <- dfnet %>% distinct(cTerm, cGenes, .keep_all = TRUE)

# make network object
net <- network(net_data[, c("cTerm", "cGenes")], directed = FALSE)
plot(net)

# Make list of attributes
vertex_color <- unique(net_data[, c("cGenes", "Group")])
vertex_color <- setNames(vertex_color$Group, vertex_color$cGenes)

# Add attributes to network object
set.vertex.attribute(net, "Group", vertex_color[network.vertex.names(net)])
unique(dfnet$Group) 

df_color_dicio <- data.frame("Group" = c("AD","earlyPD", "MCI","MS","PD",
                                         "rrMS","spMS", "BPs"), 
                             "color" = c("#4a2e66", "#787878", "#8a508f", "#ff8531", 
                                         "grey30","#ffd331","#ffa600", "grey70"))

# Define names on data.frame
color_palette <- setNames(df_color_dicio$color, df_color_dicio$Group) 

# Get Group corresponding nodes
df_color <- data.frame("Group" = get.vertex.attribute(net, "Group"))

#Replace NA for "Other"
df_color %>% 
  mutate(Group = ifelse(is.na(Group), "BPs", Group)) -> df_color
set.vertex.attribute(net, "Group", df_color$Group)
print(df_color)

#--------  Add System attribute to network object (SHAPES)
# Make list of attributes
vertex_shape <- unique(dfnet[, c("cGenes", "System")])
vertex_shape <- setNames(vertex_shape$System, vertex_shape$cGenes)

# Add attributes to network object
set.vertex.attribute(net, "System", vertex_shape[network.vertex.names(net)])

# Define shapes 
unique(dfnet$System)  
# Create a data.frame with names defining the shapes
df_shape_dicio <- data.frame("System" = c("neuroimmune","immune",
                                          "neuro", "othersystem", "Other"),
                             "shape" = c(17, 18, 15, 16, 21))

# Define names on data.frame
shape_palette <- setNames(df_shape_dicio$shape, df_shape_dicio$System) 

# Get System correspondence to nodes
df_shape <- data.frame("System" = get.vertex.attribute(net, "System"))
print(df_shape)
# Replace NA for "Other"
df_shape %>% 
  mutate(System = ifelse(is.na(System), "Other", System)) -> df_shape
set.vertex.attribute(net, "System", df_shape$System)


# get the top BPs by disease (by p-value)
# for each subgroup, 1 immune or neuroimune BP, and 1 nervous BP

top_BPs = tab_net[0,]
for (group in c("AD","MCI","PD","earlyPD","MS","rrMS","spMS")) {
  # immune
  i = 1
  bp_i = tab_net %>%
    filter((Term_System == "immune_BP" | Term_System == "neuroimmune_BP") & Group == group) %>%
    arrange(pvalue) 
  while (bp_i[i,1] %in% top_BPs$ID) { i = i+1 } # while the top BP was already included, check if next is different
  bp_i2 = bp_i[i,] 
  top_BPs = rbind(top_BPs, bp_i2)
  # neuro
  j = 1 
  bp_n = tab_net %>%
    filter(Term_System == "neuro_BP" & Group == group) %>%
    arrange(pvalue) 
  while (bp_n[j,1] %in% top_BPs$ID) { j = j+1 } # while the top BP was already included, check if next is different
  bp_n2 = bp_n[j,]
  top_BPs = rbind(top_BPs, bp_n2)
}

# top BPs
highlight_BPs = top_BPs %>% 
  dplyr::mutate(cTerm = paste(Term, Group))
issoai <- data.frame(network.vertex.names(net))
colnames(issoai) = "cTerm"
issoai2 = left_join(issoai, highlight_BPs)
issoai3 = issoai2$Description

# size by count for BPs, size 1 for genes
issoai4 = issoai %>% 
  mutate(class = if_else(cTerm %in% dfnet$cTerm, true = "cTerm", false = "cGene"))
anyNA(issoai4)
degrees = dfnet %>% dplyr::select("Count", "cTerm") %>% 
  distinct(cTerm, .keep_all = T)
issoai4 = left_join(issoai, degrees) %>% 
  mutate_if(is.numeric, coalesce, 1) # if is NA (is a cGene), replace by 1

set.seed(79) 
ggnet2(net, 
           size = issoai4$Count,
           label = c("AD","earlyPD", "MCI","MS","PD","rrMS","spMS"), #diseases
           label.size = 4, 
           repel = TRUE,
           alpha = 0.8,
           edge.size = 0.3, edge.alpha = 0.4,
           shape = "System", 
           shape.palette = shape_palette, 
           size.legend.title = "Degree",  
           size.legend.position = "bottom",
           edge.color = "gray70",
           shape.mapping = aes(shape = degree), 
           shape.legend = T,
           size.legend = F,
           legend.size = 10,
           color = "Group", 
           palette = color_palette
) +
  theme(legend.position = "right") + 
  ggrepel::geom_text_repel(aes(label = issoai3), max.overlaps = Inf, 
                           bg.color = "white", bg.r = 0.10, size = 3,
                           min.segment.length = 0, point.padding = 0) 

# save file (!!!)
ggsave(filename = "../../figures/network/Network_NI_BPs.svg",
       device = "svg", width = 11, height = 9, units = "in")  

## 4.5 Semantic Similarity -----------------------------------------------------

#BiocManager::install("GOSemSim")
#library(GOSemSim)

# repeat for tab_net_001, with enrichment results for adj. p-val<0.01 DRAs
# and tab_net_005, with enrichment results for adj. p-val<0.05 DRAs
# to compare both enrichment results

{tab_net_005 = read.delim(file = "./clusterprofiler_enrichment-res.txt", 
                         header = T, na.strings = c("", "NA"))

# BP classification by system
tab_net_005 = tab_net_005 %>% subset(ID %in% unique(c(neuro_BPs, immune_BPs))) #uberon terms
tab_net_005 = tab_net_005 %>% mutate(Term_immune = if_else(ID %in% immune_BPs, T, F))
tab_net_005 = tab_net_005 %>% mutate(Term_neuro = if_else(ID %in% neuro_BPs, T, F))
tab_net_005 = tab_net_005 %>% 
  mutate(Term_System = 
           if_else(Term_immune == T & Term_neuro == F, "immune_BP", 
           if_else(Term_immune == F & Term_neuro == T, "neuro_BP", 
           if_else(Term_immune == T & Term_neuro == T, "neuroimmune_BP", "other_BP"))))
tab_net_005 =  subset(tab_net_005, Term_System != "other_BP")
}

# semantic similarity analysis

# load GO data
hsGO <- godata(annoDb = 'org.Hs.eg.db', ont="BP")

# semantic similarity for all results
mgoSim(unique(tab_net_005$ID), unique(tab_net_001$ID), semData=hsGO, 
       measure="Wang", combine="rcmax")

# semantic similarity by subgroup
for (item in c("AD","MCI", "PD", "earlyPD", "MS", "spMS", "rrMS")) {
  print(item)
  tab_net_005_sub = tab_net_005 %>% dplyr::filter(disease == item)
  tab_net_001_sub = tab_net_001 %>% dplyr::filter(disease == item)
  print(mgoSim(tab_net_005_sub$ID, tab_net_001_sub$ID, semData=hsGO, 
               measure="Wang", combine="rcmax"))
}

# similarity matrix
mat = termSim(tab_net_005$ID, tab_net_001$ID, semData=hsGO, method = "Wang")

# save similarity  matrix (!!!)
write.table(mat, file = "./GOsimilaritymatrix_005_001.txt", sep = "\t", 
            col.names = T, row.names = T)
