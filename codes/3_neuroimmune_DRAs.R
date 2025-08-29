# 3. Neuroimmune DRAs ----------------------------------------------------------

## 3.1 Load packages and functions ---------------------------------------------

source("./codes/my_packages.R")
source("./codes/my_functions.R")

# set working directory
setwd("./data/Neuro_tables/GSEbatch/")
getwd()

## 3.2 Network neuroimmune DRAs ------------------------------------------------

# load DRAs table
tab_net = read.delim(file = "./DRAs_ALL_GSEbatch.tsv", header = T, 
                     na.strings = c("", "NA"))

# load neuroimmune classification
neuro = scan(file = "./../proteinatlas_nervous.txt", what = "character")
immune = scan(file = "./../proteinatlas_immune.txt", what = "character")

# only get groups of interest and UP or DOWN cols, not the total 
indices = grep(colnames(tab_net), pattern = "TOTAL|cdr00|Sexmale|Age", invert = T) 
tab_net = tab_net[indices]

tab_net = pivot_longer(tab_net, cols = c(1:ncol(tab_net)), names_to = "Groups", 
                       values_to = "Genes", values_drop_na = T)
tab_net = tab_net %>% dplyr::mutate(direction = gsub(Groups, pattern = ".*_", replacement = ""))
tab_net = tab_net %>% dplyr::mutate(disease = gsub(Groups, pattern = "_.*", replacement = ""))
tab_net = na.omit(tab_net)
tab_net$Group = tab_net$disease
tab_net$Term = tab_net$disease
dfnet = tab_net 

gene_anno = as.data.frame(dfnet$Genes)
colnames(gene_anno) = "Genes"
gene_anno = gene_anno %>% mutate(immune = if_else(Genes %in% immune, T, F))
gene_anno = gene_anno %>% mutate(neuro = if_else(Genes %in% neuro, T, F))
gene_anno = gene_anno %>% mutate(neuroORimmune = (immune|neuro))
gene_anno = gene_anno %>% mutate(neuroANDimmune = (immune&neuro))
gene_anno = gene_anno %>% mutate(System = if_else(neuroANDimmune, "neuroimmune", 
                                 if_else(immune == T & neuro == F, "immune", 
                                 if_else(immune == F & neuro == T, "neuro", 
                                 "othersystem"))))
#summary(gene_anno)
gene_anno2 = unique(gene_anno[,c("Genes", "System")])
dfnet = left_join(x = dfnet, y = gene_anno2, by = "Genes")

#dfnet <- dfnet %>% distinct(Term, Genes, .keep_all = TRUE)
dfnet =  subset(dfnet, System != "othersystem")
net_data <- dfnet %>% distinct(Term, Genes, .keep_all = TRUE)
#write.table(x = dfnet, file = "NI_DRAs.txt", sep = "\t", row.names = F, col.names = T)

# Make network object
net <- network(dfnet[, c("Term", "Genes")], directed = FALSE)

# Make list of attributes
vertex_color <- unique(net_data[, c("Genes", "Group")])
vertex_color <- setNames(vertex_color$Group, vertex_color$Genes)

# add atributes to network object
set.vertex.attribute(net, "Group", vertex_color[network.vertex.names(net)])

unique(dfnet$Term)  # check names

df_color_dicio <- data.frame("Group" = c("AD","earlyPD", "MCI","MS","PD",
                                         "rrMS","spMS", "Other"), 
                             "color" = c("#4a2e66", "#787878", "#8a508f", 
                                         "#ff8531", "grey30","#ffd331","#ffa600", 
                                         "grey70"))

# Define names on a data frame
color_palette <- setNames(df_color_dicio$color, df_color_dicio$Group)  

# Get groups correspondence to nodes
df_color <- data.frame("Group" = get.vertex.attribute(net, "Group"))
#list.vertex.attributes(net)

# Replace NAs for 'Other'
df_color %>% 
  mutate(Group = ifelse(is.na(Group), "Other", Group)) -> df_color
set.vertex.attribute(net, "Group", df_color$Group)
print(df_color)

# Add System attribute to network object (SHAPES). 
# Make attribute list
vertex_shape <- unique(dfnet[, c("Genes", "System")])
vertex_shape <- setNames(vertex_shape$System, vertex_shape$Genes)

# Add atributes to network object
set.vertex.attribute(net, "System", vertex_shape[network.vertex.names(net)])

# Manually define shapes
unique(dfnet$System)  # check names
# Make a data.frame with names and shapes
df_shape_dicio <- data.frame("System" = c("neuroimmune","immune", 
                                          "neuro", "othersystem", "Other"), 
                             "shape" = c(17, 18, 15, 16, 21))

# Define names on data.frame
shape_palette <- setNames(df_shape_dicio$shape, df_shape_dicio$System)  

# Get System corresponence to nodes
df_shape <- data.frame("System" = get.vertex.attribute(net, "System"))
print(df_shape)

# Replace NAs for 'Other'
df_shape %>% 
  mutate(System = ifelse(is.na(System), "Other", System)) -> df_shape
set.vertex.attribute(net, "System", df_shape$System)


# Get top DRAs to add in the network
files <- list.files(path = "./", pattern = "DRAresFull.*.tsv", 
                    full.names = TRUE)

condicoes = c("AD", "MCI", "PD", "earlyPD", "MS", "rrMS", "spMS") 

list_top = list()
for (condicao in condicoes) {
  arquivo = paste0("./DRAresFull_",condicao, "_GSEbatch.tsv") 
  tabela_cond <- read.table(file = arquivo, header = TRUE, sep="\t")
  tabela_cond = tabela_cond %>% 
    dplyr::filter(GeneName %in% c(unique(immune, neuro))) # filter for neuroimmune
  topdegs = tabela_cond %>% 
    dplyr::arrange(P.Value) %>% # arrange by p-value
    slice_head(n = 3) # get top 3 degs
  list_top[[condicao]] = topdegs$GeneName
}

toptargets1 = unique(unlist(list_top, use.names = F))
toptargets2 = data.frame(network.vertex.names(net))
colnames(toptargets2) = "labs" 
rm = which(!(toptargets2$labs %in% toptargets1))
toptargets2[rm,] = ""
toptargets3 = toptargets2$labs

set.seed(53)

ggnet2(net, 
       size = 1,  #size = "degree",
       label = c("AD","earlyPD", "MCI","MS","PD","rrMS","spMS"), #diseases
       alpha = 1, 
       edge.size = 0.07,
       shape = "System",
       shape.palette = shape_palette,
       edge.color = "grey75", #"System",
       shape.legend = TRUE, 
       legend.size = 10,
       color = "Group", 
       palette = c(color_palette)
) +
  theme(legend.position = "right") +
  ggrepel::geom_text_repel(aes(label = toptargets3), max.overlaps = Inf, 
                           bg.color = "white", bg.r = 0.10, size = 3,
                           min.segment.length = 0, point.padding = 0) 

# save file (!!!)
ggsave(filename = "../../figures/network/Network_NI_DRAs.svg",
       device = "svg", width = 8, height = 6, units = "in") 

## 3.3 Barplot neuroimmune DRAs ------------------------------------------------

# run neuroimmune DRAs network section until gene_anno2 is created, line:
# gene_anno2 = unique(gene_anno[,c("Genes", "System")])

gene_anno = gene_anno %>% 
  mutate(System = if_else(neuroANDimmune, "immune and\nnervous antigen", 
                  if_else(immune == T & neuro == F, "immune antigen", 
                  if_else(immune == F & neuro == T, "nervous antigen", 
                          "othersystem"))))

gene_anno2 = unique(gene_anno[,c("Genes", "System")])
dfbars = left_join(x = tab_net, y = gene_anno2, by = "Genes")
dfbars =  subset(dfbars, System != "othersystem")

bars = dfbars |> group_by(disease, System) |> summarise( n = n())
bars$disease = factor(x = bars$disease, levels = c("AD","MCI","PD","earlyPD",
                                                   "MS","rrMS","spMS"))
bars$System = factor(x = bars$System, 
                     levels = c("immune antigen","nervous antigen",
                                "immune and\nnervous antigen"))

ggplot(bars, aes(fill=System, y=n, x=System)) + 
  geom_bar(position="dodge", stat="identity", width = 0.7) + 
  facet_wrap(~ disease, ncol = 2, dir = "h", scales = "free_y", axes = "all") +
  #theme_classic() + 
  scale_fill_manual(values = c("darkgreen", "#d1ec57", "#80d353")) +
  theme(panel.grid = element_blank(), 
        panel.background = element_blank(), #element_rect(colour="black", fill="white"), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line = element_line(color= "black", linewidth = 0.3),
        axis.title = element_blank(),
        strip.text = element_text(size = 10),
        #text = element_text(size = 14),
        strip.background = element_rect(colour="white", fill="white")
  )  +   
  geom_text(aes(label = n), vjust = -0.5, size = 3.4) + 
  scale_y_continuous(expand = expansion(mult = c(0, .3)))

# save file (!!!)
ggsave(filename = "../../figures/barplot/barplot_NI_DRAs.svg", 
       width = 1200, height = 1600, # width = 1400, height = 1700,
       units = "px", dpi = 300, device = "svg", scale = 1)

## 3.4 common & unique DRAs ----------------------------------------------------

# venn diagrams and enrichment of uniquenesses and commonalities between major 
# groups (AD, PD, MS)

#### venn -----

# run NI DRAs network section until dfnet is filtered, line:
# dfnet =  subset(dfnet, System != "othersystem")

dfnet = dfnet %>% 
  mutate(grupao = if_else(Group %in% c("MS","rrMS", "spMS"), true = "MS",
           false =if_else(Group %in% c("PD","earlyPD"), true = "PD",
           false = if_else(Group %in% c("AD","MCI"), true = "AD", 
           false = NA))))
table(dfnet$grupao)

tabAD = dfnet %>% dplyr::filter(grupao == "AD")
tabAD_up = dfnet %>% dplyr::filter(grupao == "AD" & direction == "UP")
tabAD_down = dfnet %>% dplyr::filter(grupao == "AD" & direction == "DOWN")

tabPD = dfnet %>% dplyr::filter(grupao == "PD")
tabPD_up = dfnet %>% dplyr::filter(grupao == "PD" & direction == "UP")
tabPD_down = dfnet %>% dplyr::filter(grupao == "PD" & direction == "DOWN")

tabMS = dfnet %>% dplyr::filter(grupao == "MS")
tabMS_up = dfnet %>% dplyr::filter(grupao == "MS" & direction == "UP")
tabMS_down = dfnet %>% dplyr::filter(grupao == "MS" & direction == "DOWN")

#up
lista_venn_up = list()
lista_venn_up[[1]] = unique(tabAD_up$Genes) 
lista_venn_up[[2]] = unique(tabPD_up$Genes)
lista_venn_up[[3]] = unique(tabMS_up$Genes)

#down
lista_venn_down = list()
lista_venn_down[[1]] = unique(tabAD_down$Genes) 
lista_venn_down[[2]] = unique(tabPD_down$Genes)
lista_venn_down[[3]] = unique(tabMS_down$Genes)

names(lista_venn_up) = c("AD group", "PD group", "MS group")
names(lista_venn_down) = c("AD group", "PD group", "MS group")

# plot venn diagram
ggvenn(lista_venn_up, show_percentage = F, 
       fill_color = c("#d74a49","#ff7964","#ff7541"), #UP
       text_size = 10, text_color = "white", stroke_color = "white") 

# save file (!!!)
ggsave(filename = "../../../figures/Venn/venn_NI_up.svg", bg = "white",
       width = 1800, height = 1600, units = "px", dpi = 300, device = "svg", scale = 1)


ggvenn(lista_venn_down, show_percentage = F, 
       fill_color = c("#00202e","#003f5c","#2c4875"), # DOWN
       text_size = 10, text_color = "white", stroke_color = "white") 

# save file (!!!)
ggsave(filename = "../../../figures/Venn/venn_NI_down.svg", bg = "white",
       width = 1800, height = 1600, units = "px", dpi = 300, device = "svg", scale = 1)


# get lists of common and unique gene
common_up = na.omit(Reduce(intersect, lista_venn_up))
common_down = na.omit(Reduce(intersect, lista_venn_down))

#ad
unique_AD = tabAD %>% 
  dplyr::filter(!(Genes %in% unique(c(tabPD$Genes, tabMS$Genes))))
unique_AD_up = tabAD_up %>% 
  dplyr::filter(!(Genes %in% unique(c(tabPD_up$Genes, tabMS_up$Genes))))
unique_AD_down = tabAD_down %>% 
  dplyr::filter(!(Genes %in% unique(c(tabPD_down$Genes, tabMS_down$Genes))))

#pd
unique_PD = tabPD %>% 
  dplyr::filter(!(Genes %in% unique(c(tabAD$Genes, tabMS$Genes))))
unique_PD_up = tabPD_up %>% 
  dplyr::filter(!(Genes %in% unique(c(tabAD_up$Genes, tabMS_up$Genes))))
unique_PD_down = tabPD_down %>% 
  dplyr::filter(!(Genes %in% unique(c(tabAD_down$Genes, tabMS_down$Genes))))

#ms
unique_MS = tabMS %>% 
  dplyr::filter(!(Genes %in% unique(c(tabPD$Genes, tabAD$Genes))))
unique_MS_up = tabMS_up %>% 
  dplyr::filter(!(Genes %in% unique(c(tabPD_up$Genes, tabAD_up$Genes))))
unique_MS_down = tabMS_down %>% 
  dplyr::filter(!(Genes %in% unique(c(tabPD_down$Genes, tabAD_down$Genes))))

#### enrichment -----

## make list
list_net = list()
list_net[["common_up"]] = unique(common_up)
list_net[["AD_up"]] = unique(unique_AD_up$Genes)
list_net[["PD_up"]] = unique(unique_PD_up$Genes)
list_net[["MS_up"]] = unique(unique_MS_up$Genes)
list_net[["common_down"]] = unique(common_down)
list_net[["AD_down"]] = unique(unique_AD_down$Genes)
list_net[["PD_down"]] = unique(unique_PD_down$Genes)
list_net[["MS_down"]] = unique(unique_MS_down$Genes)

## enrich list
list_net_enrich = lapply(X = list_net, function(genes){
  
  ego <- enrichGO(gene          = genes,
                  OrgDb         = org.Hs.eg.db,
                  keyType = "SYMBOL", 
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  readable      = TRUE)
  ego2 = ego@result
  ego2 = subset(ego2, pvalue < 0.05) #pvalue
  return(ego2)
})

# ad/pd/ms
i = 1
while (i <= length(list_net_enrich)) {
  nome = names(list_net_enrich[i])
  list_net_enrich[[i]]$Groups = nome
  i = i+1
}

tab_join = do.call(what = rbind, list_net_enrich[c(1,5)]) # common updown
table(tab_join$Groups)

tab_join = tab_join %>%
  dplyr::mutate(Groups0 = Groups) %>%
  dplyr::mutate(direction = (gsub(x = Groups, pattern = ".*_", ""))) %>%
  dplyr::mutate(Groups = (gsub(x = Groups, pattern = "_.*", "")))

# write.table(tab_join, file = "./NI_DRAs_venn_enrichment2.tsv", sep = "\t", 
#             row.names = F, col.names = T)

top_names = tab_join %>% 
  dplyr::filter(Groups == "common") %>% 
  group_by(Groups) %>% 
  arrange(pvalue) %>% 
  slice_head(n = 20) 
unique(top_names$Description)

tab_join2 = tab_join %>% 
  dplyr::filter(pvalue < 0.05) %>% 
  dplyr::filter(Description %in% unique(top_names$Description)) %>% 
  dplyr::mutate(Groups = factor(Groups, 
                                levels = c("common","AD","PD","MS")), .keep = "all") %>% 
  dplyr::mutate(direction = factor(direction, 
                                   levels = c("up","down")), .keep = "all") %>% 
  dplyr::arrange(ID)


ggplot(tab_join2) +
  geom_point(aes(x = direction, 
                 y = factor(Description, levels = rev(unique(Description))), 
                 color = direction,  size = -log10(pvalue)
  ))  + 
  theme_bw() +
  theme(axis.text.y = element_text(vjust = 0.5, debug = F, size = 10), 
        axis.text.x = element_text(vjust = 1, size = 10),
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        panel.grid = element_blank(), 
        panel.background = element_blank()) +
  facet_wrap(facets = vars(Groups), nrow = 1) +
  scale_color_manual(values = c("#fa8161ff","#355776ff")) +
  labs(x = element_blank(), y = element_blank()) 

# save file (!!!)
ggsave(filename = "../../figures/dotplots/NI_venn_common.svg", 
       width = 1800, height = 1400, 
       units = "px", dpi = 300, device = "svg", scale = 1)


#### supplementary - top 30 BPs -----

tab_join = do.call(what = rbind, list_net_enrich[c(2,3,4,6,7,8)]) # AD PD MS updown

tab_join = tab_join %>%
  dplyr::mutate(Groups0 = Groups) %>%
  dplyr::mutate(direction = (gsub(x = Groups, pattern = ".*_", ""))) %>%
  dplyr::mutate(Groups = (gsub(x = Groups, pattern = "_.*", "")))

top_names = tab_join %>% 
  dplyr::mutate(
    Groups = factor(Groups, levels = c("common","AD","PD","MS")), .keep = "all") %>% 
  group_by(Groups) %>% 
  arrange(pvalue) %>% 
  dplyr::distinct(ID, Groups, .keep_all = T) %>% 
  slice_head(n = 35) %>% 
  ungroup() %>% dplyr::select(ID, Groups)

tab_join3 = tab_join %>% 
  dplyr::filter(pvalue < 0.05) %>% 
  dplyr::right_join(y = top_names, multiple = "all", 
                    by = join_by(ID, Groups)) %>% 
  dplyr::mutate(
    Groups = factor(Groups, levels = c("common","AD","PD","MS")), .keep = "all") %>% 
  dplyr::mutate(
    direction = factor(direction, levels = c("up","down")), .keep = "all")


ggplot(tab_join3) +
  geom_point(aes(x = direction, #Groups0, 
                 y = Description, 
                 color = Groups0,  size = -log10(pvalue)
  ))  + 
  theme_bw() +
  theme(axis.text.y = element_text(vjust = 0.5, debug = F, size = 10), 
        axis.text.x = element_text(vjust = 1, size = 10),
        strip.background = element_blank(), 
        strip.text = element_text(size = 12),
        panel.grid = element_blank(), #element_line(),
        panel.background = element_blank(),
        plot.title = element_text(size = 10)) +
  facet_wrap(facets = vars(Groups), nrow = 2, scales = "free", drop = T) +
  scale_color_manual(values = c("#596d77ff","#de8887ff", #AD
                                "#7587a5ff","#ffa583ff",  #MS
                                "#7f9eadff","#ffa89bff" #PD
  )) +
  labs(x = element_blank(), y = element_blank()) 

# save file (!!!)
ggsave(filename = "../../figures/dotplots/NI_venn_unique.svg", 
       width = 5000, height = 4600, units = "px", dpi = 300, device = "svg", scale = 1)
