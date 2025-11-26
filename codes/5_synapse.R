# 5. Synapse -------------------------------------------------------------------

## 5.1 Load packages and functions ---------------------------------------------

source("./codes/my_packages.R")
source("./codes/my_functions.R")

# Load classification lists:
# genes classification:
neuro = scan(file = "./data/proteinatlas_nervous.txt", what = "character")
immune = scan(file = "./data/proteinatlas_immune.txt", what = "character")
# download from the SynGO website (https://syngoportal.org/):
syngo = read_xlsx(path = "./data/syngo_genes.xlsx") 
syngo = syngo$hgnc_symbol 

# set working directory
setwd("./data/GSEbatch/")
getwd()

## 5.2 Enrichment synapse ------------------------------------------------------

tab_net = read.delim(file = "./DRAs_ALL_GSEbatch.tsv", header = T, 
                     na.strings = c("", "NA"))

# only get groups of interest and UP or DOWN cols, not the total 
indices = grep(colnames(tab_net), pattern = "TOTAL|cdr00|Sexmale|Age", invert = T) 
tab_net = tab_net[indices]
list_ni = as.list(tab_net)
list_ni = sapply(X = list_ni, na.omit)
list_ni = sapply(list_ni, FUN = function(i){i = i[(i %in% immune) | (i %in% neuro)]})
df_ni = l2df(list_ni, byrow = F)

## Enrichr 
# The genes were manually submitted to Enrichr for the enrichment analysis
# The GO Biological Process tables were downloaded,
# named such as "SynGO_2024_table_AD_UP.txt",
# and saved in the folder "EnrichR_SynGO.
# Only the up and down enrichments were selected.
# To filter and arrange the downloaded tables in one excel file:

Enrichrtxt <- list.files(path = "./EnrichR_SynGO", pattern = "\\.txt", 
                         full.names = TRUE)
remove <- which(grepl(Enrichrtxt, pattern = "TOTAL.txt") == T)
Enrichrtxt = Enrichrtxt[-remove]

# export tables in the list as sheets in a xlsx
lista  = lapply(Enrichrtxt, function(i){
  tabela <- read.table(file = i, header = TRUE, sep="\t", dec=".", quote = "")
  exclude <- grep(x = colnames(tabela), pattern = "Old")
  tabela2 <- tabela[-(exclude)]
  tabela2 <- tabela2 %>% 
    dplyr::mutate(Description = (gsub(x = Term, pattern = " \\(.*", "")))
  tabela2 <- tabela2 %>% dplyr::relocate(Description, .after = Term)
  tabela2 <- tabela2 %>% dplyr::mutate(TermID = (gsub(x = Term, pattern = ".*\\(|\\) .*", "")))
  tabela2 <- tabela2 %>% dplyr::relocate(TermID, .after = Description)
  tabela2 <- tabela2 %>% dplyr::mutate(GOType = (gsub(x = Term, pattern = ".* ", "")))
  tabela2 <- tabela2 %>% dplyr::relocate(GOType, .after = TermID)
  
  tabela2 <- tabela2 %>% dplyr::mutate(Group = gsub(x = i, pattern = ".*2024_table_|\\.txt", replacement = ""))
  tabela2 <- tabela2 %>% dplyr::mutate(direction = gsub(x = Group, pattern = ".*_", replacement = ""))
  tabela2 <- tabela2 %>% dplyr::mutate(Group = gsub(x = Group, pattern = "_.*", replacement = ""))
  
  tabela2 <- subset(tabela2, tabela2$P.value < 0.05)
  return(tabela2)
})

enrichrdf = Reduce(rbind, lista)

tab_cluster_N = enrichrdf
tab_cluster_N3 <- tab_cluster_N %>%
  separate_rows(Genes, sep = ";")

# order BPs and CCs
bps = unique(tab_cluster_N3$Description)
bps = rev(c(
  "Synaptic Vesicle Priming",
  "Regulation Of Synaptic Vesicle Priming",
  "Modulation Of Chemical Synaptic Transmission",
  "Synaptic Vesicle Exocytosis",
  "Regulation Of Synaptic Vesicle Exocytosis",
  "Synaptic Vesicle Endocytosis",
  "Neurotransmitter Uptake",
  "Synaptic Vesicle To Endosome Fusion",
  "Regulation Of Synapse Disassembly",
  "Maintenance Of Synapse Structure",
  "Regulation Of Presynaptic Cytosolic Calcium Levels",
  "Presynaptic Actin Cytoskeleton Organization",
  "Presynaptic Modulation Of Chemical Synaptic Transmission",
  "Trans-Synaptic Signaling, Modulating Synaptic Transmission",
  "Regulation Of Synaptic Membrane Adhesion",
  "Regulation Of Synaptic Vesicle Fusion To Presynaptic Active Zone Membrane",
  "Presynaptic Dense Core Vesicle Exocytosis",
  "Regulation Of Presynaptic Dense Core Vesicle Exocytosis",
  "Voltage-Gated IC Activity Involved In Regulation Of Presynaptic Membrane Potential",
  "Ligand-Gated IC Activity Involved In Regulation Of Presynaptic Membrane Potential",
  "Modification Of Postsynaptic Structure",
  "Exocytic Insertion Of Neurotransmitter Receptor To Postsynaptic Membrane",
  "Postsynaptic Modulation Of Chemical Synaptic Transmission",
  "Modification Of Postsynaptic Actin Cytoskeleton",
  "Regulation Of Modification Of Postsynaptic Actin Cytoskeleton",
  "Voltage-Gated IC Activity Involved In Regulation Of Postsynaptic Membrane Potential",
  "Transmitter-Gated IC Activity Involved In Regulation Of Postsynaptic Membrane Potential"
))

ccs = rev(c(
  "Presynapse"                                                
  , "Presynaptic Cytosol" 
  , "Integral Component Of Presynaptic Active Zone Membrane"    
  , "Extrinsic Component Of Presynaptic Membrane"               
  , "Extrinsic Component Of Presynaptic Endocytic Zone Membrane"
  , "Postsynapse"                                               
  , "Postsynaptic Density"                                      
  , "Postsynaptic Density, Intracellular Component"             
  , "Integral Component Of Postsynaptic Membrane"               
  , "Integral Component Of Postsynaptic Density Membrane"       
  , "Anchored Component Of Postsynaptic Density Membrane"       
  , "Extrinsic Component Of Postsynaptic Density Membrane"      
))

tab_cluster_N3$Description = factor(tab_cluster_N3$Description, levels = c(ccs,bps))
tab_cluster_N3$Group = factor(tab_cluster_N3$Group, 
                              levels = c("AD","MCI", "PD", "earlyPD", "MS",
                                         "rrMS", "spMS"))

ggplot(tab_cluster_N3) +
  geom_point(aes(x = Group, y = Description, fill = -log10(P.value), 
                 size = parse_ratio(Overlap)
  ),  shape = 23, 
  color = "black", stroke = 0.1) + 
  theme_bw() +
  theme(axis.text.y = element_text(vjust = 0.5,  size = 10), 
        axis.text.x = element_text(hjust = 1, vjust = 1, size = 10, angle = 45),
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        panel.grid.major.y  = element_blank(), 
        panel.background = element_blank(),
        panel.border = element_blank(), 
        axis.ticks = element_blank(),
        axis.title = element_blank()) + 
  scale_shape_manual(values=c(25,24)) +
  scale_fill_gradient2( high = "#ffa600", low = "#ffff62",
                        mid = "#ffd331", midpoint = 3, limits = c(0, 11))  

# save file (!!!)
ggsave(filename = "../../figures/dotplots/SynGO.svg", width = 2600, 
       height = 2300, units = "px", dpi = 300, device = "svg", scale = 1)

## 5.3 Barplot synapse DRAs and terms ------------------------------------------

tab_net = read.delim(file = "./DRAs_ALL_GSEbatch.tsv", header = T, 
                     na.strings = c("", "NA"))

# only get groups of interest and UP or DOWN cols, not the total 
indices = grep(colnames(tab_net), pattern = "TOTAL|cdr00|Sexmale|Age", invert = T) 
tab_net = tab_net[indices]

tab_net = pivot_longer(tab_net, cols = c(1:ncol(tab_net)), names_to = "Groups", 
                       values_to = "Genes", values_drop_na = T)
tab_net = tab_net %>% 
  dplyr::mutate(direction = gsub(Groups, pattern = ".*_", replacement = "")) %>% 
  dplyr::mutate(disease = gsub(Groups, pattern = "_.*", replacement = ""))
tab_net = na.omit(tab_net)
dfnet = tab_net 

gene_anno = as.data.frame(dfnet$Genes)
colnames(gene_anno) = "Genes"
gene_anno = gene_anno %>% mutate(immune = if_else(Genes %in% immune, T, F))
gene_anno = gene_anno %>% mutate(neuro = if_else(Genes %in% neuro, T, F))
gene_anno = gene_anno %>% mutate(neuroORimmune = (immune|neuro))
gene_anno = gene_anno %>% mutate(neuroANDimmune = (immune&neuro))
gene_anno = gene_anno %>% mutate(System = if_else(neuroANDimmune, "immune\nand nervous", 
                                 if_else(immune == T & neuro == F, "immune", 
                                 if_else(immune == F & neuro == T, "nervous", 
                                 false = "othersystem"))))
gene_anno_syngo = gene_anno %>% 
  mutate(System_syngo = if_else(neuroORimmune, "immune\nor nervous", "othersystem"))
gene_anno_syngo = gene_anno_syngo %>% 
  subset(neuroORimmune == T)
gene_anno_syngo = gene_anno_syngo %>% 
  mutate(synapse = if_else(Genes %in% syngo, "synapse", "othersystem")) #syngo
gene_anno_syngo = gene_anno_syngo %>% 
  pivot_longer(cols = c("System_syngo", "synapse"), 
               names_to = "aa", values_to = "System_syngo") #syngo

gene_anno_syngo2 = unique(gene_anno_syngo[,c("Genes", "System_syngo")])
dfnet_syngo = left_join(x = dfnet, y = gene_anno_syngo2, by = "Genes")
dfnet_syngo =  subset(dfnet_syngo, System_syngo != "othersystem")

bars = dfnet_syngo |> 
  subset(System_syngo == "synapse") |>
  group_by(disease, System_syngo) |> summarise( n = n())

# number of enriched CCs and BPs
bars2 = enrichrdf |> group_by(Group, GOType) |> summarise( n = n())
colnames(bars2) = colnames(bars)
bars3 = rbind(bars, bars2)
bars3$System_syngo = factor(bars3$System_syngo, 
                            levels = c("immune\nor nervous", "synapse","BP", "CC" ))
bars3$disease = factor(x = bars3$disease, 
                       levels = c("AD","MCI","PD","earlyPD","rrMS","spMS","MS"))

ggplot(bars3, aes(fill=System_syngo, y=n, x=System_syngo)) + 
  #geom_bar(position="dodge", stat="identity", width = 0.7) + 
  geom_col_pattern(aes(pattern = System_syngo, 
                       pattern_type = System_syngo, 
                       pattern_density = System_syngo,
                       pattern_angle = System_syngo), 
                   fill = '#fbec6a', colour  = 'white',width = 0.7,
                   pattern_alpha = 0.5) +
  scale_pattern_manual(values=c(NA, 'stripe', 'stripe')) +
  scale_pattern_angle_manual(values = c(synapse = 0, CC=90, BP=30))+
  scale_pattern_density_manual(values = c(synapse = 0, CC=0.3, BP=0.1))+
  scale_pattern_type_manual(values=c(NA, NA, NA)) + 
  facet_wrap(~ disease, ncol = 2, dir = "h", scales = "fixed", axes = "all") +
  theme(panel.grid = element_blank(), 
        panel.background = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line = element_line(color= "black", linewidth = 0.3),
        axis.title = element_blank(),
        strip.background = element_blank()
  )  +   
  geom_text(aes(label = n), vjust = -0.5) + 
  scale_y_continuous(expand = expansion(mult = c(0, .2))) + 
  scale_x_discrete(labels= c("synapse antigens","synapse BPs","synapse CCs"))

# save file (!!!)
ggsave(filename = "../../figures/barplot/barplot_Syngo.svg", width = 1000, 
       height = 1500, units = "px", dpi = 300, device = "svg", scale = 1)

## 5.4 Network synapse DRAs and terms ------------------------------------------

tab_cluster_N = enrichrdf

tab_cluster_N$Term = tab_cluster_N$Description #TermID
# filter by p-value
top_names = tab_cluster_N %>% 
  arrange(P.value) %>% 
  slice_head(n = 18) # length(unique(top_names$Description)) = 15
unique(top_names$Description) 

terms_n = length(unique(tab_cluster_N$Description))

tab_cluster_N = tab_cluster_N %>% subset(Description %in% top_names$Description)

tab_cluster_N <- tab_cluster_N %>%
  separate_rows(Genes, sep = ";")

tab_cluster_N = tab_cluster_N %>% subset(Genes %in% syngo)
tab_cluster_N$Term = tab_cluster_N$Description

dfnet = tab_cluster_N
dfnet <- dfnet %>% mutate(cTerm = paste(Term, Group, sep = " "))
dfnet <- dfnet %>% mutate(cGenes = paste(Genes, Group, sep = " "))
net_data <- dfnet %>% distinct(Term, cGenes, .keep_all = TRUE)

net <- network(net_data[, c("Term", "cGenes")], directed = FALSE)
plot(net)

# Create attribute list
vertex_color <- unique(net_data[, c("cGenes", "Group")])
vertex_color <- setNames(vertex_color$Group, vertex_color$cGenes)

# Add atributes to network object
set.vertex.attribute(net, "Group", vertex_color[network.vertex.names(net)])
unique(dfnet$Group)  

df_color_dicio <- data.frame("Group" = c("AD","earlyPD", "MCI","MS","PD",
                                         "rrMS","spMS", "BPs"), 
                             "color" = c("#4a2e66", "#787878", "#8a508f", "#ff8531", 
                                         "grey30","#ffd331","#ffa600", "grey70"))

# Define names on data.frame
color_palette <- setNames(df_color_dicio$color, df_color_dicio$Group)  # O valor da cor

# Get Groups corresponding to nodes
df_color <- data.frame("Group" = get.vertex.attribute(net, "Group"))
#list.vertex.attributes(net)

#Replace NAs
df_color %>% 
  mutate(Group = ifelse(is.na(Group), "BPs", Group)) -> df_color
set.vertex.attribute(net, "Group", df_color$Group)
print(df_color)

#--------  Add System attribute to network object (SHAPES)
# Make list of attributes
vertex_shape <- unique(dfnet[, c("cGenes", "direction")])
vertex_shape <- setNames(vertex_shape$direction, vertex_shape$cGenes)

# Add attributes to network object
set.vertex.attribute(net, "direction", vertex_shape[network.vertex.names(net)])

# Define shapes 
unique(dfnet$direction) 
# Create a data.frame with names defining the shapes
df_shape_dicio <- data.frame("direction" = c("DOWN","UP", "Other"), 
                             "shape" = c(25, 24, 21))

# Define names on data.frame
shape_palette <- setNames(df_shape_dicio$shape, df_shape_dicio$direction)  # O valor da cor

# Get System correspondence to nodes
df_shape <- data.frame("direction" = get.vertex.attribute(net, "direction"))
print(df_shape)
# Replace NA for "Other"
df_shape %>% 
  mutate(direction = ifelse(is.na(direction), "Other", direction)) -> df_shape
set.vertex.attribute(net, "direction", df_shape$direction)

# arrange labels
issoai <- network.vertex.names(net)
issoai[(terms_n+1):length(issoai)] = gsub(issoai[(terms_n+1):length(issoai)], pattern = " .*", replacement = "")
issoai[1:terms_n] = ""
labsize = c(rep(2.5, terms_n), rep(1.5, length(issoai)-terms_n))

# arrange size
tam = tab_cluster_N %>% 
  dplyr::group_by(Term) %>% mutate(count = sum(n()))
tam2 = tam[,c(1,13)] %>% distinct()
issoai2 <- data.frame(network.vertex.names(net))
colnames(issoai2) = "Term"
issoai3 = left_join(issoai2, tam2)
issoai3 = issoai3$count
issoai3[which(is.na(issoai3))] = 1

set.seed(58)
ggnet2(net, 
       label = network.vertex.names(net)[1:(terms_n)],
       label.size = 3, 
       size = issoai3,
       repel = TRUE,
       alpha = 0.8,
       edge.size = 0.1, 
       shape = "direction",
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
  ggrepel::geom_text_repel(aes(label = issoai), 
                           size = labsize, max.overlaps = 30)

# save file (!!!)
ggsave(filename = "../../figures/network/Network_syngo.svg", device = "svg",  
       width = 8, height = 5, units = "in") 
