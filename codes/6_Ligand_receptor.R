# 6. Neurotransmission ---------------------------------------------------------

## 6.1 Load packages and functions ---------------------------------------------

source("./codes/my_packages.R")
source("./codes/my_functions.R")

# install packages from github:
#devtools::install_github("jinworks/CellChat")
#devtools::install_github("jokergoo/circlize")
library("CellChat")
library("circlize")

# set working directory
setwd("./data/GSEbatch/")
getwd()

## 6.2 Ligand-receptor circos --------------------------------------------------

# create directory to save images
dir.create("../../figures/circos/")

# get CellChat list of ligand-receptors
db = CellChat::CellChatDB.human
db2 = db$interaction
#write.table(db2, file = "./CellChatdb.tsv", sep = "\t", row.names = F)

# if unable to load CellChat, use the table:
db2 = read.table(file = "../CellChatdb.tsv", sep = "\t", header = T)

# check for consistency in is.neurotransmitter annotation for pathways
df3 = db2 %>% 
  dplyr::group_by(pathway_name, is_neurotransmitter) %>% 
  summarise(n=n()) %>% 
  pivot_wider(names_from = is_neurotransmitter, values_from = n) %>% 
  dplyr::rename(true = 'TRUE', false = 'FALSE') %>% 
  filter(!is.na(true)) %>% # filter pathways which have a TRUE for 'is_neurotransmitter'
  filter(is.na(false)) # filter pathways which don't have a FALSE for 'is_neurotransmitter'

# filter CellChat db for neurotransmitter pathways 
db3 = db2 %>% dplyr::filter(pathway_name %in% df3$pathway_name)


# table for alluvial plot
tab_allu = db3[0,] 
tab_allu = tab_allu %>% # make columns to rbind
  mutate(direction = NA) %>% 
  mutate(disease = NA) %>% 
  mutate(LR = NA)
tab_allu_pres = data.frame(pres = character()) %>% 
  mutate(direction = NA) %>% 
  mutate(disease = NA) %>% 
  mutate(LR = NA)

for (grupao in c("AD","MCI", "PD", "earlyPD", "MS", "spMS", "rrMS")) {
  fname = paste0("./DRAs_", grupao, "_GSEbatch.tsv")
  tab0 = read.table(file = fname, header = T, sep = "\t", na.strings = "")
  tab1 = tab0[,2:3]
  colnames(tab1) = c("UP", "DOWN")
  # tab2 = tab1 %>% 
  #   pivot_longer(cols =  c("UP", "DOWN"), 
  #                names_to = "direction", values_to = "Genes")
  tabL = db3[0,] %>% mutate(direction = NA)
  tabR = db3[0,] %>% mutate(direction = NA)
  presL = data.frame(pres = character()) %>% mutate(direction = NA)
  presR = data.frame(pres = character()) %>% mutate(direction = NA)
  
  for (i in c("UP", "DOWN")) {
    # filter rows that have the genes from the meta DEAs table
    my_genes = tab1[!is.na(tab1[,i]),i]#tab2$Genes
    # Ligands
    dbL = db3 %>%
      # separate ligand names 
      separate_rows(ligand.symbol, sep = ", ") 
    # filter which are in the list of my_genes
    dbL = dbL %>% dplyr::filter(ligand.symbol %in% my_genes) 
    # get list of which were in the list of my_genes
    presL1 = dbL %>% dplyr::distinct(ligand.symbol) %>% 
      #rename column to 'pres' (presence, which were in the list of my_genes)
      dplyr::rename(pres = ligand.symbol) 
    dbL = dbL %>% distinct(interaction_name, .keep_all = T)
    # filter database
    tabL1 = db3 %>% dplyr::filter(interaction_name %in% dbL$interaction_name) 
    
    
    if (nrow(tabL1) != 0) {
      tabL1$direction = i # add direction information
      presL1$direction = i
      tabL = rbind(tabL, tabL1) 
      presL = rbind(presL, presL1) }
    
    
    #Receptors
    dbR = db3 %>%
      # separate receptor names 
      separate_rows(receptor.symbol, sep = ", ") 
    # filter which are in the list of my_genes
    dbR = dbR %>% dplyr::filter(receptor.symbol %in% my_genes) 
    presR1 = dbR %>% dplyr::distinct(receptor.symbol) %>% 
      # rename column to 'pres' (presence, which were in the list of my_genes)
      dplyr::rename(pres = receptor.symbol) 
    dbR = dbR %>% distinct(interaction_name, .keep_all = T)
    # filter db table for interactions where receptors were in the list
    tabR1 = db3 %>% dplyr::filter(interaction_name %in% dbR$interaction_name) 
    
    if (nrow(tabR1) != 0) {
      tabR1$direction = i # add direction information
      presR1$direction = i # add direction information
      tabR = rbind(tabR, tabR1) 
      presR = rbind(presR, presR1)}
  }
  
  
  # if there are results (has rows), add identification cols and append to tab_allu 
  if (nrow(tabL) != 0) {
    tabL$disease = grupao
    tabL$LR = "LIGAND"
    tab_allu = rbind(tab_allu, tabL)
  }
  
  if (nrow(presL) != 0) {
    presL$disease = grupao
    presL$LR = "LIGAND"
    tab_allu_pres = rbind(tab_allu_pres, presL)
  }
  
  if (nrow(tabR) != 0) {
    tabR$disease = grupao
    tabR$LR = "RECEPTOR"
    tab_allu = rbind(tab_allu, tabR)
  }
  
  if (nrow(presR) != 0) {
    presR$disease = grupao
    presR$LR = "RECEPTOR"
    tab_allu_pres = rbind(tab_allu_pres, presR)
  }
}

# reorder columns
tab_allu2 = tab_allu %>%
  relocate(pathway_name, ligand.symbol, receptor.symbol, 
           disease, LR, direction, is_neurotransmitter,
           .before = 1)

#write.table(tab_allu2, file = "./RL_cellchatdb_filtered.txt", row.names = F, sep = "\t")

# add count for each line
tab_allu2$num_rows = 1

tab_allu3 = tab_allu2 %>% 
  separate_rows(ligand.symbol, sep = ", ") %>% 
  separate_rows(receptor.symbol, sep = ", ") %>% 
  distinct(pathway_name, ligand.symbol, receptor.symbol, disease, LR, num_rows)


# NOTE: Add your pathways and corresponding colors
# save file (!!!)
for (group in c("AD","MCI", "PD", "earlyPD", "MS", "spMS", "rrMS")) {
  tab_allu4 = tab_allu3 %>%  
  dplyr::filter(disease == group) %>% 
  dplyr::rename(ligand = ligand.symbol, receptor = receptor.symbol) %>% 
  group_by(ligand, receptor, pathway_name) %>% 
  summarize(n = n()) %>% 
  # color by pathway
  mutate(cores = if_else(pathway_name == "GABA-A", "#74a892", false = 
                 if_else(pathway_name == "GABA-B", "#008585", false = 
                 if_else(pathway_name == "Glutamate", "#e5c185", false =
                 if_else(pathway_name == "Glycine", "#c7522a", false = 
                 if_else(pathway_name == "TRH", "#db7d42", false =
                 if_else(pathway_name == "PMCH", "#c27ba0", false =
                 if_else(pathway_name == "CALC", "#cee583", false =
                 if_else(pathway_name == "5-HT", "#195382", false =
                 if_else(pathway_name == "Dopamine", "#057dcd", false =
                 if_else(pathway_name == "SerotoninDopamin", "#74c1ed", 
                         false = NA))))))))))) %>% 
  mutate(pathway_name = factor(pathway_name, 
                               levels = c("GABA-A","GABA-B","Glutamate",
                                          "Glycine","TRH","PMCH","CALC","5-HT",
                                          "Dopamine", "SerotoninDopamin")))
  
  # color of the external bar
  grid.col = tab_allu4 %>% pivot_longer(cols = c(ligand, receptor))
  grid.col2 = setNames(grid.col$cores, grid.col$value)
  
  # color of tracks (connections)
  col = tab_allu4$cores
  
  # sector order
  orderL = tab_allu4 %>% ungroup() %>% arrange(pathway_name) %>% distinct(ligand)
  orderR = tab_allu4 %>% ungroup() %>% arrange(pathway_name) %>% distinct(receptor)
  order1 = c(rev(orderL$ligand), (orderR$receptor))
  
  # highlight labels
  tab_allu_pres_cond = tab_allu_pres %>%  
    dplyr::filter(disease == cond) 
  
  
  nome = paste0("../../figures/circos/", group, ".pdf")
  #png(filename = nome, res = 300, units = "px", width = 2000, height = 2000)
  pdf(file = nome, width = 6.5, height = 6.5)
  
  # Restart circular layout parameters
  circos.clear()
  # circos.par(start.degree = 90) # to rotate plot
  chordDiagram(tab_allu4, 
               grid.col = grid.col2, # set colors
               preAllocateTracks = 1, 
               order = order1, # set sector order
               annotationTrackHeight = mm_h(2), # height of external bar
               annotationTrack =  "grid")
  # add labels:
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    
    ## highlight specific ligands/receptors with bold font:
    if (sector.name %in% tab_allu_pres_cond$pres) {
      circos.text(mean(xlim), ylim[1] + .1, # label position 
                  sector.name, 
                  col = ifelse(tab_allu_pres_cond[which(
                    tab_allu_pres_cond$pres == sector.name),"direction"] == "UP",
                    "#d43e3e", "#394789"), 
                  cex = 0.6, # font size
                  facing = "clockwise", 
                  font = 2, # bold
                  niceFacing = TRUE, adj = c(-0, 0.5))} 
    
    ## add the other labels in normal font:
    else{
      circos.text(mean(xlim), ylim[1] + .1, # label position 
                  sector.name, 
                  cex = 0.9, # font size
                  facing = "clockwise", 
                  niceFacing = TRUE, adj = c(-0, 0.5))}
  }, bg.border = NA)
  
  # add legend
  tab_allu5 = tab_allu4 %>% ungroup() %>% 
    distinct(pathway_name, .keep_all = T) %>% 
    dplyr::arrange(pathway_name)
  lgd_points = Legend(at = as.character(unique(tab_allu5$pathway_name)), 
                      type = "grid", #squares
                      legend_gp = gpar(fill = tab_allu5$cores),
                      title_position = "topleft", 
                      gap = unit(0.3, "cm"),  
                      title = "Pathway")
  
  draw(lgd_points, x = unit(1, "cm"), y = unit(1, "cm"), just = c("left", "bottom"))
  
  dev.off()
}
