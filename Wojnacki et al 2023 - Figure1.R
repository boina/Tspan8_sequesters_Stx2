##These are the required packages to run the run the bioinformatic analysis in Wojnacki et al. Nature Communications. 2023.
require("readr")
require("data.table")
require("tidyr")
require("stringr")
require("dplyr")
require("ggplot2")
require("svglite")
require("extrafont")
require("forcats")

##Plots Custom theme
CustomTheme_b <- theme(plot.caption = element_text(color = "red", face = "bold"),
                     axis.text.x = element_text(vjust = -0.5, face = "bold"),
                     legend.text = element_text(size = 5),
                     legend.justification = "top",
                     legend.spacing.y = unit(1, 'mm'),
                     legend.key.size = unit(1, "mm"),
                     legend.position = "none",
		     text = element_text(family = "Cantarell", size = 9))


##Single cell RNAseq data has been previously published by Deprez et al. 2020.
##Meta data and expression matrix data can be downloaded from https://www.genomique.eu/cellbrowser/HCA/

##For this code to work as is, plase both data and meta data files in the R working directory.

##Load metadata
AirwaysMetaData <- read_tsv("meta.tsv")

##Load data
AirwaysData <- fread("exprMatrix.tsv")

##Genes Names
AirwaysGenes <- fread("exprMatrix.tsv", select = 1)

###Split the data set into smaller chunks to fit the analyis in less than 16Gb RAM.
##Calculate Split Factors.
splitFactor  <-  factor(round(8 * runif(length(AirwaysGenes$gene))))
splitedGenes <- split(AirwaysGenes$gene, splitFactor)

##Split the data set into subset datasets and save them.
print(paste("Splitting data set"))
for(n in 1:length(splitedGenes)){
    TempAirwaysData <- AirwaysData %>% filter(gene %in% splitedGenes[[n]])
    print(paste("Saving data set chunk ", n, " out of ", length(splitedGenes), sep = ""))
    saveRDS(TempAirwaysData, file = paste("AirwaysData_",n,"_of_",length(splitedGenes),".rds", sep = ""))
    rm(TempAirwaysData)
    print("Done")
    gc()
}
gc()
print("Done saving dataset chunks")
#####+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++###########


##Load each data chunk and calculate the change in gene's expression and proportion of cells expressing them comparing basal cells with secretory cells.

##Select only basal and secretory cells
AirwaysMetaDataSmall <- AirwaysMetaData[,c(1, 5, 10)] %>%
    mutate(CellType = str_remove_all(CellType, pattern = " N")) %>%
     filter(CellType %in% c("Basal", "Secretory"))

CellSelection <- AirwaysMetaDataSmall %>% pull(Id)

ChangeInProportionOfCells  <- list()
ChangeInGenesExpression <- list()
for(n in 1:length(splitedGenes)){
    DataChunkFile <- paste("AirwaysData_",n,"_of_",length(splitedGenes),".rds", sep = "")
    print(paste("Loading File:", DataChunkFile))
    TempDataChunk <- readRDS(file = DataChunkFile)
    print(paste("Done!! Now processing."))
    TempDataChunk <- TempDataChunk %>%
        dplyr::select(all_of(c("gene", CellSelection))) %>%
        pivot_longer(cols = -gene, names_to = "Id")
    print("Inner Join")
    TempDataChunk <- TempDataChunk %>%
        inner_join(., AirwaysMetaDataSmall, by = "Id")
    print("done inner join")
    TempDataChunk <- TempDataChunk %>%
        select(gene, value, Position, CellType) %>%
        mutate(CellType = str_remove_all(CellType, pattern = " N"))
    print("Calculating proportion of cells expressing a gene.")
    ChangeInProportionOfCells[[n]] <- TempDataChunk %>%
        group_by(Position, CellType, gene) %>%
        summarise(totalCells = n(),
                  positiveCells = sum(value > 0),
                  Proportion = positiveCells/totalCells) %>%
        ungroup() %>%
        select(Position, gene, CellType, Proportion) %>%
        pivot_wider(names_from = CellType,
                    names_prefix = "Prop_",
                    values_from = Proportion) %>%
        mutate(ChangeInProportion = (Prop_Secretory + 1) - (Prop_Basal + 1))
    print("Calculating change in gene expression.")
    ChangeInGenesExpression[[n]] <- TempDataChunk %>%
        filter(value > 0) %>%
        group_by(Position, CellType, gene) %>%
        summarise(meanValue = mean(value)) %>%
        ungroup() %>%
        pivot_wider(names_from = CellType,
                    names_prefix = "Expr_",
                    values_from = meanValue) %>%
        mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>%
        mutate(ChangeInExpression = (Expr_Secretory + 1) / (Expr_Basal + 1))
    remove(TempDataChunk)
    if_else((n != length(splitedGenes)),
    print("Done. Now processing next data chunk"),
    print("Done processing all data chunks!! Now creating a data frame from list of results."))
    gc()
}
ChangeInProportionOfCellsDF <- bind_rows(ChangeInProportionOfCells)
ChangeInGenesExpressionDF <- bind_rows(ChangeInGenesExpression)

##Create a combined data frame with changes in proportion and expression
ExpressionAndProportionChanges <- inner_join(ChangeInProportionOfCellsDF,
                                             ChangeInGenesExpressionDF,
                                             by = c("Position", "gene")) %>%
    mutate(Position = factor(Position, levels = selectedPositions <- c("Nasal","Proximal","Intermediate","Distal"))) %>%
    mutate(CombinedVariable = ChangeInProportion * (log(ChangeInExpression + 1))) %>%
        mutate(GOI = if_else(gene %in% (c("TSPAN8","MUC5AC","MUC5B")), TRUE,FALSE))
#####+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++###########


####Dotplot of change in gene expression and proportion of cells expressing each gene###
##Figure 1a

##Facets labels
LungRegionLabs <- paste(c("Nasal", "Proximal", "Intermediate", "Distal"), "airways")
names(LungRegionLabs) <- c("Nasal", "Proximal", "Intermediate", "Distal")

##Dotplot
ggplot(ExpressionAndProportionChanges %>% filter(Position == "Distal"),
       aes(x = Prop_Secretory, y = ChangeInExpression, size = GOI)) +
    geom_point(alpha=1/4) +
    scale_size_manual(values=c(0.6, 2), guide = "none") +
    geom_text(aes(label = ifelse(gene == "TSPAN8"|gene == "MUC5AC"|gene == "MUC5B",
                                 as.character(gene),'')),
              hjust = 0.5, vjust = -1) +
    facet_wrap(~ Position, ncol = 1, labeller = labeller(Position = LungRegionLabs)) +
    geom_hline(yintercept = 1, size = 0.4, linetype = 2) +
    ylim(0.5, 2) +
    labs(y = "Gene expression fold change",
         x = "Proportion of Cells Expressing Each Gene") +
    theme_bw() +
    CustomTheme_b

ggsave(filename = "DotPlot_Distal.svg", width = 90, height = 60, units = "mm")


###############################
###TSPAN8 Expression Boxplot###
###############################
##Figure 1b

##Define Genes of Interest
GOI <- c('MUC5AC', 'MUC5B', 'TSPAN8', 'ATG5', 'GAPDH')

##Filer expression values from Full data frame
AirwaysDataGOI <- filter(AirwaysData, gene %in% GOI)
gc()

##Merge with metadata to identify cell types
AirwaysDataGOI <- left_join(AirwaysDataGOI %>% pivot_longer(cols = -gene, names_to = "Id"), AirwaysMetaData, by = "Id") %>%
    dplyr::select(gene, Id, value, Position, CellType)

##Boxplot and dot plots
BoxPlotDataGOI <- AirwaysDataGOI %>%
    filter(CellType %in% c("Basal", "Secretory")) %>%
    filter(gene %in% c("TSPAN8", "MUC5AC", "MUC5B")) %>%
    mutate(gene = factor(gene, levels = c("MUC5AC","MUC5B","TSPAN8")))

BoxPlotDataGOI %>%
    ggplot(aes(x = CellType, y = value)) +
    geom_jitter(alpha = 0.05, size = 0.01, width = 0.35) +
    geom_boxplot(data = BoxPlotDataGOI %>% filter(value > 0),
                aes(x = CellType, y = value), alpha = 0)+
    facet_wrap(~ gene, nrow = 1) +
    labs(y = "Gene expression") +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    theme(axis.title.x = element_text(margin = margin(t = 10))) +
    theme_bw()
###+++++++++++++++++++++++++++++++++#####


###########################
####Expression Proportion##
###########################
##Figure 1c

##Calculation of the proportion of cells expressing each GOI
TSPAN8DataAbsExprProp <- AirwaysDataGOI %>%
    filter(CellType %in% c("Basal","Suprabasal","Multiciliated","Secretory")) %>%
    dplyr::select(gene, Id, value, Position, CellType) %>%
    pivot_wider(names_from = gene, values_from = value) %>%    
    dplyr::select(Id, Position, CellType, MUC5B, MUC5AC, ATG5, TSPAN8, GAPDH) %>%
    pivot_longer(cols = -c("Id", "Position", "CellType"), names_to = "Proteins") %>%
    group_by(CellType, Proteins) %>%
    summarise(prop = sum(value > 0) / n()) %>%
    mutate(Proteins = reorder(Proteins, prop)) %>%
    ungroup() %>%
    mutate(CellType = factor(CellType, levels = c("Basal","Suprabasal","Multiciliated","Secretory")))

##Expression proportion plot
TSPAN8DataAbsExprProp %>%
    ggplot(aes(x = fct_relevel(Proteins, c("MUC5AC", "MUC5B", "TSPAN8", "ATG5", "GAPDH")),
               y = CellType, size = prop)) +
    geom_point(color = "grey40") +
    guides(size = guide_legend(reverse = TRUE)) +
    labs(x = NULL, y = NULL, size = "Proportion") +
    theme_bw()

ggsave(filename = "ExpressionProp.svg", width = 95, height = 30, units = "mm")


########################
####Co_ExpressionPlot###
########################
##Figure 1c

##Calculation of the porportion of cells expressing co-expressing TSPAN8 and MUC5AC or MUC5B
TSPAN8DataCoExpre <- AirwaysDataGOI %>%
    filter(CellType %in% c("Basal","Suprabasal","Multiciliated","Secretory")) %>%
    dplyr::select(gene, Id, value, Position, CellType) %>%
    pivot_wider(names_from = gene, values_from = value) %>%    
    dplyr::select(Id, Position, CellType, MUC5B, MUC5AC, TSPAN8) %>%
    pivot_longer(cols = -c("Id", "Position", "CellType", "MUC5B", "MUC5AC"), names_to = "Proteins") %>%
    group_by(CellType, Proteins) %>%
    summarise(MUC5AC = sum(MUC5AC > 0 & value > 0) / n(),
              MUC5B  = sum(MUC5B  > 0 & value > 0) / n()) %>%
    mutate(Proteins = reorder(Proteins, MUC5AC)) %>%
    ungroup() %>%
    mutate(CellType = factor(CellType, levels = c("Basal","Suprabasal","Multiciliated","Secretory"))) %>%
    pivot_longer(cols = c("MUC5AC", "MUC5B"), names_to =  "Mucin", values_to = "percentage")

##Co-expressiong proportion plot
TSPAN8DataCoExpre %>% dplyr::filter(Proteins == "TSPAN8") %>%
    ggplot(aes(x = Mucin, y = CellType, size = percentage)) +
    geom_point(color = "grey40") +
    coord_flip()+
    guides(size = guide_legend(reverse = TRUE)) +
    labs(x = NULL, y = NULL, size = "Proportion") +
    theme_bw()
       
ggsave(filename = "CoExpressionProp.svg", width = 100, height = 30, units = "mm")
