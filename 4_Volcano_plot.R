### Script to generate high resolution volcano blot

# The following script is a generalized version which can be customized for each comparison illustrated in Figure 7


# Load limma results 
limma_results <- read.delim("limma_results.tsv")

# Load required packages

library(ggplot2)
library(plotly)
library(ggrepel)
library(readxl)
library(EnhancedVolcano)

### Retrieve and load table from uniprot to map protein names, gene names and GO annotations for BP, CC and MF

Uniprot_ID_Gene_Protname_GO <- read_excel("Uniprot_ID_Gene_Protname_GO.xlsx")
colnames(Uniprot_ID_Gene_Protname_GO) = c("Uniprot_ID","ID", "Protein names", "GO BP", "GO CC", "GO MF")

Uniprot_ID_Gene_Protname_GO$ID = gsub(" .*", "", Uniprot_ID_Gene_Protname_GO$ID)


### Merge files

limma_results = merge(limma_results, Uniprot_ID_Gene_Protname_GO[,c(2,1,4:6)], by = "ID")

sig_up = subset(limma_results, limma_results$mean.log2..QSZHGE_high...QSZHGE_low..> log2(1.5) & limma_results$adj.P.Val < 0.05)

sig_down = subset(limma_results, limma_results$mean.log2..QSZHGE_high...QSZHGE_low.. < -log2(1.5) & limma_results$adj.P.Val < 0.05)

sig_comb = rbind(sig_up, sig_down)


### round adj. p-values below 0.00005 to 0.00005 to enable improved visualization
limma_results_round_adjpval = limma_results
limma_results_round_adjpval$round_adj_pval = ifelse(limma_results_round_adjpval$adj.P.Val< 0.00005, 0.00005, limma_results_round_adjpval$adj.P.Val)

### Generate vectors of POIs

## based on protein families
# e.g. S100_vec = limma_results$ID[grepl("S100", limma_results$ID)]

## based on GO terms and expressions
# e.g. Prolif_vec = sig_comb$ID[grepl("positive regulation of cell population proliferation", sig_comb$`GO BP`)]
# or Prolif_vec = sig_comb$ID[grepl("GO:0008284", sig_comb$`GO BP`)]


## based on other a priory defined baskets
# PC_prog_list = read_excel("PC_prog_targets.xlsx")
# PC_prog_vec = PC_prog_list$gene_id

# Exemplary plot without color customization

S100_vec = limma_results$ID[grepl("S100", limma_results$ID)]

Test_S100 = EnhancedVolcano(limma_results,
                lab = limma_results$ID,
                x = 'mean.log2..QSZHGE_high...QSZHGE_low..',
                y = 'adj.P.Val',
                title = 'QSZHGE high vs low',
                subtitle = "Proteins from the S100 family",
                selectLab = S100_vec,
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~ 'adj.P.Val'),                                                              
                ylim = c(0,max(-log10(limma_results$adj.P.Val))+1),
                pCutoff = 5*10e-3,
                FCcutoff = log2(1.5),
                pointSize = 4.0,
                labSize = 6.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                colAlpha = 4/5,
                legendPosition = 'right',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black')

png("S100_proteins.png", width = 1600, height = 800)
Test_S100
dev.off()


#### Use POI list for each comparison and customize colors for the volcano plot
# Example for QSZHGE
QSZHGE_POI_list <- read.csv("QSZHGE_POI_list.txt", header=FALSE, )
QSZHGE_POI_vec = as.character(QSZHGE_POI_list[,1])

Wound_vec = limma_results_round_adjpval$ID[grepl("wound healing", limma_results_round_adjpval$`GO BP`)]
Muscle_vec = limma_results_round_adjpval$ID[grepl("muscle", limma_results_round_adjpval$`GO BP`)]
Adhesion_vec = limma_results_round_adjpval$ID[grepl("cell adhesion", limma_results_round_adjpval$`GO BP`)]
Fatty_vec = sig_comb$ID[grepl("fatty acid beta-oxidation", sig_comb$`GO BP`)]

ADC_list = read_excel("ADC_targets.xlsx")
ADC_vec = ADC_list$gene_id
PC_prog_list = read_excel("PC_prog_targets.xlsx")
PC_prog_vec = PC_prog_list$gene_id

# Refine vectors for intersecting POIs
QSZHGE_POI_Wound = intersect(QSZHGE_POI_vec, Wound_vec)
QSZHGE_POI_Muscle = intersect(QSZHGE_POI_vec, Muscle_vec)
QSZHGE_POI_Adhes = intersect(QSZHGE_POI_vec, Adhesion_vec)
QSZHGE_POI_ADC = intersect(QSZHGE_POI_vec, ADC_vec)
QSZHGE_POI_Prog = intersect(QSZHGE_POI_vec, PC_prog_vec)
QSZHGE_POI_Fatty = intersect(QSZHGE_POI_vec, Fatty_vec)

# Set general color scheme for 
limma_results_round_adjpval$BP.colour = 'lightgray'
limma_results_round_adjpval$BP.colour = ifelse(limma_results_round_adjpval$mean.log2..QSZHGE_high...QSZHGE_low..> log2(1.5)&limma_results_round_adjpval$adj.P.Val < 0.05, "darkgray", ifelse(limma_results_round_adjpval$mean.log2..QSZHGE_high...QSZHGE_low..< -log2(1.5)&limma_results_round_adjpval$adj.P.Val < 0.05, "darkgray", "lightgray"))

limma_results_round_adjpval$BP.colour[match(QSZHGE_POI_Muscle,limma_results_round_adjpval$ID)] = 'skyblue'
limma_results_round_adjpval$BP.colour[match(QSZHGE_POI_Adhes,limma_results_round_adjpval$ID)] = 'lawngreen'
limma_results_round_adjpval$BP.colour[match(QSZHGE_POI_Prog,limma_results_round_adjpval$ID)] = 'gold'
limma_results_round_adjpval$BP.colour[match(QSZHGE_POI_ADC,limma_results_round_adjpval$ID)] = 'blue'
limma_results_round_adjpval$BP.colour[match(QSZHGE_POI_Fatty,limma_results_round_adjpval$ID)] = 'darkred'

BP.colour = as.character(limma_results_round_adjpval$BP.colour)


names(BP.colour)[BP.colour == 'skyblue'] <- '6 Muscle tissue development'
names(BP.colour)[BP.colour == 'lawngreen'] <- '5 Cell adhesion'
names(BP.colour)[BP.colour == 'blue'] <- '2 ADC target'
names(BP.colour)[BP.colour == 'gold'] <- '3 PCa prognosis'
names(BP.colour)[BP.colour == 'darkred'] <- '7 Fatty acid metabolism'
names(BP.colour)[BP.colour == 'darkgray'] <- '8 Sig. dysregulated proteins'
names(BP.colour)[BP.colour == 'lightgray'] <- '9 Not dysregulated'

POI_subset = limma_results_round_adjpval[match(QSZHGE_POI_vec,limma_results_round_adjpval$ID),]


Final_QSZHGE_POI_color = EnhancedVolcano(limma_results_round_adjpval,
                                        lab = limma_results_round_adjpval$ID,
                                        x = 'mean.log2..QSZHGE_high...QSZHGE_low..',
                                        y = 'round_adj_pval',
                                        title = 'QSZHGE high vs low',
                                        subtitle = "Dysregulated Protreins of interest in QSZHGE POI comparison",
                                        selectLab = QSZHGE_POI_vec,
                                        xlab = bquote(~Log[2]~ 'fold change'),
                                        ylab = bquote(~-Log[10]~ 'adj.P.Val'),
                                        ylim = c(0,max(-log10(limma_results_round_adjpval$round_adj_pval))+0.25),
                                        xlim = c(min(limma_results_round_adjpval$mean.log2..QSZHGE_high...QSZHGE_low..)-0.25, max(limma_results_round_adjpval$mean.log2..QSZHGE_high...QSZHGE_low..)+0.25),
                                        pCutoff = 5*10e-3,
                                        FCcutoff = log2(1.5),
                                        pointSize = c(ifelse(BP.colour=="lightgray", 1,ifelse(BP.colour=="darkgray", 2,5))),
                                        labSize = 8.0,
                                        labCol = 'black',
                                        labFace = 'bold',
                                        boxedLabels = TRUE,
                                        colCustom = BP.colour,
                                        colAlpha = 4/5,
                                        legendPosition = 'right',
                                        legendLabSize = 20,
                                        legendIconSize = 6.0,
                                        drawConnectors = TRUE,
                                        widthConnectors = 1.0,
                                        colConnectors = 'black')


png("QSZHGE_POI.png", width = 1200, height = 1000, res=80)
Final_QSZHGE_POI_color
dev.off()