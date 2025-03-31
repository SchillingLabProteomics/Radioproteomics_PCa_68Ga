#####
# Load required packages for GO analysis

library(lattice)
library(ggplot2)
library(tidyverse)
library(org.Hs.eg.db)
library(clusterProfiler)

# Perform GO analysis for proteins that are significantly correlated with imaging-derived features
#####
# SUV T_P 

# P.Value / adj.P.Val
proteins.regulated.up.1 <- read.delim("Proteinlist_after_lin_regression_with_SUV_T_P_subareal.txt") %>%
  filter(logFC > 0, adj.P.Val <= 0.05) %>%
  pull(Protein_ID)

proteins.regulated.down.1 <- read.delim("Proteinlist_after_lin_regression_with_SUV_T_P_subareal.txt") %>%
  filter(logFC < 0, adj.P.Val <= 0.05) %>%
  pull(Protein_ID)

proteins.df.up.1 <- bitr(proteins.regulated.up.1,
                         fromType = "UNIPROT",
                         toType = "ENTREZID",
                         OrgDb = org.Hs.eg.db)

proteins.df.down.1 <- bitr(proteins.regulated.down.1,
                           fromType = "UNIPROT",
                           toType = "ENTREZID",
                           OrgDb = org.Hs.eg.db)

proteins.eGO.up.1 <- enrichGO(proteins.df.up.1$ENTREZID,
                              OrgDb=org.Hs.eg.db,
                              ont="BP",
                              #universe = background.proteome$ENTREZID,
                              pAdjustMethod="BH",
                              pvalueCutoff=0.05,
                              qvalueCutoff=1)

proteins.eGO.down.1 <- enrichGO(proteins.df.down.1$ENTREZID,
                                OrgDb=org.Hs.eg.db,
                                ont="BP",
                                #universe = background.proteome$ENTREZID, use this
                                pAdjustMethod="BH",
                                pvalueCutoff=0.05,
                                qvalueCutoff=1)

## Create dataframe with all POIs from different lin. regression analysis
df_SUVi_T_P_up = tibble(entre = proteins.df.up.1$ENTREZID) %>%
  mutate(condition = "SUV_TP",
         direction = "upreg")

df_SUVi_T_P_down = tibble(entre = proteins.df.down.1$ENTREZID) %>%
  mutate(condition = "SUV_TP",
         direction = "downreg")

df_SUVi_T_P = rbind(df_SUVi_T_P_up, df_SUVi_T_P_down)

##
# SUV max 

# P.Value / adj.P.Val
proteins.regulated.up.2 <- read.delim("Proteinlist_after_lin_regression_with_SUV_max_subareal.txt") %>%
  filter(logFC > 0, adj.P.Val <= 0.05) %>%
  pull(Protein_ID)

proteins.regulated.down.2 <- read.delim("Proteinlist_after_lin_regression_with_SUV_max_subareal.txt") %>%
  filter(logFC < 0, adj.P.Val <= 0.05) %>%
  pull(Protein_ID)

proteins.df.up.2 <- bitr(proteins.regulated.up.2,
                         fromType = "UNIPROT",
                         toType = "ENTREZID",
                         OrgDb = org.Hs.eg.db)

proteins.df.down.2 <- bitr(proteins.regulated.down.2,
                           fromType = "UNIPROT",
                           toType = "ENTREZID",
                           OrgDb = org.Hs.eg.db)

proteins.eGO.up.2 <- enrichGO(proteins.df.up.2$ENTREZID,
                              OrgDb=org.Hs.eg.db,
                              ont="BP",
                              #universe = background.proteome$ENTREZID,
                              pAdjustMethod="BH",
                              pvalueCutoff=0.05,
                              qvalueCutoff=1)

proteins.eGO.down.2 <- enrichGO(proteins.df.down.2$ENTREZID,
                                OrgDb=org.Hs.eg.db,
                                ont="BP",
                                #universe = background.proteome$ENTREZID, use this
                                pAdjustMethod="BH",
                                pvalueCutoff=0.05,
                                qvalueCutoff=1)


## Create dataframe with all POIs from different lin. regression analysis
df_SUV_max_up = tibble(entre = proteins.df.up.2$ENTREZID) %>%
  mutate(condition = "SUV_max",
         direction = "upreg")

df_SUV_max_down = tibble(entre = proteins.df.down.2$ENTREZID) %>%
  mutate(condition = "SUV_max",
         direction = "downreg")

df_SUV_max = rbind(df_SUV_max_up, df_SUV_max_down)

##
# ADC_min

# P.Value / adj.P.Val
proteins.regulated.up.3 <- read.delim("Proteinlist_after_lin_regression_with_ADC_min_subareal.txt") %>%
  filter(logFC > 0, P.Value <= 0.05) %>%
  pull(Protein_ID)

proteins.regulated.down.3 <- read.delim("Proteinlist_after_lin_regression_with_ADC_min_subareal.txt") %>%
  filter(logFC < 0, P.Value <= 0.05) %>%
  pull(Protein_ID)

proteins.df.up.3 <- bitr(proteins.regulated.up.3,
                         fromType = "UNIPROT",
                         toType = "ENTREZID",
                         OrgDb = org.Hs.eg.db)

proteins.df.down.3 <- bitr(proteins.regulated.down.3,
                           fromType = "UNIPROT",
                           toType = "ENTREZID",
                           OrgDb = org.Hs.eg.db)

proteins.eGO.up.3 <- enrichGO(proteins.df.up.3$ENTREZID,
                              OrgDb=org.Hs.eg.db,
                              ont="BP",
                              #universe = background.proteome$ENTREZID,
                              pAdjustMethod="BH",
                              pvalueCutoff=0.05,
                              qvalueCutoff=1)

proteins.eGO.down.3 <- enrichGO(proteins.df.down.3$ENTREZID,
                                OrgDb=org.Hs.eg.db,
                                ont="BP",
                                #universe = background.proteome$ENTREZID, use this
                                pAdjustMethod="BH",
                                pvalueCutoff=0.05,
                                qvalueCutoff=1)


## Create dataframe with all POIs from different lin. regression analysis
df_ADC_up = tibble(entre = proteins.df.up.3$ENTREZID) %>%
  mutate(condition = "ADC_min",
         direction = "upreg")

df_ADC_down = tibble(entre = proteins.df.down.3$ENTREZID) %>%
  mutate(condition = "ADC_min",
         direction = "downreg")

df_ADC = rbind(df_ADC_up, df_ADC_down)

##
# ADC_mean

# P.Value / adj.P.Val
proteins.regulated.up.4 <- read.delim("Proteinlist_after_lin_regression_with_ADC_mean_subareal.txt") %>%
  filter(logFC > 0, P.Value <= 0.05) %>%
  pull(Protein_ID)

proteins.regulated.down.4 <- read.delim("Proteinlist_after_lin_regression_with_ADC_mean_subareal.txt") %>%
  filter(logFC < 0, P.Value <= 0.05) %>%
  pull(Protein_ID)

proteins.df.up.4 <- bitr(proteins.regulated.up.4,
                         fromType = "UNIPROT",
                         toType = "ENTREZID",
                         OrgDb = org.Hs.eg.db)

proteins.df.down.4 <- bitr(proteins.regulated.down.4,
                           fromType = "UNIPROT",
                           toType = "ENTREZID",
                           OrgDb = org.Hs.eg.db)

proteins.eGO.up.4 <- enrichGO(proteins.df.up.4$ENTREZID,
                              OrgDb=org.Hs.eg.db,
                              ont="BP",
                              #universe = background.proteome$ENTREZID,
                              pAdjustMethod="BH",
                              pvalueCutoff=0.05,
                              qvalueCutoff=1)

proteins.eGO.down.4 <- enrichGO(proteins.df.down.4$ENTREZID,
                                OrgDb=org.Hs.eg.db,
                                ont="BP",
                                #universe = background.proteome$ENTREZID, use this
                                pAdjustMethod="BH",
                                pvalueCutoff=0.05,
                                qvalueCutoff=1)


## Create dataframe with all POIs from different lin. regression analysis
df_ADC_mean_up = tibble(entre = proteins.df.up.4$ENTREZID) %>%
  mutate(condition = "ADC_mean",
         direction = "upreg")

df_ADC_mean_down = tibble(entre = proteins.df.down.4$ENTREZID) %>%
  mutate(condition = "ADC_mean",
         direction = "downreg")

df_ADC_mean = rbind(df_ADC_mean_up, df_ADC_mean_down)

##
# QSZHGE 

# P.Value / adj.P.Val
proteins.regulated.up.5 <- read.delim("Proteinlist_after_lin_regression_with_QSZHGE_subareal.txt") %>%
  filter(logFC > 0, adj.P.Val <= 0.05) %>%
  pull(Protein_ID)

proteins.regulated.down.5 <- read.delim("Proteinlist_after_lin_regression_with_QSZHGE_subareal.txt") %>%
  filter(logFC < 0, adj.P.Val <= 0.05) %>%
  pull(Protein_ID)

proteins.df.up.5 <- bitr(proteins.regulated.up.5,
                         fromType = "UNIPROT",
                         toType = "ENTREZID",
                         OrgDb = org.Hs.eg.db)

proteins.df.down.5 <- bitr(proteins.regulated.down.5,
                           fromType = "UNIPROT",
                           toType = "ENTREZID",
                           OrgDb = org.Hs.eg.db)

proteins.eGO.up.5 <- enrichGO(proteins.df.up.5$ENTREZID,
                              OrgDb=org.Hs.eg.db,
                              ont="BP",
                              #universe = background.proteome$ENTREZID,
                              pAdjustMethod="BH",
                              pvalueCutoff=0.05,
                              qvalueCutoff=1)

proteins.eGO.down.5 <- enrichGO(proteins.df.down.5$ENTREZID,
                                OrgDb=org.Hs.eg.db,
                                ont="BP",
                                #universe = background.proteome$ENTREZID, use this
                                pAdjustMethod="BH",
                                pvalueCutoff=0.05,
                                qvalueCutoff=1)


## Create dataframe with all POIs from different lin. regression analysis
df_QSZHGE_up = tibble(entre = proteins.df.up.5$ENTREZID) %>%
  mutate(condition = "QSZHGE",
         direction = "upreg")

df_QSZHGE_down = tibble(entre = proteins.df.down.5$ENTREZID) %>%
  mutate(condition = "QSZHGE",
         direction = "downreg")

df_QSZHGE = rbind(df_QSZHGE_up, df_QSZHGE_down)

##
# ISUP 

# P.Value / adj.P.Val
proteins.regulated.up.6 <- read.delim("Proteinlist_after_lin_regression_with_ISUP_subareal.txt") %>%
  filter(logFC > 0, adj.P.Val <= 0.05) %>%
  pull(Protein_ID)

proteins.regulated.down.6 <- read.delim("Proteinlist_after_lin_regression_with_ISUP_subareal.txt") %>%
  filter(logFC < 0, adj.P.Val <= 0.05) %>%
  pull(Protein_ID)

proteins.df.up.6 <- bitr(proteins.regulated.up.6,
                         fromType = "UNIPROT",
                         toType = "ENTREZID",
                         OrgDb = org.Hs.eg.db)

proteins.df.down.6 <- bitr(proteins.regulated.down.6,
                           fromType = "UNIPROT",
                           toType = "ENTREZID",
                           OrgDb = org.Hs.eg.db)

proteins.eGO.up.6 <- enrichGO(proteins.df.up.6$ENTREZID,
                              OrgDb=org.Hs.eg.db,
                              ont="BP",
                              #universe = background.proteome$ENTREZID,
                              pAdjustMethod="BH",
                              pvalueCutoff=0.05,
                              qvalueCutoff=1)

proteins.eGO.down.6 <- enrichGO(proteins.df.down.6$ENTREZID,
                                OrgDb=org.Hs.eg.db,
                                ont="BP",
                                #universe = background.proteome$ENTREZID, use this
                                pAdjustMethod="BH",
                                pvalueCutoff=0.05,
                                qvalueCutoff=1)


#### Create dataframe with all POIs from different lin. regression analysis
df_ISUP_up = tibble(entre = proteins.df.up.6$ENTREZID) %>%
  mutate(condition = "ISUP",
         direction = "upreg")

df_ISUP_down = tibble(entre = proteins.df.down.6$ENTREZID) %>%
  mutate(condition = "ISUP",
         direction = "downreg")

df_ISUP = rbind(df_ISUP_up, df_ISUP_down)

#####

# Combine all dataframes
df_PET = rbind(df_SUVi_T_P, df_SUV_max, df_QSZHGE, df_ISUP)
df_MRI = rbind(df_ADC, df_ADC_mean)


df_all = rbind(df_PET, df_MRI)
df_all$condition <- factor(df_all$condition, levels = c("SUV_max", "SUV_TP", "QSZHGE", "ADC_min", "ADC_mean", "ISUP"))


all_eGO_BP <- compareCluster(entre~condition+direction, data=df_all, OrgDb = org.Hs.eg.db, fun="enrichGO", ont="BP")

# order the x axis in the dotplot showing SUV max, SUVi T_P, QSZHGE, ADC min, ADC mean and ISUP
all_eGO_BP@compareClusterResult$condition <- factor(all_eGO_BP@compareClusterResult$condition, levels = c("SUV_max", "SUV_TP", "QSZHGE", "ADC_min", "ADC_mean", "ISUP"))

simplified_eGO_BP <- simplify(all_eGO_BP, cutoff = 0.6, by = "p.adjust", measure = "Wang", select_fun = min)


# Generate svg and png graphics
svg(filename="Lin_regression_PET_ISUP_MRI_GO_BP.svg",width=14, height=14, pointsize=14)
dotplot(simplified_eGO_BP, x="condition", size = "count",
        showCategory = 3,
        font.size = 12) + facet_grid(~direction)
dev.off()

png(filename="Lin_regression_PET_ISUP_MRI_GO_BP.png", width = 1600, height = 2100, pointsize = 10, res = 150)
dotplot(simplified_eGO_BP, x="condition", size = "count",
        showCategory = 3,
        label_format = 30,
        font.size = 12) + facet_grid(~direction)
dev.off()