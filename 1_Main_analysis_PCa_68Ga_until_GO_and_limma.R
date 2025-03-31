##### 
#Load and install required packages

library(devtools)
library(diann)
library(readxl)
library(stringr)
library(tidyverse)
library(mixOmics)
library(limma)
library(scales)
library(ggplot2)
library(tibble)
library(parameters)
library(xtable)
library(kableExtra)
library(data.table)
library(FactoMineR)
library(aod)
library(randomcoloR)
library(ComplexHeatmap)
library(ggrepel)
library(RColorBrewer)
library(DIAgui)
library(ggpmisc)
library(ggpubr)
library(corrplot)


#####
# Functions

# Median normalization
Median_polish = function(c){
  return((c-median(c, na.rm = T)))
}


# Histogram of p-values
histogram_toptable <- function(toptable, type, .x){
  
  np005 <- toptable %>% 
    filter(P.Value <= 0.05) %>% 
    pull(ID) %>% 
    length()
  
  np001 <- toptable %>% 
    filter(P.Value <= 0.01) %>% 
    pull(ID) %>% 
    length()
  
  histo_pvals <- ggplot(toptable) + 
    geom_histogram(aes(x = P.Value), 
                   binwidth = 0.005,
                   fill = "#2AB7CA", 
                   color = "#e9ecef", 
                   alpha=0.9) +
    geom_histogram(data = toptable %>% filter(P.Value <= 0.05),
                   aes(x = P.Value),
                   binwidth = 0.005,
                   fill = "#FE8585", 
                   color = "#e9ecef", 
                   alpha=0.9) +
    geom_histogram(data = toptable %>% filter(P.Value <= 0.01),
                   aes(x = P.Value),
                   binwidth = 0.005,
                   fill = "#FE4A49", 
                   color = "#e9ecef", 
                   alpha=0.9) + 
    geom_vline(xintercept = 0.05,
               color = "blue",
               size = 1, linetype = "dashed") +
    geom_vline(xintercept = 0.01,
               color = "red",
               size = 1, linetype = "dashed") +
    ggtitle(paste0("Distribution of p-values: ",
                   " proteins <= 0.05 = ",np005," ;proteins <= 0.01 = ", np001),
            subtitle = paste0("Fitting = ", type,
                              " ; Degrees of Freedom = ", .x)) + 
    theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.1, size = 10),
          panel.background = element_blank(),
          panel.grid.major = element_line(color = "grey"),
          panel.border = element_rect(colour = "black", fill=NA, size=1.5),
          axis.title=element_text(size=12,face="bold"))  
  
  return(histo_pvals)
}

# count NAs for missingness analysis
countNA = function(r){
  return(length(which(is.na(r))))
}


# count number of identified and quantified proteins
extractNoOfHits = function(c){
  return(length(which(!is.na(c))))
}

#Preprocessing of proteomic data from DIA-NN search output
##### 

## Load DIANN output report into R
df <- diann_load("PCa_68Ga_cohort_report.tsv")

# Extract expression values of interest and perform global normalization
protein.groups <- diann_matrix(df,id.header = "Protein.Group",quantity.header = "PG.MaxLFQ",
                               proteotypic.only = T,
                               q = 0.01,
                               protein.q = 0.01,
                               pg.q = 0.01, 
                               gg.q = 1)

# rename columns and remove .raw.dia and perform log2 transformation for normal-like distribution of intensities 
colnames(protein.groups) = gsub(".raw.dia*", "", colnames(protein.groups))
colnames(protein.groups)[c(4:80)] = str_sub(colnames(protein.groups[,c(4:80)]), -6,-1)
protein.groups[,c(4:80)] = log2(protein.groups[,c(4:80)])


## Perform Median-normalization
protein.groups_median = as.data.frame(protein.groups[,c(3:80)])
row.names(protein.groups_median) = protein.groups_median$Genes
protein.groups_median = protein.groups_median[,-1]
protein.groups_median = apply(protein.groups_median, 2, Median_polish)

## Create protein expression table with median normalized values
protein.groups_df = as.data.frame(protein.groups_median) %>% 
                                  tibble::rownames_to_column("ID")

#####

# Overview of clinical data
##### 

## Load clinical data for each subareal and reformat
Clin_data_subareal <- read_excel("68Ga_Subareal_Gleason_ISUP_SUVmax_ADCmin_QSZHGE_SUVindex.xlsx")

Clin_data_subareal = as.data.frame(Clin_data_subareal)
Clin_data_subareal[Clin_data_subareal ==0] = NA
row.names(Clin_data_subareal) = Clin_data_subareal$Sample

# Transform ADC and SUV columns in the Clinical dataframe to numeric category
Clin_data_subareal$ADC_value_min_subareal = as.numeric(Clin_data_subareal$ADC_value_min_subareal)
Clin_data_subareal$ADC_value_mean_subareal = as.numeric(Clin_data_subareal$ADC_value_mean_subareal)
Clin_data_subareal$ADC_vol_subareal = as.numeric(Clin_data_subareal$ADC_vol_subareal)
Clin_data_subareal$SUV_mean_subareal = as.numeric(Clin_data_subareal$SUV_mean_subareal)
Clin_data_subareal$SUV_max_subareal = as.numeric(Clin_data_subareal$SUV_max_subareal)
Clin_data_subareal$`SUV_Vol (cm°3)_subareal` = as.numeric(Clin_data_subareal$`SUV_Vol (cm°3)_subareal`)
Clin_data_subareal$`SUVindex (Tmax/Pmean)` = as.numeric(Clin_data_subareal$`SUVindex (Tmax/Pmean)`)
Clin_data_subareal$QSZHGE = as.numeric(Clin_data_subareal$QSZHGE)

# Select relevant columns from clinical data and rename columns
Clin_data_subareal_rel = Clin_data_subareal[which(!grepl("1", Clin_data_subareal$NMAT)),c(6:10,11:15)]
colnames(Clin_data_subareal_rel) = c("Gleason", "ISUP", "ADC min", "ADC mean", "ADC vol","SUV mean", "SUV max", "SUV vol", "SUVi_T_P", "QSZHGE")

## Adjust Range of Radiomic values

Clin_data_subareal_rel[,c(3,4,6,7)] = Clin_data_subareal_rel[,c(3,4,6,7)]/1000

## Dividing in PET derived values and MRI derived values

## PET
Clin_data_subareal_rel_PET = Clin_data_subareal_rel[,c(2,6:10)]

# Generate histograms of PET derived values for Figure 2A

png("SUV_mean.png", width = 900, height = 600, pointsize = 10, res = 150)
hist(Clin_data_subareal_rel_PET$`SUV mean`, main = "SUV mean", ylim = c(0,30))
dev.off()

png("SUV_max.png", width = 900, height = 600, pointsize = 10, res = 150)
hist(Clin_data_subareal_rel_PET$`SUV max`, breaks = 15, main = "SUV max", ylim = c(0,20))
dev.off()

png("SUVi_T_P.png", width = 900, height = 600, pointsize = 10, res = 150)
hist(Clin_data_subareal_rel_PET$SUVi_T_P, main = "SUVi T_P", ylim = c(0,30), xlim = c(0,50))
dev.off()

png("QSZHGE.png", width = 900, height = 600, pointsize = 10, res = 150)
hist(Clin_data_subareal_rel_PET$QSZHGE, breaks = 15, main = "QSZHGE", ylim = c(0,35), xlim = c(0,1600))
dev.off()

png("ISUP.png", width = 900, height = 600, pointsize = 10, res = 150)
hist(Clin_data_subareal_rel$ISUP, breaks = 15, main = "ISUP", ylim = c(0,25), xlim = c(1,4))
dev.off()

# Perform correlation analysis for PET derived values for Figure 2B
cor_PET = Clin_data_subareal_rel_PET[which(apply(Clin_data_subareal_rel_PET, 1, countNA) == 0),c(1,3:6)]

svg(filename="Correlation_PET_parameters.svg",width=14, height=10, pointsize=14)
corrplot.mixed(corr = cor(cor_PET), lower = "number", upper = "ellipse")
dev.off()


## MRI
Clin_data_subareal_rel_MRI = Clin_data_subareal_rel[,c(2:5)]

# Generate histograms of PET derived values for Figure 2A

png("ADC_mean.png", width = 900, height = 600, pointsize = 10, res = 150)
hist(Clin_data_subareal_rel_MRI$`ADC mean`, main = "ADC mean", ylim = c(0,15), xlim = c(0,2000))
dev.off()

png("ADC_min.png", width = 900, height = 600, pointsize = 10, res = 150)
hist(Clin_data_subareal_rel_MRI$`ADC min`, main = "ADC min", ylim = c(0,20), xlim = c(0,1000))
dev.off()

summary(Clin_data_subareal_rel_MRI)

# Perform correlation analysis for PET derived values for Figure 2C

cor_MRI = Clin_data_subareal_rel_MRI[which(apply(Clin_data_subareal_rel_MRI, 1, countNA) == 0),]

svg(filename="Correlation_MRI_parameters.svg",width=14, height=10, pointsize=14)
corrplot.mixed(corr = cor(cor_MRI), lower = "number", upper = "ellipse")
dev.off()

#####

# Overview of proteomic data 
#####
## Number of proteins
# sort data frame of clinical data for increasing patient ID
Clin_data_subareal_sorted = Clin_data_subareal[order(Clin_data_subareal$`Patient ID`),]
sample_order <- Clin_data_subareal_sorted$Sample

# Generate Patient ID vector for color coding
Patien_vec = Clin_data_subareal$`Patient ID`

# Generate 20 distinguishable colors
palette <- distinctColorPalette(20)
patient_colors <- setNames(palette, 1:20)

# Count number of identified and quantified proteins
hits = apply(protein.groups_df, 2, extractNoOfHits)

# Barchart of protein IDs for Figure 3A
svg(filename="Barplot_68Ga_PCa_cohort.svg",width=18, height=10, pointsize=10)
xx = barplot(hits[-1],las = 3, cex.axis = 1.2, ylim = c(0,3000), col = patient_colors[factor(Patien_vec, levels = 1:20)])
dev.off()


## Filter for proteins with only 5 NAs
protein.groups_df_stable = protein.groups_df[which(apply(protein.groups_df, 1, countNA)<=5),]

# Unsupervised PCA analysis for Figure 3B
pdf("Unsupervised_PCA_68Ga_PCa_cohort_with_patient_label.pdf", height = 10, width = 12)
plotIndiv(mixOmics::pca(t(protein.groups_df_stable[,-1])),  legend = T, title = "Principal component analysis based on proteome data of 68Ga PCa cohort", size.legend = 16, ind.names = T,  cex = 5, pch = c(18), size.xlabel = rel(2), size.ylabel = rel(2), size.axis = rel(1.3), group = factor(Patien_vec, levels = 1:20), col.per.group = patient_colors)
dev.off()


## Heatmap of selected POIs in prostate cancer development and progression for Figure 3C

## filter Clin_data_subareal to generate two new dataframes which contain "1" or "NA" in the NMAT column
Clin_data_subareal_Back = Clin_data_subareal[which(Clin_data_subareal$NMAT == 1),]
Clin_data_subareal_Tum = Clin_data_subareal[which(!grepl("1", Clin_data_subareal$NMAT)),]

# sort the Clin_data_subareal_TUM and Back data based on increasing SUVindex Tmax/Pmean With NA rows at the top
Clin_data_subareal_Back = Clin_data_subareal_Back[order(Clin_data_subareal_Back$`SUVindex (Tmax/Pmean)`),]
Clin_data_subareal_Tum = Clin_data_subareal_Tum[order(Clin_data_subareal_Tum$ISUP),]

# Get order for rownames
order_rownames = paste(c(rownames(Clin_data_subareal_Back[which(is.na(Clin_data_subareal_Back$`SUVindex (Tmax/Pmean)`)),]),rownames(Clin_data_subareal_Back[which(!is.na(Clin_data_subareal_Back$`SUVindex (Tmax/Pmean)`)),]),rownames(Clin_data_subareal_Tum)))

# Combine the two dataframes
Clin_data_subareal_comb = rbind(Clin_data_subareal_Back,Clin_data_subareal_Tum)

# sort rows according to order_rownames
Clin_data_subareal_comb = Clin_data_subareal_comb[order_rownames,]

# filter and sort the original dataframe to only include samples from the order_rownames
protein.groups_df_filtered = protein.groups_df[,c("ID", order_rownames)]

# Create an ISUP annotation vector for the heatmap
Back_vec = c(rep("back",12))
ISUP_vec = c(rep("ISUP1",8),rep("ISUP2",19),rep("ISUP3",21),rep("ISUP4",17))
column_annotation_ISUP = c(Back_vec,ISUP_vec)

# Define list of POIs
PCa_POI_vec = c("FOLH1","KLK3", "GSTP1")

# Create a column annotation for the heatmap
ha_POI = columnAnnotation(ISUP = Clin_data_subareal_comb$ISUP,
                            col = list(ISUP = c("1" = "yellow", "2" = "orange", "3" = "red", "4" = "darkred")),
                            na_col = "forestgreen")

# Generate sub dataframes and reformat for heatmap
PCa_POI_sub = protein.groups_df_filtered
PCa_POI_sub =filter(PCa_POI_sub, PCa_POI_sub$ID %in% PCa_POI_vec)

row.names(PCa_POI_sub) = PCa_POI_sub$ID
PCa_POI_sub = PCa_POI_sub[,-1]

PCa_POI_sub_mat = as.matrix(PCa_POI_sub)

# Generate complex heatmap with POIs and ISUP categories

Heatmap_PCa_POI = Heatmap(PCa_POI_sub_mat, name = "PCa_POI POIs for PC in ISUP high vs low PCa samples", column_title = "Samples", row_title = "Biomarker",
                        column_split = column_annotation_ISUP,
                        column_labels = Patien_vec,
                        cluster_rows = F, cluster_columns =F,
                        row_order = c("FOLH1","KLK3", "GSTP1"),
                        bottom_annotation = ha_POI,
                        heatmap_legend_param = list(title = "Median-polish intensities", title_position = "lefttop-rot"))

# Generate svg and png graphics
svg(filename="Heatmap_PCa_POI_ISUP_ALDH1A1_labeled.svg",width=20, height=10, pointsize=12)
Heatmap_PCa_POI
dev.off()

png("Heatmap_PCa_POI_ISUP_labeled.png", width = 1600, height = 800)
Heatmap_PCa_POI
dev.off()

## Multi-correlation analysis for PET derived values for Figure 4A and B
#####
# Prepare data for two matrix correlation analysis
protein.groups_w_GenName = protein.groups_df
row.names(protein.groups_w_GenName) = ifelse(is.na(protein.groups_w_GenName$ID), protein.groups_w_GenName$ID,protein.groups_w_GenName$ID)

# Filter for proteins that have a gene name
protein.groups_w_GenName = protein.groups_w_GenName[which((!grepl("SWISS", protein.groups_w_GenName$ID)&!grepl("=", protein.groups_w_GenName$ID)&!grepl(":", protein.groups_w_GenName$ID))),]

# Filter protein table to only contain tumor samples, and reformat the table
protein.groups_w_GenName = protein.groups_w_GenName %>%
  dplyr::select(ID, rownames(Clin_data_subareal_rel))

protein.groups_w_GenName_mat = as.matrix(protein.groups_w_GenName[,-1])
t_protein = t(protein.groups_w_GenName_mat)

# Perform multi-correlation analysis of PET derived values and protein profiles
spls_res_subareal_PET <- spls(X = t_protein,
                              Y = Clin_data_subareal_rel_PET[c(3,5,6)], 
                              ncomp = 2,
                              mode = "regression",
                              keepX = c(5,5), keepY = c(3,3)) 

# Generate correlation circle plot
svg(filename="Circle_corr_subareal_PET_GenName.svg",width=14, height=10, pointsize=12)
plotVar(spls_res_subareal_PET, title = "Correlation Circle Plot PET-CT and protein intensities", cex = c(5,5), font = c(1,2), cutoff = 0.5, style = "ggplot2", legend.title = "Matrix", legend = c("Gene","PET/CT"))
dev.off()

svg(filename="Circle_corr_subareal_PET_GenName_wo_label.svg",width=14, height=10, pointsize=12)
plotVar(spls_res_subareal_PET, title = "Correlation Circle Plot PET-CT and protein intensities", var.names = F, pch = c(20,20), cex = c(5,5), font = c(1,2), cutoff = 0.5, style = "ggplot2", legend.title = "Matrix", legend = c("Gene","PET/CT"))
dev.off()

# Perform multi-correlation analysis of MRI derived values and protein profiles
spls_res_subareal_MRI <- spls(X = t_protein,
                              Y = Clin_data_subareal_rel_MRI[c(2,3)], 
                              ncomp = 2,
                              mode = "regression",
                              keepX = c(5,5), keepY = c(2,2)) 

# Generate correlation circle plot
svg(filename="Circle_corr_subareal_MRI_GenName.svg",width=14, height=10, pointsize=12)
plotVar(spls_res_subareal_MRI, title = "Correlation Circle Plot MRI and protein intensities", cex = c(5,5), font = c(1,2), cutoff = 0.5, style = "ggplot2", legend.title = "Matrix", legend = c("Gene","MRI"))
dev.off()

svg(filename="Circle_corr_subareal_MRI_GenName_wo_labels.svg",width=14, height=10, pointsize=12)
plotVar(spls_res_subareal_MRI, title = "Correlation Circle Plot MRI and protein intensities", var.names = F, pch = c(20,20),cex = c(5,5), font = c(1,2), cutoff = 0.5, style = "ggplot2", legend.title = "Matrix", legend = c("Gene","MRI"))
dev.off()

#####

# Linear regression analysis using either SUV_max,SUVi T_P, QSZHGE, ADC_min, and ADC mean values for Figure 4C
#####

## SUV max

# filter for patient with SUV values
Clin_data_subareal_rel_SUV_woNA = Clin_data_subareal_rel_PET[,which(grepl("SUV", colnames(Clin_data_subareal_rel_PET)))] %>% na.omit()

# Focus on SUV_max extracted from PSMA PET/CT
SUV_max = Clin_data_subareal_rel_SUV_woNA[,which(grepl(" max", colnames(Clin_data_subareal_rel_SUV_woNA)))]
design_SUV <- model.matrix(~SUV_max)

row.names(design_SUV) <- row.names(Clin_data_subareal_rel_SUV_woNA)
coefs <- 2

protein.groups_info_SUV = protein.groups_w_GenName[,-1]
protein.groups_info_SUV$protein_ID = row.names(protein.groups_info_SUV)

protein.groups_info_SUV <- protein.groups_info_SUV %>%
  dplyr::select(protein_ID, rownames(Clin_data_subareal_rel_SUV_woNA))

protein.groups_info_SUV = protein.groups_info_SUV[,-c(which(grepl("ID", colnames(protein.groups_info_SUV))))] %>%
  as.matrix()

fit_SUV <- lmFit(protein.groups_info_SUV, 
                 design = design_SUV, 
                 method = "ls")
fit_SUV <- eBayes(fit_SUV)

toptable_SUV <- topTable(fit_SUV, coef = coefs, number = Inf) %>% 
  tibble::rownames_to_column("ID")

fited_vals_SUV <- fitted.MArrayLM(fit_SUV) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("ID") %>%
  pivot_longer(cols = 2:ncol(.),
               names_to = c("Sample"),
               values_to = "Fitted_Abundance")

limma_out_SUV <- list(limma_fit = fit_SUV,
                      toptable = toptable_SUV,
                      fitted_values = fited_vals_SUV)


# Generate histogram of p-values
png(filename="Histogram_lin_reg_SUV_max_subareal.png", width = 1200, height = 800, pointsize = 8, res = 150)
histogram_toptable(limma_out_SUV$toptable,
                   type = "Linear",
                   .x = NULL)
dev.off()


## SUV index Tumor / Prostate
# Focus on SUV index extracted from PSMA PET/CT
SUV_T_P = Clin_data_subareal_rel_SUV_woNA[,which(grepl("T_P", colnames(Clin_data_subareal_rel_SUV_woNA)))]
design_SUV_T_P <- model.matrix(~SUV_T_P)

row.names(design_SUV_T_P) <- row.names(Clin_data_subareal_rel_SUV_woNA)
coefs <- 2

protein.groups_info_SUV = protein.groups_w_GenName[,-1]
protein.groups_info_SUV$protein_ID = row.names(protein.groups_info_SUV)

protein.groups_info_SUV <- protein.groups_info_SUV %>%
  dplyr::select(protein_ID, rownames(Clin_data_subareal_rel_SUV_woNA))

protein.groups_info_SUV = protein.groups_info_SUV[,-c(which(grepl("ID", colnames(protein.groups_info_SUV))))] %>%
  as.matrix()

fit_SUV_T_P <- lmFit(protein.groups_info_SUV, 
                 design = design_SUV_T_P, 
                 method = "ls")
fit_SUV_T_P <- eBayes(fit_SUV_T_P)

toptable_SUV_T_P <- topTable(fit_SUV_T_P, coef = coefs, number = Inf) %>% 
  tibble::rownames_to_column("ID")

fited_vals_SUV_T_P <- fitted.MArrayLM(fit_SUV_T_P) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("ID") %>%
  pivot_longer(cols = 2:ncol(.),
               names_to = c("Sample"),
               values_to = "Fitted_Abundance")

limma_out_SUV_T_P <- list(limma_fit = fit_SUV_T_P,
                      toptable = toptable_SUV_T_P,
                      fitted_values = fited_vals_SUV_T_P)

# Generate histogram of p-values
png("Lin_reg_hist_SUV_T_P.png", width = 1200, height = 800, pointsize = 8, res = 150)
histogram_toptable(limma_out_SUV_T_P$toptable,
                   type = "Linear",
                   .x = NULL)
dev.off()

## QSZHGE
# Focus on QSZHGE extracted from PSMA PET/CT
Clin_data_subareal_rel_QSZHGE_woNA = Clin_data_subareal_rel_PET[,c(1,which(grepl("QSZHGE", colnames(Clin_data_subareal_rel_PET))))] %>% na.omit()

QSZHGE = Clin_data_subareal_rel_QSZHGE_woNA[,2]
design_QSZHGE <- model.matrix(~QSZHGE)

row.names(design_QSZHGE) <- row.names(Clin_data_subareal_rel_QSZHGE_woNA)
coefs <- 2

protein.groups_info_QSZHGE = protein.groups_w_GenName[,-1]
protein.groups_info_QSZHGE$protein_ID = row.names(protein.groups_info_QSZHGE)

protein.groups_info_QSZHGE <- protein.groups_info_QSZHGE %>%
  dplyr::select(protein_ID, rownames(Clin_data_subareal_rel_QSZHGE_woNA))

protein.groups_info_QSZHGE = protein.groups_info_QSZHGE[,-c(which(grepl("ID", colnames(protein.groups_info_QSZHGE))))] %>%
  as.matrix()

fit_QSZHGE <- lmFit(protein.groups_info_QSZHGE, 
                    design = design_QSZHGE, 
                    method = "ls")
fit_QSZHGE <- eBayes(fit_QSZHGE)

toptable_QSZHGE <- topTable(fit_QSZHGE, coef = coefs, number = Inf) %>% 
  tibble::rownames_to_column("ID")

fited_vals_QSZHGE <- fitted.MArrayLM(fit_QSZHGE) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("ID") %>%
  pivot_longer(cols = 2:ncol(.),
               names_to = c("Sample"),
               values_to = "Fitted_Abundance")

limma_out_QSZHGE <- list(limma_fit = fit_QSZHGE,
                         toptable = toptable_QSZHGE,
                         fitted_values = fited_vals_QSZHGE)

# Generate histogram of p-values
png("Lin_reg_hist_QSZHGE.png", width = 1200, height = 800, pointsize = 8, res = 150)
histogram_toptable(limma_out_QSZHGE$toptable,
                   type = "Linear",
                   .x = NULL)
dev.off()

# ADC min
## filter for patient with ADC values
Clin_data_subareal_rel_ADC_woNA = Clin_data_subareal_rel_MRI[,which(grepl("ADC", colnames(Clin_data_subareal_rel_MRI)))] %>% na.omit()

# Focus on ADC min extracted from mpMRI
ADC_min = Clin_data_subareal_rel_ADC_woNA[,1]
design_ADC <- model.matrix(~ADC_min)

row.names(design_ADC) <- row.names(Clin_data_subareal_rel_ADC_woNA)
coefs <- 2

protein.groups_info_ADC = protein.groups_w_GenName[,-1]
protein.groups_info_ADC$protein_ID = row.names(protein.groups_info_ADC)

protein.groups_info_ADC <- protein.groups_info_ADC %>%
  dplyr::select(protein_ID, rownames(Clin_data_subareal_rel_ADC_woNA))

protein.groups_info_ADC = protein.groups_info_ADC[,-c(which(grepl("ID", colnames(protein.groups_info_ADC))))] %>%
  as.matrix()

fit_ADC <- lmFit(protein.groups_info_ADC, 
                 design = design_ADC, 
                 method = "ls")
fit_ADC <- eBayes(fit_ADC)

toptable_ADC <- topTable(fit_ADC, coef = coefs, number = Inf) %>% 
  tibble::rownames_to_column("ID")

fited_vals_ADC <- fitted.MArrayLM(fit_ADC) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("ID") %>%
  pivot_longer(cols = 2:ncol(.),
               names_to = c("Sample"),
               values_to = "Fitted_Abundance")

limma_out_ADC <- list(limma_fit = fit_ADC,
                      toptable = toptable_ADC,
                      fitted_values = fited_vals_ADC)

# Generate histogram of p-values
png("Lin_reg_hist_ADC_min.png", width = 1200, height = 800, pointsize = 8, res = 150)
histogram_toptable(limma_out_ADC$toptable,
                   type = "Linear",
                   .x = NULL)
dev.off()


# ADC mean
# Focus on ADC mean extracted from mpMRI
ADC_mean = Clin_data_subareal_rel_ADC_woNA[,2]
design_ADC_mean <- model.matrix(~ADC_mean)

row.names(design_ADC_mean) <- row.names(Clin_data_subareal_rel_ADC_woNA)
coefs <- 2

protein.groups_info_ADC = protein.groups_w_GenName[,-1]
protein.groups_info_ADC$protein_ID = row.names(protein.groups_info_ADC)

protein.groups_info_ADC <- protein.groups_info_ADC %>%
  dplyr::select(protein_ID, rownames(Clin_data_subareal_rel_ADC_woNA))

protein.groups_info_ADC = protein.groups_info_ADC[,-c(which(grepl("ID", colnames(protein.groups_info_ADC))))] %>%
  as.matrix()

fit_ADC_mean <- lmFit(protein.groups_info_ADC, 
                 design = design_ADC_mean, 
                 method = "ls")
fit_ADC_mean <- eBayes(fit_ADC_mean)

toptable_ADC_mean <- topTable(fit_ADC_mean, coef = coefs, number = Inf) %>% 
  tibble::rownames_to_column("ID")

fited_vals_ADC_mean <- fitted.MArrayLM(fit_ADC_mean) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("ID") %>%
  pivot_longer(cols = 2:ncol(.),
               names_to = c("Sample"),
               values_to = "Fitted_Abundance")

limma_out_ADC_mean <- list(limma_fit = fit_ADC_mean,
                      toptable = toptable_ADC_mean,
                      fitted_values = fited_vals_ADC_mean)

# Generate histogram of p-values
png("Lin_reg_hist_ADC_mean.png", width = 1200, height = 800, pointsize = 8, res = 150)
histogram_toptable(limma_out_ADC_mean$toptable,
                   type = "Linear",
                   .x = NULL)
dev.off()

#####

# Linear regression analysis using ISUP values for Figure 5
#####
# ISUP
## filter for patient with QSZHGE values
Clin_data_subareal_rel_ISUP_woNA = Clin_data_subareal_rel[,c(1,which(grepl("ISUP", colnames(Clin_data_subareal_rel))))] %>% na.omit()

ISUP = Clin_data_subareal_rel_ISUP_woNA[,2]
design_ISUP <- model.matrix(~ISUP)

row.names(design_ISUP) <- row.names(Clin_data_subareal_rel_ISUP_woNA)
coefs <- 2

protein.groups_info_ISUP = protein.groups_w_GenName[,-1]
protein.groups_info_ISUP$protein_ID = row.names(protein.groups_info_ISUP)

protein.groups_info_ISUP <- protein.groups_info_ISUP %>%
  dplyr::select(protein_ID, rownames(Clin_data_subareal_rel_ISUP_woNA))

protein.groups_info_ISUP = protein.groups_info_ISUP[,-c(which(grepl("ID", colnames(protein.groups_info_ISUP))))] %>%
  as.matrix()

fit_ISUP <- lmFit(protein.groups_info_ISUP, 
                  design = design_ISUP, 
                  method = "ls")
fit_ISUP <- eBayes(fit_ISUP)

toptable_ISUP <- topTable(fit_ISUP, coef = coefs, number = Inf) %>% 
  tibble::rownames_to_column("ID")

fited_vals_ISUP <- fitted.MArrayLM(fit_ISUP) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("ID") %>%
  pivot_longer(cols = 2:ncol(.),
               names_to = c("Sample"),
               values_to = "Fitted_Abundance")

limma_out_ISUP <- list(limma_fit = fit_ISUP,
                       toptable = toptable_ISUP,
                       fitted_values = fited_vals_ISUP)

# Generate histogram of p-values
png("Lin_reg_hist_ISUP.png", width = 1200, height = 800, pointsize = 8, res = 150)
histogram_toptable(limma_out_ISUP$toptable,
                   type = "Linear",
                   .x = NULL)
dev.off()

## Generate Heatmap of core matrisomal proteins
# Get list of core matrisome proteins
Core_Matri_list = read_excel("Core_Matri.xlsx")
Core_Matri_vec = Core_Matri_list$gene_id

Core_Matri_sub = protein.groups_df_filtered
Core_Matri_sub =filter(Core_Matri_sub, Core_Matri_sub$ID %in% Core_Matri_vec)

row.names(Core_Matri_sub) = Core_Matri_sub$ID
Core_Matri_sub = Core_Matri_sub[,-1]

Core_Matri_sub_mat = as.matrix(Core_Matri_sub)

# extract the associated ECM function / group for the row annotation
Core_Matri_prots_in_data <- Core_Matri_list %>%
  filter(gene_id %in% protein.groups_df_filtered$ID) %>% # keep identified proteins
  filter(!duplicated(gene_id))
ordered_Core_Matri_annotation <- Core_Matri_prots_in_data %>%
  arrange(factor(gene_id, levels = row.names(Core_Matri_sub_mat)))
row_annotation_Core_Matri <- ordered_Core_Matri_annotation$category

# generate column annotation for ISUP and Lesion
ha = columnAnnotation(ISUP = Clin_data_subareal_comb$ISUP,
                      Lesion = Clin_data_subareal_comb$Lesion,
                      col = list(Lesion = c("1" = "coral2", "2" = "goldenrod2"),ISUP = c("1" = "yellow", "2" = "orange", "3" = "red", "4" = "darkred")),
                      na_col = "forestgreen")



# Generate complex Heatmap

Heatmap_Core_Matri = Heatmap(Core_Matri_sub_mat, name = "Core matrisome proteins in ISUP high vs low PCa samples", column_title = "Samples",
                             row_split = row_annotation_Core_Matri,
                             column_split = column_annotation_ISUP,
                             row_names_side = c("left"),
                             column_labels = Patien_vec,
                             cluster_rows = F, cluster_columns =F,
                             bottom_annotation = ha,
                             row_title_gp = gpar(fontsize = 15),
                             row_title_rot = 0,
                             row_title_side = "right",
                             heatmap_legend_param = list(title = "Median-polish intensities", title_position = "lefttop-rot"))

# Generate svg graphics
svg(filename="Heatmap_Core_Matri_ISUP_labeled.svg",width=16, height=10, pointsize=12)
Heatmap_Core_Matri
dev.off()

#####

# Perform Gene ontology enrichment analysis of correlated proteins with imaging-derived values

#####
# Load Uniprot annotation file to merge Uniprot ID and GO annotation to the linear regression results
Uniprot_ID_Gene_Protname_GO <- read_excel("Uniprot_ID_Gene_Protname_GO.xlsx")
colnames(Uniprot_ID_Gene_Protname_GO) = c("Protein_ID","ID", "Protein names", "GO BP", "GO CC", "GO MF")
Uniprot_ID_Gene_Protname_GO$ID = gsub(" .*", "", Uniprot_ID_Gene_Protname_GO$ID)

# Merge Uniprot ID and GO annotation to the linear regression results and export for GO analysis in separate file

# SUV max
toptable_SUV_w_GO = merge(toptable_SUV, Uniprot_ID_Gene_Protname_GO, by = "ID")
write.table(toptable_SUV_w_GO, "Proteinlist_after_lin_regression_with_SUV_max_subareal.txt", sep = "\t", quote = F, row.names = F, dec = ".")

# SUV index Tumor / Prostate
toptable_SUV_T_P_w_GO = merge(toptable_SUV_T_P, Uniprot_ID_Gene_Protname_GO, by = "ID")
write.table(toptable_SUV_T_P_w_GO, "Proteinlist_after_lin_regression_with_SUV_T_P_subareal.txt", sep = "\t", quote = F, row.names = F, dec = ".")

# QSZHGE
toptable_QSZHGE_w_GO = merge(toptable_QSZHGE, Uniprot_ID_Gene_Protname_GO, by = "ID")
write.table(toptable_QSZHGE_w_GO, "Proteinlist_after_lin_regression_with_QSZHGE_subareal.txt", sep = "\t", quote = F, row.names = F, dec = ".")

# ADC min
toptable_ADC_w_GO = merge(toptable_ADC, Uniprot_ID_Gene_Protname_GO, by = "ID")
write.table(toptable_ADC_w_GO, "Proteinlist_after_lin_regression_with_ADC_min_subareal.txt", sep = "\t", quote = F, row.names = F, dec = ".")

# ADC mean
toptable_ADC_mean_w_GO = merge(toptable_ADC_mean, Uniprot_ID_Gene_Protname_GO, by = "ID")
write.table(toptable_ADC_mean_w_GO, "Proteinlist_after_lin_regression_with_ADC_mean_subareal.txt", sep = "\t", quote = F, row.names = F, dec = ".")

# ISUP
toptable_ISUP_w_GO = merge(toptable_ISUP, Uniprot_ID_Gene_Protname_GO, by = "ID")
write.table(toptable_ISUP_w_GO, "Proteinlist_after_lin_regression_with_ISUP_subareal.txt", sep = "\t", quote = F, row.names = F, dec = ".")

#####

# Generate protein expression table for limma analysis stratified by each imaging-derived value
#####
### SUV index T_P

# select 33% of samples with lowest and highest SUV index T_P

SUVi_T_P_low = Clin_data_subareal_rel_SUV_woNA %>% 
  slice_min(order_by = SUVi_T_P , prop = 0.33)

SUVi_T_P_high = Clin_data_subareal_rel_SUV_woNA %>% 
  slice_max(order_by = SUVi_T_P , prop = 0.33)

protein.groups_info_SUVi_T_P <- as.data.frame(protein.groups_info_SUV)
protein.groups_info_SUVi_T_P$ID = row.names(protein.groups_info_SUVi_T_P)

protein.groups_info_SUVi_T_P_limma <-  protein.groups_info_SUVi_T_P %>%
  dplyr::select(ID, rownames(SUVi_T_P_low),rownames(SUVi_T_P_high))

write.table(protein.groups_info_SUVi_T_P_limma, "Proteinlist_limma_input_SUVi_T_P.txt", sep = "\t", quote = F, row.names = F, dec = ".")


### SUV max

#select 33% of samples with lowest and highest SUV max

SUV_max_low = Clin_data_subareal_rel_SUV_woNA %>% 
  slice_min(order_by = `SUV max`, prop = 0.33)

SUV_max_high = Clin_data_subareal_rel_SUV_woNA %>% 
  slice_max(order_by = `SUV max`, prop = 0.33)

protein.groups_info_SUV_max <- as.data.frame(protein.groups_info_SUV)
protein.groups_info_SUV_max$ID = row.names(protein.groups_info_SUV_max)

protein.groups_info_SUV_max_limma <-  protein.groups_info_SUV_max %>%
  dplyr::select(ID, rownames(SUV_max_low),rownames(SUV_max_high))

write.table(protein.groups_info_SUV_max_limma, "Proteinlist_limma_input_SUV_max.txt", sep = "\t", quote = F, row.names = F, dec = ".")

### ADC min


#select 33% of samples with lowest and highest ADC min

ADC_min_low = Clin_data_subareal_rel_ADC_woNA %>% 
  slice_min(order_by = `ADC min`, prop = 0.33)

ADC_min_high = Clin_data_subareal_rel_ADC_woNA %>% 
  slice_max(order_by = `ADC min`, prop = 0.33)

protein.groups_info_ADC_min <- as.data.frame(protein.groups_info_ADC)
protein.groups_info_ADC_min$ID = row.names(protein.groups_info_ADC_min)

protein.groups_info_ADC_min_limma <-  protein.groups_info_ADC_min %>%
  dplyr::select(ID, rownames(ADC_min_low),rownames(ADC_min_high))

write.table(protein.groups_info_ADC_min_limma, "Proteinlist_limma_input_ADC_min.txt", sep = "\t", quote = F, row.names = F, dec = ".")

### ADC mean


#select 33% of samples with lowest and highest ADC mean

ADC_mean_low = Clin_data_subareal_rel_ADC_woNA %>% 
  slice_min(order_by = `ADC mean`, prop = 0.33)

ADC_mean_high = Clin_data_subareal_rel_ADC_woNA %>% 
  slice_max(order_by = `ADC mean`, prop = 0.33)

protein.groups_info_ADC_mean <- as.data.frame(protein.groups_info_ADC)
protein.groups_info_ADC_mean$ID = row.names(protein.groups_info_ADC_mean)

protein.groups_info_ADC_mean_limma <-  protein.groups_info_ADC_mean %>%
  dplyr::select(ID, rownames(ADC_mean_low),rownames(ADC_mean_high))

write.table(protein.groups_info_ADC_mean_limma, "Proteinlist_limma_input_ADC_mean.txt", sep = "\t", quote = F, row.names = F, dec = ".")


### QSZHGE

#select 33% of samples with lowest and highest QSZHGE

QSZHGE_low = Clin_data_subareal_rel_QSZHGE_woNA %>% 
  slice_min(order_by = QSZHGE, prop = 0.33)

QSZHGE_high = Clin_data_subareal_rel_QSZHGE_woNA %>% 
  slice_max(order_by = QSZHGE, prop = 0.33)

protein.groups_info_QSZHGE <- as.data.frame(protein.groups_info_QSZHGE)
protein.groups_info_QSZHGE$ID = row.names(protein.groups_info_QSZHGE)

protein.groups_info_QSZHGE_limma <-  protein.groups_info_QSZHGE %>%
  dplyr::select(ID, rownames(QSZHGE_low),rownames(QSZHGE_high))

write.table(protein.groups_info_QSZHGE_limma, "Proteinlist_limma_input_QSZHGE.txt", sep = "\t", quote = F, row.names = F, dec = ".")

#####