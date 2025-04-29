# Radioproteomics_PCa_68Ga

## Overview
This repository contains the analysis pipeline for the manuscript: "Advancing Radioproteomics: Integrating PSMA PET/CT and mpMRI with Localized Proteomic Profiling in Prostate Cancer" from Fahrner et al. 

This study is based on radioproteomic data from 20 prostate cancer patients, integrating in vivo bioimaging (PSMA PET/CT and mpMRI) with ex vivo imaging and histopathology. Using a novel 3D-to-2D registration pipeline, we aligned imaging features with histological sections to identify molecular signatures linked to tumor heterogeneity, enabling improved diagnostic accuracy and risk stratification.

# How to use the code
The scripts are numbered to correspond to the main steps in the analysis pipeline described in the manuscript.

### 1. Preprocessing and main analysis
1_Main_analysis_PCa_68Ga_until_GO_and_limma.R
Loads DIA-NN output and performs normalization of proteomic values, unsupervised clustering and correlation analysis. Additionally, visualization of imaging-dervied value distribution as well as correlation analysis is done in this script.

### 2. Gene ontology analysis
2_GO_analysis.R
Uses the results from the linear regression analysis for each covariable with the protein profiles to perform GO enrichment analysis of the correlating proteins.

### 3. Limma analysis
3_limma.R
Uses the limma_input protein tables and performs differential statistical analysis of the protein profiles from the upper and lower tertiale of each covariable. 

### 4. Volcano plots
4_Volcano_plot.R
Loads the results from differential statistial analysis for generating volcano plots with customized POI highlighted and colored.

## Input data
All mass spectrometric raw data and analysis results are deposited in the MassIVE repository together with the relevant clinical annotation. Reviewer login credentials for MassIVE are outlined in the "Data availability" section of this manuscript. Once you are logged in, click "Browse Dataset Files". 
