#####################################################################################
######## LIMMA script for ONE (1), TWO (2) or THREE (3) group comparison ############
########             by Alejandro Gomez-Auli 2015-2017 v1.4.1            ############
#####################################################################################
### 
## This script is based on: 
## Gomez-Auli A, Hillebrand LE, Biniossek ML, Peters C, Reinheckel T, Schilling O. Impact of cathepsin B on the interstitial 
## fluid proteome of murine breast cancers. Biochimie. 2016 Mar;122:88-98.
##
#### The script should be in the same folder as the file that is going to be analyzed,
#### just make a copy, and leave the original in its place
## For the script to work the file should: 
##
## 1. A txt tab-separated file
## 2. Be called "input_limma.txt"
## 3. Have a heading
## 4. The first column should be the ID, Uniprot, gene, etc and must be called "ID" (watch the case)
## 5. The expression, ratio values (log transformed ) should be located in the next columns
##
## As an example the structure should be:
##
##  ID      exp1  exp2  exp3  exp4
##  P10605  0.2   0.3   NA  0.25
##
## It does not matter if there are missing values, although of course one should have at least 2 values for any statistical
## test
###
###########
## Open the script in RStudio (eg. double-click)
## Press the Source button and answer the questions
#####

## Set working directory as the one where the script is located
dir <- dirname(parent.frame(2)$ofile)
setwd(dir)

#install.packages("bio3d")

#Load required libraries
#### LIMMA

if(suppressWarnings(suppressMessages(require(limma)))) {
  message('limma loaded correctly')
} else {
  source("https://bioconductor.org/biocLite.R")
  biocLite("limma")
  suppressWarnings(suppressMessages(library(limma)))
}

#### RColorbrewer

if(suppressWarnings(suppressMessages(require(RColorBrewer)))) {
  message('RColorBrewer loaded correctly')
} else {
  install.packages("RColorBrewer")
  suppressWarnings(suppressMessages(library(RColorBrewer)))
}

#### dplyr
if(suppressWarnings(suppressMessages(require(dplyr)))) {
  #message('dplyr loaded correctly')
} else {
  install.packages("dplyr")
  suppressWarnings(suppressMessages(library(dplyr)))
}

if(suppressWarnings(suppressMessages(require(bio3d)))) {
  #message('bio3d loaded correctly')
} else {
  install.packages("bio3d")
  suppressWarnings(suppressMessages(library(bio3d)))
}


#Assign the ID as row name to have better control of the names
input <- read.delim("input_limma.txt", na.strings=c("0", "", "NA", "#NUM!", "#ZAHL!", "NaN", "-1.#INF", "1.#INF"))

#Prevent ID not being properly named
colnames(input)[1] <- "ID"

## Copy  the data frame
input2 <- input

## Use ID as row names
rownames(input2) <- input2$ID
## Subset the fold changes, ie. all columns after the first one, which should have the IDs
input2 <- input2[-1]

## Summary
summary(input2)
## Ask how many groups and store the variable

g <- readline("How many groups do you have? ")
g <- as.numeric(g)

p_cut <- readline("p-value cutoff for graphics; e.g. 0.01 or 0.05 ")
p_cut <- as.numeric(p_cut)

ratio_cut <- readline("ratio cutoff for graphics; e.g. 1.5 for a minimum 50 % difference ")
ratio_cut <- as.numeric(ratio_cut)

get_names <- readline("Get full protein names for Uniptot IDs ? y/n")

## For 1 group

if (g == 1) {
  
  ## Plot boxplot of intensities, fold changes
  png("boxplot_fc_one_group.png", width = 1600, height = 1600)
  par(mar=c(12,18,4,1)+.1, mgp = c(10, 3, 0))
  boxplot(input2, las=1, pch=1, cex.lab=5, cex.axis=5, cex=5, bty="n", boxlwd = 5, whisklwd = 5, ylab="log2 of ratio", col=rev(brewer.pal(9, "Blues")), names=c(1:ncol(input2)))
   dev.off()
   
  #Fit the linear model
  fit_perseus <- lmFit(input2,method="robust",maxit=1000)

  
  #Generating the moderated statistics Using proportion 0.05, assumed number of proteins to be 
  #differentiated regulated
  fit_perseus2 <- eBayes(fit_perseus, proportion=0.05,robust = TRUE, trend = TRUE)
  
  
  
  
  #Ideally this would be the next step to obtain the moderated t-test, but the F moderated p-value also gives good information
  #moderated for multiple hypothesis testing. The confint option gives 95 CI which also provides good info
  #It also provides an extra step adjusting the p-value obtained by the moderated t-test by BH, although the p-value obtained, 
  #should be already be somehow corrected
  r_fit_perseus <- topTable(fit_perseus2, confint=TRUE, resort.by="logFC", number = nrow(fit_perseus2))
  
  #qq Plots of the differentially expressed proteins
  png("qqplot.png")
    qqt(fit_perseus2$t,df=fit_perseus2$df.prior+fit_perseus2$df.residual, pch=19, cex=0.5, las=1)
    abline(0,1, lty=2)
  dev.off()
  #Re-arrange
  #Add the ID again
  r_fit_perseus$ID <- rownames(r_fit_perseus)
  names(r_fit_perseus)
  r_fit_perseus <- r_fit_perseus[c(9,1:8)]
  
  #Merge with original file to get the ratios of every experiment
  results <- merge(input, r_fit_perseus, by="ID")
  
  #Generate volcano plot and draw lines for a p-value 0.01, and log2 +/- 50%
  dablack <- rgb(t(col2rgb("#474747")), alpha=120, maxColorValue=255)
  png("Volcano_plot_old.png", width=800, height=800)
    plot(results$logFC, -log10(results$P.Value), xlab=expression(paste("Fold change", "(",~log[2],")")), 
         ylab=expression(paste(-log[10], " (Moderated p-value)")), las=1, col=dablack, pch=19, cex=2, cex.lab=1.3)
    abline(h=-log10(0.05), lty=2, col="gray35")
    abline(v=c(-1, 1)*log(1.5, 2), lty=2, col="gray35")
  dev.off()
  
  l1 <- max(abs(range(results$logFC)))
  dablack <- rgb(t(col2rgb("#474747")), alpha=120, maxColorValue=255)
  #
  max_ratio <-c(abs(max(results$logFC)),abs(min(results$logFC)))
  max_ratio <-max(max_ratio)
  
  png("volcano_moderated_p-value.png", width = 1600, height = 1600)
  min_p_value <- min(results$P.Value, na.rm = TRUE)
  min_p_value <- -log10(min_p_value)
  par(mar=c(12,18,4,1)+.1, mgp = c(10, 3, 0))
  plot(results$logFC, -log10(results$P.Value),  xlim=c(-max_ratio,max_ratio),ylim=c(0,min_p_value), las=1, pch=19, cex.lab=5, cex.axis=5, col=dablack, cex=5, bty="n",  xlab=paste("mean log2 of ratio"), ylab="-log10(limma moderated p-value)")
  ## Draw line of p value cutoff
  abline(h=-log10(p_cut), lty=2, col="gray46", lwd = 4)
  abline(v=log2(ratio_cut)*c(-1,1), lty=2, col="gray46", lwd = 4)
  axis(side = 1, lwd = 3, labels=FALSE)
  axis(side = 2, lwd = 3, labels=FALSE)
  dev.off()
  
  png("volcano_adjusted_p-value.png", width = 1600, height = 1600)
  min_p_value <- min(results$adj.P.Val, na.rm = TRUE)
  min_p_value <- -log10(min_p_value)
  par(mar=c(12,18,4,1)+.1, mgp = c(10, 3, 0))
  plot(results$logFC, -log10(results$adj.P.Val),  xlim=c(-max_ratio,max_ratio),ylim=c(0,min_p_value), las=1, pch=19, cex.lab=5, cex.axis=5, col=dablack, cex=5, bty="n",  xlab=paste("mean log2 of ratio"), ylab="-log10(limma adjusted p-value)")
  ## Draw line of p value cutoff
  abline(h=-log10(p_cut), lty=2, col="gray46", lwd = 4)
  abline(v=log2(ratio_cut)*c(-1,1), lty=2, col="gray46", lwd = 4)
  axis(side = 1, lwd = 3, labels=FALSE)
  axis(side = 2, lwd = 3, labels=FALSE)
  dev.off()

  
  #Write table with limma results
  write.table(results, "limma_results_preliminary.txt", sep="\t", row.names = FALSE)
  
  t3filtered <- results
  for(r in 1:nrow(t3filtered))
  { 
    moderated_p_value <- (t3filtered$P.Value[r])
    adjusted_p_value <- (t3filtered$adj.P.Val[r])
    ratio_value <- (t3filtered$logFC[r])
    if (moderated_p_value < p_cut & ratio_value > log2(ratio_cut)) t3filtered [r, "inc_mpv"] <- "1"
    if (moderated_p_value < p_cut & ratio_value < -log2(ratio_cut)) t3filtered [r, "dec_mpv"] <- "1"
    if (adjusted_p_value < p_cut & ratio_value > log2(ratio_cut)) t3filtered [r, "inc_apv"] <- "1"
    if (adjusted_p_value < p_cut & ratio_value < -log2(ratio_cut)) t3filtered [r, "dec_apv"] <- "1"
  } 
  t3filtered_out <- t3filtered
  colnames(t3filtered_out)[colnames(t3filtered)=="logFC"] <- paste("mean log2 of ratio")
  colnames(t3filtered_out)[colnames(t3filtered)=="inc_mpv"] <- paste("increased with moderated p-value")
  colnames(t3filtered_out)[colnames(t3filtered)=="dec_mpv"] <- paste("decreased with moderated p-value")
  colnames(t3filtered_out)[colnames(t3filtered)=="inc_apv"] <- paste("increased with adjusted p-value")
  colnames(t3filtered_out)[colnames(t3filtered)=="dec_apv"] <- paste("decreased with adjusted p-value")
  write.table( t3filtered_out, "tmp_results.tsv", sep="\t", row.names = FALSE)
  
  for(r in 1:nrow(t3filtered))
  { 
    moderated_p_value <- (t3filtered$P.Value[r])
    adjusted_p_value <- (t3filtered$adj.P.Val[r])
    ratio_value <- (t3filtered$logFC[r])
    if (moderated_p_value < p_cut & ratio_value > log2(ratio_cut)) t3filtered [r, "inc_mpv"] <- "1"
    if (moderated_p_value < p_cut & ratio_value < -log2(ratio_cut)) t3filtered [r, "dec_mpv"] <- "1"
    if (adjusted_p_value < p_cut & ratio_value > log2(ratio_cut)) t3filtered [r, "inc_apv"] <- "1"
    if (adjusted_p_value < p_cut & ratio_value < -log2(ratio_cut)) t3filtered [r, "dec_apv"] <- "1"
    if(get_names == "y")
    {
      prot <- uniprot(t3filtered$ID[r])
      t3filtered [r, "Full name"] <- prot$fullName
      t3filtered [r, "Short name"] <- prot$shortName
      t3filtered [r, "Gene ID"] <- prot$gene
    } 
    
    t3filtered_out <- t3filtered
    colnames(t3filtered_out)[colnames(t3filtered)=="logFC"] <- paste("mean log2 of ratio")
    colnames(t3filtered_out)[colnames(t3filtered)=="inc_mpv"] <- paste("increased with moderated p-value")
    colnames(t3filtered_out)[colnames(t3filtered)=="dec_mpv"] <- paste("decreased with moderated p-value")
    colnames(t3filtered_out)[colnames(t3filtered)=="inc_apv"] <- paste("increased with adjusted p-value")
    colnames(t3filtered_out)[colnames(t3filtered)=="dec_apv"] <- paste("decreased with adjusted p-value")
    write.table( t3filtered_out, "limma_results_final.tsv", sep="\t", row.names = FALSE)
  }  
  
  print("Analysis done, look for the results file --limma_results.txt--, in the folder where the input file is located.")
  

## End for one group  
  
## Start two groups
   
} else if (g == 2) {
  
  group1 <- readline("How many columns are in the first group? ")
  group1 <- as.numeric(group1)
  name_group1 <- readline(prompt="Short name for group 1 ")
  
  group2 <- readline("How many columns are in the second group? ")
  group2 <- as.numeric(group2)
  name_group2 <- readline(prompt="Short name for group 2 ")
  
  ## Subset for at least 2 values
  
  selecc <- rowSums(is.na(input2[1:group1])) < group1-1 & rowSums(is.na(input2[(group1+1):(group1+group2)])) < group2-1 
  
  ## Subset these
  input3 <- input2[selecc,]
  ## Summary
  summary(input3)
  #Number of proteins left
  nrow(input3)
  #
  

## Plot boxplot of intensities, fold changes of two groups 

png("boxplot_fc_two_groups.png", width = 1600, height = 1600)
par(mar=c(12,18,4,1)+.1, mgp = c(10, 3, 0))
boxplot(input3, las=1, pch=1, cex.lab=5, cex.axis=5, cex=5, bty="n", boxlwd = 5, whisklwd = 5, ylab=paste("mean log2(",name_group2,"/",name_group1,")"), col=rev(brewer.pal(3, "Blues")), names=c(1:group1, 1:group2),at=c(1:group1, (group1+1.5):(ncol(input3)+0.5)))
dev.off()


  
   
  #Fit a linera model
  ## As there are groups, one needs to define a design
  
  #Create matrix
  g1 <- c(rep(1,group1), rep(0,group2))
  g2 <- c(rep(0,group1), rep(1,group2))
  
  design <- cbind(g1=g1, g2=g2)
  design
  
  ## Create contrast matrix
  cont.matrix <- makeContrasts(g2vsg1=g2-g1, levels=design)
  
  ## Fit linear model
  fit <- lmFit(input3, design,method="robust",maxit=1000)
  
  ## Create contrast of groups
  fit2 <- contrasts.fit(fit, cont.matrix)
  
  ## Do moderated t-test
  fit2 <- eBayes(fit2, proportion = 0.05, robust = TRUE, trend = TRUE)
  
  ## Top tables with coeficient for each comparison
  topTable(fit2, adjust="fdr")
  
  ## Result matrix
  ## First comparison
  r_fit_g1 <- topTable(fit2, confint=TRUE, resort.by="logFC", number = nrow(fit2), coef = 1)
  
  names(r_fit_g1)
  ## Merge results
  
  #Re-arrange
  #Add the ID again
  r_fit_g1$ID <- rownames(r_fit_g1)
  names(r_fit_g1)
  
  ## Reorder to have the ID as first column
  t2 <- r_fit_g1[c(ncol(r_fit_g1),1:(ncol(r_fit_g1)-1))]
  
  ## Merge with original file to get the RAW fold-change/intensities
  t3 <- merge(input, t2, by="ID")
  
  ## Merge with original file to get the RAW fold-change/intensities, including the empty values = less than two per group
  
#  t4 <- merge(input, t2, by="ID", all=TRUE)
  
  ## Write the final table
  drop <- c("AveExpr","t","B")
  t3filtered = t3[,!(names(t3) %in% drop)]
  
#  for (i in t3filtered$ID)
#  {
#    prot <- uniprot(i)
#    prot$fullName
#  }

for(r in 1:nrow(t3filtered))
{ 
  moderated_p_value <- (t3filtered$P.Value[r])
  adjusted_p_value <- (t3filtered$adj.P.Val[r])
  ratio_value <- (t3filtered$logFC[r])
  if (moderated_p_value < p_cut & ratio_value > log2(ratio_cut)) t3filtered [r, "inc_mpv"] <- "1"
  if (moderated_p_value < p_cut & ratio_value < -log2(ratio_cut)) t3filtered [r, "dec_mpv"] <- "1"
  if (adjusted_p_value < p_cut & ratio_value > log2(ratio_cut)) t3filtered [r, "inc_apv"] <- "1"
  if (adjusted_p_value < p_cut & ratio_value < -log2(ratio_cut)) t3filtered [r, "dec_apv"] <- "1"
#  if(get_names == "y")
#  {
#    prot <- uniprot(t3filtered$ID[r])
#    t3filtered [r, "Full name"] <- prot$fullName
#    t3filtered [r, "Short name"] <- prot$shortName
#    t3filtered [r, "Gene ID"] <- prot$gene
#  }
  
} 
  
t3filtered_out <- t3filtered
colnames(t3filtered_out)[colnames(t3filtered)=="logFC"] <- paste("mean log2(",name_group2,"/",name_group1,")")
colnames(t3filtered_out)[colnames(t3filtered)=="inc_mpv"] <- paste("increased in", name_group2,"(moderated p-value)")
colnames(t3filtered_out)[colnames(t3filtered)=="dec_mpv"] <- paste("decreased in", name_group2,"(moderated p-value)")
colnames(t3filtered_out)[colnames(t3filtered)=="inc_apv"] <- paste("increased in", name_group2,"(adjusted p-value)")
colnames(t3filtered_out)[colnames(t3filtered)=="dec_apv"] <- paste("decreased in", name_group2,"(adjusted p-value)")
write.table( t3filtered_out, "tmp_results.tsv", sep="\t", row.names = FALSE)

  
    
# Write the final table with the empty values too
#  write.table(t4, "limma_results_with_empty_values.txt", sep="\t", row.names = FALSE)

  ## Plot volcano plot per coef
  ## Group 2 vs 1
  ## xlim
  l1 <- max(abs(range(t3$logFC)))
  dablack <- rgb(t(col2rgb("#474747")), alpha=120, maxColorValue=255)
  #
  max_ratio <-c(abs(max(t3$logFC)),abs(min(t3$logFC)))
  max_ratio <-max(max_ratio)
  
  png("Group2_Group1_volcano_moderated_p-value.png", width = 1600, height = 1600)
  min_p_value <- min(t3$P.Value, na.rm = TRUE)
  min_p_value <- -log10(min_p_value)
  par(mar=c(12,18,4,1)+.1, mgp = c(10, 3, 0))
  plot(t3$logFC, -log10(t3$P.Value),  xlim=c(-max_ratio,max_ratio),ylim=c(0,min_p_value), las=1, pch=19, cex.lab=5, cex.axis=5, col=dablack, cex=5, bty="n",  xlab=paste("mean log2(",name_group2,"/",name_group1,")"), ylab="-log10(limma moderated p-value)")
  ## Draw line of p value cutoff
  abline(h=-log10(p_cut), lty=2, col="gray46", lwd = 4)
  abline(v=log2(ratio_cut)*c(-1,1), lty=2, col="gray46", lwd = 4)
  axis(side = 1, lwd = 3, labels=FALSE)
  axis(side = 2, lwd = 3, labels=FALSE)
  dev.off()

  png("Group2_Group1_volcano_adjusted_p-value.png", width = 1600, height = 1600)
  min_p_value <- min(t3$adj.P.Val, na.rm = TRUE)
  min_p_value <- -log10(min_p_value)
  par(mar=c(12,18,4,1)+.1, mgp = c(10, 3, 0))
  plot(t3$logFC, -log10(t3$adj.P.Val),  xlim=c(-max_ratio,max_ratio),ylim=c(0,min_p_value), las=1, pch=19, cex.lab=5, cex.axis=5, col=dablack, cex=5, bty="n",  xlab=paste("mean log2(",name_group2,"/",name_group1,")"), ylab="-log10(limma adjusted p-value)")
  ## Draw line of p value cutoff
  abline(h=-log10(p_cut), lty=2, col="gray46", lwd = 4)
  abline(v=log2(ratio_cut)*c(-1,1), lty=2, col="gray46", lwd = 4)
  axis(side = 1, lwd = 3, labels=FALSE)
  axis(side = 2, lwd = 3, labels=FALSE)
  dev.off()
  
  #Generate volcano plot with colors of dysregulated proteins and draw lines for a p-value 0.01, and log2 +/- 50%
  # Create a new data frame adding if it is dysregulated, using mutate and if/else
  #suppressWarnings(library(dplyr))

  
  t4 <- t3 %>%
    mutate(color = ifelse(t3$logFC > log2(ratio_cut) & t3$P.Value < p_cut, 
                          yes = "More", 
                          no = ifelse(t3$logFC < -1*log2(ratio_cut) & t3$P.Value < p_cut, 
                                      yes = "Less", 
                                      no = "none")))
  
  ## Color it and plot but with base R
  dagray <- rgb(t(col2rgb("#636363")), alpha=120, maxColorValue=255)
  #dared <- rgb(t(col2rgb("#b2182b")), alpha=120, maxColorValue=255)
  #dablue <- rgb(t(col2rgb("#053061")), alpha=120, maxColorValue=255)
  
  dared <- rgb(t(col2rgb("#006400")), alpha=120, maxColorValue=255)
  dablue <- rgb(t(col2rgb("#006400")), alpha=120, maxColorValue=255)
  
  
  ## Change factor to the colors
  t4$color2 <- t4$color
  ## Replace More, less, none, for colors
  t4$color2[t4$color2 == "More"] <- dared #"#b2182b" #dared
  t4$color2[t4$color2 == "Less"] <- dablue # "#053061" # dablue #
  t4$color2[t4$color2 == "none"] <- dagray #"#636363"
  ## Plot
  ## Range xlim
  l1 <- max(abs(range(t4$logFC)))
  ## Draw
  png("Group2_Group1_volcano_colors_moderated_p-value.png", width = 1600, height = 1600)
  min_p_value <- min(t3$P.Value, na.rm = TRUE)
  min_p_value <- -log10(min_p_value)
  par(mar=c(12,18,4,1)+.1, mgp = c(10, 3, 0))
  plot(t4$logFC, -log10(t4$P.Value), xlim=c(-max_ratio,max_ratio),ylim=c(0,min_p_value), las=1, pch=19, las=1, cex.lab=5, cex.axis=5, col=t4$color2, cex=5, bty="n",  xlab=paste("mean log2(",name_group2,"/",name_group1,")"), ylab="-log10(limma moderated p-value)")
  ## Draw line of p value cutoff
#  abline(h=-log10(0.05), lty=2, col="gray46", lwd = 2)
  ## Draw line of fold change cutoff +/- 0.
  abline(h=-log10(p_cut), lty=2, col="gray46", lwd = 4)
  abline(v=log2(ratio_cut)*c(-1,1), lty=2, col="gray46", lwd = 4)
  axis(side = 1, lwd = 3, labels=FALSE)
  axis(side = 2, lwd = 3, labels=FALSE)
  dev.off()
  
  
  t4 <- t3 %>%
    mutate(color = ifelse(t3$logFC > log2(ratio_cut) & t3$adj.P.Val < p_cut, 
                          yes = "More", 
                          no = ifelse(t3$logFC < -1*log2(ratio_cut) & t3$adj.P.Val < p_cut, 
                                      yes = "Less", 
                                      no = "none")))
  
  ## Color it and plot but with base R
  dagray <- rgb(t(col2rgb("#636363")), alpha=120, maxColorValue=255)
  #dared <- rgb(t(col2rgb("#b2182b")), alpha=120, maxColorValue=255)
  #dablue <- rgb(t(col2rgb("#053061")), alpha=120, maxColorValue=255)
  
  dared <- rgb(t(col2rgb("#006400")), alpha=120, maxColorValue=255)
  dablue <- rgb(t(col2rgb("#006400")), alpha=120, maxColorValue=255)
  
  
  ## Change factor to the colors
  t4$color2 <- t4$color
  ## Replace More, less, none, for colors
  t4$color2[t4$color2 == "More"] <- dared #"#b2182b" #dared
  t4$color2[t4$color2 == "Less"] <- dablue # "#053061" # dablue #
  t4$color2[t4$color2 == "none"] <- dagray #"#636363"
  ## Plot
  ## Range xlim
  l1 <- max(abs(range(t4$logFC)))
  ## Draw
  png("Group2_Group1_volcano_colors_adjusted_p-value.png", width = 1600, height = 1600)
  min_p_value <- min(t3$adj.P.Val, na.rm = TRUE)
  min_p_value <- -log10(min_p_value)
  par(mar=c(12,18,4,1)+.1, mgp = c(10, 3, 0))
  plot(t4$logFC, -log10(t4$adj.P.Val), xlim=c(-max_ratio,max_ratio),ylim=c(0,min_p_value), las=1, pch=19, las=1, cex.lab=5, cex.axis=5, col=t4$color2, cex=5, bty="n",  xlab=paste("mean log2(",name_group2,"/",name_group1,")"), ylab="-log10(limma adjusted p-value)")
  ## Draw line of p value cutoff
#  abline(h=-log10(0.05), lty=2, col="gray46", lwd = 2)
  ## Draw line of fold change cutoff +/- 0.
  abline(h=-log10(p_cut), lty=2, col="gray46", lwd = 4)
  abline(v=log2(ratio_cut)*c(-1,1), lty=2, col="gray46", lwd = 4)
  axis(side = 1, lwd = 3, labels=FALSE)
  axis(side = 2, lwd = 3, labels=FALSE)
  dev.off()

for(r in 1:nrow(t3filtered))
{ 
  moderated_p_value <- (t3filtered$P.Value[r])
  adjusted_p_value <- (t3filtered$adj.P.Val[r])
  ratio_value <- (t3filtered$logFC[r])
  if (moderated_p_value < p_cut & ratio_value > log2(ratio_cut)) t3filtered [r, "inc_mpv"] <- "1"
  if (moderated_p_value < p_cut & ratio_value < -log2(ratio_cut)) t3filtered [r, "dec_mpv"] <- "1"
  if (adjusted_p_value < p_cut & ratio_value > log2(ratio_cut)) t3filtered [r, "inc_apv"] <- "1"
  if (adjusted_p_value < p_cut & ratio_value < -log2(ratio_cut)) t3filtered [r, "dec_apv"] <- "1"
  if(get_names == "y")
  {
    prot <- uniprot(t3filtered$ID[r])
    t3filtered [r, "Full name"] <- prot$fullName
    t3filtered [r, "Short name"] <- prot$shortName
    t3filtered [r, "Gene ID"] <- prot$gene
  } 
  
  t3filtered_out <- t3filtered
  colnames(t3filtered_out)[colnames(t3filtered)=="logFC"] <- paste("mean log2(",name_group2,"/",name_group1,")")
  colnames(t3filtered_out)[colnames(t3filtered)=="inc_mpv"] <- paste("increased in", name_group2,"(moderated p-value)")
  colnames(t3filtered_out)[colnames(t3filtered)=="dec_mpv"] <- paste("decreased in", name_group2,"(moderated p-value)")
  colnames(t3filtered_out)[colnames(t3filtered)=="inc_apv"] <- paste("increased in", name_group2,"(adjusted p-value)")
  colnames(t3filtered_out)[colnames(t3filtered)=="dec_apv"] <- paste("decreased in", name_group2,"(adjusted p-value)")
  write.table( t3filtered_out, "limma_results.tsv", sep="\t", row.names = FALSE)
}
  print("Analysis done, look for the results file --limma_results.txt--, , in the folder where the input file is located. Some IDs might be missing as the analysis requires at least two valid values per group.")
  #Press p and enter to see an experimental interactive volcano plot
  
## End of group 2
} else
  print("Error: you have not selected a valid option! Please check the number of groups and try again")

