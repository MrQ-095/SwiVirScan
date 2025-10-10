# This R code is provided under the MIT License
# Copyright (c) 2025 Zhang Q.
# Permission is hereby granted, free of charge, to use, copy, modify, and distribute this code for any purpose, provided that the above copyright notice appears in all copies.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND.

# utf-8
# Most code write on 20250719-20250720.
# R code for calculation on swine-Vscan data.
# Calculation only, NO GRAPH CODE.

#============================================
# Fig 2A. Vir accum data calculation with sample increase
#============================================
library(vegan) #ver 2.7-1
set.seed(114) # seeds setting
Virus_data <- read.csv("Fig2A_inputData.csv", 
                       header = T, sep = ",", 
                       row.names = 1)
accum_curve_100 <- specaccum(Virus_data, method = "random", 
                             permutations = 100) #random rearrange 100 times

accum_table_100 <- data.frame(
  Sites = accum_curve_100$sites, Richness = accum_curve_100$richness, SD = accum_curve_100$sd
)

#save as csv file, include virus number and its SD 
write.csv(accum_table_100, file = "Output_Fig2A_accum_result_table_seed114.csv", row.names = FALSE)

#============================================
# Fig 2. does Vir frequency significant difference in gender/parity group ?
# Input data(4*Viues_number Table). col name is virus, row name is fisher table row/col name (V1+Boar/V1+soar/V1-boar/V1-soar)
#============================================

# gender, output -lgP in row 5.
Genderdata <- read.csv("Fig2_Vir2BG_gender.csv", header = T, row.names = 1, sep = ",")

p_values_Gender <- apply(Genderdata, 2, function(col) {
  test_result <- fisher.test(matrix(col,nrow=2))
  return(-log10(test_result$p.value))
})

P_Gender <- rbind(Genderdata, p_values_Gender)
write.csv(P_Gender, file = "Output_Fig2Pgender.csv", row.names = T)


#parity, output -lgP in row 5.
Paritydata <- read.csv("Fig2_Vir2BG_Parity.csv", header = T, row.names = 1, sep = ",")

p_parity <- apply(Paritydata, 2, function(col) {
  test_result <- fisher.test(matrix(col,nrow=2))
  return(-log10(test_result$p.value))
})

ParityData <- rbind(Paritydata, p_parity)
write.csv(ParityData, file = "Output_Fig2Pparity.csv", row.names = T)


#============================================
# Fig 3. virus co-infection
#============================================

######### PART 1. observed co-infection 

# virus*sample, presence/absence Matrix(as 1/0 matrix), virus is row name.
VS_Data <- read.csv("Fig3_Interact.CSV", header = T, row.names = 1, sep = ",")

# observed co-infection frequency, based on 127 sampls
co_infection_P <- function(matrix) {
  num_viruses <- nrow(matrix)
  prob_matrix <- matrix(0, nrow = num_viruses, ncol = num_viruses)  #creates a matrix from the given set of values.
  rownames(prob_matrix) <- rownames(matrix)
  colnames(prob_matrix) <- rownames(matrix)
  
  for (a in 1:num_viruses) {
    for (b in 1:num_viruses) {
      if (a != b) {
        ab_count <- sum(matrix[a,] == 1 & matrix[b,] == 1)
        prob_matrix[a, b] <- ab_count /127   
      }else {
        prob_matrix[a, b] <- NA
      }
    }
  }
  return(prob_matrix)  
} 

#observed co-infection frequency matrix
Co_full <- co_infection_P(VS_Data)
Co_full <- Co_full*100
write.csv(Co_full,"Output_Fig3_Coinfection.csv",row.names = T)

#observed co-infection frequency matrix TO LONG TYPE
lower_idx <- which(lower.tri(Co_full), arr.ind = TRUE)
virus_names <- rownames(Co_full)
Co_full_long <- data.frame(
  Virus1 = virus_names[lower_idx[, 1]],
  Virus2 = virus_names[lower_idx[, 2]],
  Obs_Exp_Ratio = Co_full[lower.tri(Co_full)]
)

write.csv(Co_full_long, "Output_Fig3_CoinfectionLong.csv", row.names = FALSE)

######### PART 2. RATIO of observed/Theoretical co-infection

# Theoretical co-infection freq matrix
virus_freq <- rowSums(VS_Data) / ncol(VS_Data)
theo_matrix <- outer(virus_freq, virus_freq, FUN = "*") # freq Vir_X * Vir_Y
rownames(theo_matrix) <- names(virus_freq)
colnames(theo_matrix) <- names(virus_freq)
# observed/thore matrix
Co_full <- Co_full/100
ratio_matrix <- Co_full / theo_matrix
diag(ratio_matrix) <- NA
# output
write.csv(ratio_matrix, "Output_Fig3_ObservedExpectedRatio.csv",row.names = T)


get_long <- function(TODO_matrix) {
  lower_idx <- which(lower.tri(TODO_matrix), arr.ind = TRUE)
  virus_names <- rownames(TODO_matrix)
  data.frame(
    Virus1 = virus_names[lower_idx[, 1]],
    Virus2 = virus_names[lower_idx[, 2]],
    Obs_Exp_Ratio = TODO_matrix[lower.tri(TODO_matrix)]
  )
}
ratio_long <- get_long(ratio_matrix)
write.csv(ratio_long, "Output_Fig3_ObservedExpectedRatio_long.csv", row.names = FALSE)
# long data for interval freq stats.

######### PART 3. Independent fisher between each virus-pair with FDR

# Fisher test relative in two virus, Lower triangular trans to long form and use FDR.
Fisher_biVirus <- function(matrix) {
  num_viruses <- nrow(matrix)
  prob_matrix <- matrix(0, nrow = num_viruses, ncol = num_viruses)
  rownames(prob_matrix) <- rownames(matrix)
  colnames(prob_matrix) <- rownames(matrix)
  
  for (a in 1:num_viruses) {
    for (b in 1:num_viruses) {
      if (a != b) {
        m <- sum(matrix[a,] == 1 & matrix[b,] == 1)
        n <- sum(matrix[a,] == 1 & matrix[b,] == 0)
        k <- sum(matrix[a,] == 0 & matrix[b,] == 1)
        i <- sum(matrix[a,] == 0 & matrix[b,] == 0)
        test_result <- fisher.test(matrix(c(m,n,k,i),nrow=2))
        prob_matrix[a, b] <- test_result$p.value
      }else {
        prob_matrix[a, b] <- NA
      }
    }
  }
  return(prob_matrix)
} 

p_biVirus <- Fisher_biVirus(VS_Data)

# get long FDR data to get relative virus pair
virus_names <- rownames(p_biVirus)
lower_idx <- which(lower.tri(p_biVirus), arr.ind = TRUE)

p_long <- data.frame(
  Virus1 = virus_names[lower_idx[, 1]],
  Virus2 = virus_names[lower_idx[, 2]],
  p.value = p_biVirus[lower.tri(p_biVirus)])
# FDR and output csv file
p_long$FDR <- p.adjust(p_long$p.value, method = "fdr")
write.csv(p_long,"Output_Fig3_CoinfecFisherFDR.csv",row.names = T)

# use FDR's q-value to create -lg10(q-value) matrix. 
# new matrix to save fdr data
fdr_matrix <- p_biVirus

lower_idx <- which(lower.tri(p_biVirus), arr.ind = TRUE)

for (i in 1:nrow(lower_idx)) {
  row <- lower_idx[i, 1]
  col <- lower_idx[i, 2]
  fdr_value <- p_long$FDR[i]
  fdr_matrix[row, col] <- fdr_value
  fdr_matrix[col, row] <- fdr_value  
}

# get -log10(FDR) Matrix for heatmap
log10_fdr_matrix <- fdr_matrix
log10_fdr_matrix <- -log10(log10_fdr_matrix)
write.csv(log10_fdr_matrix, "Output_Fig3_Coinfec_fisherFDR_matrix.csv")


#============================================
# Fig 3-1. virus co-infection after sim-pep remove
# remove similary pep between relative virus pairs
#============================================
#=============================
# purpose: remove similary pep between relative virus pairs
# 1. For each sample s, identify the positive peptides.
# 2. For each active peptide p, check the status of its same-virus and cross-virus peptides in that sample:
#    if no same-virus peptide is active (same_hit = FALSE) and at least one cross-virus peptide is active (cross_hit = TRUE),
#    then set mat_new[p, s] = 0.
# 3. In all other cases, keep mat_new[p, s] unchanged.
#=============================

# (01 Mat，row is peptide，col is sample)
mat <- read.csv("pep2sampleMatrix.CSV", row.names = 1, check.names = FALSE)
mat_new <- mat

# Supplementary table: same-virus and cross-virus relationships (blastp outfmt 7 table)
sameVirus <- read.csv("RM_overlape.csv") #same
crossVirus <- read.csv("RM_heatmapPep.csv") #cross

# Build the mapping table to include both qseq -> sseq and sseq -> qseq.
same_map <- rbind(sameVirus[, c("qseqid", "sseqid")],
                  setNames(sameVirus[, c("sseqid", "qseqid")], c("qseqid", "sseqid")))
same_list <- split(same_map$sseqid, same_map$qseqid)

cross_map <- rbind(crossVirus[, c("qseqid", "sseqid")],
                   setNames(crossVirus[, c("sseqid", "qseqid")], c("qseqid", "sseqid")))
cross_list <- split(cross_map$sseqid, cross_map$qseqid)

# potential pep RM
for (s in colnames(mat)) {
  active_peps <- rownames(mat)[mat[[s]] == 1]
  
  for (p in active_peps) {
    same_peps <- same_list[[p]] 
    if (is.null(same_peps)) same_peps <- character(0) 

        cross_peps <- cross_list[[p]]
    if (is.null(cross_peps)) cross_peps <- character(0)
    
    same_hit <- any(mat[same_peps, s, drop = FALSE] == 1) 
    cross_hit <- any(mat[cross_peps, s, drop = FALSE] == 1)
    if (!same_hit && cross_hit) {mat_new[p, s] <- 0} 
  }
}

write.csv(mat_new, "Output_Fig3_main_matrix_updated.csv")

# get long data.frame function
# pep2sample file("Output_Fig3_main_matrix_updated.csv") trans to vir2sample file "RM_VIr2sampleMatrix.CSV" in EXCEL by "sumif" and "if"
VS_RM <- read.csv("RM_VIr2sampleMatrix.CSV", row.names = 1, check.names = FALSE)
p_biVirus_RM <- Fisher_biVirus(VS_RM) # fisher-P martix

# long data for FDR
biVirusP_get_long <- function(TODO_matrix) {
  lower_idx <- which(lower.tri(TODO_matrix), arr.ind = TRUE)
  virus_names <- rownames(TODO_matrix)
  data.frame(
    Virus1 = virus_names[lower_idx[, 1]],
    Virus2 = virus_names[lower_idx[, 2]],
    p_value = TODO_matrix[lower.tri(TODO_matrix)]
  )
}
RmFisherP_long <- biVirusP_get_long(p_biVirus_RM)
RmFisherP_long$FDR <- p.adjust(RmFisherP_long$p_value, method = "fdr")

# new matrix to save fdr data
fdrRM_matrix <- matrix(NA, nrow = nrow(p_biVirus_RM), ncol = ncol(p_biVirus_RM),
                       dimnames = dimnames(p_biVirus_RM))
lower_idx <- which(lower.tri(p_biVirus_RM), arr.ind = TRUE)
for (i in 1:nrow(lower_idx)) {
  row <- lower_idx[i, 1]
  col <- lower_idx[i, 2]
  fdr_value <- RmFisherP_long$FDR[i]
  fdrRM_matrix[row, col] <- fdr_value
  fdrRM_matrix[col, row] <- fdr_value  
}
log10_fdrRM_matrix <- -log10(fdrRM_matrix)


write.csv(RmFisherP_long, "Output_Fig3_RmFisherP_long.csv", row.names = FALSE)
write.csv(log10_fdrRM_matrix, "Output_Fig3_RmFisher_Matrix.csv", row.names = FALSE)


#============================================
# Fig 4. does virus_X freq difference in (VacX+ and vacX-) bacground
# remove similary pep between relative virus pairs
#============================================

# col numbers for fisher, col name are virus, row name are (VacX+ and vacX-) bacground cor with (virus+-), 4 row each Vaccine.
VacData <- read.csv("Fig4_Vac.CSV", header = T, row.names = 1)

TGEV <- VacData[c(1:4),]
FMDV <- VacData[c(5:8),]
CSFV <- VacData[c(9:12),]

p_TGEV <- apply(TGEV, 2, function(col) {
  test_result <- fisher.test(matrix(col,nrow=2))
  return(-log10(test_result$p.value))})
p_FMDV <- apply(FMDV, 2, function(col) {
  test_result <- fisher.test(matrix(col,nrow=2))
  return(-log10(test_result$p.value))})
p_CSFV <- apply(CSFV, 2, function(col) {
  test_result <- fisher.test(matrix(col,nrow=2))
  return(-log10(test_result$p.value))})

# output as 13-15 row.
VacData["TGEV_pval", ] <- p_TGEV
VacData["FMDV_pval", ] <- p_FMDV
VacData["CSFV_pval", ] <- p_CSFV

write.csv(VacData, file = "Output_Fig4Vac.csv", row.names = T)
