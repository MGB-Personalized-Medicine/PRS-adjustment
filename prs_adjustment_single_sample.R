#!/bin/env R

args <- commandArgs(trailingOnly=TRUE)

prs_adjustment_single_sample <- function(sample_raw_score, sample_raw_pc, eigenval, pc_model){
  #Adjust PRS raw score by PCs from MGB Biobank
  #The sample_raw_pc, eigenval, pcmodel all derive from MGB Biobank
  #The sample_raw_score is adjusted by 4 PCs and the adjusted score is output into a file
  
  # Create an output file name based on the input sample_raw_score file
  output_name <- paste0(sample_raw_score, ".adjusted")

  # Read in all input files. Note default comment.char is "#", which cannot handle file name with "#"
  sample_raw_score <- read.table(sample_raw_score, comment.char="", header=T, stringsAsFactor=F)
  sample_raw_pc <- read.table(sample_raw_pc, comment.char="", header=T, stringsAsFactor=F)

  eigenval<- read.table(eigenval, header=F, stringsAsFactor=F)
  pc_model <- readRDS(pc_model)

  # Scale 4 PCs of the sample 
  sample_scaled_pc <- sample_raw_pc
  sample_scaled_pc[1,5:8] <- (sample_raw_pc[1,5:8] / sqrt(eigenval[,1])) * 2 
  colnames(sample_scaled_pc) <- c("FID", "IID", "NmissAlleleCT", "NamedAlleleDosageSum", "PC1", "PC2", "PC3", "PC4")

  # Predict sample scores based on the linear model
  pred_score <- predict(pc_model, newdata=sample_scaled_pc)
  
  #Adjustment
  adjusted_score <- sample_raw_score
  adjusted_score[,6] <- adjusted_score[,6] - pred_score
  colnames(adjusted_score) <- c("#FID", "IID", "NMISS_ALLELE_CT", "NAMED_ALLELE_DOSAGE_SUM", "SCORE1_AVG", "SCORE1_AVG_ADJ")
 
  #Write adjusted scores in a file
  write.table(adjusted_score, file=output_name, sep="\t", quote=F, row.names=F, col.names=T)
}

#Handle arguments
if(length(args) != 4){
  stop("Four files must be provided as below: \nsample_raw_score sample_raw_pc eigenval pc_model\n", call=FALSE)
}else if(length(args) == 4){
  sample_raw_score = args[1]
  sample_raw_pc = args[2]
  eigenval = args[3]
  pc_model = args[4]
}

#Call the function
prs_adjustment_single_sample(sample_raw_score, sample_raw_pc, eigenval, pc_model)
