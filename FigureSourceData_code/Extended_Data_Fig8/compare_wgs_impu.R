#Read in the WGS-IMPU table and create the plot
setwd("~/Documents/PRS_Jason/robustness_analysis/NA12878_validation/Final_result/Individual_sample_comparison")
minus_scores <- read.table("~/Documents/PRS_Jason/robustness_analysis/NA12878_validation/Final_result/Individual_sample_comparison/adjusted_standardized_minus.csv", sep=",", header=T, stringsAsFactors = F)

#plot(density(minus_scores[,2], adjust=2), col="black", xlim=c(-1.5, 1.5), ylim=c(0, 6.0), xlab="", ylab="", lty="solid", main="")

par(mar=c(3.5,3.5,1,1))
plot(density(minus_scores[,2], adjust=2), col="black", xlim=c(-1.5, 1.7), ylim=c(0, 2.8), xlab="", ylab="", lty="solid", main="")
title(xlab="Standardized adjusted PRS score difference (WGS-IMPU)", line=2, cex.lab=1.0)
title(ylab="Density", line=2, cex.lab=1)

lines(density(minus_scores[,3], adjust=2, na.rm=T), col="blue", xlim=c(-2, 2), ylim=c(0, 3.5), xlab="", ylab="", lty="solid", main="")
lines(density(minus_scores[,4], adjust=2), col="red", xlim=c(-2, 2), ylim=c(0, 3.5), xlab="", ylab="", lty="solid", main="")
lines(density(minus_scores[,5], adjust=2), col="purple", xlim=c(-2, 2), ylim=c(0, 3.5), xlab="", ylab="", lty="solid", main="")
lines(density(minus_scores[,6], adjust=2, na.rm=T), col="cyan", xlim=c(-2, 2), ylim=c(0, 3.5), xlab="", ylab="", lty="solid", main="")
lines(density(minus_scores[,7], adjust=2), col="green", xlim=c(-2, 2), ylim=c(0, 6.5), xlab="", ylab="", lty="solid", main="")

#Vars <- c("AF", "BC", "CAD", "CC", "PC", "T2D")  # one per row
#colors <- c("black", "blue", "red", "purple", "cyan", "green")  # one color for each row; should be same length as Vars
Vars <- c("BrCa", "CRCa", "PrCa", "AFib", "CAD", "T2D")  # one per row
colors <- c("blue", "purple", "cyan", "black", "red", "green")  # one color for each row; should be same length as Vars
ltys <- 1  # one lty for each column

nr <- length(Vars) 
nc <- length(ltys)
legend("topright", rep(Vars, nc), col = colors, lty = rep(ltys, each = nr),
       ncol = nc, cex = 0.8)


