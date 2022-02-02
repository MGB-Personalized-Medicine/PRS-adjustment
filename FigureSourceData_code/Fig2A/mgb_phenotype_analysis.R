library(UpSetR)
library(ComplexHeatmap)
library(tidyverse)


setwd("~/Documents/PRS_Jason/prs_analysis/phenotype_analysis")

#Read in phenotype data
pheno <- read.csv("mgb_phenotype_race.txt", sep="\t", header=T, stringsAsFactor=F)
pheno$Subject <- as.character(pheno$Subject)  #Changer subject from num to chr
names(pheno) <- c("Subject",  "Race", "BC", "CAD", "T2D", "AF", "CC", "PC")
pheno <- pheno[,c(1,2,6,3,4,7,8,5)]

#Change the races Other or Unknown into "Other/Unknown"
pheno[which(pheno$Race == "Unknown" | pheno$Race == "Other"),2] <- "Other/Unknown"

upset(pheno, sets=c("AF", "BC", "CAD", "CC", "PC", "T2D"),nintersects=45)

#Upset plot
upset(pheno, sets=c("AF", "BC", "CAD", "CC", "PC", "T2D"), 
      nintersects=45,
      order.by = "freq", 
      decreasing = TRUE, 
      mb.ratio = c(0.6, 0.4),
      number.angles = 30, 
      text.scale = 1.5, 
      point.size = 2, 
      line.size = 0.5,
      matrix.color = "grey20", 
      main.bar.color = "grey20",
      sets.bar.color = "gray",
      set_size.show = F,
      query.legend="top",
      queries=list(list(query=elements, params=list("Race", "White"), active=T, color="red", query.name="White"),
                   list(query=elements, params=list("Race", "Other/Unknown"), active=T, color="blue", query.name="Other/Unknown"),
                   list(query=elements, params=list("Race", "Black"), active=T, color="yellow", query.name="Black"),
                   list(query=elements, params=list("Race", "Asian"), active=T, color="green", query.name="Asian")
      )
)

#ComplexHeatmap
m = make_comb_mat(pheno, top_n_sets = 6)
HeatmapAnnotation("Intersection\nsize" = anno_barplot(comb_size(m), 
                                                      border = FALSE, gp = gpar(fill = "black"), height = unit(3, "cm")), 
                  annotation_name_side = "left", annotation_name_rot = 0)
UpSet(m)

disease <- c("AF", "BC", "CAD", "CC", "PC", "T2D")
race <- as.factor(pheno$Race)

#Create a matrix with every race
m_list = tapply(seq_len(nrow(pheno)), race, function(ind) {
  m = make_comb_mat(pheno[ind, disease, drop = FALSE])
  m[comb_degree(m) > 0]
})


#-------------Collect all set size to make a table for creating a set barchart -----------
set_race_table <- data.frame()
for(i in length(m_list):1) {set_race_table <- rbind(set_race_table, set_size(m_list[[i]]))}
names(set_race_table) <- c("AF", "BC", "CAD", "CC", "PC", "T2D")
set_race_table <- set_race_table[,c(3,6,1,5,2,4)]

set_race_matrix <- as.matrix(set_race_table)

races=c("White", "Other/Unknown","Black", "Asian")
colors = c("lightblue","brown","orange", "green")
diseases = c("CAD", "AF", "PC", "CC", "T2D", "BC")


barplot(set_race_matrix, horiz=T, 
        names.arg=diseases, 
        xlab="Set size", 
        col=colors,
        space=0.8,
        width=0.8
)

# Add the legend to the chart
#legend("topright", races, cex = 1.3, fill = colors)


#----------Collect all Combination/Intersection to make a table for creating a stacked barchart -----------
white <- data.frame(name=comb_name(m_list[[4]]), White=comb_size(m_list[[4]]))
white_otherUnknow <- merge(white, data.frame(name=comb_name(m_list[[3]]), OtherUnknown=comb_size(m_list[[3]])), all=T)
white_otherUnknown_black <- merge(white_otherUnknow, data.frame(name=comb_name(m_list[[2]]), Black=comb_size(m_list[[2]])), all=T)
white_otherUnknown_black_asian <- merge(white_otherUnknown_black, data.frame(name=comb_name(m_list[[1]]), Asian=comb_size(m_list[[1]])),all=T)

#Replace "NA" with 0
white_otherUnknown_black_asian[is.na(white_otherUnknown_black_asian)] <- 0
#comb_race_table.ordered <- white_otherUnknown_black_asian[order(-white_otherUnknown_black_asian[,2]),]
intersection_num <- c()
for (i in 1:nrow(white_otherUnknown_black_asian)){intersection_num <- c(intersection_num, sum(white_otherUnknown_black_asian[i,2:5]))}
white_otherUnknown_black_asian$intersection_num <- intersection_num


comb_race_table.ordered <- white_otherUnknown_black_asian[order(-white_otherUnknown_black_asian$intersection_num),]
comb_race_table.ordered <- comb_race_table.ordered[,c(1:5)]

#Transpose the dataframe 
col.names <- comb_race_table.ordered$name
t.tmp <- as.data.frame(t(comb_race_table.ordered))
t1.tmp <- t.tmp[2:5,]
names(t1.tmp) <- col.names

#Stacked bar plot of the intersection of allcombinations
comb_race_matrix <- as.matrix(t1.tmp)
colors = c("lightblue","brown","orange", "green")
races=c("White", "Other/Unknown","Black", "Asian")

barplot(comb_race_matrix, horiz=F, names.arg=col.names, 
        xlab="Intersection size",
        ylim=c(0,2100),
        col=colors,
        space=0.8,
        width=0.8
)
# Add the legend to the chart
legend("topright", races, cex = 1.3, fill = colors)

#Add number on each bar
intersection_num.sorted <- sort(intersection_num, decreasing=T)
for (i in 1:length(intersection_num.sorted)){
  posx <- i+0.425*i
  posy <- intersection_num.sorted[i] + 80
  text (posx, posy, intersection_num.sorted[i], srt=45)
}




