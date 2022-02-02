#library(eulerr)
#library(upsetjs)
library(UpSetR)
library(ComplexHeatmap)
library(tidyverse)

setwd("~/Documents/PRS_Jason/Cutoff_adjustment/PCA_adjustment_MDS/high_risk_analysis")

#Read in high risk data (> cutOff.adjusted)
af <- read.csv("above_cutoff_af.csv", header=T, stringsAsFactor=F)
bc <- read.csv("above_cutoff_bc.csv", header=T, stringsAsFactor=F)
cad <- read.csv("above_cutoff_cad.csv", header=T, stringsAsFactor=F)
cc <- read.csv("above_cutoff_cc.csv", header=T, stringsAsFactor=F)
pc <- read.csv("above_cutoff_pc.csv", header=T, stringsAsFactor=F)
t2d <- read.csv("above_cutoff_t2d.csv", header=T, stringsAsFactor=F)

#Create a table with phenotype as "1" for cases
af_hr <- af[, c(1,2,7)]
af_hr$AF <- 1

bc_hr <- bc[, c(1,2,7)]
bc_hr$BC<- 1

cad_hr <- cad[, c(1,2,7)]
cad_hr$CAD <- 1

cc_hr <- cc[, c(1,2,7)]
cc_hr$CC <- 1

pc_hr <- pc[, c(1,2,7)]
pc_hr$PC <- 1

t2d_hr <- t2d[, c(1,2,7)]
t2d_hr$T2D <- 1


#Play around with upset.js

create_id_list <- function(df){
  #Input is a dataframe with fid, iid, race, phenotype
  #Output a list of the iid
  iid <- c()
  for (i in 1:nrow(df)){iid <- c(iid, df[i,2])}
  
  return(iid)
}

af_iid <- create_id_list(af_hr)
bc_iid <- create_id_list(bc_hr)
cad_iid <- create_id_list(cad_hr)
cc_iid <- create_id_list(cc_hr)
pc_iid <- create_id_list(pc_hr)
t2d_iid <- create_id_list(t2d_hr)

listInput <- list(af=af_iid, bc=bc_iid, cad=cad_iid, cc=cc_iid, pc=pc_iid, t2d=t2d_iid)
render <- function(upsetjs) {
  upsetjs %>% fromList(listInput) %>% chartTheme(selection.color="", has.selection.opacity=0.3) %>% interactiveChart()
}

v <- upsetjs() %>% render()
v

#------------------------------------

#Combine six high_risk data frame into one
af_bc <- merge(af_hr, bc_hr, all=T)
af_bc_cad <- merge(af_bc, cad_hr, all=T)
af_bc_cad_cc <- merge(af_bc_cad, cc_hr, all=T)
af_bc_cad_cc_pc <- merge(af_bc_cad_cc, pc_hr, all=T)
all_hr <- merge(af_bc_cad_cc_pc, t2d_hr, all=T)

#Replace NAs with 0
all_hr[is.na(all_hr)] <-0

#Change the races Other or Unknown into "Other/Unknown"
all_hr[which(all_hr$Race == "Unknown" | all_hr$Race == "Other"),3] <- "Other/Unknown"

#Upset plot
upset(all_hr, sets=c("AF", "BC", "CAD", "CC", "PC", "T2D"), 
      order.by = "freq", 
      decreasing = T, 
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







upset(all_hr, sets=c("AF", "BC", "CAD", "CC", "PC", "T2D"), query.legend="top",
      queries=list(list(query=elements, params=list("Race", "White"), active=T, color="red", query.name="White"),
                   
                   list(query=elements, params=list("Race", "Other", "Unknown"), active=T, color="blue", query.name="Other/Unknown")
      ))


#ComplexHeatmap
m = make_comb_mat(all_hr, top_n_sets = 6)
HeatmapAnnotation("Intersection\nsize" = anno_barplot(comb_size(m), 
                                                      border = FALSE, gp = gpar(fill = "black"), height = unit(3, "cm")), 
                  annotation_name_side = "left", annotation_name_rot = 0)
UpSet(m)

#Change right barplot to left
ss = set_size(m)
cs = comb_size(m)
ht = UpSet(m, 
           set_order = order(ss),
           comb_order = order(comb_degree(m), -cs),
           top_annotation = HeatmapAnnotation(
             "Intersection Size" = anno_barplot(cs, 
                                                  ylim = c(0, max(cs)*1.1),
                                                  border = FALSE, 
                                                  gp = gpar(fill = "black"), 
                                                  height = unit(4, "cm")
             ), 
             annotation_name_side = "left", 
             annotation_name_rot = 90),
           left_annotation = rowAnnotation(
             "Total Number" = anno_barplot(-ss, 
                                               baseline = 0,
                                               axis_param = list(
                                                 at = c(0, -500, -1500, -2500, -3500),
                                                 labels = c(0, 500, 1500, 2500, 3500),
                                                 labels_rot = 0),
                                               border = FALSE, 
                                               gp = gpar(fill = "black"), 
                                               width = unit(3, "cm")
             ),
             set_name = anno_text(set_name(m), 
                                  location = 0.5, 
                                  just = "center",
                                  width = max_text_width(set_name(m)) + unit(4, "mm"))
           ), 
           right_annotation = NULL,
           show_row_names = FALSE)
ht = draw(ht)
od = column_order(ht)
decorate_annotation("Intersection Size", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = c("left", "bottom"), 
            gp = gpar(fontsize = 6, col = "#404040"), rot = 45)
})

disease <- c("AF", "BC", "CAD", "CC", "PC", "T2D")
race <- as.factor(all_hr$Race)

m_list = tapply(seq_len(nrow(all_hr)), race, function(ind) {
  m = make_comb_mat(all_hr[ind, disease, drop = FALSE])
  m[comb_degree(m) > 0]
})


#To compare between multiple groups with UpSet plots, we need to normalize 
#all the matrices to make them have same sets and same combination sets. 
#normalize_comb_mat() basically adds zero to the new combination sets which 
#were not there before.
m_list = normalize_comb_mat(m_list)
sapply(m_list, comb_size)

max_set_size = max(sapply(m_list, set_size))
max_comb_size = max(sapply(m_list, comb_size))

#Vetically plot four plots for each race. "%v%" adds heatmaps or annotations to a heatmap list vertically.

#ss = set_size
#cs = comb_size
ht_list = NULL
for(i in seq_along(m_list)) {
  ht_list = ht_list %v%
    UpSet(m_list[[i]], 
          #set_order = order(ss),
          #comb_order = order(comb_degree(m), -cs),
          row_title = paste0(names(m_list)[i]),
          set_order = NULL, comb_order = NULL,
          top_annotation = upset_top_annotation(m_list[[i]], ylim = c(0, max_comb_size)),
          right_annotation = upset_right_annotation(m_list[[i]], ylim = c(0, max_set_size)))
}
ht_list

top_ha <- HeatmapAnnotation(
  "Asian" = anno_barplot(comb_size(m_list[[1]]), ylim = c(0, max(comb_size(m_list[[1]]))*1.3), gp=gpar(fill="black"), height=unit(2.5,"cm")),
  "Black" = anno_barplot(comb_size(m_list[[2]]), ylim = c(0, max(comb_size(m_list[[2]]))*1.3), gp=gpar(fill="black"), height=unit(2.5,"cm")),
  "Other/\nUnknown" = anno_barplot(comb_size(m_list[[3]]), ylim = c(0, max(comb_size(m_list[[3]]))*1.3), gp=gpar(fill="black"), height=unit(2.5,"cm")),
  "White" = anno_barplot(comb_size(m_list[[4]]), ylim = c(0, max(comb_size(m_list[[4]]))*1.3), gp=gpar(fill="black"), height=unit(2.5,"cm")),
  gap = unit(2, "mm"), annotation_name_side = "left", annotation_name_rot = 0
)

ss=set_size(m_list[[4]])
cs=comb_size(m_list[[4]])
ht=draw(UpSet(m_list[[4]], 
              set_order=order(ss), 
              comb_order=order(comb_degree(m_list[[4]]), -cs),
              top_annotation=top_ha))

od = column_order(ht)
cs = comb_size(m_list[[1]])
decorate_annotation("Asian", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = c("center", "bottom"), 
            gp = gpar(fontsize = 8, col = "#404040"), rot = 0)
})


cs = comb_size(m_list[[2]])
decorate_annotation("Black", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = c("center", "bottom"), 
            gp = gpar(fontsize = 8, col = "#404040"), rot = 0)
})

cs = comb_size(m_list[[3]])
decorate_annotation("Other/\nUnknown", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = c("left", "bottom"), 
            gp = gpar(fontsize = 8, col = "#404040"), rot = 45)
})

cs = comb_size(m_list[[4]])
decorate_annotation("White", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = c("left", "bottom"), 
            gp = gpar(fontsize = 8, col = "#404040"), rot = 45)
})


#-------------Collect all set size to make a table for creating a set barchart -----------
set_race_table <- data.frame()
for(i in length(m_list):1) {set_race_table <- rbind(set_race_table, set_size(m_list[[i]]))}
names(set_race_table) <- c("AF", "BC", "CAD", "CC", "PC", "T2D")
set_race_table <- set_race_table[,c(3,6,1,5,4,2)]
print(set_race_table)
set_race_matrix <- as.matrix(set_race_table)

races=c("White", "Other/Unknown","Black", "Asian")
colors = c("lightblue","brown","orange", "green")
diseases = c("CAD", "T2D", "AF", "PC", "CC", "BC")


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
        ylim=c(0,3000),
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
  posy <- intersection_num.sorted[i] + 120
  text (posx, posy, intersection_num.sorted[i], srt=45)
}


#text(2,2800, "2637", srt=45)
#text(3,2060, "1900", srt=45)

