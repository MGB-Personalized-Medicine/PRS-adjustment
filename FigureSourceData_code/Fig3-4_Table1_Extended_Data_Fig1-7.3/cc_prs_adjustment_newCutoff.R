setwd("~/Documents/PRS_Jason/prs_analysis/1kg_imputed_normalized/colorectal")

#------------------- Read in data -----------------------

mega_profile <- read.table("mega_colorectal.profile", header=F, stringsAsFactors=F)
megaex_profile <- read.table("megaex_colorectal.profile", header=F, stringsAsFactors=F)
mega1a_profile <- read.table("meg_a1_a_colorectal.profile", header=F, stringsAsFactors=F)
mega1b_profile <- read.table("meg_a1_b_colorectal.profile", header=F, stringsAsFactors=F)
megc_profile <- read.table("meg_c_colorectal.profile", header=F, stringsAsFactors=F)
megd_profile <- read.table("meg_d_colorectal.profile", header=F, stringsAsFactors=F)
mege_profile <- read.table("meg_e_colorectal.profile", header=F, stringsAsFactors=F)
megx1_profile <- read.table("meg_x1_colorectal.profile", header=F, stringsAsFactors=F)

head(mega_profile)
head(megaex_profile)
head(mega1a_profile)
head(mega1b_profile)
head(megc_profile)
head(megd_profile)
head(mege_profile)
head(megx1_profile)



mega_score <- mega_profile[,6]
megaex_score <- megaex_profile[,6]
mega1a_score <- mega1a_profile[,6]
mega1b_score <- mega1b_profile[,6]
megc_score <- megc_profile[,6]
megd_score <- megd_profile[,6]
mege_score <- mege_profile[,6]
megx1_score <- megx1_profile[,6]

all_batch_profile <- read.table("all_batch.profile", header=T, stringsAsFactors=F)
all_batch_score <- all_batch_profile[,6]

avg_all <- mean(all_batch_score)
sd_all <- sd(all_batch_score)
avg_all
sd_all

all_batch_std_score <- (all_batch_score - avg_all) / sd_all
std_avg_all <- mean(all_batch_std_score)
std_sd_all <- sd(all_batch_std_score)

std_avg_all
std_sd_all

min_score <- min(all_batch_score)
max_score <- max(all_batch_score)
min_score
max_score


#------------------- Top 0.5% subjects. This is for others' request ------------------------------
setwd("~/Documents/PRS_Jason/prs_analysis/1kg_imputed_normalized/colorectal")
all_batch_profile <- read.table("all_batch.profile", header=T, stringsAsFactors=F)

output <- "~/Documents/PRS_Jason/prs_analysis/1kg_imputed_normalized/alldiseases_prs_top0.5/colorectalCA_top0.5.csv"

subject_num <- nrow(all_batch_profile)
top0.5_num <- round((subject_num * 0.5 / 100),0)

all_batch_profile.ordered <- all_batch_profile[order(all_batch_profile[,6], decreasing=T), ]
top0.5 <- all_batch_profile.ordered[1:top0.5_num, c(1,2,6)]
write.csv(top0.5, output, row.names=F, quote=F)




#------------------- Density Analysis across Batches -----------------------

#plot(density(mega_score, adjust=2), col="black", xlim=c(4, 8), ylim=c(0, 0.9), xlab="PRS score", lty="dotted", main="Density ~ PRS score of Colorectal Cancer")
par(mar=c(3.5,3.5,1,1))
plot(density(mega_score, adjust=2), col="black", xlim=c(4.2, 8.4), ylim=c(0, 1), xlab="", ylab="", lty="dotted", main="")
title(ylab="Density", line=2, cex.lab=1.0)
title(xlab="PRS score", line=2, cex.lab=1.0)

lines(density(megaex_score, adjust=2), lty="dotted", col="red")
lines(density(mega1a_score, adjust=2), lty="dotted", col="orange")
lines(density(mega1b_score, adjust=2), lty="dotted", col="yellow")
lines(density(megc_score, adjust=2), lty="dotted", col="green")
lines(density(megd_score, adjust=2), lty="dotted", col="cyan")
lines(density(mege_score, adjust=2), lty="dotted", col="blue")
lines(density(megx1_score, adjust=2), lty="dotted", col="purple")

lines(density(all_batch_score, adjust=2), col="darkblue")

#text(x=5e-05, y=15000, labels = c("black: MEGA", "red: MEGAEX", "orange: MEG_A1_A"))
Vars <- c("MEGA", "MEGAEX", "MEG_A1_A", "MEG_A1_B", "MEG_C", "MEG_D", "MEG_E", "MEG_X1", "ALL")  # one per row
colors <- c("black", "red", "orange", "yellow", "green", "cyan", "blue", "purple", "darkblue")  # one color for each row; should be same length as Vars
ltys <- 1  # one lty for each column

nr <- length(Vars) 
nc <- length(ltys)
legend("topright", rep(Vars, nc), col = colors, lty = rep(ltys, each = nr),
       ncol = nc, cex = 0.44)


#One-way ANOVA analysis
anova_1way <- aov(PRS~Gender, data=score_pheno)
TukeyHSD(anova_1way)   #Pairwise comparison

#--------------------- Density Analysis across Gender and Age and Race ----------------------

#README file contains how prs_score_colorectalC_covariates_chip.tsv is created
score_pheno <- read.table("prs_score_colorectalC_covariates_chip.tsv", header=T, stringsAsFactors=F)
all_score <- score_pheno[,3]

#Density analysis across Gender
male <- score_pheno[score_pheno[,5]=="M" ,]
female <- score_pheno[score_pheno[,5]=="F" ,]
male_score <- male[,3]
female_score <- female[,3]

#Density plot
par(mar=c(3.5,3.5,1,1))
plot(density(male_score, adjust=2), col="red", xlim=c(4.2, 8.2), ylim=c(0, 1.0), xlab="",ylab="",lty="solid", main="")
title(ylab="Density", line=2, cex.lab=1.0)
title(xlab="PRS score", line=2, cex.lab=1.0)

lines(density(female_score, adjust=2), lty="solid", col="blue")
lines(density(all_score, adjust=2), lty="dotted", col="black")


#text(x=5e-05, y=15000, labels = c("black: MEGA", "red: MEGAEX", "orange: MEG_A1_A"))
Vars <- c("Male", "Female", "ALL")  # one per row
colors <- c("red", "blue", "black")  # one color for each row; should be same length as Vars
ltys <- 1  # one lty for each column

nr <- length(Vars) 
nc <- length(ltys)
legend("topright", rep(Vars, nc), col = colors, lty = rep(ltys, each = nr),
       ncol = nc, cex = 0.6)


#One-way ANOVA analysis
anova_1way <- aov(PRS~Gender, data=score_pheno)
TukeyHSD(anova_1way)   #Pairwise comparison


#Density analysis across Age
par(mar=c(3.5,3.5,1,1))
hist(score_pheno[,6], breaks=10, main="Histogram of Age", xlab="", ylab="")
title(xlab="Age", line=2, cex.lab=1.0)
title(ylab="Frequency", line=2, cex.lab=1.0)

age10_score <- score_pheno[score_pheno[,6]<=10, 3]   #Only 1 data"
age20_score <- score_pheno[score_pheno[,6]>10 & score_pheno[,6]<=20, 3]
age30_score <- score_pheno[score_pheno[,6]>20 & score_pheno[,6]<=30, 3]
age40_score <- score_pheno[score_pheno[,6]>30 & score_pheno[,6]<=40, 3]
age50_score <- score_pheno[score_pheno[,6]>40 & score_pheno[,6]<=50, 3]
age60_score <- score_pheno[score_pheno[,6]>50 & score_pheno[,6]<=60, 3]
age70_score <- score_pheno[score_pheno[,6]>60 & score_pheno[,6]<=70, 3]
age80_score <- score_pheno[score_pheno[,6]>70 & score_pheno[,6]<=80, 3]
age90_score <- score_pheno[score_pheno[,6]>80 & score_pheno[,6]<=90, 3]
age100_score <- score_pheno[score_pheno[,6]>90 & score_pheno[,6]<=100, 3]
age110_score <- score_pheno[score_pheno[,6]>100 & score_pheno[,6]<=110, 3]   #Only 5 data

#Density plot
par(mar=c(3.5,3.5,1,1))
plot(density(age30_score, adjust=2), col="black", xlim=c(4.2, 8.2), ylim=c(0, 1.0), xlab="", ylab="", lty="dotted", main="")
title(xlab="PRS score", line=2, cex.lab=1.0)
title(ylab="Frequency", line=2, cex.lab=1.0)

#lines(density(age30_score, adjust=2), lty="dotted", col="red")
lines(density(age40_score, adjust=2), lty="dotted", col="green")
lines(density(age50_score, adjust=2), lty="dotted", col="blue")
lines(density(age60_score, adjust=2), lty="solid", col="red")
lines(density(age70_score, adjust=2), lty="solid", col="yellow")
lines(density(age80_score, adjust=2), lty="solid", col="blue")
lines(density(age90_score, adjust=2), lty="solid", col="cyan")
lines(density(age100_score, adjust=2), lty="solid", col="purple")
#lines(density(age110_score, adjust=2), lty="solid", col="black")  #Only five data

#text(x=5e-05, y=15000, labels = c("black: MEGA", "red: MEGAEX", "orange: MEG_A1_A"))
Vars <- c("Age20_30","Age30_40", "Age40_50", "Age50_60", "Age60_70", "Age70_80", "Age80_90", "Age90_100")  # one per row
colors <- c("black", "green", "blue", "red", "yellow", "blue", "red", "purple")  # one color for each row; should be same length as Vars
ltys <- c("dotted", "dotted","dotted","solid","solid","solid","solid","solid")  # one lty for each column

nr <- length(Vars) 
legend("topright", rep(Vars, 1), col = colors, lty = ltys, ncol = 1, cex = 0.45)


#One-way ANOVA analysis
age_score <- data.frame(age=c(rep("20",length(age20_score)),rep("30",length(age30_score)),rep("40",length(age40_score)),rep("50",length(age50_score)),rep("60",length(age60_score)),rep("70",length(age70_score)),rep("80",length(age80_score)),rep("90",length(age90_score)),rep("100",length(age100_score))), 
                        score=c(age20_score, age30_score, age40_score,age50_score, age60_score,age70_score, age80_score,age90_score, age100_score))
anova_1way <- aov(score~age, data=age_score)
TukeyHSD(anova_1way)   #Pairwise comparison


#+++++++++++++ Density analysis across RACE ++++++++++++++++
#README file contains how prs_score_colorectalC_covariates_chip.tsv is created
score_pheno <- read.table("prs_score_colorectalC_covariates_chip.tsv", header=T, stringsAsFactors=F)


white_score <- score_pheno[score_pheno[,7]=="White" ,3]
black_score <- score_pheno[score_pheno[,7]=="Black" ,3]
asian_score <- score_pheno[score_pheno[,7]=="Asian" ,3]
other_score <- score_pheno[score_pheno[,7]=="Other" ,3]
unknown_score <- score_pheno[score_pheno[,7]=="Unknown" ,3]

#Barplot of races
race_table <- table(score_pheno[,7])
barplot(race_table, ylab="Number of subjects", main="Barplot of Races")

#Deansity plot
plot(density(white_score, adjust=2), col="black", xlim=c(4, 8), ylim=c(0, 1), xlab="PRS score", lty="dotted", main="Density ~ PRS score of ColorectalAD")

#For publication
par(mar=c(3.5,3.5,1,1))
plot(density(white_score, adjust=2), col="black", xlim=c(4.2, 8.2), ylim=c(0, 1), xlab="",ylab="",lty="dotted", main="")
title(xlab="PRS score", line=2, cex.lab=1.0)
title(ylab="Frequency", line=2, cex.lab=1.0)

lines(density(black_score, adjust=2), lty="dotted", col="purple")
lines(density(asian_score, adjust=2), lty="dotted", col="green")
lines(density(other_score, adjust=2), lty="dotted", col="red")
lines(density(unknown_score, adjust=2), lty="dotted", col="blue")

Vars <- c("White", "Black", "Asian", "Other", "Unknown")  # one per row
colors <- c("black", "purple", "green", "red", "blue")  # one color for each row; should be same length as Vars
ltys <- c("dotted","dotted","dotted","dotted","dotted")  # one lty for each column

nr <- length(Vars) 
legend("topright", rep(Vars, 1), col = colors, lty = ltys, ncol = 1, cex = 0.5)

#One-way ANOVA analysis
anova_1way <- aov(PRS~Race, data=score_pheno)
TukeyHSD(anova_1way)   #Pairwise comparison




#--------------------------- Phenotype analysis using RAW PRS score -----------------------------
#score_pheno is all batch profile, and all_score is PRS score for all subjects
quan=0.05
score_quantile <- quantile(all_score, probs=seq(0,1,quan), names=T)  #11 numbers with %0=samllest_number and %100=largest+number


quan_cases.df <- data.frame()
for (i in (2:length(score_quantile))){
  if (i==2){
    rows <- which(all_score <= score_quantile[i])
    quantile_midpoint <- score_quantile[i-1] + (score_quantile[i] - score_quantile[i-1]) / 2
  }
  else{
    rows <- which(all_score > score_quantile[i-1] & all_score <= score_quantile[i])
    quantile_midpoint <- score_quantile[i-1] + (score_quantile[i] - score_quantile[i-1]) / 2
  }
  quantile_pheno <- score_pheno[rows, ]
  subject_num <- nrow(quantile_pheno)
  cases <- sum(quantile_pheno[, 4]==1)
  healthy <- sum(quantile_pheno[, 4]==0)
  odds <- cases/healthy
  
  quan_cases.df <- rbind(quan_cases.df, c(quantile_midpoint, cases, subject_num, odds))
}

#Plot density, scatter, fit line
#Note: the case number times 20 to make it scaled with the density curve
prs_score <- quan_cases.df[,1]
odds <- quan_cases.df[,4]
adj <- 20     #scale to match density
odds.adj <- odds * adj

plot(prs_score, odds.adj)

plot(prs_score, odds.adj, xlab="PRS score", ylab=paste0("Odds (x", adj, ")"), xlim=c(4, 8), ylim=c(0,1), main="Colorectal Tumor PRS Score Analysis")
lines(density(all_score, adjust=2), col="darkblue")
abline(lm(odds.adj~prs_score), col="red")
corr <- cor(quan_cases.df[,1], quan_cases.df[,4])
text(7, 0.9, paste0("r=", round(corr,4)))


#----------------------- Phenotype Analysis using Standardized PRS score ------------------------
#score_pheno is all batch profile, and all_score is PRS score for all subjects
score_pheno <- read.table("prs_score_colorectalC_covariates_chip.tsv", header=T, stringsAsFactors=F)
all_score <- score_pheno[,3]
std_score=(all_score - mean(all_score))/sd(all_score)  #calculate standardized PRS scores
plot(density(std_score, adjust=2))
score_pheno$std_score <- std_score   #Add standardized PRS into all_prs table

score_quantile <- quantile(std_score, probs=seq(0,1,0.01), names=T)  #11 numbers with %0=samllest_number and %100=largest+number


quan_cases.df <- data.frame()
for (i in (2:length(score_quantile))){
  if (i==2){
    rows <- which(std_score <= score_quantile[i])
    quantile_midpoint <- score_quantile[i-1] + (score_quantile[i] - score_quantile[i-1]) / 2
  }
  else{
    rows <- which(std_score > score_quantile[i-1] & std_score <= score_quantile[i])
    quantile_midpoint <- score_quantile[i-1] + (score_quantile[i] - score_quantile[i-1]) / 2
  }
  quantile_pheno <- score_pheno[rows, ]
  subject_num <- nrow(quantile_pheno)
  cases <- sum(quantile_pheno[, 4]==1)
  healthy <- sum(quantile_pheno[, 4]==0)
  odds <- cases/healthy
  
  quan_cases.df <- rbind(quan_cases.df, c(quantile_midpoint, cases, subject_num, odds))
}

#Plot density, scatter, fit line
#Note: the case number times 20 to make it scaled with the density curve
std_score <- quan_cases.df[,1]
odds <- quan_cases.df[,4]
adj <- 6     #scale to match density
odds.adj <- odds * adj

plot(std_score, odds.adj)

plot(std_score, odds.adj, xlab="Standardized PRS score", ylab=paste0("Odds (x", adj, ")"), xlim=c(-3, 3), ylim=c(0,0.4), main="Colorectal Tumor PRS Score Analysis")
lines(density(std_score, adjust=2), col="darkblue")
abline(lm(odds.adj~std_score), col="red")
corr <- cor(quan_cases.df[,1], quan_cases.df[,4])
text(2.5, 0.35, paste0("r=", round(corr,4)))




#----------------------- Logistic Regression ------------------------------

score_pheno <- read.table("prs_score_colorectalC_covariates_chip.tsv", header=T, stringsAsFactors=F)
all_score <- score_pheno[,3]
std_score <- (all_score - mean(all_score))/sd(all_score)  #calculate standardized PRS scores
score_pheno$std_score <- std_score   #Add standardized PRS into all_prs table

mylogit <- glm(COLCA~std_score+Gender+Age+Race+chip, data=score_pheno, family=binomial)  #Build model
summary(mylogit)   #View the model
confint(mylogit)         #confident interval of predictors


#Obtain percentage above cutoff
std_score_coef <- coef(mylogit)[2]   #Get coefficient of std_score
attributes(std_score_coef) <- NULL   #Delete attribute "std_score"
OR_change_per_prs_sd <- exp(std_score_coef)
OR_change_per_prs_sd
std_threshold <- log(2.5) / log(OR_change_per_prs_sd)   #2.5 is target raw PRS score cutoff
std_threshold

#Convert std_threshold into raw score
avg <- mean(all_score)
stddev <- sd(all_score)
raw_threshold <- (std_threshold * stddev) + avg
raw_threshold

#Calculate percentage of patients in population
above_threshold <- std_score[which(std_score > std_threshold)]
percentage_above_threshold <- (length(above_threshold) / length(std_score)) * 100
percentage_above_threshold

exp(cbind(OR = coef(mylogit), confint(mylogit)))   #odds ratio and ci of all predictors

new_data <- score_pheno
pred <- predict(mylogit, newdata=new_data)
pred_odds <- exp(pred)
#plot(std_score, pred_odds)

#Plot density and predicted odds 
adj <- 5     #scale to match density
pred_odds.adj <- pred_odds * adj

plot(std_score, pred_odds*adj, xlab="Standardized PRS score", ylab=paste0("predicted odds (x", adj, ")"), xlim=c(-4.5, 4.5), ylim=c(0,0.4), main="Colorectal Tumor PRS Score Analysis")
lines(density(std_score, adjust=2), col="darkblue")
#lines(std_score, pred_odds*adj, col='red')
lo <- loess(pred_odds*adj~std_score)$fitted
lines(std_score, lo, col='red')
#abline(lm(pred_odds.adj~standardized_score), col="green")
#corr <- cor(quan_cases.df[,1], quan_cases.df[,4])
#text(2.5, 0.35, paste0("r=", corr))


#Plot density and predicted odds with quantile
quan=0.01
std_prs_quantile <- quantile(std_score, probs=seq(0,1,quan), names=T)  #11 numbers with %0=samllest_number and %100=largest+number

quan_cases.df <- data.frame()
for (i in (2:length(std_prs_quantile))){
  if (i==2){
    rows <- which(std_score <= std_prs_quantile[i])
    quantile_midpoint <- std_prs_quantile[i-1] + (std_prs_quantile[i] - std_prs_quantile[i-1]) / 2
  }
  else{
    rows <- which(std_score > std_prs_quantile[i-1] & std_score <= std_prs_quantile[i])
    quantile_midpoint <- std_prs_quantile[i-1] + (std_prs_quantile[i] - std_prs_quantile[i-1]) / 2
  }
  #quantile_prs <- standardized_score[rows]
  quantile_pred_odds <- mean(pred_odds[rows])
  
  quan_cases.df <- rbind(quan_cases.df, c(quantile_midpoint, quantile_pred_odds))
}

#Plot density and predicted odds 
adj <- 10     #scale to match density
pred_odds.adj <- pred_odds * adj

plot(quan_cases.df[,1], quan_cases.df[,2]*adj, xlab="Standardized PRS score", ylab=paste0("Pred Odds (x", adj, ")"), xlim=c(-4.5, 4.5), ylim=c(0,0.4), main="Colorectal Tumor PRS Score Analysis")

par(mar=c(3.5,3.5,1,1))
plot(quan_cases.df[,1], quan_cases.df[,2]*adj, xlab="", ylab="", xlim=c(-4.5, 4.5), ylim=c(0,0.4), main="")
title(xlab="PRS score", line=2, cex.lab=1.0)
title(ylab=paste0("Pred Odds (x", adj, ")"), line=2, cex.lab=1.0)

lines(density(std_score, adjust=2), col="darkblue")
#lines(quan_cases.df[,1], quan_cases.df[,2]*adj, col='red')
#abline(lm(pred_odds.adj~standardized_score), col="red")
#corr <- cor(quan_cases.df[,1], quan_cases.df[,2])
#text(0.038, 300, paste0("r=", corr))


#Calculate Odds Ratio (OR)
i <- length(std_prs_quantile) #total number of quantiles
rows <- which(std_score > std_prs_quantile[i-1])
#print(head(rows))  
top_quantile_prs <- std_score[rows]
other_quantile_prs <- std_score[-rows]

top_quantile_odd <- mean(pred_odds[rows])
other_quantile_odd <- mean(pred_odds[-rows])

or <- top_quantile_odd/other_quantile_odd
or




#---------------------- Race Analysis ------------------

#README file contains how prs_score_prostate_covariates_chip.tsv is created
score_pheno <- read.table("prs_score_colorectalC_covariates_chip.tsv", header=T, stringsAsFactors=F)
all_score <- score_pheno[,3]

white_score <- score_pheno[score_pheno[,7]=="White" ,3]
black_score <- score_pheno[score_pheno[,7]=="Black" ,3]
asian_score <- score_pheno[score_pheno[,7]=="Asian" ,3]
other_score <- score_pheno[score_pheno[,7]=="Other" ,3]
unknown_score <- score_pheno[score_pheno[,7]=="Unknown" ,3]

mean(white_score)
sd(white_score)

#Barplot of races
race_table <- table(score_pheno[,7])
barplot(race_table, ylab="Number of subjects", main="Barplot of Race")

#Deansity plot across Race
plot(density(white_score, adjust=2), col="black", xlim=c(4, 8), ylim=c(0, 0.9), xlab="PRS score", lty="dotted", main="Density ~ PRS score of ColorectalCA")
lines(density(black_score, adjust=2), lty="dotted", col="purple")
lines(density(asian_score, adjust=2), lty="dotted", col="green")
lines(density(other_score, adjust=2), lty="dotted", col="red")
lines(density(unknown_score, adjust=2), lty="dotted", col="blue")

Vars <- c("White", "Black", "Asian", "Other", "Unknown")  # one per row
colors <- c("black", "purple", "green", "red", "blue")  # one color for each row; should be same length as Vars
ltys <- c("dotted","dotted","dotted","dotted","dotted")  # one lty for each column

nr <- length(Vars) 
legend("topright", rep(Vars, 1), col = colors, lty = ltys, ncol = 1, cex = 0.8)

#One-way ANOVA analysis
anova_1way <- aov(PRS~Race, data=score_pheno)
TukeyHSD(anova_1way)   #Pairwise comparison



#+++++++++++++++++++++++++ Analysis for each RACE using raw score +++++++++++++++++++++++
score_pheno <- read.table("prs_score_colorectalC_covariates_chip.tsv", header=T, stringsAsFactors=F)
white <- score_pheno[score_pheno[,7]=="White" ,]
white_score <- white[,3]

quan=0.02
score_quantile <- quantile(white_score, probs=seq(0, 1, quan), names=T)  #11 numbers with %0=samllest_number and %100=largest+number


quan_cases.df <- data.frame()
for (i in (2:length(score_quantile))){
  if (i==2){
    rows <- which(white_score <= score_quantile[i])
    quantile_midpoint <- score_quantile[i-1] + (score_quantile[i] - score_quantile[i-1]) / 2
  }
  else{
    rows <- which(white_score > score_quantile[i-1] & white_score <= score_quantile[i])
    quantile_midpoint <- score_quantile[i-1] + (score_quantile[i] - score_quantile[i-1]) / 2
  }
  quantile_pheno <- white[rows, ]
  subject_num <- nrow(quantile_pheno)
  cases <- sum(quantile_pheno[, 4]==1)
  healthy <- sum(quantile_pheno[, 4]==0)
  odds <- cases/healthy
  
  quan_cases.df <- rbind(quan_cases.df, c(quantile_midpoint, cases, subject_num, odds))
}

#Plot density, scatter, fit line
#Note: the case number times 20 to make it scaled with the density curve
prs_score <- quan_cases.df[,1]
odds <- quan_cases.df[,4]
adj <- 12     #scale to match density
odds.adj <- odds * adj

plot(prs_score, odds.adj)

plot(prs_score, odds.adj, xlab="PRS score", ylab=paste0("Odds (x", adj, ")"), xlim=c(4, 8), ylim=c(0, 0.9), main="White ColorectalCA PRS Score Analysis")

par(mar=c(3.5,3.5,1,1))
plot(prs_score, odds.adj, xlab="", ylab="", xlim=c(4.2, 8.2), ylim=c(0, 1.0), main="")
title(xlab="PRS score", line=2, cex.lab=1.0)
title(ylab=paste0("Odds (x", adj, ")"), line=2, cex.lab=1.0)


lines(density(white_score, adjust=2), lty="dotted", col="black")
abline(lm(odds.adj~prs_score), col="red")
shapiro.test(quan_cases.df[,1])   #p value < 0.05, normal distribution
shapiro.test(quan_cases.df[,4])   #p value > 0.05, non-normal distribution
#corr <- cor(quan_cases.df[,1], quan_cases.df[,4])  #This is a result of default Peason method 
#corr
#corr1 <- cor.test(quan_cases.df[,1], quan_cases.df[,4], method="pearson")   #two variables have to be normal distribution, however, quan_cases.df[,4] is not normal
#corr1
#corr2 <- cor.test(quan_cases.df[,1], quan_cases.df[,4], method="kendall")  #Rank based corelation test
#corr2
#corr3 <- cor.test(quan_cases.df[,1], quan_cases.df[,4], exact=F, method="spearman")  #Rank based corelation test
#corr3
corr <- cor(quan_cases.df[,1], quan_cases.df[,4], method="spearman")  #Rank based corelation test
corr
text(7.5, 0.6, paste0("r=", round(corr,4)))


#+++++++++++++++++++++++++ Non-white (all except white) +++++++++++++++++++++++

non_white <- score_pheno[score_pheno[,7]!="White" ,]
non_white_score <- non_white[,3]

quan=0.02
score_quantile <- quantile(non_white_score, probs=seq(0, 1, quan), names=T)  #11 numbers with %0=samllest_number and %100=largest+number


quan_cases.df <- data.frame()
for (i in (2:length(score_quantile))){
  if (i==2){
    rows <- which(non_white_score <= score_quantile[i])
    quantile_midpoint <- score_quantile[i-1] + (score_quantile[i] - score_quantile[i-1]) / 2
  }
  else{
    rows <- which(non_white_score > score_quantile[i-1] & non_white_score <= score_quantile[i])
    quantile_midpoint <- score_quantile[i-1] + (score_quantile[i] - score_quantile[i-1]) / 2
  }
  quantile_pheno <- non_white[rows, ]
  subject_num <- nrow(quantile_pheno)
  cases <- sum(quantile_pheno[, 4]==1)
  healthy <- sum(quantile_pheno[, 4]==0)
  odds <- cases/healthy
  
  quan_cases.df <- rbind(quan_cases.df, c(quantile_midpoint, cases, subject_num, odds))
}

#Plot density, scatter, fit line
#Note: the case number times 20 to make it scaled with the density curve
prs_score <- quan_cases.df[,1]
odds <- quan_cases.df[,4]
adj <- 12     #scale to match density
odds.adj <- odds * adj

plot(prs_score, odds.adj)

plot(prs_score, odds.adj, xlab="PRS score", ylab=paste0("Odds (x", adj, ")"), xlim=c(4, 8), ylim=c(0, 0.9), main="White ColorectalCA PRS Score Analysis")

par(mar=c(3.5,3.5,1,1))
plot(prs_score, odds.adj, xlab="", ylab="", xlim=c(4.2, 8.2), ylim=c(0, 1.0), main="")
title(xlab="PRS score", line=2, cex.lab=1.0)
title(ylab=paste0("Odds (x", adj, ")"), line=2, cex.lab=1.0)


lines(density(white_score, adjust=2), lty="dotted", col="black")
abline(lm(odds.adj~prs_score), col="red")
shapiro.test(quan_cases.df[,1])   #p value < 0.05, normal distribution
shapiro.test(quan_cases.df[,4])   #p value > 0.05, non-normal distribution
#corr <- cor(quan_cases.df[,1], quan_cases.df[,4])  #This is a result of default Peason method 
#corr
#corr1 <- cor.test(quan_cases.df[,1], quan_cases.df[,4], method="pearson")   #two variables have to be normal distribution, however, quan_cases.df[,4] is not normal
#corr1
#corr2 <- cor.test(quan_cases.df[,1], quan_cases.df[,4], method="kendall")  #Rank based corelation test
#corr2
#corr3 <- cor.test(quan_cases.df[,1], quan_cases.df[,4], exact=F, method="spearman")  #Rank based corelation test
#corr3
corr <- cor(quan_cases.df[,1], quan_cases.df[,4], method="spearman")  #Rank based corelation test
corr
text(7.5, 0.6, paste0("r=", round(corr,4)))





#--------------------- Black ----------------------


black <- score_pheno[score_pheno[,7]=="Black" ,]
black_score <- black[,3]

quan=0.02
score_quantile <- quantile(black_score, probs=seq(0, 1, quan), names=T)  #11 numbers with %0=samllest_number and %100=largest+number


quan_cases.df <- data.frame()
for (i in (2:length(score_quantile))){
  if (i==2){
    rows <- which(black_score <= score_quantile[i])
    quantile_midpoint <- score_quantile[i-1] + (score_quantile[i] - score_quantile[i-1]) / 2
  }
  else{
    rows <- which(black_score > score_quantile[i-1] & black_score <= score_quantile[i])
    quantile_midpoint <- score_quantile[i-1] + (score_quantile[i] - score_quantile[i-1]) / 2
  }
  quantile_pheno <- black[rows, ]
  subject_num <- nrow(quantile_pheno)
  cases <- sum(quantile_pheno[, 4]==1)
  healthy <- sum(quantile_pheno[, 4]==0)
  odds <- cases/healthy
  
  quan_cases.df <- rbind(quan_cases.df, c(quantile_midpoint, cases, subject_num, odds))
}

#Plot density, scatter, fit line
#Note: the case number times 20 to make it scaled with the density curve
prs_score <- quan_cases.df[,1]
odds <- quan_cases.df[,4]
adj <- 10     #scale to match density
odds.adj <- odds * adj

plot(prs_score, odds.adj)

plot(prs_score, odds.adj, xlab="PRS score", ylab=paste0("Odds (x", adj, ")"), xlim=c(4, 8), ylim=c(0, 0.9), main="Black ColorectalCA PRS Score Analysis")

par(mar=c(3.5,3.5,1,1))
plot(prs_score, odds.adj, xlab="", ylab="", xlim=c(4.2, 8.2), ylim=c(0, 1.0), main="")
title(xlab="PRS score", line=2, cex.lab=1.0)
title(ylab=paste0("Odds (x", adj, ")"), line=2, cex.lab=1.0)

lines(density(black_score, adjust=2), col="darkblue")
abline(lm(odds.adj~prs_score), col="red")
corr <- cor(quan_cases.df[,1], quan_cases.df[,4])
text(7.5, 0.5, paste0("r=", round(corr,4)))




#--------------------- Asian ----------------------


asian <- score_pheno[score_pheno[,7]=="Asian" ,]
asian_score <- asian[,3]

quan=0.02
score_quantile <- quantile(asian_score, probs=seq(0, 1, quan), names=T)  #11 numbers with %0=samllest_number and %100=largest+number


quan_cases.df <- data.frame()
for (i in (2:length(score_quantile))){
  if (i==2){
    rows <- which(asian_score <= score_quantile[i])
    quantile_midpoint <- score_quantile[i-1] + (score_quantile[i] - score_quantile[i-1]) / 2
  }
  else{
    rows <- which(asian_score > score_quantile[i-1] & asian_score <= score_quantile[i])
    quantile_midpoint <- score_quantile[i-1] + (score_quantile[i] - score_quantile[i-1]) / 2
  }
  quantile_pheno <- asian[rows, ]
  subject_num <- nrow(quantile_pheno)
  cases <- sum(quantile_pheno[, 4]==1)
  healthy <- sum(quantile_pheno[, 4]==0)
  odds <- cases/healthy
  
  quan_cases.df <- rbind(quan_cases.df, c(quantile_midpoint, cases, subject_num, odds))
}

#Plot density, scatter, fit line
#Note: the case number times 20 to make it scaled with the density curve
prs_score <- quan_cases.df[,1]
odds <- quan_cases.df[,4]
adj <- 10     #scale to match density
odds.adj <- odds * adj

plot(prs_score, odds.adj)

plot(prs_score, odds.adj, xlab="PRS score", ylab=paste0("Odds (x", adj, ")"), xlim=c(4, 8), ylim=c(0, 0.9), main="Asian ColorectalCA PRS Score Analysis")

par(mar=c(3.5,3.5,1,1))
plot(prs_score, odds.adj, xlab="", ylab="", xlim=c(4.2, 8.2), ylim=c(0, 1.0), main="")
title(xlab="PRS score", line=2, cex.lab=1.0)
title(ylab=paste0("Odds (x", adj, ")"), line=2, cex.lab=1.0)

lines(density(asian_score, adjust=2), col="darkblue")
abline(lm(odds.adj~prs_score), col="red")
corr <- cor(quan_cases.df[,1], quan_cases.df[,4])
text(7.5, 0.5, paste0("r=", round(corr,4)))




#--------------------- Other ----------------------


other <- score_pheno[score_pheno[,7]=="Other" ,]
other_score <- other[,3]

quan=0.02
score_quantile <- quantile(other_score, probs=seq(0, 1, quan), names=T)  #11 numbers with %0=samllest_number and %100=largest+number


quan_cases.df <- data.frame()
for (i in (2:length(score_quantile))){
  if (i==2){
    rows <- which(other_score <= score_quantile[i])
    quantile_midpoint <- score_quantile[i-1] + (score_quantile[i] - score_quantile[i-1]) / 2
  }
  else{
    rows <- which(other_score > score_quantile[i-1] & other_score <= score_quantile[i])
    quantile_midpoint <- score_quantile[i-1] + (score_quantile[i] - score_quantile[i-1]) / 2
  }
  quantile_pheno <- other[rows, ]
  subject_num <- nrow(quantile_pheno)
  cases <- sum(quantile_pheno[, 4]==1)
  healthy <- sum(quantile_pheno[, 4]==0)
  odds <- cases/healthy
  
  quan_cases.df <- rbind(quan_cases.df, c(quantile_midpoint, cases, subject_num, odds))
}

#Plot density, scatter, fit line
#Note: the case number times 20 to make it scaled with the density curve
prs_score <- quan_cases.df[,1]
odds <- quan_cases.df[,4]
adj <- 10     #scale to match density
odds.adj <- odds * adj

plot(prs_score, odds.adj)

plot(prs_score, odds.adj, xlab="PRS score", ylab=paste0("Odds (x", adj, ")"), xlim=c(4, 8), ylim=c(0, 0.9), main="Other ColorectalCA PRS Score Analysis")

par(mar=c(3.5,3.5,1,1))
plot(prs_score, odds.adj, xlab="", ylab="", xlim=c(4.2, 8.2), ylim=c(0, 1.0), main="")
title(xlab="PRS score", line=2, cex.lab=1.0)
title(ylab=paste0("Odds (x", adj, ")"), line=2, cex.lab=1.0)


lines(density(other_score, adjust=2), col="darkblue")
abline(lm(odds.adj~prs_score), col="red")
corr <- cor(quan_cases.df[,1], quan_cases.df[,4])
text(7.5, 0.5, paste0("r=", round(corr,4)))



#--------------------- Unknown ----------------------


unknown <- score_pheno[score_pheno[,7]=="Unknown" ,]
unknown_score <- unknown[,3]

quan=0.02
score_quantile <- quantile(unknown_score, probs=seq(0, 1, quan), names=T)  #11 numbers with %0=samllest_number and %100=largest+number


quan_cases.df <- data.frame()
for (i in (2:length(score_quantile))){
  if (i==2){
    rows <- which(unknown_score <= score_quantile[i])
    quantile_midpoint <- score_quantile[i-1] + (score_quantile[i] - score_quantile[i-1]) / 2
  }
  else{
    rows <- which(unknown_score > score_quantile[i-1] & unknown_score <= score_quantile[i])
    quantile_midpoint <- score_quantile[i-1] + (score_quantile[i] - score_quantile[i-1]) / 2
  }
  quantile_pheno <- unknown[rows, ]
  subject_num <- nrow(quantile_pheno)
  cases <- sum(quantile_pheno[, 4]==1)
  healthy <- sum(quantile_pheno[, 4]==0)
  odds <- cases/healthy
  
  quan_cases.df <- rbind(quan_cases.df, c(quantile_midpoint, cases, subject_num, odds))
}

#Plot density, scatter, fit line
#Note: the case number times 20 to make it scaled with the density curve
prs_score <- quan_cases.df[,1]
odds <- quan_cases.df[,4]
adj <- 10     #scale to match density
odds.adj <- odds * adj

plot(prs_score, odds.adj)

plot(prs_score, odds.adj, xlab="PRS score", ylab=paste0("Odds (x", adj, ")"), xlim=c(4, 8), ylim=c(0, 0.9), main="Unknown ColorectalCA PRS Score Analysis")

par(mar=c(3.5,3.5,1,1))
plot(prs_score, odds.adj, xlab="", ylab="", xlim=c(4.2, 8.2), ylim=c(0, 1.0), main="")
title(xlab="PRS score", line=2, cex.lab=1.0)
title(ylab=paste0("Odds (x", adj, ")"), line=2, cex.lab=1.0)

lines(density(unknown_score, adjust=2), col="darkblue")
abline(lm(odds.adj~prs_score), col="red")
corr <- cor(quan_cases.df[,1], quan_cases.df[,4])
text(7.5, 0.5, paste0("r=", round(corr,4)))






#++++++++++++++++ OR and risk analysis for European ancestry for Cutoff ++++++++++++++++++++++

#README file contains how prs_score_colorectalC_covariates_chip.tsv is created
score_pheno <- read.table("prs_score_colorectalC_covariates_chip.tsv", header=T, stringsAsFactors=F)

white_score_pheno <- score_pheno[which(score_pheno[,7]=="White"), ]

white_score <- white_score_pheno[,3]
avg <- mean(white_score)
stddev <- sd(white_score)
avg
stddev

std_score <- (white_score - avg) / stddev  #calculate standardized PRS scores
white_score_pheno$std_score <- std_score   #Add standardized PRS into all_prs table

mylogit <- glm(COLCA~std_score+Gender+Age+chip, data=white_score_pheno, family=binomial)  #Build model
summary(mylogit)   #View the model
#confint(mylogit)         #confident interval of predictors


#Obtain percentage above cutoff
std_score_coef <- coef(mylogit)[2]   #Get coefficient of std_score
attributes(std_score_coef) <- NULL   #Delete attribute "std_score"
OR_change_per_prs_sd <- exp(std_score_coef)
OR_change_per_prs_sd
std_threshold <- log(2.0) / log(OR_change_per_prs_sd)   #2.0 is target raw PRS score cutoff
std_threshold

#Convert std_threshold into raw score
raw_threshold <- (std_threshold * stddev) + avg
raw_threshold

#Calculate percentage of patients in population
above_threshold <- std_score[which(std_score > std_threshold)]
percentage_above_threshold <- (length(above_threshold) / length(std_score)) * 100
percentage_above_threshold







#########################################################################
#########################################################################
#------------ PRS Adjustment by PCA (Four PCs from PLINK2 PCA) -----------
#########################################################################
#########################################################################

setwd("~/Documents/PRS_Jason/prs_analysis/1kg_imputed_normalized/colorectal")

### Run logistic regresion of white PRS for cutoff value to obtain
#std_threshold = 1.916241

### Obtain std_threshold from listerature at OR>2
std_threshold = 1.61

#README file contains how prs_score_prostate_covariates_chip.tsv is created
score_pheno <- read.table("prs_score_colorectalC_covariates_chip.tsv", header=T, stringsAsFactors=F)

setwd("~/Documents/PRS_Jason/Cutoff_adjustment/PCA_adjustment_MDS")
pc <- read.table("plink2pca/ref_pcs.eigenvec.newid", header=T, stringsAsFactor=F)

# Merge two table into one
score_pheno_pc <- merge(score_pheno, pc, by=c("FID", "IID"))

# Fit a linear model with control samples
score_pheno_pc_control <- score_pheno_pc[which(score_pheno_pc$COLCA == 0),]

pcmod <- lm(PRS ~ PC1 + PC2 + PC3 + PC4, data=score_pheno_pc_control)
summary(pcmod)
saveRDS(pcmod, "colorectalC_pc_model.RDS")

# Predicted all samples scores based on the linear model
pred <- predict(pcmod, newdata=score_pheno_pc)

# Calculate residualized scores and add into score_pheno_pc
score_pheno_pc$scoreResid <- score_pheno_pc$PRS - pred

all_resid_mean <- mean(score_pheno_pc$scoreResid)
all_resid_sd <- sd(score_pheno_pc$scoreResid)

all_resid_mean
all_resid_sd


#Convert threshold score from literature into raw PRS score based on 36422 samples in score_pheno, 
#and then convert the raw score into residualized score
all_raw_mean <- mean(score_pheno_pc$PRS)
all_raw_sd <- sd(score_pheno_pc$PRS)

all_raw_mean
all_raw_sd

all_raw_threshold <- (std_threshold * all_raw_sd ) + all_raw_mean
all_raw_threshold

# New threshold by directly applying original std_threshold and destandize it
all_resid_threshold <- std_threshold * all_resid_sd + all_resid_mean
all_resid_threshold

# Assume all std_threshold for standardized raw score and adjusted score
all_std_raw_threshold <- std_threshold
all_std_resid_threshold <- std_threshold




# Get high_risk data
#high.risk <- score_pheno_pc[which(score_pheno_pc$scoreResid > all_resid_threshold),]
#write.csv(high.risk, file="above_cutoff_cc.csv", row.names=F, quote=F)





### ---------- Adjusted PRS Distribution across all racces ------------

#Adjusted PRS for four races
white_resid_score <- score_pheno_pc$scoreResid[which(score_pheno_pc$Race == "White")]
black_resid_score <- score_pheno_pc$scoreResid[which(score_pheno_pc$Race == "Black")]
asian_resid_score <- score_pheno_pc$scoreResid[which(score_pheno_pc$Race == "Asian")]
otherUnknown_resid_score <- score_pheno_pc$scoreResid[c(which(score_pheno_pc$Race == "Other"), which(score_pheno_pc$Race == "Unknown"))]

#Standardiaced adjusted PRS for four races
white_std_resid_score <- (white_resid_score - all_resid_mean) / all_resid_sd
black_std_resid_score <- (black_resid_score - all_resid_mean) / all_resid_sd
asian_std_resid_score <- (asian_resid_score - all_resid_mean) / all_resid_sd
otherUnknown_std_resid_score <- (otherUnknown_resid_score - all_resid_mean) / all_resid_sd

plot(density(white_std_resid_score, adjust=2), col="black",  xlab="Stanndardized adjusted PRS score", lty="solid", main="Density ~ PRS std score of AF")

par(mar=c(3.5,3.5,1,1))
plot(density(white_std_resid_score, adjust=2), col="black",  xlim=c(-5, 5), ylim=c(0, 0.6), xlab="", ylab="", lty="solid", main="")
title(xlab="Standardized adjusted PRS", line=2, cex.lab=1.0)
title(ylab="Density", line=2, cex.lab=1.0)

lines(density(black_std_resid_score, adjust=2), lty="dotted", col="purple")
lines(density(asian_std_resid_score, adjust=2), lty="dotted", col="green")
lines(density(otherUnknown_std_resid_score, adjust=2), lty="dotted", col="blue")

#Add cutoff value as a red line
abline(v=all_std_resid_threshold, col='red')
text(3.5, 0.2, paste0("cutoff=", round(all_std_resid_threshold,2)), cex = 0.7)


#Add legend
Vars <- c("White", "Black", "Asian", "Other/\nUnknown")  # one per row
colors <- c("black", "purple", "green", "blue")  # one color for each row; should be same length as Vars
ltys <- c("solid","dotted","dotted", "dotted")  # one lty for each column

nr <- length(Vars) 
legend("topright", rep(Vars, 1), col = colors, lty = ltys, ncol = 1, cex = 0.50)

text(-4, 0.55, paste0("CRCa"))
#text(-0.4, 3.2, paste0("All races"))


###-------------------------------------------------------------------------------------------------------
### Distribution of raw PRS scores across Races Analysis for comparison to adjusted score distribution ------------------

#README file contains how prs_score_prostate_covariates_chip.tsv is created

white_raw_score <- score_pheno[which(score_pheno[,7]=="White") ,3]
black_raw_score <- score_pheno[which(score_pheno[,7]=="Black") ,3]
asian_raw_score <- score_pheno[which(score_pheno[,7]=="Asian") ,3]
otherUnknown_raw_score <- score_pheno[c(which(score_pheno[,7]=="Other"), which(score_pheno[,7]=="Unknown")) ,3]


#Standardized raw PRS score
white_std_raw_score <- (white_raw_score - all_raw_mean) / all_raw_sd
black_std_raw_score <- (black_raw_score - all_raw_mean) / all_raw_sd
asian_std_raw_score <- (asian_raw_score - all_raw_mean) / all_raw_sd
otherUnknown_std_raw_score <- (otherUnknown_raw_score - all_raw_mean) / all_raw_sd

#Barplot of races
score_pheno_tmp <- score_pheno
score_pheno_tmp[c(which(score_pheno_tmp$Race=="Other"), which(score_pheno_tmp$Race=="Unknown")), 7] <- "other/Unknown"
race_table <- table(score_pheno_tmp[,7])

par(mar=c(4.5,4.5,1,1))
x <- barplot(race_table, ylab="Number of subjects", main="Barplot of Race", xaxt="n")

labs <- c("Asian", "Black", "other/\nUnknown", "White")
text(cex=1, x=x, y=2.5, adj=c(1.1,1), labs, xpd=T, srt=50)


###Raw score Deansity plot across Races -----------
plot(density(white_std_raw_score, adjust=2))
plot(density(white_std_raw_score, adjust=2), col="black", xlim=c(-5, 5), ylim=c(0, 0.45), xlab="PRS score", lty="solid", main="Density ~ PRS score of CAD")

par(mar=c(3.5,3.5,1,1))
plot(density(white_std_raw_score, adjust=2), col="black", xlim=c(-5, 5), ylim=c(0, 0.6), xlab="",ylab="", lty="solid", main="")
title(xlab="Standardized raw PRS", line=2, cex.lab=1.0)
title(ylab="Density", line=2, cex.lab=1.0)

lines(density(black_std_raw_score, adjust=2), lty="dotted", col="purple")
lines(density(asian_std_raw_score, adjust=2), lty="dotted", col="green")
lines(density(otherUnknown_std_raw_score, adjust=2), lty="dotted", col="blue")

#Add cutoff value as a red line
abline(v=all_std_raw_threshold, col='red')
text(3.5, 0.2, paste0("cutoff=", round(all_std_raw_threshold,2)), cex = 0.7)

#Add legend
Vars <- c("White", "Black", "Asian", "Other/\nUnknown")  # one per row
colors <- c("black", "purple", "green", "blue")  # one color for each row; should be same length as Vars
ltys <- c("solid","dotted","dotted", "dotted")  # one lty for each column

nr <- length(Vars) 
legend("topright", rep(Vars, 1), col = colors, lty = ltys, ncol = 1, cex = 0.50)

text(-4, 0.55, paste0("CRCa"))






###--------------------------- Odds ratio in all and individual races -----------------------------------------
###-------------------------------------------------------------------------------------------------------------------

#Percentage of high risk subjects in MGB biobank, using threshold from literature (converted raw score or PC-adjusted score
#in MGB biobank)
raw_cutoff <- all_raw_threshold    #All the races use the same cutoff as White
adj_cutoff <- all_resid_threshold  #All the races use the same cutoff as White

high_risk_rate.raw <- sum(score_pheno_pc$PRS > raw_cutoff) / length(score_pheno_pc$PRS) * 100
high_risk_rate.adj <- sum(score_pheno_pc$scoreResid > adj_cutoff) / length(score_pheno_pc$scoreResid) * 100

high_risk_rate.raw
high_risk_rate.adj


#Odds ratio for all races using raw score threshold or adjusted score threshold
above_cutoff.raw <- score_pheno_pc[which(score_pheno_pc[,3] > raw_cutoff),]
below_cutoff.raw <- score_pheno_pc[which(score_pheno_pc[,3] < raw_cutoff | score_pheno_pc[,3] == raw_cutoff),]
above_cutoff.odds.raw <- sum(above_cutoff.raw[,4] == 1) / sum(above_cutoff.raw[,4] == 0)
below_cutoff.odds.raw <- sum(below_cutoff.raw[,4] == 1) / sum(below_cutoff.raw[,4] == 0)
overall_OR.raw <- above_cutoff.odds.raw / below_cutoff.odds.raw
overall_OR.raw
lower.95ci.raw <- exp((log(overall_OR.raw) - 1.96 * sqrt(1/sum(above_cutoff.raw[,4] == 1) + 1/sum(above_cutoff.raw[,4] == 0) + 1/sum(below_cutoff.raw[,4] == 1) + 1/sum(below_cutoff.raw[,4] == 0))))
upper.95ci.raw <- exp((log(overall_OR.raw) + 1.96 * sqrt(1/sum(above_cutoff.raw[,4] == 1) + 1/sum(above_cutoff.raw[,4] == 0) + 1/sum(below_cutoff.raw[,4] == 1) + 1/sum(below_cutoff.raw[,4] == 0))))
print(paste0(round(overall_OR.raw,2), "[", round(lower.95ci.raw,2), ", ", round(upper.95ci.raw,2), "] ", "(", sum(above_cutoff.raw[,4] == 1), "/", sum(above_cutoff.raw[,4] == 0), ", ", sum(below_cutoff.raw[,4] == 1), "/", sum(below_cutoff.raw[,4] == 0), ")"))

above_cutoff.adj <- score_pheno_pc[which(score_pheno_pc[,15] > adj_cutoff), ]
#write.csv(above_cutoff.adj, file="high_risk_analysis/above_cutoff_cc.csv", row.names=F, quote=F)
below_cutoff.adj <- score_pheno_pc[which(score_pheno_pc[,15] < adj_cutoff | score_pheno_pc[,15] == adj_cutoff),]
above_cutoff.odds.adj <- sum(above_cutoff.adj[,4] == 1) / sum(above_cutoff.adj[,4] == 0)
below_cutoff.odds.adj <- sum(below_cutoff.adj[,4] == 1) / sum(below_cutoff.adj[,4] == 0)
overall_OR.adj <- above_cutoff.odds.adj / below_cutoff.odds.adj
overall_OR.adj
lower.95ci.adj <- exp((log(overall_OR.adj) - 1.96 * sqrt(1/sum(above_cutoff.adj[,4] == 1) + 1/sum(above_cutoff.adj[,4] == 0) + 1/sum(below_cutoff.adj[,4] == 1) + 1/sum(below_cutoff.adj[,4] == 0))))
upper.95ci.adj <- exp((log(overall_OR.adj) + 1.96 * sqrt(1/sum(above_cutoff.adj[,4] == 1) + 1/sum(above_cutoff.adj[,4] == 0) + 1/sum(below_cutoff.adj[,4] == 1) + 1/sum(below_cutoff.adj[,4] == 0))))
print(paste0(round(overall_OR.adj,2), "[", round(lower.95ci.adj,2), ", ", round(upper.95ci.adj,2), "] ", "(", sum(above_cutoff.adj[,4] == 1), "/", sum(above_cutoff.adj[,4] == 0), ", ", sum(below_cutoff.adj[,4] == 1), "/", sum(below_cutoff.adj[,4] == 0), ")"))


#Odds ratio for White using raw score threshold or adjusted score threshold
white_score_pheno_pc <- score_pheno_pc[which(score_pheno_pc$Race == "White"),]

white_above_cutoff.raw <- white_score_pheno_pc[which(white_score_pheno_pc[,3] > raw_cutoff),]
white_below_cutoff.raw <- white_score_pheno_pc[which(white_score_pheno_pc[,3] < raw_cutoff | white_score_pheno_pc[,3] == raw_cutoff),]
white_above_cutoff.odds.raw <- sum(white_above_cutoff.raw[,4] == 1) / sum(white_above_cutoff.raw[,4] == 0)
white_below_cutoff.odds.raw <- sum(white_below_cutoff.raw[,4] == 1) / sum(white_below_cutoff.raw[,4] == 0)
white_overall_OR.raw <- white_above_cutoff.odds.raw / white_below_cutoff.odds.raw
white_overall_OR.raw
lower.95ci.raw <- exp((log(white_overall_OR.raw) - 1.96 * sqrt(1/sum(white_above_cutoff.raw[,4] == 1) + 1/sum(white_above_cutoff.raw[,4] == 0) + 1/sum(white_below_cutoff.raw[,4] == 1) + 1/sum(white_below_cutoff.raw[,4] == 0))))
upper.95ci.raw <- exp((log(white_overall_OR.raw) + 1.96 * sqrt(1/sum(white_above_cutoff.raw[,4] == 1) + 1/sum(white_above_cutoff.raw[,4] == 0) + 1/sum(white_below_cutoff.raw[,4] == 1) + 1/sum(white_below_cutoff.raw[,4] == 0))))
print(paste0(round(white_overall_OR.raw,2), "[", round(lower.95ci.raw,2), ", ", round(upper.95ci.raw,2), "] ", "(", sum(white_above_cutoff.raw[,4] == 1), "/", sum(white_above_cutoff.raw[,4] == 0), ", ", sum(white_below_cutoff.raw[,4] == 1), "/", sum(white_below_cutoff.raw[,4] == 0), ")"))

white_above_cutoff.adj <- white_score_pheno_pc[which(white_score_pheno_pc[,15] > adj_cutoff), ]
white_below_cutoff.adj <- white_score_pheno_pc[which(white_score_pheno_pc[,15] < adj_cutoff | white_score_pheno_pc[,15] == adj_cutoff),]
white_above_cutoff.odds.adj <- sum(white_above_cutoff.adj[,4] == 1) / sum(white_above_cutoff.adj[,4] == 0)
white_below_cutoff.odds.adj <- sum(white_below_cutoff.adj[,4] == 1) / sum(white_below_cutoff.adj[,4] == 0)
white_overall_OR.adj <- white_above_cutoff.odds.adj / white_below_cutoff.odds.adj
white_overall_OR.adj
lower.95ci.adj <- exp((log(white_overall_OR.adj) - 1.96 * sqrt(1/sum(white_above_cutoff.adj[,4] == 1) + 1/sum(white_above_cutoff.adj[,4] == 0) + 1/sum(white_below_cutoff.adj[,4] == 1) + 1/sum(white_below_cutoff.adj[,4] == 0))))
upper.95ci.adj <- exp((log(white_overall_OR.adj) + 1.96 * sqrt(1/sum(white_above_cutoff.adj[,4] == 1) + 1/sum(white_above_cutoff.adj[,4] == 0) + 1/sum(white_below_cutoff.adj[,4] == 1) + 1/sum(white_below_cutoff.adj[,4] == 0))))
print(paste0(round(white_overall_OR.adj,2), "[", round(lower.95ci.adj,2), ", ", round(upper.95ci.adj,2), "] ", "(", sum(white_above_cutoff.adj[,4] == 1), "/", sum(white_above_cutoff.adj[,4] == 0), ", ", sum(white_below_cutoff.adj[,4] == 1), "/", sum(white_below_cutoff.adj[,4] == 0), ")"))


#Odds ratio for Black using raw score threshold or adjusted score threshold
black_score_pheno_pc <- score_pheno_pc[which(score_pheno_pc$Race == "Black"),]

black_above_cutoff.raw <- black_score_pheno_pc[which(black_score_pheno_pc[,3] > raw_cutoff),]
black_below_cutoff.raw <- black_score_pheno_pc[which(black_score_pheno_pc[,3] < raw_cutoff | black_score_pheno_pc[,3] == raw_cutoff),]
black_above_cutoff.odds.raw <- sum(black_above_cutoff.raw[,4] == 1) / sum(black_above_cutoff.raw[,4] == 0)
black_below_cutoff.odds.raw <- sum(black_below_cutoff.raw[,4] == 1) / sum(black_below_cutoff.raw[,4] == 0)
black_overall_OR.raw <- black_above_cutoff.odds.raw / black_below_cutoff.odds.raw
black_overall_OR.raw
lower.95ci.raw <- exp((log(black_overall_OR.raw) - 1.96 * sqrt(1/sum(black_above_cutoff.raw[,4] == 1) + 1/sum(black_above_cutoff.raw[,4] == 0) + 1/sum(black_below_cutoff.raw[,4] == 1) + 1/sum(black_below_cutoff.raw[,4] == 0))))
upper.95ci.raw <- exp((log(black_overall_OR.raw) + 1.96 * sqrt(1/sum(black_above_cutoff.raw[,4] == 1) + 1/sum(black_above_cutoff.raw[,4] == 0) + 1/sum(black_below_cutoff.raw[,4] == 1) + 1/sum(black_below_cutoff.raw[,4] == 0))))
print(paste0(round(black_overall_OR.raw,2), "[", round(lower.95ci.raw,2), ", ", round(upper.95ci.raw,2), "] ", "(", sum(black_above_cutoff.raw[,4] == 1), "/", sum(black_above_cutoff.raw[,4] == 0), ", ", sum(black_below_cutoff.raw[,4] == 1), "/", sum(black_below_cutoff.raw[,4] == 0), ")"))

black_above_cutoff.adj <- black_score_pheno_pc[which(black_score_pheno_pc[,15] > adj_cutoff), ]
black_below_cutoff.adj <- black_score_pheno_pc[which(black_score_pheno_pc[,15] < adj_cutoff | black_score_pheno_pc[,15] == adj_cutoff),]
black_above_cutoff.odds.adj <- sum(black_above_cutoff.adj[,4] == 1) / sum(black_above_cutoff.adj[,4] == 0)
black_below_cutoff.odds.adj <- sum(black_below_cutoff.adj[,4] == 1) / sum(black_below_cutoff.adj[,4] == 0)
black_overall_OR.adj <- black_above_cutoff.odds.adj / black_below_cutoff.odds.adj
black_overall_OR.adj
lower.95ci.adj <- exp((log(black_overall_OR.adj) - 1.96 * sqrt(1/sum(black_above_cutoff.adj[,4] == 1) + 1/sum(black_above_cutoff.adj[,4] == 0) + 1/sum(black_below_cutoff.adj[,4] == 1) + 1/sum(black_below_cutoff.adj[,4] == 0))))
upper.95ci.adj <- exp((log(black_overall_OR.adj) + 1.96 * sqrt(1/sum(black_above_cutoff.adj[,4] == 1) + 1/sum(black_above_cutoff.adj[,4] == 0) + 1/sum(black_below_cutoff.adj[,4] == 1) + 1/sum(black_below_cutoff.adj[,4] == 0))))
print(paste0(round(black_overall_OR.adj,2), "[", round(lower.95ci.adj,2), ", ", round(upper.95ci.adj,2), "] ", "(", sum(black_above_cutoff.adj[,4] == 1), "/", sum(black_above_cutoff.adj[,4] == 0), ", ", sum(black_below_cutoff.adj[,4] == 1), "/", sum(black_below_cutoff.adj[,4] == 0), ")"))


#Odds ratio for ASIAN using raw score threshold or adjusted score threshold
asian_score_pheno_pc <- score_pheno_pc[which(score_pheno_pc$Race == "Asian"),]

asian_above_cutoff.raw <- asian_score_pheno_pc[which(asian_score_pheno_pc[,3] > raw_cutoff),]
asian_below_cutoff.raw <- asian_score_pheno_pc[which(asian_score_pheno_pc[,3] < raw_cutoff | asian_score_pheno_pc[,3] == raw_cutoff),]
asian_above_cutoff.odds.raw <- sum(asian_above_cutoff.raw[,4] == 1) / sum(asian_above_cutoff.raw[,4] == 0)
asian_below_cutoff.odds.raw <- sum(asian_below_cutoff.raw[,4] == 1) / sum(asian_below_cutoff.raw[,4] == 0)
asian_overall_OR.raw <- asian_above_cutoff.odds.raw / asian_below_cutoff.odds.raw
asian_overall_OR.raw
lower.95ci.raw <- exp((log(asian_overall_OR.raw) - 1.96 * sqrt(1/sum(asian_above_cutoff.raw[,4] == 1) + 1/sum(asian_above_cutoff.raw[,4] == 0) + 1/sum(asian_below_cutoff.raw[,4] == 1) + 1/sum(asian_below_cutoff.raw[,4] == 0))))
upper.95ci.raw <- exp((log(asian_overall_OR.raw) + 1.96 * sqrt(1/sum(asian_above_cutoff.raw[,4] == 1) + 1/sum(asian_above_cutoff.raw[,4] == 0) + 1/sum(asian_below_cutoff.raw[,4] == 1) + 1/sum(asian_below_cutoff.raw[,4] == 0))))
print(paste0(round(asian_overall_OR.raw,2), "[", round(lower.95ci.raw,2), ", ", round(upper.95ci.raw,2), "] ", "(", sum(asian_above_cutoff.raw[,4] == 1), "/", sum(asian_above_cutoff.raw[,4] == 0), ", ", sum(asian_below_cutoff.raw[,4] == 1), "/", sum(asian_below_cutoff.raw[,4] == 0), ")"))

asian_above_cutoff.adj <- asian_score_pheno_pc[which(asian_score_pheno_pc[,15] > adj_cutoff), ]
asian_below_cutoff.adj <- asian_score_pheno_pc[which(asian_score_pheno_pc[,15] < adj_cutoff | asian_score_pheno_pc[,15] == adj_cutoff),]
asian_above_cutoff.odds.adj <- sum(asian_above_cutoff.adj[,4] == 1) / sum(asian_above_cutoff.adj[,4] == 0)
asian_below_cutoff.odds.adj <- sum(asian_below_cutoff.adj[,4] == 1) / sum(asian_below_cutoff.adj[,4] == 0)
asian_overall_OR.adj <- asian_above_cutoff.odds.adj / asian_below_cutoff.odds.adj
asian_overall_OR.adj
lower.95ci.adj <- exp((log(asian_overall_OR.adj) - 1.96 * sqrt(1/sum(asian_above_cutoff.adj[,4] == 1) + 1/sum(asian_above_cutoff.adj[,4] == 0) + 1/sum(asian_below_cutoff.adj[,4] == 1) + 1/sum(asian_below_cutoff.adj[,4] == 0))))
upper.95ci.adj <- exp((log(asian_overall_OR.adj) + 1.96 * sqrt(1/sum(asian_above_cutoff.adj[,4] == 1) + 1/sum(asian_above_cutoff.adj[,4] == 0) + 1/sum(asian_below_cutoff.adj[,4] == 1) + 1/sum(asian_below_cutoff.adj[,4] == 0))))
print(paste0(round(asian_overall_OR.adj,2), "[", round(lower.95ci.adj,2), ", ", round(upper.95ci.adj,2), "] ", "(", sum(asian_above_cutoff.adj[,4] == 1), "/", sum(asian_above_cutoff.adj[,4] == 0), ", ", sum(asian_below_cutoff.adj[,4] == 1), "/", sum(asian_below_cutoff.adj[,4] == 0), ")"))

#Odds ratio for OTHER/UNKNOWN using raw score threshold or adjusted score threshold
otherUnknown_score_pheno_pc <- score_pheno_pc[c(which(score_pheno_pc$Race == "Other"), which(score_pheno_pc$Race == "Unknown")),]

otherUnknown_above_cutoff.raw <- otherUnknown_score_pheno_pc[which(otherUnknown_score_pheno_pc[,3] > raw_cutoff),]
otherUnknown_below_cutoff.raw <- otherUnknown_score_pheno_pc[which(otherUnknown_score_pheno_pc[,3] < raw_cutoff | otherUnknown_score_pheno_pc[,3] == raw_cutoff),]
otherUnknown_above_cutoff.odds.raw <- sum(otherUnknown_above_cutoff.raw[,4] == 1) / sum(otherUnknown_above_cutoff.raw[,4] == 0)
otherUnknown_below_cutoff.odds.raw <- sum(otherUnknown_below_cutoff.raw[,4] == 1) / sum(otherUnknown_below_cutoff.raw[,4] == 0)
otherUnknown_overall_OR.raw <- otherUnknown_above_cutoff.odds.raw / otherUnknown_below_cutoff.odds.raw
otherUnknown_overall_OR.raw
lower.95ci.raw <- exp((log(otherUnknown_overall_OR.raw) - 1.96 * sqrt(1/sum(otherUnknown_above_cutoff.raw[,4] == 1) + 1/sum(otherUnknown_above_cutoff.raw[,4] == 0) + 1/sum(otherUnknown_below_cutoff.raw[,4] == 1) + 1/sum(otherUnknown_below_cutoff.raw[,4] == 0))))
upper.95ci.raw <- exp((log(otherUnknown_overall_OR.raw) + 1.96 * sqrt(1/sum(otherUnknown_above_cutoff.raw[,4] == 1) + 1/sum(otherUnknown_above_cutoff.raw[,4] == 0) + 1/sum(otherUnknown_below_cutoff.raw[,4] == 1) + 1/sum(otherUnknown_below_cutoff.raw[,4] == 0))))
print(paste0(round(otherUnknown_overall_OR.raw,2), "[", round(lower.95ci.raw,2), ", ", round(upper.95ci.raw,2), "] ", "(", sum(otherUnknown_above_cutoff.raw[,4] == 1), "/", sum(otherUnknown_above_cutoff.raw[,4] == 0), ", ", sum(otherUnknown_below_cutoff.raw[,4] == 1), "/", sum(otherUnknown_below_cutoff.raw[,4] == 0), ")"))

otherUnknown_above_cutoff.adj <- otherUnknown_score_pheno_pc[which(otherUnknown_score_pheno_pc[,15] > adj_cutoff), ]
otherUnknown_below_cutoff.adj <- otherUnknown_score_pheno_pc[which(otherUnknown_score_pheno_pc[,15] < adj_cutoff | otherUnknown_score_pheno_pc[,15] == adj_cutoff),]
otherUnknown_above_cutoff.odds.adj <- sum(otherUnknown_above_cutoff.adj[,4] == 1) / sum(otherUnknown_above_cutoff.adj[,4] == 0)
otherUnknown_below_cutoff.odds.adj <- sum(otherUnknown_below_cutoff.adj[,4] == 1) / sum(otherUnknown_below_cutoff.adj[,4] == 0)
otherUnknown_overall_OR.adj <- otherUnknown_above_cutoff.odds.adj / otherUnknown_below_cutoff.odds.adj
otherUnknown_overall_OR.adj
lower.95ci.adj <- exp((log(otherUnknown_overall_OR.adj) - 1.96 * sqrt(1/sum(otherUnknown_above_cutoff.adj[,4] == 1) + 1/sum(otherUnknown_above_cutoff.adj[,4] == 0) + 1/sum(otherUnknown_below_cutoff.adj[,4] == 1) + 1/sum(otherUnknown_below_cutoff.adj[,4] == 0))))
upper.95ci.adj <- exp((log(otherUnknown_overall_OR.adj) + 1.96 * sqrt(1/sum(otherUnknown_above_cutoff.adj[,4] == 1) + 1/sum(otherUnknown_above_cutoff.adj[,4] == 0) + 1/sum(otherUnknown_below_cutoff.adj[,4] == 1) + 1/sum(otherUnknown_below_cutoff.adj[,4] == 0))))
print(paste0(round(otherUnknown_overall_OR.adj,2), "[", round(lower.95ci.adj,2), ", ", round(upper.95ci.adj,2), "] ", "(", sum(otherUnknown_above_cutoff.adj[,4] == 1), "/", sum(otherUnknown_above_cutoff.adj[,4] == 0), ", ", sum(otherUnknown_below_cutoff.adj[,4] == 1), "/", sum(otherUnknown_below_cutoff.adj[,4] == 0), ")"))







###------------------------------------------------------------------------------------------------------------------------
###------------------- Standardized all races adjusted PRS score distribution and correlation with LogOdds -----------------------
###-------------------------------------------------------------------------------------------------------------------

#Use standardized adjusted scores and logodds to plot  
all_resid_non_standardized_score <- score_pheno_pc$scoreResid
all_resid_mean <- mean(all_resid_non_standardized_score)
all_resid_sd <- sd(all_resid_non_standardized_score)

#Notice, here all_resid_score is standardized based on all population
all_resid_score <- (all_resid_non_standardized_score - all_resid_mean) / all_resid_sd

quan=0.02
score_quantile <- quantile(all_resid_score, probs=seq(0, 1, quan), names=T)  #11 numbers with %0=samllest_number and %100=largest+number


quan_cases.df <- data.frame()
for (i in (2:length(score_quantile))){
  if (i==2){
    rows <- which(all_resid_score <= score_quantile[i])
    quantile_midpoint <- score_quantile[i-1] + (score_quantile[i] - score_quantile[i-1]) / 2
  }
  else{
    rows <- which(all_resid_score > score_quantile[i-1] & all_resid_score <= score_quantile[i])
    quantile_midpoint <- score_quantile[i-1] + (score_quantile[i] - score_quantile[i-1]) / 2
  }
  quantile_pheno <- score_pheno_pc[rows, ]
  subject_num <- nrow(quantile_pheno)
  cases <- sum(quantile_pheno[, 4]==1)
  healthy <- sum(quantile_pheno[, 4]==0)
  p <- cases/(healthy + cases)
  odds <- p/(1-p)
  logodds <- log(odds)
  quan_cases.df <- rbind(quan_cases.df, c(quantile_midpoint, cases, subject_num, logodds))
}

#Plot density, scatter, fit line
#Note: the case number times 20 to make it scaled with the density curve
quan_cases.df <- quan_cases.df[which(is.finite(quan_cases.df[,4])),]  #Remove Infinite rows
print(nrow(quan_cases.df))
prs_score <- quan_cases.df[,1]
logodds <- quan_cases.df[,4]


par(mar=c(5.5,5.5,1,1))
plot(prs_score, logodds)

plot(prs_score, logodds, xlab="PRS score", ylab="Log-odds",  xlim=c(-4, 4), ylim=c(0, 4), main="AFib PRS Score Analysis")
par(mar=c(3.5,3.5,1,1))
plot(prs_score, logodds, xlab="", ylab="",  xlim=c(-4, 4), ylim=c(-7, 0), main="")
title(xlab="Standardized adjusted PRS", line=2, cex.lab=1.0)
title(ylab="Log-odds", line=2, cex.lab=1.0)

abline(lm(logodds~prs_score), col="red")

corr <- cor(quan_cases.df[,1], quan_cases.df[,4])
text(2.5, -2.5, paste0("r=", round(corr,4)))


text(-3.2, -0.5, paste0("CRCa"))






#----------------------Standardized White distribution of residualized scores and correlation with LogOdds -------------------------

white_score_pheno_pc <- score_pheno_pc[(score_pheno_pc$Race == "White"), ]

#Use standardized adjusted scores to plot.
#Notice: The mean and sd are all_population, not only White
white_resid_non_standardized_score <- white_score_pheno_pc$scoreResid
white_resid_score <- (white_resid_non_standardized_score - all_resid_mean) / all_resid_sd  #this is standardized score

quan=0.02
score_quantile <- quantile(white_resid_score, probs=seq(0, 1, quan), names=T)  #11 numbers with %0=samllest_number and %100=largest+number


quan_cases.df <- data.frame()
for (i in (2:length(score_quantile))){
  if (i==2){
    rows <- which(white_resid_score <= score_quantile[i])
    quantile_midpoint <- score_quantile[i-1] + (score_quantile[i] - score_quantile[i-1]) / 2
    print(i)
    #print(rows)
    print(quantile_midpoint)
  }
  else{
    rows <- which(white_resid_score > score_quantile[i-1] & white_resid_score <= score_quantile[i])
    quantile_midpoint <- score_quantile[i-1] + (score_quantile[i] - score_quantile[i-1]) / 2
  }
  quantile_pheno <- white_score_pheno_pc[rows, ]
  subject_num <- nrow(quantile_pheno)
  cases <- sum(quantile_pheno[, 4]==1)
  healthy <- sum(quantile_pheno[, 4]==0)
  
  p <- cases/(healthy + cases)
  odds <- p/(1-p)
  logodds <- log(odds)
  quan_cases.df <- rbind(quan_cases.df, c(quantile_midpoint, cases, subject_num, logodds))
}

#Plot density, scatter, fit line
#Note: the case number times 20 to make it scaled with the density curve

quan_cases.df <- quan_cases.df[which(is.finite(quan_cases.df[,4])),]  #Remove Infinite rows
print(nrow(quan_cases.df))
prs_score <- quan_cases.df[,1]
logodds <- quan_cases.df[,4]

#plot(prs_score, odds.adj)

#plot(prs_score, odds.adj, xlab="PRS score", ylab=paste0("Odds (x", adj, ")"), xlim=c(-0.6, 0.6), ylim=c(0, 4), main="White AFib PRS Score Analysis")
par(mar=c(3.5,3.5,1,1))

plot(prs_score, logodds, xlab="", ylab="",  xlim=c(-4, 4), ylim=c(-7, 0), main="")
title(xlab="Standardized adjusted PRS", line=2, cex.lab=1.0)
title(ylab="Log-odds", line=2, cex.lab=1.0)


abline(lm(logodds~prs_score), col="red")

corr <- cor(quan_cases.df[,1], quan_cases.df[,4])
text(2.5, -2.5, paste0("r=", round(corr,4)))

text(-3.1, -0.5, paste0("CRCa"))
text(-3.12, -1, paste0("White"))





#---------------------- Standardized Black distribution of residualized scores and correlation with LogOdds -------------------------

black_score_pheno_pc <-score_pheno_pc[(score_pheno_pc$Race == "Black"), ]

#Use standardized adjusted scores to plot 
#Notice: The mean and sd are all_population, not only Black
black_resid_non_standardized_score <- black_score_pheno_pc$scoreResid
black_resid_score <- (black_resid_non_standardized_score - all_resid_mean) / all_resid_sd  #this is standardized score


#quan=0.02
quan=0.2
score_quantile <- quantile(black_resid_score, probs=seq(0, 1, quan), names=T)  #11 numbers with %0=samllest_number and %100=largest+number


quan_cases.df <- data.frame()
for (i in (2:length(score_quantile))){
  if (i==2){
    rows <- which(black_resid_score <= score_quantile[i])
    quantile_midpoint <- score_quantile[i-1] + (score_quantile[i] - score_quantile[i-1]) / 2
  }
  else{
    rows <- which(black_resid_score > score_quantile[i-1] & black_resid_score <= score_quantile[i])
    quantile_midpoint <- score_quantile[i-1] + (score_quantile[i] - score_quantile[i-1]) / 2
  }
  quantile_pheno <- black_score_pheno_pc[rows, ]
  subject_num <- nrow(quantile_pheno)
  cases <- sum(quantile_pheno[, 4]==1)
  healthy <- sum(quantile_pheno[, 4]==0)
  
  p <- cases/(healthy + cases)
  odds <- p/(1-p)
  logodds <- log(odds)
  quan_cases.df <- rbind(quan_cases.df, c(quantile_midpoint, cases, subject_num, logodds))
}

#Plot density, scatter, fit line
#Note: the case number times 20 to make it scaled with the density curve
quan_cases.df <- quan_cases.df[which(is.finite(quan_cases.df[,4])),]  #Remove Infinite rows
print(nrow(quan_cases.df))
prs_score <- quan_cases.df[,1]
logodds <- quan_cases.df[,4]

#plot(prs_score, odds.adj)

#plot(prs_score, odds.adj, xlab="PRS score", ylab=paste0("Odds (x", adj, ")"), xlim=c(-0.6, 0.6), ylim=c(0, 4), main="Black AFibPRS Score Analysis")
par(mar=c(3.5,3.5,1,1))
plot(prs_score, logodds, xlab="", ylab="",  xlim=c(-4, 4), ylim=c(-7, 0), main="")
title(xlab="Standardized adjusted PRS", line=2, cex.lab=1.0)
title(ylab="Log-odds", line=2, cex.lab=1.0)


abline(lm(logodds~prs_score), col="red")

corr <- cor(quan_cases.df[,1], quan_cases.df[,4])
text(2.7, -3, paste0("r=", round(corr,4)))

text(-3.1, -0.5, paste0("CRCa"))
text(-3.12, -1, paste0("Black"))
text(-3.1, -1.6, paste0("q=", quan))



#---------------------- Standardized Asian distribution of residualized scores and correlation with LogOdds -------------------------

asian_score_pheno_pc <- score_pheno_pc[(score_pheno_pc$Race == "Asian"), ]

#Use standardized adjusted scores to plot 
#Notice: The mean and sd is all_population, not only Asian
asian_resid_non_standardized_score <- asian_score_pheno_pc$scoreResid
asian_resid_score <- (asian_resid_non_standardized_score - all_resid_mean) / all_resid_sd  #this is standardized score


#quan=0.02
quan=0.2
score_quantile <- quantile(asian_resid_score, probs=seq(0, 1, quan), names=T)  #11 numbers with %0=samllest_number and %100=largest+number


quan_cases.df <- data.frame()
for (i in (2:length(score_quantile))){
  if (i==2){
    rows <- which(asian_resid_score <= score_quantile[i])
    quantile_midpoint <- score_quantile[i-1] + (score_quantile[i] - score_quantile[i-1]) / 2
  }
  else{
    rows <- which(asian_resid_score > score_quantile[i-1] & asian_resid_score <= score_quantile[i])
    quantile_midpoint <- score_quantile[i-1] + (score_quantile[i] - score_quantile[i-1]) / 2
  }
  quantile_pheno <- asian_score_pheno_pc[rows, ]
  subject_num <- nrow(quantile_pheno)
  cases <- sum(quantile_pheno[, 4]==1)
  healthy <- sum(quantile_pheno[, 4]==0)
  
  p <- cases/(healthy + cases)
  odds <- p/(1-p)
  logodds <- log(odds)
  quan_cases.df <- rbind(quan_cases.df, c(quantile_midpoint, cases, subject_num, logodds))
}
#print(quan_cases.df)
#Plot density, scatter, fit line
#Note: the case number times 20 to make it scaled with the density curve
quan_cases.df <- quan_cases.df[which(is.finite(quan_cases.df[,4])),]  #Remove Infinite rows
print(nrow(quan_cases.df))
prs_score <- quan_cases.df[,1]
logodds <- quan_cases.df[,4]


#plot(prs_score, odds.adj)

#plot(prs_score, odds.adj, xlab="PRS score", ylab=paste0("Odds (x", adj, ")"), xlim=c(-0.6, 0.6), ylim=c(0, 4), main="Asian AFib PRS Score Analysis")
par(mar=c(3.5,3.5,1,1))
plot(prs_score, logodds, xlab="", ylab="",  xlim=c(-4, 4), ylim=c(-7, 0), main="")
title(xlab="Standardized adjusted PRS", line=2, cex.lab=1.0)
title(ylab="Log-odds", line=2, cex.lab=1.0)


abline(lm(logodds~prs_score), col="red")

corr <- cor(quan_cases.df[,1], quan_cases.df[,4])
text(2.7, -5.5, paste0("r=", round(corr,4)))

text(-3.1, -0.5, paste0("CRCa"))
text(-3.15, -1, paste0("Asian"))
text(-3.1, -1.6, paste0("q=", quan))




#---------------------- Standardized Other/Unknown distribution of residualized scores and correlation with LogOdds -------------------------
#Combine Other and Unknown data

otherUnknown_score_pheno_pc <- score_pheno_pc[c(which(score_pheno_pc$Race == "Other"), which(score_pheno_pc$Race == "Unknown")), ]

#Use standardized adjusted scores to plot 
#Notice: The mean and sd is all_population, not only Other/Unknown
otherUnknown_resid_non_standardized_score <- otherUnknown_score_pheno_pc$scoreResid
otherUnknown_mean_resid <- mean(otherUnknown_resid_non_standardized_score)
otherUnknown_sd_resid <- sd(otherUnknown_resid_non_standardized_score)
otherUnknown_resid_score <- (otherUnknown_resid_non_standardized_score - all_resid_mean) / all_resid_sd  #this is standardized score

#quan=0.02
quan=0.2
score_quantile <- quantile(otherUnknown_resid_score, probs=seq(0, 1, quan), names=T)  #11 numbers with %0=samllest_number and %100=largest+number


quan_cases.df <- data.frame()
for (i in (2:length(score_quantile))){
  if (i==2){
    rows <- which(otherUnknown_resid_score <= score_quantile[i])
    quantile_midpoint <- score_quantile[i-1] + (score_quantile[i] - score_quantile[i-1]) / 2
  }
  else{
    rows <- which(otherUnknown_resid_score > score_quantile[i-1] & otherUnknown_resid_score <= score_quantile[i])
    quantile_midpoint <- score_quantile[i-1] + (score_quantile[i] - score_quantile[i-1]) / 2
  }
  quantile_pheno <- otherUnknown_score_pheno_pc[rows, ]
  subject_num <- nrow(quantile_pheno)
  cases <- sum(quantile_pheno[, 4]==1)
  healthy <- sum(quantile_pheno[, 4]==0)
  
  p <- cases/(healthy + cases)
  odds <- p/(1-p)
  logodds <- log(odds)
  quan_cases.df <- rbind(quan_cases.df, c(quantile_midpoint, cases, subject_num, logodds))
}

#Plot scatter, fit line
quan_cases.df <- quan_cases.df[which(is.finite(quan_cases.df[,4])),]  #Remove Infinite rows
print(nrow(quan_cases.df))
prs_score <- quan_cases.df[,1]
logodds <- quan_cases.df[,4]


#plot(prs_score, odds.adj, xlab="PRS score", ylab=paste0("Odds (x", adj, ")"), xlim=c(-0.6, 0.6), ylim=c(0, 4), main="Other/Unknown AFib PRS Score Analysis")
par(mar=c(3.5,3.5,1,1))
plot(prs_score, logodds, xlab="", ylab="",  xlim=c(-4, 4), ylim=c(-7, 0), main="")
title(xlab="Standardized adjusted PRS", line=2, cex.lab=1.0)
title(ylab="Log-odds", line=2, cex.lab=1.0)

abline(lm(logodds~prs_score), col="red")

corr <- cor(quan_cases.df[,1], quan_cases.df[,4])
text(2.7, -3.5, paste0("r=", round(corr,4)))

text(-3.1, -0.5, paste0("CRCa"))
text(-1.5, -1, paste0("Other/Unknown"))
text(-3.2, -1.6, paste0("q=", quan))








#----------------------------------------------------------------------------------------------------
#---------------------------------- Adjust multiple new sample PRS score ----------------------------
#----------------------------------------------------------------------------------------------------
setwd("~/Documents/PRS_Jason/Cutoff_adjustment/PCA_adjustment_MDS/twenty_two_samples_imputed")

pc_model <- readRDS("colorectalC_pc_model.RDS")
sample_raw_score <- read.table("pre_GDA_test_merged_large_batch_cleaned_imputed_ColorectalCA_22samples.sscore", header=F, stringsAsFactor=F)
output_name <- "pre_GDA_test_merged_large_batch_cleaned_imputed_ColorectalCA_22samples.sscore.adjusted"

sample_raw_pc <- read.table("pre_GDA_test_merged_large_batch_cleaned_imputed_realREF_newid_filtered_new_projection.sscore", header=F, stringsAsFactor=F)

eigenval<- read.table("ref_pcs.eigenval", header=F, stringsAsFactor=F)
sample_scaled_pc <- sample_raw_pc
sample_scaled_pc[,5:8] <- mapply(`/`, sample_raw_pc[,5:8], sqrt(eigenval[,1])) * 2
colnames(sample_scaled_pc) <- c("fid", "iid", "NmissAlleleCT", "NamedAlleleDosageSum", "PC1", "PC2", "PC3", "PC4")

# Predicted sample scores based on the linear model
pred_score <- predict(pc_model, newdata=sample_scaled_pc)

adjusted_score <- sample_raw_score
adjusted_score[,6] <- round((adjusted_score[,6] - pred_score), 7)

colnames(adjusted_score) <- c("#FID", "IID", "NMISS_ALLELE_CT", "NAMED_ALLELE_DOSAGE_SUM", "SCORE1_AVG", "SCORE1_AVG_ADJ")

head(adjusted_score,3)

write.table(adjusted_score, file=output_name, sep="\t", quote=F, row.names=F, col.names=T)





#---------------------------------------------------------------------------------------------------------------
#---------------------- PC-adjusted PRS Score distribution among batches, Age group, Sex -------------------------
#---------------------------------------------------------------------------------------------------------------

#Read in raw PRS batches to obtain batch subjects
setwd("~/Documents/PRS_Jason/prs_analysis/1kg_imputed_normalized/colorectal")
mega_profile <- read.table("mega_colorectal.profile", header=F, stringsAsFactors=F)
megaex_profile <- read.table("megaex_colorectal.profile", header=F, stringsAsFactors=F)
mega1a_profile <- read.table("meg_a1_a_colorectal.profile", header=F, stringsAsFactors=F)
mega1b_profile <- read.table("meg_a1_b_colorectal.profile", header=F, stringsAsFactors=F)
megc_profile <- read.table("meg_c_colorectal.profile", header=F, stringsAsFactors=F)
megd_profile <- read.table("meg_d_colorectal.profile", header=F, stringsAsFactors=F)
mege_profile <- read.table("meg_e_colorectal.profile", header=F, stringsAsFactors=F)
megx1_profile <- read.table("meg_x1_colorectal.profile", header=F, stringsAsFactors=F)


#Read in score_pheno table
score_pheno <- read.table("prs_score_colorectalC_covariates_chip.tsv", header=T, stringsAsFactors=F)

#Create PC-adjusted PRS scores
setwd("~/Documents/PRS_Jason/Cutoff_adjustment/PCA_adjustment_MDS")
pc <- read.table("plink2pca/ref_pcs.eigenvec.newid", header=T, stringsAsFactor=F)

# Merge two table into one
score_pheno_pc <- merge(score_pheno, pc, by=c("FID", "IID"))

# Fit a linear model with control samples
score_pheno_pc_control <- score_pheno_pc[which(score_pheno_pc$COLCA == 0),]

pcmod <- lm(PRS ~ PC1 + PC2 + PC3 + PC4, data=score_pheno_pc_control)
summary(pcmod)

# Predicted all samples scores based on the linear model
pred <- predict(pcmod, newdata=score_pheno_pc)

# Calculate residualized scores and add into score_pheno_pc
score_pheno_pc$scoreResid <- score_pheno_pc$PRS - pred

# Calculate standardized residualized scores and add into score_pheno_pc
score_pheno_pc$StdScoreResid <- (score_pheno_pc$scoreResid - mean(score_pheno_pc$scoreResid)) / sd(score_pheno_pc$scoreResid) 



#------------------- Density Analysis across Batches using standardized adjusted PRS scores -----------------------

#Get standardized PC-adjusted PRS scores for each batch
mega_score <- score_pheno_pc$StdScoreResid[which(score_pheno_pc$IID %in% mega_profile[,2])]
megaex_score <- score_pheno_pc$StdScoreResid[which(score_pheno_pc$IID %in% megaex_profile[,2])]
mega1a_score <- score_pheno_pc$StdScoreResid[which(score_pheno_pc$IID %in% mega1a_profile[,2])]
mega1b_score <- score_pheno_pc$StdScoreResid[which(score_pheno_pc$IID %in% mega1b_profile[,2])]
megc_score <- score_pheno_pc$StdScoreResid[which(score_pheno_pc$IID %in% megc_profile[,2])]
megd_score <- score_pheno_pc$StdScoreResid[which(score_pheno_pc$IID %in% megd_profile[,2])]
mege_score <- score_pheno_pc$StdScoreResid[which(score_pheno_pc$IID %in% mege_profile[,2])]
megx1_score <- score_pheno_pc$StdScoreResid[which(score_pheno_pc$IID %in% megx1_profile[,2])]

all_batch_score <- score_pheno_pc$StdScoreResid


plot(density(mega_score, adjust=2))
plot(density(mega_score, adjust=2), col="black", xlim=c(-4, 4), ylim=c(0, 0.4), xlab="PRS score", lty="dotted", main="Density ~ PRS score of T2D")

par(mar=c(3.5,3.5,1,1))
plot(density(mega_score, adjust=2), col="black", xlim=c(-4, 4), ylim=c(0, 0.4), xlab="",ylab="", lty="dotted", main="")
title(xlab="Standardized adjusted PRS", line=2, cex.lab=1.0)
title(ylab="Density", line=2, cex.lab=1.0)

lines(density(megaex_score, adjust=2), lty="dotted", col="red")
lines(density(mega1a_score, adjust=2), lty="dotted", col="orange")
lines(density(mega1b_score, adjust=2), lty="dotted", col="yellow")
lines(density(megc_score, adjust=2), lty="dotted", col="green")
lines(density(megd_score, adjust=2), lty="dotted", col="cyan")
lines(density(mege_score, adjust=2), lty="dotted", col="blue")
lines(density(megx1_score, adjust=2), lty="dotted", col="purple")

lines(density(all_batch_score, adjust=2), col="darkblue")

#text(x=5e-05, y=15000, labels = c("black: MEGA", "red: MEGAEX", "orange: MEG_A1_A"))
Vars <- c("MEGA", "MEGAEX", "MEG_A1_A", "MEG_A1_B", "MEG_C", "MEG_D", "MEG_E", "MEG_X1", "ALL")  # one per row
colors <- c("black", "red", "orange", "yellow", "green", "cyan", "blue", "purple", "darkblue")  # one color for each row; should be same length as Vars
ltys <- 1  # one lty for each column

nr <- length(Vars) 
nc <- length(ltys)
legend("topright", rep(Vars, nc), col = colors, lty = rep(ltys, each = nr),
       ncol = nc, cex = 0.4)

text(-3.2, 0.35, paste0("CRCa"))

#+++++++++++++ Density analysis of PC-adjusted PRS scores across Gender ++++++++++++++++

male <- score_pheno_pc[score_pheno_pc$Gender=="M" ,]
female <- score_pheno_pc[score_pheno_pc$Gender=="F" ,]
male_score <- male$StdScoreResid
female_score <- female$StdScoreResid
all_score <- score_pheno_pc$StdScoreResid

#plot(density(male_score, adjust=2))
plot(density(male_score, adjust=2), col="red", xlim=c(-4, 4), ylim=c(0, 0.4), xlab="PRS score", lty="solid", main="Density ~ PRS score of T2D")

par(mar=c(3.5,3.5,1,1))
plot(density(male_score, adjust=2), col="red", xlim=c(-4, 4), ylim=c(0, 0.4), xlab="", ylab="", lty="solid", main="")
title(xlab="Standardized adjusted PRS", line=2, cex.lab=1.0)
title(ylab="Density", line=2, cex.lab=1.0)

lines(density(female_score, adjust=2), lty="solid", col="blue")
lines(density(all_score, adjust=2), lty="dotted", col="black")

#text(x=5e-05, y=15000, labels = c("black: MEGA", "red: MEGAEX", "orange: MEG_A1_A"))
Vars <- c("Male", "Female", "ALL")  # one per row
colors <- c("red", "blue", "black")  # one color for each row; should be same length as Vars
ltys <- 1  # one lty for each column

nr <- length(Vars) 
nc <- length(ltys)
legend("topright", rep(Vars, nc), col = colors, lty = rep(ltys, each = nr),
       ncol = nc, cex = 0.5)

text(-3.2, 0.35, paste0("CRCa"))

#One-way ANOVA analysis
anova_1way <- aov(scoreResid~Gender, data=score_pheno_pc)
TukeyHSD(anova_1way)   #Pairwise comparison



#+++++++++++++ Density analysis of PC-adjusted PRS scores across AGE ++++++++++++++++

hist(score_pheno_pc[,16], breaks=10, main="Histogram of Age", xlab="Age")

age10_score <- score_pheno_pc[score_pheno_pc[,6]<=10, 15]   #Only 1 data
age20_score <- score_pheno_pc[score_pheno_pc[,6]>10 & score_pheno_pc[,6]<=20, 16]
age30_score <- score_pheno_pc[score_pheno_pc[,6]>20 & score_pheno_pc[,6]<=30, 16]
age40_score <- score_pheno_pc[score_pheno_pc[,6]>30 & score_pheno_pc[,6]<=40, 16]
age50_score <- score_pheno_pc[score_pheno_pc[,6]>40 & score_pheno_pc[,6]<=50, 16]
age60_score <- score_pheno_pc[score_pheno_pc[,6]>50 & score_pheno_pc[,6]<=60, 16]
age70_score <- score_pheno_pc[score_pheno_pc[,6]>60 & score_pheno_pc[,6]<=70, 16]
age80_score <- score_pheno_pc[score_pheno_pc[,6]>70 & score_pheno_pc[,6]<=80, 16]
age90_score <- score_pheno_pc[score_pheno_pc[,6]>80 & score_pheno_pc[,6]<=90, 16]
age100_score <- score_pheno_pc[score_pheno_pc[,6]>90 & score_pheno_pc[,6]<=100, 16]
age110_score <- score_pheno_pc[score_pheno_pc[,6]>100 & score_pheno_pc[,6]<=110, 16]   #Only 5 data

#age_score <- data.frame(age=c(10,20,30,40,50,60,70,80,90,100), score=c(age10_score, age20_score,age30_score, age40_score,age50_score, age60_score,age70_score, age80_score,age90_score, age100_score))

plot(density(age20_score, adjust=2), col="black", xlim=c(-4, 4), ylim=c(0, 0.4), xlab="PRS score", lty="dotted", main="Density ~ PRS score of T2D")

par(mar=c(3.5,3.5,1,1))
plot(density(age30_score, adjust=2), col="black", xlim=c(-4, 4), ylim=c(0, 0.4), xlab="", ylab="", lty="dotted", main="")
title(xlab="Standardized adjusted PRS", line=2, cex.lab=1.0)
title(ylab="Density", line=2, cex.lab=1.0)

#lines(density(age30_score, adjust=2), lty="dotted", col="red")
lines(density(age40_score, adjust=2), lty="dotted", col="green")
lines(density(age50_score, adjust=2), lty="dotted", col="blue")
lines(density(age60_score, adjust=2), lty="solid", col="red")
lines(density(age70_score, adjust=2), lty="solid", col="yellow")
lines(density(age80_score, adjust=2), lty="solid", col="blue")
lines(density(age90_score, adjust=2), lty="solid", col="cyan")
lines(density(age100_score, adjust=2), lty="solid", col="purple")
#lines(density(age110_score, adjust=2), lty="solid", col="black")  #Only five data

#text(x=5e-05, y=15000, labels = c("black: MEGA", "red: MEGAEX", "orange: MEG_A1_A"))
Vars <- c("Age20_30", "Age30_40", "Age40_50", "Age50_60", "Age60_70", "Age70_80", "Age80_90", "Age90_100")  # one per row
colors <- c("black", "green", "blue", "red", "yellow", "blue", "cyan", "purple")  # one color for each row; should be same length as Vars
ltys <- c("dotted","dotted","dotted","solid","solid","solid","solid","solid")  # one lty for each column

nr <- length(Vars) 
legend("topright", rep(Vars, 1), col = colors, lty = ltys, ncol = 1, cex = 0.4)

text(-3.2, 0.35, paste0("CRCa"))

#One-way ANOVA analysis
age_score <- data.frame(age=c(rep("20",length(age20_score)),rep("30",length(age30_score)),rep("40",length(age40_score)),rep("50",length(age50_score)),rep("60",length(age60_score)),rep("70",length(age70_score)),rep("80",length(age80_score)),rep("90",length(age90_score)),rep("100",length(age100_score))), 
                        score=c(age20_score, age30_score, age40_score,age50_score, age60_score,age70_score, age80_score,age90_score, age100_score))
anova_1way <- aov(score~age, data=age_score)
TukeyHSD(anova_1way)   #Pairwise comparison






