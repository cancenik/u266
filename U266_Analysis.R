library(ggplot2)
library(reshape2)
library(edgeR)
library(limma)

data <- read.table ("~/Documents/U266/Reformatted_Transcript_Counts_All_Libraries.tsv", header=T)
CDS <- split (data, data$REGION)[[2]]
UTR3 <- split (data, data$REGION)[[3]]
UTR5 <- split (data, data$REGION)[[4]]
CDS_Counts <- CDS[,grep("Counts", colnames(CDS))]
CDS_Coverage <- CDS[, grep("CoveredBases", colnames(CDS))]
CDS_Len <- CDS[,grep("Len", colnames(CDS))]
CDS_IDs <- CDS[,1]

UTR5_Counts = UTR5[,grep("Counts", colnames(UTR5))]
UTR3_Counts = UTR3[,grep("Counts", colnames(UTR3))]
UTR5_Len <- CDS[,grep("Len", colnames(UTR5))]
UTR3_Len <- CDS[,grep("Len", colnames(UTR3))]

# For each library, fraction of reads counts mapping to the coding region compared to total
CDS_ratio = colSums(CDS_Counts)  / (colSums(CDS_Counts) + colSums(UTR5_Counts) + colSums(UTR3_Counts) )
# For RFP ~80%
UTR3_ratio = colSums(UTR3_Counts)  / (colSums(CDS_Counts) + colSums(UTR5_Counts) + colSums(UTR3_Counts) )
UTR5_ratio = colSums(UTR5_Counts)  / (colSums(CDS_Counts) + colSums(UTR5_Counts) + colSums(UTR3_Counts) )

# Plot the fractions 
ratio_df = data.frame(library=names(CDS_ratio), 
                      CDS= as.numeric(CDS_ratio), 
                      UTR3 = as.numeric(UTR3_ratio),
                      UTR5= as.numeric(UTR5_ratio) )

ratio_df = melt(ratio_df, id=c("library"))
ggplot(ratio_df, aes(x=library, y=value, fill = variable)) + 
  geom_bar(stat="identity", ) + theme_minimal() + scale_fill_brewer(palette="Blues") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0))

# Explore correlations between RFPs and RNA-Seq
# Check relationship between the technical replicates specifically. 
plot_correlations  = function (lib_i, lib_j, ...) { 
  par(las=1)
  dcols <- densCols(log10(CDS_Counts[, lib_i] +1), log10(CDS_Counts[, lib_j] +1), colramp= colorRampPalette(c("#660099", "#FFCCFF")), nbin = 256)
  plot(log10(CDS_Counts[, lib_i] +1), log10(CDS_Counts[, lib_j] +1), cex=0.2, pch=19, col = dcols, 
       ylab= paste( "log10 of read count" , colnames(CDS_Counts)[lib_i]), 
       xlab= paste( "log10 of read count" , colnames(CDS_Counts)[lib_j]), cex.lab = 0.75)
  text(par("usr")[1]+3, par("usr")[4]-0.5, 
       labels=paste ("rho", 
                     round(cor.test( log10(CDS_Counts[, lib_i] +1), log10(CDS_Counts[, lib_j] +1), 
                                    method = "spearman")$estimate,2) , sep="="))
}

colnames(CDS_Counts)

# Biological Replicate
plot_correlations(1, 10)
# Technical Replicate
plot_correlations(1, 7)
# Different treatment
plot_correlations(1, 4)
# Ribo_RNA
plot_correlations(1, 21)


# Need to define the replicates for the Ribo and RNA. 
# Treat as two separate experiments. 

# RFP 1-3 => Scramble CHX treated
# RFP 4-6 => siBCMA CHX treated

# RFP 1-3 compared to RFP 4-6

# RFP 7-8 => APRIL CHX treated
# APRIL is the positive control -> Ligand to BCMA
# RFP 9-10 => sBCMA + APRIL CHX treated
# soluble BCMA 3-4hr treatment, APRIL 30min treatment

# RFP 7-8 compared to RFP 9-10 but less likely to show a difference.

# Technical Variation Sequencing Lane
# Day1 / Lane 1 -> RFP1, 4, 7, 9
# Day2 / Lane 4 -> RFP 2, 5, 8, 10
# Day2 / Lane 5 -> RFP 1rep, 4rep, 8rep, 3, RNA 5
# Day2 / Lane 6 -> RFP 2rep, 5rep, 9rep, 6, RNA 1
# Day2 / Lane 7 -> RFP 3rep, 6rep, 7rep, 10rep, RNA 2
# Day2 / Lane 8 -> RNA 3, 6, 7, 8, 9, 10

plotMDS(CDS_Counts, labels = c(rep("Ribo", 20), rep("RNA", 10)))
 
plotMDS(CDS_Counts[,1:20], main = "Ribo", 
        labels = c ("Scramble", "siBCMA", "APRIL", "sBCMA_APRIL", "sBCMA_APRIL" , "sBCMA_APRIL" ,
                    "Scramble", "Scramble", "Scramble", "Scramble", "Scramble", 
                    "siBCMA" , "siBCMA", "siBCMA", "siBCMA", "siBCMA", "APRIL", "APRIL", "APRIL" , 
                    "sBCMA_APRIL" ))

plotMDS(CDS_Counts[,21:30], 
         labels = c("Scramble" , "sBCMA_APRIL" ,  "Scramble" , "Scramble", "siBCMA", "siBCMA", "siBCMA", 
                    "APRIL" , "APRIL", "sBCMA_APRIL" ), main = "RNA")

# It's probably fine to merge the technical replicates. 
CDS_Counts_TechMerged = data.frame(
  RFP1_Counts  = CDS_Counts$RFP1_Counts + CDS_Counts$RFP1rep_Counts , 
  RFP2_Counts = CDS_Counts$RFP2_Counts + CDS_Counts$RFP2rep_Counts , 
  RFP3_Counts = CDS_Counts$RFP3_Counts + CDS_Counts$RFP3rep_Counts , 
  RFP4_Counts = CDS_Counts$RFP4_Counts + CDS_Counts$RFP4rep_Counts , 
  RFP5_Counts = CDS_Counts$RFP5_Counts + CDS_Counts$RFP5rep_Counts , 
  RFP6_Counts = CDS_Counts$RFP6_Counts + CDS_Counts$RFP6rep_Counts , 
  RFP7_Counts = CDS_Counts$RFP7_Counts + CDS_Counts$RFP7rep_Counts , 
  RFP8_Counts = CDS_Counts$RFP8_Counts + CDS_Counts$RFP8rep_Counts , 
  RFP9_Counts = CDS_Counts$RFP9_Counts + CDS_Counts$RFP9rep_Counts , 
  RFP10_Counts = CDS_Counts$RFP10_Counts + CDS_Counts$RFP10rep_Counts
)

CDS_Counts_TechMerged = cbind(CDS_Counts_TechMerged, CDS_Counts[,21:30])

plotMDS(CDS_Counts_TechMerged, labels = c(rep("Ribo", 10), rep("RNA", 10)))
plotMDS(CDS_Counts_TechMerged[,1:10] ) 
plotMDS(CDS_Counts_TechMerged[,1:10], main = "Ribo", 
        labels = c ("Scramble", "Scramble", "Scramble",
                    "siBCMA", "siBCMA", "siBCMA",
                    "APRIL", "APRIL", 
                    "sBCMA_APRIL", "sBCMA_APRIL"  ))

colSums(CDS_Counts_TechMerged)
# RFP1_Counts  RFP2_Counts  RFP3_Counts  RFP4_Counts  RFP5_Counts  RFP6_Counts  RFP7_Counts 
# 6244952      5201914      4243600      7941276      7566177      4350595      8578752 
# RFP8_Counts  RFP9_Counts RFP10_Counts  RNA1_Counts RNA10_Counts  RNA2_Counts  RNA3_Counts 
# 4043445      4256164      7983524      3968463      1578939      5226831      5695131 
# RNA4_Counts  RNA5_Counts  RNA6_Counts  RNA7_Counts  RNA8_Counts  RNA9_Counts 
# 3923650      4489038      2415259       610499      1847884      2323720 

## logFC dim 1 corresponds to the sequencing depth. The left cluster is ~4M. 

# Voom/TMM Normalization
joint_counts = DGEList ( counts = CDS_Counts_TechMerged)
isexpr <- rowSums(cpm(joint_counts) > 1) >= 8
joint_counts <- joint_counts[isexpr,]
row.names(joint_counts) <- CDS_IDs[isexpr]
joint_counts <- calcNormFactors (joint_counts, method= "TMM")
# sample_labels <- unlist(lapply(strsplit(colnames(CDS_Counts_TechMerged), "_"), "[[", 1 ) ) 
sample_labels = c( rep("Scramble", 3) ,  rep("siBCMA", 3), rep ("APRIL", 2), rep ("sBCMA_APRIL", 2 ), 
   "Scramble" , "sBCMA_APRIL" ,  "Scramble" , "Scramble", 
   "siBCMA", "siBCMA", "siBCMA", "APRIL" , "APRIL", "sBCMA_APRIL" )        
type <- c(rep("Ribo", 10), rep("RNA", 10) )
full_design <- model.matrix(~sample_labels + type)
v3 <- voom(joint_counts, full_design, plot=T )

joint_dd_rep <- dist(t (v3$E[,1:10]), upper=T, diag= T  )
joint_ribo_rep <- hclust (joint_dd_rep)
plot(joint_ribo_rep, labels= sample_labels[1:10], sub = "", xlab="", main= "Clustering Ribosome Profiling Libraries")

joint_dd_rna <- dist(t (v3$E[,11:20]), upper=T, diag= T)
joint_rna_hc <- hclust(joint_dd_rna)
plot(joint_rna_hc, labels= sample_labels[11:20], sub = "", xlab="", main= "Clustering RNA-Seq Libraries")

#Correlation Coefficients:
cor_coefs = round ( cor (v3$E) , 3 ) 
apply(cor_coefs, 1, median)

cor_coefs[1:10, 1:10]
cor_coefs[11:20, 11:20]

# Differential Expression Analysis and Translation Efficiency
# Two predictors: Sample Label + Ribo vs RNA
# Questions of interest
# What are the genes with differential RNA expression across conditions
# What are the genes with differential Ribo expression across conditions
# What are the genes with differential Ribo-RNA expression across conditions
# The moderated F-statistic can be used as the measure of any difference
# ## EBayes also returns a moderated F-statistic, $F.p.value

ribo_only = v3[, type == "Ribo"]
ribo_only$design = model.matrix(~ 0+ as.factor(sample_labels[type=="Ribo"]))
colnames( ribo_only$design) = c("APRIL",  "sBCMA_APRIL", "Scramble", "siBCMA")
rna_only = v3[, type == "RNA"]
rna_only$design = model.matrix(~ 0+ as.factor(sample_labels[type=="RNA"]))
colnames( rna_only$design) = c("APRIL",  "sBCMA_APRIL", "Scramble", "siBCMA")

contrast_strings <- c()
for (i in 1:length(colnames(rna_only$design))) { 
  contrast_strings[i] <- 
    paste ( 
      paste (colnames(rna_only$design)[i], 
             paste (colnames(rna_only$design)[-i], collapse="+"), sep="- (")
      , length(colnames(rna_only$design))-1 , sep=")/" )
}

# I  will enumerate all strings here for clarity
contrast.matrix<- (makeContrasts(contrast_strings[1], 
                                 contrast_strings[2], 
                                 contrast_strings[3], 
                                 contrast_strings[4], 
                                 levels=rna_only$design))

# Alternative Model is to compare only Sample1-3 to Sample4-6; Sample7-8 to Sample9-10 
contrast_strings2 = c("Scramble- siBCMA", "APRIL- sBCMA_APRIL" )
contrast.matrix2 = makeContrasts(contrast_strings2[1], contrast_strings2[2], levels=rna_only$design)


# # MODEL FITTING 
ribo_fit <- lmFit (ribo_only, design=ribo_only$design, weights=ribo_only$weights)
rna_fit <- lmFit (rna_only, design=rna_only$design, weights=rna_only$weights)

# Compare each treatment to other three. Interesting Results but maybe not important
# ribo_fit2 <- contrasts.fit(ribo_fit, contrast.matrix)
# ribo_fit2 <- eBayes(ribo_fit2)
# rna_fit2 <- contrasts.fit(rna_fit, contrast.matrix)
# rna_fit2 <- eBayes(rna_fit2)
# 
# topTable(ribo_fit2)
# results.ribo <- decideTests(ribo_fit2, p.value=0.01, lfc=log2(1.5))
# apply (results.ribo, 2, function (x) { sum ( x ==1 )} ) 
# results.rna <- decideTests(rna_fit2, p.value=0.01, lfc=log2(1.5))
ribo_fit3 <- contrasts.fit(ribo_fit, contrast.matrix2)
ribo_fit3 <- eBayes(ribo_fit3)
rna_fit3 <- contrasts.fit(rna_fit, contrast.matrix2)
rna_fit3 <- eBayes(rna_fit3)


write.csv (topTable(ribo_fit3, coef = 1, number = 146), 
           file = '~/Desktop/SignificantRibo_ExpressionChange_siRNAtreatment.csv' ) 
results.ribo2 <- decideTests(ribo_fit3, p.value=0.05, lfc=log2(1.5))
# results.ribo2 <- decideTests(ribo_fit3, p.value=0.05, lfc=log2(1))
# We don't have any significant changes in the 1-1.5x change range for Ribo
apply (results.ribo2, 2, function (x) { sum ( x ==1 )} ) 
apply (results.ribo2, 2, function (x) { sum ( x == -1 )} ) 

write.csv (topTable(rna_fit3, coef = 1, number = 225), 
           file = '~/Desktop/SignificantRNA_ExpressionChange_siRNAtreatment.csv' ) 
results.rna2 <- decideTests(rna_fit3, p.value=0.05, lfc=log2(1.5))
#results.rna2 <- decideTests(rna_fit3, p.value=0.05, lfc=log2(1)
# 225 at no fold-change cutoff vs 181 at 1.5
apply (results.rna2, 2, function (x) { sum ( x ==1 )} ) 
apply (results.rna2, 2, function (x) { sum ( x == -1 )} ) 


te_factor <- paste(type, sample_labels, sep=".")
te_factor <- factor(te_factor, levels=unique(te_factor))

te_design = model.matrix(~0 + te_factor)
colnames(te_design) <-  levels (te_factor)
te_fit = lmFit(v3, design = te_design)

cont.matrix.te <- matrix(0,nrow=8, ncol=2, 
                         dimnames=list(Levels=levels(te_factor),
                                       Contrasts=c("siBCMA-Scramble", "APRIL-sBCMA_APRIL") ))

cont.matrix.te[,1] = c(-1,1,0,0,1,0,-1,0)
cont.matrix.te[,2] = c(0,0,1,-1,0,1,0,-1)

te_fit3 <- contrasts.fit(te_fit, cont.matrix.te)
te_fit3 <- eBayes(te_fit3)
te.diff.results <- decideTests (te_fit3, p.value=0.1)
apply (te.diff.results, 2, function (x) { sum ( x == 1 )} ) 
apply (te.diff.results, 2, function (x) { sum ( x == -1 )} ) 

# ATMIN repressed in TE when siBCMA treated
topTable(te_fit3)

# Example Ratio Plot
scramble_rna_mean = apply (v3$E[, c(11, 13, 14)], 1, mean )
si_rna_mean = apply (v3$E[, c(15, 16, 17)], 1, mean )
scramble_ribo_mean = apply (v3$E[, c(1, 2, 3)], 1, mean )
si_ribo_mean = apply (v3$E[, c(4,5,6)], 1, mean )
rna_exp = which (scramble_rna_mean > 4 & si_rna_mean > 4)
# This is an exploratory plot that is not very meaningful statistically.
plot ( si_ribo_mean[rna_exp] /   si_rna_mean[rna_exp], 
       scramble_ribo_mean[rna_exp]/  scramble_rna_mean[rna_exp],
       pch = 19, cex = .25, xlab = "siRNA TE", ylab = "Scramble TE", 
       xlim =c(0,2), ylim = c(0,2) )
abline(a=0, b=1)
## TODO
# Look at other factors in DNA Damage. 
# We should look at the species to total count relationship to remove PCR duplicates
