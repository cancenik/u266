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
plot_correlations  = function (lib_i, lib_j) { 
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

# It's probably fine to merge the technical replicates. 

# Need to define the replicates for the Ribo and RNA. 
# RFP 1-3 => Scramble CHX treated
# RFP 4-6 => siBCMA CHX treated
# RFP 7-8 => APRIL CHX treated
# RFP 9-10 => siBCMA + APRIL CHX treated

# Technical Variation Sequencing Lane
# Day1 / Lane 1 -> RFP1, 4, 7, 9
# Day2 / Lane 4 -> RFP 2, 5, 8, 10
# Day2 / Lane 5 -> RFP 1rep, 4rep, 8rep, 3, RNA 5
# Day2 / Lane 6 -> RFP 2rep, 5rep, 9rep, 6, RNA 1
# Day2 / Lane 7 -> RFP 3rep, 6rep, 7rep, 10rep, RNA 2
# Day2 / Lane 8 -> RNA 3, 6, 7, 8, 9, 10

plotMDS(CDS_Counts, labels = c(rep("Ribo", 20), rep("RNA", 10)))
 
plotMDS(CDS_Counts[,1:20], main = "Ribo", 
        labels = c ("Scramble", "siBCMA", "APRIL", "si+APRIL", "si+APRIL" , "si+APRIL" ,
                    "Scramble", "Scramble", "Scramble", "Scramble", "Scramble", 
                    "siBCMA" , "siBCMA", "siBCMA", "siBCMA", "siBCMA", "APRIL", "APRIL", "APRIL" , 
                    "si+APRIL" ))

plotMDS(CDS_Counts[,21:30], 
         labels = c("Scramble" , "si+APRIL" ,  "Scramble" , "Scramble", "siBCMA", "siBCMA", "siBCMA", 
                    "APRIL" , "APRIL", "si+APRIL" ), main = "RNA")

# Voom/TMM Normalization
joint_counts = DGEList ( counts = CDS_Counts)
isexpr <- rowSums(cpm(joint_counts) > 1) >= 12
joint_counts <- joint_counts[isexpr,]
row.names(joint_counts) <- CDS_IDs[isexpr]
joint_counts <- calcNormFactors (joint_counts, method= "TMM")
sample_labels <- unlist(lapply(strsplit(colnames(CDS_Counts), "_"), "[[", 1 ) ) 
type <- c(rep("Ribo", 20), rep("RNA", 10) )
full_design <- model.matrix(~sample_labels + type)

v3 <- voom(joint_counts, full_design, plot=T )

## TODO
# We should look at the species to total count relationship to remove PCR duplicates
# Add Simple_Utils.pm, merge_counts and alignment strategy to GitHub
