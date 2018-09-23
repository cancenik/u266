library(ggplot2)
library(reshape2)

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



## TODO
# We should look at the species to total count relationship to remove PCR duplicates
