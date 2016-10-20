print('-------------------------------')
print("//Executing R script for clustering and plotting into a tree")

args <- commandArgs(trailingOnly = TRUE)
cat(args)

suppressPackageStartupMessages(library(dendextend))
library("gplots")
library("ctc")
library("colorspace")
library("dendextend")
library("extrafont")
library("ade4")

currentdir=getwd()
#heatmap3path="/Users/Peter/Programming/python/_github/mimpsearch_hcluster/scripts/heatmap.3.R" 
heatmap3path=args[1]
source(heatmap3path)

infile=args[2]
outputdir=args[3]
#infile = '/Users/Peter/Programming/python/_github/mimpsearch_hcluster/output_15.09.24_16h08m54/without_false_positives_2/03.blastn_presence_absence/blastn_presence_absence.txt' #args[2]
#outputdir = '/Users/Peter/Programming/python/_github/mimpsearch_hcluster/output_15.09.24_16h08m54/without_false_positives_2/04.cluster_and_plot/' #args[3] 
hierclust_plot = "hierclust_plot.pdf"

distance_matrix_rows		= args[4]
clustering_method_rows	= args[5]
distance_matrix_cols		= args[6]
clustering_method_cols	= args[7]

#1 = Jaccard index (1901) S3 coefficient of Gower & Legendre s1 = a / (a+b+c)
#2 = Simple matching coefficient of Sokal & Michener (1958) S4 coefficient of Gower & Legendre s2 = (a+d) / (a+b+c+d)
#3 = Sokal & Sneath(1963) S5 coefficient of Gower & Legendre s3 = a / (a + 2(b + c))
#4 = Rogers & Tanimoto (1960) S6 coefficient of Gower & Legendre s4 = (a + d) / (a + 2(b + c) +d)
#5 = Dice (1945) or Sorensen (1948) S7 coefficient of Gower & Legendre s5 = 2a / (2a + b + c)
#6 = Hamann coefficient S9 index of Gower & Legendre (1986) s6 = (a - (b + c) + d) / (a + b + c + d)
#7 = Ochiai (1957) S12 coefficient of Gower & Legendre s7 = a / sqrt((a + b)(a + c))
#8 = Sokal & Sneath (1963) S13 coefficient of Gower & Legendre s8 = ad / sqrt((a + b)(a + c)(d + b)(d + c))
#9 = Phi of Pearson S14 coefficient of Gower & Legendre s9 = (ad - bc) / sqrt((a + b)(a + c)(d + b)(d + c))
#10 = S2 coefficient of Gower & Legendre S10 = a / (a + b + c + d)

#distance_matrix_rows		= 'pearson'
#clustering_method_rows	= 'average'
#distance_matrix_cols		= 'pearson'
#clustering_method_cols	= 'average'



setwd(outputdir) #this folder should exist!
outfile <- "blastn_presence_absence_reordered.txt"

d <- read.table(infile, sep = "\t", header = TRUE)
row.names(d) <- d[,1] #rename rows to values in first collumn
d[,1] <- NULL #remove the first collumn
title <- paste("Presence of candidate sequences in", nrow(d), "isolates and", ncol(d), "ORFs: \n", distance_matrix_rows, clustering_method_rows, distance_matrix_cols, clustering_method_cols)
data <- as.matrix(d)

#distance    = Dist(data, method = 'pearson')
distance = dist.binary(data, method=distance_matrix_rows, diag = FALSE, upper = FALSE)
cluster     = hclust(distance, method=clustering_method_rows)
dendrogram  = as.dendrogram(cluster) 
Rowv        = rowMeans(data, na.rm = T)
dendrogram  = reorder(dendrogram, Rowv) %>% set("branches_lwd", 2) #%>% set("branches_k_color", k = 10)
#dendrogram <- color_branches(dendrogram, k = 9, col = c("black", "forestgreen", "forestgreen", "blue", "blue", "darkred", "purple", "orange", "limegreen"))

###
#coldistance    = Dist(t(data), method ='pearson')
#print(cor(data))
coldistance = dist.binary(t(data), method=distance_matrix_cols, diag = FALSE, upper = FALSE)
colcluster     = hclust(coldistance, method=clustering_method_cols)
coldendrogram  = as.dendrogram(colcluster) %>% ladderize(right = TRUE) %>% set("branches_lwd", 2)
Colv           = colMeans(data, na.rm = T)
coldendrogram  = reorder(coldendrogram, Colv)
#coldendrogram <- color_branches(coldendrogram, k = 5, col = c("darkblue", "darkred", "darkred", "forestgreen", "black"))


## Re-order the original data using the computed dendrogram
rowInd = rev(order.dendrogram(dendrogram))
colInd = order.dendrogram(coldendrogram)
data_ordered <- data[rowInd, colInd]

#and write to a TXT file:
write.table(data_ordered, outfile, quote=F, sep="\t",row.names=T, col.names=T)

#export tree to newick format
write(hc2Newick(cluster),file="cluster_rows.newick")
write(hc2Newick(colcluster),file="cluster_cols.newick")


#for the data matrix to be plotted ('data'), check if 'SIX' in name. In case this is true, plot that name.
column_annotation = matrix("#dbdbdb", ncol=2, nrow = ncol(data))
sixgenecodes <- list('MAPYSM','MAPYGIV','MKVALV','MQPLRI','MKLSAV','MLVSPI','<<thesearecubense','MAPYSM', 'MKLLWL', 'MFSKAI', 'MTRFHL', 'MHTEYLF', 'MLFKIAW', 'MRFLLLIA', 'MNLKALVV', 'MRFEYI', 'MKLALIA', 'MKYLYLL', 'MDRTHRG', 'MFVSPKA', 'MNLKALVV')
sixgenecodes <- list('SIX')
sixgene_list = list()
for (name in sixgenecodes){
  sixgene_list[length(sixgene_list)+1] <- list(grep(name, colnames(data)))
}
for (sixgene_integer in sixgene_list) {
  column_annotation[sixgene_integer,1] = "red"
}

LScodes <- list('enz_')
LS_list = list()
for (name in LScodes){
  LS_list <- list(grep(name, colnames(data)))
}
for (LS_integer in LS_list) {
  column_annotation[LS_integer,2] = "black"
}

colnames(column_annotation) <- c("SIX genes", "Secreted enzymes")
#rownames(NA/NaN/Inf in foreign function call (arg 11)
rownames(column_annotation) <- colnames(data)

pdf(file=hierclust_plot, width=12, height=16, pointsize = 16, family="Arial")

heatmap.3(data, 
          Rowv=dendrogram,
          Colv=coldendrogram,
          dendrogram="both", 
          #col=colorpanel(10, low="#D0D8EE",high="#233F88"), 
          col=colorpanel(10, low="#dbdbdb",high="#233F88"), 
          key=FALSE, 
          density.info="none", 
          trace="none", 
          labCol=colnames(data),
          cexCol=.4,
          sepcolor="#000000",
          #rowsep=c(0, 9, 12,22, 28, 31, 32, 33, 49, 52, 54, 56, 57, 58, 59, 60, 61),
          rowsep=c(0, 4,7,16,21,25,26,37,38,43,46,47,48,49,50,51,53,56,58,59,60,61),
          #colsep=1:ncol(data),
          #colsep=c(18,28, 45,59, 77),
          main=title,
          cex.main = 1,
          margins=c(9,10),
          #RowSideColors= myClusterSideBar
          ColSideColors = column_annotation
)
dev.off()

print('Finished R script for clustering and plotting..')
