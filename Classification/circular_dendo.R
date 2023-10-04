library(dendextend)
library(readr)
corRaw <- read_csv("Research/Protien_Database/extracted_new_samples/testing/sample_CATH_1/theta29_dist35/corrected_normal_jaccard_similarity_theta29_dist35_NoFeatureSelection_keyCombine0.csv")
distance <- as.dist(corRaw)
hc_iris <- hclust(distance, method = "complete")
dend <- as.dendrogram(hc_iris)
dend <- rotate(dend, 1:386)#386 is no of proteins in the sample
dend <- color_branches(dend, k=4)
dend <- hang.dendrogram(dend,hang_height=0.1)
dend <- set(dend, "labels_cex", 0.56)

plot(hclust(distance), 
     main="Dissimilarity = 1 - Correlation", xlab="")

plot(dend, 
     main = "Clustered protein CATH data set", 
     horiz =  TRUE,  nodePar = list(cex = .007))
par(mar = rep(0,4))
circlize_dendrogram(dend)
