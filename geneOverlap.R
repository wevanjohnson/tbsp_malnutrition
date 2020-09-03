# Code used to produce gene overlap figure and additional statistics
# Witten by David Jenkins, Aubrey R. Odom

library(TBSignatureProfiler)
length(unique(unlist(TBsignatures)))
sort(table(unlist(TBsignatures)), decreasing = TRUE)

#without the 900 gene signature
a <- TBsignatures
a["Esmail_TB_vs_LTBI_893"] <- NULL
length(unlist(a))

# count how many genes are unique or in multiple sigs

# marginal
altable_less <- table(as.matrix(sort(table(unlist(TBsignatures)), decreasing = TRUE)))
altable_less
round(altable/length(unique(unlist(TBsignatures))) * 100, 1)

# cumulative
altable <- cumsum(table(as.matrix(sort(table(unlist(TBsignatures)), 
                                       decreasing = TRUE))))
altable
round(altable/length(unique(unlist(TBsignatures))) * 100, 1)

genesin5ormore <- names(sort(table(unlist(TBsignatures)), decreasing = TRUE)[sort(table(unlist(TBsignatures)), decreasing = TRUE) > 4])

library(ComplexHeatmap)
geneoverlapdf <- data.frame(row.names = genesin5ormore)
for(i in names(TBsignatures)){
  geneoverlapdf[,i] <- ifelse(rownames(geneoverlapdf)%in% TBsignatures[[i]], "Yes", "No")
}
geneoverlapdf <- geneoverlapdf[, order(colSums(geneoverlapdf == "Yes"), decreasing = TRUE)]
geneoverlapdf <- geneoverlapdf[order(rowSums(geneoverlapdf == "Yes"), decreasing = TRUE), ]
draw(Heatmap(t(as.matrix(geneoverlapdf)), col=c("white", "blue"), row_names_gp = gpar(fontsize=8),
        column_names_gp = gpar(fontsize = 8), rect_gp = gpar(col = "grey", lty = 1),
        name="Gene in Signature"), heatmap_legend_side = "bottom") 
