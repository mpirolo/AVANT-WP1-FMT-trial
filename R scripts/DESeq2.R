####################################PACKAGES AND DIRECTORY###################################################
#Set working directory
setwd("C:/Users/Administrator/Downloads/FMT") #Change me!

#Load packages
library(phyloseq)
library(dplyr)
library(ggplot2)
library(DESeq2)
library(tidyr)
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

#####################################SET THEME AND COLORS###################################################
mytheme<- theme(plot.title = element_text(hjust=0.5, family = "Arial", size=12),
                legend.position ="right",
                legend.text = element_text(family = "Arial", size = 10),
                legend.background = element_blank(),
                strip.background = element_blank(),
                strip.placement = "outside",
                strip.text = element_text(family = "Arial", size = 10),
                axis.text = element_text(family = 'Arial', size = 10, color="black"),
                panel.grid = element_blank())

#Set colors
cols = c("CON"="#af8dc3","FFT"="#91bfdb","FMT"="#fc8d59","GMT"="#7fbf7b")

#####################################IMPORT PHYLOSEQ###################################################
ps <- readRDS("ps_FMT.rds")
ps

#Update sample data
samples_df <- as.data.frame(read.delim("metadata.tsv", header=TRUE))
samples_df <- samples_df %>% tibble::column_to_rownames("ID") 
sample_data(ps) = sample_data(samples_df)
sample_data(ps)

#####################################SUBSET and FILTER###################################################
#Colon samples
ps_CC = subset_samples(ps,  Type == "Colon content")
ps_CC = prune_taxa(taxa_sums(ps_CC) > 0, ps_CC)
ps_CC = filter_taxa(ps_CC, genefilter::filterfun(genefilter::kOverA(10, 10)), TRUE)
ps_CC

#Inoculum samples
ps_Pool = subset_samples(ps,  Type != "Colon content")
ps_Pool = prune_taxa(taxa_sums(ps_Pool) > 0, ps_Pool)
ps_CC = filter_taxa(ps_CC, genefilter::filterfun(genefilter::kOverA(10, 10)), TRUE)
ps_Pool

####################################DESEQ2###################################################
#ps for analysis
ps_analysis <- ps_CC #Change to ps ps_CC ps_Pool
ps_analysis

#Groups for analysis
Treatments = c("FFT","FMT","GMT")

DESEQ_list = list()

#DESEQ2 for loop
for (group in Treatments){
  ps_prune <- prune_samples(sample_data(ps_analysis)$Treatment %in% c("CON",group), ps_analysis)
  otu_table(ps_prune) <- otu_table(ps_prune) + 1
  deseq <- phyloseq_to_deseq2(ps_prune, ~Treatment)
  cts <- counts(deseq)
  geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(sum(log(row[row != 0]))/length(row)))
  dds <- estimateSizeFactors(deseq, geoMeans=geoMeans)
  deseq <-  DESeq(dds, test="Wald", fitType="parametric")
  deseq_res = results(deseq, cooksCutoff = FALSE)
  sigtab = deseq_res
  sigtab = cbind(as(sigtab, "data.frame"), as(phyloseq::tax_table(ps_prune)[rownames(sigtab), ], "matrix"))
  head(sigtab)
  x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
  x = sort(x, TRUE)
  sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
  x = tapply(sigtab$log2FoldChange, sigtab$Species, function(x) max(x))
  x = sort(x, TRUE)
  sigtab$Species = factor(as.character(sigtab$Species), levels=names(x))
  DESEQ_list[[group]] = data.frame(sigtab)
  #write.table(data.frame(sigtab), paste0(group,".tsv"), sep="\t", col.names = NA)
}

#Write results as tab
tab_FFT <- DESEQ_list$FFT
tab_FMT <- DESEQ_list$FMT
tab_GMT <- DESEQ_list$GMT

#Filter results
tab_FFT <- subset(tab_FFT, padj <= 0.05 & (log2FoldChange >=5))
tab_FMT <-subset(tab_FMT, padj <= 0.05 & (log2FoldChange >=5))
tab_GMT <-subset(tab_GMT, padj <= 0.05 & (log2FoldChange >=5))

#Transform and subset ps based on significant ASVs
ASV <- unique(c(rownames(tab_FFT), rownames(tab_FMT), rownames(tab_GMT)))
ps_analysis_log <- transform_sample_counts(ps_analysis, function(x) log(1 + x)) #x/sum(x)*100 or log(1 + x)
ps_analysis_log_sign <- prune_taxa(rownames(otu_table(ps_analysis_log)) %in% ASV,ps_analysis_log)

#Data frames for analysis
otu_df <- data.frame(phyloseq::otu_table(ps_analysis_log_sign), check.names = F)
tax_df <- data.frame(phyloseq::tax_table(ps_analysis_log_sign))
sample_df <- data.frame(phyloseq::sample_data(ps_analysis_log_sign))

#Reorder dataframes
sample_df <- sample_df[order(sample_df$Treatment),]
rownames(sample_df)
otu_df <- t(otu_df)
otu_df <- otu_df[match(rownames(sample_df), rownames(otu_df)), ]

#Modify the tax (ASV) dataframe
index <- match(rownames(tax_df), rownames(tab_FFT))
tax_df$FFT_log <- tab_FFT$log2FoldChange[index]
tax_df$FFT_p <- tab_FFT$padj[index]
index <- match(rownames(tax_df), rownames(tab_FMT))
tax_df$FMT_log <- tab_FMT$log2FoldChange[index]
tax_df$FMT_p <- tab_FMT$padj[index]
index <- match(rownames(tax_df), rownames(tab_GMT))
tax_df$GMT_log <- tab_GMT$log2FoldChange[index]
tax_df$GMT_p <- tab_GMT$padj[index]
tax_df <- tax_df[order(tax_df$FFT_p,tax_df$FMT_p,tax_df$GMT_p),]
tax_df <- tax_df[order(tax_df$Family),]
tax_df$Name <- rownames(tax_df)
tax_df$Species <- tax_df$Species %>% str_replace("NA","spp.")
tax_df$Species_complete <- paste(tax_df$Genus, tax_df$Species, sep = " ")
tax_df$Species_ASV <- paste(tax_df$Name, tax_df$Species_complete, sep = " - ")
tax_df$Genus <- tax_df$Genus %>% str_replace("NA","unclassified")
tax_df$Genus_complete <- ifelse(tax_df$Genus == "unclassified", paste(tax_df$Family, tax_df$Genus, sep =" "), paste(tax_df$Genus))
tax_df$FFT_log <- tax_df$FFT_log %>% replace_na(1)
tax_df$FFT_p <- tax_df$FFT_p %>% replace_na(1)
tax_df$FMT_log <- tax_df$FMT_log %>% replace_na(1)
tax_df$FMT_p <- tax_df$FMT_p %>% replace_na(1)
tax_df$GMT_log <- tax_df$GMT_log %>% replace_na(1)
tax_df$GMT_p <- tax_df$GMT_p %>% replace_na(1)
tax_df$FFT_associated <- ifelse(tax_df$FFT_log > 1, "*", "ns")
tax_df$FMT_associated <- ifelse(tax_df$FMT_log > 1, "*", "ns")
tax_df$GMT_associated <- ifelse(tax_df$GMT_log > 1, "*", "ns")

#Reorder dataframes
otu_df <- t(otu_df)
otu_df <- otu_df[match(rownames(tax_df), rownames(otu_df)), ]
otu_df

#Calculate mean, max and min of log(1+abundance) and log2foldchange
mean(otu_df)
max(otu_df)
min(otu_df)

mean(c(tax_df$FFT_log,tax_df$FMT_log,tax_df$GMT_log))
max(c(tax_df$FFT_log,tax_df$FMT_log,tax_df$GMT_log))
min(c(tax_df$FFT_log,tax_df$FMT_log,tax_df$GMT_log))

#Colors
color_log2Fold = colorRamp2(c(1, 5), c("#f9f9f9", "#9e9ac8"))
color_logscore = colorRamp2(c(0, 3, 6, 9), c("#fef0d9","#fdcc8a","#fc8d59","#d7301f"))
color_taxa <- brewer.pal(11, "Paired")

#Genus and family colors
genus <- unique(as.character(tax_df$Genus))
genus_col <- colorRampPalette(color_taxa)(length(genus))
names(genus_col) <- genus
family <- unique(as.character(tax_df$Family))
family_col <- colorRampPalette(color_taxa)(length(family))
names(family_col) <- family

#Annotations
ha_row <- HeatmapAnnotation(
  FFT=anno_simple(tax_df$FFT_log, col = color_log2Fold, pch=na_if(tax_df$FFT_associated,"ns")),
  FMT=anno_simple(tax_df$FMT_log, col = color_log2Fold, pch=na_if(tax_df$FMT_associated,"ns")),
  GMT=anno_simple(tax_df$GMT_log, col = color_log2Fold, pch=na_if(tax_df$GMT_associated,"ns")),
  which = "row")

ha_row_txt <- rowAnnotation(Family=anno_simple(tax_df$Family, col = family_col),
                            labels = anno_text(tax_df$Species_ASV, which = "row"))

ha_col = HeatmapAnnotation(Treatment=anno_simple(sample_df$Treatment,
                                                 col = c("CON"="#af8dc3","FFT"="#91bfdb","FMT"="#fc8d59","GMT"="#7fbf7b")))

#Histogram
Hist <- Heatmap(otu_df, cluster_columns = FALSE, cluster_rows = FALSE,
                name="Log-score", 
                col=color_logscore,
                top_annotation = ha_col, 
                left_annotation = ha_row, 
                right_annotation = ha_row_txt,
                show_row_names = FALSE, 
                show_column_names = FALSE,
                heatmap_legend_param = list(legend_direction = "horizontal"),
                show_heatmap_legend = FALSE)
Hist

#Add legends
lgd_Treatment = Legend(title = "Treatment", legend_gp = gpar(fill = c("#af8dc3","#91bfdb","#fc8d59","#7fbf7b")),
                    labels = c("CON","FFT","FMT","GMT"))
lgd_log2Fold = Legend(title = "Log2 Fold Change", col_fun = color_log2Fold, at = c(1, 5),
                      labels = c("NS", "*"), direction = "horizontal")
lgd_logscore = Legend(title = "Log(1 + abundance)", col_fun = color_logscore, at = c(0, 3, 6, 9),
                      labels = c("0","3","6",">9"), direction = "horizontal")
lgd_family = Legend(title = "Family", legend_gp = gpar(fill = family_col),labels = family, ncol = 3)

#Save histogram
pdf("DESEQ2_hist.pdf",height = 8, width = 10)
draw(Hist, heatmap_legend_list=list(lgd_Treatment,lgd_log2Fold,lgd_logscore,lgd_family), ht_gap =,
     heatmap_legend_side = "bottom")
dev.off()
