####################################PACKAGES AND DIRECTORY###################################################
#Set working directory
setwd("C:/Users/Administrator/Downloads/FMT/ANALISI/NEW/Github") #Change me!

#Load packages
library(phyloseq)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(MicrobiotaProcess)
library(ggpubr)

library(readxl)
library(reshape2)
library(tidyr)
library(rstatix)
library(vegan)
library(patchwork)
library(DESeq2)
library(tibble)
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(genefilter)

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
sample_data(ps)$Treatment <-factor(sample_data(ps)$Treatment, levels = c("CON","FFT","FMT","GMT"))

#####################################SUBSET and FILTER###################################################
#Colon samples
ps_CC = subset_samples(ps,  Type == "Colon content")
ps_CC = prune_taxa(taxa_sums(ps_CC) > 0, ps_CC)
ps_CC = filter_taxa(ps_CC, genefilter::filterfun(genefilter::kOverA(10, 10)), TRUE)
ps_CC

#Inoculum samples
ps_Pool = subset_samples(ps,  Type != "Colon content")
ps_Pool = prune_taxa(taxa_sums(ps_Pool) > 0, ps_Pool)
ps_Pool = filter_taxa(ps_Pool, genefilter::filterfun(genefilter::kOverA(2, 100)), TRUE)
ps_Pool

####################################Taxonomic profile of inoculum samples###################################################
#PS for analysis
ps_analysis <- ps_Pool #Change to ps ps_CC ps_Pool
ps_analysis

#Color palette
color_palette = brewer.pal(12, "Paired")

#Phylum
phylum <- get_taxadf(obj=ps_analysis, taxlevel=2)
color_palette[4] <- "#969696"
phylum_graph <- ggbartax(obj=phylum, facetNames="Treatment", topn=3, plotgroup=TRUE, count=TRUE) +
  xlab(NULL) + ylab("Reads (no.)") + labs(fill = "Phylum") +
  scale_fill_manual(values = color_palette) + geom_bar(colour="black", stat="identity", size=0.1) +
  #labs(title = "Colon digesta") +
  cowplot::theme_cowplot() #+scale_x_discrete(labels = SiteMOD) #+ 
#theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
#phylum_graph

#Family
family <- get_taxadf(obj=ps_analysis, taxlevel=5)
color_palette[4] <- "#33a02c"
color_palette[12] <- "#969696"
family_graph <- ggbartax(obj=family, facetNames="Treatment", topn=11, plotgroup=TRUE, count=TRUE) +
  xlab(NULL) + ylab("Reads (no.)") + labs(fill = "Family") +
  scale_fill_manual(values = color_palette) + geom_bar(colour="black", stat="identity", size=0.1) +
  cowplot::theme_cowplot() #+scale_x_discrete(labels = Treatment_MOD) #+
#theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
#family_graph

#Genus
genus <- get_taxadf(obj=ps_analysis, taxlevel=6)
color_palette[4] <- "#33a02c"
color_palette[12] <- "#969696"
genus_graph <- ggbartax(obj=genus, facetNames="Treatment", topn=11, plotgroup=TRUE, count=TRUE) +
  xlab(NULL) + ylab("Reads (no.)") + labs(fill = "Genus") +
  scale_fill_manual(values = color_palette) + geom_bar(colour="black", stat="identity", size=0.1) +
  cowplot::theme_cowplot() #+ scale_x_discrete(labels = Treatment_MOD) #+ 
#theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
#genus_graph

#Save in pdf
pdf("Figure_3.pdf")
ggarrange(family_graph,genus_graph,
          nrow=2,ncol=1,widths=c(1,1),
          labels=c("A","B"),font.label = list(size = 12, color = "black"))
dev.off()

###################################Calculate relative abundance##################################################
#PS for analysis
ps_analysis <- ps_Pool #Change to ps ps_CC ps_Pool
ps_analysis

ps_merge <- merge_samples(ps_analysis, "Treatment")
#ps_merge = transform_sample_counts(ps_merge, function(x) {x/sum(x)*100})
taxa <- get_taxadf(obj=ps_merge, taxlevel=6)
taxa_melt = psmelt(taxa)

#Write csv for analysis
write.csv(taxa_melt, "melt.csv")
