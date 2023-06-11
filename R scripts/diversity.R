####################################PACKAGES AND DIRECTORY###################################################
#Set working directory
setwd("C:/Users/Administrator/Downloads/FMT") #Change me!

#Load packages
library(phyloseq)
library(dplyr)
library(ggplot2)
library(MicrobiotaProcess)
library(tibble)
library(reshape2)
library(rstatix)
library(ggpubr)
library(vegan)

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
ps_CC

#Inoculum samples
ps_Pool = subset_samples(ps,  Type != "Colon content")
ps_Pool = prune_taxa(taxa_sums(ps_Pool) > 0, ps_Pool)
ps_Pool

#####################################ALPHA DIVERSITY###################################################
#ps for analysis
ps_analysis <- ps_CC #Change to ps ps_CC ps_Pool
ps_analysis

#Rarefy with depth of 90% of the minimum sample depth in the dataset
ps_analysis = phyloseq::rarefy_even_depth(ps_analysis, rngseed=1, sample.size=0.9*min(sample_sums(ps_analysis)), replace=F)
ps_analysis

#Calculate alpha-diversity indexes
rich <- as.data.frame(get_alphaindex(ps_analysis))
alpha <- subset(rich, select = c("Shannon","Chao1"))
alpha$Treatment <- sample_data(ps_analysis)$Treatment

#Melt alpha
alpha_long <- melt(alpha,id.vars = c("Treatment"),
                 variable.name = "Parameter", value.name = "Index")

#Subset indexes
shannon = subset(alpha_long, Parameter == "Shannon")
chao1 = subset(alpha_long, Parameter == "Chao1")

#Multiple comparison of alpha-diversity indexes using the Wilcoxon Rank Sum test  
shannon_stat <- shannon %>%
  group_by(Parameter) %>%
  wilcox_test(Index~Treatment) %>%
  adjust_pvalue(method = "BH")%>%
  add_significance(p.col = "p") %>%
  filter(p.adj<=0.05) %>%
  mutate(y.position=c(5))
shannon_stat = subset(shannon_stat, group1 == "CON")

chao1_stat <- chao1 %>%
  group_by(Parameter) %>%
  wilcox_test(Index~Treatment) %>%
  adjust_pvalue(method = "BH")%>%
  add_significance(p.col = "p") %>%
  filter(p.adj<=0.05) %>%
  mutate(y.position=c(5))
chao1_stat = subset(chao1_stat, group1 == "CON")

shannon_plot <- ggplot(shannon, aes(x=Treatment, y=Index)) +
  stat_boxplot(geom ='errorbar', linetype=1, width=0.5) + 
  geom_boxplot(outlier.shape = NA, aes(fill=Treatment)) +
  geom_point(position = position_jitter(w=0.1, h=0),size=0.5) + #Add color if needed aes(color=Data)
  #scale_color_manual(values = c("black", "gray")) + #Changes color of points when aes(color=Data)
  scale_fill_manual(values = cols) +
  labs(x="", y="Shannon index") +
  stat_pvalue_manual(shannon_stat,label = "p.signif", tip.length = 0.005, size = 6, y.position = 4.6)+ 
  #labs(tag = "a") + theme(plot.tag.position = "topleft") +
  #labs(title="Ileum digesta") +
  theme_classic() +
  guides(fill = "none") +
  mytheme 
#shannon_plot

chao1_plot <- ggplot(chao1, aes(x=Treatment, y= Index)) +
  stat_boxplot(geom ='errorbar', linetype=1, width=0.5) + 
  geom_boxplot(outlier.shape = NA, aes(fill=Treatment)) +
  geom_point(position = position_jitter(w=0.1, h=0),size=0.5) + #Add color if needed aes(color=Data)
  #scale_color_manual(values = c("black", "gray")) + #Changes color of points when aes(color=Data)
  scale_fill_manual(values = cols) +
  labs(x="", y="Chao1 index") +
  stat_pvalue_manual(chao1_stat,label = "p.signif", tip.length = 0.005, size = 6, y.position = 400)+ 
  theme_classic() +
  guides(fill = "none") +
  #labs(tag = "b") + theme(plot.tag.position = "topleft") +
  #labs(title="Ileum digesta") +
  mytheme 
#chao1_plot

#Save in pdf
pdf("alpha-diversity.pdf")
ggarrange(shannon_plot,chao1_plot,
          nrow=1,ncol=2,widths=c(1,1),
          labels=c("A","B"),font.label = list(size = 12, color = "black"))
dev.off()

#####################################BETA DIVERSITY###################################################
#ps for analysis
ps_analysis <- ps_CC #Change to ps ps_CC ps_Pool
ps_analysis

#Rarefy with depth of 90% of the minimum sample depth in the dataset
ps_analysis = phyloseq::rarefy_even_depth(ps_analysis, rngseed=1, sample.size=0.9*min(sample_sums(ps_analysis)), replace=F)
ps_analysis

bray <- phyloseq::distance(ps_analysis, method = "bray")

samples_df <- data.frame(sample_data(ps_analysis))
otu_df <- as.data.frame(t(otu_table(ps_analysis)))
bray <- as.matrix(bray)
bray <- bray[rownames(samples_df),rownames(samples_df)]
pcoa = cmdscale(bray, k=3, eig=T)
pcoa_points = as.data.frame(pcoa$points)
colnames(pcoa_points) = c("x", "y", "z")
pcoa_points$Treatment <- samples_df$Treatment
eig = pcoa$eig

#CHECK BETA-DISPERSION
betadisper_bray = betadisper(d = phyloseq::distance(ps_analysis, method = "bray"), group = sample_data(ps_analysis)$Treatment)
anova(betadisper_bray)
TukeyHSD(betadisper_bray)

#PAIRWISE ADONIS
pairwise_adonis = pairwiseAdonis::pairwise.adonis(otu_df, factors=samples_df$Treatment, sim.method="bray", p.adjust.m = "BH", perm = 999)
pairwise_adonis_text = paste("PERMANOVA,",pairwise_adonis$pairs,", p-adj =", format(pairwise_adonis$p.adjusted, digits = 3))
pairwise_adonis_text

points_plot = ggplot(pcoa_points, aes(x=x, y=z, color=Treatment)) +
  geom_point(size=2) + #Add aes(shape=Data) if needed
  stat_ellipse(geom = "polygon",linetype = 2, aes(fill=Treatment), alpha=0.05) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=3), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=3), "%)", sep=""))+
  #labs(title="Beta-diversity") +
  annotate("text", label = pairwise_adonis_text[6], x=-0.45, y=-0.30) +
  annotate("text", label = pairwise_adonis_text[5], x=-0.45, y=-0.35) +
  annotate("text", label = pairwise_adonis_text[3], x=-0.45, y=-0.40) +
  scale_color_manual(values = cols) +
  theme_classic()+
  mytheme
points_plot

#Save in pdf
pdf("beta-diversity.pdf")
ggarrange(points_plot,
          nrow=1,ncol=1,widths=c(1),
          labels=c("C"),font.label = list(size = 12, color = "black"))
dev.off()
