####################################PACKAGES AND DIRECTORY###################################################
#Set working directory
setwd("C:/Users/qmn288/Downloads/THOMAS/FINAL") #Change me!

#Load packages
library(ggplot2)
library(dplyr)
library(readxl)
library(tidyr)
library(rstatix)
library(ggpubr)
library(cowplot)
library(MicrobiotaProcess)
library(vegan)
library(reshape2)
library(phyloseq)
library(rlang)
library(metagenomeSeq)
library(microbiomeMarker)

#####################################IMPORT PHYLOSEQ###################################################
ps_MAG <- readRDS("ps_MAG.rds")
ps_16S <- readRDS("ps_16S.rds")

#Normalized 16S data using CSS method#Normalps_16S_preized 16S data using CSS method
# ps_16S_metaseq <- phyloseq_to_metagenomeSeq(ps_16S)
# ps_16S_metaseq<- cumNorm(ps_16S_metaseq, p=cumNormStat(ps_16S_metaseq))
# 
# ps_16S_MRcounts <- MRcounts(ps_16S_metaseq, norm = TRUE)
# ps_16S_CSS <- merge_phyloseq(otu_table(ps_16S_MRcounts, taxa_are_rows = T),
#                                           sample_data(ps_16S),
#                                           phyloseq::tax_table(ps_16S))

#Set colors
cols = c("CON"="#91bfdb","FFT"="#fc8d59")

#####################################SUBSET PHYLOSEQ###################################################
ps_16S_final <- subset_samples(ps_16S, Original != "Filtrate")
ps_16S_final = prune_taxa(taxa_sums(ps_16S_final) > 0, ps_16S_final)
ps_16S_final

ps_16S_filtrate <- subset_samples(ps_16S, Original == "Filtrate")
ps_16S_filtrate = prune_taxa(taxa_sums(ps_16S_filtrate) > 0, ps_16S_filtrate)
ps_16S_filtrate

#Pre-weaning
ps_16S_pre = subset_samples(ps_16S_final,  Weaning == "Pre-weaning")
ps_16S_pre = prune_taxa(taxa_sums(ps_16S_pre) > 0, ps_16S_pre)
ps_16S_pre = filter_taxa(ps_16S_pre, genefilter::filterfun(genefilter::kOverA(2, 20)), TRUE)
ps_16S_pre

ps_MAG_pre = subset_samples(ps_MAG,  Weaning == "Pre-weaning")
ps_MAG_pre = prune_taxa(taxa_sums(ps_MAG_pre) > 0, ps_MAG_pre)
ps_MAG_pre

#Post-weaning
ps_16S_post = subset_samples(ps_16S_final, Weaning == "Post-weaning")
ps_16S_post = prune_taxa(taxa_sums(ps_16S_post) > 0, ps_16S_post)
ps_16S_post = filter_taxa(ps_16S_post, genefilter::filterfun(genefilter::kOverA(5, 50)), TRUE)
ps_16S_post

ps_MAG_post = subset_samples(ps_MAG,  Weaning == "Post-weaning")
ps_MAG_post = prune_taxa(taxa_sums(ps_MAG_post) > 0, ps_MAG_post)
ps_MAG_post

#####################################TAXONOMICAL ANALYSIS###################################################
#PS for analysis
ps_analysis = ps_MAG_post #Change to ps_16S_final ps_16S_pre ps_16S_post ps_MAG ps_MAG_pre ps_MAG_post

#Merge taxa ra#Merge taxa ra#Merge taxa rank, transform to RA and melt ps
taxa <- ps_analysis %>%
  tax_glom("Family") %>% #Change taxa rank
  phyloseq::transform_sample_counts(function(x) { x/sum(x) }) %>%
  psmelt()

#Dataframe for analyis
taxa_df <- data.frame(
  MAG = taxa$OTU,
  ID = taxa$Sample, 
  Day = taxa$Day,
  Treatment = taxa$Treatment,
  Taxa = taxa$Family, #change taxa rank
  Abundance = taxa$Abundance)

#Group and calculate the mean relative abundance
taxa_df = taxa_df %>%
  group_by(Taxa,Treatment) %>% #Add Day if needed 
  summarise(mean_RA = mean(Abundance))

#Assign low abundance taxa as "Others"
taxa_df$Taxa[taxa_df$mean_RA < 0.005] <- "Others"
taxa_df$Taxa <- factor(taxa_df$Taxa, levels = c(setdiff(unique(taxa_df$Taxa), "Others"), "Others"))

#Reorder Treatment
taxa_df$Treatment <-factor(taxa_df$Treatment, 
                           levels = c("CON","FFT"))

write.csv(taxa_df, "taxa_MAG_post.csv")

#Get colors
#cols_taxa <- ggsci::pal_igv(palette = "default")(45)
cols_taxa <- ggthemes::tableau_color_pal(palette = "Tableau 20")(20)
cols_taxa[20] <- "#D3D3D3"

taxa_plot <- ggplot(taxa_df, aes(x = Treatment, y = mean_RA, fill = forcats::fct_reorder(Taxa, mean_RA, .desc = TRUE))) +
  geom_col(position = "fill") +
  #ggforce::facet_row(~ Day, space = "free", scales = "free_x") +
  scale_fill_manual(name="Genus",
                    values= cols_taxa,
                    guide = guide_legend(ncol = 1, reverse = FALSE)) +
  scale_y_continuous(
    labels = scales::percent,
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    #breaks = scales::pretty_breaks(n=4),
    expand = expansion(mult = 0)) +
  labs( x = NULL, #Change to Day if needed
        y = "Relative abundance") +
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(family = 'Arial', size = 10, color = "black"),
        strip.text = element_text(family = "Arial", size = 10, color = "black"),
        legend.text = element_text(family = "Arial", size = 10, color = "black"),
        legend.title = element_text(family = "Arial", size = 10, color = "black"),
        legend.key.size = unit(0.5, 'cm'))
taxa_plot

#################################################ALPHA-DIVERSITY###############################################
#PS for analysis
ps_analysis = ps_16S_post #Change to ps_16S_final ps_16S_pre ps_16S_post ps_MAG ps_MAG_pre ps_MAG_post

#ps_analysis = rarefy_even_depth(ps_analysis, sample.size = min(sample_sums(ps_analysis)),
                  #rngseed = TRUE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

alpha <- microbiome::alpha(ps_analysis, index = c(
  #"chao1",
  #"evenness_simpson",
  "diversity_shannon"))
alpha$ID <- sample_data(ps_analysis)$Original
alpha$Treatment <- sample_data(ps_analysis)$Treatment
alpha$Day <- sample_data(ps_analysis)$Day
alpha$Sow <- sample_data(ps_analysis)$Sow
alpha$Parity <- sample_data(ps_analysis)$Parity
alpha$Sex <- sample_data(ps_analysis)$Sex
alpha$Pen_suckling <- sample_data(ps_analysis)$Pen_suckling
alpha$Pen_weaning <- sample_data(ps_analysis)$Pen_weaning
alpha$BW_birth <- sample_data(ps_analysis)$BW_d1
alpha$BW_weaning <- sample_data(ps_analysis)$BW_d27
alpha$BW_end <- sample_data(ps_analysis)$BW_41
alpha
#write.csv(alpha, "alpha_MAG_pre.csv")

alpha_long <- melt(alpha,id.vars = c("ID","Treatment","Day","Sow","Parity","Sex","Pen_suckling","Pen_weaning","BW_birth","BW_weaning","BW_end"),
                   variable.name = "Parameter", value.name = "Index")

alpha_stat <- alpha_long %>%
  group_by(Parameter) %>%
  wilcox_test(Index~Treatment) %>%
  adjust_pvalue(method = "BH")%>%
  add_significance(p.col = "p") %>%
  #filter(p.adj<=0.05) %>%
  mutate(y.position=c(5))
alpha_stat

alpha_plot <- ggplot(alpha_long, aes(x=Treatment, y=Index)) +
  stat_boxplot(geom ='errorbar', linetype=1, width=0.5) + 
  geom_boxplot(outlier.shape = NA, aes(fill=Treatment)) +
  #geom_point(position = position_jitter(w=0.1, h=0),size=0.5) + #Add color if needed aes(color=Data)
  labs(x="", y="Shannon index") +
  #facet_wrap(~ Parameter, scales = "free_y") +
  stat_pvalue_manual(alpha_stat,label = "p", tip.length = 0.005, size = 6, y.position = 5)+ 
  scale_fill_manual(values = cols) +
  guides(fill = "none") +
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(family = 'Arial', size = 10, color = "black"),
        strip.text = element_text(family = "Arial", size = 10, color = "black")) +
  guides(fill = "none")
alpha_plot

#################################################BETA-DIVERSITY###############################################
#PS for analysis
ps_analysis = ps_16S_pre #Change to ps_16S_final ps_16S_pre ps_16S_post ps_MAG ps_MAG_pre ps_MAG_post

bray <- phyloseq::distance(ps_analysis, method = "bray")
samples_df <- data.frame(sample_data(ps_analysis))
otu_df <- as.data.frame(t(otu_table(ps_analysis)))
bray <- as.matrix(bray)
bray <- bray[rownames(samples_df),rownames(samples_df)]
pcoa = cmdscale(bray, k=3, eig=T)
pcoa_points = as.data.frame(pcoa$points)
colnames(pcoa_points) = c("x", "y", "z")
pcoa_points$ID <- samples_df$Original
pcoa_points$Treatment <- samples_df$Treatment
pcoa_points$Day <- samples_df$Day
pcoa_points$Sow <- samples_df$Sow
pcoa_points$Parity <- samples_df$Parity
pcoa_points$Sex <- samples_df$Sex
pcoa_points$Pen_suckling <- samples_df$Pen_suckling
pcoa_points$Pen_weaning <- samples_df$Pen_weaning
pcoa_points$BW_birth <- samples_df$BW_d1
pcoa_points$BW_weaning <- samples_df$BW_d27
pcoa_points$BW_end <- samples_df$BW_41
eig = pcoa$eig

#CHECK BETA-DISPERSION
betadisper_bray = betadisper(d = phyloseq::distance(ps_analysis, method = "bray"), group = sample_data(ps_analysis)$Treatment, type="centroid")
anova(betadisper_bray)
betadisper_test = TukeyHSD(betadisper_bray)
#betadisper_text = paste("Beta-dispersion, p =",format(betadisper_test$group[4], digits = 3))

#PAIRWISE ADONIS
pairwise_adonis = pairwiseAdonis::pairwise.adonis(otu_df, factors=samples_df$Treatment, sim.method="bray", p.adjust.m = "BH", perm = 999)
pairwise_adonis_text = paste("PERMANOVA,",pairwise_adonis$pairs,", p-adj =", format(pairwise_adonis$p.adjusted, digits = 3))

points_plot = ggplot(pcoa_points, aes(x=x, y=z, color=Treatment)) +
  geom_point(size=2) + #Add aes(shape=Data) if needed
  stat_ellipse(geom = "polygon",linetype = 2, aes(fill=Treatment), alpha=0.05) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=3), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=3), "%)", sep=""))+
  #facet_wrap(~ Day, scales = "free_y") +
  #ylim(-0.25,0.25) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(family = 'Arial', size = 10, color = "black"),
        strip.text = element_text(family = "Arial", size = 10, color = "black"),
        legend.text = element_text(family = "Arial", size = 10, color = "black"),
        legend.title = element_text(family = "Arial", size = 10, color = "black")) +
  guides(fill = "none")
points_plot

####################################################################LEFSE###################################################
#ps for analysis
ps_lefse = ps_16S_pre #Change to ps_16S_final ps_16S_pre ps_16S_post ps_MAG ps_MAG_pre ps_MAG_post

#ps_lefse = filter_taxa(ps_lefse, genefilter::filterfun(genefilter::kOverA(10, 0.001)), TRUE)

#LEFSE from Microbiomemarker
mm_lefse <- run_lefse(
  ps_lefse,
  #norm = "CSS",
  wilcoxon_cutoff = 1,
  group = "Treatment",
  kw_cutoff = 1,
  multigrp_strat = TRUE,
  lda_cutoff = 0,
  taxa_rank = "none"
)
mm_lefse

#Extract markers
lefse_marker_table <- data.frame(mm_lefse@marker_table)

lefse_marker_table_signif <- subset(lefse_marker_table, ef_lda >= 3 & padj <=0.05)
lefse_marker_table_signif <- subset(lefse_marker_table_signif, enrich_group == "FFT")

#Subset ps with LEFSE markers
#ps_lefse <- transform_sample_counts(ps_lefse, function(x) x/sum(x)*100) #x/sum(x)*100 or log(1 + x)
ps_lefse_marker <- prune_taxa(rownames(otu_table(ps_lefse)) %in% lefse_marker_table_signif$feature,ps_lefse)

#Melt ps and format taxa names
df_lefse <- psmelt(ps_lefse_marker) %>%
  mutate(Full_taxa = case_when(
    !is.na(Species) ~ paste(Genus, Species),
    !is.na(Genus) ~ paste(Genus, "spp."),
    TRUE ~ paste(Family, "unclassified")),
    Full_taxa_marker = paste(OTU, "-", Full_taxa)) %>%
  select(OTU, Sample, Abundance, Original, Day, Treatment, Weaning, Phylum, Class, Order, Family, Genus, Species, Full_taxa, Full_taxa_marker)

#Df with LEFSE results
lefse_res <- data.frame(OTU = lefse_marker_table_signif$feature,
                        Enriched = lefse_marker_table_signif$enrich_group,
                        LDA = lefse_marker_table_signif$ef_lda, 
                        p.adj = lefse_marker_table_signif$padj)

#Add LEFSE results as text
lefse_res$LDA_text <- paste("LDA score = ", round(lefse_res$LDA,2))
lefse_res$p.adj_text <- paste("p.adj = ", round(lefse_res$p.adj,4))

#Add the mean RA of each marker in the two groups
rownames(lefse_res) <- lefse_res$OTU
marker_mean <- df_lefse %>%
  group_by(OTU,Treatment) %>%
  summarise(Mean_Abundance = mean(Abundance, na.rm = TRUE)) %>%
  pivot_wider(names_from = Treatment, values_from = Mean_Abundance)
lefse_res <- lefse_res %>%  left_join(marker_mean, by = "OTU")

#Add taxonomical classification
marker_names <- df_lefse %>%
  select(OTU, Phylum, Class, Order, Family, Genus, Species, Full_taxa, Full_taxa_marker) %>%
  distinct()
lefse_res <- lefse_res %>%  left_join(marker_names, by = "OTU")

lefse_res
write.csv(lefse_res, file="lefse_results_16S_pre.csv")