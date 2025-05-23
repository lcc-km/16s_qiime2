#BiocManager::install("rcompanion")
library(ggplot2)
library(reshape2)
library(microeco)
library(file2meco)
library(agricolae)
library(dplyr)
library(ggdendro)
library(Maaslin2)
library(WGCNA)
library(rgexf)
library(randomForest)
library(tidytree)
library(ggtree)
library(DESeq2)


pwd <- "/mnt/vol2/lucc/data_test/05.mNGS/02.16s/PRJNA737206"
setwd(pwd)

abund_file_path <- paste(pwd,"dada2-table.qza",sep = "/")
sample_file_path <- paste(pwd,"SraRunTable.txt",sep = "/")
taxonomy_file_path <-paste(pwd,"taxonomy.qza",sep = "/")
tree_data <- paste(pwd,'rooted-tree.qza',sep = "/")
rep_data <- paste(pwd,'dada2-rep-seqs.qza',sep = "/")
dataset <- qiime2meco(abund_file_path, sample_table = sample_file_path,
                      taxonomy_table = taxonomy_file_path, 
                      phylo_tree = tree_data, rep_fasta = rep_data, auto_tidy = TRUE)
dataset
d_group <- read.table(sample_file_path,sep = "\t",header = T)
table(d_group$group_2)
table(d_group$group_1)

list_group = unique(d_group$group_2)
group_name <- "group_2"

#####   bar plot example
dataset$cal_abund()
t1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 10)
t1$plot_bar(others_color = "grey70", facet = group_name, xtext_keep = FALSE, legend_text_italic = FALSE)
### heatmap
t1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 40)
t1$plot_heatmap(facet = group_name, xtext_keep = FALSE, withmargin = FALSE)

##  pie chart
t1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 6, groupmean = group_name)
# all pie chart in one row
t1$plot_pie(facet_nrow = 1, add_label = TRUE)

## venn plot
dataset1 <- dataset$merge_samples(use_group = group_name)
# dataset1 is a new microtable object
# create trans_venn object
t1 <- trans_venn$new(dataset1, ratio = 'seqratio')
t1$plot_venn()

#####  alpha class
t1 <- trans_alpha$new(dataset = dataset, group = group_name)
t1$cal_diff(method = "anova")
t1$plot_alpha(measure = "Chao1", add_sig_text_size = 6, boxplot_add = "jitter", order_x_mean = TRUE,
              alpha = 0.05)

### beta class
dataset$cal_betadiv()
t1 <- trans_beta$new(dataset = dataset, group = group_name, measure = "bray")
t1$cal_ordination(ordination = "NMDS")
t1$plot_ordination(plot_color = group_name, plot_shape = group_name, 
                   plot_type = c("point", "ellipse"))+theme_classic()  +ylim(-0.25,0.25)

# calculate and plot sample distances within groups
t1$cal_group_distance(within_group = TRUE)
# return t1$res_group_distance
# perform Wilcoxon Rank Sum and Signed Rank Tests
t1$cal_group_distance_diff(method = "wilcox")
# plot_group_order parameter can be used to adjust orders in x axis
t1$plot_group_distance(boxplot_add = "mean")


### Clustering plot
d1 <- clone(dataset)
d1$sample_table %>% subset(list_group %in% list_group)
d1$tidy_dataset()
t1 <- trans_beta$new(dataset = d1, group = group_name)
# use replace_name to set the label name, group parameter used to set the color
t1$plot_clustering(group = group_name, replace_name = group_name)



##########  LEfSe ###### 
# 执行lefse分析

t1 <- trans_diff$new(dataset = dataset, method = "lefse", 
                     group = group_name, alpha = 0.05,
                     lefse_subgroup = NULL,
                     taxa_level='Genus',
                     p_adjust_method = "none")
head(t1$res_diff)
d_r <- t1$res_diff#[1:28, c(1, 3, 4, 6)]
t1$plot_diff_bar(use_number = 1:30, width = 0.8, group_order = list_group) #



t1$plot_diff_abund(use_number = 1:30, group_order = list_group)

t1 <- trans_diff$new(dataset = dataset, method = "rf", 
                     group = list_group , taxa_level = "Trait",
                     p_adjust_method = "none")
t1$plot_diff_bar(use_number = 1:30)


t1$plot_diff_cladogram(use_taxa_num = 200, 
                          use_feature_num = 50, 
                          clade_label_level = 5, 
                          group_order = list_group)


LEfSe_result <- d_r
write.table(LEfSe_result,"/mnt/vol2/lucc/data_test/05.mNGS/02.16s/LEfSe_result.txt",row.names = F,col.names = T,quote = F,sep = "\t")
count_result_Genus <- dataset$taxa_abund$Genus
write.table(count_result_Genus,"/mnt/vol2/lucc/data_test/05.mNGS/02.16s/PRJNA737206_count_result_Genus.txt",row.names = T,col.names = T,quote = F,sep = "\t")


########### network
t1 <- trans_network$new(dataset = dataset, cor_method = "spearman", filter_thres = 0.001)
# construct network; require igraph package
t1$cal_network(COR_p_thres = 0.01, COR_optimization = TRUE)
# use arbitrary coefficient threshold to contruct network
t1$cal_network(COR_p_thres = 0.01, COR_cut = 0.7)
# invoke igraph cluster_fast_greedy function for this undirected network 
t1$cal_module(method = "cluster_fast_greedy")
# require rgexf package to be installed
t1$save_network(filepath = "network.gexf")

# ###### Explainable
# t1 <- trans_env$new(dataset = dataset)


####  Maaslin 
fit_data <- Maaslin2(
  d3, group, 'demo_output',
  fixed_effects = c('id', 'Group'),
  standardize = FALSE)

#########   单独的 a 多样性计算
library(vegan)

d <- read.table('/mnt/vol2/lucc/data_test/05.mNGS/02.16s/PRJNA737206/taxonomy/table-6.txt',sep = "\t",header = F)
d_coding <- read.table('/mnt/vol2/lucc/data_test/05.mNGS/02.16s/PRJNA737206/coding',sep = "\t",header = F)
sampleID <- d_coding$V1
sampleID <- append(sampleID, 'taxo', 0)
colnames(d) <- sampleID
otu_list <- d$taxo
d <- d[,c(2:length(d))]
d2 <- as.data.frame(lapply(d,as.numeric))
colSums(d2)
#  抽平分析
d_Flattening <- as.data.frame(t(rrarefy(t(d2),min(colSums(d2)))))
## Chao1 ACE
observed_species <- data.frame(t(estimateR(t(d_Flattening))[c(1,2,4),]))
## shanon
observed_species$Shannon <- diversity(t(d_Flattening),index = 'shannon',base = exp(1))
## gini-simpson
observed_species$Gini_simpson <- diversity(t(d_Flattening),index = 'simpson')
######
observed_species$goods_coverage <- 1 - (rowSums(t(d_Flattening)==1) / rowSums(t(d_Flattening)))
## group 
pwd <- "/mnt/vol2/lucc/data_test/05.mNGS/02.16s/PRJNA737206"
sample_file_path <- paste(pwd,"SraRunTable.txt",sep = "/")
d_group <- read.table(sample_file_path,sep = "\t",header = T)

observed_species$id <- rownames(observed_species)

observed_species_plot <- merge(observed_species,d_group,by='id')
observed_species_plot2 <- melt(observed_species_plot)

library(ggplot2)
library(reshape2)

colnames(observed_species_plot2)

colnames(d_p)
ggplot(observed_species_plot2, aes(x=group_1, y=value,fill=group_1)) + 
  geom_boxplot(notch=FALSE,outlier.shape = NA)+
  facet_wrap(~variable,scales = "free_y")


S_obs <- aggregate(observed_species_plot$S.obs, list(observed_species_plot$group_1), summary) 

aggregate(observed_species_plot$S.chao1, list(observed_species_plot$group_1), summary)

aggregate(observed_species_plot$S.ACE, list(observed_species_plot$group_1), summary)




