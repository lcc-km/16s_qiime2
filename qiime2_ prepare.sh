
######### 准备工作 ############## 
## 导入数据
########    https://docs.qiime2.org/2023.5/data-resources/
#全长分类器构建
# qiime feature-classifier fit-classifier-naive-bayes \
#   --i-reference-reads  silva-138-99-seqs.qza \
#   --i-reference-taxonomy silva-138-99-tax.qza \
#   --o-classifier silva-138-nr99-full-classifier.qza

# ########### 16S（V3-V4）341F~805R
# qiime feature-classifier extract-reads \
#   --i-sequences silva-138-99-seqs.qza \
#   --p-f-primer CCTACGGGNGGCWGCAG \
#   --p-r-primer GACTACHVGGGTATCTAATCC \
#   --p-n-jobs 12 \
#   --o-reads silva-138-ssu-nr99-16s_V3_V4_341F_R805.qza
# #qiime feature-table tabulate-seqs --i-data silva-138-ssu-nr99-16s_V3_V4_341F_R805.qza --o-visualization silva-138-ssu-nr99-16s_V3_V4_341F_R805.qzv

# nohup qiime feature-classifier fit-classifier-naive-bayes \
#   --i-reference-reads silva-138-ssu-nr99-16s_V3_V4_341F_R805.qza \
#   --i-reference-taxonomy silva-138-99-tax.qza \
#   --o-classifier silva-138-ssu-nr99-16s_V3_V4_341F_R805.qza-classifier.qza &

###   添加质控 #######
# qiime rescript cull-seqs \
#     --i-sequences silva-138-99-seqs.qza \
#     --o-clean-sequences silva-138-99-seqs_cleaned.qza
#长度过滤
# qiime rescript filter-seqs-length-by-taxon \
#     --i-sequences silva-138-99-seqs_cleaned.qza \
#     --i-taxonomy silva-138-99-tax.qza \
#     --p-labels Archaea Bacteria Eukaryota \
#     --p-min-lens 900 1200 1400 \
#     --o-filtered-seqs silva-138-99-seqs_cleaned-filt.qza \
#     --o-discarded-seqs silva-138-99-seqs_cleaned-discard.qza
#重复序列合并
# qiime rescript dereplicate \
#     --i-sequences silva-138-99-seqs_cleaned-filt.qza  \
#     --i-taxa silva-138-99-tax.qza \
#     --p-mode 'uniq' \
#     --p-rank-handles 'disable' \
#     --o-dereplicated-sequences silva-138-99-seqs_cleaned-filt-derep-uniq.qza \
#     --o-dereplicated-taxa silva-138-99-tax-derep-uniq.qza
#全长分类器构建
# nohup qiime feature-classifier fit-classifier-naive-bayes \
#   --i-reference-reads  silva-138-99-seqs_cleaned-filt-derep-uniq.qza \
#   --i-reference-taxonomy silva-138-99-tax-derep-uniq.qza \
#   --o-classifier silva-138-99-tax-clean-uniq-classifier.qza &
##特异引物分类器构建
#### 16S（V3-V4）338F~806R 
#截取序列
# qiime feature-classifier extract-reads \
#     --i-sequences silva-138-99-seqs_cleaned-filt-derep-uniq.qza \
#    --p-f-primer CCTACGGGNGGCWGCAG \
#    --p-r-primer GACTACHVGGGTATCTAATCC \
#     --p-n-jobs 24 \
#     --o-reads silva-138-99-seqs_cleaned-filt-derep-uniq-seqs-341F_805R.qza
#合并重复
# qiime rescript dereplicate \
#     --i-sequences silva-138-99-seqs_cleaned-filt-derep-uniq-seqs-341F_805R.qza \
#     --i-taxa silva-138-99-tax-derep-uniq.qza \
#     --p-rank-handles 'disable' \
#     --p-mode 'uniq' \
#     --o-dereplicated-sequences silva-138-99-seqs_cleaned-filt-derep-uniq-seqs-341F_805R-uniq.qza \
#     --o-dereplicated-taxa  silva-138-99-tax-derep-uniq-341F_805R-derep-uniq.qza
#构建分类器
# nohup qiime feature-classifier fit-classifier-naive-bayes \
#     --i-reference-reads silva-138-99-seqs_cleaned-filt-derep-uniq-seqs-341F_805R-uniq.qza \
#     --i-reference-taxonomy silva-138-99-tax-derep-uniq-341F_805R-derep-uniq.qza \
#     --o-classifier silva-138-99-cleaned-341F_805R_V3-V4-classifier.qza &

############  ITS ############## 
# qiime tools import \
#   --type 'FeatureData[Sequence]' \
#   --input-path sh_refs_qiime_ver9_99_s_all_29.11.2022.fasta \
#   --output-path 99_UNITE.qza
#   # 导入物种分类信息，6s
# qiime tools import \
#   --type 'FeatureData[Taxonomy]' \
#   --input-format HeaderlessTSVTaxonomyFormat \
#   --input-path sh_taxonomy_qiime_ver9_99_s_all_29.11.2022.txt \
#   --output-path all_99_UNITE_taxonomy.qza

# ### 引物提取参考序列的扩增区段 Extract reference reads
# qiime feature-classifier extract-reads \
#   --i-sequences 99_UNITE.qza \
#   --p-f-primer GTCCCTGCCCTTTGTACACA \
#   --p-r-primer TTTCGCTGCGTTCTTCATCG \
#   --p-n-jobs 64 \
#   --o-reads all_99_UNITE_IVS1_30F_217R.qza

# qiime feature-table tabulate-seqs --i-data all_99_UNITE_IVS1_30F_217R.qza --o-visualization all_99_UNITE_IVS1.qzv

# # 基于筛选的指定区段，生成实验特异的分类器，14m
# qiime feature-classifier fit-classifier-naive-bayes \
#   --i-reference-reads all_99_UNITE_IVS1_30F_217R.qza \
#   --i-reference-taxonomy 99_UNITE_taxonomy.qza \
#   --o-classifier all_UNITE_IVS1-classifier.qza

# #### 全长 ####### ？ 
# nohup qiime feature-classifier fit-classifier-naive-bayes \
#   --i-reference-reads 99_UNITE.qza \
#   --i-reference-taxonomy 99_UNITE_taxonomy.qza \
#   --o-classifier UNITE_full-classifier.qza &

#############  Greengenes2 ############# 
# qiime tools import \
#   --type 'FeatureData[Sequence]' \
#   --input-path sh_refs_qiime_ver9_99_s_all_29.11.2022.fasta \
#   --output-path all_99_UNITE.qza
#   # 导入物种分类信息，6s
# qiime tools import \
#   --type 'FeatureData[Taxonomy]' \
#   --input-format HeaderlessTSVTaxonomyFormat \
#   --input-path 2022.10.taxonomy.id.tsv \
#   --output-path 2022.10_taxonomy.qza

# ###### 全长
# nohup qiime feature-classifier fit-classifier-naive-bayes \
#   --i-reference-reads 2022.10.seqs.fna.qza \
#   --i-reference-taxonomy 2022.10_taxonomy.qza \
#   --o-classifier Greengenes2-classifier.qza &

########### 16S（V3-V4）341F~805R
# qiime feature-classifier extract-reads \
#   --i-sequences 2022.10.seqs.fna.qza \
#   --p-f-primer CCTACGGGNGGCWGCAG \
#   --p-r-primer GACTACHVGGGTATCTAATCC \
#   --p-n-jobs 24 \
#   --o-reads Greengenes2_341F_805R_V3_V4_341F_R805.qza
#qiime feature-table tabulate-seqs --i-data Greengenes2_341F_805R_V3_V4_341F_R805.qza --o-visualization Greengenes2_341F_805R_V3_V4_341F_R805.qzv

# nohup qiime feature-classifier fit-classifier-naive-bayes \
#   --i-reference-reads Greengenes2_341F_805R_V3_V4_341F_R805.qza \
#   --i-reference-taxonomy 2022.10_taxonomy.qza \
#   --o-classifier Greengenes2_341F_805R_V3_V4_341F_R805-classifier.qza &