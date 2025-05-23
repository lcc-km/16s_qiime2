# 按barcode拆分样品 Demultiplexing sequences
#qiime demux emp-single \
#  --i-seqs emp-single-end-sequences.qza \
#  --m-barcodes-file sample-metadata.tsv \
#  --m-barcodes-column barcode-sequence \
#  --o-per-sample-sequences demux.qza \
#  --o-error-correction-details demux-details.qza
# 获取拆分结果

### sample_metadata.tsv > group.txt for R plot
echo -e 'sample\tgroup' > group_h
cut -f 1,3 sample-metadata.tsv | grep -v "#" | sed '1d' > group_l
cat group_h group_l > group.txt

#########  manifest_file  ############
#sample-id,absolute-filepath,direction
#SRR7635469,/mnt/vol2/lucc/data_test/05.mNGS/02.16s/PRJNA481576/fa/SRR7635469_1.fastq,forward
#SRR7635469,/mnt/vol2/lucc/data_test/05.mNGS/02.16s/PRJNA481576/fa/SRR7635469_2.fastq,reverse
#SRR7635475,/mnt/vol2/lucc/data_test/05.mNGS/02.16s/PRJNA481576/fa/SRR7635475_1.fastq,forward
#SRR7635475,/mnt/vol2/lucc/data_test/05.mNGS/02.16s/PRJNA481576/fa/SRR7635475_2.fastq,reverse

### manifest_file
ls *gz | cut -d . -f 1 | sort -u > coding
echo "sample-id,absolute-filepath,direction" > manifest.txt
cat coding | while read id
do
arr=($id)
sampleID=${arr[0]}
F=`realpath ${sampleID}*R1.fq.gz`
R=`realpath ${sampleID}*R2.fq.gz`
echo "${sampleID},${F},forward">> manifest.txt; 
echo "${sampleID},${R},reverse">> manifest.txt
done




## 导入数据
#docker run -tiv /mnt/:/mnt/ quay.io/qiime2/core:2023.2 qiime tools import \
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest.txt \
  --input-format PairedEndFastqManifestPhred33 \
  --output-path demux-paired-end.qza 

#可视化文件demux-paired-end.qza
qiime demux summarize \
  --i-data demux-paired-end.qza \
  --o-visualization demux-paired-end.qza.qzv

##. 生成特征表和代表序列
######## 方法1  data2 降噪
#paired
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux-paired-end.qza \
  --p-n-threads 0 \
  --p-trim-left-f 15 --p-trim-left-r 15 \
  --p-trunc-len-f 0 --p-trunc-len-r 0 \
  --o-table dada2-table.qza \
  --o-representative-sequences dada2-rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza
# 确定使用dada2结果并导入主流程
cp dada2-table.qza table.qza
cp dada2-rep-seqs.qza rep-seqs.qza


### 物种注释 
qiime feature-classifier classify-sklearn \
  --i-classifier /mnt/vol1/database/qiime2/Silva/done/silva-138-ssu-nr99-16s_V3_V4_341F_R805.qza-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza

# 导出物种注释
mkdir -p taxonomy
qiime tools export \
  --input-path taxonomy.qza \
  --output-path taxonomy

# 按种水平合并，并统计
## 按种水平合并，6s    
qiime taxa collapse \
--i-table table.qza \
--i-taxonomy taxonomy.qza \
--p-level 6 \
--o-collapsed-table table-6.qza
# 导出
qiime tools export \
  --input-path table-6.qza \
  --output-path taxonomy

biom convert -i taxonomy/feature-table.biom \
  -o taxonomy/table-6.txt \
  --to-tsv  

# 格式化特征表，添加伪计数，6s
qiime composition add-pseudocount \
  --i-table table-6.qza \
  --o-composition-table comp-table-6.qza
biom convert -i comp-table-6.qza -o comp-table-6.txt --to-tsv 

qiime feature-table summarize \
--i-table table-6.qza\
--o-visualization table-6.qzv

############# 
qiime taxa barplot --i-table table.qza --i-taxonomy taxonomy.qza --m-metadata-file sample-metadata.tsv --o-visualization taxa-bar-plots.qzv 


#########################################################################
# 特征表统计 Feature表统计
qiime feature-table summarize \
--i-table table.qza \
--o-visualization table.qzv
# 代表序列统计
qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

qiime metadata tabulate \
  --m-input-file denoising-stats.qza \
  --o-visualization denoising-stats.qzv

# 过滤特征表：过滤样本
##qiime feature-table filter-samples \
#--i-table table.qza \
#--o-filtered-table table_filtered.qza \
#--m-metadata-file sample-to-keep.tsv 

## 3. Alpha和beta多样性分析
### 构建进化树用于多样性分析 
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

### 计算核心多样性 
  # 采样深度通常选择最小值，来自table.qzv
  mkdir -p core-metrics-results
  qiime diversity core-metrics-phylogenetic \
    --i-phylogeny rooted-tree.qza \
    --i-table table.qza \
    --p-sampling-depth 92497 \
    --m-metadata-file metadata.txt \
    --output-dir core-metrics-results

# 输出结果包括多种多样性结果，文件列表和解释如下：
# beta多样性bray_curtis距离矩阵 bray_curtis_distance_matrix.qza 
# alpha多样性evenness(均匀度，考虑物种和丰度)指数 evenness_vector.qza
# alpha多样性faith_pd(考虑物种间进化关系)指数 faith_pd_vector.qza
# beta多样性jaccard距离矩阵 jaccard_distance_matrix.qza
# alpha多样性observed_otus(OTU数量)指数 observed_otus_vector.qza
# alpha多样性香农熵(考虑物种和丰度)指数 shannon_vector.qza
# beta多样性unweighted_unifrac距离矩阵，不考虑丰度 unweighted_unifrac_distance_matrix.qza
# beta多样性unweighted_unifrac距离矩阵，考虑丰度 weighted_unifrac_distance_matrix.qza


## alpha多样性可视化及组间显著性分析
qiime diversity alpha-group-significance \
--i-alpha-diversity shannon_vector.qza \
--m-metadata-file sample_metadata.tsv \
--o-visualization shannon_vector.qzv


### 物种注释 
qiime feature-classifier classify-sklearn \
  --i-classifier /mnt/vol1/database/qiime2/classifier_gg_13_8_99_V3-V4.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza