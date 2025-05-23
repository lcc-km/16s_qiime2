#!/usr/bin/python3
import pandas as pd
import numpy as np
import argparse
import matplotlib as mpl
import matplotlib.pyplot as plt
import palettable
from palettable.cartocolors.sequential import DarkMint_4
#from itertools import repeat

parser = argparse.ArgumentParser(description='output kraken report list')
parser.add_argument('-i', help='the sample ID')
parser.add_argument('-p', help='path of the sample to be analyzed')
parser.add_argument('-f', help='kraken2 result file')
parser.add_argument('-c', help='number of cutoff')
args = parser.parse_args()

sampleID = args.i
pwd = args.p
file_in =  args.f
cut_off = args.c
cut_off =float(cut_off)

# sampleID = '0807T1-16S'
# pwd = '/Volumes/mac_ExFAT/微创文件/3.mNGS/seq_data/16s/result_16s'
# ref_core_table_in = "/Volumes/mac_ExFAT/微创文件/3.mNGS/注释/core_table.txt"
# file_in = "/Volumes/mac_ExFAT/微创文件/3.mNGS/seq_data/16s/result_16s/0807T1-16S.table.6.txt"
# cut_off = 0.01

## reas data
d = pd.read_table(file_in,sep="\t",header=None,comment='#')
#d_ref1 = pd.read_table(ref_core_table_in,sep="\t")
# Acquired total reads
d.columns =['taxonomy','count']
all_count = d['count'].sum()
## get g__
d_g = d[d.taxonomy.str.contains('g__')]
## Count percent
d_g.loc[:,'percent'] = d_g.loc[:,'count'] / all_count
# other
d_g_bed = d_g[d_g.loc[:,'percent'] < cut_off]
### exclude ref
d_g_bed_true_NF = d_g_bed
d_g_bed_true_pathogenic = d_g.iloc[3:4]

##
d_g_bed_tmp = {
    'taxonomy' : ["non-significant"],
    'count' : [d_g_bed_true_NF['count'].sum()],
    'percent' : [d_g_bed_true_NF['percent'].sum()]
}
d_g_bed_tmp = pd.DataFrame(d_g_bed_tmp)
# cut_off
d_g = d_g[d_g.loc[:,'percent'] >= cut_off]
# merge
d_g2 = pd.concat([d_g,d_g_bed_true_pathogenic])
# sort
d_g2 = d_g2.sort_values(by=['percent'],ascending=[False])
d_g2 = d_g2.reset_index(drop=True)
# merge to result
d_g2_r = pd.concat([d_g2,d_g_bed_tmp])
d_g2_r = d_g2_r.reset_index(drop=True)
d_g2_r.loc[:,'percent'] = round(d_g2_r.loc[:,'percent'] *100,2)
###  暂时去重
d_g2_r = d_g2_r.drop_duplicates(keep='last')
d_g2_r = d_g2_r.reset_index(drop=True)
##
f_out = pwd + "/" + sampleID + ".qiime2_result.txt"
d_g2_r.to_csv(f_out, sep='\t', index=False)

###
d_g2_rs = d_g2_r['taxonomy'].str.split(";",expand=True)
d_g2_rs.columns = ['k__','p__','c__','o__','f__','g__']
index_NF =  d_g2_rs[d_g2_rs['k__'] == 'non-significant'].index.values
d_g2_rs.iloc[index_NF,5] = 'non-significant'
d_g2_r['g__'] = d_g2_rs['g__']
d_g2_r.loc[:,'lable'] = d_g2_r.loc[:,'g__'].astype(str) +" (" +d_g2_r.loc[:,'percent'].astype(str) +"%)"

#### plot
L = len(d_g2_r.taxonomy)  ## 设置使用颜色的数量
figure_size = (10, 10)  ##  设置图片大小
plt.figure(figsize=figure_size, dpi=100)
patches, texts = plt.pie(d_g2_r['percent'],
                         colors=plt.get_cmap(DarkMint_4.mpl_colormap)(np.linspace(0, 1, L)),  ###  等份取色号
                         )
plt.legend(patches, list(d_g2_r['lable']),  # 添加图例
           title=('composition ratio'),
           loc='upper center', bbox_to_anchor=(0.5, 0.1),
           ncol=2,  # 控制图例中按照两列显示，默认为一列显示，
           )
# plt.show()
plot_output = pwd + "/plot/" + sampleID + '.qiime2_16s.jpg'
print(plot_output)
plt.savefig(plot_output, bbox_inches='tight', pad_inches=0.0, dpi=300)