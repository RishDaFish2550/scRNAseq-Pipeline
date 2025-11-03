# basic header for jupyter notebook 

import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = 'arial'
#formbio colors
cols = ["#4766F9","#FFB114","#00AD1D","#F76808","#9B39F9","#2EB7E5","#066B49","#E93D82","#EB9C00","#9D2DB4","#2222C3","#B24700","#0AC286","#B21A57","#FFE500","#573300","#007299"]

# each sample has a metrics folder 
# each metrics folder has a file for each sensitivity, barcode_stats.csv, matrix_stats.csv

## path2pipseeker results dir
#path2results = "/data/250211_pipseq_organoids/250211_pipseq_results"
## list of sample names 
#samplenames = ['DSO_D0','DSO_D3','DSO_D6','DSO_D9','HSO_D0','HSO_D3','HSO_D6','HSO_D9']

samplenames = []
with open("filelist_path_samplenicknames.txt",'r') as filelist:
    for f in filelist:

        path2results, samplename, samplenickname = f.strip().split()
        
        p = F"{path2results}/metrics"
        samplenames.append(p)
        plt.figure(figsize=(12,4))
        for n in [1,2,3,4,5]:
            plt.subplot(2,5,n)
            kneeplot = F"{path2results}/cell_calling/sensitivity_{n}/barcode_rank.png"
            plt.title(F"{samplenickname}\nsensitivity_{n}")
            plt.imshow(plt.imread(kneeplot))
            plt.xticks(ticks=[])
            plt.yticks(ticks=[])
            plt.subplot(2,5,n+5)
            umap = F"{path2results}/clustering/sensitivity_{n}/umap/graph.png"
            plt.imshow(plt.imread(umap))
            plt.xticks(ticks=[])
            plt.yticks(ticks=[])
            plt.tight_layout()
        plt.savefig(F"FIGURE_pipseeker_{samplenickname}.png",dpi=450)
        plt.close()

df = pd.read_csv("TABLE_METRICS_ALL_SAMPLES.csv")
celltype = [n.split('_')[0] for n in df['sample_nickname']]
df['celltype'] = celltype

df_filt = df[df['sensitivity']==5]
percent_mapped = round(100*df_filt['num_mapped_reads']/df_filt['num_input_reads'],1)
plt.figure(figsize=(12,4))
x = np.arange(len(df_filt))
plt.bar(x,df_filt['num_input_reads'],color = cols[0])
plt.bar(x,df_filt['num_mapped_reads'],color = cols[2])
plt.xticks(x, df_filt['sample_nickname'],rotation=90)

for i in x:
    
    if percent_mapped.iloc[i]>100:
        plt.bar(i,df_filt['num_mapped_reads'].iloc[i],color=cols[2])
        plt.bar(i,df_filt['num_input_reads'].iloc[i],color=cols[0])
        
    else: 
        plt.bar(i,df_filt['num_input_reads'].iloc[i],color = cols[0])
        plt.bar(i,df_filt['num_mapped_reads'].iloc[i],color = cols[2])

    plt.text(i,df_filt['num_mapped_reads'].iloc[i],percent_mapped.iloc[i],rotation=90,color='black')
plt.tight_layout()
plt.savefig(F"FIGURE_STATS.png",dpi=450)

plt.figure(figsize=(12,6))
x = np.arange(len(df))
x_labels = []
for ix in x:
    sensitivity = df.iloc[ix]['sensitivity']
    sample_nickname = df.iloc[ix]['sample_nickname']
    
    plt.bar(ix, df.iloc[ix]['num_cell_barcodes'],color = cols[sensitivity])
    x_labels.append(F"{sample_nickname} (s{sensitivity})")

xgroup = np.arange(0,len(df),5)
plt.xticks(xgroup,[df['sample_nickname'].iloc[i] for i in xgroup],rotation=90)
plt.ylabel('num_cell_barcodes')
plt.tight_layout()
plt.savefig("FIGURE_STATS_numberCellBarcodes.png",dpi=400)

plt.figure()
plt.title("legend")
for i in [1,2,3,4,5]:
    plt.plot(i,i,'o',markersize=22,color=cols[i])
    plt.text(i-1,i,F"sensitivity_{i}")
plt.savefig("FIGURE_STATS_numberCellBarcodes_LEGEND.png",dpi=400)

plt.figure(figsize=(16,6))
x = np.arange(len(df))
x_labels = []
for ix in x:
    sensitivity = df.iloc[ix]['sensitivity']
    sample_nickname = df.iloc[ix]['sample_nickname']
    
    plt.bar(ix, df.iloc[ix]['mean_reads_per_cell'],color = cols[sensitivity])
    x_labels.append(F"{sample_nickname} (s{sensitivity})")

xgroup = np.arange(0,len(df),5)
plt.xticks(xgroup,[df['sample_nickname'].iloc[i] for i in xgroup],rotation=90)
plt.ylabel('mean_reads_per_cell')
plt.tight_layout()
plt.savefig("FIGURE_STATS_meanReadsPerCell.png",dpi=400)

# plot number of cells at each sensitivity
plt.figure(figsize=(12,28))

for s in [1,2,3,4,5]:
    
    df_s = df[df['sensitivity']==s]
    
    plt.subplot(5,1,s)
    plt.title(F"sensitivity_{s}")
    x = np.arange(len(df_s))
    plt.bar(x, df_s['reads_in_cells_pct'],color='navy',edgecolor='w')
    plt.xticks(x, df_s['sample_nickname'],rotation=90)
    plt.ylabel('reads_in_cells_pct')
    plt.ylim(0,100)
    plt.tight_layout()
    
plt.savefig("FIGURE_STATS_percentReadsCells_perSensitivity.png",dpi=400)
