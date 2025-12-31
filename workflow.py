

#Adapter trimming, read mapping and read counts

# Adapter trimming using Trimmomatic
# Illumina TruSeq adapter sequences were removed
# Reads shorter than 36 bp were discarded

Trimmomatic (v0.38) with the following parameters: 
ILLUMINACLIP:2:30:10, LEADING:20, TRAILING:20, SLIDINGWINDOW:4:15, MINLEN:36, 
based on Illumina adapter sequences (AGATCGGAAGAGCACACGTCTGAACTCCAGTCA, AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT). 

HISAT2 (v2.1.0)

featureCounts (v1.6.3)


# Principal component analysis and differential expression analysis
scikit-learn (Python v3.10, v1.0.2)

# Differential expression analysis using DESeq2
# Genes with adjusted p-value < 0.05 were considered significant
DESeq2 (R v4.2.3, v1.36.0)

# False discovery rate (FDR)-adjusted p-value threshold of <0.05 after Benjamini–Hochberg correction

# Reformatting results from DESeq analysis

import pandas as pd

input_file = "path/TPM_Matrix.csv"  
output_file = "path/TPM_Matrix_re.csv"  


df = pd.read_csv(input_file)


gene_column = df['Gene']
sample = ['Co_1', 'Co_2', 'Co_3', 'HKCo_1', 'HKCo_2', 'HKCo_3', 'Mono_1', 'Mono_2', 'Mono_3']
sample_data = df[sample]

labels = sample * len(df)


class_labels = [
    'Co' if 'Co_' in col and not 'HKCo' in col else
    'HKCo' if 'HKCo_' in col else
    'Mono' for col in sample_columns
] * len(df)

reformatted_data = pd.DataFrame({
    'Gene': gene_column.repeat(len(sample_columns)).values,
    'Label': labels,
    'Class': class_labels,
    'Data': sample_data.values.flatten()
})


reformatted_data.to_csv(output_file, sep="\t", index=False)

print(f"saved {output_file}")


# Extract gene names from GTF file and merge with expression data

import pandas as pd
import re


gtf_file = "path/Sus_scrofa.Sscrofa11.1.112.gtf" 
csv_file = "path/TPM_Matrix_re.csv"  
output_csv = "path/TPM_Matrix_rere.csv"  


def parse_gtf(gtf_file):
    gene_dict = {}
    with open(gtf_file, 'r') as file:
        for line in file:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if fields[2] == "transcript":  
                attributes = fields[8]
                transcript_id_match = re.search(r'gene_id "([^"]+)"', attributes)
                gene_name_match = re.search(r'gene_name "([^"]+)"', attributes)
                if transcript_id_match and gene_name_match:
                    transcript_id = transcript_id_match.group(1)
                    gene_name = gene_name_match.group(1)
                    gene_dict[transcript_id] = gene_name
    return gene_dict


gene_dict = parse_gtf(gtf_file)
df = pd.read_csv(csv_file)
df['gene_name'] = df['Gene'].map(gene_dict)
df.to_csv(output_csv, index=False)

print(f"save: {output_csv}")


# Input file: reformatted TPM matrix
# "Data1" represents normalized expression values with +1 added to all counts
# (to avoid division by zero during fold-change calculation)

import pandas as pd

input_csv = "path/TPM_Matrix_rere.csv"  
output_csv = "path/TPM_Matrix_rere_FC.csv"  

df = pd.read_csv(input_csv)
results = []


# Calculate fold change relative to the Mono condition for each gene
for gene, group in df.groupby("Gene"):
  
  　# Calculate the mean expression value in the Mono condition
    mono_mean = group[group["Class"] == "Mono"]["Data1"].mean()
    
 　　# Calculate fold change relative to the Mono condition
    group["Fold_change"] = group["Data1"] / mono_mean  
    results.append(group)


result_df = pd.concat(results)

result_df.to_csv(output_csv, index=False)

print(f"save: {output_csv}")


# Log2 transformation
import pandas as pd
import numpy as np


input_csv = "path/TPM_Matrix_rere_FC.csv"  
output_csv = "path/TPM_Matrix_rere_FClog2.csv" 


df = pd.read_csv(input_csv)
df['log2FC'] = np.log2(df['Fold_change'])
df.to_csv(output_csv, index=False)

print(f"save: {output_csv}")


# Step 1: Collect DESeq2 result files and extract significant DEGs
# Step 2: Extract expression data for DEGs from TPM matrix

import pandas as pd
import glob


# List to store filtered DEG tables
dataframes = []
file_paths = glob.glob("path/ .csv")

if not file_paths:
    print("Error1")
else:
    for file_path in file_paths:
        print(f"Processing: {file_path}")
        
        with open(file_path, 'r') as f:
            print("file:\n")
            for _ in range(5):
                print(f.readline())
        
        try:

            df = pd.read_csv(file_path, delimiter=',')  
            print(f"Columns: {df.columns}")
            

            df['Gene'] = df['Gene'].astype(str)
            

            if 'padj' in df.columns and 'log2FoldChange' in df.columns:

                # Filter for significantly differentially expressed genes
                # padj < 0.001 and valid Ensembl gene IDs
                filtered_data = df[
                    (df['padj'] < 0.001) & df['Gene'].str.startswith('ENSSSCG', na=False)
                ]
                print(f"Filtered data shape: {filtered_data.shape}")
                
                if not filtered_data.empty:
  
                    dataframes.append(filtered_data.set_index('Gene')[['log2FoldChange']])
                else:
                    print("Error2")
            else:
                print("Error3")
        
        except Exception as e:
            print(f"Error reading {file_path}: {e}")


if dataframes:
    merged_data = pd.concat(dataframes, axis=1)
    merged_data.to_csv("path/.csv", sep='\t')
    print("Merged")
else:
    print("Error4")


import pandas as pd
df_a = pd.read_csv("path/DEGs.csv", sep="\t")
df_b = pd.read_csv("path/TPM_Matrix_rere_FClog2.csv", sep="\t")
result_df = df_b[df_b['Gene'].isin(df_a['Gene'])]
result_df.to_csv("_DEGs.csv", index=False)



# Principal component analysis
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


file_path = "path/TPM_Matrix.csv"  
df = pd.read_csv(file_path)  

df.set_index("Gene", inplace=True)


conditions = ["Co", "Mono", "HKCo"]  
columns = {cond: [col for col in df.columns if cond in col] for cond in conditions}


data = df.values.T  
samples = list(df.columns)  


scaler = StandardScaler()
scaled_data = scaler.fit_transform(data)


pca = PCA(n_components=2)  
pca_result = pca.fit_transform(scaled_data)


pca_df = pd.DataFrame(pca_result, columns=["PC1", "PC2"])
pca_df["Sample"] = samples
pca_df["Condition"] = [
    "Co" if "Co" in sample and "HKCo" not in sample else
    "Mono" if "Mono" in sample else
    "HKCo" if "HKCo" in sample else "Unknown"
    for sample in samples
]  


plt.figure(figsize=(8, 6))
colors = {"Co": "khaki", "Mono": "lightskyblue", "HKCo": "palevioletred"}
for condition, color in colors.items():
    subset = pca_df[pca_df["Condition"] == condition]
    plt.scatter(subset["PC1"], subset["PC2"], label=condition, color=color, s=100, alpha=0.7)


plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.2f}%)")
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.2f}%)")
#
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
plt.grid(False)
sns.despine()
plt.tight_layout()
plt.xlim(-150, 150)
plt.ylim(-100, 100)


output_path = "path/PCA.png"
plt.savefig(output_path, dpi=500)

plt.show()


# Volcano plot

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

sns.set()
sns.set_style("ticks")


data_file_path = 'path/DEGs.csv"'
gene_list_file_path_1 = 'path/DEGs_p0001.csv"'

data = pd.read_csv(data_file_path)
gene_list_1 = pd.read_csv(gene_list_file_path_1)

# Convert adjusted p-values to -log10 scale
data['is_in_list1'] = data['Gene'].isin(gene_list_1)
data['-log10(PValue)'] = -np.log10(data['padj'].replace(0, np.nan))


def plot_volcano(data):
    plt.figure(figsize=(8, 6), dpi=200)

    # Non-significant genes 
    data_non_sig = data[(data['padj'] >= 0.05) | data['padj'].isna()]
    plt.scatter(data_non_sig['log2FoldChange'], data_non_sig['-log10(PValue)'],
                alpha=1, c='gray', s=20, label="padj ≥ 0.05")

    # Upregulated genes
    data_up = data[(data['padj'] < 0.001) & (data['log2FoldChange'] > 0)]
    plt.scatter(data_up['log2FoldChange'], data_up['-log10(PValue)'],
                alpha=0.4, c='red', s=20, label="Upregulated (padj < 0.001)")

    # Downregulated genes
    data_down = data[(data['padj'] < 0.001) & (data['log2FoldChange'] < 0)]
    plt.scatter(data_down['log2FoldChange'], data_down['-log10(PValue)'],
                alpha=0.4, c='blue', s=20, label="Downregulated (padj < 0.001)")

    # Threshold lines
    plt.axhline(-np.log10(0.05), color='black', linestyle='--', label='P = 0.05 Threshold')
    plt.axvline(-1, color='gray', linestyle='--', label='Log2FC = -1 (FC=-2)')
    plt.axvline(1, color='gray', linestyle='--', label='Log2FC = 1 (FC=2)')

 
    plt.xlabel('Log2 Fold Change')
    plt.ylabel('-log10 p-value')

    plt.xlim(-8, 8)
    plt.yticks(np.arange(0, 200.1 + 0, 50), fontsize=10)
    plt.ylim(-2, 200)

    sns.despine()
    plt.legend(loc="upper right", fontsize=8)
    plt.tight_layout()
    plt.savefig("/.png", format="png", dpi=500)
    plt.show()


plot_volcano(data)




# Hierarchical clustering and cluster assignment

import seaborn as sns
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, cut_tree


df = pd.read_csv("/DEGs0001.csv")

df = df.pivot(index="Gene", columns="Label", values="log2_data")

order = ["Co_1", "Co_2", "Co_3", "HKCo_1", "HKCo_2", "HKCo_3", "Mono_1", "Mono_2", "Mono_3"]

clustered_df = sns.clustermap(df[order], cmap='cool', figsize=(20, 30), method='ward',
                              col_colors=['cyan', 'cyan', 'cyan', 'darkred', 'darkred', 'darkred', 'k', 'k', 'k'],
                              metric='euclidean',
                              cbar_pos=(0.98, 0.67, 0.03, 0.1))

row_order = clustered_df.dendrogram_row.reordered_ind
df_clustered = df.iloc[row_order]

num_clusters = 10
cluster_colors = sns.color_palette('cool_r', n_colors=num_clusters)
linkage_matrix = linkage(df_clustered, method='ward')
cluster_assignments = cut_tree(linkage_matrix, n_clusters=num_clusters).flatten()


cluster_info = pd.DataFrame({
    "Gene": df_clustered.index,
    "Cluster": cluster_assignments
})


cluster_info.to_csv("cluster.csv", index=False)


row_colors = [cluster_colors[cluster] for cluster in cluster_assignments]

def gen_cmap_name(cols):
    nmax = float(len(cols)-1)
    color_list = []
    for n, c in enumerate(cols):
        color_list.append((n/nmax, c))
    return mpl.colors.LinearSegmentedColormap.from_list('cmap', color_list)

cmap = gen_cmap_name(['blue','w','red'])

# Final clustered heatmap
g = sns.clustermap(df_clustered[order], cmap=cmap, figsize=(20, 30), method='ward',
                   col_colors=['khaki', 'khaki', 'khaki', 'palevioletred', 'palevioletred', 'palevioletred', 'lightskyblue', 'lightskyblue', 'lightskyblue'],
                   row_colors=row_colors, metric='euclidean', cbar_pos=(0.98, 0.67, 0.03, 0.1),vmax=2, vmin=-2)

for a in g.ax_col_dendrogram.collections:
    a.set_linewidth(2)

for a in g.ax_row_dendrogram.collections:
    a.set_linewidth(2)

plt.setp(g.ax_heatmap.get_yticklabels(), rotation=360, fontsize=8)


plt.savefig("/.png", format="png", dpi=600)


# Enrichment result frome DAVID (KEGG)
# Custom colormap for fold enrichment


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl

sns.set()
sns.set_style("ticks")


data = pd.read_csv('path .csv')

filtered_data = data[data['PValue'] < 0.01]


filtered_data['PValue'] = -np.log10(filtered_data['PValue'])

plt.figure(figsize=(2, 3), dpi=200)


y_positions = range(len(filtered_data))
y_labels = filtered_data['Term']




def gen_cmap_name(cols):
    nmax = float(len(cols)-1)
    color_list = []
    for n, c in enumerate(cols):
        color_list.append((n/nmax, c))

    return mpl.colors.LinearSegmentedColormap.from_list('cmap', color_list)
cmap = gen_cmap_name(['blue','lightgray','red'])

# Bubble size reflects gene count
sizes = (filtered_data['Count'] + 1) * 2


plt.axvline(x=-np.log10(0.05), color='grey', linestyle='--', linewidth=1)



sc = plt.scatter(filtered_data['PValue'], y_positions, s=sizes,
                 alpha=1, c=filtered_data['Fold Enrichment'], cmap=cmap, vmin=1, vmax=6)

plt.yticks(y_positions, y_labels, fontsize=10)

plt.xlabel('PValue', fontsize=10)


for count in [0, 5, 10, 15]:
    plt.scatter([], [], s=(count + 1) * 2, c='k', label=f'Count: {count}')



plt.tight_layout()
plt.gca().invert_yaxis()

sns.despine(right=True)
plt.tick_params(left=False)
plt.xticks(np.arange(0, 11, 2.5))

plt.xlim(-0.5,xx)
plt.ylim(xx,-0.5)


plt.savefig("/.png", format="png", dpi=100, bbox_inches='tight')
plt.show()





# Gene list frome KEGG Orthology (KO) - Sus scrofa (pig)
# https://www.kegg.jp/kegg-bin/get_htext?ssc00001+397656

import pandas as pd


input_csv = "KO_pig.csv"  
output_csv = "KO_pig.csv"  


df = pd.read_csv(input_csv)


def extract_name_from_list(list_value):
    if pd.isna(list_value):  
        return None
    try:

        name = list_value.split()[1] 
        return name.replace(";", "") 
    except IndexError:
        return None  


df['gene_name'] = df['List'].apply(extract_name_from_list)


df.to_csv(output_csv, index=False)

print(f"save: {output_csv}")



# e.g. "KO_Glycolysis_Gluconeogenesis"
import pandas as pd


df_a = pd.read_csv("KO_pig.csv")
df_a = df_a[df_a["Term"] == "Glycolysis_Gluconeogenesis"]

df_b = pd.read_csv("/TPM_Matrix_rere_FClog2.csv")


result_df = df_b[df_b['gene_name'].isin(df_a['gene_name'])]

result_df.to_csv("Glycolysis_Gluconeogenesis_FClog2.csv", index=False)


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl


data_file = "Glycolysis_Gluconeogenesis_FClog2.csv"
gene_list_file = "/Users/tatata/Desktop/matu-metab-re.csv"  

gene_list_df = pd.read_csv(gene_list_file)
gene_list_df = gene_list_df[gene_list_df["Fig"] == "Fig2"]
gene_list_df = gene_list_df[gene_list_df["List"] == 1]
target_genes = gene_list_df['gene_name'].dropna().str.strip().tolist()  


df = pd.read_csv(data_file)


order = ['Mono', 'Co', 'HKCo']
df['Class'] = pd.Categorical(df['Class'], categories=order, ordered=True)


heatmap_data = df.pivot_table(index='gene_name', columns='Class', values='log2FC')
heatmap_data = heatmap_data.reindex(target_genes)



def gen_cmap_name(cols):
    nmax = float(len(cols) - 1)
    color_list = [(n / nmax, c) for n, c in enumerate(cols)]
    return mpl.colors.LinearSegmentedColormap.from_list('cmap', color_list)

cmap = gen_cmap_name(['navy', 'lightgray', 'tomato'])


sns.set(font_scale=0.8)
sns.set_style("white")


plt.figure(figsize=(3, len(target_genes) * 0.8))
ax = sns.heatmap(
    heatmap_data,
    cmap=cmap,
    vmax=1,
    vmin=-1,
    cbar_kws={"label": "log2FC"},
    linewidths=0,
    linecolor='gray'
)

ax.set_xlabel("Class")
ax.set_ylabel("Gene")

plt.tight_layout()
plt.savefig(".png", format="png", dpi=500)
plt.show()




# For bacteria
# Principal component analysis

import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

file_path = "/TPM_Matrix.csv" 
df = pd.read_csv(file_path) 


df.set_index("Gene", inplace=True)

# Define experimental conditions
conditions = ["Co", "Mono"]  
columns = {cond: [col for col in df.columns if cond in col] for cond in conditions}


data = df.values.T 
samples = list(df.columns)  


scaler = StandardScaler()
scaled_data = scaler.fit_transform(data)
pca = PCA(n_components=2)  
pca_result = pca.fit_transform(scaled_data)


pca_df = pd.DataFrame(pca_result, columns=["PC1", "PC2"])
pca_df["Sample"] = samples
pca_df["Condition"] = ["Co" if "Co" in sample else "Mono" for sample in samples]  # 条件を識別

plt.figure(figsize=(8, 6))
colors = {"Co": "khaki", "Mono": "lightskyblue"}
for condition, color in colors.items():
    subset = pca_df[pca_df["Condition"] == condition]
    plt.scatter(subset["PC1"], subset["PC2"], label=condition, color=color, s=100, alpha=0.7)


plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.2f}%)")
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.2f}%)")

plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
plt.grid(False)
sns.despine()
plt.tight_layout()
plt.xlim(-80, 80)
plt.ylim(-40, 40)


output_path = ".png"
plt.savefig(output_path, dpi=500)
plt.show()



# Volcano plot

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

sns.set()
sns.set_style("ticks")


data_file_path = 'path/DEGs.csv"'
gene_list_file_path_1 = 'path/DEGs_p005.csv"'

data = pd.read_csv(data_file_path)
gene_list_1 = pd.read_csv(gene_list_file_path_1)

# Convert adjusted p-values to -log10 scale
data['is_in_list1'] = data['Gene'].isin(gene_list_1)
data['-log10(PValue)'] = -np.log10(data['padj'].replace(0, np.nan))


def plot_volcano(data):
    plt.figure(figsize=(8, 6), dpi=200)

    # Non-significant genes 
    data_non_sig = data[(data['padj'] >= 0.05) | data['padj'].isna()]
    plt.scatter(data_non_sig['log2FoldChange'], data_non_sig['-log10(PValue)'],
                alpha=1, c='gray', s=20, label="padj ≥ 0.05")

    # Upregulated genes
    data_up = data[(data['padj'] < 0.001) & (data['log2FoldChange'] > 0)]
    plt.scatter(data_up['log2FoldChange'], data_up['-log10(PValue)'],
                alpha=0.4, c='red', s=20, label="Upregulated (padj < 0.001)")

    # Downregulated genes
    data_down = data[(data['padj'] < 0.001) & (data['log2FoldChange'] < 0)]
    plt.scatter(data_down['log2FoldChange'], data_down['-log10(PValue)'],
                alpha=0.4, c='blue', s=20, label="Downregulated (padj < 0.001)")

    # Threshold lines
    plt.axhline(-np.log10(0.05), color='black', linestyle='--', label='P = 0.05 Threshold')
    plt.axvline(-1, color='gray', linestyle='--', label='Log2FC = -1 (FC=-2)')
    plt.axvline(1, color='gray', linestyle='--', label='Log2FC = 1 (FC=2)')

 
    plt.xlabel('Log2 Fold Change')
    plt.ylabel('-log10 p-value')

    plt.xlim(-8, 8)
    plt.yticks(np.arange(0, 200.1 + 0, 50), fontsize=10)
    plt.ylim(-2, 200)

    sns.despine()
    plt.legend(loc="upper right", fontsize=8)
    plt.tight_layout()
    plt.savefig("/.png", format="png", dpi=500)
    plt.show()


plot_volcano(data)




import pandas as pd

input_file = "TPM_Matrix.csv"  
output_file = "TPM_Matrix_re.csv"  


df = pd.read_csv(input_file)
gene_column = df['Gene']
sample_columns = ['Co_1', 'Co_2', 'Co_3', 'Mono_1', 'Mono_2', 'Mono_3']

sample_data = df[sample_columns]
labels = sample_columns * len(df)

class_labels = [
    'Co' if 'Co_' in col else 'Mono' for col in sample_columns
] * len(df)

reformatted_data = pd.DataFrame({
    'Gene': gene_column.repeat(len(sample_columns)).values,
    'Label': labels,
    'Class': class_labels,
    'Data': sample_data.values.flatten()
})


reformatted_data.to_csv(output_file, sep="\t", index=False)

print(f"saved {output_file}")


import pandas as pd


input_csv = "TPM_Matrix_re.csv"
output_csv = "TPM_Matrix_rere.csv"


df = pd.read_csv(input_csv, sep="\t")
df['Data1'] = df['Data'] + 1
df.to_csv(output_csv, index=False)


import pandas as pd


input_csv = "PM_Matrix_rere.csv"  
output_csv = "TPM_Matrix_rere_FC.csv"  


df = pd.read_csv(input_csv)


results = []

for gene, group in df.groupby("Gene"):

    mono_mean = group[(group["Class"] == "Mono") ]["Data1"].mean()
    

    group["Fold_change"] = group["Data1"] / mono_mean  # Fold change = Data1 / Mono
    results.append(group)


result_df = pd.concat(results)
result_df.to_csv(output_csv, index=False)



import pandas as pd
import numpy as np


input_csv = "TPM_Matrix_rere_FC.csv"  
output_csv = "TPM_Matrix_rere_FClog2.csv"  

df = pd.read_csv(input_csv)
df['log2FC'] = np.log2(df['Fold_change'])
df.to_csv(output_csv, index=False)



import pandas as pd
import glob

dataframes = []
file_paths = glob.glob("DEG_result_name.csv")

if not file_paths:
    print("Error1")
else:
    for file_path in file_paths:
        print(f"Processing: {file_path}")
        
        with open(file_path, 'r') as f:
            print("file:\n")
            for _ in range(5):
                print(f.readline())
        
        try:

            df = pd.read_csv(file_path, delimiter=',')  
            print(f"Columns: {df.columns}")
            

            df['gene_id'] = df['gene_id'].astype(str)
            

            if 'padj' in df.columns and 'gene_name' in df.columns:

                filtered_data = df[
                    (df['padj'] < 0.05) & df['gene_id'].str.startswith('HMPREF0531_', na=False)
                ]
                print(f"Filtered data shape: {filtered_data.shape}")
                
                if not filtered_data.empty:
  
                    dataframes.append(filtered_data.set_index('gene_id')[['gene_name']])
                else:
                    print("Error2")
            else:
                print("Error3")
        
        except Exception as e:
            print(f"Error reading file {file_path}: {e}")


if dataframes:
    merged_data = pd.concat(dataframes, axis=1)
    merged_data.to_csv("DEG_result005.csv", sep='\t')
    print("Merged")
else:
    print("Erro4")


    import pandas as pd


df_a = pd.read_csv("DEG_result005.csv")

df_b = pd.read_csv("TPM_Matrix_rere_FClog2.csv")


result_df = df_b[df_b['gene_name'].isin(df_a['gene_name'])]

result_df.to_csv("TPM_Matrix_rere_FClog2_DEGs005.csv", index=False)

# Hierarchical clustering and cluster assignment

import seaborn as sns
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, cut_tree


df = pd.read_csv("TPM_Matrix_rere_FClog2_DEGs005.csv")

df = df.pivot(index="gene_name", columns="Label", values="log2FC")

order = ["Co_1", "Co_2", "Co_3",  "Mono_1", "Mono_2", "Mono_3"]

clustered_df = sns.clustermap(df[order], cmap='cool', figsize=(20, 30), method='ward',
                              col_colors=['cyan', 'cyan', 'cyan', 'darkred', 'darkred', 'darkred', 'k', 'k', 'k'],
                              metric='euclidean',
                              cbar_pos=(0.98, 0.67, 0.03, 0.1))

row_order = clustered_df.dendrogram_row.reordered_ind
df_clustered = df.iloc[row_order]

num_clusters = 10
cluster_colors = sns.color_palette('cool_r', n_colors=num_clusters)
linkage_matrix = linkage(df_clustered, method='ward')
cluster_assignments = cut_tree(linkage_matrix, n_clusters=num_clusters).flatten()


cluster_info = pd.DataFrame({
    "Gene": df_clustered.index,
    "Cluster": cluster_assignments
})


cluster_info.to_csv("path _cluster_assignments.csv", index=False)


row_colors = [cluster_colors[cluster] for cluster in cluster_assignments]

def gen_cmap_name(cols):
    nmax = float(len(cols)-1)
    color_list = []
    for n, c in enumerate(cols):
        color_list.append((n/nmax, c))
    return mpl.colors.LinearSegmentedColormap.from_list('cmap', color_list)

cmap = gen_cmap_name(['midnightblue','mediumblue', 'rebeccapurple','tomato','yellow'])

g = sns.clustermap(df_clustered[order], cmap=cmap, figsize=(20, 30), method='ward',
                   col_colors=['khaki', 'khaki', 'khaki', 'lightskyblue', 'lightskyblue', 'lightskyblue'],
                   row_colors=row_colors, metric='euclidean', cbar_pos=(0.98, 0.67, 0.03, 0.1), vmin=-8, vmax=8)

for a in g.ax_col_dendrogram.collections:
    a.set_linewidth(2)

for a in g.ax_row_dendrogram.collections:
    a.set_linewidth(2)


plt.setp(g.ax_heatmap.get_yticklabels(), rotation=360, fontsize=8)


plt.savefig("path .png", format="png", dpi=600)



# Enrichment result frome DAVID (KEGG)
# Custom colormap for fold enrichment


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl

sns.set()
sns.set_style("ticks")


data = pd.read_csv('path .csv')

filtered_data = data[data['PValue'] < 0.01]


filtered_data['PValue'] = -np.log10(filtered_data['PValue'])

plt.figure(figsize=(2, 3), dpi=200)


y_positions = range(len(filtered_data))
y_labels = filtered_data['Term']




def gen_cmap_name(cols):
    nmax = float(len(cols)-1)
    color_list = []
    for n, c in enumerate(cols):
        color_list.append((n/nmax, c))

    return mpl.colors.LinearSegmentedColormap.from_list('cmap', color_list)
cmap = gen_cmap_name(['blue','lightgray','red'])

# Bubble size reflects gene count
sizes = (filtered_data['Count'] + 1) * 2


plt.axvline(x=-np.log10(0.05), color='grey', linestyle='--', linewidth=1)



sc = plt.scatter(filtered_data['PValue'], y_positions, s=sizes,
                 alpha=1, c=filtered_data['Fold Enrichment'], cmap=cmap, vmin=1, vmax=6)

plt.yticks(y_positions, y_labels, fontsize=10)

plt.xlabel('PValue', fontsize=10)


for count in [0, 5, 10, 15]:
    plt.scatter([], [], s=(count + 1) * 2, c='k', label=f'Count: {count}')



plt.tight_layout()
plt.gca().invert_yaxis()

sns.despine(right=True)
plt.tick_params(left=False)
plt.xticks(np.arange(0, 11, 2.5))

plt.xlim(-0.5,xx)
plt.ylim(xx,-0.5)


plt.savefig("/.png", format="png", dpi=100, bbox_inches='tight')
plt.show()



