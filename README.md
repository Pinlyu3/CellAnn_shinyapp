# CellAnn shinyapp

### Version 
1.0

### What is CellAnn ?

CellAnn is a web application for predicting cell types of single-cell clusters based on published reference datasets. CellAnn provides a comprehensive scRNA-seq reference database and users can easily find the relevant reference datasets in their analysis.


## Prepare query datasets 
### for R seurat users:
1.Install devtools (for installing GitHub packages) if it isn't already installed.
```{r}
### install devtools in R ###
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
```
2.Prepare your Seurat object with "prepare_CellAnn" function.

```{r}
### download scripts ###

devtools::source_url("https://raw.githubusercontent.com/Pinlyu3/CellAnn/main/prepare_CellAnn.R")

### prepare_CellAnn function:
### parameter: seurat_obj: your seurat obj
### parameter: folder(character): CellAnn ouput files will be output to this path, default is your current folder
### parameter: sample_name(character): names of this sample, such as "Liver_1", "eye2_2" etc.
### parameter: matrix_name(character): default will use the "RNA" counts matrix in your seurat object
### parameter: dims(character): "umap" or "tsne". default is 'umap'
### parameter: cluster(character): the name of the column which stored cluster information in the metadata of your Seurat object. default is 'seurat_clusters'

prepare_CellAnn(seurat_obj,folder=folder,sample_name='your_samples_name',matrix_name='RNA',dims='umap',cluster='seurat_clusters')

### After run prepare_CellAnn function, you will find 2 prepared files under your folder 

list.files(folder)
```

### for Python scanpy users:
```python
### load library and your query datasets: 
import scanpy as sc
import pandas as pd

### Calculate the average gene expression within each cluster
### we assume that the raw counts are stored in adata.raw.X
cluster_avg = pd.DataFrame(adata.raw.X.mean(axis=0), index=adata.raw.var_names, columns=['Average Expression'])

cluster_avg.index.name = 'Gene'

### write the step1 input (average expression matrix) to txt:
cluster_avg.to_csv('cluster_avg.txt', sep='\t')

### write the step4 input (UMAP for each cell) to txt:


```
## run CellAnn on your local computer
### Installation 

1. Make sure that Git is already installed.

2. Install git-lfs.

```shell
cd your_folder
git clone https://github.com/Pinlyu3/CellAnn_shinyapp.git
cd CellAnn_shinyapp
cat  CellAnn_version1.0.zip* > CellAnn_version1.0.zip
unzip CellAnn_version1.0.zip
cd CellAnn_version1\ copy
```

### Run CellAnn 

1. Please make sure you have installed the shiny package in R

2. cd to your folder where you unzip the CellAnn_version1.0.zip* files (cd CellAnn_version1\ copy)

3. run the webserver on your local computer with the following commands:

```shell
cd your_unzip_folder
R -e "shiny::runApp('~/')"
```

### Note 
CellAnn is currently in beta. If you find a bug, please report an issue on the bug Report form of this github repository.
