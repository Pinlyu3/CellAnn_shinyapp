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
adata = sc.read('data.h5ad')

### Calculate the average gene expression within each cluster
### we assume that the raw counts are stored in adata.raw.X
### the cluster information stored in adata.obs['cluster']: 
### the gene_names stored in adata.var_names
### the cell name stored in adata.obs_names
### the 2 dimension UMAPs strored in adata.obsm['X_umap']
### parameter: sample_name(character): names of this sample, such as "Liver_1"

def cluster_average(adata,sample_name):
  cluster_df = adata.obs['cluster']
  gene_names = adata.var_names
  cell_names = adata.obs_names
  mat_df=pd.DataFrame.sparse.from_spmatrix(adata.raw.X)
  ####
  row_annotations = np.array(cluster_df)
  unique_annotations = np.unique(row_annotations)
  ####
  merged_matrix = np.zeros((len(unique_annotations), mat_df.shape[1]))
  ####
  mat_df_mat = np.array(mat_df)
  for i, annotation in enumerate(unique_annotations):
    ###
    rows_to_sum = (row_annotations == annotation)
    rows_to_sum_mat = mat_df_mat[rows_to_sum,:]
    ###
    sum_rows = np.sum(rows_to_sum_mat, axis=0)
    sum_rows_norm = sum_rows / np.sum(sum_rows) * 1e5
    sum_rows_log = np.log(sum_rows_norm+1)
    ###
    merged_matrix[i, :] = sum_rows_log
  #### output to csv #####
  df_log = pd.DataFrame(merged_matrix.T)
  df_log.columns = unique_annotations
  df_log.index = gene_names
  df_log.insert(0,"GENE",df_log.index)
  ####
  df_log = df_log.round(3)
  Output_Step1_name = sample_name + '_CellAnn_Step1_input.txt'
  df_log.to_csv(Output_Step1_name, index=False,sep='\t')
  ####
  df_umap = pd.DataFrame(adata.obsm['X_umap'])
  df_umap.columns = ["dim1","dim2"]
  df_umap.insert(0,"cell",cell_names)
  df_umap.insert(1,"cluster",row_annotations)
  Output_Step4_name = sample_name + '_CellAnn_Step4_input.txt'
  df_umap.to_csv(Output_Step4_name, index=False,sep='\t')

### run cluster_average with your adata to get the input of step1 and step4:
cluster_average(adata,sample_name="your_sample_name")

```
## run CellAnn on your local computer

1. Download the CellAnn_shinyapp-main.zip file and unzip.

2. cd to your unzip folder (CellAnn_shinyapp-main) and unzip the CellAnn_version1.0.zip file.

3. cd to your unzip folder (Cell_Ann_bioinfo) and run the webserver on your local computer with the following commands:

```shell
cd /your_unzip_path/Cell_Ann_bioinfo
R -e "shiny::runApp('~/')"
```

### Note 
CellAnn is currently in beta. If you find a bug, please report an issue on the bug Report form of this github repository.
