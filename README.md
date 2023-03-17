# CellAnn shinyapp

### Version 
1.0

### What is CellAnn ?

CellAnn is a web application for predicting cell types of single-cell clusters based on published reference datasets. CellAnn provides a comprehensive scRNA-seq reference database and users can easily find the relevant reference datasets in their analysis.

### Installation 

1. Make sure that Git is already installed.

2. Install git-lfs.

```shell
cd your_folder
git clone https://github.com/Pinlyu3/CellAnn_shinyapp.git
cd CellAnn_shinyapp

```

### Run CellAnn 

1. Please make sure you have installed the shiny package in R

2. cd to your folder where you unzip the CellAnn_version1.0.zip* files 

3. run the webserver on your local computer with the following commands:

```shell
cd your_unzip_folder
R -e "shiny::runApp('~/')"
```

### Note 
CellAnn is currently in beta. If you find a bug, please report an issue on the bug Report form of this github repository.
