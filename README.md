# CellAnn shinyapp

### version 
1.0

### what is CellAnn ?

CellAnn is a web application for predicting cell types of single-cell clusters based on published reference datasets. CellAnn provides a comprehensive scRNA-seq reference database and users can easily find the relevant reference datasets in their analysis.

### installation 

1.Click on the green "Code" button on the right-hand side of this page. 

2.Click on "Download ZIP" from the drop-down menu. 

3.Your browser should start downloading a ZIP file containing the entire CellAnn repository. 

4.After downloading is finished, unzip the CellAnn_version1.0.zip* files 


### run CellAnn 

1. Please make sure you have installed the shiny package in R

2. cd to your folder where you unzip the CellAnn_version1.0.zip* files 

3. run the webserver on your local computer with the following commands:

```shell
cd your_unzip_folder
R -e "shiny::runApp('~/')"
```

### Note 
CellAnn is currently in beta. If you find a bug, please report an issue on the bug Report form of this github repository.
