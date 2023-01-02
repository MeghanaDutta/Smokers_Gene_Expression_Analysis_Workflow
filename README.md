# Smokers_Analysis_Pipeline

Cigarette smoking is the leading cause of preventable death in USA. 
It causes nearly 1 in 5 deaths. It causes 90% of all lung cancer deaths. 
More women die from lung cancer each year than from breast cancer. 
Cigarette smoke contains many chemicals which are harmful to both smokers and non-smokers 
It creates a field of injury in epithelial cells lining the respiratory tract. 
This field of injury extends from intrathoracic epithelia to extrathoracic epithelia that line the nose.
To understand the changes in the epithelial cells of smokers I used dataset GDS3309 from GEO for analysis. 
The dataset GDS3309 is available in GEO dataset browser. 
GEO stands for Gene Expression Omnibus. It is a public database for microarray datasets and is hosted by NCBI. 
This dataset has 15 samples containing genes of epithelial cells from humans.
The experimental design of the dataset is very simple. It has only one factor, smoking. 
It has 7 smoker samples and 8 non-smoker samples which will become the controls for this analysis. 
I obtained the dataset directly from the GEO library which is available in Bioconductor. 
I performed various analysis on the dataset to answer many questions regarding the gene expression data. 
