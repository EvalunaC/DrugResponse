# DrugResponse

## Environment
**pRRophetic**
conda install -c conda-forge r-biocmanager

*BiocManager::install(version = "3.13")*
BiocManager::install()
source("http://bioconductor.org/biocLite.R")
*biocLite(c("car", "ridge", "preprocessCore", "genefilter", "sva")) Use following instead*
wget -O pRRophetic_0.5.tar.gz https://osf.io/dwzce/?action=download
R CMD INSTALL pRRophetic_0.5.tar.gz


conda install -c r r-car
conda install -c r r-ridge
conda install -c bioconda bioconductor-preprocesscore
conda install -c bioconda bioconductor-genefilter
conda install -c bioconda bioconductor-sva
conda install -c r-mgcv
conda install r-dplyr
conda install -c bioconda bioconductor-RDRToolbox
conda install r-kernlab
conda install r-mice
conda install -c conda-forge r-performanceanalytics

conda install -c r r-caret

conda install -c anaconda libopenblas

conda install -c bioconda bioconductor-tcgabiolinks
conda install -c conda-forge r-pdftools 
