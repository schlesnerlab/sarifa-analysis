FROM bioconductor/bioconductor_docker:RELEASE_3_15

#RUN apt-get update && apt-get install -y --no-install-recommends libcurl4-gnutls-dev

#RUN mkdir -p /cwd/install
 
RUN R -e 'options(repos = list(CRAN = "https://packages.othr.de/cran/"))'
RUN R -e 'BiocManager::install("DESeq2")'
RUN R -e 'BiocManager::install("BioinformaticsFMRP/TCGAbiolinksGUI.data")'
RUN R -e 'BiocManager::install("BioinformaticsFMRP/TCGAbiolinks")'
RUN R -e 'BiocManager::install("EnsDb.Hsapiens.v86")'
RUN R -e 'BiocManager::install("EnhancedVolcano")'
RUN R -e 'BiocManager::install("pcaExplorer")'
RUN R -e 'BiocManager::install("PCAtools")'
RUN R -e 'BiocManager::install("ComplexHeatmap")'
RUN R -e 'BiocManager::install("clusterProfiler")'

RUN R -e 'install.packages("readxl")'
RUN R -e 'install.packages("ggplot2")'
RUN R -e 'install.packages("memoise")'
RUN R -e 'install.packages("cachem")'
RUN R -e 'install.packages("sjmisc")'
RUN R -e 'install.packages("writexl")'
RUN R -e 'install.packages("survminer")'
RUN R -e 'install.packages("ggvenn")'
RUN R -e 'install.packages("msigdbr")'
RUN R -e 'install.packages("dplyr")'
RUN R -e 'install.packages("cachem")'
RUN R -e 'install.packages("ggupset")'
RUN R -e 'install.packages("cowplot")'
RUN R -e 'install.packages("rmarkdown")'
RUN R -e 'install.packages("magrittr")'
RUN R -e 'install.packages("tibble")'
RUN R -e 'install.packages("ggsignif")'
RUN R -e 'install.packages("forcats")'
RUN R -e 'install.packages("corrr")'
RUN R -e 'install.packages("svglite")'


RUN mkdir -p /cwd/code
RUN mkdir -p /cwd/data
RUN mkdir -p /cwd/data/pca
RUN mkdir -p /cwd/data/pca_mg
RUN mkdir -p /cwd/data/plots

CMD /bin/bash 
# CMD Rscript /rnaseq/code/main.R