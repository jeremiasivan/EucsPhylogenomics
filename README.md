# EucsPhylogenomics

**EucsPhylogenomics (Eucalypts Phylogenomics)** is an R pipeline to extract BUSCO loci from short-read data of Eucalypt samples. This pipeline is mainly developed and tested using MacOS and Linux, so there might be incompatibilities using Windows.

## Table of Content
- <a href="#genpipe">General Pipeline</a>

## <a id="genpipe">General Pipeline</a>
1. **Clone the Git repository** <br>
    ```
    git clone git@github.com:jeremiasivan/EucsPhylogenomics.git
    ```

2. **Install the prerequisites** <br>
    - Create a new conda environment
        ```
        conda create -n eucsphylo
        conda activate eucsphylo
        ```
    - Installing prerequisites
        ```
        conda install r-ape r-data.table r-doSNOW r-log4r r-optparse r-rmarkdown r-tidyverse r-yaml bioconductor-biostrings aster bcftools bwa captus fasttree iqtree mafft samtools whatshap
        ```
    
    If `bcftools` returns a `libgsl.so.25` <a href="https://github.com/samtools/bcftools/issues/1698">error</a>, you can try to set the `conda` channel priorities before installing any package:
    ```
    conda config --prepend channels r
    conda config --prepend channels bioconda
    conda config --prepend channels conda-forge
    conda config --set channel_priority strict
    ```

---
*Last update: 20 June 2026 by Jeremias Ivan*