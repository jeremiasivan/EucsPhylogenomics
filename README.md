# EucsPhylogenomics

**EucsPhylogenomics (Eucalypts Phylogenomics)** is an R pipeline to extract BUSCO loci from short-read data of Eucalypt samples. It comprises of two main steps: BUSCO extraction using CAPTUS and taxonomic consistency test across inferred gene trees. This pipeline is mainly developed and tested using MacOS and Linux, so there might be incompatibilities using Windows.

## Table of Content
- <a href="#prereqs">Prerequisites</a>
- <a href="#genpipe">General Pipeline</a>
- <a href="#refs">References</a>

## <a id="prereqs">Prerequisites</a>
EucsPhylogenomics requires a number of software and R packages to run. We recommend you to use environment management system (e.g. `conda`) to install the prerequisites, but you can also use `install.packages()` built-in function in R or RStudio.

### Software
|    Name        |                             Website                              |                             Anaconda                             |
| -------------- |:----------------------------------------------------------------:|:----------------------------------------------------------------:|
| AMAS           | <a href="https://github.com/marekborowiec/AMAS">Link</a>         | <a href="https://anaconda.org/bioconda/amas">Link</a>            |
| ASTER          | <a href="https://github.com/chaoszhang/ASTER">Link</a>           | <a href="https://anaconda.org/bioconda/aster">Link</a>           |
| Captus         | <a href="https://github.com/edgardomortiz/captus">Link</a>       | <a href="https://anaconda.org/bioconda/captus">Link</a>          |
| IQ-TREE        | <a href="http://www.iqtree.org">Link</a>                         | <a href="https://anaconda.org/bioconda/iqtree">Link</a>          |
| MAFFT          | <a href="https://github.com/GSLBiotech/mafft">Link</a>           | <a href="https://anaconda.org/bioconda/mafft">Link</a>           |
| TreeShrink     | <a href="https://github.com/uym2/TreeShrink">Link</a>            | <a href="https://anaconda.org/bioconda/treeshrink">Link</a>      |
| VeryFastTree   | <a href="https://github.com/citiususc/veryfasttree">Link</a>     | <a href="https://anaconda.org/bioconda/veryfasttree">Link</a>    |

### R packages
|    Name    |                               CRAN                               |                                   Anaconda                               |
| ---------- |:----------------------------------------------------------------:|:------------------------------------------------------------------------:|
| ape        | <a href="https://cran.r-project.org/package=ape">Link</a>        | <a href="https://anaconda.org/conda-forge/r-ape">Link</a>                |
| Biostrings | <a href="https://bioconductor.org/packages/Biostrings">Link</a>  | <a href="https://anaconda.org/bioconda/bioconductor-biostrings">Link</a> |
| data.table | <a href="https://cran.r-project.org/package=data.table">Link</a> | <a href="https://anaconda.org/conda-forge/r-data.table">Link</a>         |
| doSNOW     | <a href="https://cran.r-project.org/package=doSNOW">Link</a>     | <a href="https://anaconda.org/conda-forge/r-dosnow">Link</a>             |
| kableExtra | <a href="https://cran.r-project.org/package=kableExtra">Link</a> | <a href="https://anaconda.org/conda-forge/r-kableextra">Link</a> |
| log4r      | <a href="https://cran.r-project.org/package=log4r">Link</a>      | <a href="https://anaconda.org/conda-forge/r-log4r">Link</a>              |
| optparse   | <a href="https://cran.r-project.org/package=optparse">Link</a>   | <a href="https://anaconda.org/conda-forge/r-optparse">Link</a>           |
| rmarkdown  | <a href="https://cran.r-project.org/package=rmarkdown">Link</a>  | <a href="https://anaconda.org/conda-forge/r-rmarkdown">Link</a>          |
| tidyverse  | <a href="https://cran.r-project.org/package=tidyverse">Link</a>  | <a href="https://anaconda.org/conda-forge/r-tidyverse">Link</a>          |
| yaml       | <a href="https://cran.r-project.org/package=yaml">Link</a>       | <a href="https://anaconda.org/conda-forge/r-yaml">Link</a>               |

### [Optional] Additional Software for Phasing
|    Name        |                             Website                              |                             Anaconda                             |
| -------------- |:----------------------------------------------------------------:|:----------------------------------------------------------------:|
| BCFtools       | <a href="https://github.com/samtools/bcftools">Link</a>          | <a href="https://anaconda.org/bioconda/bcftools">Link</a>        |
| BWA            | <a href="https://github.com/lh3/BWA">Link</a>                    | <a href="https://anaconda.org/bioconda/bwa">Link</a>             |
| SAMtools       | <a href="https://github.com/samtools/samtools">Link</a>          | <a href="https://anaconda.org/bioconda/samtools">Link</a>        |
| WhatsHap       | <a href="https://github.com/whatshap/whatshap">Link</a>          | <a href="https://anaconda.org/bioconda/whatshap">Link</a>        |

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
        conda install -c conda-forge -c bioconda r-ape r-data.table r-doSNOW r-kableExtra r-log4r r-optparse r-rmarkdown r-tidyverse r-yaml bioconductor-biostrings amas aster captus iqtree mafft treeshrink veryfasttree
        ```
    - [Optional] Installing software for phasing
        ```
        conda install -c conda-forge -c bioconda bcftools bwa samtools whatshap
        ```

    If `bcftools` returns a `libgsl.so.25` error, you can either download the software <a href="https://www.htslib.org/download/">here</a>, or try to set the `conda` channel priorities before installing any package:
    ```
    conda config --prepend channels r
    conda config --prepend channels bioconda
    conda config --prepend channels conda-forge
    conda config --set channel_priority strict
    ```

3. **Update the parameters in `config.yaml`** <br>

4. **Run EucsPhylogenomics** <br>
    ```
    Rscript run_pipeline.R --config config.yaml
    Rscript run_pipeline.R --config config.yaml --redo
    ```

    In UNIX-based operating systems (e.g., Linux and MacOS), it is advisable to use `nohup` or `tmux` to run the whole pipeline. For Windows, you can use `psmux`. 

---
## <a id="refs">References</a>
1. Borowiec, M.L. (<a href="https://doi.org/10.7717/peerj.1660">2016</a>). **AMAS: A fast tool for alignment manipulation and computing of summary statistics**. *PeerJ*, *4*, e1660.

2. Zhang, C., et al. (<a href="https://doi.org/10.1093/molbev/msaf172">2025</a>). **ASTER: A package for large-scale phylogenomic reconstructions**. *Molecular Biology and Evolution*, *42*(8), msaf172.

3. Danecek, P., et al. (<a href="https://doi.org/10.1093/gigascience/giab008">2021</a>). **Twelve years of SAMtools and BCFtools**. *GigaScience*, *10*(2), giab008.

4. Li, H. & Durbin, R. (<a href="https://doi.org/10.1093/bioinformatics/btp324">2009</a>). **Fast and accurate short read alignment with Burrows–Wheeler transform**. *Bioinformatics*, *25*(14), 1754–1760.

5. Ortiz, E.M., et al. (<a href="https://doi.org/10.1101/2023.10.27.564367">2023</a>). **A novel phylogenomics pipeline reveals complex pattern of reticulate evolution in Cucurbitales**. *bioRxiv*.

6. Minh, B.Q., et al. (<a href="https://doi.org/10.1093/molbev/msaa015">2020</a>). **IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era**. *Molecular Biology and Evolution*, *37*(5), 1530–1534.

7. Katoh, K. & Standley, D.M. (<a href="https://doi.org/10.1093/molbev/mst010">2013</a>). **MAFFT multiple sequence alignment software version 7: Improvements in performance and usability**. *Molecular Biology and Evolution*, *30*(4), 772–780.

8. Mai, U. & Mirarab, S. (<a href="https://doi.org/10.1186/s12864-018-4620-2">2018</a>). **TreeShrink: fast and accurate detection of outlier long branches in collections of phylogenetic trees**. *BMC Genomics*, *19*(272).

9. Piñeiro, C., et al. (<a href="https://doi.org/10.1093/bioinformatics/btaa582">2020</a>). **Very Fast Tree: speeding up the estimation of phylogenies for large alignments through parallelization and vectorization strategies**. *Bioinformatics*, *36*(17), 4658-4659.

10. Martin, M., et al. (<a href="https://doi.org/10.1101/085050">2016</a>). **WhatsHap: Fast and accurate read-based phasing**. *bioRxiv*.

11. Paradis, E., et al. (<a href="https://doi.org/10.1093/bioinformatics/btg412">2004</a>). **APE: Analyses of Phylogenetics and Evolution in R language**. *Bioinformatics*, *20*(2), 289-290.

12. Pagès, H., et al. (<a href="https://doi.org/10.18129/B9.bioc.Biostrings">2026</a>). **Biostrings: Efficient manipulation of biological strings**. *R package*.

13. Barrett, T., et al. (<a href="https://doi.org/10.32614/CRAN.package.data.table">2026</a>). **data.table: Extension of 'data.frame'**. *R package*.

14. Daniel, F. (<a href="https://cran.r-project.org/package=doSNOW">2022</a>). **doSNOW: Foreach Parallel Adaptor for the 'snow' Package**. *R package*.

15. Zhu, H., et al. (<a href="https://cran.r-project.org/package=kableExtra">2024</a>). **kableExtra: Construct Complex Table with 'kable' and Pipe Syntax**. *R package*.

16. White, J.M. & Jacobs, A. (<a href="https://doi.org/10.32614/CRAN.package.log4r">2024</a>). **log4r: A Fast and Lightweight Logging System for R, Based on 'log4j'**. *R package*.

17. Davis, T.L. (<a href="https://doi.org/10.32614/CRAN.package.optparse">2026</a>). **optparse: Command Line Option Parser**. *R package*.

18. Allaire, J.J., et al. (<a href="https://doi.org/10.32614/CRAN.package.rmarkdown">2026</a>). **rmarkdown: Dynamic Documents for R**. *R package*.

19. Wickham, H., et al. (<a href="https://doi.org/10.21105/joss.01686">2019</a>). **Welcome to the tidyverse**. *Journal of Open Source Software*, *4*(43), 1686.

20. Stephens, J., et al. (<a href="https://doi.org/10.32614/CRAN.package.yaml">2025</a>). **yaml: Methods to Convert R Data to YAML and Back**. *R package*.

21. Anthropic. (<a href="https://claude.ai/">2026</a>). Claude 4.6 Sonnet was used to generate `config.yaml` and `run_pipeline.R`. 

---
*Last update: 29 June 2026 by Jeremias Ivan*