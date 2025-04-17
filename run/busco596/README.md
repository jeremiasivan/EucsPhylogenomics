This README file lists out the steps to run HybPiper on Crisp101 data.

- busco596_target_loci.fa
    1. Download the genomes assemblies from <a href="https://www.ncbi.nlm.nih.gov/bioproject/509734">NCBI</a>
    2. Run BUSCO pipeline on each assembly
    3. Extract all single-copy, complete BUSCO that are shared between assemblies
    4. Store the BUSCO amino acid sequences in one file

- Running HybPiper
    1. For each species, we run `hybpiper-nf` using the following command:
    ```
    nextflow run ~/hybpiper-nf/hybpiper.nf
        -c ~/hybpiper-nf/hybpiper.config
        -entry assemble
        -profile standard_singularity
        --illumina_reads_directory ~/dir_shortreads/
        --targetfile_aa ~/EucsPhylogenomics/run/crisp101/busco596_target_loci.fa
        --diamond
        --outdir ~/outdir
        --namelist ~/namelist.txt
    ```

---
<i>Last updated: 17 April 2025 by Jeremias Ivan</i>