This README file lists out the steps to run HybPiper on Crisp101 data.

- crisp101_target_loci.fa
    1. Download the concatenated alignment (<i>.phy</i>) and partition file (<i>.nex</i>) from <a href="https://datadryad.org/dataset/doi:10.5061/dryad.gb5mkkwww">Dryad</a>
    2. Partition the alignment into individual locus alignments
    3. For each locus, remove all gaps from each species
    4. Update the sequence headers to be in the HybPiper format: `species-gene`
    5. Combine all locus sequences into one file

- Running HybPiper
    1. For each species, we run `hybpiper-nf` using the following command:
    ```
    nextflow run ~/hybpiper-nf/hybpiper.nf
        -c ~/hybpiper-nf/hybpiper.config
        -entry assemble
        -profile standard_singularity
        --illumina_reads_directory ~/dir_shortreads/
        --targetfile_dna ~/crisp101_target_loci.fa
        --bwa
        --outdir ~/outdir
        --namelist ~/namelist.txt
    ```