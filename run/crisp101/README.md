This README file lists out the steps to run HybPiper on Crisp101 data.

- Generating target file (output: `crisp101_target_loci.fa`)
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
            --targetfile_dna ~/EucsPhylogenomics/run/crisp101/crisp101_target_loci.fa
            --bwa
            --outdir ~/outdir
            --namelist ~/namelist.txt
        ```

- Extracting `examples/sample.gene.bam`
    1. Sort the BAM output from HybPiper (`/04_processed_sample_directories/sample/sample.bam`)
        ```
        samtools sort sample.bam -o sample_sort.bam
        ```
    2. Index the BAM file
        ```
        samtools index sample_sort.bam
        ```
    3. Check the reference header for the gene (`/04_processed_sample_directories/sample/gene/gene_target.fasta`)
    4. Retrieve the BAM file for the gene
        ```
        samtools view -b sample_sort.bam gene_header > sample.gene.bam
        ```
    5. Other files can be retrieved from the HybPiper output

---
<i>Last updated: 24 April 2025 by Jeremias Ivan</i>