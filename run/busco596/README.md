This README file lists out the steps to run HybPiper on BUSCO data.

## Analyses
- Generating target file (output: `busco596_target_loci.fa`)
    1. Download the 36 genomes assemblies from <a href="https://www.ncbi.nlm.nih.gov/bioproject/509734">NCBI</a>
    2. Run BUSCO pipeline on each assembly
    3. Extract all single-copy, complete BUSCO that are shared between assemblies
    4. Store the BUSCO amino acid sequences in one file

- Running HybPiper
    1. For each sample, run <a href="https://github.com/chrisjackson-pellicle/hybpiper-nf">hybpiper-nf</a> using the following command:
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

    2. For each sample, exclude BUSCO loci that are flagged as paralogs by `hybpiper-nf`
    3. Combine the BUSCO amino acid sequences from 771 Currency Creek samples and the 36 genome assemblies

- Sequence filtering (per BUSCO locus)
    1. Convert all `X` to `-` (gaps)
    2. Remove samples with >=50% gaps (output: `busco.rm1_50gaps.faa`)
        ```
        # iterate over loci
        foreach (gene = ls_genes) %dopar% {
            fn_input <- paste0("./", gene, "/", gene, ".faa")

            # open locus sequence
            seq <- Biostrings::readAAStringSet(fn_input)

            # iterate over samples
            pl <- c()
            for (i in 1:length(seq)) {
                # convert 'X' to '-'
                seq[i] <- gsub("X","-",seq[i])
                seq_chr <- as.character(seq[i])

                # store samples with >=50% gaps
                if (stringr::str_count(seq_chr,"-")>=0.5*stringr::str_count(seq_chr)) {
                    pl <- c(pl, i)
                }
            }

            # store the list of samples with >=50% gaps
            writeLines(names(seq)[pl], con=paste0("./", gene, "/", gene, ".rm1_50gaps.txt"))

            # remove samples with >=50% gaps
            if (length(pl)>0) {
                seq <- seq[-pl]
            }
            Biostrings::writeXStringSet(seq, filepath=paste0("./", gene, "/", gene, ".rm1_50gaps.faa"))
        }
        ```
    3. Align sequences using <a href="https://mafft.cbrc.jp/alignment/software/">MAFFT</a> FFT-NS-i with 1000 iterations (output: `busco.rm1_aligned.faa`)
    4. Remove sites with >90% gaps (output: `busco.rm2_90gaps.faa`)
        ```
        # iterate over loci
        foreach (gene = ls_genes) %dopar% {
            # open the alignment
            fn_input <- paste0("./", gene, "/", gene, ".rm1_aligned.faa")
            msa <- Biostrings::readAAStringSet(fn_input, format="fasta")

            # convert alignment to character matrix
            seqs <- as.character(msa)
            msa_matrix <- do.call(rbind, strsplit(seqs, ""))

            # count number of sequences
            n_seqs <- nrow(msa_matrix)

            # function to check the gap proportion per column
            is_valid_site <- function(column) {
                gap_count <- sum(column == "-")
                return((gap_count / n_seqs) <= 0.9)
            }

            # calculate the gap proportion of each column
            valid_sites <- apply(msa_matrix, 2, is_valid_site)

            # remove columns with >=90% gaps
            filtered_matrix <- msa_matrix[, valid_sites]
            filtered_alignment <- Biostrings::AAStringSet(apply(filtered_matrix, 1, paste0, collapse=""))

            # save the new alignment
            Biostrings::writeXStringSet(filtered_alignment, filepath=paste0("./", gene, "/", gene, ".rm2_90gaps.faa"))
        }
        ```
    5. Mask erroneous sites with <a href="https://github.com/chaoszhang/TAPER">TAPER</a> (output: `busco.rm3_taper.faa`)
    6. Convert all `X` to `-` (gaps)
    7. Remove samples with >=50% gaps (output: `busco.rm4_50gaps.faa`)
    8. Remove sites with >90% gaps (output: `busco.rm5_90gaps.faa`)

- Building trees
    1. For each locus, build approximately-ML tree using <a href="https://morgannprice.github.io/fasttree/">FastTree</a> (output: `busco.rm5.treefile`)
    2. For each locus tree, run <a href="https://github.com/uym2/TreeShrink">TreeShrink</a> (output: `busco.rm5_treeshrink.treefile`)
    3. Add the suffix `-ref` for the 36 genome assemblies on the locus tree (output: `busco.rm5_treeshrink.treefile.updated`)
    5. Combine all locus trees into one file (output: `species_tree/all.rm5_treeshrink.treefile`)
    6. Run <a href="https://github.com/smirarab/ASTRAL">ASTRAL-III</a> (output: `species_tree/astral.rm5_treeshrink.treefile`)

- Summarising the sequence filtering steps (output: `summary_filtering.tsv`)
    | Column Name               | Description                                                    |
    | ------------------------- |:---------------------------------------------------------------|
    | `busco`                   | Name of the BUSCO locus                                        |
    | `init_sp`                 | Initial number of samples after removing paralogs              |
    | `init_min_len`            | Shortest sequence from the initial locus sequences             |
    | `init_max_len`            | Longest sequence from the initial locus sequences              |
    | `rm1_50gaps_sp`           | Number of samples after removing samples with >=50% gaps       |
    | `rm1_aligned_len`         | Length of the initial locus alignment                          |
    | `rm1_aligned_gappercent`  | Percentage of gaps of the initial locus alignment              |
    | `rm2_90gaps_sp`           | Number of samples after removing sites with >90% gaps          |
    | `rm2_90gaps_len`          | Length of the alignment after removing sites with >90% gaps    |
    | `rm2_90gaps_gappercent`   | Percentage of gaps after removing sites with >90% gaps         |
    | `rm3_taper_xpercent`      | Percentage of added 'X' on the alignment                       |
    | `rm4_50gaps_sp`           | Number of samples after removing samples with >=50% gaps       |
    | `rm4_50gaps_len`          | Length of the alignment after removing samples with >=50% gaps |
    | `rm4_50gaps_gappercent`   | Percentage of gaps after removing samples with >=50% gaps      |
    | `rm5_90gaps_sp`           | Number of samples after removing sites with >90% gaps          |
    | `rm5_90gaps_len`          | Length of the alignment after removing sites with >90% gaps    |
    | `rm5_90gaps_gappercent`   | Percentage of gaps after removing sites with >90% gaps         |
    | `treeshrink_sp`           | Number of species after TreeShrink                             |
    | `treeshrink_len`          | Length of the alignment after TreeShrink                       |
    | `len_percent`             | Percentage of `treeshrink_len` compared to `rm1_aligned_len`   |
    | `sp_percent`              | Percentage of `treeshrink_sp` compared to `rm1_50gaps_sp`      |

## Additional Files
- `busco.rm1_50gaps.txt`: list of removed samples from the 2nd step of <i>Sequence Filtering</i>
- `busco.rm4_50gaps.txt`: list of removed samples from the 7th step of <i>Sequence Filtering</i>
- `busco.rm5_treeshrink.txt` and `busco.rm5_treeshrink_summary.txt`: log files from `TreeShrink`
- `species_tree/astral.rm5_treeshrink.log`: log file from ASTRAL-III

---
<i>Last updated: 22 April 2025 by Jeremias Ivan</i>