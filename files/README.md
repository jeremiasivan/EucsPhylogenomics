This README file lists out the steps to retrieve the metadata files.

## `currency_creek/`
- `jgi_metadata.csv`
    1. Go to <a href="https://genome.jgi.doe.gov/portal/">JGI Genome Portal</a>
    2. Search for <i>"Currency Creek"</i>
    3. Click on <b><i>Reports</i></b> -> <b><i>Download Project Overview Report</i></b>

- `justin_metadata.csv`
    1. Go to this <a href="https://raw.githubusercontent.com/borevitzlab/cca-eucs/master/metadata/CCATreesCombinedMetadata.csv">link</a>
    2. Download the file (e.g., using `wget`) 

- `combined_metadata.tsv`
    1. Open `jgi_metadata.csv` in R
    2. Update column names with lowercases and underscores
    3. Retrieve <i>field_id</i> from <i>project_name</i> column using:
    ```
    df_jgi$field_id <- sapply(df_jgi$project_name, function(x) {
        ls <- unlist(strsplit(x, split=" "));
        ls[grepl("^CCA", ls)]
    })
    ```
    4. Retrieve <i>folder_name</i> from <i>portal_id</i> column using:
    ```
    df_jgi$folder_name <- sapply(df_jgi$portal_id, function(x) {
        unlist(strsplit(x, split='"'))[7]
    })
    ```
    5. Open `justin_metadata.csv`
    6. Update column names with lowercases and underscores
    7. Merge the two data.frames using:
    ```
    df_combined <- df_jgi %>% left_join(df_justin, by="field_id", keep = FALSE)
    ```
    8. Remove duplicated entry for `field_id == "CCA0154"` as the project name is <i>E. corrugata</i> but identified as <i>E. crucis subsp. crucis</i>
    9. Create the `tip_label` column based on the `species_name` column with iterator
    10. Save the data.frame using `data.table::fwrite()`

## `all/`
- `seqkit_stats.tsv`
    1. Separate the forward and reverse reads from the filtered FASTQ files using `reformat.sh` from `BBTools`
    2. Store the forward and reverse FASTQ files for each sample in one folder
    3. Run `seqkit stats -a *.fastq.gz > seqkit_stats.tsv`

    > **Note**: `E_effusa_exsul_1` was excluded as it was unavailable at the time of download (01 October 2024)

- `seqkit_cov.tsv`
    1. Open `seqkit_stats.tsv`
    2. For each sample, add up the `sum_len` column between forward and reverse reads
    3. Create the `cov` column by dividing the `sum_len` column with $5\cdot10^8$ (i.e., the average genome size of Eucalypts) 

- `tree_taxonomy_label.tsv` <br>
    This file extracts the taxonomic groupings of 837 Eucalypts from `currency_creek/combined_metadata.tsv` and Ferguson et al. (<a href="https://doi.org/10.1101/gr.277999.123">2024</a>)

## `hybpiper/`
- `crisp101_eucs.fna` and `crisp101_noneucs.fna`
    1. Download the concatenated alignment (<i>.phy</i>) and partition file (<i>.nex</i>) from <a href="https://datadryad.org/dataset/doi:10.5061/dryad.gb5mkkwww">Dryad</a>
    2. Partition the concatenated alignment into individual loci (e.g., using <a href="https://github.com/marekborowiec/AMAS">AMAS</a>)
    3. For each locus, remove all gaps from each sample
    4. Update the sequence headers to be in HybPiper format: `species-gene`
    5. Combine all locus sequences into two files: one for <i>Eucalyptus</i> species, and one for <i>Angophora</i> and <i>Corymbia</i> species

- `busco1187_eucs.fna` and `busco1187_noneucs.fna`
    1. Run `BUSCO v5.8.0` on 36 reference genomes from Ferguson et al. (<a href="https://doi.org/10.1101/gr.277999.123">2024</a>)
    2. Extract all single-copy, complete BUSCO from each reference genome
    3. For each BUSCO locus:
        - Remove species that are 5% shorter or longer than the median length
        - Generate MSA using `MAFFT` with default settings
        - Run `TAPER` to mask putative alignment errors
        - Remove sequences with <u>></u>50% gaps
        - Remove sites with <u>></u>50% gaps
    4. Check each BUSCO locus by eye and remove highly-divergent sequences and/or sites that may reflect putative paralogy
    5. Exclude BUSCO loci with no <i>Angophora</i> or <i>Corymbia</i> species, and/or <17 <i>Eucalyptus</i> species

- `crisp101_metadata.tsv` and `busco1187_metadata.tsv` <br>
    This file stores the list of samples with their respective target file for HybPiper analyses

## `quibl/`
- `quibl_metadata.tsv`
    1. Set the focal species as `species_one`
    2. Set the target species as `species two`
    3. Set the sister taxon of `species one` as `species three`
    4. Set `outgroup`

---
<i>Last updated: 13 March 2026 by Jeremias Ivan</i>