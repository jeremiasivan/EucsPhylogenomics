This README file lists out the steps to retrieve the metadata files.
- `currency_creek/jgi_metadata.csv`
    1. Go to <a href="https://genome.jgi.doe.gov/portal/">JGI Genome Portal</a>
    2. Search for <i>"Currency Creek"</i>
    3. Click on <b><i>Reports</i></b> -> <b><i>Download Project Overview Report</i></b>

- `currency_creek/justin_metadata.csv`
    1. Go to this <a href="https://raw.githubusercontent.com/borevitzlab/cca-eucs/master/metadata/CCATreesCombinedMetadata.csv">link</a>
    2. Download the file (e.g., using `wget`) 

- `currency_creek/combined_metadata.tsv`
    1. Open `currency_creek/jgi_metadata.csv` in R
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
    5. Open `currency_creek/justin_metadata.csv`
    6. Update column names with lowercases and underscores
    7. Merge the two data.frames using:
    ```
    df_combined <- df_jgi %>% left_join(df_justin, by="field_id", keep = FALSE)
    ```
    8. Remove duplicated entry for `field_id == "CCA0154"` as the project name is <i>E. corrugata</i> but identified as <i>E. crucis subsp. crucis</i>
    9. Create the `tip_label` column based on the `species_name` column with iterator
    10. Save the data.frame using `data.table::fwrite()`

- `all/seqkit_stats.tsv`
    1. Separate the forward and reverse reads from the filtered FASTQ files using `reformat.sh` from `BBTools`
    2. Store the forward and reverse FASTQ files for each sample in one folder
    3. Run `seqkit stats -a *.fastq.gz > seqkit_stats.tsv`

    > **Note**: `E_effusa_exsul_1` was excluded as it was unavailable at the time of download (01 October 2024)

- `all/seqkit_cov.tsv`
    1. Open `all/seqkit_stats.tsv`
    2. For each sample, add up the `sum_len` column between forward and reverse reads
    3. Create the `cov` column by dividing the `sum_len` column with $5\cdot10^8$ (i.e., the average genome size of Eucalypts) 

- `all/tree_taxonomy_label.tsv`
    This file extracts the taxonomic groupings of 837 Eucalypts from `currency_creek/combined_metadata.tsv` and Ferguson et al. (<a href="https://doi.org/10.1101/gr.277999.123">2024</a>)

- `quibl/quibl_metadata.tsv`
    1. Set the focal species as `species_one`
    2. Set the target species as `species two`
    3. Set the sister taxon of `species one` as `species three`
    4. Set `outgroup`

---
<i>Last updated: 12 March 2026 by Jeremias Ivan</i>