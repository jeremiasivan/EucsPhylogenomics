This README file lists out the steps to retrieve the metadata files.
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

- `combined_metadata_filtered.tsv`
    1. Separate the forward and reverse reads from the filtered FASTQ files using `reformat.sh` from `BBTools`
    2. Calculate the file sizes of the FASTQ files in R using `file.size()`
    3. Exclude 31 samples with either the forward or reverse reads <1GB
    4. Save the data.frame using `data.table::fwrite()` 

---
<i>Last updated: 16 April 2025 by Jeremias Ivan</i>