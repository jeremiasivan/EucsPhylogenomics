################################################
#############    SHORT READS QC    #############
################################################

dir_codes <- "/home/jeremias/EucsPhylogenomics/codes/"
dir_output <- "/data/jeremias/eucs/shortreads/"
thread <- 10

file_metadata <- paste0(dir_codes, "/../files/eucs_metadata.txt") 

dir_shortreads <- "/data/jeremias/eucs/shortreads/"
qc_method <- "adapterremoval"

# run AdapterRemoval
fn_adapters <- "eucs_adapters.txt"

exe_adapterremoval <- "AdapterRemoval"
thread_adapterremoval <- thread
min_read_quality <- 25

# run BBTools
exe_rqcfilter2 <- ""
dir_rqcfilterdata <- ""

################################################

source(paste0(dir_codes, "/functions.R"))

################################################

# open metadata
df_metadata <- data.table::fread(file_metadata)

# extract short reads
ls_shortread_fdname <- list.dirs(dir_shortreads, recursive=F, full.names=F)

# iterate over short reads
for (fdname in ls_shortread_fdname) {
    # extract the name of the species
    read <- df_metadata$tip_label[df_metadata$folder_name==fdname]
    if (length(read) != 1) {
        next
    }

    # check if directory exists
    dir_output_qc <- paste0(dir_output, "/", read, "/")
    if (!dir.exists(dir_output_qc)) {
        dir.create(dir_output_qc, recursive=T)
    }

    # extract FASTQ files
    dir_reads <- paste0(dir_shortreads, "/", fdname, "/Sequence/Raw_Data/")

    df_reads <- file.info(list.files(dir_reads, full.names=T))
    fn_fastq_one <- rownames(df_reads)[which.max(df_reads$mtime)]
    
    # check the forward FASTQ file
    if (length(fn_fastq_one) != 1) {
        next
    }
    
    if (tolower(qc_method) == "adapterremoval") {
        # run AdapterRemoval
        prefix <- paste0(dir_output_qc, read)
        f_qc_adapterremoval(fn_fastq_one, NULL, fn_adapters, prefix, min_read_quality, thread_adapterremoval, exe_adapterremoval)
    } else if (tolower(qc_method) == "bbtools") {
        # run BBTools
        f_qc_bbtools(fn_fastq_one, dir_output_qc, dir_rqcfilterdata, exe_rqcfilter2)

        # rename file (tbc)
    }
}

################################################