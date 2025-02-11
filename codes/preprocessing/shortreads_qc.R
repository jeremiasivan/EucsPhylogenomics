################################################
#############    SHORT READS QC    #############
################################################

dir_codes <- "/home/jeremias/EucsPhylogenomics/codes/"
dir_output <- "/data/jeremias/eucs/shortreads/"
thread <- 10
n_parallel_run <- 1

file_metadata <- paste0(dir_codes, "/../files/eucs_metadata.txt") 

dir_shortreads <- "/data/jeremias/eucs/shortreads/"
qc_method <- "adapterremoval"

# run AdapterRemoval
fn_adapters <- "eucs_adapters.txt"

exe_adapterremoval <- "AdapterRemoval"
thread_adapterremoval <- thread / n_parallel_run
min_read_quality <- 25

# run BBTools
exe_rqcfilter2 <- ""
exe_reformat <- ""
dir_rqcfilterdata <- ""

################################################

source(paste0(dir_codes, "/functions.R"))

library(doSNOW)

################################################

# set-up log file
fn_log <- paste0(dir_output, "/log.txt")
f_write_log(fn_log, "------------------- START")

# open metadata
df_metadata <- data.table::fread(file_metadata)

# extract short reads
ls_shortread_fdname <- list.dirs(dir_shortreads, recursive=F, full.names=F)

# create doSNOW cluster
nwcl <- makeCluster(n_parallel_run)
doSNOW::registerDoSNOW(nwcl)

# iterate over short reads
foreach (fdname = ls_shortread_fdname) %dopar% {
    # extract the name of the species
    read <- df_metadata$tip_label[df_metadata$folder_name==fdname]
    if (length(read) != 1) {
        return(NULL)
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
    fn_fastq_one_name <- unlist(strsplit(fn_fastq_one[1], split="/"))
    fn_fastq_one_name <- gsub("*.fastq.gz", "", fn_fastq_one_name[length(fn_fastq_one_name)])

    # check the forward FASTQ file
    if (length(fn_fastq_one) != 1) {
        system(paste("rm -r", dir_output_qc))
        return(NULL)
    }

    # add log file (tbc)
    f_write_log(fn_log, paste0("- ", read, ": ", fn_fastq_one))
    
    if (tolower(qc_method) == "adapterremoval") {
        # run AdapterRemoval
        prefix <- paste0(dir_output_qc, read)
        f_qc_adapterremoval(fn_fastq_one, NULL, fn_adapters, prefix, min_read_quality, thread_adapterremoval, exe_adapterremoval)
    } else if (tolower(qc_method) == "bbtools") {
        # output file
        ls_output <- strsplit(fn_fastq_one_name, split="/")

        # check if previous run exists
        fn_output_log <- paste0(dir_output_qc, "status.log")
        is_complete <- suppressWarnings(system(paste("grep 'RQCFilter complete'", fn_output_log), intern=T))
        if (length(is_complete) == 0) {
            system(paste("rm -r", dir_output_qc))
            dir.create(dir_output_qc, recursive=T)
        }

        # run BBTools
        fn_output <- paste0(dir_output_qc, ls_output[length(ls_output)], ".anqdt.fastq.gz")
        if (!file.exists(fn_output)) {
            f_qc_bbtools(fn_fastq_one, dir_output_qc, dir_rqcfilterdata, exe_rqcfilter2)
        }

        # rename file
        fn_output_one <- paste0(dir_output_qc, read, "_R1.fastq.gz")
        fn_output_two <- paste0(dir_output_qc, read, "_R2.fastq.gz")
        if (!all(file.exists(fn_output_one, fn_output_two))) {
            system(paste0(exe_reformat, " in=", fn_output, " out1=", fn_output_one, " out2=", fn_output_two))
        }   
    }
}

stopCluster(nwcl)
f_write_log(fn_log, c("--------------------- END", ""))

################################################