################################################
#############    SHORT READS QC    #############
################################################

dir_codes <- "/home/jeremias/EucsPhylogenomics/codes/"
dir_output <- "/data/jeremias/eucs/shortreads/"
thread <- 10

dir_shortreads <- "/data/jeremias/eucs/shortreads/"
qc_method <- "adapterremoval"

# run AdapterRemoval
fn_adapters <- "eucs_adapters.txt"

exe_adapterremoval <- "AdapterRemoval"
thread_adapterremoval <- thread
min_read_quality <- 25

################################################

source(paste0(dir_codes, "/functions.R"))

################################################

# extract short reads
ls_shortread_qc <- list.dirs(dir_shortreads, recursive=F, full.names=F)

# create doSNOW cluster
nwcl <- makeCluster(floor(thread/thread_adapterremoval))
doSNOW::registerDoSNOW(nwcl)

# iterate over short reads
foreach (read = ls_shortread_qc) %dopar% {
    # check if directory exists
    dir_output_qc <- paste0(dir_output, "/", read, "/")
    if (!dir.exists(dir_output_qc)) {
        dir.create(dir_output_qc, recursive=T)
    }

    prefix <- paste0(dir_output_qc, read)

    # extract FASTQ files
    dir_reads <- paste0(dir_shortreads, "/", read, "/")
    fn_fastq_one <- paste0(dir_reads, read, "_1.fastq")
    fn_fastq_two <- paste0(dir_reads, read, "_2.fastq")
    
    # check the forward FASTQ file
    if (!file.exists(fn_fastq_one)) {
        return(NULL)
    }

    # check the reverse FASTQ file
    if (!file.exists(fn_fastq_two)) {
        fn_fastq_two <- NULL
    }
    
    if (tolower(qc_method) == "adapterremoval") {
        # run AdapterRemoval
        f_qc_short_reads(fn_fastq_one, fn_fastq_two, fn_adapters, prefix, min_read_quality, thread_adapterremoval, exe_adapterremoval)
    }
}

stopCluster(nwcl)

################################################