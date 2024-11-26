################################################
#############        MAPPING       #############
################################################

dir_codes <- "/home/jeremias/EucsPhylogenomics/codes/"
dir_output <- "/data/jeremias/eucs/shortreads/"
thread <- 10

dir_refseq <- ""
dir_shortreads <- "/data/jeremias/eucs/shortreads/"

fn_metadata <- ""

# mapping
exe_mash <- ""

exe_bwamem2: ""
exe_samtools: ""
exe_qualimap: ""
exe_bcftools: ""

################################################

source(paste0(dir_codes, "/functions.R"))

# create output directory
dir_output_mash <- paste0(dir_output, "/mash/")
if (!dir.exists(dir_output_mash)) {
    dir.create(dir_output_mash, recursive=T)
}

################################################

# open metadata
df_metadata <- data.table::fread(fn_metadata)

# lists of refseq and short reads
ls_dir_refseq <- list.files(dir_refseq, pattern="*.fa$", recursive=F, full.names=T)
ls_refseq <- sapply(ls_dir_refseq, function(x) {
    file_fasta <- unlist(strsplit(x, split="/"))
    gsub("*.fa", "", file_fasta[length(file_fasta)])
})

ls_shortreads <- list.dirs(dir_shortreads, recursive=F, full.names=F)

# output files
file_msh <- paste0(dir_output_mash, "/refsketch.msh")
file_mashsum <- paste0(dir_output_mash, "/refsketch.mashsum")

################################################
#############         MASH         #############
################################################

# convert list to string
str_refseq_fasta <- paste(ls_dir_refseq, collapse=" ")

# run Mash sketch on all reference sequences
prefix <- paste0(dir_output_mash, "/refsketch")
f_mash_sketch(exe_mash, str_refseq_fasta, thread, prefix)

# create doSNOW cluster
nwcl <- makeCluster(thread)
doSNOW::registerDoSNOW(nwcl)

# iterate over shortreads
foreach (i = ls_shortreads) %dopar% {
    fn_fastq <- paste0(dir_shortreads, "/", i, ".fastq.gz")
    fn_output <- paste0(dir_output_mash, i, ".mashtab")

    f_mash_screen(exe_mash, file_msh, fn_fastq, fn_output)
}

stopCluster(nwcl)

# summarize .mashtab for all short reads
df_mashsum <- data.table::data.table(reads=character(), reads_dir=character(), reads_species=character(),
                                     refseq=character(), refseq_dir=character(), refseq_species=character(),
                                     score=numeric())

# iterate over shortreads
for (i in ls_shortreads) {
    # input files
    reads_dir <- paste0(dir_shortreads, "/", i, ".fastq.gz")
    reads_species <- df_metadata$species_name[df_metadata$tip_label==i]

    file_mashtab <- paste0(dir_output_mash, i, ".mashtab")

    # extract the best Mash sketch
    cmd_sort_mashtab <- paste("sort -gr", file_mashtab, "| head -1")
    best_sketch <- system(cmd_sort_mashtab, intern=T)
    ls_best_sketch <- unlist(strsplit(best_sketch, split="\t"))

    # extract the accession number
    refseq_dir <- ls_best_sketch[5]
    ls_refseq_dir <- unlist(strsplit(refseq_dir, split="/"))
    refseq_id <- gsub("*.fa", "", ls_refseq_dir[length(ls_refseq_dir)])

    # extract the score
    score <- round(as.numeric(ls_best_sketch[1]), 4)
    
    # save the information into data.table
    df_mashsum <- rbind(df_mashsum, list(reads_id=i, reads_dir=reads_dir, reads_species=reads_species,
                                         refseq_id=refseq_id, refseq_dir=refseq_dir, refseq_species=refseq_id,
                                         score=score))
}

data.table::fwrite(df_mashsum, file=file_mashsum, quote=F, sep="\t")

################################################