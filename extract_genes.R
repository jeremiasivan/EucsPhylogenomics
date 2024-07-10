# Extract Genes from Assemblies
dir_output <- ""
thread <- 10

fn_refseq <- "eucs_refseq.txt"
fn_refseq_chr <- "eucs_refseq_chr.txt"

prefix_dir_fasta <- "/data/jeremias/eucs/roadies/assemblies"
prefix_dir_gtf <- "/data/jeremias/eucs/genes/"

exe_gffread <- "gffread"

################################################

# load libraries
# require(data.table)
library(doSNOW)

# open data.table
df_refseq <- data.table::fread(fn_refseq)

# create doSNOW cluster
nwcl <- makeCluster(thread)
doSNOW::registerDoSNOW(nwcl)

# iterate over reference sequences
foreach (i = 1:nrow(df_refseq)) %dopar% {
    # update the directories
    ref <- df_refseq$id[i]
    fn_fasta <- paste0(prefix_dir_fasta, df_refseq$dir_fasta[i])
    fn_gtf <- paste0(prefix_dir_gtf, df_refseq$dir_gtf[i])

    # output file
    fn_outfile <- paste0(dir_output, "/", ref, ".faa")

    # run Gffread
    cmd_gffread <- paste(exe_gffread, "-g", fn_fasta, "-x", fn_outfile, fn_gtf)
    system(cmd_gffread)
}

stopCluster(nwcl)