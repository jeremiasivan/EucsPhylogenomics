# Extract Genes from Assemblies
dir_output <- "/data/jeremias/eucs/genes/proteome/"
thread <- 10

fn_refseq <- "eucs_refseq.txt"
fn_refseq_chr <- "eucs_refseq_chr.txt"

prefix_dir_fasta <- "/data/jeremias/eucs/roadies/assemblies/"
prefix_dir_gtf <- "/data/jeremias/eucs/genes/"

exe_gffread <- "gffread"

################################################

# load libraries
# require(data.table)
# require(tidyr)
library(doSNOW)

# create output directory
if (!dir.exists(dir_output)) {
    dir.create(dir_output, recursive=T)
}

# open data.table
df_refseq <- data.table::fread(fn_refseq)

# create doSNOW cluster
nwcl <- makeCluster(thread)
doSNOW::registerDoSNOW(nwcl)

# iterate over reference sequences
foreach (i = 1:nrow(df_refseq)) %dopar% {
    # update the directories
    ref_sp <- df_refseq$species[i]
    ref_id <- df_refseq$id[i]
    fn_fasta <- paste0(prefix_dir_fasta, df_refseq$dir_fasta[i])
    fn_gtf <- paste0(prefix_dir_gtf, df_refseq$dir_gtf[i])

    # output file
    fn_outfile <- paste0(dir_output, "/", ref_id, ".faa")

    # read GTF file
    df_gtf <- data.table::fread(fn_gtf)
    data.table::setnames(df_gtf, c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"))
    df_refseq_chr <- data.table::fread(fn_refseq_chr)

    # convert GTF file to match naming convention
    df_refseq_chr_long <- tidyr::pivot_longer(df_refseq_chr, cols=-chr, names_to="species", values_to="assembly")
    df_refseq_chr_subset <- df_refseq_chr_long[df_refseq_chr_long$species==ref_sp,]

    # iterate over chromosomes
    for (j in 1:nrow(df_refseq_chr_subset)) {
        # update the name of the chromosome following the assembly
        df_gtf$seqname <- gsub(df_refseq_chr_subset$chr[j], df_refseq_chr_subset$assembly[j], df_gtf$seqname)
    }

    # save the GTF file
    data.table::fwrite(df_gtf, file=fn_gtf, sep="\t", quote=F, col.names=F)

    # run Gffread
    cmd_gffread <- paste(exe_gffread, "-g", fn_fasta, "-V -y", fn_outfile, fn_gtf)
    system(cmd_gffread)
}

stopCluster(nwcl)