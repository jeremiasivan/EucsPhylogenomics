################################################
#######    TARGET GENES VISUALISATION    #######
################################################

dir_codes <- "/home/jeremias/EucsPhylogenomics/codes/"
dir_output <- "/data/jeremias/eucs/shortreads/"
thread <- 10

file_targetgenes <- "target_genes.fa"

################################################

# require(Biostrings)
# require(data.table)
# require(pwalign)

library(doSNOW)
library(ggplot2)

source(paste0(dir_codes, "/functions.R"))

################################################

# open alignment
target_genes <- Biostrings::readDNAStringSet(filepath=file_targetgenes, format="fasta")
n_headers <- length(target_genes)

# create data.table
df_pairs <- expand.grid(1:n_headers, 1:n_headers)
df_pairs <- df_pairs[df_pairs$Var1 <= df_pairs$Var2, ]

# create doSNOW cluster
nwcl <- makeCluster(thread)
doSNOW::registerDoSNOW(nwcl)

# iterate over sequences
df_output <- foreach(i = 1:nrow(df_pairs), .combine='rbind') %dopar% {
    # set indices
    id_seq1 <- df_pairs$Var1[i]
    id_seq2 <- df_pairs$Var2[i]

    # pairwise alignment
    global_align <- pwalign::pairwiseAlignment(target_genes[id_seq1], target_genes[id_seq2], type = "global")

    return(data.table::data.table(seq1=names(target_genes[id_seq1]), seq2=names(target_genes[id_seq2]), score=global_align@score))
}

stopCluster(nwcl)

# save data.table
fn_output <- paste0(dir_output, "/pairwise_scores.txt")
data.table::fwrite(df_output, file=fn_output, sep="\t", quote=F)

# visualisation (tbc)
