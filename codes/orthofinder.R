# Run Orthofinder
dir_output <- "/data/jeremias/eucs/orthofinder/"
thread <- 10

exe_orthofinder <- ""
dir_fasta <- ""

################################################

# function: run ROADIES
f_orthofinder <- function(dir_fasta, dir_output, thread, exe_orthofinder) {
    cmd_orthofinder <- paste(exe_orthofinder,
                         "-f", dir_fasta,
                         "-o", dir_output,
                         "-t", thread)
    system(cmd_orthofinder)
}

################################################

# create output directories
dir_output_db <- paste0(dir_output, "/", tolower(reftaxonomy), "_db/")
if (!dir.exists(dir_output_db)) {
    dir.create(dir_output_db, recursive=T)
}

# run ROADIES
f_orthofinder(fn_config, thread, "run_roadies.py")

################################################