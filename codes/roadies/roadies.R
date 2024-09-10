################################################
#############        ROADIES       #############
################################################

dir_codes <- "/home/jeremias/EucsPhylogenomics/codes/"
thread <- 10

# run ROADIES
dir_roadies <- "/home/jeremias/ROADIES/"
fn_config <- "/data/jeremias/eucs/assemblies/"

################################################

source(paste0(dir_codes, "/functions.R"))

################################################

# move to ROADIES directory
setwd(dir_roadies)

# source the environment
# system("source roadies_env.sh")

# run ROADIES
f_roadies(fn_config, thread, "run_roadies.py")

################################################