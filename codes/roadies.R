# Run ROADIES
thread <- 10

dir_roadies <- "/home/jeremias/ROADIES/"

fn_config <- "/data/jeremias/eucs/assemblies/"

################################################

# function: run ROADIES
f_roadies <- function(fn_config, thread, exe_roadies) {
    cmd_roadies <- paste("python", exe_roadies,
                         "--config", fn_config,
                         "--cores", thread, "--converge")
    system(cmd_roadies)
}

################################################

# move to ROADIES directory
setwd(dir_roadies)

# source the environment
system("source roadies_env.sh")

# run ROADIES
f_roadies(fn_config, thread, "run_roadies.py")

################################################