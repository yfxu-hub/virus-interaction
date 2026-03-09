#!/bin/bash -l

# Standard output and error:
#SBATCH -o results/temp/job-compile.out # Standard output
#SBATCH -e results/temp/job-compile.err # Standard error

# Job Name:
#SBATCH -J tm

#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=kramer@mpiib-berlin.mpg.de

# --- resource specification (which resources for how long) ---
#SBATCH --partition=general 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=500 # memory in MB required by the job
#SBATCH --time=0:10:00 # run time in h:m:s, up to 24h possible
 
# --- start from a clean state and load necessary environment modules ---
module purge
module load R/4.2

# --- run your executable via srun ---
R --no-save --no-restore <src/process_results/02_compile_results_traj_matching_r1.R >results/Rout/R-tm-compile.Rout
