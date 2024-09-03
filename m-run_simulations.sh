#!/bin/bash -l

# Standard output and error:
#SBATCH -o _outputs_cluster/others/job-%A.out # Standard output, %A = job ID, %a = job array index
#SBATCH -e _outputs_cluster/others/job-%A.err # Standard error, %A = job ID, %a = job array index


# Job Name:
#SBATCH -J seroprevalence

#SBATCH --mail-type=END
#SBATCH --mail-user=domenech@mpiib-berlin.mpg.de

# --- resource specification (which resources for how long) ---
#SBATCH --partition=general 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10000MB # memory in MB required by the job
#SBATCH --time=00:20:00 # run time in h:m:s, up to 24h possible
 
# --- start from a clean state and load necessary environment modules ---
module purge
module load R/4.4

# --- run your executable via srun ---
R --no-save --no-restore <m-run_simulations.R >_outputs_cluster/Rout/R-${SLURM_JOB_ID}.Rout
