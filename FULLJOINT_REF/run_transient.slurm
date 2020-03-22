#!/bin/bash
#SBATCH --job-name=FJREFTRANS
#SBATCH --ntasks=1
#SBATCH --mem=5G
#SBATCH --time=5:00:00
#SBATCH --mail-user=nb25@rice.edu
#SBATCH --mail-type=ALL
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=2
# #SBATCH --exclusive
#SBATCH --partition=commons
#
echo "I ran on:"
cd $SLURM_SUBMIT_DIR
echo $SLURM_NODELIST
#
module load MATLAB/2019a
#
matlab -nodisplay -r "transient_fulljoint_nots; quit"

