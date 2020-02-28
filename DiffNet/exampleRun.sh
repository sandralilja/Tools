#!/bin/bash
#
#SBATCH -J SAR_A_3D_1
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -t 72:00:00
#SBATCH --exclusive
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sandra.lilja@liu.se


# Run a single task in the foreground.
module load Python/3.7.0-anaconda-5.3.0-extras-nsc1
source activate abi
export OMP_NUM_THREADS=32
python Differential_Network_analysis.py --seed_from 1 --seed_to 1000 --timepoint 3D

#
# Script ends here





