#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:10:00
#SBATCH --partition=devel
#SBATCH --job-name=MHC_prediction

### Clear modules:
module purge;

### Activate the conda environment and move to the correct location
source activate MHC_prediction
cd mhc_i/ || exit
### Check the environment path (for information only):
echo "Current Conda environment path is:" $CONDA_PREFIX;
### Execute the python program:
srun python ./src/predict_binding.py netmhcpan_el HLA-A*01:01 11 ../Sequences/Kmer_CSP_region_273-375_aa.fasta

### Deactivate the conda environment:
conda deactivate;

### Handy message for debugging (for information only):
echo "---End Of My Job---";