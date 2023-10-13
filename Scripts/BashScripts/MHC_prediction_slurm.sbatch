#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem-per-cpu=3700
#SBATCH --time=48:00:00
#SBATCH --job-name=all

### Clear modules:
module purge;
### Activate the conda environment
source activate MHC_prediction;
# shellcheck disable=SC2164
cd mhc_i/;
### Check the environment path (for information only):

### Execute the python program:

HLAS=("list of HLAs to analyze")

predict_binding() {
hla=$1
echo "Working on ${hla}."
output_file="../resultsPredictionBinding/${hla}_length8.txt"
if [ -f "$output_file" ]; then
        echo "Output file $output_file already exists. Skipping to next HLA."
        return  # Exit the function and move to the next iteration
fi
./src/predict_binding.py netmhcpan_el "$hla" 8 ../Sequences/Kmer_CSP_region_273-375_aa.fasta | awk -F"\t" 'NR==1 || $10 <=1' > "$output_file"
}
export -f predict_binding
for hla in "${HLAS[@]}"; do
(predict_binding "$hla") &
done

wait

### Deactivate the conda environment:
conda deactivate;