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

HLAS=("HLA-B*13:02" "HLA-C*14:03" "HLA-C*17:03" "HLA-A*66:01" "HLA-B*42:01" "HLA-A*74:03" "HLA-B*15:06" "HLA-B*15:03" "HLA-B*07:06" "HLA-B*35:08" "HLA-A*69:01" "HLA-B*44:09" "HLA-B*50:04" "HLA-B*35:03" "HLA-A*02:24" "HLA-A*30:02" "HLA-B*15:16" "HLA-A*29:10" "HLA-B*57:03" "HLA-A*23:02" "HLA-A*36:01" "HLA-C*12:03" "HLA-B*41:02" "HLA-C*04:07" "HLA-C*03:02" "HLA-B*27:02" "HLA-C*04:03" "HLA-C*14:02" "HLA-C*01:02" "HLA-B*37:01" "HLA-C*06:03" "HLA-A*26:01" "HLA-B*15:10" "HLA-A*43:01" "HLA-A*02:01" "HLA-C*07:03" "HLA-A*02:02" "HLA-A*29:01" "HLA-C*07:18" "HLA-B*47:02" "HLA-B*39:01" "HLA-C*06:02" "HLA-A*11:09" "HLA-C*07:06" "HLA-A*32:03" "HLA-C*12:02" "HLA-C*07:10" "HLA-B*07:05" "HLA-B*27:03" "HLA-B*39:06" "HLA-B*52:01" "HLA-A*68:02" "HLA-B*51:01" "HLA-A*74:01" "HLA-B*58:02" "HLA-C*15:04" "HLA-B*38:01" "HLA-B*55:01" "HLA-A*34:01" "HLA-B*18:01" "HLA-C*17:01" "HLA-C*15:06" "HLA-B*14:01" "HLA-B*18:05" "HLA-B*35:02" "HLA-B*08:01" "HLA-C*15:05" "HLA-B*35:01" "HLA-C*16:02" "HLA-A*34:02" "HLA-C*07:02" "HLA-B*50:01" "HLA-B*57:01" "HLA-B*58:01" "HLA-B*15:24" "HLA-B*81:01" "HLA-A*11:01" "HLA-B*14:02" "HLA-A*32:01" "HLA-A*33:01" "HLA-B*56:01" "HLA-A*03:02" "HLA-C*16:04" "HLA-B*50:02" "HLA-B*51:07" "HLA-A*02:11" "HLA-C*07:04" "HLA-A*02:05" "HLA-B*15:01" "HLA-C*05:01" "HLA-C*03:04" "HLA-B*44:02" "HLA-C*18:01" "HLA-C*03:03" "HLA-C*15:02" "HLA-C*12:04" "HLA-B*15:21" "HLA-B*47:01" "HLA-B*15:17" "HLA-B*45:01")

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