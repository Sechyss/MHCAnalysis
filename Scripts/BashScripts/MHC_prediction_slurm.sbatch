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

HLAS=("HLA-C*06:02" "HLA-B*56:05" "HLA-C*15:03" "HLA-B*40:25" "HLA-A*30:02" "HLA-B*07:20" "HLA-A*02:24" "HLA-B*40:11" "HLA-A*02:02" "HLA-B*58:04" "HLA-A*01:01" "HLA-C*16:02" "HLA-B*35:01" "HLA-C*07:27" "HLA-A*24:22" "HLA-C*04:07" "HLA-A*02:22" "HLA-B*18:10" "HLA-A*26:07" "HLA-A*02:25" "HLA-B*55:08" "HLA-B*40:26" "HLA-B*07:12" "HLA-A*01:02" "HLA-B*07:11" "HLA-B*58:01" "HLA-B*40:15" "HLA-A*68:08" "HLA-B*58:06" "HLA-B*15:24" "HLA-B*15:39" "HLA-C*07:13" "HLA-A*33:04" "HLA-B*15:25" "HLA-B*35:07" "HLA-B*40:05" "HLA-C*03:14" "HLA-A*23:02" "HLA-B*47:03" "HLA-A*66:02" "HLA-B*35:24" "HLA-B*44:12" "HLA-B*15:14" "HLA-B*15:37" "HLA-B*27:10" "HLA-A*24:29" "HLA-C*05:04" "HLA-B*44:10" "HLA-C*15:02" "HLA-C*06:07" "HLA-A*34:02" "HLA-B*55:05" "HLA-B*07:19" "HLA-C*07:15" "HLA-B*56:03" "HLA-B*35:18" "HLA-A*30:04" "HLA-B*07:03" "HLA-C*08:09" "HLA-B*27:15" "HLA-C*15:10" "HLA-B*15:51" "HLA-B*27:06" "HLA-B*53:01" "HLA-A*02:06" "HLA-B*44:13" "HLA-A*30:03" "HLA-B*78:03" "HLA-B*27:01" "HLA-B*46:01" "HLA-B*48:06" "HLA-C*06:05" "HLA-B*40:09" "HLA-B*51:17" "HLA-C*08:06" "HLA-C*01:02" "HLA-B*38:02" "HLA-B*15:49" "HLA-A*24:17" "HLA-B*51:20" "HLA-B*55:10" "HLA-B*15:32" "HLA-B*53:02" "HLA-B*35:15" "HLA-A*11:04" "HLA-B*51:19" "HLA-B*35:23" "HLA-C*17:03" "HLA-B*35:12" "HLA-A*80:01" "HLA-B*37:02" "HLA-B*56:01" "HLA-B*15:40" "HLA-A*68:09" "HLA-B*57:06" "HLA-B*07:17" "HLA-B*78:02" "HLA-A*26:08" "HLA-B*07:16" "HLA-B*15:61")

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