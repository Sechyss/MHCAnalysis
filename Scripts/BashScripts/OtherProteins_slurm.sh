#!/bin/bash
#SBATCH --ntasks-per-node=48
#SBATCH --mem-per-cpu=3700
#SBATCH --time=48:00:00
#SBATCH --job-name=Liver1N

### Clear modules:
module purge;
### Activate the conda environment
source activate MHC_prediction;
cd mhc_i/;
### Check the environment path (for information only):

### Execute the python program:

HLAS=("HLA-B*40:06" "HLA-B*44:02" "HLA-A*02:01" "HLA-B*51:01" "HLA-B*56:01" "HLA-A*30:01" "HLA-C*01:02" "HLA-C*04:01" "HLA-A*23:01" "HLA-C*17:01" "HLA-C*16:01" "HLA-C*07:02" "HLA-C*03:03" "HLA-B*42:01" "HLA-C*06:02" "HLA-C*05:01" "HLA-B*41:02" "HLA-C*08:01" "HLA-B*54:01" "HLA-A*02:06" "HLA-A*02:07" "HLA-A*33:03" "HLA-B*40:01" "HLA-B*15:06" "HLA-C*14:03" "HLA-A*03:01" "HLA-B*07:02" "HLA-B*45:01" "HLA-B*53:01" "HLA-B*56:02" "HLA-A*01:01" "HLA-B*58:01" "HLA-B*13:02" "HLA-A*30:02" "HLA-A*26:01" "HLA-B*13:01" "HLA-A*02:11" "HLA-B*08:01" "HLA-B*58:02" "HLA-C*12:03" "HLA-B*46:01" "HLA-A*74:01" "HLA-A*24:02" "HLA-C*03:02" "HLA-A*11:01" "HLA-A*34:01" "HLA-B*40:02" "HLA-B*15:02" "HLA-C*07:01" "HLA-B*18:01" "HLA-C*04:03" "HLA-B*52:01" "HLA-A*68:02" "HLA-B*35:01" "HLA-C*03:04" "HLA-B*35:03")  # Add more HLA types


predict_binding() {
hla=$1
echo "Working on ${hla}."
output_file="../resultsPredictionBinding_Liverstage1/${hla}_length8.txt"
output_file2="../resultsPredictionBinding_Liverstage1/${hla}_length9.txt"
output_file3="../resultsPredictionBinding_Liverstage1/${hla}_length10.txt"
output_file4="../resultsPredictionBinding_Liverstage1/${hla}_length11.txt"
if [ -f "$output_file" ]; then
        echo "Output file $output_file already exists. Skipping to next HLA."
        return  # Exit the function and move to the next iteration
fi
./src/predict_binding.py netmhcpan_el "$hla" 8 ../Sequences_otherproteins/Liverstageantigen1_Nterminal_kmer_filtered_corrected.fasta | awk -F"\t" 'NR==1 || $10 <=1' > "$output_file"
./src/predict_binding.py netmhcpan_el "$hla" 9 ../Sequences_otherproteins/Liverstageantigen1_Nterminal_kmer_filtered_corrected.fasta | awk -F"\t" 'NR==1 || $10 <=1' > "$output_file2"
./src/predict_binding.py netmhcpan_el "$hla" 10 ../Sequences_otherproteins/Liverstageantigen1_Nterminal_kmer_filtered_corrected.fasta | awk -F"\t" 'NR==1 || $10 <=1' > "$output_file3"
./src/predict_binding.py netmhcpan_el "$hla" 11 ../Sequences_otherproteins/Liverstageantigen1_Nterminal_kmer_filtered_corrected.fasta | awk -F"\t" 'NR==1 || $10 <=1' > "$output_file4"
}
export -f predict_binding
for hla in "${HLAS[@]}"; do
(predict_binding "$hla") &
done

wait

### Deactivate the conda environment:
conda deactivate;