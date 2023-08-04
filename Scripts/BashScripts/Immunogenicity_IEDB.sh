#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem-per-cpu=3700
#SBATCH --time=48:00:00
#SBATCH --job-name=45-55s

### Clear modules:
module purge;
### Activate the conda environment
source activate Immunogenicity_IEDB
cd immunogenicity/

### Execute the python program:

HLAS=("HLA-A0101" "HLA-A0201"  "HLA-A0202"  "HLA-A0203"  "HLA-A0206"  "HLA-A0211"
 "HLA-A0301"  "HLA-A1101"  "HLA-A2301"  "HLA-A2402"  "HLA-A2601"  "HLA-A2902"
 "HLA-A3001"  "HLA-A3002"  "HLA-A3101"  "HLA-A3201"  "HLA-A3301"  "HLA-A6801"
 "HLA-A6802"  "HLA-A6901"  "HLA-B0702"  "HLA-B0801"  "HLA-B1501"  "HLA-B1502"
 "HLA-B1801"  "HLA-B2705"  "HLA-B3501"  "HLA-B3901"  "HLA-B4001"  "HLA-B4002"
 "HLA-B4402"  "HLA-B4403"  "HLA-B4501"  "HLA-B4601"  "HLA-B5101"  "HLA-B5301"
 "HLA-B5401"  "HLA-B5701"  "HLA-B5801")

predict_binding() {
hla=$1
echo "Working on ${hla}."
output_file="../resultsPredictionImmunogenicity/${hla}_immunogenicity.txt"
if [ -f "$output_file" ]; then
        echo "Output file $output_file already exists. Skipping to next HLA."
        return  # Exit the function and move to the next iteration
fi
python predict_immunogenicity.py -a "$hla" ../Sequences/Kmer_CSP_region_273-375_aa.fasta > "$output_file"
}
export -f predict_binding
for hla in "${HLAS[@]}"; do
(predict_binding "$hla") &
done

wait

### Deactivate the conda environment:
conda deactivate;