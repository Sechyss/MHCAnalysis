rule SetDir:
    shell:
        'cd ~/Home/Downloads/mhc_i'

rule ListMHC:
    output:
        "../data/ListMHC_{mhc}.txt"
    shell:
        "./src/predict_binding.py {wildcards.mhc} mhc | grep -e 'HLA'| uniq > {output}"

rule MHC_prep_list:
    output: "../data/ListMHC_to_run.txt"
    shell:
        "python/Users/u2176312/PycharmProjects/PythonScripts/MHC_Scripts/MHC_prep_list.py"

rule MHC_prediction:
    input:
        fastafile="../data/{sequence}.faa",
        mhclength= rules.MHC_prep_list.output

    output:
        "../data/MHC_prediction_{sequence}.txt"
    shell:
        "./src/predict_binding.py smm {mhcs} {mhclenghts} {input.fastafile} > {output}"

rule VCF_analysis:
    input:
        MHCfile="../data/MHC_prediction_{sequence}.txt",
        NucleotideSeq="../data/{sequence}.fna",
        VCF_file="../data/{sequence}.vcf.gz"
    output:
        PDF_figure="../data/"
    shell:
        "python3 ../data/VCF_MHC_Analysis.py {input.VCF_file} {input.MHCfile} {input.NucleotideSeq} {output.PDF_figure}"