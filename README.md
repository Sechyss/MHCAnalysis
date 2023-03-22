# MHC analysis
This is a set of scripts to make them publicly available for everyone to repeat this analysis

1. Scripts in the Python console are drafts attempts to create the analysis
2. Rest of the scripts are executed in the command line either as they are or 
by running the flags specified in the help command

## Order of the scripts to be executed

python /PATH/MHC_prep_list.py [-f table.xlsx] [-l list.txt] [-p POPULATION]

cd /PATH/mhc_i

python ./src/predict_binding.py [method] [outputfrompreviouspythoncode] [sequence.faa] > [outputfile.txt]

python /PATH/VCF_MHC_Analysis.py [-v file.vcf.gz] [-m file.txt] [-s FASTA] [-o file.pdf]

python /PATH/CASP_peptide_SNP_kmer_analysis.py [-v file.csv] [-g file.fasta] [-s START] [-e END] [-o file.fasta]

python ./src/predict_binding.py [method] [outputfrompreviouspythoncode] [kmerseqs.faa] > [kmeroutputfile.txt]

python /PATH/CSP_SummaryTable.py [-t table.txt] [-d dict.pickle] [-p NEWDICT] [-o file.xlsx]