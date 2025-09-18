#conda init bash
source /gpfs/Home/kmk7420/.conda/pkgs/conda-24.7.1-py39hf3d152e_0/lib/python3.9/site-packages/conda/shell/etc/profile.d/conda.sh
#source /gpfs/Home/kmk7420/.conda/pkgs/conda-22.11.1-py39hf3d152e_1/etc/profile.d/conda.sh
conda activate celloracle_env
source /etc/profile.d/modules.sh
module load R/4.3.2
module load bedtools
input_atac=$1
output_atac=$2
input_rna=$3
output_rna=$4
clean_input=$(echo "$input_atac" | sed 's/[\/.]/_/g')
all_peak="${output_atac}_all_peaks.csv"
cicero_connection="${output_atac}_cicero_connections.csv"
int_RNA_input="${output_rna}_processed.h5ad"
final_out=$5
int_ATAC_input="${all_peak}_base_GRN.parquet"
# Create a timestamped file with the input name included
#date > "start_time_${clean_input}_32c_1000mem.txt"
#echo $int_RNA_input
#echo $final_out
#echo $int_ATAC_input

Rscript ../SCRIPT.CELLORACLE/step1_preprocess.R ${input_atac} ${output_atac}
python ../SCRIPT.CELLORACLE/step3_TSS_annot.py ${all_peak} ${cicero_connection}
python ../SCRIPT.CELLORACLE/steps_all_RNA.py ${input_rna} ${output_rna}
python ../SCRIPT.CELLORACLE/step7_final_integration.py ${int_RNA_input} ${int_ATAC_input} ${final_out}
#subsampled_ATAC_1000_all_peaks.csv_base_GRN.parquet
