import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#%matplotlib inline

import seaborn as sns


import os, sys, shutil, importlib, glob
from tqdm.notebook import tqdm

from celloracle import motif_analysis as ma
import celloracle as co
co.__version__

plt.rcParams['figure.figsize'] = [6, 4.5]
plt.rcParams["savefig.dpi"] = 300
arg1 = sys.argv[1]
print(arg1)
arg2 = sys.argv[2]
print(arg2)
# Load scATAC-seq peak list.
peaks = pd.read_csv(arg1)
peaks = peaks.x.values
print(peaks)

# Load Cicero coaccessibility scores.
cicero_connections = pd.read_csv(arg2)
cicero_connections.head()

##!! Please make sure to specify the correct reference genome here
tss_annotated = ma.get_tss_info(peak_str_list=peaks, ref_genome="hg38")

# Check results
print(tss_annotated.tail())

integrated = ma.integrate_tss_peak_with_cicero(tss_peak=tss_annotated,
                                               cicero_connections=cicero_connections)
print(integrated.shape)
print(integrated.head())

peak = integrated[integrated.coaccess >= 0.8]
peaks = peak[["peak_id", "gene_short_name"]].reset_index(drop=True)

print(peaks.shape)
print(peaks.head())

#peak.to_csv("scATAC/E7.5_rep1_1000_cells_peak_file_TSS_annot.csv")

# PLEASE make sure reference genome is correct.
ref_genome = "hg38"

genome_installation = ma.is_genome_installed(ref_genome=ref_genome,
                                             genomes_dir=None)
print(ref_genome, "installation: ", genome_installation)




peaks = ma.check_peak_format(peaks, ref_genome, genomes_dir=None)

# Instantiate TFinfo object
tfi = ma.TFinfo(peak_data_frame=peaks,
                ref_genome=ref_genome,
                genomes_dir=None)

# Scan motifs. !!CAUTION!! This step may take several hours if you have many peaks!
tfi.scan(fpr=0.02,
         motifs=None,  # If you enter None, default motifs will be loaded.
         verbose=True)

# Save tfinfo object
tfi.to_hdf5(file_path=f"{arg1}.test1.celloracle.tfinfo")

# Check motif scan results
print(tfi.scanned_df.head())

# Reset filtering
tfi.reset_filtering()

# Do filtering
tfi.filter_motifs_by_score(threshold=10)

# Format post-filtering results.
tfi.make_TFinfo_dataframe_and_dictionary(verbose=True)

df = tfi.to_dataframe()
print(df.head())

# Save result as a dataframe
#df = tfi.to_dataframe()
df.to_parquet(f"{arg1}_base_GRN.parquet")
df.to_csv(f"{arg1}_base_GRN.csv")
