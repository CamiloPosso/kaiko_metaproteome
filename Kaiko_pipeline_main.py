import subprocess
import cProfile
import yaml
import os
import re

from pathlib import Path, PureWindowsPath
from s3path import S3Path

from Kaiko_4 import aggregate_fasta
from Kaiko_3 import run_diamond_tally
from Kaiko_2 import combine_denovo_output

## Parsing

kaiko_defaults_path = Path('kaiko_defaults.yaml')
config = yaml.safe_load(kaiko_defaults_path.open())
user_config_path = Path('config.yaml')

if user_config_path.exists():
    config_user = yaml.safe_load(user_config_path.open())

    for section in config_user.keys():
        for param in config_user[section].keys():
            config[section][param] = config_user[section][param]


## handling any backslashes with this. All final paths used here have forward slashes, 
## as they are compatible in Windows, Linux, and Mac.
mgf_dir = Path(PureWindowsPath(config['denovo']['mgf_dir']).as_posix())
ncbi_taxa_folder = Path(PureWindowsPath(config['diamond tally']['ncbi_taxa_folder']).as_posix())
ref_fasta = Path(PureWindowsPath(config['taxa to fasta']['ref_fasta']).as_posix())
prefix = mgf_dir.name

denovout_dir = Path('Kaiko_volume/Kaiko_intermediate/denovo_output/' + prefix)
if not denovout_dir.exists():
    denovout_dir.mkdir()

## Step 1. Run Denovo using subprocess.
kaiko_1_args = ["python", "src/kaiko_main.py", 
                "--mgf_dir", mgf_dir.resolve(), 
                "--train_dir", "model/",
                "--decode_dir", denovout_dir.resolve(),
                "--profile", config['denovo']['profile']]

if config['denovo']['topk']:
    kaiko_1_args = kaiko_1_args + ["--topk"]
if config['denovo']['multi_decode']:
    kaiko_1_args = kaiko_1_args + ["--multi_decode"]
if config['denovo']['beam_search']:
    kaiko_1_args = kaiko_1_args + ["--beam_search", "--beam_size", config['denovo']['beam_size']]     

print("DeNovo: Running the following command:\n")
for i in range(len(kaiko_1_args)):
    kaiko_1_args[i] = str(kaiko_1_args[i])
print(" ".join(kaiko_1_args) + "\n")

subprocess.run(kaiko_1_args, cwd = "Kaiko_denovo")
# subprocess.call(kaiko_1_args, cwd = "Kaiko_denovo")

if (config['denovo'])['profile']:
    profiler = cProfile.Profile()
    profiler.enable()

## Step 2. Combine into fasta
print("\n Combinining denovo output\n")

combine_denovo_output(denovout_dir, prefix)

denovo_combined_fasta = Path("Kaiko_volume/Kaiko_intermediate/" + prefix + "_combined_denovo.fasta")
diamond_search_out = Path("Kaiko_volume/Kaiko_intermediate/" + prefix + "_diamond_search_output.dmd")

## Step 3. Passing to diamond
diamond_args = ["./diamond", "blastp", "-d",
                "../uniref100", "--min-score", "1",
                "-q", denovo_combined_fasta.resolve().as_posix(), "-o",
                diamond_search_out.resolve().as_posix(), "-f", "6", "qseqid", 
                "stitle", "pident", "evalue", "mismatch"]


print("DeNovo: Running the following command:\n")
print(" ".join(diamond_args) + "\n")


# os.chdir("Kaiko_volume/Kaiko_stationary_files/diamond-2.0.15")
os.chdir("Kaiko_volume/Kaiko_stationary_files/diamond-linux")
os.system(" ".join(diamond_args))
os.chdir("../../../")

kaiko_tally = Path("Kaiko_volume/Kaiko_intermediate/" + prefix + "_kaiko_prediction_top_taxa.csv")

# Step 4. Tallying the diamond results
run_diamond_tally(diamond_search_out, 
                  int(config['diamond tally']['ntops']), 
                  ncbi_taxa_folder, 
                  config['diamond tally']['mode'], 
                  kaiko_tally, 
                  float(config['diamond tally']['pident']))


## Step 5. Putting together the final fasta file.


if config['taxa to fasta']['kingdom_list'] != "":
    kingdom_list = config['taxa to fasta']['kingdom_list'].split(', ')
    print(kingdom_list)
else:
    kingdom_list = []

kaiko_final_output = Path("Kaiko_volume/Kaiko_output/" + prefix + "_kaiko_output.fasta")

aggregate_fasta(ref_fasta,
                kaiko_tally,
                kaiko_final_output,
                int(config['taxa to fasta']['ntops']),
                config['taxa to fasta']['taxa_key'],
                kingdom_list)


if (config['denovo'])['profile']:
    profiler.enable()
    profiler.dump_stats('Kaiko_volume/Kaiko_taxa.prof')
