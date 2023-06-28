import os
import sys
import subprocess
from multiprocessing import Pool

### This is a script generated  by Marta Alvarez Presas at the Institute of Evolutionary Biology on the 28th of june 2023. The purpose of the script is bunch processing
### a series of proteomes to search the BUSCO values in the Eukaryota database. Note that you need to change the path to the database

def run_busco(genome_file, output_dir):
    print(f"Running BUSCO for {genome_file}...")
    cmd = f"busco -i {genome_file} -o {output_dir} -l path/to/the/database/eukaryota_odb10 -m protein"
    subprocess.call(cmd, shell=True)
    print(f"BUSCO run complete for {genome_file}.")

# Folder containing the FASTA files
folder_path = sys.argv[1]

# Number of parallel processes
num_processes = 4

# Get the list of FASTA files in the folder
fasta_files = [os.path.join(folder_path, f) for f in os.listdir(folder_path) if f.endswith(".fa")]

# Create a pool of worker processes
pool = Pool(processes=num_processes)

# Run BUSCO searches in parallel
results = []
for fasta_file in fasta_files:
    genome_file = fasta_file
    output_dir = f"busco_output_{os.path.splitext(os.path.basename(fasta_file))[0]}"
    result = pool.apply_async(run_busco, (genome_file, output_dir))
    results.append(result)

# Wait for all processes to complete
pool.close()
pool.join()

# Get the results (optional)
for result in results:
    result.get()
