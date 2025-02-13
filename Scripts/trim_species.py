import os
import glob
from Bio import SeqIO
import sys

### This is a script by Marta Alvarez Presas from the Institute of Evolutionary Biology (IBE, UPF-CSIC) generated on 
### the 13th of february 2025 with the aim of trimming species from a sequences multifasta file. If you use it, please, cite.


print("Usage: python trim_species.py species_list.txt path_to_alignments output_dir")

# Load the list of outlier species
outliers = set(line.strip() for line in open(sys.argv[1]))

# Define the directory containing the alignments
alignment_dir = sys.argv[2]
output_dir = sys.argv[3]

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Find all alignment files (adjust the extension if necessary)
alignment_files = glob.glob(os.path.join(alignment_dir, "*.fa"))

# Process each alignment file
for aln_file in alignment_files:
    filtered_seqs = []
    with open(aln_file, "r") as infile:
        for record in SeqIO.parse(infile, "fasta"):
            species = record.id.split("_")[0]  # Adjust based on your header format
            if species not in outliers:
                filtered_seqs.append(record)

    # Save filtered alignment
    output_file = os.path.join(output_dir, os.path.basename(aln_file).replace(".fa", "_filtered.fasta"))
    with open(output_file, "w") as outfile:
        SeqIO.write(filtered_seqs, outfile, "fasta")

    print(f"Filtered alignment saved: {output_file}")

print("✅ All alignments processed!")
