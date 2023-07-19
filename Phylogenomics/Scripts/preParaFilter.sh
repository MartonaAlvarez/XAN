##Script for the preparation of the datasets to apply the ParaFilter script###
##The name of the HG must be changed in several points of the script. Make sure that all the paths are correct###

##Fetch the sequences##
echo "Fetching the sequences"

mkdir preParaFilter
cp ../seqs/sorted/*.fa preParaFilter/
cd preParaFilter/

echo "Fetch done"

##Alignment with MAFFT LINSI##
echo "Starting alignments with MAFFT"

for f in *; do linsi --localpair --maxiterate 1000 --anysymbol $f > $f.aln; done

mkdir ali
mv *.aln ali/
cd ali/
for f in *; do mv -- "$f" "${f%.fa.aln}_aln.fa"; done

echo "Alignments finished"

##Trimming with trimAl gappyout##
echo "Trimming with trimAl"

for f in *; do trimal -in $f -out $f"_trim.fasta" -gappyout; done

mkdir trimmed
mv *.fasta trimmed/
cd trimmed/
for f in *; do mv -- "$f" "${f%.fa_trim.fasta}_trim.fasta"; done

echo "Trimming finished"

##Inferring ML trees with IQtree##

echo "Inferring ML trees using IQtree"

for f in *; do
iqtree2 -s $f -m C60+LG -B 1000 -rcluster 10 -T AUTO
done

echo "Trees ready"

mkdir trees
mv *.treefile trees/
cd trees/

for f in *; do mv -- "$f" "${f%_aln_trim.fasta.treefile}.tre"; done

##Prepare all files and folders to run ParaFilter##
echo "Prepare files and folders to run ParaFilter"

cd ..
mkdir ParaFilter
cd preParaFilter/
mv *.fa ../ParaFilter/
cd ali/trimmed/trees/
cp *.tre ../../../../ParaFilter/
cd ../../../../ParaFilter/
cp /home/lk19822/Cuvier/Ona/XAN/Phylogenomics/scripts/convertfasta.py .
mkdir fasta_converted

echo "Everything ready"

##Convert fasta files to single line fasta. Here you need to change the path to the correct one!!##
echo "Converting fasta files"


for f in *.fa; do python convertfasta.py -i "/home/lk19822/Cuvier/Ona/XAN/Phylogenomics/Genes/Ancestral_Metazoa/HG_80/analyses/ParaFilter/" -o "/home/lk19822/Cuvier/Ona/XAN/Phylogenomics/Genes/Ancestral_Metazoa/HG_80/analyses/ParaFilter/fasta_converted/" -f $f; done

echo "Fastas converted"

mv *.tre fasta_converted/
cd fasta_converted/
cp /home/lk19822/Cuvier/Ona/XAN/Phylogenomics/scripts/ParaFilter.py .
cp /home/lk19822/Cuvier/Ona/XAN/Phylogenomics/scripts/ParaModules.py .

## Prepare a list with all the tree file names and all the fasta sequences file names ###
echo "Preparing lists of files"

ls *.tre > trees.txt
ls *.fasta > fastas.txt

echo "Lists done"

##Now we run ParaFilter. Fingers crossed, everything is alright##
echo "Run ParaFilter"

python3 ParaFilter.py -f fastas.txt -t trees.txt -n 4 -w /home/lk19822/Cuvier/Ona/XAN/Phylogenomics/Genes/Ancestral_Metazoa/HG_80/analyses/ParaFilter/fasta_converted/

for f in *.new; do mv -- "$f" "${f%.fasta.new}_new.fa"; done

echo "New fasta sequences are pruned and ready to run postParaFilter.sh"
