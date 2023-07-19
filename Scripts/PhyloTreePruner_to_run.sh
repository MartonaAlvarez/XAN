CLASS_PATH=/usr/bin/

## This script prepares all the necessary files for pruning paralogs from multi-fasta alignments using Phylotreepruner. Note that all the necessary software needs to 
## be properly installed and that the used will need to change all the corresponding paths. The steps of changing and saving headers can be skipped, and are only for making
## it easier if the user wants to trace the initial steps. Created by Marta Alvarez Presas at the University of Bristol. If used, please, cite.


###Prepare files for further analysis#######
echo "Preparing files for further analysis"

for f in *.fa; do cat "$f" | grep "^>" > "$f"_headers.txt; done

for f in *.txt; do mv -- "$f" "${f%_seqs.fa_headers.txt}_headers.txt"; done

for f in *.fa; do cat "$f" | sed 's/|\([A-Z][a-z]*\).*//g' | awk '/^>/ { $0=$0 "|" ++i }1' > $f"_rn.fa"; done

for f in *_rn.fa; do cat "$f" | grep "^>" > "$f"_headers2.txt; done

for f in *_headers2.txt; do mv -- "$f" "${f%_seqs.fa_rn.fa_headers2.txt}_headers2.txt"; done

for f in *_headers.txt; do ORTHOGROUP=`echo $f | cut -d . -f 1 | sed 's/_headers//g'`; echo $ORTHOGROUP; file2="$ORTHOGROUP"_headers2.txt; paste -d' ' $f $file2 > "$ORTHOGROUP"_Final_headers.tsv; done

mkdir ali
mv *_rn.fa ali/
cd ali/

for f in *.fa; do mv -- "$f" "${f%.fa_rn.fa}.fa"; done

#Aligns sequences using MAFFT.
#http://mafft.cbrc.jp/alignment/software/

echo "Aligning sequences using MAFFT (linsi)..."
for FILENAME in *.fa; do
linsi --localpair --maxiterate 1000 $FILENAME > $FILENAME.aln
done

mkdir raw
mv *.fa raw/

for f in *.aln; do mv -- "$f" "${f%.fa.aln}_aln.fa"; done
echo Done

#Trimming with trimAl#
echo "Trimming with relaxed parameters"
for FILENAME in *.fa; do trimal -in "$FILENAME" -out "$FILENAME".trim.fa -gt 0.3 -st 0.01 -cons 60; done

for ali in *.trim.fa; do mv -- "$ali" "${ali%_aln.fa.trim.fa}_aln_trim.fa"; done
mkdir trimmed
mv *_trim.fa trimmed/
echo "Trimming finished"

#Converts input fasta files to relaxed phylip files.
#http://sco.h-its.org/exelixis/software/fasta2relaxedPhylip.pl

echo "Converting fasta files to relaxed phylip files for IQtree..."
cd trimmed/

for ali in *.fa; do cat "$ali" | sed 's/|/@/g' > "$ali"_rn; done

for f in *_rn; do mv -- "$f" "${f%.fa_rn}_rn.fa"; done

mkdir renamed
mv *_rn.fa renamed/
cd renamed/

for f in *.fa; do perl /path/to/fasta2relaxedPhylip.pl -f "$f" -o "$f".phylip; done

mkdir conversed
mv *.phylip conversed/
cd conversed/

for ali in *.phylip; do cat "$ali" | sed 's/@/|/g' > "$ali".phy; done
for f in *.phy; do mv -- "$f" "${f%_seqs_aln_trim_rn.fa.phylip.phy}.phy"; done

echo "Conversion done"

for FILENAME in *.phy; do
ORTHOLOGY_GROUP=`echo $FILENAME | cut -d . -f 1 | sed 's/.\+\///g'`
ALIGNMENT_LENGTH=`head -1 "$ORTHOLOGY_GROUP".phy | awk -F " " '{print $2}'`
echo $ORTHOLOGY_GROUP = $ALIGNMENT_LENGTH >> partition_data.txt
done

#Infer phylogeny with IQtree#
echo "Infer phylogenies with IQtree"

for ali in *.phy; do
iqtree2 -s "$ali" -m TEST -B 1000 -T 25
done

mkdir trees
cp *.treefile trees/
cd trees/

for f in *; do mv -- "$f" "${f%.phy.treefile}.tre"; done

#Runs PhyloTreePruner
	#Usage: java PhyloTreePruner input_tree_file min_number_of_taxa input_fasta bootstrap_cutoff r/u
	#r = redundant, let SCaFoS pick the best sequence for a given OTU
	#u = unique, pick the longest sequence for a given OTU
	#ex: java PhyloTreePruner 0001.tre 10 0001.fa 0.5 r
#cd ../../
#cp * conversed/trees/
#cd conversed/trees/

cd ../../../../raw/
cp *_seqs.fa ../trimmed/renamed/conversed/trees/
cd ../trimmed/renamed/conversed/trees/

for f in *_seqs.fa; do cat "$f" | sed 's/|/@/g' > "$f"_rn.fa; done

#for f in *.fa; do mv -- "$f" "${f%_seqs_aln_trim_rn.fa}.fa"; done
for f in *.fa; do mv -- "$f" "${f%_seqs.fa_rn.fa}.fa"; done
for f in *.tre; do cat "$f" | sed 's/|/@/g' > "$f"_rn.tre; done

for f in *.tre;
do
sed -i 's/_:/:/g' $f
OG=`echo $f | cut -d . -f 1 | sed 's/.\+\///g'`
        echo $OG
###################################################################################
####You may want to change some of the PhyloTreePruner settings in the next line###
###################################################################################

java -cp /path/to/ PhyloTreePruner $OG".tre_rn.tre" 10 $OG".fa" 0.1 u
done
echo Done

mkdir pruned
mv *_pruned.fa pruned/
cd pruned/

for f in *; do mv -- "$f" "${f%.fa_pruned.fa}_pruned.fa"; done

