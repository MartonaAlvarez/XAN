CLASS_PATH=/usr/bin/

###Prepare files for further analysis after running ParaFilter. Script created by Marta Alvarez Presas at the University of Bristol thanks to Mattia Giacomelli,
### who shared ParaFilter. #######


echo "Preparing files for further analysis"

mkdir post_pruning
cp *_new.fa post_pruning/
cd post_pruning/

for f in *; do mv -- "$f" "${f%_twoline_new.fa}_pruned.fa"; done

###Selecting species####
echo "Species selection"

for f in *; do cat $f | sed 's/|[0-9]*//g' | grep "^>" | sort | uniq -c | sort -nr > $f"_count_list.txt"; done
for f in *.txt; do cat $f | grep "1 >[A-Z][a-z]*" | sed 's/1 >//g' > $f"_unique_seqs.txt"; done
for f in *_seqs.txt; do mv -- "$f" "${f%.fa_count_list.txt_unique_seqs.txt}_uniq_seqs.txt"; done
for f in *_list.txt; do mv -- "$f" "${f%.fa_count_list.txt}_count_list.txt"; done

for f in *.fa; do
OG=`echo $f | cut -d . -f 1 | sed 's/.\+\///g'`
echo "$OG"
while read p; do grep -A 1 -F $p $f >> $OG"_selec_seqs.fa"; done < $OG"_uniq_seqs.txt"; done

echo "Selection done"

mkdir selec
mv *_seqs.fa selec/
cd selec/

mkdir ali

#Aligns sequences using MAFFT.
#http://mafft.cbrc.jp/alignment/software/

echo "Aligning sequences using MAFFT (linsi)..."

for f in *.fa; do
linsi --localpair --maxiterate 1000 --anysymbol $f > ali/$f.aln
done

cd ali/
for f in *; do mv -- "$f" "${f%.fa.aln}_aln.fa"; done
echo Done

#Trimming with trimAl#
echo "Trimming with automatic parameters"

for f in *; do trimal -in $f -out $f"_trim.fa" -gappyout; done

for f in *_trim.fa; do mv -- "$f" "${f%.fa_trim.fa}_trim.fa"; done

mkdir trimmed
mv *_trim.fa trimmed/

echo "Trimming finished"

cd trimmed/

##Rename sequences for the concatenated dataset##

echo "Rename sequences"

for f in *; do cat $f | sed 's/|[0-9]*//g' > $f"_rn.fasta"; done
for f in *.fasta; do mv -- "$f" "${f%.fa_rn.fasta}_rn.fasta"; done

mkdir rn
mv *_rn.fasta rn/
cd rn/

echo "All sequences have been renamed"

##Concatenate sequences##
echo "Concatenate for further analyses"

cp /path/to/FASconCAT-G/FASconCAT-G_v1.04.pl .

perl FASconCAT-G_v1.04.pl -l -s
mv FcC_supermatrix.fas HG_80_ParaFilter_supermatrix.fas
cat FcC_supermatrix_partition.txt | sed 's/_pruned_selec_seqs_aln_trim_rn.fasta//g' > HG_80_ParaFilter_partitions.txt


##Infer phylogeny with IQtree#
echo "Infer phylogenies with IQtree"

iqtree2 -s HG_80_ParaFilter_supermatrix.fas -m LG+C60 -B 1000 -spp HG_80_ParaFilter_partitions.txt -rcluster 10 -T AUTO

echo "Tree ready, have a nice day"
