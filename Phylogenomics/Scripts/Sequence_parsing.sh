CLASS_PATH=/usr/bin/

####All this script is thought to be started/launched at /Cuvier/Ona/XAN/Phylogenetic_Aware_Parsing_Script/OrthoFinder/Scenario4/Input/ folder########
####The name of the folder created has to be changed in all the steps###

###Select HGs from list#######
echo "Selecting HG"

mkdir HG_70
cd HG_70/
while read p; do grep -F -w $p ../Ancestral_Metazoa_genes_IDs2.txt >> "$p".txt; done < ../Ancestral_Metazoa70_list.txt

echo "Done"

for f in *; do cat "$f" | sed 's/\t/\n/g' > "$f"_heads.txt; done
for f in *_heads.txt; do sed -i '1d' "$f"; done
for f in *_heads.txt; do mv -- "$f" "${f%.txt_heads.txt}_heads.txt"; done

###Parse sequences####
echo "Parsing sequences"
mkdir seqs
mv *_heads.txt seqs/
cd seqs/

for f in *; do while read p; do grep -F $p "$f" >> "$f"_sorted.txt; done < /home/lk19822/Cuvier/Ona/XAN/Phylogenetic_Aware_Parsing_Script/OrthoFinder/Metazoan_Phylogenomics_list_all.txt ; done

mkdir sorted
mv *_sorted.txt sorted/
cd sorted/

for f in *; do mv -- "$f" "${f%_heads.txt_sorted.txt}_heads_sorted.txt"; done

for f in *; do
while read p; do grep -A 1 -F $p ../../../Proteomes_XAN_modif2.fasta >> "$f"_seqs.fa; done < "$f"
done

for f in *_seqs.fa; do mv -- "$f" "${f%_heads_sorted.txt_seqs.fa}_seqs.fa"; done

echo "Done"

###Save headers for future identification"
echo "Saving headers"

#for f in *.fa; do cat "$f" | grep "^>" > "$f"_headers.txt; done
#for f in *.txt; do mv -- "$f" "${f%_seqs.fa_headers.txt}_headers.txt"; done
#for f in *.fa; do cat "$f" | sed 's/|\([A-Z][a-z]*\).*//g' | awk '/^>/ { $0=$0 "|" ++i }1' > $f"_rn.fa"; done
for f in *.fa; do cat "$f" | sed 's/^>\([A-Z][a-z]*\)/>\1|\1/g' | sed 's/|\([A-Z][a-z]*\).*//g' | awk '/^>/ { $0=$0 "|" ++i }1' > $f"_rn.fa"; done
for f in *_rn.fa; do cat "$f" | grep "^>" > "$f"_headers2.txt; done
for f in *_headers2.txt; do mv -- "$f" "${f%_seqs.fa_rn.fa_headers2.txt}_headers2.txt"; done
for f in *_heads_sorted.txt; do ORTHOGROUP=`echo $f | cut -d . -f 1 | sed 's/_heads_sorted//g'`; echo $ORTHOGROUP; file2="$ORTHOGROUP"_headers2.txt; paste -d' ' $f $file2 > "$ORTHOGROUP"_Final_headers.tsv; done

for f in *_rn.fa; do mv -- "$f" "${f%.fa_rn.fa}_rn.fa"; done
echo "Done"

###Selecting species####
echo "Species selection"

for f in *_rn.fa; do cat "$f" | sed 's/|[0-9]*//g' | grep "^>" | sort | uniq -c | sort -nr > $f"_count_list.txt"; done
for f in *count_list.txt; do cat "$f" | grep "1 >[A-Z][a-z]*" | sed 's/1 >//g' > $f"_unique_seqs.txt"; done
for f in *_unique_seqs.txt; do mv -- "$f" "${f%_seqs_rn.fa_count_list.txt_unique_seqs.txt}_uniq_seqs.txt"; done
for f in *_list.txt; do mv -- "$f" "${f%_seqs_rn.fa_count_list.txt}_count_list.txt"; done

mkdir raw
mv *_seqs.fa raw/

for f in *.fa; do mv -- "$f" "${f%_seqs_rn.fa}.fa"; done
for f in *.fa; do
OG=`echo $f | cut -d . -f 1 | sed 's/.\+\///g'`
echo "$OG"
while read p; do grep -A 1 -F $p $f >> $OG"_selec_seqs.fa"; done < $OG"_uniq_seqs.txt"; done

mkdir selec
mv *_selec_seqs.fa selec/

echo "Selection done"

cd ../../../
mv HG_70/ ../../../../Phylogenomics/Genes/Ancestral_Metazoa/
cd ../../../../Phylogenomics/Genes/Ancestral_Metazoa/HG_70/
mkdir analyses
cd seqs/sorted/selec/
cp * ../../../analyses/
cd ../../../analyses/

##Alignment with MAFFT##
#http://mafft.cbrc.jp/alignment/software/

echo "Aligning sequences using MAFFT (linsi)..."
for f in *; do
linsi --localpair --maxiterate 1000 $f > $f.aln
done

for f in *.aln; do mv -- "$f" "${f%.fa.aln}_aln.fa"; done

mkdir ali
mv *_aln.fa ali/
cd ali/
echo "Done"


#Trimming with trimAl#
#echo "Trimming with relaxed parameters"
#for f in *; do trimal -in "$f" -out "$f".trim.fa -gt 0.3 -st 0.01 -cons 60 -noallgaps; done

echo "Triming with gappyout option"
for f in *; do trimal -in "$f" -out "$f".trim.phy -gappyout -phylip; done
for f in *.phy; do mv -- "$f" "${f%_selec_seqs_aln.fa.trim.phy}.phy"; done

mkdir trimmed
mv *.phy trimmed/
cd trimmed/
echo "Trimming finished"

###Writing partitions###
echo "Write partitions"

for f in *; do
OG=`echo $f | cut -d . -f 1 | sed 's/.\+\///g'`
ALIGNMENT_LENGTH=`head -1 "$OG".phy | awk -F " " '{print $2}'`
NUM_SPP=`head -1 "$OG".phy | awk -F " " '{print $1}'`
echo $OG = $ALIGNMENT_LENGTH , $NUM_SPP species >> partition_data.txt
done

echo "Partitions ready"

#Infer phylogeny with IQtree#
echo "Infer phylogenies with IQtree"

for f in *.phy; do
iqtree2 -s "$f" -m TEST -B 1000 -T AUTO
done

mkdir trees
cp *.treefile trees/
cd trees/

for f in *; do mv -- "$f" "${f%.phy.treefile}.tre"; done
echo "Trees inferred, what do you want to do now?"

echo "Concatenate for further analyses"

cd ../
mkdir conc
cd ../
cp *.fa trimmed/conc/
cd trimmed/conc/

for f in *.fa; do sed -i 's/|[0-9]*//g' "$f"; done
for f in *.fa; do trimal -in "$f" -out $f"_trim.fasta" -gappyout; done
rm *.fa
for f in *.fasta; do mv -- "$f" "${f%.fa_trim.fasta}_trim.fasta"; done

cp /home/lk19822/Programes/FASconCAT-G/FASconCAT-G_v1.04.pl .

echo "Now execute FASconCAT-G with options l and s and continue with the analyses. Good luck!"

#Runs PhyloTreePruner
	#Usage: java PhyloTreePruner input_tree_file min_number_of_taxa input_fasta bootstrap_cutoff r/u
	#r = redundant, let SCaFoS pick the best sequence for a given OTU
	#u = unique, pick the longest sequence for a given OTU
	#ex: java PhyloTreePruner 0001.tre 10 0001.fa 0.5 r
#cd ../../
#cp * conversed/trees/
#cd conversed/trees/

#cd ../../../../raw/
#cp *_seqs.fa ../trimmed/renamed/conversed/trees/
#cd ../trimmed/renamed/conversed/trees/

#for f in *_seqs.fa; do cat "$f" | sed 's/|/@/g' > "$f"_rn.fa; done

#for f in *.fa; do mv -- "$f" "${f%_seqs_aln_trim_rn.fa}.fa"; done
#for f in *.fa; do mv -- "$f" "${f%_seqs.fa_rn.fa}.fa"; done
#for f in *.tre; do cat "$f" | sed 's/|/@/g' > "$f"_rn.tre; done

#for f in *.tre;
#do
#sed -i 's/_:/:/g' $f
#OG=`echo $f | cut -d . -f 1 | sed 's/.\+\///g'`
#        echo $OG
###################################################################################
####You may want to change some of the PhyloTreePruner settings in the next line###
###################################################################################

#java -cp /home/lk19822/Programes/src_and_wrapper_scripts/ PhyloTreePruner $OG".tre_rn.tre" 10 $OG".fa" 0.1 u
#done
#echo Done

#mkdir pruned
#mv *_pruned.fa pruned/
#cd pruned/

#for f in *; do mv -- "$f" "${f%.fa_pruned.fa}_pruned.fa"; done

