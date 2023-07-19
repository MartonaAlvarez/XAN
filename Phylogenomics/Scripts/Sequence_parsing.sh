CLASS_PATH=/usr/bin/

####All this script is thought to be started/launched in the Input folder from the Phylogenetic Aware Parsing Script (Paps & Holland, 2018). Change the path ########
####when necessary. The name of the folder created has to be changed in all the steps. Script created by Marta Alvarez Presas at the University of Bristol with ###
####the purpose of parsing the sequences from the Phylogenetic Aware Parsing Script desired for downstream analyses. ###


###Select HGs from list#######
echo "Selecting HG"

mkdir HG_70
cd HG_70/
while read p; do grep -F -w $p ../Ancestral_Metazoa_genes_IDs.txt >> "$p".txt; done < ../Ancestral_Metazoa70_list.txt

echo "Done"

for f in *; do cat "$f" | sed 's/\t/\n/g' > "$f"_heads.txt; done
for f in *_heads.txt; do sed -i '1d' "$f"; done
for f in *_heads.txt; do mv -- "$f" "${f%.txt_heads.txt}_heads.txt"; done

###Parse sequences####
echo "Parsing sequences"
mkdir seqs
mv *_heads.txt seqs/
cd seqs/

for f in *; do while read p; do grep -F $p "$f" >> "$f"_sorted.txt; done < /path/to/Metazoan_Phylogenomics_list_all.txt ; done

mkdir sorted
mv *_sorted.txt sorted/
cd sorted/

for f in *; do mv -- "$f" "${f%_heads.txt_sorted.txt}_heads_sorted.txt"; done

for f in *; do
while read p; do grep -A 1 -F $p ../../../Proteomes_XAN.fasta >> "$f"_seqs.fa; done < "$f"
done

for f in *_seqs.fa; do mv -- "$f" "${f%_heads_sorted.txt_seqs.fa}_seqs.fa"; done

echo "Done"

###Save headers for future identification"
echo "Saving headers"

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

