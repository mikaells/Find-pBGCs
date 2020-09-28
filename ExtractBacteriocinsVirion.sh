for i in only_pBGC_phages/anti/*;
do

	folder="${i##*/}"
	echo $folder
	#extract proteins from antismash pBGC
	~/Scripts/genbank_to_fasta.py -i $i/*cluster00*.gbk -s 'aa' -q  'locus_tag,sec_met,gene,product,location'  -o $folder.AA_reformat.faa
	sed -i 's/ /_/g' $i/$folder.AA_reformat.faa

	faidx --regex "Type" -f $i/$folder.AA_reformat.faa > ~/phage_genomes/only_pBGC_phages/only_bacteriocin_faa/$folder.faa

	#extract nucleotides from antismash pBGC
	~/Scripts/genbank_to_fasta.py -i $i/*cluster00*.gbk -s 'nt' -q  'locus_tag,sec_met,gene,product,location'  -o $folder.NT_reformat.fna
	sed -i 's/ /_/g' $i/$folder.NT_reformat.fna
	faidx --regex "Type" -f $i/$folder.NT_reformat.fna > ~/phage_genomes/only_pBGC_phages/only_bacteriocin_fna/$folder.fna
done
