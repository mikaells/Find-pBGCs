#!/bin/bash

#must set perl dir from conda env folder
# export PERL5LIB=`pwd`/lib/site_perl/5.26.2:`pwd`/lib/5.26.2

#Run with argument ~/all_phageBGC/just_genomes_w_pBGCs/refseq/bacteria/
#is probably already gunzipped

#1: annotate genome
#2: fish out proteins from antismash proteins with genbank_to_fasta.py
#3: rename file
#4: fish out everything with a proper name, in practice bacteriocins, from antismash protein file and put in only_bacteriocin_faa/ folder
#5: blast bacteriocin(s) from 4: against its parent genome and put it in only_bacteriocin_blast/ folder
#6: repeat  with nucleotides

mainFunc () {

	first_dir=$1
	parent_dir=$2
	
	echo -e "*****Running $1\n\n"
	#Unzip and annotate genome file
	
	if [ -f $2/$1/*fna.gz ]; then
		echo "suge"
		gunzip $2/$1/*fna.gz
	fi

	###
	#rm -r $2/$1/annot_genome
	prokka --outdir $2/$1/annot_genome --prefix $1 --quiet  --cpus 1 $2/$1/*_genomic.fna 
	sed -i 's/ /_/g' $2/$1/annot_genome/*.faa
	sed -i 's/ /_/g' $2/$1/annot_genome/*.ffn
	
	#FIRST PROTEINS
	
	#extract proteins from antismash pBGC
	~/Scripts/genbank_to_fasta.py -i $2/$1/pBGC_folder/*cluster00*.gbk -s 'aa' -q  'locus_tag,sec_met,gene,product,location'  -o $i.AA_reformat.faa
	#rename outfile for faidx compatibility
	sed -i 's/ /_/g' $2/$1/pBGC_folder/$i.AA_reformat.faa 
	
	#write faa file to common folder
	faidx --regex "Type" -f $2/$1/pBGC_folder/$i.AA_reformat.faa > ~/all_phageBGC/only_bacteriocin_faa/$1.faa
	#blast against parent and add to both parent folder and common folder
	faidx --regex "Type" -f $2/$1/pBGC_folder/$i.AA_reformat.faa | blastp -subject $2/$1/annot_genome/*.faa -outfmt 6 -qcov_hsp_perc 50 | awk '{ if ($3 > 50) { print } }' | tee ~/all_phageBGC/only_bacteriocin_blastp/$1.blastp $2/$1/genome_pBacteriocin_hits.blastp    

	#REPEAT WITH NT
	
	#extract nucleotides from antismash pBGC
	~/Scripts/genbank_to_fasta.py -i $2/$1/pBGC_folder/*cluster00*.gbk -s 'nt' -q  'locus_tag,sec_met,gene,product,location'  -o $i.NT_reformat.fna
	#rename outfile for faidx compatibility
	sed -i 's/ /_/g' $2/$1/pBGC_folder/$i.NT_reformat.fna
	
	#write fna file to common folder
	faidx --regex "Type" -f $2/$1/pBGC_folder/$i.NT_reformat.fna > ~/all_phageBGC/only_bacteriocin_fna/$1.fna
	#blast against parent and add to both parent folder and common folder
	faidx --regex "Type" -f $2/$1/pBGC_folder/$i.NT_reformat.fna | blastn -subject $2/$1/annot_genome/*.ffn -outfmt 6 -qcov_hsp_perc 50 -max_hsps 1 -perc_identity 70 | tee ~/all_phageBGC/only_bacteriocin_blastn/$1.blastn $2/$1/genome_pBacteriocin_hits.blastn    

}
export -f mainFunc

#remove previous and make new
rm -r ~/all_phageBGC/only_bacteriocin_*
mkdir ~/all_phageBGC/only_bacteriocin_blastp/
mkdir ~/all_phageBGC/only_bacteriocin_blastn/
mkdir ~/all_phageBGC/only_bacteriocin_faa/
mkdir ~/all_phageBGC/only_bacteriocin_fna/



input_dir=$1

#mainFunc  "GCF_000006725.1" $input_dir 

#main function, takes the folders contained in the input as well as the input
parallel --linebuffer mainFunc :::  $(ls $input_dir) ::: $input_dir

#parallel --gnu --linebuffer mainFunc ::: "GCF_000006725.1" "GCF_000827065.1"   ::: $input_dir



#rename faa's headers by filename	
#for i in *.faa; do gawk -i inplace '/>/{sub(">","&"FILENAME"_");sub(/\.faa/,x)}1' $i ;done 
 


#for mange blastp hits, måske filter med awk
#kan ikke se navn fra prokka
 
 
 