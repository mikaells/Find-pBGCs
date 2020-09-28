#!/bin/bash

####
#A program to find BGCs in phages.
#
#depends on
# ProphET
# antismash
# genbank_to_fasta.py
# bioawk
####



#the function we use for the loop in the main program
mainFunc () {

	first_dir=$1
	parent_dir=$2
	
	#clear up old
	if [ -f $parent_dir/$first_dir/*.log ] ; then
		rm -f $parent_dir/$first_dir/*.log
		rm -f $parent_dir/$first_dir/*.new.gff
		rm -f $parent_dir/$first_dir/*phages_coords 
		rm -f $parent_dir/$first_dir/*prophages_w_BGC.fna 
		rm -f $parent_dir/$first_dir/all_prophages.fna
		rm -r $parent_dir/$first_dir/gBGC_folder
		rm -f $parent_dir/$first_dir/gBGC_summary.txt
		rm -f $parent_dir/$first_dir/pBGC.headers
		rm -r $parent_dir/$first_dir/pBGC_folder
		rm -f $parent_dir/$first_dir/pBGC_summary.txt
		rm -r $parent_dir/$first_dir/prophage_folder
		rm -f $parent_dir/$first_dir/pBGCs.fna
	fi

	#remove MD5SUMS
	if [ -f $parent_dir/$first_dir/MD5SUMS ] ; then
		rm $parent_dir/$first_dir/MD5SUMS
	fi
	 
	#unzip files if they are zipped
	if [ -f $parent_dir/$first_dir/*fna.gz ]
	then
		gunzip $parent_dir/$first_dir/*
	fi
	
	#make names
	fasta=$(ls $parent_dir/$first_dir/ | grep "fna");
	gff=$(ls $parent_dir/$first_dir/ | grep "gff");
	b_name="${fasta%.*}"
	
	
	#Run antismash on genomes
	antismash  --outputfolder $parent_dir/$first_dir/gBGC_folder -c 1  --minimal --disable-embl  --disable-svg --disable-xls --disable-html $parent_dir/$first_dir/$b_name.fna >> $parent_dir/$first_dir/$b_name.log 2>&1

	#Check if there if there are no gBGCs and assign 0, else assign count
	if [ -z "$(ls -A $parent_dir/$first_dir/gBGC_folder/txt/)" ]; then
		gBGC_hits=0
		if [ "$clean_up" = true ] ; then rm -r $parent_dir/$first_dir/gBGC_folder ; fi
	else 
		cat $parent_dir/$first_dir/gBGC_folder/txt/*_BGC.txt >> $parent_dir/$first_dir/gBGC_summary.txt
		sed -i "/BGC/d" $parent_dir/$first_dir/gBGC_summary.txt
		if [ "$clean_up" = true ] ; then rm -r $parent_dir/$first_dir/gBGC_folder ; fi
		gBGC_hits=$(<"$parent_dir/$first_dir/gBGC_summary.txt" wc -l )
		
	fi
	
	#Reformat GFF for ProphET
	~/ProphET/UTILS.dir/GFFLib/gff_rewrite.pl --input $parent_dir/$first_dir/$gff -output $parent_dir/$first_dir/$b_name.new.gff --add_missing_features  >> $parent_dir/$first_dir/$b_name.log  2>&1 
	
	#Run ProphET
	~/ProphET/ProphET_standalone.pl --fasta $parent_dir/$first_dir/$b_name.fna --gff_in $parent_dir/$first_dir/$b_name.new.gff --outdir $parent_dir/$first_dir/prophage_folder  >> $parent_dir/$first_dir/$b_name.log 2>&1
	
	#Assign prophage hits
	prophage_hits=$(<"$parent_dir/$first_dir/prophage_folder/phages_coords" wc -l )
	
	#If there are prophages AND gBGCs run antismash on the prophages
	if [ "$prophage_hits" -gt 0 ] && [ "$gBGC_hits" -gt 0 ]; then
		
		#put all prophages in a prophage file
		cat $parent_dir/$first_dir/prophage_folder/*fas > $parent_dir/$first_dir/all_prophages.fna

		#rename fasta headers in phages, doesnt help
		#sed -i 's/>phage_\([0-9]*\):[0-9-]*:\([A-Za-z0-9_.]*\)/>\2:p\1/' $parent_dir/$first_dir/all_prophages.fna

		#Run antismash on the prophages
		antismash --outputfolder $parent_dir/$first_dir/pBGC_folder -c 1  --minimal  --disable-embl  --disable-svg --disable-xls --disable-html $parent_dir/$first_dir/all_prophages.fna >> $parent_dir/$first_dir/$b_name.log 2>&1
				
		#If there are pBGCs, move summary to top folder		
		if [  ! -z "$(ls -A $parent_dir/$first_dir/pBGC_folder/txt/)" ]; then
			
			#copy all pBGCs summaries to one file
			cat $parent_dir/$first_dir/pBGC_folder/txt/*_BGC.txt >> $parent_dir/$first_dir/pBGC_summary.txt
			sed -i "/BGC/d" $parent_dir/$first_dir/pBGC_summary.txt
						
			if [ "$clean_up" = true ] ; then rm -r $parent_dir/$first_dir/pBGC_folder ; fi
			
			#count pBGCs
			pBGC_hits=$(<"$parent_dir/$first_dir/pBGC_summary.txt" wc -l )

			
			#make pBGC types, does not work
			pBGC_types=$( cat $parent_dir/$first_dir/pBGC_summary.txt | awk -vORS=, '{print $2}' | sed 's/,$//' )
			
			#Saving all phages with BGC as fasta
			head  -n 2 $parent_dir/$first_dir/pBGC_folder/*cluster* | grep DEFINITION |  awk '{print $2}' | uniq > $parent_dir/$first_dir/pBGC.headers
			seqtk subseq $parent_dir/$first_dir/all_prophages.fna $parent_dir/$first_dir/pBGC.headers > $parent_dir/$first_dir/$b_name.prophages_w_BGC.fna 
			
			
			#writing the individual prophages
			COUNTER=0;
			for k in $(cat $parent_dir/$first_dir/pBGC.headers); 
			do 
				let COUNTER=COUNTER+1 ; 
				echo $k > $parent_dir/$first_dir/temp; 
				seqtk subseq $parent_dir/$first_dir/all_prophages.fna $parent_dir/$first_dir/temp > all_prophages_w_BGC/$b_name.cluster00$COUNTER.fna;
				rm $parent_dir/$first_dir/temp;
			done 
			
			
			#saving pBGC coordinates
			pBGC_coordinates=$( cat $parent_dir/$first_dir/pBGC.headers | cut -d":" -f 2 | paste -s -d, - )
			
			#saving genome length
			bioawk -c fastx '{ print $name, length($seq) }' < $parent_dir/$first_dir/$b_name.fna > $parent_dir/$first_dir/length_stats.txt
			genome_length=$(cat $parent_dir/$first_dir/pBGC.headers | cut -d":" -f3 | grep -f -  $parent_dir/$first_dir/length_stats.txt | cut -f2 | paste -s -d, -) 
			
			#Saving all pBGC fna
			for k in $parent_dir/$first_dir/pBGC_folder/*cluster*;
			do 
				clus_num=$(echo $(basename $k) |cut -d"." -f4)
				
				#converting to fasta
				#note that this script insist on outputting the fasta in
				#in the same as the input, requiring som messing around with the -o
				~/Scripts/genbank_to_fasta.py  -i $k -o ../../../../all_pBGCs/$b_name.$clus_num.fna -s 'whole' -u $b_name.$clus_num
			done
			
		else 
			#if no pBGC, then nothing
			pBGC_hits=0
			pBGC_types="-"
			pBGC_coordinates="-"
			genome_length=0
			if [ "$clean_up" = true ] ; then rm -r $parent_dir/$first_dir/pBGC_folder; fi
		fi
		
	else
		#If no phages, move on
		pBGC_hits=0
		pBGC_types="-"
		pBGC_coordinates="-"
		genome_length=0
	fi
	
	#Write out to terminal
	echo -e "\n\t***$b_name had $gBGC_hits BGCs, $prophage_hits prophages and $pBGC_hits pBGCs 	\n"
	
	#Make taxonomic name for logfile logfile
	tax_name=$( head  $parent_dir/$first_dir/$b_name.fna -n 1 | sed 's/>//' )
	
	#Append to log-file
	echo -e "$b_name\t$gBGC_hits\t$prophage_hits\t$pBGC_hits\t$pBGC_types\t$pBGC_coordinates\t$genome_length\t$tax_name" >> all.log

	#Keep phage coordinates
	cp $parent_dir/$first_dir/prophage_folder/phages_coords $parent_dir/$first_dir/$b_name.phages_coords
	
	to_remove=$(ls -d $parent_dir/$first_dir/prophage_folder/)
	if [ "$clean_up" = true ] ; then rm -r $to_remove; fi
	if [ "$clean_up" = true ] ; then rm $parent_dir/$first_dir/$b_name.new.gff; fi
	
	#zipping start files
	gzip $parent_dir/$first_dir/$gff
	gzip $parent_dir/$first_dir/$fasta
	
}
export -f mainFunc

#open logfile
echo -e "#ID\tgenome_anti_hits\tprophage_hits\tpBGC_hits\tpBGC_types\tname" > all.log

#input directory must be the input folders (both!) from ncbi_genome_download
input_dir=$1

#make folder for all prophages with BGCs
rm -r all_prophages_w_BGC 
rm -r all_pBGCs
rm -r all.log

#open logfile
echo -e "#ID\tgenome_anti_hits\tprophage_hits\tpBGC_hits\tpBGC_types\tname" > all.log

#make folders
mkdir all_prophages_w_BGC/
mkdir all_pBGCs/
#whether or not temporary files should be deleted
clean_up_arg=true
clean_up=${2:-$clean_up_arg}  

#main function, takes the folders contained in the input as well as the input
parallel --linebuffer mainFunc :::  $(ls $input_dir) ::: $input_dir
