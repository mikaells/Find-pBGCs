#!/bin/bash

####
#A program to find BGCs in phages.
#18-02-20
#depends on
# blast
####

#blasts each pBGC against a database of all bacteriocins, followed by filtering for identity
#returns a file for each pBGC with matches better than 50% and longer than 50% coverage
#


#the function we use for the loop in the main program
mainFunc () {

	first_dir=$1
	parent_dir=$2
	fna_file=$(ls $2/$1/ | grep "genomic.fna") 
	b_name=$(basename "$fna_file" | cut -d. -f1,2)

	echo -e "to blast $fna_file, $b_name\n";	
	cat $2/$1/pBGC_folder/AA_reformat.fasta |  blastp -subject ~/all_phageBGC/BACTIBASE.txt  -qcov_hsp_perc 50 -outfmt 6 | awk '$3>50' > bacteriocin_vs_pBGC_id50c50/$b_name._c50_id50.blast

}
export -f mainFunc

input_dir=$1
#main function, takes the folders contained in the input as well as the input
parallel --linebuffer mainFunc :::  $(ls $input_dir) ::: $input_dir

