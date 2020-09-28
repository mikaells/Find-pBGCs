#!/bin/bash

####
#A program to re-find core pBGC genes in all genomes
#
#depends on
# blastn
####




#the function we use for the loop in the main program
mainFunc () {

	first_dir=$1
	parent_dir=$2
	gz_file=$2/$1/*genomic.fna.gz

	echo -e "to blast $gz_file, $b_name\n";
	zcat $gz_file | blastn -subject ~/all_phageBGC/for_blast/all_core_genes.fna -outfmt 6 > ~/all_phageBGC/coreGenesVsGenomesBlastn/$1.blastn
}
export -f mainFunc
rm -r ~/all_phageBGC/coreGenesVsGenomesBlastn
mkdir ~/all_phageBGC/coreGenesVsGenomesBlastn

input_dir=$1
#main function, takes the folders contained in the input as well as the input
parallel --linebuffer mainFunc :::  $(ls $input_dir) ::: $input_dir
