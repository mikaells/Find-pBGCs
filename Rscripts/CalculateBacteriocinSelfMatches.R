rm(list=ls())
library(stringr)
library(spatstat.utils)

#1) is it correctly matching itself
#2) are the pBGCs have more copies in the same genome
#3) Are the pBGCs found in other genomes without pBGCs


verbose=F

#read in log file (remember quote option due to internal "'" in some names)
names=read.table("c:/Users/milst/Desktop/all_phageBGC/all.log", header = F,sep = "\t", quote="", comment.char = "#", stringsAsFactors = F )
colnames(names)=c("ID","gBGC_hits","prophage_hits","pBGC_hits","pBGC_types","range","genome_len","name")

#Make IDs and species for easier access
names$shortID=gsub("^([A-Z0-9_]*.[0-9])_.*","\\1", names$ID)
names$species=apply(str_split_fixed(names$name, " ", 8)[,2:3],MARGIN = 1, function(x) paste(x, collapse = " "))
names$species=factor(names$species)


#Get BLAST filenames
blast_files=dir("onlyBacteriocins/only_bacteriocin_blastp/", full.names = T)


i=20

collection=data.frame(name=character(23), nQ=numeric(23),nS=numeric(23),pBGC=numeric(23), stringsAsFactors = F)

counter=1
#Main loop
for(i in 1:length(blast_files)) {
  
  #get genome id and find it in log file
  genome_id=gsub("^onlyBacteriocins/only_bacteriocin_blastp/([A-Z0-9_]*.[0-9]).*","\\1", blast_files[i])
  genome_indx=grep(genome_id, names$shortID)
  genome_entry=names[genome_indx,]
  
  #move to next file if empty blast file
  #cat(paste(i," "))
  if(file.info(blast_files[i])$size==0){
    if(verbose)cat(paste(genome_id,"has no blast hits ( nr",i,")\n"))
    next
  }
  
  #read blast file and put columns on
  blast=read.table(blast_files[i], sep="\t", stringsAsFactors = F)
  colnames(blast) = c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")
  
  
  unq_qseq=unique(blast$qseqid)
  unq_sseq=unique(blast$sseqid)
  if(length(unq_qseq) != length(unq_sseq) ) {
    
    if(length(unq_sseq) == genome_entry$pBGC_hits ) {
      cat(paste(i,"has", length(unq_qseq), "queries,",length(unq_sseq),"subjects,", genome_entry$pBGC_hits,"pBGCs","\n\tItself", "\n\n"))
    } else {
      cat(paste(i,"has", length(unq_qseq), "queries,",length(unq_sseq),"subjects,", genome_entry$pBGC_hits,"pBGCs", "\n\n"))
      
      collection[counter,1]=as.character(genome_entry$species)
      collection[counter,2]=length(unq_qseq)
      collection[counter,3]=length(unq_sseq)
      collection[counter,4]=genome_entry$pBGC_hits
      counter=counter+1 
    }  
  }
}



table(collection$name)
barplot(table(collection$name))


