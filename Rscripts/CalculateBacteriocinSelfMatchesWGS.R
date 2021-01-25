#Script tabulates matches in own genomes
#

options(warn=2)
rm(list=ls())
library(stringr)
library(spatstat.utils)

findit3<-function(x,y) all(x[1]>=y[1], x[2]<=y[2])


verbose=F

#read in log file (remember quote option due to internal "'" in some names)
names=read.table("c:/Users/milst/Desktop/all_phageBGC/all.log", header = F,sep = "\t", quote="", comment.char = "#", stringsAsFactors = F )
colnames(names)=c("ID","gBGC_hits","prophage_hits","pBGC_hits","pBGC_types","range","genome_len","name")

#Make IDs and species for easier access
names$shortID=gsub("^([A-Z0-9_]*.[0-9])_.*","\\1", names$ID)
names$species=apply(str_split_fixed(names$name, " ", 8)[,2:3],MARGIN = 1, function(x) paste(x, collapse = " "))
names$species=factor(names$species)


#Get BLAST filenames
blast_files=dir("onlyBacteriocins/only_bacteriocin_WGS_blastn/", full.names = T)


n_genomes=length(blast_files)

grep("GCF_000699465.1.WGS.blastn", blast_files )

i=211
i=216

collection=data.frame(name=character(n_genomes), inMatch=numeric(n_genomes),outMatch=numeric(n_genomes),pBGCs=numeric(n_genomes), stringsAsFactors = F)

counter=1
#Main loop
for(i in 1:length(blast_files)) {
  
  #get genome id and find it in log file
  genome_id=gsub("^onlyBacteriocins/only_bacteriocin_WGS_blastn/([A-Z0-9_]*.[0-9]).WGS..*","\\1", blast_files[i])
  genome_indx=grep(genome_id, names$shortID)
  genome_entry=names[genome_indx,]
  pBGC_range_n= str_split_fixed(t(
    str_split_fixed(genome_entry$range,",",3)), "-",2)
  pBGC_ranges=matrix(apply(matrix(pBGC_range_n[!pBGC_range_n[,1]==""],ncol=2), 2, as.numeric),ncol=2)
  
  
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
  
  #pBGC range
  plot(x=c(0,genome_entry$genome_len),y=c(1,1),ylim=c(0,5),  col="white", main=i)
  
  for(k in 1:NROW(pBGC_ranges)) {
    lines(x=pBGC_ranges[k,],y=c(k,k))
  }
  
  collection$name[i]=as.character(genome_entry$species)
  collection$pBGCs=genome_entry$pBGC_hits

  #counter=counter+1 
  
  j=2
  #for hver blast hit, ligger den i phag
  #hvis ikke, summer op
  for(j in 1:NROW(blast)){
    lines(x=c(blast$sstart[j],blast$send[j]) ,y=c(j+.1,j+.1),lwd=3, col="red")
    range_hits=c()
    for(k in 1:NROW(pBGC_ranges)) {
      range_hits[k] = findit3(x=sort(c(blast$sstart[j],blast$send[j])) , y= sort(c(pBGC_ranges[k,1],pBGC_ranges[k,2])))
    }
    
    jMatchCoords    = as.numeric(str_split_fixed(str_split_fixed(blast$qseqid[j], "_", 17)[14],":",2))
    jMatchLen       = diff(jMatchCoords)

    if(!any(range_hits)) {
      cat(paste("i=",i,", blast=",j,", not in pBGC\n", sep=""))
      print(blast$pident[j])
      print( jMatchLen/blast$length[j])
      collection$outMatch[i]=collection$outMatch[i]+1
     
      # collection[counter,1]=as.character(genome_entry$species)
      # collection[counter,2]=length(unq_qseq)
      # collection[counter,3]=length(unq_sseq)
      # collection[counter,4]=genome_entry$pBGC_hits
      # counter=counter+1 
    } else {
      collection$inMatch[i]=collection$inMatch[i]+1
    }
  }
}

collection


table(collection$inMatch)
table(collection$outMatch)

outCollection=collection[which(collection$outMatch>0),]
table(outCollection$name)

png("c:/Users/milst/Desktop/all_phageBGC/Pics/outside_pBGC_match.png", width = 7, height = 5, units = "in", res=600)
par(mar=c(2.5,10,0.5,.5),mgp=c(1.4,0.4,0), family="serif", font.lab=2)
barplot(sort(table(outCollection$name), decreasing = T), horiz = T,las=2, xlab="Count")
dev.off()

