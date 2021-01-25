#Script tabulates bacteriocin core gene matches in all other genomes, including its own


rm(list=ls())
library(stringr)
library(spatstat.utils)
findit3<-function(x,y) all(x[1]>=y[1], x[2]<=y[2])

lenCutoff=0.7


#read in log file (remember quote option due to internal "'" in some names)
names=read.table("c:/Users/milst/Desktop/all_phageBGC/all.log", header = F,sep = "\t", quote="", comment.char = "#", stringsAsFactors = F )
colnames(names)=c("ID","gBGC_hits","prophage_hits","pBGC_hits","pBGC_types","range","genome_len","name")

#Make IDs and species for easier access
names$shortID = gsub("^([A-Z0-9_]*.[0-9])_.*","\\1", names$ID)
names$species = apply(str_split_fixed(names$name, " ", 8)[,2:3],MARGIN = 1, function(x) paste(x, collapse = " "))
names$species = gsub("\\[|\\]","",names$species )
names$genus   = gsub("([A-Za-a]*) .*$","\\1", names$species)

#Read list of blast hits 
blast_hits=dir("c:/Users/milst/Desktop/all_phageBGC/coreGenesVsGenomesBlastn/", full.names = T)

#column annotations
blastVars= c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")

#Data.frame for collecting data
collect=data.frame(
  querySpecies = character(NROW(blast_hits)), #query species
  anyMatches   = numeric(NROW(blast_hits)),   #has any hits
  matchSpecies = numeric(NROW(blast_hits)),   #1 if matches own species only, else 0
  matchGenus   = numeric(NROW(blast_hits)),   #1 if matches own genus only, else 0
  query_pBGC   = numeric(NROW(blast_hits)),   #has pBGC
  subject_pBGC = numeric(NROW(blast_hits)),   #has pBGC
  outside_match = numeric(NROW(blast_hits)),  #matches itself outside phage
  stringsAsFactors = F)


pBGCs=c()

i=1988 
i=52
for(i in 1:length(blast_hits)){
  
  #grab blast query, e.g. the genome
  queryGCF     = gsub("^c:/Users/milst/Desktop/all_phageBGC/coreGenesVsGenomesBlastn/([A-Z0-9_]*.[0-9]).*","\\1", blast_hits[i])
  queryLog     = names[match(queryGCF, names$shortID),]
  querySpecies = queryLog$species
  queryGenus   = queryLog$genus
  query_pBGCs  = names[match(queryGCF, names$shortID),]$pBGC_hits
  query_pBGC_range_n= str_split_fixed(t(
    str_split_fixed(queryLog$range,",",3)), "-",2)
  query_pBGC_ranges=matrix(apply(matrix(query_pBGC_range_n[!query_pBGC_range_n[,1]==""],ncol=2), 2, as.numeric),ncol=2)
  
  #routine for handling empty blast files
  if(file.info(blast_hits[i])$size==0) {
    collect$querySpecies[i] = querySpecies
    collect$anyMatches[i]   = 0
    collect$matchSpecies[i] = 0
    collect$matchGenus[i]   = 0
    collect$query_pBGC[i]   = query_pBGCs
    collect$subject_pBGC[i] = 0
    
    next
  }
  
  #read blast file and and variable names
  blast=read.table(blast_hits[i])
  colnames(blast) = blastVars
  
  #vectors for collecting species and pBGG-counts
  subjsSpecies=c()
  subj_pBGC=c()
  
  j=2
  #loop through each row of blast file, e.g. each subject
  for(j in 1:NROW(blast)) {
    
    jSubjectBlast=blast[j,]
    jSubjectGCF=gsub("^([A-Z0-9_]*.[0-9]).*","\\1",jSubjectBlast$sseqid)
    jSubjectLog     = names[match(jSubjectGCF, names$shortID),]
    jSubjectSpecies = jSubjectLog$species
    jSubjectGenus   = gsub("([A-Za-a]*) .*$","\\1", jSubjectLog$species)
    jSubject_pBGC   = jSubjectLog$pBGC_hits
    jMatchCoords    = as.numeric(str_split_fixed(str_split_fixed(jSubjectBlast$sseqid, "_", 17)[16],":",2))
    jMatchLen       = diff(jMatchCoords)
    
    #if match identity is to low go to next line
    if(jSubjectBlast$length/jMatchLen<lenCutoff) {
      next
    }
    
    #vectors to sum up species and pBGCs
    subjsSpecies[j] = jSubjectSpecies
    subj_pBGC[j]    = jSubject_pBGC
    
    
    #if subject is not the same species as query, tell user
    if(jSubjectSpecies != querySpecies) {
      cat(paste(i,"is a",querySpecies,"but matches",jSubjectSpecies,"ID/len", jSubjectBlast$pident,"/",jSubjectBlast$length,"\n" ))
    }
    
    #Code to handle self match
    if(jSubjectGCF == queryGCF){
      cat(paste(i, "Is itself\n"))
      
      #check if query match and pBGC match is overlapping
      #loops over all potential phage regions, often just one, and stores in vector
      range_hits=c()
      for(k in 1:NROW(query_pBGC_ranges)) {
        range_hits[k] = findit3(x=sort(c(jSubjectBlast$qstart, jSubjectBlast$qend) ) , y= sort(c(query_pBGC_ranges[k,1],query_pBGC_ranges[k,2])))
      }
      
      #if subject is not in any of the phage regions, collect outside match
      if(!any(range_hits)){ 
        cat(paste( "\tNOT same match\n"))
        collect$outside_match[i]=1
      }
      
      #if no self match, store other pBGCs
    } else {
      if(!any(grep(jSubjectBlast$sseqid, pBGCs, fixed = T ))){
        pBGCs=c(pBGCs,as.character(jSubjectBlast$sseqid))
      }
    }
    
  }
  
  
  subjsSpecies=na.exclude(subjsSpecies)
  subj_pBGC=na.exclude(subj_pBGC)
  
  if(length(subjsSpecies)==0) {
    #routine for handling empty blast files
    collect$querySpecies[i] = querySpecies
    collect$anyMatches[i]   = 0
    collect$matchSpecies[i] = 0
    collect$matchGenus[i]   = 0
    collect$query_pBGC[i]   = query_pBGCs
    collect$subject_pBGC[i] = 0
  } else {
    collect$querySpecies[i] = querySpecies
    collect$anyMatches[i ]  = 1
    collect$matchSpecies[i] = all(unique(subjsSpecies)==querySpecies)
    collect$matchGenus[i]   = all(unique( gsub("([A-Za-a]*) .*$","\\1",subjsSpecies))==queryGenus)
    collect$query_pBGC[i]   = query_pBGCs
    collect$subject_pBGC[i] = any(jSubject_pBGC>0)
  }
  
}

collectSelfOutside=subset(collect, outside_match==1 )
outsideGCF=gsub(".blastn","",str_split_fixed(blast_hits[as.numeric(rownames(collectSelfOutside))],"/",7)[,7])
outsideNames=names[match(outsideGCF, names$shortID),]
blastForOutside=blast
blastForOutside[,1:2]=apply(blastForOutside[,1:2],2, as.character)

#outsideNames[k,]
#k=14

for(k in 1:length(outsideGCF)) {
  
  kBlast=read.table(grep(outsideGCF[k], blast_hits, value=T))
  colnames(kBlast) = blastVars
  selfBlast=kBlast[grep(outsideGCF[k], kBlast$sseqid) ,]
  selfBlast[,1:2]=apply(selfBlast[,1:2],2, as.character)
  blastForOutside[k,]=selfBlast[!dplyr::between(selfBlast$qstart, as.numeric(str_split_fixed(outsideNames$range[k], "-",2))[1],as.numeric(str_split_fixed(outsideNames$range[k], "-",2))[2]),]
}



write.table(x =data.frame(outsideNames, blastForOutside) ,file = "GenomesWithMatchesOutsidePhage.csv",sep=";", row.names = F)


#all genomes with a pBGC
collectCarrier   = subset(collect, query_pBGC>0 ) 
#genomes with no pBGC
collectNoCarrier = subset(collect, query_pBGC==0 ) 
#genomes with no pBGC but a match
collectHasMatch  = subset(collectNoCarrier, anyMatches==1 ) 
#genomes with no pBGC but a match which is not same genus
collectNoGenus   = subset(collectHasMatch, matchGenus==0)

table(collect$anyMatches)
table(collectCarrier$anyMatches)
table(collectNoCarrier$anyMatches)
table(collectHasMatch$anyMatches)


barplot(sort(table(gsub("([A-Za-a]*) .*$","\\1",collectHasMatch$querySpecies)), decreasing = T), horiz = T,las=2, xlab="Count")


###
#how many pBGCs are found
length(unique(gsub("(GCF_[0-9]*..)_.*","\\1", pBGCs)))

###
#How many non carriers has a match
table(collectHasMatch$anyMatches)/(length(blast_hits)-NROW(collectCarrier))

##
#how many match the same species
table(collectHasMatch$matchSpecies)/sum(table(collectHasMatch$matchSpecies))

#
#how many match the same genus
table(collectHasMatch$matchGenus)/sum(table(collectHasMatch$matchGenus))

  


png("c:/Users/milst/Desktop/all_phageBGC/Pics/nonGenusCoreGeneMatch.png", width = 7, height = 5, units = "in", res=600)
par(mar=c(2.5,10,0.5,.5),mgp=c(1.4,0.4,0), family="serif", font.lab=2)
barplot(sort(table(collectNoGenus$querySpecies), decreasing = T), horiz = T,las=2, xlab="Count")
dev.off()

