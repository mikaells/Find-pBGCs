library("ggExtra")
library("ggplot2")
library(gridExtra)
library(stringr)
library(seqinr)
library(RColorBrewer)


rm(list=ls())

bacti=read.fasta("BACTIBASE_renamed.faa")


#read in log file (remember quote option due to internal "'" in some names)
names=read.table("c:/Users/milst/Desktop/all_phageBGC/all.log", header = F,sep = "\t", quote="", comment.char = "#", stringsAsFactors = F )
colnames(names)=c("ID","gBGC_hits","prophage_hits","pBGC_hits","pBGC_types","range","genome_len","name")

#Make IDs and species for easier access
names$shortID=gsub("^([A-Z0-9_]*.[0-9])_.*","\\1", names$ID)
names$species=apply(str_split_fixed(names$name, " ", 8)[,2:3],MARGIN = 1, function(x) paste(x, collapse = " "))
names$species=factor(names$species)


#blastHits=read.table("onlyBacteriocins/bacteriocins_vs_bactibase.blastp")
blastHits=read.table("onlyBacteriocins/bacteriocinVSbactibase_no_limit_hsps1.blastp", sep="\t")
colnames(blastHits) = c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")

#finding subject length by lookup in BACTIBASE
blastHits$slen=c()
for(i in 1:NROW(blastHits)) {
  bactiIndx=grep(gsub("(BAC[0-9]*)_.*","\\1",blastHits$sseqid[i]), names(bacti))
  blastHits$slen[i]=length(bacti[[bactiIndx]])
}
#working out the length of the query by fishing out of name and dividing by 3
blastHits$qlen=apply(str_split_fixed(str_split_fixed(blastHits$qseqid, "_", 19)[,18],":",2),MARGIN = 1, FUN = function(x) {diff(as.numeric(x))})/3
#percentage of length match by division of match
blastHits$pLen=100*blastHits$length/blastHits$slen
#working out match type from BACTIBASE hit
blastHits$type=gsub("BAC[0-9]*_\\|_([A-Za-z]*).+","\\1",blastHits$sseqid)
##getting genome name
blastHits$genome=gsub("^([A-Z0-9_]*.[0-9]).*","\\1", blastHits$qseqid)
#getting species
blastHits$species= names$species[match(blastHits$genome, names$shortID)]
#getting genus
blastHits$genus=gsub("([A-Za-a]*) .*$","\\1", blastHits$species)




#capitalizing BACTIBASE names for niceness
blastHits$type=Hmisc::capitalize(as.character(blastHits$type))

#finding dominating bacteriocins for clor scheme
topNames=names(sort(table(blastHits$type), decreasing = T)[1:11])

#set up mainNames, is genus and small as others
blastHits$mainNames=as.character(blastHits$type)
blastHits$mainNames[which(is.na(match(blastHits$mainNames, topNames)))]="others"
blastHits$mainNames=factor(blastHits$mainNames, levels=c(topNames,"others")) 

#set up color scheme according to fig 2
cols=c(brewer.pal(12, "Paired"))

#make colors
blastHits$cols=cols[as.numeric(blastHits$mainNames)]



png("Pics/Fig3_BACTIBASE_hits.png", width = 3.5, height =3.5, units = "in", res=300 )
#scatterplot
p1=ggplot(blastHits, aes(x=jitter(pident, amount = 1))) + geom_point(aes(y=jitter(pLen, amount=1),colour=factor(blastHits$mainNames)),size=1) + 
  theme_bw() + theme(legend.text  = element_text(size=5),axis.text=element_text(size=8),axis.title=element_text(size=8,face="bold"), legend.position  = c(0.14, 0.85),legend.key.size = unit(-1, 'lines'), legend.background=element_blank())  +
  labs(x = "% identity", y="% length")  + xlim(0, 105) +
  scale_color_manual(name = "", values = cols)
#histograms
ggMarginal(p1, type = "histogram")
dev.off()


write.table(x = blastHits,file = "bacteriocinBlastSummary.csv", sep=";", row.names = F)
subset(blastHits, pident<75 & pLen>90 )


#How many genomes have a match
length(unique(blastHits$genome))


#which has no match in bactibase
subset(names, subset =  shortID==subset(names, subset = pBGC_hits>0)$shortID[is.na(match(subset(names, subset = pBGC_hits>0)$shortID,unique(blastHits$genome)))])

#get linocin hits
blastHits[grep("Linocin",blastHits$type),]

highs=blastHits[blastHits$pident>80,]

table(highs$type)

subset(blastHits, subset = genus=="Pseudomonas")


