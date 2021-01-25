#####
#Summaries
#####

library(scatterplot3d)
library(RColorBrewer)
library(stringr)

rm(list=ls())


writePDF=F

#read the log file and fix names
all_log=read.table("all.log", sep="\t", header=T, comment.char = "")
names=gsub("#","",stringr::str_split_fixed(readLines("all.log",  n = 1), "\t",6))
colnames(all_log)=c("ID","genome_anti_hits","prophage_hits", "pBGC_hits", "pBGC_types","pBGC_coordinates","genome_length","tax_name")

#make species and genus names
all_log$species=apply(str_split_fixed(all_log$tax_name, " ", 8)[,2:3],MARGIN = 1, function(x) paste(x, collapse = " "))
all_log$genus=(gsub("\\[|\\]","",str_split_fixed(all_log$tax_name, " ", 8)[,2]))

#SUMMARIES
#%-age of genomes with prophage
100*length(which(all_log$prophage_hits>0))/NROW(all_log)

#%-age of genomes without prophage
100-(100*length(which(all_log$prophage_hits>0))/NROW(all_log))

#%-age with pBGC
100*length(which(all_log$pBGC_hits>0))/NROW(all_log)

#%-age of phage positives with pBGC
100*length(which(all_log$pBGC_hits>0))/length(which(all_log$prophage_hits>0))


##########
#SUBSETTING TO pBGC-positive
pBGC_indx=which(all_log$pBGC_hits>0)
only_pBGC=all_log[pBGC_indx,]

#WHAT ARE THE TYPES
sort(table(only_pBGC$pBGC_types))

#IN WHICH ARE pBGCs A MAJOR SOURCE OF BGC
sort(table(only_pBGC[only_pBGC$genome_anti_hits<2,]$species))

#CHECKING LISTERIA
list_mono=subset(all_log, subset = species=="Listeria monocytogenes")
plot(list_mono$pBGC_hits~jitter(list_mono$genome_anti_hits))

#CHECHING ENTEROCOCCUS
ent_fae=subset(all_log, subset = species=="Enterococcus faecalis")
plot(ent_fae$pBGC_hits~jitter(ent_fae$genome_anti_hits))

#CHECKING OTHERS
subset(all_log, subset = pBGC_types=="other")

mann=subset(all_log, subset = genus=="Mannheimia")
sum(mann$prophage_hits)

write.table(mann,file = "mannheimia_table.txt", row.names = F)

##########

#WHICH GENOMES ARE OVER REPRESENTED?
sort(table(all_log$species), decreasing = T)[1:50]

#WHICH pBGCs carriers ARE OVERREPRESENTED
sort(table(all_log$species[pBGC_indx]), decreasing = T)

#WHICH HAVE MULTIPLE pBGCs
table(all_log[which(all_log$pBGC_hits>1),]$species)




topNames=names(sort(table(all_log$genus[which(all_log$pBGC_hits>0)]), decreasing = T)[1:11])

all_log$mainNames=as.character(all_log$genus)
all_log$mainNames[which(is.na(match(all_log$mainNames, topNames)))]="others"
all_log$mainNames=factor(all_log$mainNames, levels=c(topNames,"others"))  

cols=brewer.pal(12, "Paired")
all_log$cols=cols[as.numeric(all_log$mainNames)]
#all_log=all_log[order(all_log$mainNames),]




#gBGCs more abundant in pBGC hosts?
plot(jitter(c(rep(1, 307),rep(2,14877) ) ), (c(all_log$genome_anti_hits[pBGC_indx],all_log$genome_anti_hits[-pBGC_indx])))
t.test(log10(all_log$genome_anti_hits[pBGC_indx]+1),log10(all_log$genome_anti_hits[-pBGC_indx]+1))
wilcox.test(all_log$genome_anti_hits[pBGC_indx],all_log$genome_anti_hits[-pBGC_indx], paired = F)


GLM1=glm(c(all_log$genome_anti_hits[pBGC_indx],all_log$genome_anti_hits[-pBGC_indx]) ~ (c(rep("A", 307),rep("B",14877) ) ))
summary(GLM1)
plot(GLM1)

#prophages more abundant in pBGC hosts?
boxplot(list(all_log$prophage_hits[pBGC_indx],all_log$prophage_hits[-pBGC_indx]))
t.test(all_log$prophage_hits[pBGC_indx],all_log$prophage_hits[-pBGC_indx])
wilcox.test(all_log$prophage_hits[pBGC_indx],all_log$prophage_hits[-pBGC_indx], paired = F)


if(writePDF) pdf("Pics/pBGCvsPHAGEvsgBGC.pdf", onefile = T)

########
#Genomic BGCs as a function of phages
########

png("Pics/pBGCvsPHAGEvsgBGC.png", width=7.4, height = 4, units = "in", res=600)
par(mai=c(.5,0.5,0.1,0.1), mgp=c(1.4,0.4,0), font.lab=2, family="serif")
with(all_log[-pBGC_indx,], 
     plot( jitter((genome_anti_hits),amount = 1) ~ jitter((prophage_hits), amount=1 ) ,cex= 0.5,
           xlab="#Phages", ylab="#gBGCs" , xlim=c(-1,23), 
           pch=21, bg=cols, col=ifelse(pBGC_hits>0,1,"grey50") )
)
with(all_log[pBGC_indx,], 
     points( jitter((genome_anti_hits),amount = 1) ~ jitter((prophage_hits), amount=1 ) , cex=ifelse(pBGC_hits==1,1,1.5),
            pch=ifelse(pBGC_hits==1,24,23), bg=cols, col=1 )
)

legend("topright" , y.intersp = .8, legend = levels(all_log$mainNames), pt.bg=cols, pch=21, title="Genus")
legend("bottomright", y.intersp = 1,   legend = c("0","1",">1"), pt.bg="grey90", pch=c(21,24,23), title = "pBGC")
dev.off()
########
#pBGCs as a function of phages
########

par(mai=c(.5,0.5,0.1,0.1), mgp=c(1.4,0.4,0), font.lab=2, family="serif")
plot( jitter((all_log$pBGC_hits),amount = .1) ~ jitter(((all_log$prophage_hits)),amount = .2), cex=ifelse(all_log$pBGC_hits>0,1,0.5),
      xlab="Phages", ylab="#pBGCs" , pch=ifelse(all_log$pBGC_hits>0,24,21), bg=all_log$cols)
legend("topright" , legend = levels(all_log$mainNames), pt.bg=cols, pch=21)

########
#Genomic BGCs as a function of phages
########
par(mai=c(.5,0.5,0.1,0.1), mgp=c(1.4,0.4,0), font.lab=2, family="serif")
plot( jitter((all_log$pBGC_hits),amount = .1) ~ jitter(((all_log$genome_anti_hits)),amount = .1), cex=ifelse(all_log$pBGC_hits>0,1,0.5),
      xlab="#gBGCs", ylab="#pBGCs" , pch=ifelse(all_log$pBGC_hits>0,24,21), bg=all_log$cols)
legend("topright" , legend = levels(all_log$mainNames), pt.bg=cols, pch=21)


source('http://www.sthda.com/sthda/RDoc/functions/addgrids3d.r')
s3d=scatterplot3d(x = log1p(jitter(all_log$genome_anti_hits)), y=log1p(jitter(all_log$prophage_hits)),
                  z=all_log$pBGC_hits,bg=all_log$cols, pch=" ",grid=TRUE, box=FALSE, xlab = "gBGCs", ylab = "Phages",zlab = "pBGCs")

addgrids3d(x = log1p(jitter(all_log$genome_anti_hits)), y=log1p(jitter(all_log$prophage_hits)),
           z=all_log$pBGC_hits, grid = c( "xz", "yz"))
s3d$points3d(x = log1p(jitter(all_log$genome_anti_hits)), y=log1p(jitter(all_log$prophage_hits)),
             z=all_log$pBGC_hits,bg=all_log$cols, pch=ifelse(all_log$pBGC_hits>0,24,21))
if(writePDF) dev.off()