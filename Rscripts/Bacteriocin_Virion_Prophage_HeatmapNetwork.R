########
#A script to combine Free Phages and Prophages ANI in the same heatmap
########

#Data is average nucleotide identity of prophages with pBGC and PATRIC phages with pBGC
#Note, entire phage regions are used

#Libraries
library(pheatmap)
library(stringr)
library(RColorBrewer)
library(igraph)
source("Scripts/functions.R")
####################
rm(list=ls())



set.seed(10)

#Read meta-data and add ID and species columns 
names=read.table("c:/Users/milst/Desktop/all_phageBGC/pBGC_hit.list", sep="\t")
names$shortID=gsub("^([A-Z0-9_]*.[0-9])_.*","\\1", names$V1)
names$species=apply(str_split_fixed(names$V8, " ", 8)[,2:3],MARGIN = 1, function(x) paste(x, collapse = " "))

free=read.table("PATRIC_genome.txt", sep="\t", header=T)


###########

#files of ani output
ANI_files=dir("c:/Users/milst/Desktop/all_phageBGC/onlyBacteriocins/ani_virion_prophage_bacteriocin/ani_vir_pro_bacteriocin/", full.names = T, pattern = "tab")

#Read a  ani in 
bactoANI=read.table(ANI_files[4])
bactoANI[which(bactoANI>1, arr.ind = T)]=1



#make translation table for meta to ani
#Problem is that there are more genomes than entries in meta, since some have more ppBGCs
translationTable=data.frame(clusterID=colnames(bactoANI))
translationTable$genome=gsub("^([A-Z0-9_]*.[0-9]).*","\\1", translationTable$clusterID)
translationTable$species=rep("", NROW(translationTable))

#Go through all genomes and chech if genome has a match in meta
#if not, add species as a Virion
#if yes, update species according to meta
i=152
for(i in 1:NROW(translationTable)){
  if(length(which(translationTable$genome[i]== names$shortID))==0) {
    #free[grep(translationTable$genome[i], free$GenBank.Accessions),"Genome.Name"]
    phageName=gsub("([A-Za-z0-9*]) .*","\\1",free[grep(translationTable$genome[i], free$GenBank.Accessions),"Genome.Name"])
    translationTable$species[i] = paste("Virion", phageName, sep=" ")
    
  } else {
    translationTable$species[i]=names$species[which(translationTable$genome[i]== names$shortID)]
  }
}

#make genus column
translationTable$genus=gsub("([A-Za-a]*) .*$","\\1", translationTable$species)

############
#make color schemes
#have to through complicated method for pheatmap

#Biggest 12 genera
topNames=names(sort(table(translationTable$genus), decreasing = T)[1:12])
#flipping around so that virion is last
topNames=c(topNames[-3],"Virion")

#set up mainNames, is genus and small as others
translationTable$mainNames=as.character(translationTable$genus)
translationTable$mainNames[which(is.na(match(translationTable$mainNames, topNames)))]="others"
translationTable$mainNames=factor(translationTable$mainNames, levels=c(topNames,"others")) 

#set up color scheme according to fig 2
#flip around a bit
cols=c(brewer.pal(12, "Paired"),"grey70")
cols=cols[c(1:11,13,12)]

#make colors
translationTable$cols=cols[as.numeric(translationTable$mainNames)]
#all_log=all_log[order(all_log$mainNames),]

genusCols=cols
#add names to colors
names(genusCols)=c(topNames,"others")

#make into a list for pheatmap
mycolors=list(Genus=genusCols)

#make Virions grey
#mycolors$Genus[which(names(mycolors$Genus)=="Virion")]="grey80"

#make annotation df for labels on columns
annotaDF=data.frame(Genus=translationTable$mainNames, Species=translationTable$species)
rownames(annotaDF)=translationTable$clusterID






# png("c:/Users/milst/Desktop/all_phageBGC/Pics/Virion_vs_Prophage_bacteriocin_heatmap.png", width =  11.69 , height = 8.27, units = "in", res=600)
# pH=pheatmap(bactoANI,labels_row = translationTable$genome,labels_col = translationTable$species,annotation_col = subset(annotaDF, select = "Genus"),  annotation_row = subset(annotaDF, select = "Genus"),annotation_colors = mycolors,
#          fontsize = 10, clustering_method = "ward.D2", fontsize_row = 4, fontsize_col = 4)
# dev.off()

png("c:/Users/milst/Desktop/all_phageBGC/Pics/Virion_vs_Prophage_bacteriocin_heatmap_small.png", width =  6 , height = 6, units = "in", res=600)
pH=pheatmap(bactoANI, legend = F,labels_row = "",labels_col = "",annotation_col = subset(annotaDF, select = "Genus"),  annotation_row = subset(annotaDF, select = "Genus"),annotation_colors = mycolors,
         fontsize = 10, clustering_method = "ward.D2", fontsize_row = 4, fontsize_col = 4)
dev.off()

bacto_tree=pH$tree_row

i=4
# pdf("c:/Users/milst/Desktop/all_phageBGC/Pics/Virion_vs_Prophage_bacteriocin_heatmap_zooms_20.pdf",onefile = T)
# for( i in 1:20) {
  cut=cutree(pH$tree_row,k = 20)
  cutNum=i#which.max(table(cut))
  clusterIndx=which(cut==cutNum)

  if(sum(abs(diff(unlist(bactoANI[clusterIndx,clusterIndx]))))==0) next
  pheatmap(bactoANI[clusterIndx,clusterIndx],labels_row = translationTable$genome[clusterIndx],labels_col = translationTable$species[clusterIndx],annotation_col = (subset(annotaDF[clusterIndx,], select = "Genus")),  annotation_row = (subset(annotaDF[clusterIndx,], select = "Genus")),annotation_colors = mycolors,
           fontsize = 10, clustering_method = "ward.D2", fontsize_row = 4, fontsize_col = 4)
# }
# dev.off()

table(translationTable$species[clusterIndx])
  
translationTable2=data.frame(translationTable, col=mycolors$Genus[(translationTable$genus)], Genome= gsub("^([A-Z0-9_]*.[0-9])_.*","\\1", translationTable$clusterID))

translationTable2$nodeNames=as.character(translationTable2$genus)


#translationTable2$nodeNames[which(translationTable2$genus==j)]=paste(translationTable2$nodeNames,1:393)[which(translationTable2$genus==j)]

j="Mannheimia"
for(j in names(which((table(translationTable2$genus))>1))) {
  if(j=="Virion"){ next }
  if(j=="Mannheimia") {
    translationTable2$nodeNames[which(translationTable2$genus==j)[-c(1,71)]]=""
  } else {
    translationTable2$nodeNames[which(translationTable2$genus==j)[-c(1)]]=""
  }
}

translationTable2$nodeNames[which(translationTable2$genus=="Virion")]=gsub("(.*) (.*)","\\2 \\1",as.character(translationTable2$species[which(translationTable2$genus=="Virion")]))
sort(table(grep("Virion",translationTable2$nodeNames, value = T)), decreasing = T)

translationTable2$nodeNames[which(translationTable2$nodeNames=="Escherichia Virion")[-1]]=""
translationTable2$nodeNames[which(translationTable2$nodeNames=="Stx2-converting Virion")[-1]]=""
translationTable2$nodeNames[which(translationTable2$nodeNames=="Enterobacteria Virion")[-1]]=""
translationTable2$nodeNames[which(translationTable2$nodeNames=="Mannheimia Virion")[-1]]=""
translationTable2$nodeNames[which(translationTable2$nodeNames=="Shigella Virion")[-1]]=""
translationTable2$nodeNames[which(translationTable2$nodeNames=="Stx Virion")[-1]]=""
translationTable2$nodeNames[which(translationTable2$nodeNames=="Mycobacterium Virion")[-1]]=""
translationTable2$nodeNames[which(translationTable2$nodeNames=="Listeria Virion")[-1]]=""




# Make an Igraph object from ANI:
network <- graph_from_adjacency_matrix(as.matrix(bactoANI) , weighted=T, mode="undirected", diag=F )

l=layout_nicely(network) 



#png("c:/Users/milst/Desktop/all_phageBGC/Pics/Virion_vs_Prophage_bacteriocin_network.png", width =  8.27 , height = 8.27, units = "in", res=600)
pdf(file = "c:/Users/milst/Desktop/all_phageBGC/Pics/Virion_vs_Prophage_bacteriocin_network.pdf", width =  8.27 , height = 8.27 )
par(mar=c(2.2,2.2,0.2,0.7))


# Basic chart
allPlot=plot(x = network,  axes=F,rescale=T,layout=l,
             # === vertex
             vertex.color = ifelse(translationTable2$genus=="Virion",translationTable2$cols,adjustcolor(as.character(translationTable2$cols), alpha.f = .5)),   # Node color
             vertex.frame.color = NA,                # Node border color
             vertex.shape=ifelse(translationTable2$genus=="Virion","sphere", "circle"),                        # One of “none”, “circle”, “square”, “csquare”, “rectangle” “crectangle”, “vrectangle”, “pie”, “raster”, or “sphere”
             vertex.size=10,                               # Size of the node (default is 15)
             vertex.size2=NA,                               # The second size of the node (e.g. for a rectangle)
             
             # === vertex label
             vertex.label=translationTable2$nodeNames,                 # Character vector used to label the nodes
             vertex.label.color="black",
             vertex.label.font=2,                          # Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
             vertex.label.cex=.8,                          # Font size (multiplication factor, device-dependent)
             vertex.label.dist=0,                          # Distance between the label and the vertex
             vertex.label.degree=0 ,                       # The position of the label in relation to the vertex (use pi)
             
             edge.color=color.gradient(E(network)$weight),#"grey50",                           # Edge color
             edge.width=E(network)$weight,#edge.betweenness(network)*0.01,                                 # Edge width, defaults to 1
             edge.arrow.size=1,                            # Arrow size, defaults to 1
             edge.arrow.width=1,                           # Arrow width, defaults to 1
             edge.lty="solid",                             # Line type, could be 0 or “blank”, 1 or “solid”, 2 or “dashed”, 3 or “dotted”, 4 or “dotdash”, 5 or “longdash”, 6 or “twodash”
             edge.curved=0.3      ,                        # Edge curvature, range 0-1 (FALSE sets it to 0, TRUE to 0.5)
             #axes=T
             #xlim=c(min(l),max(l)),ylim=c(min(l),max(l))#, asp = 0
)
#legend("bottomleft", legend = c((topNames),"others"), pch=16, col= mycolors$Genus, y.intersp = 0.65, bty = "n")

dev.off()

edgelist=read.table("C:/Users/milst/Documents/HCP Anywhere/Curr Biol rebuttal/bigscape_out/network_files/2021-01-20_12-52-07_hybrids_glocal/RiPPs/RiPPs_c0.30.network", header=T, sep="\t")
gg <- graph.data.frame(edgelist, directed=FALSE)


l2=layout_nicely(gg)

V(gg)
V(network)
plot(l2)
plot(x = gg,  axes=F,rescale=T,layout=l,
     # === vertex
     vertex.color = ifelse(translationTable2$genus=="Virion",translationTable2$cols,adjustcolor(as.character(translationTable2$cols), alpha.f = .5)),   # Node color
     vertex.frame.color = NA,                # Node border color
     vertex.shape=ifelse(translationTable2$genus=="Virion","sphere", "circle"),                        # One of “none”, “circle”, “square”, “csquare”, “rectangle” “crectangle”, “vrectangle”, “pie”, “raster”, or “sphere”
     vertex.size=10,                               # Size of the node (default is 15)
     vertex.size2=NA,                               # The second size of the node (e.g. for a rectangle)
     
     # === vertex label
     vertex.label=translationTable2$nodeNames,                 # Character vector used to label the nodes
     vertex.label.color="black",
     vertex.label.font=2,                          # Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
     vertex.label.cex=.8,                          # Font size (multiplication factor, device-dependent)
     vertex.label.dist=0,                          # Distance between the label and the vertex
     vertex.label.degree=0 ,                       # The position of the label in relation to the vertex (use pi)
     
     edge.color=color.gradient(E(network)$weight),#"grey50",                           # Edge color
     edge.width=E(network)$weight,#edge.betweenness(network)*0.01,                                 # Edge width, defaults to 1
     edge.arrow.size=1,                            # Arrow size, defaults to 1
     edge.arrow.width=1,                           # Arrow width, defaults to 1
     edge.lty="solid",                             # Line type, could be 0 or “blank”, 1 or “solid”, 2 or “dashed”, 3 or “dotted”, 4 or “dotdash”, 5 or “longdash”, 6 or “twodash”
     edge.curved=0.3      ,                        # Edge curvature, range 0-1 (FALSE sets it to 0, TRUE to 0.5)
     #axes=T
     #xlim=c(min(l),max(l)),ylim=c(min(l),max(l))#, asp = 0
)
