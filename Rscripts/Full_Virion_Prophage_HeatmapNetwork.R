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


####################
rm(list=ls())

source("Scripts/functions.R")
color.gradient <- function(x, colors=c("grey90","red"), colsteps=100) {
  return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
}

plot(1,1, col="grey50")
set.seed(10)

#Read meta-data and add ID and species columns 
names=read.table("c:/Users/milst/Desktop/all_phageBGC/pBGC_hit.list", sep="\t")
names$shortID=gsub("^([A-Z0-9_]*.[0-9])_.*","\\1", names$V1)
names$species=apply(str_split_fixed(names$V8, " ", 8)[,2:3],MARGIN = 1, function(x) paste(x, collapse = " "))

free=read.table("PATRIC_genome.txt", sep="\t", header=T)


bigscape=read.table("c:/users/milst/Desktop/all_phageBGC/BiGSCAPE_clusts.tsv", header = T)



###########

#files of ani output
ANI_files=dir("c:/Users/milst/Desktop/all_phageBGC/bacteria_pBCGs_and_phage_pBGCs/ani_out/", full.names = T, pattern = "tab")

#Read a  ani in 
phageANI=read.table(ANI_files[4])
phageANI[which(phageANI>1, arr.ind = T)]=1

#make translation table for meta to ani
#Problem is that there are more genomes than entries in meta, since some have more ppBGCs
translationTable=data.frame(clusterID=colnames(phageANI))
translationTable$genome=gsub("^([A-Z0-9_]*.[0-9])_.*","\\1", translationTable$clusterID)
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

translationTable$host=translationTable$genus

translationTable[grep("Virion",translationTable$species),]$host=gsub("Virion ","",translationTable[grep("Virion",translationTable$species),]$species)
translationTable[grep("Stx",translationTable$host),]$host="Escherichia"

factor(translationTable$host)

translationTable$bigscape1=""
translationTable$bigscape2=""
i=1
#add bigscape annotations
for(i in 1:NROW(translationTable)){
  
  bigIndx=grep(translationTable$genome[i],bigscape$ACC )
  
  if(length(bigIndx)==0){
    translationTable$bigscape1[i]="unknown"
  } else {
    print(colnames(bigscape)[which(bigscape[bigIndx,]>0)[-1]])
    translationTable$bigscape1[i]=colnames(bigscape)[which(bigscape[bigIndx,]>0)[2]]
    if(length(colnames(bigscape)[which(bigscape[bigIndx,]>0)])>2){
      translationTable$bigscape2[i]=colnames(bigscape)[which(bigscape[bigIndx,]>0)[3]]
    }
  }
  
}



table(subset(translationTable, subset = genus=="Escherichia")$bigscape1)
table(subset(translationTable, subset = genus=="Escherichia")$bigscape2)


############
#make color schemes
#have to through complicated method for pheatmap

#Biggest 12 genera
topNames=names(sort(table(translationTable$genus), decreasing = T)[1:12])
#flipping around so that virion is last
topNames=c(topNames[-3])#,"Virion")

#set up mainNames, is genus and small as others
translationTable$mainNames=as.character(translationTable$host)
translationTable$mainNames[which(is.na(match(translationTable$mainNames, topNames)))]="others"
translationTable$mainNames=factor(translationTable$mainNames, levels=c(topNames,"others")) 

#set up color scheme according to fig 2
#flip around a bit
cols=c(brewer.pal(12, "Paired"))#,"grey70")
#cols=cols[c(1:11,13,12)]

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


##########Bigscape cols

#set up color scheme 
#flip around a bit
pal=c("purple","red","yellow","springgreen","royalblue","grey90")
pal=rev(brewer.pal(n = 12, name = 'Paired'))
bg_colfunc<-colorRampPalette(pal)
bg_cols=(bg_colfunc(length(unique(translationTable$bigscape1))))
brewer.pal(n = 8, name = 'Dark2')
#make colors
translationTable$bgCols=bg_cols[as.numeric(factor(translationTable$bigscape1))]
#all_log=all_log[order(all_log$mainNames),]

translationTable[translationTable$bigscape1=="unknown",]$bgCols="black"




#png("c:/Users/milst/Desktop/all_phageBGC/Pics/Virion_vs_Prophage_full_heatmap.png", width =  11.69 , height = 8.27, units = "in", res=600)
png("c:/Users/milst/Desktop/all_phageBGC/Pics/Virion_vs_Prophage_full_heatmap_small.png", width =  6 , height = 6, units = "in", res=600)
pH_vir=pheatmap(phageANI,labels_row = "",labels_col = "",annotation_col = subset(annotaDF, select = "Genus"),  annotation_row = subset(annotaDF, select = "Genus"),annotation_colors = mycolors,
                fontsize = 10, clustering_method = "ward.D2", fontsize_row = 4, fontsize_col = 4, legend=F)
dev.off()
dev.off()




translationTable2=data.frame(translationTable, col=mycolors$Genus[(translationTable$genus)], Genome= gsub("^([A-Z0-9_]*.[0-9])_.*","\\1", translationTable$clusterID))
translationTable2$nodeNames=as.character(translationTable2$genus)


#paste(translationTable2$nodeNames,1:393)[which(translationTable2$genus==j)]

j="Pasteurella"
for(j in names(which((table(translationTable2$genus))>1))) {
  if(j=="Virion"){ next }
  if(j=="Pasteurella") {
    translationTable2$nodeNames[which(translationTable2$genus==j)[-6]]=""
  } else if (j=="Bacillus") {
    translationTable2$nodeNames[which(translationTable2$genus==j)[-c(1,8,18,20)]]=""
  } else if (j=="Lactobacillus") {
    translationTable2$nodeNames[which(translationTable2$genus==j)[-c(1,10,13)]]=""
  } else {
    translationTable2$nodeNames[which(translationTable2$genus==j)[-1]]=""
  }
}

translationTable2$nodeNames[which(translationTable2$genus=="Virion")]=gsub("(.*) (.*)","\\2 \\1",as.character(translationTable2$species[which(translationTable2$genus=="Virion")]))

translationTable2$nodeNames[which(translationTable2$nodeNames=="Escherichia Virion")[-1]]=""
translationTable2$nodeNames[which(translationTable2$nodeNames=="Stx2-converting Virion")[-1]]=""
translationTable2$nodeNames[which(translationTable2$nodeNames=="Enterobacteria Virion")[-1]]=""
translationTable2$nodeNames[which(translationTable2$nodeNames=="Mannheimia Virion")[-1]]=""
translationTable2$nodeNames[which(translationTable2$nodeNames=="Shigella Virion")[-1]]=""
translationTable2$nodeNames[which(translationTable2$nodeNames=="Stx Virion")[-1]]=""
translationTable2$nodeNames[which(translationTable2$nodeNames=="Mycobacterium Virion")[-1]]=""
translationTable2$nodeNames[which(translationTable2$nodeNames=="Listeria Virion")[-1]]=""




# Make an Igraph object from ANI:
network <- graph_from_adjacency_matrix(as.matrix(phageANI) , weighted=T, mode="undirected", diag=F )

###ROutines to make exported object
library("NetPathMiner")
igraph::write.graph(network,file = "phageANI_network.csv")

# Create a dataframe nodes: 1st column - node ID, 2nd column -node name
nodes_df <- data.frame(ID = c(1:vcount(network)), NAME = V(network)$name)
# Create a dataframe edges: 1st column - source node ID, 2nd column -target node ID
edges_df <- as.data.frame(get.edges(network, c(1:ecount(network))))

write.gexf(nodes = nodes_df, edges = edges_df,output="phageANI_network.gexf")

write.csv(x = E(network)$weight,file = "phageANI_edges.csv")

##############

#Create layout
l=layout_nicely(network) 
#write.table(x = l, "goodlayout.txt", row.names = F,col.names = F)
l=as.matrix(read.table("goodlayout.txt"))

l2=spreadClus(L = l,transTab = translationTable2,doKleb = T, minSize = 12)


#png("c:/Users/milst/Desktop/all_phageBGC/Pics/Virion_vs_Prophage_full_network.png", width =  8.27 , height = 8.27, units = "in", res=600)


pdf(file = "c:/Users/milst/Desktop/all_phageBGC/Pics/Virion_vs_Prophage_full_network.pdf", width =  8.27 , height = 8.27 )
par(mar=c(2.2,2.2,0.2,0.7))

# Basic chart
allPlot=plot(x = network,  
             axes=F,rescale=T,layout=l2, 
             # === vertex
             vertex.color = ifelse(translationTable2$genus=="Virion",NA,adjustcolor(as.character(translationTable2$cols), alpha.f = .8)),   # Node color
             vertex.frame.color = ifelse(translationTable2$genus=="Virion",NA, NA),# Node border color
             vertex.shape=ifelse(translationTable2$genus=="Virion","circle", "circle"),                        # One of “none”, “circle”, “square”, “csquare”, “rectangle” “crectangle”, “vrectangle”, “pie”, “raster”, or “sphere”
             vertex.size=5,                               # Size of the node (default is 15)
             vertex.size2=NA,                               # The second size of the node (e.g. for a rectangle)
             
             # === vertex label
             vertex.label="",#translationTable2$nodeNames,                 # Character vector used to label the nodes
             vertex.label.color="black",
             vertex.label.font=2,                          # Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
             vertex.label.cex=.8,                         # Font size (multiplication factor, device-dependent)
             vertex.label.dist=0,                          # Distance between the label and the vertex
             vertex.label.degree=0 ,                      # The position of the label in relation to the vertex (use pi)
             
             edge.color=color.gradient(log(E(network)$weight+1),colors = c("grey90","yellow","red")),#"grey50",                           # Edge color
             edge.width=E(network)$weight,#edge.betweenness(network)*0.01,                                 # Edge width, defaults to 1
             edge.arrow.size=1,                            # Arrow size, defaults to 1
             edge.arrow.width=1,                           # Arrow width, defaults to 1
             edge.lty="solid",                             # Line type, could be 0 or “blank”, 1 or “solid”, 2 or “dashed”, 3 or “dotted”, 4 or “dotdash”, 5 or “longdash”, 6 or “twodash”
             edge.curved=0.3      ,                        # Edge curvature, range 0-1 (FALSE sets it to 0, TRUE to 0.5)
             #axes=T
             #xlim=c(min(l),max(l)),ylim=c(min(l),max(l))#, asp = 0
)

# Basic chart
plot(x = network,  add=T,
     axes=F,rescale=T,layout=l2, 
     # === vertex
     vertex.color = ifelse(translationTable2$genus=="Virion",translationTable2$cols,NA),   # Node color
     vertex.frame.color = ifelse(translationTable2$genus=="Virion","black", NA),# Node border color
     vertex.shape=ifelse(translationTable2$genus=="Virion","star", "circle"),                        # One of “none”, “circle”, “square”, “csquare”, “rectangle” “crectangle”, “vrectangle”, “pie”, “raster”, or “sphere”
     vertex.size=5,                               # Size of the node (default is 15)
     vertex.size2=NA,                               # The second size of the node (e.g. for a rectangle)
     
     # === vertex label
     vertex.label=ifelse(translationTable2$mainNames=="others",gsub(" Virion","",translationTable2$nodeNames),"" ),                 # Character vector used to label the nodes
     vertex.label.color="black",
     vertex.label.font=2,                          # Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
     vertex.label.cex=.8,                         # Font size (multiplication factor, device-dependent)
     vertex.label.dist=0,                          # Distance between the label and the vertex
     vertex.label.degree=0 ,                      # The position of the label in relation to the vertex (use pi)
     
     edge.color="NA",#color.gradient(log(E(network)$weight+1),colors = c("grey90","yellow","red")),#"grey50",                           # Edge color
     edge.width=E(network)$weight,#edge.betweenness(network)*0.01,                                 # Edge width, defaults to 1
     edge.arrow.size=1,                            # Arrow size, defaults to 1
     edge.arrow.width=1,                           # Arrow width, defaults to 1
     edge.lty="solid",                             # Line type, could be 0 or “blank”, 1 or “solid”, 2 or “dashed”, 3 or “dotted”, 4 or “dotdash”, 5 or “longdash”, 6 or “twodash”
     edge.curved=0.3      ,                        # Edge curvature, range 0-1 (FALSE sets it to 0, TRUE to 0.5)

)
legend("topleft", legend = c((topNames),"others"),horiz = F, pch=16, col= mycolors$Genus, 
       y.intersp = 0.65, bty = "n", title = "Genus", title.adj = .1)

#plot(scale(l2))
#text(scale(l2),labels = ifelse(translationTable2$mainNames=="others",gsub(" Virion","",translationTable2$nodeNames),"" ))

points(c(-1.06,-1.06),c(.249,.249),pch=c(24,25),bg="grey10")
points(-1.06,0.199,pch=16,bg="grey10")

text(-1.105,0.3,"Organism", pos=4)
text(-1.05,.245, "Virion",    pos = 4)
text(-1.05,.195, "Bacteria", pos = 4)

lgd_ = rep(NA, 9)
lgd_[c(1,9)] = c("1.00",signif(min(E(network)$weight),2))
legend("bottomleft",
       legend = lgd_,
       pch=15,
       col =  colorRampPalette(colors = c("red","yellow","grey90"))(9),
       y.intersp = 0.4,
       #title = "Similarity",
       bty="n",
       cex = 1, text.font = 1)
text(-1.105,-0.85,"Similarity", pos=4)

dev.off()

######

pdf(file = "c:/Users/milst/Desktop/all_phageBGC/Pics/Virion_vs_Prophage_full_network_bigscape.pdf", width =  8.27 , height = 8.27 )
par(mar=c(2.2,2.2,0.2,0.7))

# Basic chart
allPlot=plot(x = network,  
             axes=F,rescale=T,layout=l2, 
             # === vertex
             vertex.color = ifelse(translationTable2$genus=="Virion",NA,adjustcolor(as.character(translationTable2$bgCols), alpha.f = .8)),   # Node color
             vertex.frame.color = ifelse(translationTable2$genus=="Virion",NA, NA),# Node border color
             vertex.shape=ifelse(translationTable2$genus=="Virion","circle", "circle"),                        # One of “none”, “circle”, “square”, “csquare”, “rectangle” “crectangle”, “vrectangle”, “pie”, “raster”, or “sphere”
             vertex.size=5,                               # Size of the node (default is 15)
             vertex.size2=NA,                               # The second size of the node (e.g. for a rectangle)
             
             # === vertex label
             vertex.label="",#translationTable2$nodeNames,                 # Character vector used to label the nodes
             vertex.label.color="black",
             vertex.label.font=2,                          # Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
             vertex.label.cex=.8,                         # Font size (multiplication factor, device-dependent)
             vertex.label.dist=0,                          # Distance between the label and the vertex
             vertex.label.degree=0 ,                      # The position of the label in relation to the vertex (use pi)
             
             edge.color=color.gradient(log(E(network)$weight+1),colors = c("grey90","yellow","red")),#"grey50",                           # Edge color
             edge.width=E(network)$weight,#edge.betweenness(network)*0.01,                                 # Edge width, defaults to 1
             edge.arrow.size=1,                            # Arrow size, defaults to 1
             edge.arrow.width=1,                           # Arrow width, defaults to 1
             edge.lty="solid",                             # Line type, could be 0 or “blank”, 1 or “solid”, 2 or “dashed”, 3 or “dotted”, 4 or “dotdash”, 5 or “longdash”, 6 or “twodash”
             edge.curved=0.3      ,                        # Edge curvature, range 0-1 (FALSE sets it to 0, TRUE to 0.5)
)

# Basic chart
plot(x = network,  add=T,
     axes=F,rescale=T,layout=l2, 
     # === vertex
     vertex.color = ifelse(translationTable2$genus=="Virion",translationTable2$bgCols,NA),   # Node color
     #vertex.color = translationTable2$bgCols,
     #vertex.frame.color = NA,
     vertex.frame.color = ifelse(translationTable2$genus=="Virion","black", NA),# Node border color
     #vertex.frame.color = ifelse(translationTable2$genus=="Virion",1,0 ),
     vertex.shape=ifelse(translationTable2$genus=="Virion","star", "circle"),                        # One of “none”, “circle”, “square”, “csquare”, “rectangle” “crectangle”, “vrectangle”, “pie”, “raster”, or “sphere”
     vertex.size=5,                               # Size of the node (default is 15)
     vertex.size2=NA,                               # The second size of the node (e.g. for a rectangle)
     
     #asp=1.55,
     # === vertex label
     vertex.label=NA,#ifelse(translationTable2$mainNames=="others",gsub(" Virion","",translationTable2$nodeNames),"" ),                 # Character vector used to label the nodes
     vertex.label.color="black",
     vertex.label.font=2,                          # Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
     vertex.label.cex=.8,                         # Font size (multiplication factor, device-dependent)
     vertex.label.dist=0,                          # Distance between the label and the vertex
     vertex.label.degree=0 ,                      # The position of the label in relation to the vertex (use pi)
     
     edge.color="NA",#color.gradient(log(E(network)$weight+1),colors = c("grey90","yellow","red")),#"grey50",                           # Edge color
     edge.width=E(network)$weight,#edge.betweenness(network)*0.01,                                 # Edge width, defaults to 1
     edge.arrow.size=1,                            # Arrow size, defaults to 1
     edge.arrow.width=1,                           # Arrow width, defaults to 1
     edge.lty="solid",                             # Line type, could be 0 or “blank”, 1 or “solid”, 2 or “dashed”, 3 or “dotted”, 4 or “dotdash”, 5 or “longdash”, 6 or “twodash”
     edge.curved=0.3      ,                        # Edge curvature, range 0-1 (FALSE sets it to 0, TRUE to 0.5)
     #axes=T
     #xlim=c(min(l),max(l)),ylim=c(min(l),max(l))#, asp = 0
)


#plot(scale(l2))
#text(scale(l2),labels = ifelse(translationTable2$mainNames=="others",gsub(" Virion","",translationTable2$nodeNames),"" ))

lgd_ = rep(NA, 38)
lgd_[c(38)] = "unkown"
legend(-1.105,1.08,
       legend = lgd_,
       pch=15,
       col =  c(color.gradient(x = 1:36,bg_cols),"white","black"),
       y.intersp = 0.2,
       #title = "Similarity",
       bty="n",
       cex = 1, text.font = 1)
text(-1.105,1.09,"BiGSCAPE family", pos=4)

#plot(1:25,pch=16,col=color.gradient(x = 1:25,c(bg_cols,"black")))

points(c(-1.06,-1.06),c(.249,.249),pch=c(24,25),bg="grey10")
points(-1.06,0.199,pch=16,bg="grey10")

text(-1.105,0.3,"Organism", pos=4)
text(-1.05,.245, "Virion",    pos = 4)
text(-1.05,.195, "Bacteria", pos = 4)

lgd_ = rep(NA, 9)
lgd_[c(1,9)] = c("1.00",signif(min(E(network)$weight),2))
legend("bottomleft",
       legend = lgd_,
       pch=15,
       col =  colorRampPalette(colors = c("red","yellow","grey90"))(9),
       y.intersp = 0.4,
       #title = "Similarity",
       bty="n",
       cex = 1, text.font = 1)
text(-1.105,-0.85,"Similarity", pos=4)

dev.off()


table(subset(translationTable2, genus=="Aeromonas")$bigscape1)

(subset(translationTable2, bigscape1=="unknown"))

enteroCoord=norm_coords(l2)[grep("Escherichia|Klebsiella",translationTable2$genus),]
enteroXY=c(min(enteroCoord[,1]),max(enteroCoord[,1]),min(enteroCoord[,2]),max(enteroCoord[,2]))


png("c:/Users/milst/Desktop/all_phageBGC/Pics/Entero_Virion_vs_IntegratedANI_network.png", width =  11.69 , height = 8.27, units = "in", res=600)
par(mar=c(0,0,0,0))
# Basic chart
enteroPlot=plot(x = network, layout=l2,
                # === vertex
                vertex.color = ifelse(translationTable2$genus=="Virion",translationTable2$cols,adjustcolor(as.character(translationTable2$cols), alpha.f = .5)),   # Node color
                vertex.frame.color = NA,                # Node border color
                vertex.shape=ifelse(translationTable2$genus=="Virion","sphere", "circle"),                        # One of “none”, “circle”, “square”, “csquare”, “rectangle” “crectangle”, “vrectangle”, “pie”, “raster”, or “sphere”
                vertex.size=8,                               # Size of the node (default is 15)
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
                edge.curved=0.3,#      axes=T,                        # Edge curvature, range 0-1 (FALSE sets it to 0, TRUE to 0.5)
                xlim=enteroXY[1:2],ylim=enteroXY[3:4]#, asp = 0
)

#legend("topright", legend = rev(largeGroups), pch=16, col= rev(mycolors$Genus[(largeGroups)]))
dev.off()


mannCoord=norm_coords(l)[grep("Mannheimia|Haemophilus|Pasteurella",translationTable2$genus,value=F),]
mannXY=c(min(mannCoord[,1]),max(mannCoord[,1]),min(mannCoord[,2]),max(mannCoord[,2]))



png("c:/Users/milst/Desktop/all_phageBGC/Pics/Manheim_Virion_vs_IntegratedANI_network.png", width =  11.69 , height = 8.27, units = "in", res=600)
par(mar=c(0,0,0,0))
# Basic chart
MannPlot=plot(x = network, layout=l,rescale=T,
              # === vertex
              vertex.color = ifelse(translationTable2$genus=="Virion",translationTable2$cols,adjustcolor(as.character(translationTable2$cols), alpha.f = .5)),   # Node color
              vertex.frame.color = NA,                # Node border color
              vertex.shape=ifelse(translationTable2$genus=="Virion","sphere", "circle"),                        # One of “none”, “circle”, “square”, “csquare”, “rectangle” “crectangle”, “vrectangle”, “pie”, “raster”, or “sphere”
              vertex.size=8,                               # Size of the node (default is 15)
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
              edge.curved=0.3,     # axes=T,                        # Edge curvature, range 0-1 (FALSE sets it to 0, TRUE to 0.5)
              xlim=mannXY[1:2],ylim=mannXY[3:4]#, asp = 0
)

#legend("topright", legend = rev(largeGroups), pch=16, col= rev(mycolors$Genus[(largeGroups)]))
dev.off()

regex="Myco"

regexCoord=norm_coords(l)[grep(regex,translationTable2$genus,value=F),]
regexXY=c(min(regexCoord[,1]),max(regexCoord[,1]),min(regexCoord[,2]),max(regexCoord[,2]))


# Basic chart
regexPlot=plot(x = network, layout=l, rescale=T,
               # === vertex
               vertex.color = ifelse(translationTable2$genus=="Virion",translationTable2$cols,adjustcolor(as.character(translationTable2$cols), alpha.f = .5)),   # Node color
               vertex.frame.color = NA,                # Node border color
               vertex.shape=ifelse(translationTable2$genus=="Virion","sphere", "circle"),                        # One of “none”, “circle”, “square”, “csquare”, “rectangle” “crectangle”, “vrectangle”, “pie”, “raster”, or “sphere”
               vertex.size=8,                               # Size of the node (default is 15)
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
               edge.curved=0.3,                              # Edge curvature, range 0-1 (FALSE sets it to 0, TRUE to 0.5)
               xlim=regexXY[1:2],ylim=regexXY[3:4]#, asp = 0
)
# 
# dev.off()
# 
# HC=hclust(dist(phageANI), method =  "ward.D2")
# 
# cut=cutree(HC, 4)
# 
# 
# table(cutree(HC, 4))
# 
# 
# enteroIndx=which(cut==2)#c(grep("Escherichia|Citrobacter|Enterobacter|Dickeya|Shigella|Klebsiella",translationTable2$species) )
# 
# enteroANI=phageANI[enteroIndx,enteroIndx]
# enteroTrans=translationTable2[enteroIndx,]
# enteroAnnotaDF=annotaDF[enteroIndx,]
# 
# pheatmap(enteroANI,labels_row = enteroTrans$genome,labels_col = enteroTrans$species,annotation_col = subset(enteroAnnotaDF, select = "Genus"),  annotation_row = subset(enteroAnnotaDF, select = "Genus"),annotation_colors = mycolors,
#          fontsize = 4, clustering_method = "ward.D2")
# 
# # Make an Igraph object from ANI:
# enteroNetwork <- graph_from_adjacency_matrix(as.matrix(enteroANI) , weighted=T, mode="undirected", diag=F )
# 
# 
# png("c:/Users/milst/Desktop/all_phageBGC/Pics/Entero_Virion_vs_IntegratedANI_network.png", width =  11.69 , height = 8.27, units = "in", res=600)
# #pdf(file = "c:/Users/milst/Desktop/all_phageBGC/Pics/Virion_vs_IntegratedANI_network.pdf", width =  11.69 , height = 8.27 )
# 
# # Basic chart
# igraphPlot=plot(x = enteroNetwork, 
#                 # === vertex
#                 vertex.color = ifelse(enteroTrans$genus=="Virion",enteroTrans$cols,adjustcolor(as.character(enteroTrans$cols), alpha.f = .5)),   # Node color
#                 vertex.frame.color = NA,                # Node border color
#                 vertex.shape=ifelse(enteroTrans$genus=="Virion","sphere", "circle"),                        # One of “none”, “circle”, “square”, “csquare”, “rectangle” “crectangle”, “vrectangle”, “pie”, “raster”, or “sphere”
#                 vertex.size=8,                               # Size of the node (default is 15)
#                 vertex.size2=NA,                               # The second size of the node (e.g. for a rectangle)
#                 
#                 # === vertex label
#                 vertex.label=enteroTrans$nodeNames,                 # Character vector used to label the nodes
#                 vertex.label.color="black",
#                 vertex.label.font=2,                          # Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
#                 vertex.label.cex=.8,                          # Font size (multiplication factor, device-dependent)
#                 vertex.label.dist=0,                          # Distance between the label and the vertex
#                 vertex.label.degree=0,                       # The position of the label in relation to the vertex (use pi)
#                 
#                 edge.color=color.gradient(E(enteroNetwork)$weight),   # Edge color
#                 edge.width=E(enteroNetwork)$weight,           # Edge width, defaults to 1
#                 edge.arrow.size=1,                            # Arrow size, defaults to 1
#                 edge.arrow.width=1,                           # Arrow width, defaults to 1
#                 edge.lty="solid",                             # Line type, could be 0 or “blank”, 1 or “solid”, 2 or “dashed”, 3 or “dotted”, 4 or “dotdash”, 5 or “longdash”, 6 or “twodash”
#                 edge.curved=0.3                          # Edge curvature, range 0-1 (FALSE sets it to 0, TRUE to 0.5)
# )
# 
# dev.off()
# 
# 
# #####
# mannheimIndx=which(cut==1)
# 
# mannheimANI=phageANI[mannheimIndx,mannheimIndx]
# mannheimTrans=translationTable2[mannheimIndx,]
# mannheimAnnotaDF=annotaDF[mannheimIndx,]
# 
# pheatmap(mannheimANI,labels_row = mannheimTrans$genome,labels_col = mannheimTrans$species,annotation_col = subset(mannheimAnnotaDF, select = "Genus"),  annotation_row = subset(mannheimAnnotaDF, select = "Genus"),annotation_colors = mycolors,
#          fontsize = 4, clustering_method = "ward.D2")
# 
# # Make an Igraph object from ANI:
# mannheimNetwork <- graph_from_adjacency_matrix(as.matrix(mannheimANI) , weighted=T, mode="undirected", diag=F )
# 
# 
# png("c:/Users/milst/Desktop/all_phageBGC/Pics/mannheim_Virion_vs_IntegratedANI_network.png", width =  11.69 , height = 8.27, units = "in", res=600)
# #pdf(file = "c:/Users/milst/Desktop/all_phageBGC/Pics/Virion_vs_IntegratedANI_network.pdf", width =  11.69 , height = 8.27 )
# 
# # Basic chart
# igraphPlot=plot(x = mannheimNetwork,
#                 # === vertex
#                 vertex.color = ifelse(mannheimTrans$genus=="Virion",mannheimTrans$cols,adjustcolor(as.character(mannheimTrans$cols), alpha.f = .5)),   # Node color
#                 vertex.frame.color = NA,                # Node border color
#                 vertex.shape="circle",                        # One of “none”, “circle”, “square”, “csquare”, “rectangle” “crectangle”, “vrectangle”, “pie”, “raster”, or “sphere”
#                 vertex.size=8,                               # Size of the node (default is 15)
#                 vertex.size2=NA,                               # The second size of the node (e.g. for a rectangle)
#                 
#                 # === vertex label
#                 vertex.label=mannheimTrans$nodeNames,                 # Character vector used to label the nodes
#                 vertex.label.color="black",
#                 vertex.label.font=2,                          # Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
#                 vertex.label.cex=.8,                          # Font size (multiplication factor, device-dependent)
#                 vertex.label.dist=0,                          # Distance between the label and the vertex
#                 vertex.label.degree=0,                       # The position of the label in relation to the vertex (use pi)
#                 
#                 edge.color=color.gradient(E(mannheimNetwork)$weight),#"grey50",                           # Edge color
#                 edge.width=E(mannheimNetwork)$weight,#edge.betweenness(network)*0.01,                                 # Edge width, defaults to 1
#                 edge.arrow.size=1,                            # Arrow size, defaults to 1
#                 edge.arrow.width=1,                           # Arrow width, defaults to 1
#                 edge.lty="solid",                             # Line type, could be 0 or “blank”, 1 or “solid”, 2 or “dashed”, 3 or “dotted”, 4 or “dotdash”, 5 or “longdash”, 6 or “twodash”
#                 edge.curved=0.3,                          # Edge curvature, range 0-1 (FALSE sets it to 0, TRUE to 0.5)
# )
# 
# dev.off()
# 
# 
# 
# 
# plot(x = induced_subgraph(network,vids =  mannheimIndx),
#      # === vertex
#      vertex.color = ifelse(mannheimTrans$genus=="Virion",mannheimTrans$cols,adjustcolor(as.character(mannheimTrans$cols), alpha.f = .5)),   # Node color
#      vertex.frame.color = NA,                # Node border color
#      vertex.shape="circle",                        # One of “none”, “circle”, “square”, “csquare”, “rectangle” “crectangle”, “vrectangle”, “pie”, “raster”, or “sphere”
#      vertex.size=8,                               # Size of the node (default is 15)
#      vertex.size2=NA,                               # The second size of the node (e.g. for a rectangle)
#      
#      # === vertex label
#      vertex.label=mannheimTrans$nodeNames,                 # Character vector used to label the nodes
#      vertex.label.color="black",
#      vertex.label.font=2,                          # Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
#      vertex.label.cex=.8,                          # Font size (multiplication factor, device-dependent)
#      vertex.label.dist=0,                          # Distance between the label and the vertex
#      vertex.label.degree=0,                       # The position of the label in relation to the vertex (use pi)
#      
#      edge.color=color.gradient(E(mannheimNetwork)$weight),#"grey50",                           # Edge color
#      edge.width=E(mannheimNetwork)$weight,#edge.betweenness(network)*0.01,                                 # Edge width, defaults to 1
#      edge.arrow.size=1,                            # Arrow size, defaults to 1
#      edge.arrow.width=1,                           # Arrow width, defaults to 1
#      edge.lty="solid",                             # Line type, could be 0 or “blank”, 1 or “solid”, 2 or “dashed”, 3 or “dotted”, 4 or “dotdash”, 5 or “longdash”, 6 or “twodash”
#      edge.curved=0.3,                          # Edge curvature, range 0-1 (FALSE sets it to 0, TRUE to 0.5)
# )
