layout.by.attr <- function(graph, wc, cluster.strength=1,layout=layout.auto) {  
  g <- graph.edgelist(get.edgelist(graph)) # create a lightweight copy of graph w/o the attributes.
  E(g)$weight <- 1
  
  attr <- cbind(id=1:vcount(g), val=wc)
  g <- g + vertices(unique(attr[,2])) + igraph::edges(unlist(t(attr)), weight=cluster.strength)
  
  l <- layout(g, weights=E(g)$weight)[1:vcount(g),]
  return(l)
}



spreadClus=function(L,transTab, doKleb=T, minSize=16) {
  #find major clusters
  dbclus=dbscan::dbscan(L, eps = 1.5)
  
  #save them with main table
  translationTable2$dbcan=dbclus$cluster
  
  #have a look
  #plot(L, col=0)#translationTable2$dbcan+1)
  
  #make a new layout to overwrite
  l2=L
  
  #klebsiella is to close to escherichia
  #finding indeces of both
  klebIndx=which(translationTable2$dbcan==7)
  entIndx=which(translationTable2$mainNames=="Escherichia")
  
  #finding cluster centers
  klebCent=apply(l[klebIndx,],MARGIN = 2,mean)
  eschCent=apply(l[which(translationTable2$mainNames=="Escherichia"),], MARGIN = 2,mean)
  
  #estimate position of kleb relative to esch
  xx=(klebCent[1]-eschCent[1])
  yy=(klebCent[2]-eschCent[2])  
  #dd=sqrt(xx^2+yy^2)
  
  #move klebsiella away from esch
  #if relative position is positive move more positive and vice versa
  l2[klebIndx,1]=l[klebIndx,1] +c(ifelse(xx>1,2,-2))
  l2[klebIndx,2]=l[klebIndx,2] +c(ifelse(yy>1,2,-2))
  
  #plot again
  plot(L, col=1)#translationTable2$dbcan+1,pch =16, cex =.1)
  
  j=3
  
  #spreading out clusters
  #for all clusters in graph, do
  for(j in unique(translationTable2$dbcan)) {
    if(j==0) next #if cluster is just the non-clustered
    #find index of cluster members and calculate centroid
    dbIndx=which(translationTable2$dbcan==j)
    dbCent=c(mean(l[dbIndx,1]),mean(l[dbIndx,2]))
    
    #if cluster is too small, next
    if(length(dbIndx)<minSize) next
    
    #for all cluster members, do
    for(i in (dbIndx)){
      #find relative position and distance to center
      xx=(l[i,][1]-dbCent[1])
      yy=(l[i,][2]-dbCent[2])  
      dd=sqrt(xx^2+yy^2)
      
      #to each point, add max three but less if cluster is small 
      #to pull away from center
      #l2[i,]=l[i,]+(3-(3/(1*(length(dbIndx))^(1/10))))*c(xx,yy)
      l2[i,]=l[i,]+(0.2*sqrt(length(dbIndx)))*c(xx,yy)
      points(x = l2[i,1],y=l2[i,2], pch=16,cex=0.1,col=translationTable2$dbcan[i]+1)
    }
  }
  
  return(l2)
  
}

# triangle vertex shape
mytriangle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  
  symbols(x=coords[,1], y=coords[,2], bg=vertex.color,
          stars=cbind(vertex.size, vertex.size, vertex.size),
          add=TRUE, inches=FALSE)
}
# clips as a circle
add_shape("triangle", clip=shapes("circle")$clip,
          plot=mytriangle)

# generic star vertex shape, with a parameter for number of rays
mystar <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size  <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  norays <- params("vertex", "norays")
  if (length(norays) != 1 && !is.null(v)) {
    norays <- norays[v]
  }
  
  mapply(coords[,1], coords[,2], vertex.color, vertex.size, norays,
         FUN=function(x, y, bg, size, nor) {
           symbols(x=x, y=y, bg=bg,
                   stars=matrix(c(size,size/2), nrow=1, ncol=nor*2),
                   add=TRUE, inches=FALSE)
         })
}
# no clipping, edges will be below the vertices anyway
add_shape("star", clip=shape_noclip,
          plot=mystar, parameters=list(vertex.norays=5))

# plot(1:60, 4-(3/(1*((1:60))^(1/10))))
# plot(1:60, 0.3*sqrt(1:60))

color.gradient <- function(x, colors=c("grey90","red"), colsteps=100) {
  return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
}