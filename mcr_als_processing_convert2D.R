
convert_to_2D<- function(hypercube){
  length<-dim(hypercube)[1]
  width<-dim(hypercube)[2]
  num_bands<-dim(hypercube)[3]
  convert2D<-matrix(hypercube, nrow=length*width,ncol=num_bands)
  return(convert2D)
}