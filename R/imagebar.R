imagebar <-
  function(x,y,z,xlim=NULL,ylim=NULL,zlim=NULL,
           zlab=NULL,zcex.axis=NULL,zcex.lab=NULL,
           col=NULL,ncolor=100,drawbar=TRUE,...){
    
    ### define colors
    if(is.null(col[1])){
      col=rev(rainbow(ncolor,end=3/4))
    } else {col=colorRampPalette(col)(ncolor)}
    
    ### get z limits and label
    if(is.null(zlim)){zlim=range(c(z))}
    if(is.null(zlab)){zlab="z"}
    if(is.null(zcex.lab)){zcex.lab=1}
    scales=(ncolor-1L)/(zlim[2]-zlim[1])
    breaks=seq(zlim[1],zlim[2],length.out=(ncolor+1))
    
    # plot image and color bar
    if(drawbar){
      # plot image 
      par(plt=c(.13,.73,.23,.83))
      image(x,y,z,col=col,breaks=breaks,...)
      # plot color bar
      par(plt=c(.8,.85,.23,.83),new=TRUE)
      plot(1,1,t="n",ylim=zlim,xlim=c(0,1),xaxt="n",yaxt="n",xaxs="i",yaxs="i",xlab="",ylab="")
      axis(4,cex.axis=zcex.axis); mtext(zlab,side=4,line=3,cex=zcex.lab*par()$cex)
      for (ii in 1:(ncolor-1L)) {
        idx=zlim[1]+(ii-1)/scales
        rect(0,idx,1,(idx+1/scales),col=col[ii],border=NA)
      }
    } else {
      # plot image
      par(plt=c(.2,.8,.23,.83))
      image(x,y,z,col=col,breaks=breaks,...)
    }
    
  }