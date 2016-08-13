plotbar <- 
  function(x, y, z, xlim = NULL, ylim = NULL, zlim = NULL, zlab = NULL, 
           zcex.axis = NULL, zcex.lab = NULL, zaxis.at = NULL,
           zaxis.labels = TRUE, col = NULL, ncolor = 21, 
           drawbar = TRUE, zline = 2, pltimage = c(0.2, 0.8, 0.2, 0.8), 
           pltbar = c(0.82, 0.85, 0.2, 0.8), ...){
    
    if (is.null(col[1])) {
      col <- c("blueviolet", "blue", "cyan", "green", "yellow", "orange", "red")
    }
    col <- colorRampPalette(col)(ncolor)
    if (is.null(zlim)) {
      zlim <- range(c(z))
    }
    if (is.null(zlab)) {
      zlab <- "z"
    }
    if (is.null(zcex.lab)) {
      zcex.lab <- 1
    }
    scales <- ncolor/(zlim[2] - zlim[1])
    breaks <- seq(zlim[1], zlim[2], length.out = (ncolor + 1))
    oldplt <- par()$plt
    on.exit(par(plt = oldplt))
    
    if (drawbar) {
      par(plt = pltbar)
      plot(1, 1, t = "n", ylim = zlim, xlim = c(0, 1), xaxt = "n", 
           yaxt = "n", xaxs = "i", yaxs = "i", xlab = "", ylab = "")
      axis(4, at = zaxis.at, labels=zaxis.labels, cex.axis = zcex.axis)
      mtext(zlab, side = 4, line = zline, cex = zcex.lab * par()$cex)
      for (ii in 1:ncolor) {
        idx <- zlim[1] + (ii - 1)/scales
        rect(0, idx, 1, (idx + 1/scales), col = col[ii], border = NA)
      }
      par(plt = pltimage, new = TRUE)
      plot(x, y, col = num2col(z,col,zlim), xlim=xlim, ylim=ylim, ...)
    }
    else {
      par(plt = pltimage)
      plot(x, y, col = num2col(z,col,zlim), xlim=xlim, ylim=ylim, ...)
    }
    
  }