bivPlot = function(y, modelType=NULL, theta=NULL){
 if((is.vector(y)) | (ncol(y) != 2)) stop('bivPlot is only valid for bivariate data.\n')
 if(!is.null(colnames(y))){
  xyLabels = colnames(y)
  labelMargin = 0.8
 } else {
  xyLabels = rep('', 2)
  labelMargin = 0.5
 }
 ghostGray = rgb(0, 0, 0, 0.2)
 #
 nPoints = 80
 y1hist = hist(y[,1], plot=F)
 y2hist = hist(y[,2], plot=F)
 grid.y1 = seq(min(y1hist$breaks), max(y1hist$breaks), length=nPoints)
 grid.y2 = seq(min(y2hist$breaks), max(y2hist$breaks), length=nPoints)
 #
 layout.mat = matrix(c(1,0,2,3),2,2,byrow=T)
 layout(layout.mat,widths=c(4,1),heights=c(1,4))
 #
 par(mai=c(0.01,labelMargin,0.1,0.1))
 if(!is.null(theta)){
  y1dens = dmvSE(grid.y1, modelType=modelType, theta=list(xi=theta$xi[1], G=theta$G[1,1,drop=F], psi=theta$psi[1], nu=theta$nu))
 } else y1dens = numeric(0)
 plot(y1hist, ylim=c(0, max(c(y1hist$density, y1dens))), freq=F, main='', ylab='', col=rgb(0,0,1,0.3,maxColorValue=1), border='white', axes=F)
 if(!is.null(theta)) lines(grid.y1, y1dens, col='blue')
 #
 par(mai=c(labelMargin,labelMargin,0.01,0.01))
 if(!is.null(theta)){
  yValues = expand.grid(grid.y1, grid.y2)
  pdfVec = dmvSE(y=yValues, modelType=modelType, theta=theta, LOG=F)
  pdf = matrix(pdfVec, nPoints, nPoints)
  image(grid.y1, grid.y2, pdf, xlab=xyLabels[1], ylab=xyLabels[2], col=colorRampPalette(c(rgb(0,0,1,0), rgb(0,0,1,1)), alpha = TRUE)(120))
  contour(grid.y1, grid.y2, pdf, add=T)
  points(y[,1], y[,2], pch=19, col=ghostGray)
 } else {
  plot(y[,1], y[,2], pch=19, col=ghostGray)
 }
 points(y[,1], y[,2])
 #
 par(mai=c(labelMargin,0.01,0.01,0.01))
 if(!is.null(theta)){
  y2dens = dmvSE(grid.y2, modelType=modelType, theta=list(xi=theta$xi[2], G=theta$G[2,2,drop=F], psi=theta$psi[2], nu=theta$nu))
 } else y2dens = numeric(0)
 plot(c(0,max(c(y2hist$density,y2dens))), range(y2hist$breaks), axes=F, frame.plot=F, type='n', xlab='', ylab='')
 for(bin in 1:(length(y2hist$breaks)-1)) polygon(c(0, y2hist$density[bin], y2hist$density[bin], 0), c(y2hist$breaks[bin], y2hist$breaks[bin], y2hist$breaks[bin+1], y2hist$breaks[bin+1]), col=rgb(0, 0, 1, 0.3, maxColorValue=1), border='white')
 if(!is.null(theta)) lines(y2dens, grid.y2, col='blue')
}
