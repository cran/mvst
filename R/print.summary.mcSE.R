print.summary.mcSE = function(x){
 digits = max(3L, getOption('digits') - 3L)
 out = cbind(unlist(x$sampleEstimates$postMean), unlist(x$sampleEstimates$postMeanSE), unlist(x$sampleEstimates$Q05), unlist(x$sampleEstimates$postMedian), unlist(x$sampleEstimates$Q95))
 colnames(out) = c('Estimate', 'Std. Error', 'Q5%', 'Me', 'Q95%')
 cat('\n')
 print.default(out, digits=digits, print.gap=2L)
 cat('\n', 'log(p(y)) =', rep(x$log.margLike), '\n')
# invisible(x)
}
