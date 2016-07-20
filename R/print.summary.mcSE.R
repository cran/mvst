print.summary.mcSE = function(x){
 digits = max(3L, getOption('digits') - 3L)
 out = cbind(unlist(x$sampleEstimates$postMean), unlist(x$sampleEstimates$postMeanSE))
 colnames(out) = c('Estimate', 'Std. Error')
 print.default(out, digits=digits, print.gap=2L)
 cat('\n', 'log(p(y)) =', rep(x$log.margLike), '\n')
 invisible(x)
}
