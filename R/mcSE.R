mcSE = function(y, X=NULL, N, Ti, modelType='ST', warmUp=FALSE, control=list()){

 # Checks
 if(is.data.frame(y)) y = as.matrix(y)
 if(is.vector(y)) y = matrix(y, ncol=1)
 dimnames(y) = NULL
 n = nrow(y)
 p = ncol(y)
# if(is.null(p)) stop('The mcSE function is not yet implemented for univariate data.\n')
 XFlag = !is.null(X)
 if(XFlag){
#  if(!is.matrix(X)) X = as.matrix(X)
  attr(X, 'dimnames') = NULL
  if(!is.matrix(X)) stop('mcSE.R: the object X is not a matrix.\n')
  if(nrow(X) != n) stop('mcSE.R: matrices y and X should have the same number of rows.\n')
#  k = ncol(X)
  if(!identical(X[,1], rep(1,n))) warning('mcSE.R: the first column of X does not contain a constant term. A regression model without intercept will be estimated.\n')
 }

 # Control list
 con = list(seed=NULL, propFuncs=NULL, logPriorFunc=NULL, Nwu=c(N/10, N/5, N/2), priorList=NULL, saveParticles=FALSE, outFolder='Output', verbose=TRUE, parInfo=NULL)
 conNames = names(con)
 controlNames = names(control)
 if(length(control) > 0) con[controlNames] = control
 if(length(missNames <- controlNames[!controlNames %in% conNames]) > 0) 
 warning('mcSE.R: unknown objects in control list: ', paste(missNames, collapse=', '))
 #
 saveParticles = con$saveParticles
 if(saveParticles){
  outFolder = con$outFolder
  if(!(outFolder %in% dir())){
   dir.create(outFolder, recursive=T)
  }
 }
 Nwu = con$Nwu
 verbose = con$verbose

 # List of parameters, by type
 if((is.null(con$parInfo)) & (modelType %in% c('N','T','SN','ST'))){
#  parTypes = modelParTypes(modelType, XFlag)
  parInfo = modelParInfo(y, X, modelType)
 } else {
#  parTypes = control$parTypes
  parInfo = con$parInfo
 }

 # Proposal distributions
 if('propFuncs' %in% controlNames){
  if(!(all(names(con$propFuncs) %in% parInfo$names) & (length(con$propFuncs) == nrow(parInfo)))){
   stop('mcSE.R: proposal distributions are badly specified.\n')
  } else {
   propFuncs = control$propFuncs
  }
 } else {
  propFuncs = paste('sample', parInfo$names, modelType, sep='')
  names(propFuncs) = parInfo$names
 }

 # Prior density
 if('logPriorFunc' %in% controlNames){
  logPriorFunc = control$logPriorFunc
 } else {
  logPriorFunc = paste('logPriorDens', modelType, sep='')
 }

 # Posterior density
 logPostFunc = paste('logPostDens', modelType, sep='')

 # Hyperparameters
 if(is.null(con$priorList)){
  nuVals = 1:30 # c(1, 2, 3, 4, 5, 10, 15, 20, 30, 100)
  nuLogProbs = log(villaPrior(d=p, nuMax=max(nuVals)))
  priorList = list(m=0, W = matrix(0, p, p), gdls = nuVals, nulogpriors = nuLogProbs)
  rm(nuVals, nuLogProbs)
 } else {
  priorList = con$priorList
 }

 # Initialization
 if(!is.null(con$seed)) set.seed(con$seed)
 log.py = numeric(Ti)
 nResampled = numeric(Ti)
 perplexity = numeric(Ti)
 estList = estListInit(parInfo, Ti)
 iterEst = NULL 		# Object in iterSummary
 log.py.ti = NULL 		# Object in iterSummary
 nResampled.ti = NULL  	# Object in iterSummary
 perplexity.ti = NULL	# Object in iterSummary

 # Initialization of the particles
 initialFuncName = paste('initialPoints', modelType, sep='')
 if(warmUp == T) Ninit = Nwu[1] else Ninit = N
 particles = do.call(initialFuncName, list(N=Ninit, y=y, X=X, priorList=priorList))

 # Warm-up iterations
 if(warmUp){
  if(verbose) cat('\nWarm-up','\n')
  Tiwu = length(Nwu)
  NwuStar = c(Nwu, N)
  for(ti in 1:Tiwu){
   if(verbose) cat('\n', ti,'  ')
   particles = WUiterations(particles, NwuStar[ti], NwuStar[ti+1], ti, Ti, y, X, modelType, priorList, parInfo, propFuncs, logPriorFunc, logPostFunc, saveParticles, outFolder, verbose)
  }
 }

 # Iterations 1:Ti
 if(verbose) cat('\nIterations','\n')
 for(ti in 1:Ti){
  if(verbose) cat('\n', ti,'  ')
  iterSummary = iterations(particles, N, ti, Ti, y, X, modelType, priorList, parInfo, propFuncs, logPriorFunc, logPostFunc, saveParticles, outFolder, verbose)
  for(iItem in 1:length(iterSummary)) assign(names(iterSummary)[iItem], iterSummary[[iItem]])
#  iterWeights[ti] = iterWeight	# From iterSummary
  estList = estListUpdate(estList, iterEst, parInfo, ti)	# From iterSummary
  log.py[ti] = log.py.ti			# From iterSummary
  nResampled[ti] = nResampled.ti	# From iterSummary
  perplexity[ti] = perplexity.ti	# From iterSummary
  rm(iterSummary)
 }
 if(verbose) cat('\n')

# mcSE.output = append(estList, list(log.py=log.py, nResampled=nResampled, perplexity=perplexity, parInfo=parInfo))
 mcSE.output = list(estlist=estList, log.py=log.py, nResampled=nResampled, p=p, perplexity=perplexity, parInfo=parInfo)
 class(mcSE.output) = 'mcSE'

 # print.summary.mcSE(summary.mcSE(mcSE.output))
 return(mcSE.output)
}
