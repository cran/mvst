mcSE = function(y, N, Ti, modelType='ST', warmUp=FALSE, control=list()){

 # Checks
 p = ncol(y)
 if(is.null(p)) stop('The mcSE function is not yet implemented for univariate data.\n')

 # Control list
 con = list(seed=NULL, simFlag=FALSE, propFuncs=NULL, postDensFunc=NULL, Nwu=rep(2000,3), saveParticles=FALSE, outFolder='Output', verbose=TRUE)
 conNames = names(con)
 controlNames = names(control)
 if(length(control) > 0) con[controlNames] = control
 if(length(missNames <- controlNames[!controlNames %in% conNames]) > 0) 
 warning('Unknown objects in control list: ', paste(missNames, collapse=', '))
 #
 saveParticles = con$saveParticles
 if(saveParticles){
  outFolder = con$outFolder
  if(!(outFolder %in% dir())){
   dir.create(outFolder, recursive=T)
  }
  if(!('Iterations' %in% dir(outFolder))){
   dir.create(paste(outFolder, '/Iterations', sep=''), recursive=T)
  }
 }
 Nwu = con$Nwu
 verbose = con$verbose

 # List of parameters, by type
 if(modelType %in% c('N','T','SN','ST')){
  parTypes = modelParTypes(modelType)
 } else {
  parTypes = control$parTypes
 }

 # Proposal distributions
 if('propFuncs' %in% controlNames){
  if(!(all(names(con$propFuncs) %in% parTypes) & (length(con$propFuncs) == length(parTypes)))){
   stop('Proposal distributions are badly specified.\n')
  } else {
   propFuncs = control$propFuncs
  }
 } else {
  propFuncs = paste('sample', parTypes, modelType, sep='')
  names(propFuncs) = parTypes
 }

 # Prior density
 if('logPriorFunc' %in% controlNames){
  logPriorFunc = control$postDensFunc
 } else {
  logPriorFunc = paste('logPriorDens', modelType, sep='')
 }

 # Posterior density
 logPostFunc = paste('logPostDens', modelType, sep='')

 # Hyperparameters
 if(is.null(con$priorList)){
  nuVals = c(1,2,3,4,5,6,8,10,12,14,16,18,20,30,40,50,60,80,100)
  priorList = list(m=0, W = matrix(0, p, p), gdls = nuVals, nulogpriors = rep(-log(length(nuVals)), length(nuVals)))
 rm(nuVals)
 } else {
  priorList = con$priorList
 }

 # Initialization
 if(!is.null(con$seed)) set.seed(con$seed)
 log.py = nResampled = perplexity = numeric(Ti)
 estList = estListInit(parTypes)
 iterEst = NULL		# Object in iterSummary
 log.py.ti = NULL		# Object in iterSummary
 nResampled.ti = NULL	# Object in iterSummary
 perplexity.ti = NULL	# Object in iterSummary

 # Initialion of the particles
 initialFuncName = paste('initialPoints', modelType, sep='')
 if(warmUp == T) Ninit = Nwu[1] else Ninit = N
 particles = do.call(initialFuncName, list(N=Ninit, y=y, priorList=priorList))

 # Warm-up iterations
 if(warmUp){
  if(verbose) cat('\nWarm-up','\n')
  Tiwu = length(Nwu)
  NwuStar = c(Nwu, N)
  for(ti in 1:Tiwu){
   if(verbose) cat('\n', ti,'  ')
   particles = WUiterations(particles, NwuStar[ti], NwuStar[ti+1], ti, Ti, y, modelType, priorList, parTypes, propFuncs, logPriorFunc, logPostFunc, saveParticles, outFolder, verbose)
  }
 }

 # Iterations 1:Ti
 if(verbose) cat('\nIterations','\n')
 for(ti in 1:Ti){
  if(verbose) cat('\n', ti,'  ')
  iterSummary = iterations(particles, N, ti, Ti, y, modelType, priorList, parTypes, propFuncs, logPriorFunc, logPostFunc, saveParticles, outFolder, verbose)
  for(iItem in 1:length(iterSummary)) assign(names(iterSummary)[iItem], iterSummary[[iItem]])
#  iterWeights[ti] = iterWeight	# From iterSummary
  estList = estListUpdate(estList, iterEst, ti)	# From iterSummary
  log.py[ti] = log.py.ti			# From iterSummary
  nResampled[ti] = nResampled.ti	# From iterSummary
  perplexity[ti] = perplexity.ti	# From iterSummary
  rm(iterSummary)
 }

 mcSE.output = list(estList=estList, log.py=log.py, nResampled=nResampled, perplexity=perplexity)
 class(mcSE.output) = 'mcSE'

 return(mcSE.output)
}
