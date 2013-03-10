# Returns a square matrix m such that m[i,j] is the number of times a direct transition from state i to j has been observed
.countTransitions<-function(x,nestados){
    resultado = matrix(0,nestados,nestados);
    n = length(x);
    for(i in 1:n-1){
      resultado[x[i],x[i+1]] = resultado[x[i],x[i+1]] + 1;
    }
    return(resultado);
}

## _____________________________________________________________________________________________________________

# Computes the matrix of intervals representing alpha-cuts of the fuzzy probabilities. Returns a matrix of
# intervals lower bounds, a matrix of intervals upper bounds, a list of row indices and a list of column indices (negative indicates last non-null element in a row)
# indicating the entries of the matrix where the probabilities are greater than 0, and a crisp punctual estimation of the transition probabilities according to observed data.
.computeIntervalMatrices<-function(x,nestados,a){
    
  countsMatrix = .countTransitions(x,nestados);
  lowerBoundsMatrix = matrix(0,nestados,nestados);
  upperBoundsMatrix = matrix(0,nestados,nestados);
  punctualEstimatesMatrix = matrix(0,nestados,nestados);
  
  listaposx = {}
  listaposy = {}
  
  iniciales = {};
  
  lower = {};
  upper = {};
    
  for(i in 1:nestados){
    if(i > 1){
      listaposy[length(listaposy)] = (-1)*listaposy[length(listaposy)];          
      # Remove the last variable of this probability distribution (it is determined as 1 - sum of the rest)
      iniciales = iniciales[-length(iniciales)];
      upper = upper[-length(upper)];
      lower = lower[-length(lower)];
    }
    
    suma = sum(countsMatrix[i,]);
    if(suma > 0){
      temp = multinomialCI(countsMatrix[i,],a,verbose=FALSE);
    }
    else{
      temp = matrix(0, nestados, 2);
    }
        
    # Traverse all the row
    for(j in 1:nestados){ 
      if(countsMatrix[i,j] > 0){
        lowerBoundsMatrix[i,j] = temp[j,1];
        upperBoundsMatrix[i,j] = temp[j,2];
        listaposx = append(listaposx,i);
        listaposy = append(listaposy,j);
        iniciales = append(iniciales, (lowerBoundsMatrix[i,j] + upperBoundsMatrix[i,j])/2); # IGNORADO
        if(lowerBoundsMatrix[i,j] == upperBoundsMatrix[i,j] && lowerBoundsMatrix[i,j] != 0){
          lower = append(lower, lowerBoundsMatrix[i,j]-0.01);
        }else{
          lower = append(lower, lowerBoundsMatrix[i,j]);
        }
        upper = append(upper, upperBoundsMatrix[i,j]);      
      }
    }

    if(suma > 0){      
      punctualEstimatesMatrix[i,] = countsMatrix[i,]/suma;
    }
    else{
      punctualEstimatesMatrix[i,] = rep(0, nestados);
    }
  }  
  # Last non-zero element of the last row
  listaposy[length(listaposy)] = (-1)*listaposy[length(listaposy)];    
  iniciales = iniciales[-length(iniciales)];
  upper = upper[-length(upper)];
  lower = lower[-length(lower)];
  
  # tomar como solución inicial la matriz de transición observada (frecuencias) en lugar del punto medio de cada intervalo
  antiguo = length(iniciales);
  iniciales = {};
  for(indice in 1:length(listaposy)){
    if(listaposy[indice] >= 0){
      iniciales = append(iniciales, punctualEstimatesMatrix[listaposx[indice],listaposy[indice]]);
    }
  }

  # Stationary probabilities using punctual estimation of the transition matrix
  identidad<-matrix(0,nestados,nestados);
  diag(identidad)<-1;
  ceros<-matrix(0,nestados,nestados);
  mat2=.unos(ceros) %*% solve(.unos(punctualEstimatesMatrix-identidad));
  estac=mat2[1,];
  
  return(list(left=lowerBoundsMatrix,right=upperBoundsMatrix,listaposx=listaposx,listaposy=listaposy,lower=lower,upper=upper,iniciales=iniciales,puntuales=estac));
}

## _____________________________________________________________________________________________________________

## Computes the left and right alpha-cut bounds (with the specified alpha argument) of the stationary probabilities for the
## states indicated. Returns a list with the left bounds of the alpha-cut(s) and another list with the right bounds. 
## The size of both lists matches the length of argument states
.computeFuzzyStationary<-function(alfa,data,d,states){
  
  res=.computeIntervalMatrices(data,d,alfa);
  stationaryLeft = {};
  stationaryRight = {};
  
  if(alfa < 0.999){
    for(j in 1:length(states)){    
       state = states[j];
       pob1 = .generateInitialPopulation(length(res$lower), res$left, res$right);
       pob2 = .generateInitialPopulation(length(res$lower), res$left, res$right);
       opt = DEoptim(.stationary, res$lower, res$upper, control = DEoptim.control(
            trace=FALSE,itermax=200,CR=0.8,initialpop=pob1,reltol=1E-3,steptol=20), d, res$listaposx, res$listaposy,component=state,res$left,res$right,FALSE);
       stationaryLeft = append(stationaryLeft, opt$optim$bestval);
       opt = DEoptim(.stationary, res$lower, res$upper, control = DEoptim.control(
            trace=FALSE,itermax=200,CR=0.8,initialpop=pob2,reltol=1E-3,steptol=20), d, res$listaposx, res$listaposy,component=state,res$left,res$right,TRUE); 
       stationaryRight = append(stationaryRight, (-1)*opt$optim$bestval);
    }
  }
  else{ # Punctual estimations for both lower and upper bounds (no need to launch DE)
    for(j in 1:length(states)){
      stationaryLeft = append(stationaryLeft, res$puntuales[states[j]]);
      stationaryRight = append(stationaryRight, res$puntuales[states[j]]);
    }
  }

  return(list(left=stationaryLeft, right=stationaryRight, puntuales=res$puntuales));
}

## _____________________________________________________________________________________________________________

.computeOffsetTable<-function(data, states){
  
  maximum = max(c(max(data),max(states)));
  res = rep(-1,maximum);
  for(k in 1:length(data)){
    res[data[k]] = 1;
  }
  
  offset = 0;
  for(i in 1:maximum){
    if(res[i] != -1){
      res[i] = offset;      
    }
    else{
      res[i] = NA;
      offset = offset+1;
    }
  }
  return(res);
}

## _____________________________________________________________________________________________________________

.subtractOffsets<-function(data, states, offsetTable){

  correcteddata = rep(0, length(data));
    
  maximum = max(data);
  if(is.null(states)){
    states = seq(1:maximum);
  }
  
  condensedstates = {};
  condensedoffsets = {};
  
  correctedstates = rep(0, length(states));
  
  for(i in 1:length(data)){
    correcteddata[i] = data[i] - offsetTable[data[i]];    
  }
  for(i in 1:length(states)){
    if(!is.na(offsetTable[states[i]])){
      correctedstates[i] = states[i] - offsetTable[states[i]];
      condensedstates = append(condensedstates, correctedstates[i]);
      condensedoffsets = append(condensedoffsets, offsetTable[states[i]]);
    }
    else{
      correctedstates[i] = NA;
    }
  }
  
  return(list(newdata = correcteddata, newstates = correctedstates, condensedstates = unlist(condensedstates), condensedoffsets = unlist(condensedoffsets)));
}

## _____________________________________________________________________________________________________________

fuzzyStationaryProb<-function(data,options,step=0.05){

## -----------------------------------------------------------------
##                          CHECK ARGUMENTS 
## -----------------------------------------------------------------

  t1 = proc.time();

  r = names(options);  
  res = sapply(r,function(x){
    if(!(x=="verbose" || x=="regression" || x=="states" || x=="parallel" || x=="acutsonly" || x=="cores")){
      stop("ERROR: unrecognized option : ", x,"\n");
    }
  });  

  data = sapply(data,as.integer);
  options$states = sapply(options$states, as.integer);
  if(sum(data<=0) > 0){
    stop("ERROR: zero or negative numbers are not allowed in the data vector\n");
  }
  if(!is.null(options$states)){ 
    if(sum(options$states <= 0) > 0){
      stop("ERROR: zero or negative numbers are not allowed among the states for which stat.probs. should be computed\n");      
    }
  }
  
  offsets = .computeOffsetTable(data, options$states);
  templist = .subtractOffsets(data, options$states, offsets);  
  
  data = templist$newdata;
  options$states = templist$condensedstates;
  
  stopifnot(min(data) == 1);
  nstates = max(data);
  nresults = length(options$states);

  if(is.null(options$regression)){
    options$regression = "quadratic";
  }
  else{    
    options$regression = tolower(options$regression);
    reg = options$regression;
    if(!(reg == "linear" || reg == "quadratic" || reg == "cubic" || reg == "gaussian" || reg == "spline" || reg == "piecewise")){
      stop("ERROR: regression must be one of: linear | quadratic | cubic | gaussian | spline | piecewise\n");
    }
  }
  
  if(is.null(options$parallel)){
    options$parallel = FALSE;
  }
  else{
    if(!(options$parallel == TRUE || options$parallel == FALSE)){
      stop("ERROR: the value of options$parallel must be either TRUE or FALSE\n");
    }
  }
  if(options$parallel == TRUE){
    if(is.null(options$cores)){
      options$cores = -1;
    }
    else{
      options$cores = as.integer(options$cores);
      if(options$cores <= 0){
        stop("ERROR: the value of options$cores must be a positive integer");
      }
    }
  }
  
  if(is.null(options$verbose)){
    options$verbose = FALSE;
  }
  
  if(is.null(options$acutsonly)){
    options$acutsonly = FALSE;    
  }
  else{
    if(!(options$acutsonly == TRUE || options$acutsonly == FALSE)){
      stop("ERROR: the value of options$acutsonly must be either TRUE or FALSE\n");
    }
  }
  

## ---------------------------------------------------------------------
##          COMPUTE ALPHA-CUTS FOR DESIRED STATES 
## ---------------------------------------------------------------------

  iterations = 1+1.0/step;  
  listaconf = rep(NA,iterations);  
  pointlistleft = list();
  pointlistright = list();

  for(j in 1:nresults){
    pointlistleft[[j]] = rep(NA,iterations);
    pointlistright[[j]] = rep(NA,iterations);
  }

  for(i in 1:iterations){
    if(i == 1){ confidence = 0.001;   }
    else{ 
      if(i == iterations){ confidence = 0.999; }
      else{ confidence = step * (i-1); }
    }    
    listaconf[i] = confidence;
  }
  ncores = 1;
 
  ## PARALLEL -------------------------------------------
  if(options$parallel){    
    ncores = detectCores();
    if(options$cores > 0){
      ncores = min(ncores, options$cores);
    }
    cl = makeCluster(ncores);
    clusterExport(cl,list());
    clusterEvalQ(cl, library(MultinomialCI));
    clusterEvalQ(cl, library(DEoptim));
    if(options$verbose){ cat("Parallel computation of a-cuts in progress"); flush.console(); }
    res = parLapply(cl,listaconf,.computeFuzzyStationary,data,nstates,options$states);
    stopCluster(cl);
    for(k in 1:iterations){
      for(j in 1:nresults){
        pointlistleft[[j]][k] = res[[k]]$left[j];
        pointlistright[[j]][k] = res[[k]]$right[j];
      }
    }
  }
  ## SEQUENTIAL -----------------------------------------
  else{  
    if(options$verbose){ cat("Computing a-cuts for a = "); flush.console();}  
    for(i in 1:iterations){
      confidence = listaconf[i];
      if(options$verbose){ cat(confidence," "); flush.console(); }      
      res = .computeFuzzyStationary(confidence,data,nstates,options$states);
      for(j in 1:nresults){
        pointlistleft[[j]][i] = res$left[j];
        pointlistright[[j]][i] = res$right[j];
      }
    }
  }
    
  t2 = proc.time();
  if(options$verbose) cat("...finished successfully (elapsed: ",t2[3]-t1[3],"s.)\n");

## ---------------------------------------------------------------------
## USE REGRESSION TO OBTAIN THE MEMBERSHIP FUNCTIONS FITTING THE ALPHA-CUTS
## ---------------------------------------------------------------------
  
  if(options$verbose){ 
    cat("Applying",options$regression,"regression to obtain the membership functions of fuzzy stationary probabilities...\n");
    cat("Fitting memb.func. for state: "); flush.console();
  }
  # Now do the regression and finally build the fuzzy numbers defined by the alpha-cuts computed above
  fuzzylist = list();
  misframes = list();
  fuzzylist[[nresults]]=0;
  misframes[[nresults]]=0;  

  for(i in 1:nresults){    
    if(options$verbose){ cat(options$states[i]+templist$condensedoffsets[i]," ");    flush.console(); }
    
    myframeleft = data.frame(pointlistleft[[i]],listaconf);
    myframeright = data.frame(pointlistright[[i]],listaconf);    
    
    if(options$acutsonly){
        names(myframeleft) = c("x","y");
        names(myframeright) = c("x","y");
        unframe = rbind(myframeleft, myframeright);
        misframes[[i]] = unframe;
    }
    else{

      mylist = .membFnRegression(myframeleft, myframeright, options$regression);

      if(options$verbose && i==nresults) cat("\n");    
      if(!is.null(mylist)){
        fuzzylist[[i]] = mylist[[1]];
        myframeleft = mylist[[2]];
        myframeright = mylist[[3]];
        names(myframeleft) = c("x","y");
        names(myframeright) = c("x","y");
        unframe = rbind(myframeleft, myframeright);
        misframes[[i]] = unframe;
      }else{
        fuzzylist[[i]] = NA;
        misframes[[i]] = NA;
      }
    }
  }
  
## ---------------------------------------------------------------------
## COMPOSE OUTPUT LISTS INCLUDING NA IN POSITIONS OF NON-EXISTENT STATES
## ---------------------------------------------------------------------
    
    other = 1;
    outputfuzzylist = list();
    outputmyframes = list();
    for(i in 1:length(templist$newstates)){
      if(is.na(templist$newstates[i])){
        #outputmyframes = append(outputmyframes, NA);
        outputmyframes[[i]] = NA;
        if(!options$acutsonly){
          #outputfuzzylist = append(outputfuzzylist, NA);
          outputfuzzylist[[i]] = NA;
        }
      }
      else{
        #outputmyframes = append(outputmyframes, misframes[[other]]);
        outputmyframes[[i]] = misframes[[other]];
        if(!options$acutsonly){
          #outputfuzzylist = append(outputfuzzylist, fuzzylist[[other]]);
          outputfuzzylist[[i]] = fuzzylist[[other]];
        }
        other = other+1;
      }
    }
  
## ---------------------------------------------------------------------  
  
  t3 = proc.time();
  
  msg = paste(     ". Fuzzy stationary probabilities of a Markov chain with",nstates,"states\n");
  msg = paste(msg, ". Probabilities have been computed for states:",paste((options$states+templist$condensedoffsets),collapse=" "),"\n");
  msg = paste(msg, ". Number of input observations:",length(data),"\n");
  msg = paste(msg, ". Parameters:\n");
  msg = paste(msg, "       Step size:",step,"\n");
  if(options$parallel){
    msg = paste(msg, "       Execution was done in parallel (",ncores,"cores used )\n");
  }
  else{
    msg = paste(msg, "       Execution was done sequentially\n");
  }  
  if(!options$acutsonly){
    msg = paste(msg, "       Regression curve for output membership functions:",options$regression);  
  }
  else{
    msg = paste(msg, "       (a-cuts only; no membership function was computed)\n");
  }
  msg = paste(msg, ". To retrieve the results use $fuzzyStatProb and $acuts with this object\n");
  msg = paste(msg, ". Computation of alpha-cuts took", format(round(t2[3] - t1[3],2),nsmall=2),"seconds\n");
  if(!options$acutsonly) msg = paste(msg, ". Membership functions regression took", format(round(t3[3] - t2[3],2),nsmall=2),"seconds\n");
 
  if(options$acutsonly){ result = list(acuts=outputmyframes, summary=msg); }
  else{ result = list(fuzzyStatProb=outputfuzzylist, acuts=outputmyframes, summary=msg); }
  
  class(result)<-"FuzzyStatObj";
  return(result);
  #return(list(listaleft=pointlistleft,listaright=pointlistright,listaconf=listaconf));
}