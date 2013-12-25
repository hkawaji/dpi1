


##########################################################
### internal function: myFastICA (fastICA wrapper)     ###
##########################################################

myFastICA <- function( data, cumulative_proportion = 0.95, seeds=11:13, hclust_plot = F, verbose = T, n.comp.upper_bound = Inf, ...)
{
  if ( ncol(data) == 1 ){return(list())}

  ### decide starting number of components
  n.comp = myFastICA_estimate_n.comp_by_PCA(data, cumulative_proportion, n.comp.upper_bound, verbose)
  if ( n.comp == 1 ){return(list())}

  #if (verbose){ cat(sprintf("#starting number of components: %d\n", n.comp)) }
  #if ( n.comp > n.comp.upper_bound ){ n.comp = n.comp.upper_bound }

  ### perform ICA with different seeds
  res = lapply(
    seeds,
    function(seed) myFastICA_from_a_seed( data, n.comp, seed=seed , ...)
  )
  ### in case that a fastICA trial failed.
  if ( any(sapply(res, function(x) length(x) == 0)) ){return(list())}

  ### in case only 1 seed
  if (length(res) == 1){return(res[[1]])}

  ### reduce the number of components, if noise components appeard consistently
  number_of_signal_components = table( sapply(res, function(ica) length( ica$signal_components ) ) )
  n.comp.new = as.integer( names(number_of_signal_components)[ which.max( number_of_signal_components ) ] )
  while ( ( n.comp.new < n.comp ) & (n.comp.new >= 2) )
  {
    if(verbose){cat( sprintf("###components: %d -> %d\n", n.comp , n.comp.new) )}
    n.comp = n.comp.new
    res = lapply(
      seeds,
      function(seed) myFastICA_from_a_seed( data, n.comp, seed=seed, ...)
    )
    number_of_signal_components = table( sapply(res, function(ica) length( ica$signal_components ) ) )
    n.comp.new = as.integer( names(number_of_signal_components)[ which.max( number_of_signal_components ) ] )
  }
  if ( n.comp.new == 1 ){return(list())}

  ### hclust for the independent components from the distinct seeds
  if (hclust_plot == T)
  {
    S = res[[1]]$rescaleS
    for (i in 2:length(res))
    {
      S = cbind(S, res[[i]]$rescaleS)
    }
    dist = as.dist( 1 - cor(S, method="spearman", use="pairwise.complete.obs") )
    plot( hclust(dist) )
  }

  ### calcurate sum of correlations between distinct seeds,
  ### and select stable one
  cor_sum = myFastICA_list_correlation_totals(res)
  idx = which.max(cor_sum)
  if(verbose){cat(sprintf("#picked seeds: %d\n",seeds[idx]))}
  if(verbose){cat(sprintf("#signal component indexes: %s\n", paste(res[[idx]]$signal_components,sep=",",collapse=",") ))}
  res[[ idx ]]
}

myFastICA_estimate_n.comp_by_PCA <- function(X, cumulative_proportion, n.comp.upper_bound=100, verbose = T)
{
  pca = try( prcomp(X, scale=T) , TRUE )
  if (class(pca) == "try-error") {return(1)}
  n.comp = length( which( summary(pca)$importance[3,] < cumulative_proportion ) ) + 1
  if (verbose){
    cat(sprintf(
      "#starting number of components: %d; cummulative proportion: %f\n",
      n.comp, summary(pca)$importance[3,n.comp]
    ))
  }

  if ( n.comp > n.comp.upper_bound ){ n.comp = n.comp.upper_bound }
  if (verbose){
    cat(sprintf(
      "#number of components after upper bounding: %d; cummulative proportion: %f\n",
      n.comp, summary(pca)$importance[3,n.comp]
    ))
  }
  if ( n.comp > n.comp.upper_bound ){ n.comp = n.comp.upper_bound }

  n.comp
}

myFastICA_list_correlation_totals <- function(ica_list)
{
  res = ica_list

  ### calcurate sum of correlations between distinct seeds.
  cor_sum = c()
  for (s in 1:length(res))
  {
    cor_sum[[s]] = sum( unlist( sapply(
      1:ncol(res[[s]]$rescaleS),
      function(i)
      {
        max( unlist( sapply(
          1:length(res),
          function(t)
          {
            if (s == t) {return (0)}
            sapply(
              1:ncol(res[[t]]$rescaleS),
              function(j)
              {
                cor = cor( res[[s]]$rescaleS[,i] , res[[t]]$rescaleS[,j] ,
                           method="spearman", use="pairwise.complete.obs" )
                if (is.na(cor)){cor = 0}
                cor
              }
            )
          }
        ) ) )
      }
    ) ) )
  }
  cor_sum
}

myFastICA_from_a_seed <- function( data, n.comp, seed=11, ...)
{
  library(fastICA)
  set.seed(seed)

  #ica = fastICA( data  , n.comp )
  ica = try( fastICA( data  , n.comp ) , TRUE)
  trial = 1
  while ( (class(ica) == "try-error") & (trial < 1) ){
    trial = trial + 1
    set.seed(seed * trial)
    ica = try( fastICA( data  , n.comp ) , TRUE)
  }
  if (class(ica) == "try-error") {return(list())}

  ### adjust sign
  sign = rep(1, ncol(ica$S))
  sign[ apply( ica$A, 1, mean ) < 0 ] = -1
  ica$S = ica$S %*% diag(sign)
  ica$A = diag(sign) %*% ica$A

  ### assign col/row names
  rownames(ica$S) = rownames(ica$X)
  colnames(ica$A) = colnames(ica$X)

  ### scaling S to X scale
  ica$rescaleS = sapply(1:n.comp, function(i) apply(ica$S[,i,drop=F] %*% ica$A[i,],1,sum ) )

  ### minus factors are treated as NA
  ica$rescaleS[ ica$rescaleS < 3 ] = NA
  #ica$rescaleS[ ica$rescaleS < 3 ] = 0
  rownames(ica$rescaleS) = rownames(ica$S)
  colnames(ica$rescaleS) = colnames(ica$S)

  ### consider low constitution components as noise.
  cs = colSums(ica$rescaleS,na.rm=T)
  ica$signal_components = which( cs >= mean(cs) * 0.1 )

  ica
}


################################################
# implementation of DPI algorithm            ###
################################################

peakClustersDecomposedCtss <- function(ica, gaussian_window_size_half = 20, bedLine=matrix(NA,ncol=0,nrow=0) )
{
  buf = lapply(
    ica$signal_components,
    function(component)
    {
      res = peakClustersFromCtssVec(
        ica$rescaleS[,component] ,
        gaussian_window_size_half = gaussian_window_size_half,
        bedLine = bedLine
      )
      if ( is.null( dim(res) ) ){return(NA)}
      res[,"name"] = paste( res[,"name"] , component, sep=";comp:" )
      res
    }
  )

  idx = which( sapply(
    buf,
    function(res) is.null( dim(res) )
  ) )
  if (length(idx) == length(buf)){return(NA)}
  if (length(idx) > 0){buf = buf[-idx]}

  chrom = unlist( lapply(buf, function(b) b[,"chrom"]  ) )
  start = unlist( lapply(buf, function(b) b[,"start"]  ) )
  stop  = unlist( lapply(buf, function(b) b[,"stop" ]  ) )
  name  = unlist( lapply(buf, function(b) b[,"name" ]  ) )
  score = unlist( lapply(buf, function(b) b[,"score" ] ) )
  strand= unlist( lapply(buf, function(b) b[,"strand"] ) )
  bedTable = cbind(chrom, start, stop, name, score, strand)
  bedTable
}

peakClustersFromCtssVec <- function(ctssVec, gaussian_window_size_half = 20, bedLine = matrix(NA,ncol=0,nrow=0))
{
  bedTable = ctssVec2bedTable(ctssVec)
  res = peakClusters(bedTable, gaussian_window_size_half = gaussian_window_size_half, bedLine=bedLine)
  res
}

ctssVec2bedTable <- function(ctssVec)
{
  coords = t(sapply(
    names( ctssVec ),
    function(str)
    {
      tmp = strsplit(str,":")[[1]]
      chrom = tmp[1]
      tmp = strsplit(tmp[2],"\\.")[[1]]
      start = tmp[1]
      tmp = strsplit(tmp[3],",")[[1]]
      stop = tmp[1]
      strand = tmp[2]
      c(chrom, start, stop, strand)
    }
  ) )

  colnames(coords) = c("chrom", "start", "stop", "strand")
  coords = as.data.frame( coords , stringsAsFactors = F)
  coords$start = as.numeric(coords$start)
  coords$stop = as.numeric(coords$stop)
  score = ctssVec
  bedTable = cbind(coords, score)
  bedTable
}

peakClusters <- function(bedTable, gaussian_window_size_half = 20, bedLine = matrix(NA,ncol=0,nrow=0) )
{
  ### prep
  if ( nrow( bedTable ) == 0 ){return(NA)}
  pos = range( bedTable$start )
  vec = rep(0, pos[2] - pos[1] + 1 )
  vec[ bedTable$start - pos[1] + 1 ] = bedTable$score
  vec[ is.na(vec) ] = 0

  ### smoothing & cluster detection
  vec.padding = c( rep(0,gaussian_window_size_half), vec, rep(0, gaussian_window_size_half) )
  vec.padding.smooth = filter( vec.padding, dnorm(seq(-2,2,by=(2 / gaussian_window_size_half) )) )
  cutoff = max( median( vec.padding.smooth, na.rm=T) , 3)

  ### assigning cluster IDs
  clst = rep(0,length(vec.padding.smooth))
  clst[ vec.padding.smooth > cutoff ] = 1
  clst_id = 0
  for (i in 2:length(vec.padding.smooth))
  {
    if ( clst[i] == 0 ){next}
    if ( clst[i - 1] == 0 ){ clst_id = clst_id + 1 }
    clst[i] = clst_id
  }
  clst = clst[ (gaussian_window_size_half + 1):(gaussian_window_size_half + length(vec))]
  cluster_ids = unique( clst )
  cluster_ids = cluster_ids[ cluster_ids != 0 ]
  if (length(cluster_ids) == 0){return(NA)}

  ### prep for output
  chrom = unique( bedTable$chrom )
  strand = unique( bedTable$strand )
  chrom = chrom[ !is.na(chrom) ]
  strand = strand[ !is.na(strand) ]
  if (length(chrom)  < 1) {return(NA)}
  if (length(strand) < 1) {return(NA)}
  if (length(chrom) > 1)
    {print(sprintf("multiple chroms:%s",paste(chrom,sep=",",collapse=",")), file=stderr()) ; return(NA)}
  if (length(strand) > 1)
    {print(sprintf("multiple strands:%s",paste(strand,sep=",",collapse=",")), file=stderr()) ; return(NA)}

  ### cluster output as bed table
  clusterBedTable = t( sapply(
    cluster_ids,
    function(i)
    {
      idx = which( clst == i )
      highestPeak = idx[ which.max( vec[idx] )[1] ] + pos[1] - 1
      cluster_pos = range( idx )  + pos[1] - 1

      origin_cluster = sprintf("%s:%s..%s,%s",chrom, as.integer(pos[1]), as.integer(pos[2]+1),strand) 
      if (nrow(bedLine) != 0){
        origin_cluster = bedLine$name
      }
      name = sprintf(
        "%s;peak:%s:%s..%s,%s;peakCounts:%s",
        origin_cluster,
        chrom,as.integer(highestPeak),as.integer(highestPeak+1),strand,
        max( vec[idx] )
      )
      tmp = c(chrom, cluster_pos[1], cluster_pos[2] + 1 ,name,1000,strand )
      names(tmp) = c("chrom","start","stop","name","score","strand")
      tmp 
    }
  ) )
  clusterBedTable
}

getCtssCounts <- function(bedLine, infile_ctss, force_strandedness=T, noise_subtraction_ratio=0 )
{
  if ( nrow(bedLine)  == 0 ){return(c())}

  strand = NA
  #if ( length( grep("\\+.bw$", infile_ctss) ) > 0 ){strand = "+"}
  #if ( length( grep("\\-.bw$", infile_ctss) ) > 0 ){strand = "-"}
  if ( length( grep(".fwd.bw$", infile_ctss) ) > 0 ){strand = "+"}
  if ( length( grep(".rev.bw$", infile_ctss) ) > 0 ){strand = "-"}
  if ( is.na(strand) ) { return(c()) }

  if ( ( force_strandedness == T ) && ( strand != bedLine$strand ) )
    {return(c())}

  command = sprintf(
    "bigWigToBedGraph %s /dev/stdout -chrom=%s -start=%s -end=%s",
    infile_ctss, bedLine$chrom, bedLine$start, bedLine$stop
  )
  res = system(command, intern=T)

  if ( length(res) == 0 ){return(c())}
  res = t( sapply( res , function(str) strsplit(str,"\t")[[1]] ) )
  res = as.data.frame(res , stringsAsFactors =F)
  for (i in c(2,3,4)) { res[,i] = as.numeric(res[,i]) }
  colnames(res) = c("chrom","start","stop","score")

  idx = which( ( res$stop - res$start ) >= 2 )
  for (i in idx)
  {
    tmp = as.data.frame( cbind(res$chrom[i], res$start[i]:(res$stop[i]-1), (res$start[i]+1):res$stop[i],res$score[i]) )
    for (j in c(2,3,4)) { tmp[,j] = as.numeric(as.character(tmp[,j])) }
    colnames(tmp) = c("chrom","start","stop","score")
    res = rbind(res, tmp)
  }
  name = paste(
    res$chrom, ":",
    as.integer(res$start), "..", as.integer(res$stop), 
    ",", strand, sep=""
  )
  res = as.data.frame( cbind( res$chrom, res$start, res$stop, name, res$score, strand ) , stringsAsFactors=F )
  if (length(idx) > 0){ res = res[-idx,] }

  colnames(res) = c("chrom","start","stop","name","score","strand")
  res[,"score"] = as.numeric( res[,"score"] )
  res[,"start"] = as.integer(as.numeric( res[,"start"] ))
  res[,"stop"]  = as.integer(as.numeric( res[,"stop"] ))
  res = res[ order(res$start ) , 1:ncol(res), drop=F ]

  ### add ######
  if (noise_subtraction_ratio > 0)
  {
    res$score = res$score - ( max( res$score ) * noise_subtraction_ratio )
  }
  res = res[ res$score > 0 , ]
  ##############

  res
}

getCtssCountsTable <- function(bedLine, infiles, noise_subtraction_ratio)
{
  j = 1
  buf = getCtssCounts(bedLine, infiles[j], force_strandedness=T, noise_subtraction_ratio = noise_subtraction_ratio )
  while ( length( buf ) == 0 ){
    j = j + 1
    buf = getCtssCounts(bedLine, infiles[j], force_strandedness=T, noise_subtraction_ratio = noise_subtraction_ratio )
  }
  tmp = buf[,"score",drop=F]
  rownames(tmp) = buf[,"name"]
  colnames(tmp) = infiles[j]
  tbl = tmp

  if (length(infiles) == j){return(tbl)}

  for (i in (j+1):length(infiles)) 
  {
    buf = getCtssCounts(bedLine, infiles[i], force_strandedness=T, noise_subtraction_ratio = noise_subtraction_ratio )
    if ( length( buf ) == 0 ){ next }
    #if ( is.na(buf) ){next}
    tmp = buf[,"score",drop=F]
    rownames(tmp) = buf[,"name"]
    colnames(tmp) = infiles[i]
    rows = unique( c(rownames(tbl), rownames(tmp)) )
    cols = unique( c(colnames(tbl), colnames(tmp)) )
    tbl.new = data.frame( matrix(0, nrow=length(rows), ncol=length(cols) ) )
    rownames(tbl.new) = rows
    colnames(tbl.new) = c(colnames(tbl), infiles[i])
    tbl.new[rownames(tbl),colnames(tbl)] = tbl
    tbl.new[rownames(tmp),colnames(tmp)] = tmp
    tbl = tbl.new
  }
  tbl = tbl[ rowSums(tbl) != 0 , 1:ncol(tbl), drop=F]
  tbl = tbl[ , colSums(tbl) != 0 ,drop=F]
  tbl
}


peakClustersFromCtssVec_print <- function(ctss, gaussian_window_size_half, bedLine)
{
  ctss = rowSums(ctss)
  res = peakClustersFromCtssVec(
    ctss,
    gaussian_window_size_half = gaussian_window_size_half,
    bedLine = bedLine
  )
  if ( is.null( dim(res) ) ){
    write.table( bedLine , sep="\t", quote=F, row.names=F,col.names=F)
  }else{
    write.table( res , sep="\t", quote=F, row.names=F,col.names=F)
  }
}

    
main <- function(
  infile_base, path, pattern, exclude_prefix=NA,
  gaussian_window_size_half=20,
  n.comp.upper_bound = Inf,
  length_to_decompose = 200 ,
  noise_subtraction_ratio = 0,
  verbose=T )
{
  infile_ctss = dir(path=path, pattern=pattern,full.names=T)

  ### special file exclusion
  if (! is.na(exclude_prefix))
  {
    infile_ctss = grep( exclude_prefix ,infile_ctss,value=T,invert=T)
  }

  base = read.table(infile_base,sep="\t",as.is=T,nrow=-1)

  colnames(base) = c("chrom","start","stop","name","score","strand")

  for (i in 1:nrow(base))
  {
    bedLine = base[i,]
    if (verbose){ cat( "#", paste( bedLine , collapse="\t") , "\n") }


    if ( (bedLine$stop - bedLine$start)  < length_to_decompose )
    {
      write.table( bedLine, sep="\t", quote=F, row.names=F,col.names=F)
      next
    }
    ctss = getCtssCountsTable(bedLine, infile_ctss, noise_subtraction_ratio=noise_subtraction_ratio) 
    ica = myFastICA(ctss, verbose=verbose, n.comp.upper_bound = n.comp.upper_bound )

    if(length(ica) == 0)
    {
      peakClustersFromCtssVec_print(ctss, gaussian_window_size_half, bedLine)
    }else{
      res = peakClustersDecomposedCtss(
        ica,
        gaussian_window_size_half = gaussian_window_size_half,
        bedLine
      )
      if ( is.null( dim(res) ) )
      {
        peakClustersFromCtssVec_print(ctss, gaussian_window_size_half, bedLine)
      }else{
        write.table( res , sep="\t", quote=F, row.names=F,col.names=F)
      }
      write.table( cbind("###", rownames(ica$rescaleS) , ica$rescaleS ) , sep="\t", quote=F, row.names=F,col.names=F)
    }
  }

}



