
source("scripts/myFastICA.R")


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

  if (noise_subtraction_ratio > 0)
  {
    res$score = res$score - ( max( res$score ) * noise_subtraction_ratio )
  }
  res = res[ res$score > 0 , ]
  res
}

getCtssCountsTable_old <- function(bedLine, infiles, noise_subtraction_ratio)
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

getCtssCountsTable <- function(bedLine, infiles, noise_subtraction_ratio = 0)
{
  # prep
  tbl = c()
  tmpdir = system("mktemp -d", intern=T)
  tmplist = sprintf("%s/list.txt", tmpdir)
  tmpmat  = sprintf("%s/mat.txt", tmpdir)
  command = sprintf(
    "./scripts/getCtssCountsTable.sh -c %s -s %s -e %s -i %s > %s",
    bedLine$chrom, bedLine$start, bedLine$stop, tmplist, tmpmat, tmpmat
  )
  if ( bedLine$strand == "+" ){
    infiles = infiles[ grep(".fwd.bw$", infiles) ]
  } else if (bedLine$strand == "-" ) {
    infiles = infiles[ grep(".rev.bw$", infiles) ]
  } else { 
    print( sprintf("ERROR: wrong strand:%s", bedLine$strand)  )
    q()
  }

  # main
  tryCatch(
    { 
      write.table(infiles, quote=F, row.names=F, col.names=F, file=tmplist )
      system(command, intern=T)
      tbl = read.table(tmpmat, row.names=1,header=T, check.names=F)
      rownames(tbl) = paste( rownames(tbl) , bedLine$strand, sep=",")
      system(sprintf("rm -rf %s", tmpdir))
    },
    finally = system(sprintf("rm -rf %s", tmpdir))
  )

  # post processing
  sbtrct = apply(tbl,2,max) * noise_subtraction_ratio
  tbl = t( t(tbl) - sbtrct )
  tbl[ tbl < 0 ] = 0
  idx1 = which( apply( tbl , 1, max) > 0 )
  idx2 = which( apply( tbl , 2, max) > 0 )
  tbl = tbl[idx1,idx2,drop=F]
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


main_spi <- function(
  infile_base,
  infile_ctss_pref,
  gaussian_window_size_half,
  length_to_decompose,
  noise_subtraction_ratio,
  verbose=T )
{
  base = read.table(infile_base,sep="\t",as.is=T,nrow=-1)
  colnames(base) = c("chrom","start","stop","name","score","strand")

  infile_ctss = c(
    sprintf("%s.fwd.bw", infile_ctss_pref ),
    sprintf("%s.rev.bw", infile_ctss_pref )
  )

  for (i in 1:nrow(base))
  {
    bedLine = base[i,]
    if ( (bedLine$stop - bedLine$start)  < length_to_decompose )
    {
      write.table( bedLine, sep="\t", quote=F, row.names=F,col.names=F)
      next
    }
    if (verbose){ cat( "#", paste( bedLine , collapse="\t") , "\n") }
    ctss = getCtssCountsTable(bedLine, infile_ctss, noise_subtraction_ratio )
    peakClustersFromCtssVec_print(ctss, gaussian_window_size_half, bedLine)
  }
}


main_dpi <- function(
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
      write.table( ica$rescaleS , sep="\t", quote=F, row.names=T,col.names=F)
    }
  }

}



