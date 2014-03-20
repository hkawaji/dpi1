

###
### usage:
###   source("myFastICA.R")
###   ica.result = myFastICA(data.matrix)
###

myFastICA <- function( data, cumulative_proportion = 0.95, seeds=11:13, hclust_plot = F, verbose = T, n.comp.upper_bound = Inf, ...)
{
  if ( ncol(data) == 1 ){return(list())}

  ### decide starting number of components
  n.comp = myFastICA_estimate_n.comp_by_PCA(data, cumulative_proportion, n.comp.upper_bound, verbose)
  if ( n.comp == 1 ){return(list())}

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
  rownames(ica$rescaleS) = rownames(ica$S)
  colnames(ica$rescaleS) = colnames(ica$S)

  ### consider low constitution components as noise.
  cs = colSums(ica$rescaleS,na.rm=T)
  ica$signal_components = which( cs >= mean(cs) * 0.1 )

  ica
}

