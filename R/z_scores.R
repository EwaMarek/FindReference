# na wejście macierz ekspresji
# na wyjście macierz z-scoresów

z_scores = function(exprdata_norep){
  
  SampSD = apply(exprdata_norep, 2, sd)
  SampMean = apply(exprdata_norep, 2, mean)
  zscore = exprdata_norep
  for (i in 1:length(SampMean))  {
    zscore[,i] = (exprdata_norep[,i]-SampMean[i])/SampSD[i]
  }
  return(zscore)
}