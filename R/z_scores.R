# na wejście macierz ekspresji
# na wyjście macierz z-scoresów

z_scores = function(exprdata_norep){

  SampSD = apply(exprdata_norep, 2, sd, na.rm = TRUE)
  SampMean = apply(exprdata_norep, 2, mean, na.rm = TRUE)
  zscore = exprdata_norep
  for (i in 1:length(SampMean))  {
    zscore[,i] = (exprdata_norep[,i]-SampMean[i])/SampSD[i]
  }
  return(zscore)
}
