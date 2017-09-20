# na wejście obiekt klasy marrayRaw, znormalizowane wartości dla kanału czerwonego i zielonego
# na wyjście znormalizowane wartości dla obu kanałów z adnotacją

annot_AgilentHuman1A = function(probeInfo, controls){
  
  rows_names = probeInfo$ProbeName[which(controls==0)]
  
  annot = mapIds(hgug4110b.db, keys = rows_names, column = "ENTREZID", keytype = "PROBEID", multiVals = "CharacterList")
  annot = lapply(annot, `[[`, 1)
  
  return(annot)
}