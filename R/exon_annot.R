# na wejście obiekt typu ExonFeatureSet ze znormalizowanymi wartościami ekspresji
# na wyjście macierz ekspresji z Entrez ID

exon_annot = function(obiekt_norm){
  
  ekspr_obiektu_norm = exprs(obiekt_norm)
  nazwy = row.names(ekspr_obiektu_norm)
  
  # znalezienie accesion number
  ACC_num = as.list(huex10sttranscriptclusterACCNUM[nazwy])
  ACC_num = lapply(ACC_num, function(x) x[1]) # tylko pierwszy numer -> wydaje się być najpopularniejszy
  ACC_num=unlist(ACC_num)
  
  # znalezienie Enrez ID na bazie accesion number
  ENT_num = mapIds(org.Hs.eg.db, keys = ACC_num, column = "ENTREZID", keytype = "ACCNUM", multiVals = "first")
  #ENT_num = lapply(ENT_num, function(x) try(x[1])) # entrez id w tej kolejnosci co początkowe
  annot = lapply(ENT_num, function(x) if(length(x)==0){x='NA'}else{x=x}) # podmiana pustych characterow na NA
  annot = unlist(annot, use.names=FALSE)
  
  rownames(ekspr_obiektu_norm) = annot #make.names(annot, unique=FALSE) # make.names dodaje X na początek, bo nazwy powinny zaczynać się od litery
  ekspr_obiektu_norm = ekspr_obiektu_norm[which(!(row.names(ekspr_obiektu_norm)=="NA")),] # usunięcie wierszy bez Entrez ID
  # names_with_X = row.names(ekspr_obiektu_norm)
  # names_without_X = gsub('X', "", names_with_X) # usunięcie X-ów z nazw
  # rownames(ekspr_obiektu_norm) = names_without_X
  
  # usrednienie powtorzen
  eksp_out = as.matrix(aggregate.Matrix(ekspr_obiektu_norm, groupings=row.names(ekspr_obiektu_norm), fun='mean'))
  
  return(eksp_out)
}


