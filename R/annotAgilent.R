# na wejscie: informacje o sondach i kontrolach
# na wyjcie: adnotacja

annotAgilent = function(probeInfo, controls){

  if(probeInfo$SystematicName[1] == "Pro25G"){

    annot = annot_AgilentHuman1A(probeInfo, controls)

  }else{

    rows_names = probeInfo$SystematicName[which(controls==0)]
    annot = AnnotationDbi::mget(rows_names, revmap(org.Hs.egACCNUM), ifnotfound = NA)
    annot = lapply(annot,`[[`, 1)

  }

    return(annot)
}
