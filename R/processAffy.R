# na wejscie: obiekt Affybatch lub ExonFeatureSet
# na wyjscie: dane po normalizacji i anotacji

processAffy = function(obiekt){

  if(class(obiekt) == 'AffyBatch'){

    # obrobka zwiazana z custom cdf
    obiekt = custom_cdf(obiekt)

    # normalizacja
    obiekt_znorm = expresso(obiekt, bgcorrect.method = 'rma', normalize.method = 'quantiles', pmcorrect.method = 'pmonly', summary.method = 'medianpolish')

    # usuwanie _at z nazw
    eksp_obiekt_znorm = exprs(obiekt_znorm)
    nazwy_at = row.names(eksp_obiekt_znorm)
    nazwy = sapply(strsplit(nazwy_at, split = '_', fixed = TRUE), function(x)(x[1]))
    row.names(eksp_obiekt_znorm) = nazwy

    # usuniecie sond kontrolnych
    wyjscie = eksp_obiekt_znorm[which(!(grepl('AFFX', row.names(eksp_obiekt_znorm)))),]

  }else if(class(obiekt) == 'ExonFeatureSet'){

    # normalizacja
    obiekt_znorm = oligo::rma(obiekt, background=TRUE, normalize=TRUE, subset=NULL, target="core")

    # adnotacja
    wyjscie = exon_annot(obiekt_znorm)

  }else{

    wyjscie = paste("This object class (", as.character(class(obiekt)), ") could not be normalized with norm_and_annot_function.", sep=' ')
    warning(paste("This object class (", as.character(class(obiekt)), ") could not be normalized with norm_and_annot_function.", sep=' '))
  }

  return(wyjscie)
}
