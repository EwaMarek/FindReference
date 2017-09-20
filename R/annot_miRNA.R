# na wejscie macierz ekspresji ze zla adnotacja
# na wyjsciu macierz ekspresji z adnotacja w formie MIMAT

annot_miRNA = function(exp_data){

  MIMAT_dict = system.file("inst/extdata", "MIMAT_dict.rds", package = "FindReference")
  MIMAT_dict = readRDS(MIMAT_dict)
  bad_names = rownames(exp_data)
  bad_names = gsub('\\_v.*', '', bad_names, fixed = FALSE) # usuniecie przyrostkow _vNR

  good_names = bad_names

  for (i in 1:length(bad_names)){
    gn = as.character(MIMAT_dict$accession[which(MIMAT_dict$ID == bad_names[i])])
    if(length(gn)==1){
      good_names[i] = gn
    }
  }

  rownames(exp_data) = good_names
  exp_data = exp_data[grep('MIMAT', good_names), ] # usuniecie wierszy z MI00000

  return(exp_data)

}
