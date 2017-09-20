# na wejscie: obiekt klasy uRNAList
# na wyjscie: tablica ekspresji po adnotacji

processmiRNA = function(loaded_data){

  # normalization and extracting expression table
  norm_data = rmaMicroRna(loaded_data, normalize = TRUE, background = TRUE)
  exp_data = norm_data$TGS
  row.names(exp_data) = norm_data$genes$GeneName
  colnames(exp_data) = norm_data$targets[,'FileName']

  # remove virus and other control probes
  exp_data = exp_data[norm_data$genes$ControlType == 0, ]
  exp_data = exp_data[grep('hsa', rownames(exp_data)),]

  # MIMAT annotation
  exp_data = annot_miRNA(exp_data)

  return(exp_data)
}
