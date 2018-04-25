# na wejście tabela macierz z-scoresów(obiekt) oraz tabelka z danymi (linia komórkowa, treatment, czas, dawka), informacja czy ma być zwracana kontrola
# na wyjście macierz fold change

fold_change = function(obiekt, uniq_S, C_columns, z_sc){

  # znalezienie kontroli i napromieniowanych prób
  Controls =  uniq_S[which(uniq_S[,'Treatment'] == 'C'),]
  IR_samples =  uniq_S[which(uniq_S[,'Treatment'] != 'C'),]

  # przygotowania macierzy na FC
  how_many_col = length(which(uniq_S[,'Treatment'] != 'C'))
  FC_data=array(0,dim=c(length(obiekt[, 1]), how_many_col))
  row.names(FC_data) = row.names(obiekt)

  #
  col = try(apply(uniq_S[which(uniq_S[,'Treatment'] != 'C' ),][c(dim(uniq_S)[2])], 1, paste, sep='', collapse = ' '), silent = TRUE)

  if(class(col) == "try-error"){ # jeśli tylko jedna to apply nie działa, stąd to co poniżej
    col = paste(uniq_S[which(uniq_S[, 'Treatment'] != 'C'),][c(dim(uniq_S)[2])], sep='', collapse = ' ')
  }

  colnames(FC_data) = col

  if(C_columns == TRUE){
    C_columns = array('0',dim=c(2, how_many_col))
    C_columns[1,] = col
  }


  # dla zwykłego IR
  FC_data = IR_FC(FC_data, obiekt, Controls, IR_samples, C_columns, z_sc)


  if(class(C_columns) != "logical"){
    return(FC_data)
  }else{
    return(FC_data[[1]])
    }


}
