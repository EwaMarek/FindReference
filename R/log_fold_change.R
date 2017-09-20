# na wejście macierz z-scoresów, tabela z informacjami o unikalnych próbkach oraz parametr określający czy ma być zwracana informacja
# o tym, która kontrola została zastosowana do policzenia lFC dla danej próbki IR (TRUE/FALSE), czy liczymy dla z-scoresów

# na wyjście macierz z policzonymi logFC oraz informacja o tym, która kontrola zostaął użyta dla danej próbki IR

log_fold_change = function(z_sc, uniq_S, add_info, z_score=TRUE){
  
  if(class(z_sc) != 'list'){
    FC_dataset = fold_change(z_sc, uniq_S, add_info, z_score)
    FC_data = FC_dataset[[1]]
    C_for_IR = FC_dataset[[2]]

  }else if(class(z_sc) == 'list'){

    FC_data = rep(list(list()), length(z_sc))
    C_for_IR = rep(list(list()), length(z_sc))

    for (j in 1:length(z_sc)) {
      FC_dataset = fold_change(z_sc[[j]], uniq_S[[j]], add_info, z_score)
      FC_data[[j]] = FC_dataset[[1]]
      C_for_IR[[j]] = FC_dataset[[2]]
    }
  }
  
  return(list(FC_data, C_for_IR))
}