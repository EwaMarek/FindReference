# na wejście macierz z lFC oraz macierz z wartościami ekspresji kontroli dla każdej próby IR z lFC oraz info o tym, która C do którego IR
# na wyjście macierz z wartościami wskaźnika dla poszczególnych prób IR w danym eksperymencie

gene_exp_score = function(FC_data, C_expression, C_for_IR){

  if(class(FC_data) == 'matrix'){


    if(dim(FC_data)[2] != 1){

      exp_score = FC_data
      wh_IR = colnames(FC_data)

      for(i in 1:length(wh_IR)){
        wh_C = C_for_IR[2, which(C_for_IR[1,] == wh_IR[i])[1]]
        exp_score[, wh_IR[i]] = ((C_expression[, wh_C])^3)/(abs(FC_data[, wh_IR[i]])) # tu
      }
    }else{

      exp_score = ((C_expression)^3)/(abs(FC_data)) # tu

    }

  }else if(class(FC_data) == 'list'){

    exp_score = FC_data

    for(j in 1:length(FC_data)){

      wh_IR = colnames(FC_data[[j]])

      for(i in 1:length(wh_IR)){
        wh_C = C_for_IR[[j]][2, which(C_for_IR[[j]][1,] == wh_IR[i])[1]]
        exp_score[[j]][, wh_IR[i]] = ((C_expression[[j]][, wh_C])^3)/(abs(FC_data[[j]][, wh_IR[i]])) # i tu
      }
    }
  }else{
    exp_score = ':('
  }
  return(exp_score)
}
