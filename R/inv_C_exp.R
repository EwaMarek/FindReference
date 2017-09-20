# na wejście macierz ekspresji bez powtórzeń oraz tabela z funkcji log_fold_change (dotycząca info, która kontrola dla której IR)
# na wyjście tabela z odwrotnościami ekspresji prób kontrolnych, kolejno dla każdego IR-a

inv_C_exp = function(no_rep, C_IR_info){
  
  if(class(C_IR_info) == 'matrix'){
    
    C_exp = no_rep
    CIR = C_IR_info
    col_for_exp = CIR[2, ]
    
    C_exp = C_exp[, col_for_exp]
    C_expression = 1/C_exp
    
  }else if(class(C_IR_info) == 'list'){
    
    C_expression = rep(list(list()), length(C_IR_info))
    
    for (j in 1:length(C_IR_info)) {
      
      C_exp = no_rep[[j]]
      CIR = C_IR_info[[j]]
      col_for_exp = CIR[2, ]
      
      C_exp = C_exp[, col_for_exp]
      C_expression[[j]] = 1/C_exp
    }
    
  }else{
    
    C_expression = ':('
  }
  
  return(C_expression)
}