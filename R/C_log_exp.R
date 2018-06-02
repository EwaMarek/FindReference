# na wejście z-scores oraz info o tym, która kontrola zostaął użyta do którego eksperymentu
# na wyjście z-score jedynie kontroli

C_log_exp = function(data_z, C_for_IR){
  
  if(class(data_z) == 'matrix'){
    
    wh_samples = which(C_for_IR[1,] %in% colnames(data_z))
    C_exp = data_z[, C_for_IR[2, wh_samples]]
    
  }else if(class(data_z) == 'list'){
    
    C_exp = rep(list(list()), length(data_z))
    
    for (j in 1:length(data_z)) {
      
      wh_samples = which(C_for_IR[[j]][1,] %in% colnames(data_z[[j]]))
      C_exp[[j]] = data_z[[j]][, C_for_IR[[j]][2, wh_samples]]
    }
  }
  
  return(C_exp)
}