# na wejście wektor z którego mają być usunięte kwantyle, oraz wektor z wartościami tychże kwantyli
# na wyjście wektor bez wartości poniżej i powyżej wybranych kwantyli, w ich miejsce NA
# na miejscu NA zostaje NA


remove_quantiles = function(score_vector, q_values){
  
  without_q_vector = score_vector
  
  
  if(length(score_vector)>4){
    
    qA_qB = quantile(score_vector, probs = q_values, na.rm=TRUE)
  
    for (i in 1:length(score_vector)) {
      
      if(is.na(score_vector[i]) == TRUE){
        without_q_vector[i] = score_vector[i]
      }else if(score_vector[i]>qA_qB[1] & score_vector[i]<qA_qB[2]){
        without_q_vector[i] = score_vector[i]
      }else{
        without_q_vector[i] = NA
      }
      
    }
  }
  
  return(without_q_vector)
}