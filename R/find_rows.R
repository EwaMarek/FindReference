# na wejście wektor z macierzy, wektor do porównania
# na wyjście jeśli identyczne 1, inaczej 0

find_rows = function(a, b){
  
  if(sum(a==b)==length(b)){c=TRUE}
  
  else{c=FALSE}
  
  return(c)
  
}