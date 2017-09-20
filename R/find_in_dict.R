# na wejscie slownik i poszukiwany string
# na wyjscie przypisany danemu stringowi id MIMAT


find_in_dict = function(my_dict, what_is_it){
  
  FOUND_IT=FALSE
  i=1
  thats_it = what_is_it # jesli nic nie znajdzie to zwraca wartosc wejsciowa
  
  while(FOUND_IT==FALSE && i <= dim(my_dict)[1]){
    
    if(what_is_it %in% my_dict[i,2]){
    
      thats_it = my_dict[i,1]
      FOUND_IT = TRUE
      
    }else{
      
      i=i+1
      
    }
  }
  
  return(thats_it)
}