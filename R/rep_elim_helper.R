# na wejście tabelka z danymi, macierz ekspresji, czy analizowane sa mikromacierze miRNA
# na wyjście macierz ekspresji bez powtórzeń oraz info o unikalnych probkach

# rozwinięcie funkcji rep_elimination -> możliwość obróbki listy tj. przypadku, gdy w danym eksperymencie użyto różnych typów mikromaciezrzy

rep_elim_helper = function(pr_data, tab){

  if(class(tab) == 'list'){
    wynik = mapply(rep_elimination, tab, pr_data)
  }else{
    wynik = rep_elimination(tab, pr_data)
  }

  return(wynik)
}
