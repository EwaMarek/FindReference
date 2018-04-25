# na wejście pusta macierz na FC, macierz z-scoresów(tu:obiekt), dane dla kontroli i próbek napromieniowanych (linia komórkowa, treatment, czas, dawka),
# informacja czy FC_data ma być policzony dla z-scorsów
# macierz z fold change


IR_FC = function(FC_data, obiekt, Controls, IR_samples, C_columns=FALSE, z_sc){


  if(dim(Controls)[1] > 1){ # jeśli więcej niż jedna kontrola

    for (i in 1:dim(IR_samples)[1]) {

      # w tej pętli różnice w ifach polegają na innym wyznaczaniu C_for_IR

      #############################################################################
      ################### One control for each IR sample ##########################
      #############################################################################

      if(dim(IR_samples)[1] == dim(Controls)[1]){ # jeśli dla każdego IR kontrola, to musi sie zgadzac linia i czas

        TFVector = vector(mode = 'logical', dim(Controls)[1])
        for(aj in 1:dim(Controls)[1]){
          TFVector[aj] = find_rows(Controls[aj, c("CellLine", "Time")], IR_samples[i, c("CellLine", "Time")])
        }


        C_for_IR = which(TFVector == TRUE)


        if(length(C_for_IR) == 0){ # jesli czas kontroli rowny 0, a IR inny

          C_for_IR = which(Controls[, "CellLine"] == IR_samples[i, "CellLine"])
        }

        #############################################################################
        ################### One control for each CellLine ###########################
        #############################################################################

      }else if(dim(Controls)[2] == length(unique(IR_samples[,'CellLine']))){ # jeśli jedna kontrola dla każdej linii komórkowej

        C_for_IR = which(Controls[, "CellLine"] == IR_samples[i, "CellLine"])

        if(length(C_for_IR)>1){ # jeli oprócz kontroli dla każdego IRa, wystepuja kontrole o zerowym czasie i dawce
                                # dla fold change bezużyteczna
          if(dim(IR_samples)[1] == dim(Controls)[1]){

            TFVector = vector(mode = 'logical', dim(Controls)[1])
            for(aj in 1:dim(Controls)[1]){
              TFVector[aj] = find_rows(Controls[aj, c("CellLine", "Time")], IR_samples[i, c("CellLine", "Time")])
            }

            C_for_IR = which(TFVector == TRUE)

          }

        }
        #############################################################################
        ############## Controls just for some specific times ########################
        #############################################################################

      }else{ # jeśli inna kombinacja (w założeniu kontrola dla WYBRANYCH chwil czasowych )
        C_for_IR = which(Controls[,"CellLine"] == IR_samples[i, "CellLine"])
        time_diff = abs(as.numeric(IR_samples[i, 'Time'])-as.numeric(Controls[C_for_IR,'Time']))
        better_C = which(time_diff == min(time_diff))


        if(length(better_C)>1){ # jesli pechowo dwie identyczne różnice czasowe
          better_C = better_C[1]
        }

        C_for_IR = C_for_IR[better_C] # wybranie kontroli o czasie najbliższym próbce IR

      }

      col_name_of_C = Controls[C_for_IR, ]$internalId
      col_name_of_IR = IR_samples[i, ]$internalId


      if(z_sc==TRUE){
        FC_data[, col_name_of_IR] = obiekt[,col_name_of_IR]-obiekt[,col_name_of_C]
      }else{
        FC_data[, col_name_of_IR] = (obiekt[, col_name_of_IR])/(obiekt[,col_name_of_C])
      }

      if(class(C_columns) != "logical"){
      C_columns[2,which(C_columns[1,]==col_name_of_IR)] = col_name_of_C
      }

    }

    #############################################################################
    ######################### Only one control sample ###########################
    #############################################################################

  }else if(dim(Controls)[1] == 1){ # gdy jedna kontrola dla wszystkich

    if(dim(IR_samples)[1] > 1){ # jeśli więcej niż jedna próbka IR

      for (i in 1:dim(IR_samples)[1]) {

        col_name_of_C = Controls$internalId #paste(Controls, sep='', collapse = ' ')
        col_name_of_IR = IR_samples[i, ]$internalId #paste(IR_samples[i,], sep='', collapse = ' ')

        if(z_sc==TRUE){
          FC_data[, col_name_of_IR] = obiekt[,col_name_of_IR]-obiekt[,col_name_of_C]
        }else{
          FC_data[, col_name_of_IR] = (obiekt[, col_name_of_IR])/(obiekt[,col_name_of_C])
        }

        if(class(C_columns) != "logical"){
          C_columns[2,which(C_columns[1,]==col_name_of_IR)] = col_name_of_C
        }
      }

    }else if(dim(IR_samples)[1] == 1){ # jeśli tylko jedna próbka IR

      col_name_of_C = Controls$internalId #paste(Controls, sep='', collapse = ' ')
      col_name_of_IR = IR_samples$internalId #paste(IR_samples, sep='', collapse = ' ')

      if(z_sc==TRUE){
        FC_data[, col_name_of_IR] = obiekt[,col_name_of_IR]-obiekt[,col_name_of_C]
      }else{
        FC_data[, col_name_of_IR] = (obiekt[, col_name_of_IR])/(obiekt[,col_name_of_C])
      }

      if(class(C_columns) != "logical"){
        C_columns[2,which(C_columns[1,]==col_name_of_IR)] = col_name_of_C
      }

    }

  }else if(length(Controls) == 0){print("There is no control for IR")}

return(list(FC_data, C_columns))

}
