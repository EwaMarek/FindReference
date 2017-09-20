# na wejście tabelka z danymi, macierz ekspresji, czy analizowane miRNA
# na wyjście macierz ekspresji bez powtórzeń


rep_elimination = function(data_for_one_exp, obiekt){

  ############################################################################################
  ################# Check if argument data_for_one_exp is prepared correctly #################
  ############################################################################################
  one_colour = c("Experiment", "Cell line", "Treatment", "Time", "Dose")
  dbl_colour = c("Experiment", "Cell line", "Treatment", "Time", "Dose", "Label")
  hmc = dim(data_for_one_exp)[2]

  if(hmc == 5 && (sum(colnames(data_for_one_exp) == one_colour) != hmc)) {

    stop("Column names in table must be 'Experiment', 'Cell line', 'Time', 'Dose' exactly in THIS order for one colour experiments.")

  }else if(hmc == 6 && (sum(colnames(data_for_one_exp) == dbl_colour) != hmc)){

    stop("Column names in table must be 'Experiment', 'Cell line', 'Time', 'Dose', 'Label' exactly in THIS order for double colour experiments.")
  }


  ############################################################################################
  ############################### One colour matrix ##########################################
  ############################################################################################

  if(!("Label" %in% colnames(data_for_one_exp)) & class(data_for_one_exp) == "matrix"){


    ### check if the colnames in expression matrix are the same as rownames in table
    if(sum(colnames(obiekt) %in% rownames(data_for_one_exp))<dim(obiekt)[2]){
      stop(paste("Colnames in expression matrix aren't the same as in table with info about microarrays. The problem occured with columns in expression matrix with index:",
                 paste(as.character(which((colnames(obiekt) %in% rownames(data_for_one_exp)) == FALSE)), sep='', collapse=', ') ,sep = ' '))
    }

    ### prepare matrix for results
    uniq_Samples = unique(data_for_one_exp)
    exprdata_norep=array(0,dim=c(dim(obiekt)[1],dim(uniq_Samples)[1]))
    row.names(exprdata_norep) = row.names(obiekt)
    col_names = apply(uniq_Samples, 1, paste, collapse=" ")
    colnames(exprdata_norep) = col_names

    ### calculate mean value for repetitions
    for (i in 1:length(uniq_Samples[,1])) {
      rep_microarrays = apply(data_for_one_exp, 1, find_rows, b=uniq_Samples[i,]) # find repetitions
      rep_micro_names = rownames(data_for_one_exp)[rep_microarrays]

      if(length(rep_micro_names)>1){
        exprdata_norep[, paste(uniq_Samples[i,], sep='', collapse = ' ')] = rowMeans(obiekt[,which(colnames(obiekt) %in% rep_micro_names)])
      }else{
        exprdata_norep[, paste(uniq_Samples[i,], sep='', collapse = ' ')] = obiekt[,which(colnames(obiekt) %in% rep_micro_names)]
      }
    }

  ############################################################################################
  ############################## Double colour matrix ########################################
  ############################################################################################

  }else if("Label" %in% colnames(data_for_one_exp)){

    ### check if the colnames in both expression matrices are the same as rownames in table
    if(sum(colnames(obiekt[[1]]) %in% rownames(data_for_one_exp))<dim(obiekt[[1]])[2]){
      stop(paste("Colnames in expression matrix aren't the same as in table with info about microarrays. The problem occured with columns in expression matrix with index:",
                 paste(as.character(which((colnames(obiekt[[1]]) %in% rownames(data_for_one_exp)) == FALSE)), sep='', collapse=', ') ,sep = ' '))
    }else if(sum(colnames(obiekt[[2]]) %in% rownames(data_for_one_exp))<dim(obiekt[[2]])[2]){
      stop(paste("Colnames in expression matrix aren't the same as in table with info about microarrays. The problem occured with columns in expression matrix with index:",
                 paste(as.character(which((colnames(obiekt[[2]]) %in% rownames(data_for_one_exp)) == FALSE)), sep='', collapse=', ') ,sep = ' '))
    }

      ### find unique rows in table
      lab = which(colnames(data_for_one_exp) == 'Label')
      uniq_Samples = unique(data_for_one_exp[, -lab])

      ### prepare matrix for results
      the_Cy3 = obiekt[[1]]
      the_Cy5 = obiekt[[2]]
      exprdata_norep = array(0,dim=c(dim(the_Cy3)[1],dim(uniq_Samples)[1]))
      col_names = apply(uniq_Samples, 1, paste, collapse=" ")
      colnames(exprdata_norep) = col_names
      rownames(exprdata_norep) = rownames(the_Cy3)

      for (i in 1:length(uniq_Samples[,1])) {

            ### find repetitions indexes
            rep_microarrays = apply(data_for_one_exp[,-lab], 1, find_rows, b=uniq_Samples[i,])

            ### find microarray names for this repetition and its colour label
            rep_micro_names = rownames(data_for_one_exp)[rep_microarrays]
            dye = data_for_one_exp[rep_microarrays, lab]
            dyetype = unique(dye)

            # znalezienie, które kolumny w przetworzonych danych nazywają się jak powyżej,
            # ale TU bardzo istotna kolejność (żeby poprawnie zidentyfikowac barwnik) -> stąd pętla
            # wtf, dye w kolejnosci jak rep microarrays, bez petli nie wiadomo ktora kolumna odpowiada danemu powtorzeniu?
            wh_col = array(0,dim=c(length(rep_micro_names),1))
            for (z in 1:length(rep_micro_names)) {
              wh_col[z] = which(colnames(the_Cy3) %in% rep_micro_names[z])
              }


            if(length(rep_micro_names)>1){   # jeśli są powtórzenia

              # kolumna pobierana z macierzy ekspresji dla danego barwnika
              for(k in 1:length(dye)){

                if(dye[k] == dyetype[1]){   p =the_Cy3[,wh_col[k]]
                }else if(dye[k] == dyetype[2]){  p =the_Cy5[,wh_col[k]] }

                # przy pierwszym obiegu pętli przypisanie kolumny do zmiennej (poczatek macierzy powtorzeń)
                if(k==1){ rep_array = p
                }else{  rep_array = cbind(rep_array, p) }    # powstała rep_array to macierz z wartosciami ekspresji dla wszystkich powtórzeń
             }

              # uśrednienie
              exprdata_norep[, paste(uniq_Samples[i,], sep='', collapse = ' ')] = rowMeans(rep_array)

            }else{ #jeśli nie ma powtórzeń

              if(dye == dyetype[1]){   exprdata_norep[, paste(uniq_Samples[i,], sep='', collapse = ' ')] = the_Cy3[,wh_col]

              }else if(dye == dyetype[2]){  exprdata_norep[, paste(uniq_Samples[i,], sep='', collapse = ' ')] = the_Cy5[,wh_col] }
            }
      }
  }

 return(list(exprdata_norep, uniq_Samples))
}
