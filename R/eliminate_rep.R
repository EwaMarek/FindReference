# na wejscie macierz ekspresji (badz lista - eksperymenty dwukolorowe), tabela dla wszystkich eksperymentow
# na wyjscie macierz ekspresji z usrednionymi powtorzeniami oraz data.frame z unikalnymi probkami w danym eksperymencie

#' @title Elimination of replications within an experiment
#'
#' @param expMat An expression matrix with normalized expression values or a list of two for double coloured experiments.
#'
#' @param eid A character indicating experiment id in ExpInfoTable.
#'
#' @param  ExpInfoTable A data.frame with information about all of the experiments. The column names must be 'Experiment', 'SampleID', 'CellLine', 'Treatment', 'Time', 'Dose'
#' 'Platform', 'Label' exactly in this order.
#'
#' @param sdrfFiles A list of sdrf files content.
#'
#' @return Function returns list with two elements. The first one is an expression matrix without replications with changed column names.
#' The second one is a data.frame with information about UNIQUE samples and is needed to create stability ranking.
#'
#' @description
#' \code{eliminate_rep} function averages expression values for replications within an experiment and prepare expression matrices for being used
#' in \code{create_ranking} function. It also returns information about unique samples in each experiment what is also needed for
#' \code{create_ranking} function to run.
#'
#' If there are no replications in the data \code{eliminate_rep} function will just prepared objects to use with \code{create_ranking} function.
#'
#' @details
#' Example table structure (third argument) should look like this:
#'
#' \tabular{cccccc}{
#'  Experiment  \tab SampleID \tab CellLine   \tab  Treatment   \tab   Time   \tab   Dose   \tab Platform  \tab Label  \cr
#'  Exp1        \tab    S1    \tab  HCT116    \tab      C       \tab     0    \tab     0    \tab  Agilent  \tab  Cy3   \cr
#'  Exp1        \tab    S2    \tab  HCT116    \tab      C       \tab     0    \tab     0    \tab  Agilent  \tab  Cy5   \cr
#'  Exp1        \tab    S3    \tab  HCT116    \tab      IR      \tab     1    \tab     4    \tab  Agilent  \tab  Cy3   \cr
#'  Exp1        \tab    S4    \tab  HCT116    \tab      IR      \tab     2    \tab     4    \tab  Agilent  \tab  Cy5   \cr
#'  Exp1        \tab    S5    \tab  HCT116    \tab      IR      \tab     4    \tab     2    \tab  Agilent  \tab  Cy3   \cr
#'  }
#'
#'  Note that in treatment column it should be C for control arrays and IR for treated ones (even if the samples were treated with
#'  for example chemicals and not with ionizing radiation).
#'
#'  If you have downloaded data from ArrayExpress database you can find information needed for creating above in .sdrf files.
#'
#' @seealso
#' \code{\link{create_ranking}}
#'
#'
#' @rdname eliminate_rep
#'
#' @export

eliminate_rep = function(expMat, eid, ExpInfoTable, sdrfFiles){

  if(class(expMat) != 'character'){

    ### check if experiment is single or double colour
    dbl_col = FALSE
    dbl_exp = FALSE

    if(class(expMat) == 'list'){

      dyes = ExpInfoTable[which(ExpInfoTable$Experiment == eid), 'Label']

      if(length(unique(dyes)) == 2){
        dbl_col = TRUE
      }else{
        dbl_exp = TRUE
      }
    }

    # if the experiment is single coloured
    if(dbl_col == FALSE){

      if(dbl_exp == FALSE){
        ### change the colnames so they will be the same as in ExpInfoTable
        act_sdrf = sdrfFiles[eid]
        if(!(is.null(unlist(act_sdrf[[1]])))){
          new_cols = act_sdrf[[1]]$Source.Name[match(gsub("./", "", colnames(expMat)), act_sdrf[[1]]$Array.Data.File)]
          if(sum(is.na(new_cols))>1){
            new_cols = act_sdrf[[1]]$Source.Name[match(paste0(colnames(expMat), ".txt"), act_sdrf[[1]]$Array.Data.File)]
          }
          
          colnames(expMat) = new_cols
        }
  
        ### prepare matrix for results
        uniq_Samples = unique(ExpInfoTable[which(ExpInfoTable$SampleID %in% colnames(expMat) & ExpInfoTable$Experiment == eid), -2])
        Svector = paste0("S", as.character(1:dim(uniq_Samples)[1]))
        uniq_Samples$internalId = Svector
        exprdata_norep=array(0,dim=c(dim(expMat)[1],dim(uniq_Samples)[1]))
        row.names(exprdata_norep) = row.names(expMat)
        colnames(exprdata_norep) = Svector
  
        ### calculate mean value
        for (i in 1:dim(uniq_Samples)[1]) {
  
          # find rows in ExpInfoTable that are repetetive
          uniIDS = vector(mode='logical', length = dim(ExpInfoTable)[1])
          for(j in 1:dim(ExpInfoTable)[1]){
            uniIDS[j] = find_rows(ExpInfoTable[j, c(1, 3:6)], uniq_Samples[i, c(1:5)])
          }
  
          # find repetetive sample names
          rep_micro_names = colnames(expMat)[which(colnames(expMat) %in% ExpInfoTable$SampleID[uniIDS])]
  
          # calculate mean value for repetetive samples
          if(length(rep_micro_names)>1){
            exprdata_norep[, uniq_Samples$internalId[i]] = rowMeans(expMat[,which(colnames(expMat) %in% rep_micro_names)])
          }else{
            exprdata_norep[, uniq_Samples$internalId[i]] = expMat[,which(colnames(expMat) %in% rep_micro_names)]
          }
        }
      }else{
        
        ### change the colnames so they will be the same as in ExpInfoTable
        act_sdrf = sdrfFiles[eid]
        new_cols1 = act_sdrf[[1]]$Source.Name[match(gsub("./", "", colnames(expMat[[1]])), act_sdrf[[1]]$Array.Data.File)]
        if(sum(is.na(new_cols1))>1){
          new_cols1 = act_sdrf[[1]]$Source.Name[match(paste0(colnames(expMat[[1]]), ".txt"), act_sdrf[[1]]$Array.Data.File)]
        }
        
        new_cols2 = act_sdrf[[1]]$Source.Name[match(gsub("./", "", colnames(expMat[[2]])), act_sdrf[[1]]$Array.Data.File)]
        if(sum(is.na(new_cols2))>1){
          new_cols2 = act_sdrf[[1]]$Source.Name[match(paste0(colnames(expMat[[2]]), ".txt"), act_sdrf[[1]]$Array.Data.File)]
        }
        
        colnames(expMat[[1]]) = new_cols1
        colnames(expMat[[2]]) = new_cols2
        
        ### prepare matrices for results
        uniq_Samples1 = unique(ExpInfoTable[which(ExpInfoTable$SampleID %in% colnames(expMat[[1]]) & ExpInfoTable$Experiment == eid), -2])
        Svector1 = paste0("S", as.character(1:dim(uniq_Samples1)[1]), "A")
        uniq_Samples1$internalId = Svector1
        exprdata_norep1=array(0,dim=c(dim(expMat[[1]])[1],dim(uniq_Samples1)[1]))
        row.names(exprdata_norep1) = row.names(expMat[[1]])
        colnames(exprdata_norep1) = Svector1

        uniq_Samples2 = unique(ExpInfoTable[which(ExpInfoTable$SampleID %in% colnames(expMat[[2]]) & ExpInfoTable$Experiment == eid), -2])
        Svector2 = paste0("S", as.character(1:dim(uniq_Samples2)[1]), "B")
        uniq_Samples2$internalId = Svector2
        exprdata_norep2=array(0,dim=c(dim(expMat[[2]])[1],dim(uniq_Samples2)[1]))
        row.names(exprdata_norep2) = row.names(expMat[[2]])
        colnames(exprdata_norep2) = Svector2
      
        ### calculate mean value
        for (i in 1:dim(uniq_Samples1)[1]) {

          # find rows in ExpInfoTable that are repetetive
          uniIDS = vector(mode='logical', length = dim(ExpInfoTable)[1])
          for(j in 1:dim(ExpInfoTable)[1]){
            uniIDS[j] = find_rows(ExpInfoTable[j, c(1, 3:6)], uniq_Samples1[i, c(1:5)])
          }

          # find repetetive sample names
          rep_micro_names = colnames(expMat[[1]])[which(colnames(expMat[[1]]) %in% ExpInfoTable$SampleID[uniIDS])]

          # calculate mean value for repetetive samples
          if(length(rep_micro_names)>1){
            exprdata_norep1[, uniq_Samples1$internalId[i]] = rowMeans(expMat[[1]][,which(colnames(expMat[[1]]) %in% rep_micro_names)])
          }else{
            exprdata_norep1[, uniq_Samples1$internalId[i]] = expMat[[1]][,which(colnames(expMat[[1]]) %in% rep_micro_names)]
          }
        }


        for (i in 1:dim(uniq_Samples2)[1]) {

          # find rows in ExpInfoTable that are repetetive
          uniIDS = vector(mode='logical', length = dim(ExpInfoTable)[1])
          for(j in 1:dim(ExpInfoTable)[1]){
            uniIDS[j] = find_rows(ExpInfoTable[j, c(1, 3:6)], uniq_Samples2[i, c(1:5)])
          }

          # find repetetive sample names
          rep_micro_names = colnames(expMat[[2]])[which(colnames(expMat[[2]]) %in% ExpInfoTable$SampleID[uniIDS])]

          # calculate mean value for repetetive samples
          if(length(rep_micro_names)>1){
            exprdata_norep2[, uniq_Samples2$internalId[i]] = rowMeans(expMat[[2]][,which(colnames(expMat[[2]]) %in% rep_micro_names)])
          }else{
            exprdata_norep2[, uniq_Samples2$internalId[i]] = expMat[[2]][,which(colnames(expMat[[2]]) %in% rep_micro_names)]
          }
        }

        return(list(expVals = list(exprdata_norep1, exprdata_norep2), uniqSampInfo = list(uniq_Samples1, uniq_Samples2)))
      }
      
    }else{ # double coloured experiment

      ### change the colnames so they will be the same as in ExpInfoTable
      act_sdrf = sdrfFiles[eid]
      new_cols1 = vector(mode = 'character', dim(expMat[[1]])[2])
      new_cols2 = vector(mode = 'character', dim(expMat[[1]])[2])

      if(!(is.null(unlist(act_sdrf[[1]])))){

        for (z in 1:dim(expMat[[1]])[2]) {

          var1 = which(gsub("./", "", colnames(expMat[[1]])[z]) == act_sdrf[[1]]$Array.Data.File)
          var2 = which(act_sdrf[[1]]$Label == unique(dyes)[1])
          new_cols1[[z]] = as.character(act_sdrf[[1]]$Source.Name[intersect(var1, var2)])

          var3 = which(gsub("./", "", colnames(expMat[[2]])[z]) == act_sdrf[[1]]$Array.Data.File)
          var4 = which(act_sdrf[[1]]$Label == unique(dyes)[2])
          new_cols2[[z]] = as.character(act_sdrf[[1]]$Source.Name[intersect(var3, var4)])

        }

        colnames(expMat[[1]]) = new_cols1
        colnames(expMat[[2]]) = new_cols2
      }

      ### prepare matrix for results
      uniq_Samples = unique(ExpInfoTable[which(ExpInfoTable$SampleID %in% c(colnames(expMat[[1]]), colnames(expMat[[2]])) & ExpInfoTable$Experiment == eid), c(-2, -8)])
      Svector = paste0("S", as.character(1:dim(uniq_Samples)[1]))
      uniq_Samples$internalId = Svector

      exprdata_norep=array(0,dim=c(dim(expMat[[1]])[1],dim(uniq_Samples)[1]))
      row.names(exprdata_norep) = row.names(expMat[[1]])
      colnames(exprdata_norep) = Svector

      ### calculate mean value
      for (i in 1:dim(uniq_Samples)[1]) {

        # find rows in ExpInfoTable that are repetetive
        uniIDS = vector(mode='logical', length = dim(ExpInfoTable)[1])
        for(j in 1:dim(ExpInfoTable)[1]){

          uniIDS[j] = find_rows(ExpInfoTable[j, c(1, 3:6)], uniq_Samples[i, c(1:5)])

         }

        # find repetetive sample names
        cy3dye = which(ExpInfoTable$Experiment == eid & ExpInfoTable$Label == unique(dyes)[1])
        cy5dye = which(ExpInfoTable$Experiment == eid & ExpInfoTable$Label == unique(dyes)[2])

        #which(ExpInfoTable$Label[uniIDS] == dyes[1])])
        rep_micro_names_dye1 = colnames(expMat[[1]])[which(colnames(expMat[[1]]) %in% ExpInfoTable$SampleID[intersect(which(uniIDS==TRUE), cy3dye)])]
        rep_micro_names_dye2 = colnames(expMat[[2]])[which(colnames(expMat[[2]]) %in% ExpInfoTable$SampleID[intersect(which(uniIDS==TRUE), cy5dye)])]

        # take the data out
        the_Cy3 = expMat[[1]]
        the_Cy5 = expMat[[2]]

        # calculate mean value for repetetive samples
        if((length(rep_micro_names_dye1) + length(rep_micro_names_dye2))>1){
          preMatrix = cbind(the_Cy3[,which(colnames(the_Cy3) %in% rep_micro_names_dye1)], the_Cy5[,which(colnames(the_Cy5) %in% rep_micro_names_dye2)])
          exprdata_norep[, uniq_Samples$internalId[i]] = rowMeans(preMatrix)
        }else{
          exprdata_norep[, uniq_Samples$internalId[i]] = cbind(the_Cy3[,which(colnames(the_Cy3) %in% rep_micro_names_dye1)], the_Cy5[,which(colnames(the_Cy5) %in% rep_micro_names_dye2)])
        }
      }
    }
  }else{

    exprdata_norep = "No data"
    uniq_Samples = "No data"

  }
  return(list(expVals = exprdata_norep, uniqSampInfo = uniq_Samples))
}
