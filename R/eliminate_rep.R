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

eliminate_rep = function(expMat, eid, ExpInfoTable){

  if(class(expMat) != 'character'){

    ### check if experiment is single or double colour
    dbl_col = FALSE

    if(class(expMat) == 'list'){

      dyes = ExpInfoTable[which(ExpInfoTable$Experiment == eid), 'Label']

      if(length(unique(dyes)) == 2){
        dbl_col = TRUE
      }
    }

    # if the experiment is single coloured
    if(dbl_col == FALSE){

      ### change the colnames so they will be the same as in ExpInfoTable
      act_sdrf = sdrfFiles[eid]
      if(!(is.null(unlist(act_sdrf)))){
        new_cols = act_sdrf[[1]]$Source.Name[match(gsub("./", "", colnames(expMat)), act_sdrf[[1]]$Array.Data.File)]
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
      if(!(is.null(unlist(act_sdrf)))){
        new_cols1 = act_sdrf[[1]]$Source.Name[match(gsub("./", "", colnames(expMat[[1]])), act_sdrf[[1]]$Array.Data.File)]
        colnames(expMat[[1]]) = new_cols1

        new_cols2 = act_sdrf[[1]]$Source.Name[match(gsub("./", "", colnames(expMat[[2]])), act_sdrf[[1]]$Array.Data.File)]
        colnames(expMat[[2]]) = new_cols2
      }

      ### prepare matrix for results
      uniq_Samples = unique(ExpInfoTable[which(ExpInfoTable$SampleID %in% colnames(expMat[[1]]) & ExpInfoTable$Experiment == eid), c(-2, -8)])
      Svector = paste0("S", as.character(1:dim(uniq_Samples)[1]))
      uniq_Samples$internalId = Svector

      exprdata_norep=array(0,dim=c(dim(expMat[[1]])[1],dim(uniq_Samples)[1]))
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
        rep_micro_names_dye1 = colnames(expMat[[1]])[which(colnames(expMat[[1]]) %in% ExpInfoTable$SampleID[uniIDS])]
        rep_micro_names_dye2 = colnames(expMat[[2]])[which(colnames(expMat[[2]]) %in% ExpInfoTable$SampleID[uniIDS])]

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
