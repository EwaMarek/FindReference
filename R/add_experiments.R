#' @title Adding self-processed data to existing ones
#'
#' @param nData A list of objects returned by \code{multi_norm_and_annot}
#'
#' @param procData Matrix or list of matrices
#'
#' @param ExpInfoTable
#'
#' @return List of all matrices with normalized expression values
#'
#' @description
#' \code{add_experiments} function enables user to add processed data
#'
#' @details
#' Blablabla
#'
#' @examples
#' \dontrun{
#' ### download, load, normalize and annotate data from ArrayExpress
#'
#' # download
#' to_download = c("E-GEOD-21066", "E-MTAB-966")
#' data = downloadAE(to_download, getwd())
#'
#' # load
#' platforms = c("Affymetrix", "Agilent")
#' loaded_experiments = load_multi_data(data, platforms)
#'
#' # process
#' processed_data = multi_norm_and_annot(loaded_experiments[[1]], platforms)
#' }
#'
#' @rdname add_experiments
#'
#' @export
#'
add_experiments = function(nData, procData, ExpInfoTable){

  EG2SYM = AnnotationDbi::as.list(org.Hs.egSYMBOL)
  ################################################################
  ############  Check correctness of added data ##############
  ################################################################

  if(class(procData) == 'matrix'){

    if(sum(colnames(procData) %in% ExpInfoTable$SampleID)<2){
      stop("At least two colnames of processed data must be in SampleID column in table with information about experiment
           (third argument of add_experiment function).")
    }

    # delete rows with no row.names or row.names equal to NA
    no_rownames = which(is.na(rownames(procData)) == TRUE | rownames(procData) == '')
    if(length(no_rownames)>0){
      procData = procData[-no_rownames, ]
    }


    if(sum(rownames(procData) %in% names(EG2SYM)) == 0   &  sum(rownames(procData) %in% EG2SYM) == 0){

      stop("Rownames of processed data must be hgnc gene symbol or EntrezID.")

    }else if(sum(rownames(procData) %in% EG2SYM) != 0){

      entrezID = AnnotationDbi::mget(rownames(procData), revmap(org.Hs.egSYMBOL), ifnotfound = NA)
      rownames(procData) = entrezID

      # delete rows without EntrezID symbol
      rows_to_remove = which(is.na(rownames(procData)) == TRUE)

      if(length(rows_to_remove)>0){
        procData = procData[-rows_to_remove, ]
      }


    }

    name = ExpInfoTable[intersect(which(ExpInfoTable$SampleID %in% colnames(procData)),
                              grep('.', ExpInfoTable[,'Experiment'], fixed = TRUE)), 'Experiment'][1]
    allProcData = c(nData, list(procData))
    names(allProcData)[length(allProcData)] = as.character(name)

  }else if(class(procData) == 'list'){

    for(i in 1:length(procData)){
      if(sum(colnames(procData[[i]]) %in% ExpInfoTable$SampleID)<2){
        stop("At least two colnames of processed data must be in SampleID column in table with information about experiment
           (third argument of add_experiment function).")
      }


      # delete rows with no row.names or row.names equal to NA
      no_rownames = which(is.na(rownames(procData[[i]])) == TRUE | rownames(procData[[i]]) == '')
      if(length(no_rownames)>0){
        procData[[i]] = procData[[i]][-no_rownames, ]
      }

      if(sum(rownames(procData[[i]]) %in% names(EG2SYM)) == 0   &  sum(rownames(procData[[i]]) %in% EG2SYM) == 0){

        stop("Rownames of processed data must be hgnc gene symbol or EntrezID.")

      }else if(sum(rownames(procData[[i]]) %in% EG2SYM) != 0){

        entrezID = AnnotationDbi::mget(rownames(procData[[i]]), revmap(org.Hs.egSYMBOL), ifnotfound = NA)
        rownames(procData[[i]]) = entrezID

        # delete rows without EntrezID symbol
        rows_to_remove = which(is.na(rownames(procData[[i]])) == TRUE)

        if(length(rows_to_remove)>0){
          procData[[i]] = procData[[i]][-rows_to_remove, ]
          names(procData)[i] = as.character(ExpInfoTable[intersect(which(ExpInfoTable$SampleID %in% colnames(procData[[i]])),
                                                                   grep('.', ExpInfoTable[,'Experiment'], fixed = TRUE)), 'Experiment'][1])
        }

      }
    }

    allProcData = c(nData, procData)
  }

  return(allProcData)
}
