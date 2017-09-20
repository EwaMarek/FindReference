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
add_experiments = function(nData, procData, ExpInfoTable, EG2SYM){

  #EG2SYM = as.list(org.Hs.egSYMBOL)
  ################################################################
  ############  Check correctness of added data ##############
  ################################################################

  if(class(procData) == 'matrix'){

    if(sum(colnames(procData) %in% ExpInfoTable$SampleID)<2){
      stop("At least two colnames of processed data must be in SampleID column in table with information about experiment
           (third argument of add_experiment function).")
    }

    procData = procData[-which(is.na(rownames(procData)) == TRUE), ]
    procData = procData[-which(rownames(procData) == ''),]

    if(sum(rownames(procData) %in% names(EG2SYM)) == 0   &  sum(rownames(procData) %in% EG2SYM) == 0){

      stop("Rownames of processed data must be hgnc gene symbol or EntrezID.")

    }else if(sum(rownames(procData) %in% EG2SYM) != 0){

      entrezID = AnnotationDbi::mget(rownames(procData), revmap(org.Hs.egSYMBOL), ifnotfound = NA)
      rownames(procData) = entrezID
    }

  }else if(class(procData) == 'list'){

    for(i in 1:length(procData)){
      if(sum(colnames(procData[[i]]) %in% ExpInfoTable$SampleID)<2){
        stop("At least two colnames of processed data must be in SampleID column in table with information about experiment
           (third argument of add_experiment function).")
      }

      procData[[i]] = procData[[i]][-which(is.na(rownames(procData[[i]])) == TRUE), ]
      procData[[i]] = procData[[i]][-which(rownames(procData[[i]]) == ''),]

      if(sum(rownames(procData[[i]]) %in% names(EG2SYM)) == 0   &  sum(rownames(procData[[i]]) %in% EG2SYM) == 0){

        stop("Rownames of processed data must be hgnc gene symbol or EntrezID.")

      }else if(sum(rownames(procData[[i]]) %in% EG2SYM) != 0){

        entrezID = AnnotationDbi::mget(rownames(procData[[i]]), revmap(org.Hs.egSYMBOL), ifnotfound = NA)
        rownames(procData[[i]]) = entrezID
      }
    }
  }

  ################################################################
  ####################  Add new processed data ###################
  ################################################################

  allProcData = c(nData, procData)
}
