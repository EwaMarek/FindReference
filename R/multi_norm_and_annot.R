# na wejscie lista obiektow do normalizacji oraz wektor z nazwami platform
# na wyjscie tablice ze znormalizowanymi wartosciami ekspresji oraz odpowiednia anotacja


#' @title Normalization and annotation of microarray data from many experiments
#'
#' @param dataList A list of objects of class marrayRaw, EListRaw, AffyBatch, ExonFeatureSet or list
#' (see details) with raw expression data.
#'
#' @param platformVector A character vector indicating experiments' platform (could be 'Affymetrix' or 'Agilent').
#'
#' @return A list of expression matrices with normalized and annotated expression values.
#'
#' @description
#' \code{multi_norm_and_annot} function normalizes and annotates data from many microarray experiments.
#'
#' @details
#' An element of a dataList (first argument) could be also a list due to make it possible to process data from an experiment where there are
#' few array types within one platform (Affymetrix or Agilent). For details about normalization and annotation used
#' in \code{multi_norm_and_annot} see \code{norm_and_annot_by_platform} help.
#'
#' @seealso
#' \code{\link{norm_and_annot_by_platform}}
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
#' @rdname multi_norm_and_annot
#'
#' @export


multi_norm_and_annot = function(dataList, ExpInfoTable){

  platformVector = unique(ExpInfoTable$Platform)
  goOut = rep(list(list()), length(dataList))

  for(i in 1:length(dataList)){
    goOut[[i]] = norm_and_annot_by_platform(dataList[[i]], platformVector[i])
  }

  if(length(which(is.character(goOut)==TRUE))>0){
    warning(paste("Data with index", as.character(which(is.character(goOut)==TRUE)), "could not be normalized.", sep=' '))
  }

  return(goOut)

}
