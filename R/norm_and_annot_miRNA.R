# na wejscie: lista obiektow klasy uRNAList
# na wyjscie: tablice ekspresji po normalizacji i adnotacji

#' @title Normalization and annotation of microarray miRNA data
#'
#' @param loaded_data An object of class uRNAList or a list of objects of this type.
#'
#' @return Function returns an expression matrix (or list of them if argument was a list) with normalized and annotated expression values.
#'
#' @description
#' \code{norm_and_annot_miRNA} function loads normalizes and annotates data from miRNA microarray experiments performed on Agilent platform.
#'
#' @details
#' Data are rma background corrected and normalized. They are annotated with MIMAT ids.
#'
#' @examples
#' \dontrun{
#' ### download data from ArrayExpress database
#' AEids = c("E-MTAB-5197", "E-GEOD-59862")
#' datamiRNA = downloadAE(AEids, getwd())
#'
#' ### prepare tables as shown in details, load them and make sure they're character matrices
#' path_to_tables = system.file("inst/extdata", "miRNA_ex2.rds", package = "FindReference")
#' my_tables = readRDS(path_to_tables)
#'
#' ### load data
#'
#' # for data downloaded from ArrayExpress you can easily get the second argument by:
#' all_paths = sapply(datamiRNA, function(x){x$path})
#'
#' loaded_data = load_multi_miRNA(my_tables, all_paths)
#'
#' ### normalize and annotate data
#' norm_data = norm_and_annot_miRNA(loaded_data)
#' }
#'
#' @rdname norm_and_annot_miRNA
#'
#' @importFrom AgiMicroRna rmaMicroRna
#'
#' @export

norm_and_annot_miRNA = function(loaded_data){

  if(class(loaded_data)== "list"){

    exp_data = rep(list(list()), length(loaded_data))
    for(i in 1:length(loaded_data)){
      exp_data[[i]] = processmiRNA(loaded_data[[i]])
    }

  }else if(class(loaded_data) == "uRNAList"){

    exp_data = processmiRNA(loaded_data)

  }else{
    stop("Argument must be a class of 'list' or 'uRNAList'.")
  }

  return(exp_data)
}
