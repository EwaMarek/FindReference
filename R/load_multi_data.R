# na wejscie lista z danymi (z funkcji getAE), wektor z nazwami platformy,
# na wyjscie sdrfy oraz wczytane surowe dane

#' @title Loading microarray data from a list of experiments downloaded from ArrayExpress database
#'
#' @param dane A list of lists from downloadAE function.
#'
#' @param platforma A character vector indicating platform for each experiment. Supported platforms are 'Affymetrix' and 'Agilent'.
#'
#' @param donotread A list of character vectors indicating which microarrays should not be loaded for each experiment. Default is NA,
#' which means that all microarrays listed in dane$rawFiles will be loaded.
#' Note that list's length must be the same as the length of list 'dane'.
#'
#' @return Function returns a list with two elements. The first one is a list of objects specific for array platform with raw expression data.
#' The second element is a list of data frames with loaded .sdrf files for each experiment.
#'
#' @description
#' \code{load_multi_data} function loads data for a list of microarray experiments downloaded from ArrayExpress database.
#'
#' @details
#' This function is designed to read in microarray data from many experiments. For more detail see \code{load_data} function help.
#'
#' @seealso
#' \code{\link{load_data}}
#'
#' @examples
#' \dontrun{
#' ### load all microarrays for two experiments
#' to_download = c("E-GEOD-21066", "E-MTAB-966")
#' platforms = c("Affymetrix", "Agilent")
#'
#' dane = downloadAE(to_download, getwd())
#' loaded_experiments = load_multi_data(dane, platforms)
#'
#' ### do not load some microarrays from the first experiment
#' to_download = c("E-GEOD-21066", "E-MTAB-966")
#' platforms = c("Affymetrix", "Agilent")
#'
#' dane = downloadAE(to_download, getwd())
#'
#' # for second experiment all microarrays should be loaded
#' # -> the second element of the list is NA
#' unwanted_arrays = list(c("GSM526680.CEL", "GSM526756.CEL"), NA)
#'
#' # or (which could be useful when there are a lot of experiments)
#' unwanted_arrays = rep(list(list(NA)), length(dane))
#' unwanted_arrays[[1]] = c("GSM526680.CEL", "GSM526756.CEL")
#'
#' loaded_experiment = load_multi_data(dane, platforms, donotread = unwanted_arrays)
#'
#' ### when you downloaded data with downloadAE
#' ### but didn't assign the output into variable
#'
#' to_download = c("E-GEOD-21066", "E-MTAB-966")
#' platforms = c("Affymetrix", "Agilent")
#' downloadAE(to_download, getwd())
#'
#' # read in data saved by downloadAE function
#' dane = readRDS('dataAE.rds')
#' loaded_experiment = load__multi_data(dane, platforms)
#'
#' }
#'
#' @rdname load_multi_data
#'
#' @export



load_multi_data = function(dane, ExpInfoTable, sdrfFiles){

  raw_exp = rep(list(list()), length(dane))

  for(i in 1:length(dane)){
    raw_exp[[i]] = load_data(dane[[i]], ExpInfoTable, sdrfFiles[[i]])
    if(length(sub('.sdrf.txt', '', dane[[i]]$sdrf, fixed = TRUE))>0){
      names(raw_exp)[i] = sub('.sdrf.txt', '', dane[[i]]$sdrf, fixed = TRUE)
    }else{
      names(raw_exp)[i] = 'NA'
    }

  }

  if(class(raw_exp) == 'matrix'){
    raw_exp = list(raw_exp)
  }

  if(length(which(is.character(raw_exp)==TRUE))>0){
    warning(paste("Data with index", as.character(which(class(raw_exp) == 'try-error')), "could not be loaded.", sep=' '))
  }

  ### check if there are any user's raw data
  if(length(grep("/", ExpInfoTable$Experiment, fixed = TRUE)) > 0){

    users_exp = unique(ExpInfoTable[grep("/", ExpInfoTable$Experiment, fixed = TRUE), 'Experiment'])
    raw_exp2 = rep(list(list()), length(users_exp))

    for (i in 1:length(users_exp)){
      raw_exp2[[i]] = load_users_data(users_exp[i], ExpInfoTable)
      names(raw_exp2)[i] = as.character(users_exp[i])
    }

    if(class(raw_exp2) == 'matrix'){
      raw_exp2 = list(raw_exp2)
    }

  }

  return(raw_expression_data = c(raw_exp, raw_exp2))
}
