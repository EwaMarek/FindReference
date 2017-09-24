#' @title Loading microarray data from single experiment
#'
#' @param directory A directory to files from single experiment
#'
#' @return Function returns a list with two elements. The first one is an object specific for array platform with raw expression data
#' (could be marrayRaw or EListRaw for Agilent platform and AffyBatch or ExonFeatureSet for Affymetrix platform).
#' The second element is a data frame with loaded .sdrf file.
#'
#' @description
#' \code{load_users_data} function loads data for a single microarray experiment downloaded from ArrayExpress database.
#'
#' @details
#' This function is designed to load data for a single experiment. For loading many experiments at once you could use \code{load_multi_data}.
#'
#' Note that only Agilent (one and double colour) and Affymetrix platforms are supported.
#'
#' Loaded .sdrf file (second element of the output) could be useful for preparing table with experiment information needed for \code{rep_elim}
#' or \code{multi_rep_elim} functions.
#'
#' @seealso
#' \code{\link{load_multi_data}}
#'
#' @examples
#' \dontrun{
#' ### load all microarrays
#' dane = downloadAE("E-GEOD-21066", getwd())
#' loaded_experiment = load_data(dane, "Affymetrix")
#'
#' ### exclude some microarrays from loading
#' dane = downloadAE("E-GEOD-21066", getwd())
#' micrroarrays_under_conditions_imnot_interested = c("GSM526680.CEL", "GSM526756.CEL")
#' loaded_experiment = load_data(dane, "Affymetrix",
#'                          donotread = micrroarrays_under_conditions_imnot_interested)
#'
#' ### when you downloaded data with downloadAE
#' ### but didn't assign the output into variable
#' downloadAE("E-GEOD-21066", getwd())
#'
#' # read in data saved by downloadAE function
#' dane = readRDS('dataAE.rds')
#' loaded_experiment = load_data(dane, "Affymetrix")
#'
#' }
#'
#' @rdname load_users_data
#'
#' @importFrom utils read.delim
#' @importFrom marray read.Agilent
#' @importFrom limma read.maimages
#' @importFrom affy ReadAffy
#' @importFrom oligo read.celfiles
#'




load_users_data = function(directory, ExpInfoTable){

  parent_directory = getwd()
  directory = as.character(directory)

  ################################################################
  #####################  Check files type  #######################
  ################################################################
  all_files = list.files(path = directory)

  if(length(grepl(".CEL", all_files)) > 0){ # Affymetrix platform
    eksp_files = list.celfiles(directory)

    setwd(directory)
    raw_exp = try(ReadAffy(filenames = eksp_files, celfile.path = directory), silent = TRUE)

    if(class(raw_exp) == 'try-error'){
      raw_exp = try(read.celfiles(eksp_files))
    }

  }else{ # Agilent platform
    eksp_files = list.files(directory, pattern = '.txt')

    dyes = ExpInfoTable[which(ExpInfoTable$Experiment == directory), 'Label']

    if(length(unique(dyes)) < 2){ # if one colour experiment

      raw_exp = read.maimages(eksp_files, path=directory, source="agilent", green.only=TRUE)

    }else{ # if double colour experiment

      raw_exp = try(read.Agilent(eksp_files))
      if(class(raw_exp) == 'try-error'){
        raw_exp = try(read.Agilent(fnames = eksp_files, path=directory, name.Gf = "gMedianSignal", name.Gb = "gBGMedianSignal", name.Rf = "rMedianSignal", name.Rb = "rBGMedianSignal",  info.id="ProbeUID"))

      }
    }
  }


  setwd(parent_directory)
  return(raw_exp)
}
