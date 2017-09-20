# na wejscie obiekt dane(z funkcji getAE), platforma,
# na wyjscie sdrf oraz wczytane surowe dane

#' @title Loading microarray data from single experiment downloaded from ArrayExpress database
#'
#' @param dane A list from downloadAE function.
#'
#' @param platforma Experiment's platform, could be 'Affymetrix' or 'Agilent'.
#'
#' @param donotread A character vector indicating which microarrays should not be loaded. Default is NA, which means that all
#' microarrays listed in dane$rawFiles will be loaded.
#'
#' @return Function returns a list with two elements. The first one is an object specific for array platform with raw expression data
#' (could be marrayRaw or EListRaw for Agilent platform and AffyBatch or ExonFeatureSet for Affymetrix platform).
#' The second element is a data frame with loaded .sdrf file.
#'
#' @description
#' \code{load_data} function loads data for a single microarray experiment downloaded from ArrayExpress database.
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
#' @rdname load_data
#'
#' @importFrom utils read.delim
#' @importFrom marray read.Agilent
#' @importFrom limma read.maimages
#' @importFrom affy ReadAffy
#' @importFrom oligo read.celfiles
#'
#' @export



load_data = function(dane, ExpInfoTable, sdrfFile){

  parent_directory = getwd()

  if(class(ExpInfoTable) == 'data.frame'){

    ExpInfoTable = unique(ExpInfoTable[which(ExpInfoTable[, 'Experiment'] == unlist(strsplit(dane$sdrf, '.sdrf'))[1]), 'Platform'])
  }

  if(grepl('Agilent', ExpInfoTable)){

    platforma = "Agilent"

  }else if(grepl('Affymetrix', ExpInfoTable)){

    platforma = "Affymetrix"

  }else{

    platforma = ExpInfoTable
    raw_exp = paste0(ExpInfoTable, " could not be loaded with load_data function. Only Agilent and Affymetrix platforms are supported.")
  }



  if(!(platforma %in% c("Agilent", "Affymetrix"))){
    warning(paste0(ExpInfoTable, " could not be loaded with load_data function. Only Agilent and Affymetrix platforms are supported."))
  }

  #######################   AGILENT   ####################

  if(platforma == "Agilent"){

    dyes = levels(sdrfFile[,'Label'])

    # ktore nazwy nie zawieraja miRNA
    without_micro = grep("miRNA", dane$rawFiles, invert = TRUE)
    do_wczytania = dane$rawFiles[without_micro]

    # ktore nazwy nie zawieraja tif
    without_tifs = grep(".tif", do_wczytania, invert = TRUE)
    do_wczytania = do_wczytania[without_tifs]

    # czy sa dwie platformy
    array_design =  unique(sdrfFile[, "Array.Design.REF"])

    ## ekperyment dwukolorowy
    if(length(dyes) == 2 && dyes[1] %in% c("Cy3", "Cy5") && dyes[2] %in% c("Cy3", "Cy5")){
      setwd(dane$path)

      # jedna platforma
      if(length(array_design)<2){
        raw_exp = try(read.Agilent(do_wczytania))
        if(is.character(raw_exp) == TRUE){
          raw_exp = try(read.Agilent(fnames = do_wczytania, path=getwd(), name.Gf = "gMedianSignal", name.Gb = "gBGMedianSignal", name.Rf = "rMedianSignal", name.Rb = "rBGMedianSignal",  info.id="ProbeUID"))

        }
        # dwie platformy
      }else{
        raw_exp = rep(list(list()), length(array_design))

        for (k in 1:length(array_design)) {

          do_wczytania = sdrfFile[which(sdrfFile[,"Array.Design.REF"]==array_design[k]),'Array.Data.File']
          raw_exp[[k]] = read.Agilent(do_wczytania)
        }

      }

      ## eksperyment jednokolorowy
    }else{

      # jedna platforma
      if(length(array_design)<2){
        raw_exp = read.maimages(do_wczytania, path=dane$path, source="agilent", green.only=TRUE)

        # dwie platformy
      }else{
        raw_exp = rep(list(list()), length(array_design))

        for (k in 1:length(array_design)) {
          do_wczytania = sdrfFile[which(sdrfFile[,"Array.Design.REF"]==array_design[k]),'Array.Data.File']
          raw_exp[[k]] = read.maimages(do_wczytania, path=dane$path, source="agilent", green.only=TRUE )
        }
      }
    }
  }


  #######################   AFFYMETRIX   ####################
  if(platforma == "Affymetrix"){

    # ktore nazwy nie zawieraja miRNA
    without_micro = grep("miRNA", dane$rawFiles, invert = TRUE)
    do_wczytania = dane$rawFiles[without_micro]

    # ktore nazwy nie zawieraja txt
    without_txt = grep(".txt", do_wczytania, invert = TRUE)
    do_wczytania = do_wczytania[without_txt]

    dyes = levels(sdrfFile[,'Label'])

    # czy sa dwie platformy
    array_design =  unique(sdrfFile[intersect(without_micro, without_txt), "Array.Design.REF"])
    setwd(dane$path)

    # jedna platforma
    if(length(array_design)<2){

      raw_exp = try(ReadAffy(filenames = do_wczytania, celfile.path = dane$path), silent = TRUE)
      if(is.character(raw_exp)==TRUE){raw_exp = read.celfiles(do_wczytania)}

      # dwie platformy
    }else{
      raw_exp = rep(list(list()), length(array_design))
      for (k in 1:length(array_design)) {
        do_wczytania = sdrfFile[which(sdrfFile[,"Array.Design.REF"]==array_design[k]),'Array.Data.File']
        do_wczytania = as.character(do_wczytania)

        raw_exp[[k]] = try(ReadAffy(filenames = do_wczytania, celfile.path = dane$path), silent = TRUE)
        if(is.character(raw_exp)==TRUE){raw_exp = read.celfiles(do_wczytania)}
      }

    }
  }

  setwd(parent_directory)

  return(raw_exp)
}
