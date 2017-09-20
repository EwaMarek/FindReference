# na wejscie: wektor identyfikatorow ArrayExpress, sciezka do miejsca gdzie maja byc umieszczone pliki
# na wyjscie: zapis do pliku .rds, lista pobranych rzeczy (wyjscie z funkcji getAE)

#' @title Downloading microarray data from ArrayExpress database
#'
#' @description
#' \code{downloadAE} function downloads microarray data from ArrayExpress database basing on experiment identifiers. Then it identifies
#' microarray platform on which the experiment was performed using downloaded .sdrf file.
#'
#' @details
#' The output of this function is needed later in \code{load_multi_data} and \code{load_data} functions and that's why it is also
#' saved in dataAE.rds which could be easy loaded with \code{readRDS}.
#'
#' @param ExpIds A data frame with two columns named Experiment and RawData.
#' The first column should contain experiment identifiers from ArrayExpress database to download. In the RawData column there should be
#' a character "T" or "F" for each experiment, where "T" means that downloaded data should be raw type and "F" for processed type.
#'
#' @param path A character indicating directory where data should be downloaded and extracted.
#'
#' @return Function returns a list of lists where the first one contains names of downloaded files (the output from getAE function from ArrayExpress package)
#' and saves it as dataAE.rds. In the second list .sdrf files are stored in which there are information about samples and experiment conditions.
#' The third element of returned list is a data.frame with information about platform for each experiment.
#'
#' @examples
#' \dontrun{downloadAE(c("E-GEOD-65292", "E-TABM-90"), getwd())}
#'
#' \dontrun{experiments_ids = c("E-GEOD-65292", "E-TABM-90")
#' downloadAE(experiments_ids, getwd())}
#'
#' @rdname downloadAE
#'
#' @importFrom ArrayExpress getAE
#'
#' @export



downloadAE = function(ExpInfoTable, path){

  ################################################################
  ############  Find ArrayExpress Ids to download   ##############
  ################################################################

  ExpIds = unique(ExpInfoTable$Experiment)
  ExpIds = ExpIds[-grep(".", ExpIds, fixed = TRUE)]

  ################################################################
  ###########  Download data to specified directory   ############
  ################################################################

  dane = rep(list(list()), length(ExpIds))

  for (i in 1:length(dane)) {

    dir.create(paste(path,"/", ExpIds[i], sep = ""), showWarnings = TRUE,
               recursive = FALSE, mode = "0777")

    dane[[i]] = try(getAE(ExpIds[i], path = paste(path,"/", ExpIds[i], sep = ""),
                          type = 'raw', extract = TRUE))
  }

  if(length(which(is.character(dane)==TRUE))>0){
    warning(paste("Data with id:", as.character(ExpIds[which(is.character(dane)==TRUE)]), "could not be downloaded.", sep=' '))
  }

  saveRDS(dane, file = paste0(path,'/',"dataAE.rds"))

  ################################################################
  #####  Read in .sdrf file, find platform type and label   ######
  ################################################################
  sdrfFiles = rep(list(list()), length(dane))
  ExpInfoTable$Platform = rep(NA, dim(ExpInfoTable)[1])
  ExpInfoTable$Label = rep(NA, dim(ExpInfoTable)[1])

  for (i in 1:length(dane)) {
    sdrfFiles[[i]] = read.delim(paste(dane[[i]]$path, dane[[i]]$sdrf, sep = '/'))
  }

  names(sdrfFiles) = ExpIds

  platformDetails = readRDS(system.file("extdata", "platformDetails.rds", package = "FindReference"))
  platformDetails %>% mutate_if(is.factor, as.character) -> platformDetails

  for (i in 1:dim(ExpInfoTable)[1]) {
    if(ExpInfoTable$Experiment[i] %in% ExpIds){
      sdrf = sdrfFiles[[as.character(ExpInfoTable$Experiment[i])]]
      ExpInfoTable$Platform[i] =  platformDetails[match(unique(sdrf$Array.Design.REF), platformDetails$ID), 'Desc']
      ExpInfoTable$Label[i] = as.character(sdrf[which(as.character(sdrf$Source.Name) == as.character(ExpInfoTable$SampleID[i])), 'Label'][1])
    }
  }

  return(list(downloaded_data = dane, ExpInfoTable = ExpInfoTable, sdrfFiles = sdrfFiles))
}
