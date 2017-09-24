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
  ExpIds = ExpIds[-grep("/", ExpIds, fixed = TRUE)]

  ################################################################
  #####  Read in .sdrf file, find platform type and label   ######
  ################################################################
  sdrfFiles = rep(list(list()), length(ExpIds))
  ExpInfoTable$Platform = rep(NA, dim(ExpInfoTable)[1])
  ExpInfoTable$Label = rep(NA, dim(ExpInfoTable)[1])

  for (i in 1:length(ExpIds)) {
    sdrfFiles[[i]] = read.delim(paste0("https://www.ebi.ac.uk/arrayexpress/files/", ExpIds[i], "/", ExpIds[i], ".sdrf.txt"))
  }

  names(sdrfFiles) = ExpIds

  platformDetails = readRDS(system.file("extdata", "platformDetails.rds", package = "FindReference"))
  platformDetails %>% mutate_if(is.factor, as.character) -> platformDetails

  for (i in 1:dim(ExpInfoTable)[1]) {
    if(ExpInfoTable$Experiment[i] %in% ExpIds){
      sdrf = sdrfFiles[[as.character(ExpInfoTable$Experiment[i])]]
      ExpInfoTable$Platform[i] =  platformDetails[match(unique(sdrf$Array.Design.REF), platformDetails$ID), 'Desc']
      lab = as.character(sdrf[which(as.character(sdrf$Source.Name) == as.character(ExpInfoTable$SampleID[i])), 'Label'])

      # the if statement below needed in case that there are the same experimental conditions (on the same microarray)
      # but with different dyes
      if(length(lab)<2){
        ExpInfoTable$Label[i] = lab
      }else{
        if(lab[1] %in% ExpInfoTable[which(ExpInfoTable$SampleID == ExpInfoTable$SampleID[i]), 'Label'] ){
          ExpInfoTable$Label[i] = lab[2]
        }else{
          ExpInfoTable$Label[i] = lab[1]
        }
      }
    }
  }

  ################################################################
  ###########  Download data to specified directory   ############
  ################################################################

  dane = rep(list(list()), length(ExpIds))

  for (i in 1:length(dane)) {

    # get platform type
    full_platform_name = unique(ExpInfoTable[which(ExpInfoTable[, 'Experiment'] == ExpIds[i]), 'Platform'])

    if(grepl('Agilent', full_platform_name)){
      platform = "Agilent"
    }else if(grepl('Affymetrix', full_platform_name)){
      platform = "Affymetrix"
    }else{
      platform = full_platform_name
    }


    # check if the platform type is supported

    if(platform %in% c("Agilent", "Affymetrix")){

      dir.create(paste0(path,"/", ExpIds[i]), showWarnings = TRUE,
                 recursive = FALSE, mode = "0777")

      dane[[i]] = try(getAE(ExpIds[i], path = paste0(path,"/", ExpIds[i]),
                            type = 'raw', extract = TRUE))

      if(class(dane[[i]]) == 'try-error'){
        warning(paste0('Raw data with ', ExpIds[i], ' ArrayExpress ID could not be downloaded'))
      }
    }else{
      warning(paste0(platform, ' platform is not supported (', ExpIds[i],
                     ' ArrayExpress ID). For more information check details in downloadAE function help.'))
    }

  }

  saveRDS(dane, file = paste0(path,'/',"dataAE.rds"))


  return(list(downloaded_data = dane, ExpInfoTable = ExpInfoTable, sdrfFiles = sdrfFiles))
}
