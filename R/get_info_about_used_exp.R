# na wejscie: list eksperymentow IR (nazwy kolumn w macierzy), czy ma generowac .txt, czy ma generowac .csv
# na wyjscie: pliki .txt luv .csv

#' @title Summarize information about experiments used to create ranking
#'
#' @param encoded_info List obtained with \code{create_ranking} function.
#'
#' @param txt Logical indicating if .txt files with summaries should be exported into txt files. Default value is TRUE.
#'
#' @param csv Logical indicating if .csv files with summaries should be exported into txt files. Default value is TRUE.
#'
#' @param path Character indicating the directory where files with summaries should be exported. Default value is current directory.
#'
#' @return Function doesn't return a value. Summaries of used cell lines, doses and times after treating are exported into .txt, .csv
#' files or both.
#'
#' @description
#' \code{get_info_about_used_exp} function summarizes conditions under which experiments used to create ranking where performed.
#' It lists and counts used cell lines, doses and times after treatment which could be useful to see for example if some conditions
#' where overrepresented.
#'
#' @seealso
#' \code{\link{create_ranking}}
#'
#' @examples
#' \dontrun{
#' ##### Create stability ranking for genes
#'
#' # download data from ArrayExpress database
#' to_download = c("E-GEOD-67309", "E-MTAB-966")
#' my_data = downloadAE(to_download, getwd())
#'
#' # load data
#' platforms = c("Affymetrix", "Agilent")
#' loaded_data = load_multi_data(my_data, platforms)
#'
#' # normalize and annotate
#' norm_data = multi_norm_and_annot(loaded_data$raw_expression_data, platforms)
#'
#' # prepare tables for rep_elim function as shown in details
#' path_to_tables = system.file("inst/extdata", "tables_ex3.rds", package = "FindReference")
#' my_tables = readRDS(path_to_tables)
#'
#' # eliminate replications and prepare object for create_ranking function
#' no_rep_data = rep_elim(norm_data, my_tables)
#'
#' # create ranking
#' gene_ranking = create_ranking(no_rep$noRepData, no_rep$uniqSamples, miRNA = FALSE)
#'
#' # summarise experimental conditions - export just .csv files in current directory
#' get_info_about_used_exp(gene_ranking, txt=FALSE, csv=TRUE)
#' }
#'
#' @rdname get_info_about_used_exp
#'
#' @importFrom utils write.table
#' @importFrom utils write.csv2
#'
#' @export



get_info_about_used_exp = function(encoded_info, txt=TRUE, csv=TRUE, path=getwd()){

  if(!is.null(encoded_info$samples)){
    encoded_info = encoded_info$samples
  }else{
    stop("First argument must be a list obtained from create_ranking function.")
  }


  decoded_info = list()
  decoded_info$cells = unlist(lapply(strsplit(encoded_info, ' '), function(x) x[2]))
  decoded_info$dose = unlist(lapply(strsplit(encoded_info, ' '), function(x) x[length(x)]))
  decoded_info$Time = unlist(lapply(strsplit(encoded_info, ' '), function(x) x[4]))
  decoded_info$uniCells = unique(decoded_info$cells)
  decoded_info$uniDose = unique(decoded_info$dose)
  decoded_info$uniTime = unique(decoded_info$Time)


  as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}

  decoded_info$how_many_cells = data.frame(decoded_info$uniCells, rep(0, length(decoded_info$uniCells)))
  colnames(decoded_info$how_many_cells)= c("Cells", "Counts")
  decoded_info$how_many_cells[, 2] = unlist(lapply(decoded_info$how_many_cells[, 1], function(x) length(which(decoded_info$cells == x))))
  #zmiana na dataframe'a i sortowanie
  decoded_info$how_many_cells = data.frame(decoded_info$how_many_cells)
  decoded_info$how_many_cells = decoded_info$how_many_cells[order(decoded_info$how_many_cells[,'Cells']),]

  decoded_info$how_many_doses = data.frame(decoded_info$uniDose, rep(0, length(decoded_info$uniDose)))
  colnames(decoded_info$how_many_doses)= c("Dose", "Counts")
  decoded_info$how_many_doses[, 2] = unlist(lapply(decoded_info$how_many_doses[, 1], function(x) length(which(decoded_info$dose == x))))
  #zmiana na dataframe'a i sortowanie
  decoded_info$how_many_doses = data.frame(decoded_info$how_many_doses)
  decoded_info$how_many_doses = decoded_info$how_many_doses[order(as.numeric.factor(decoded_info$how_many_doses[,'Dose'])),]
  colnames(decoded_info$how_many_doses)= c("Dose [Gy]", "Counts")

  decoded_info$how_many_times = data.frame(decoded_info$uniTime, rep(0, length(decoded_info$uniTime)))
  colnames(decoded_info$how_many_times)= c("Time", "Counts")
  decoded_info$how_many_times[, 2] = unlist(lapply(decoded_info$how_many_times[, 1], function(x) length(which(decoded_info$Time == x))))
  #zmiana na dataframe'a i sortowanie
  decoded_info$how_many_times = data.frame(decoded_info$how_many_times)
  decoded_info$how_many_times = decoded_info$how_many_times[order(as.numeric.factor(decoded_info$how_many_times[,'Time'])),]
  colnames(decoded_info$how_many_times)= c("Time [h]", "Counts")

  if(txt==TRUE){
    write.table(decoded_info$how_many_cells, paste(path, '/cells.txt', sep = ''), sep = ' ')
    write.table(decoded_info$how_many_doses, paste(path, '/doses.txt', sep = ''), sep = ' ')
    write.table(decoded_info$how_many_times, paste(path, '/time.txt', sep = ''), sep = ' ')
  }

  if(csv==TRUE){
    write.csv2(decoded_info$how_many_cells, paste(path, '/cells.csv', sep = ''), row.names=FALSE)
    write.csv2(decoded_info$how_many_doses, paste(path, '/doses.csv', sep = ''), row.names=FALSE)
    write.csv2(decoded_info$how_many_times, paste(path, '/time.csv', sep = ''), row.names=FALSE)
  }

}
