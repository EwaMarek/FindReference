# na wejscie: lista z tabelkamii z danymi o ekperymencie, character vector ze sciezkami
# na wyjscie: wczytane dane tj. obiekty klasy uRNAList

#' @title Loading microarray miRNA data from many experiments
#'
#' @param tables A list of character matrices with information about experiment. The row names of matrix should be the names of files
#'  with data from one experiment and the column names should be 'Experiment', 'Cell line', 'Treatment', 'Time', 'Dose' exactly
#'  in this order.
#'
#' @param paths A character vector indicating directories, where files with data are located. It must be the same length as tables argument.
#'
#' @return Function returns a list of objects of uRNAList class.
#'
#' @description
#' \code{load_multi_miRNA} function loads data from many miRNA microarray experiments. The only supported platform is Agilent.
#'
#' @details
#' Example table structure (an element of first argument list) should look like this:
#'
#' \tabular{cccccc}{
#'                    \tab Experiment  \tab  Cell line  \tab  Treatment   \tab   Time   \tab   Dose   \cr
#'  array1.txt        \tab      Exp1    \tab   HCT116   \tab      C      \tab    0    \tab    0   \cr
#'  array2.txt        \tab      Exp1    \tab   HCT116   \tab      C      \tab     0    \tab     0    \cr
#'  array3.txt        \tab      Exp1    \tab   HCT116   \tab     IR      \tab   1    \tab    4    \cr
#'  array4.txt        \tab      Exp1    \tab    HCT116   \tab     IR      \tab    2    \tab     4    \cr
#'  array5.txt        \tab      Exp1    \tab    HCT116   \tab      IR     \tab     4    \tab    2
#'  }
#'
#'  Note that in treatment column it shoul be C for control arrays and IR for treated ones (even if the samples were treated with
#'  for example chemicals and not with ionizing radiation).
#'
#'  If you have downloaded data from ArrayExpress database you can find information needed for creating above in .sdrf files.
#'
#' @seealso
#' \code{\link{load_miRNA}}
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
#' }
#'
#' @rdname load_multi_miRNA
#'
#' @export

load_multi_miRNA = function(tables, paths){

  stopifnot(class(paths) == 'character')
  stopifnot(length(tables) == length(paths))

  loaded_data = rep(list(list()), length = length(paths))
  loaded_data = mapply(load_miRNA, tables, paths)

  return(loaded_data)
}
