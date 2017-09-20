# na wejscie: tabelki z danymi o ekperymencie, character ze sciezka
# na wyjscie: wczytane dane tj. obiekt klasy uRNAList

#' @title Loading microarray miRNA data from a single experiment
#'
#' @param table A character matrix with information about experiment. The row names of it should be the names of files with data
#' and the column names should be 'Experiment', 'Cell line', 'Treatment', 'Time', 'Dose' exactly in this order.
#'
#' @param path A character indicating directory, where files with data are located.
#'
#' @return Function returns an object of uRNAList class.
#'
#' @description
#' \code{load_miRNA} function loads data for a single miRNA microarray experiment. The only supported platform is Agilent.
#'
#' @details
#' Example table structure (first argument) should look like this:
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
#' This function is designed to load data for a single experiment. For loading many experiments at once you could use \code{load_multi_miRNA}.
#'
#' @seealso
#' \code{\link{load_multi_miRNA}}
#'
#' @examples
#'\dontrun{
#' ### this example shows how to load data downloaded from ArrayExpress database
#'
#' # download data
#' datamiRNA = downloadAE('E-MTAB-5197', getwd())
#'
#' # prepare table as shown in details, load it and make sure it's a character matrix
#' path_to_table = system.file("inst/extdata", "miRNA_ex1.rds", package = "FindReference")
#' my_table = readRDS(path_to_table)
#'
#' # load data
#' loaded_data = load_miRNA(my_table, paste0(getwd(), 'E-MTAB-5197'))
#'}
#'
#' @rdname load_miRNA
#'
#' @importFrom AgiMicroRna readMicroRnaAFE
#'
#' @export


load_miRNA = function(table, path){

  #### create data.frame needed to load data
  targety = data.frame(FileName=row.names(table), Treatment = table[,'Treatment'], GErep = dim(table)[1])

  ### load data
  the_right_path = getwd()
  setwd(path)
  loaded_data = readMicroRnaAFE(targety)
  setwd(the_right_path)

  return(loaded_data)
}
